/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * MULTEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MULTEM is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef HOST_FUNCTIONS_H
#define HOST_FUNCTIONS_H

#include <thread>
#include <type_traits>
#include <algorithm>
#include <numeric>
#include <random>

#include <fftw3.h>
#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "fft2.cuh"
#include "stream.cuh"
#include "lapack.hpp"

#include "host_device_functions.cuh"

namespace multem
{

	namespace host_detail
	{
		template <class TFn, class ...TArg> 
		void matrix_iter(const Range &range, TFn &fn, TArg &...arg)
		{
			for(auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for(auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					fn(ix, iy, arg...);
				}
			}
		}

		template <class TFn, class ...TArg> 
		void vector_iter(const Range &range, TFn &fn, TArg &...arg)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				fn(ixy, arg...);
			}
		}

		template<class T>
		T atom_cost_function(const Grid<T> &grid, const Atom_Ip<T> &atom_Ip, rVector<T> M_i)
		{
			T sum = 0;

			for(auto ix0 = 0; ix0 < atom_Ip.ixn; ix0++)
			{
				for(auto iy0 = 0; iy0 < atom_Ip.iyn; iy0++)
				{
					int ix = ix0 + atom_Ip.ix0;
					int iy = iy0 + atom_Ip.iy0;

					T R2 = grid.R2(ix, iy, atom_Ip.x, atom_Ip.y);
					if (R2 < atom_Ip.R2_max)
					{
						int ixy = grid.ind_col(ix, iy);
						T M = M_i[ixy];
						T V = host_device_detail::eval_cubic_poly<T>(R2, atom_Ip);
						sum += (V-2*M)*V;
					}
				}
			}

			return sum;
		}

		template<class T>
		void subtract_atom(Stream<e_host> &stream, const Grid<T> &grid, const Atom_Ip<T> &atom_Ip, rVector<T> M_i)
		{
			for(auto ix0 = 0; ix0 < atom_Ip.ixn; ix0++)
			{
				int iyc = 0;
				for(auto iy0 = 0; iy0 < atom_Ip.iyn; iy0++)
				{
					int ix = ix0 + atom_Ip.ix0;
					int iy = iy0 + atom_Ip.iy0;

					T R2 = grid.R2(ix, iy, atom_Ip.x, atom_Ip.y);
					if (R2 < atom_Ip.R2_max)
					{
						int ixy = grid.ind_col(ix, iy);

						T V = host_device_detail::eval_cubic_poly<T>(R2, atom_Ip);

						atom_Ip.iv[iyc] = ixy;
						atom_Ip.v[iyc] = V;
						iyc++;
					}
				}

				stream.stream_mutex.lock();
				for(auto iy0 = 0; iy0 < iyc; iy0++)
				{
					M_i.V[atom_Ip.iv[iy0]] += atom_Ip.v[iy0];
				}
				stream.stream_mutex.unlock();
			}
		}

		// Linear projected potential: V and zV
		template<ePotential_Type potential_type, int charge, class T> 
		void linear_Vz(const Q1<T, e_host> &qz, Atom_Vp<T> &atom_Vp)
		{	
			for(auto iR =0; iR < c_nR; iR++)
			{
				T R2 = atom_Vp.R2[iR];
				T V0s = 0;
				T dV0s = 0;

				T a = (atom_Vp.split)?(-atom_Vp.z0h):(atom_Vp.zeh-atom_Vp.z0h);
				T b = (atom_Vp.split)?(atom_Vp.z0h):(atom_Vp.zeh+atom_Vp.z0h);
				for(auto ix = 0; ix < qz.size(); ix++)
				{
					T z = a*qz.x[ix] + b;
					T r = sqrt(z*z + R2);
					T V, dVir;
					Vr_dVrir<potential_type, charge, T>(r, atom_Vp.cl, atom_Vp.cnl, a*qz.w[ix], V, dVir);
					V0s += V;
					dV0s += dVir;
				}

				if (atom_Vp.split)
				{
					T a = atom_Vp.zeh;
					T b = atom_Vp.zeh;
					for(auto ix = 0; ix < qz.size(); ix++)
					{
						T z = a*qz.x[ix] + b;
						T r = sqrt(z*z + R2);
						T V, dVir;
						Vr_dVrir<potential_type, charge, T>(r, atom_Vp.cl, atom_Vp.cnl, a*qz.w[ix], V, dVir);
						V0s += V;
						dV0s += dVir;
					}
				}

				atom_Vp.c0[iR] = V0s; 		// V_0
				atom_Vp.c1[iR] = 0.5*dV0s; 	// dR2V0
			}
		}

		// Get Local interpolation coefficients
		template<class T> 
		void cubic_poly_coef(Atom_Vp<T> &atom_Vp)
		{
			for(auto iR = 0; iR < c_nR-1; iR++)
			{
				host_device_detail::cubic_poly_coef(iR, atom_Vp);
			}
		}

		// Cubic polynomial evaluation
		template<class T> 
		void eval_cubic_poly(Stream<e_host> &stream, const Grid<T> &grid, const Atom_Vp<T> &atom_Vp, const rVector<T> &V0g)
		{	
			for(auto ix0 = 0; ix0 < atom_Vp.ixn; ix0++)
			{
				int iyc = 0;
				for(auto iy0 = 0; iy0 < atom_Vp.iyn; iy0++)
				{
					int ix = ix0 + atom_Vp.ix0;
					int iy = iy0 + atom_Vp.iy0;

					T R2 = grid.R2(ix, iy, atom_Vp.x, atom_Vp.y);
					if (R2 < atom_Vp.R2_max)
					{
						int ixy;
						T V = host_device_detail::eval_cubic_poly(ix, iy, grid, R2, atom_Vp, ixy);
						
						atom_Vp.iv[iyc] = ixy;
						atom_Vp.v[iyc] = V;
						iyc++;
					}
				}

				stream.stream_mutex.lock();
				for(auto iy0 = 0; iy0 < iyc; iy0++)
				{
					V0g.V[atom_Vp.iv[iy0]] += atom_Vp.v[iy0];
				}
				stream.stream_mutex.unlock();
			}
		}

		template<class TGrid>
		Value_type<TGrid> Lorentz_factor(Stream<e_host> &stream, TGrid &grid, EELS<Value_type<TGrid>> &eels)
		{
			using value_type_r = Value_type<TGrid>;
			value_type_r sum = 0;

			auto thr_Lorentz_factor = [&](const Range &range)
			{
				value_type_r sum_partial = 0;
				matrix_iter(range, host_device_detail::Lorentz_factor<TGrid>, grid, eels.gc2, eels.ge2, sum_partial);

				stream.stream_mutex.lock();
				sum += sum_partial;
				stream.stream_mutex.unlock();
			};

			stream.set_n_act_stream(grid.nx);
			stream.set_grid(grid.nx, grid.ny);
			stream.exec(thr_Lorentz_factor);

			return sqrt(eels.occ)/sum;
		}
	} // host_detail

	/***************************************************************************/
	/***************************************************************************/

	// get index (with typ =0: bottom index for equal values and typ =1: upper index for equal values)
	int getIndex(int ixmin, int ixmax, double *x, int typ, double x0)
	{
		int ixmid; 
		switch(typ)
		{
			case 0:
			{
				do{
					ixmid = (ixmin + ixmax)>>1; 	// divide by 2
					if(x0 <= x[ixmid]) 
						ixmax = ixmid;
					else 
						ixmin = ixmid;
				}while ((ixmax-ixmin)>1);
			}
			break;
			case 1:
			{
				do{
					ixmid = (ixmin + ixmax)>>1; 	// divide by 2
					if(x0 < x[ixmid]) 
						ixmax = ixmid;
					else 
						ixmin = ixmid;
				}while ((ixmax-ixmin)>1);
			}
			break;
		}

		if(x0 == x[ixmax])
			return ixmax;
		else
			return ixmin;
	}

	template<class T>
	void get_bn(const T &R, const int &nR, const T &dR, const T &R_max, const bool &pbc, int &iR0, int &iRn)
	{
		int iR_0 = static_cast<int>(floor((R-R_max)/dR));
		int iR_e = static_cast<int>(ceil((R+R_max)/dR));

		if(!pbc)
		{
			auto set_Bound = [](const int &i, const int &n)->int{ return (i<0)?0:((i>=n)?n-1:i); };
			iR_0 = set_Bound(iR_0, nR);
			iR_e = set_Bound(iR_e, nR);
		}

		iR0 = iR_0;
		iRn = (iR_0 == iR_e)?0:iR_e-iR_0+1;
	}

	// symmetric coordinates(Fourier space Coordinates)
	inline
	int FSC(const int &i, const int &nh, bool shift =false)
	{
		return (shift)?((i<nh)?(i):(i-2*nh)):(i-nh);
	}

	// symmetric coordinates(Fourier space Coordinates)
	inline
	int RSC(const int &i, const int &nh, bool shift = false)
	{
		return (shift)?((i<nh)?(i+nh):(i-nh)):i;
	}

	// get one dimensional Hanning_Filter
	template<class T>
	void filter_Hanning_1D(int nx, T dx, T k, bool shift, T *fI)
	{	
		int nxh = nx/2;
		auto cx = c_2Pi/(dx*nx);

		for(auto ix = 0; ix < nx; ix++)
		{
			auto Rx = RSC(ix, nxh, shift)*dx; 
			auto f = 0.5*(1.0-cos(cx*Rx));
			fI[ix] = (f>k)?1.0:f/k;
		}
	}

	// get two dimensional Hanning_Filter
	template<class T>
	void filter_Hanning_2D(int ny, int nx, T dx, T dy, T k, bool shift, T *fI)
	{	
		int nxh = nx/2, nyh = ny/2;
		auto cx = c_2Pi/(dx*nx);
		auto cy = c_2Pi/(dy*ny);

		Vector<T, e_host> fx(nx);
		Vector<T, e_host> fy(ny);

		for(auto ix = 0; ix < fx.size(); ix++)
		{
			auto Rx = RSC(ix, nxh, shift)*dx; 
			fx[ix] = 0.5*(1.0-cos(cx*Rx));
		}

		for(auto iy = 0; iy < fy.size(); iy++)
		{
			auto Ry = RSC(iy, nyh, shift)*dy;
			fy[iy] = 0.5*(1.0-cos(cy*Ry));
		}

		for(auto ix = 0; ix < fx.size(); ix++)
		{
			for(auto iy = 0; iy < fy.size(); iy++)
			{
				auto f = fx[ix]*fy[iy];
				fI[ix*ny+iy] = (f>k)?1.0:f/k;
			}
		}
	}


	// get one dimensional Gaussian_Filter
	template<class T>
	void filter_Gaussian_1D(int nx, T dx, T Sigma, bool shift, T *fI)
	{	
		int nxh = nx/2;
		auto cx = 0.5/(Sigma*Sigma);

		for(auto ix = 0; ix < nx; ix++)
		{
			auto Rx = FSC(ix, nxh, shift)*dx; 
			fI[ix] = exp(-cx*Rx*Rx);
		}
	}

	// get two dimensional Gaussian_Filter
	template<class T>
	void filter_Gaussian_2D(int ny, int nx, T dx, T dy, T Sigma, bool shift, T *fI)
	{	
		int nxh = nx/2, nyh = ny/2;
		T Rx, Ry;
		T cx = 0.5/(Sigma*Sigma);
		T cy = 0.5/(Sigma*Sigma);

		Vector<T, e_host> fx(nx);
		Vector<T, e_host> fy(ny);

		for(auto ix = 0; ix < fx.size(); ix++)
		{
			Rx = FSC(ix, nxh, shift)*dx; 
			fx[ix] = exp(-cx*Rx*Rx);
		}

		for(auto iy = 0; iy < fy.size(); iy++)
		{
			Ry = FSC(iy, nyh, shift)*dy; 
			fy[iy] = exp(-cy*Ry*Ry);
		}

		for(auto ix = 0; ix < fx.size(); ix++)
		{
			for(auto iy = 0; iy < fy.size(); iy++)
			{
				fI[ix*ny+iy] = fx[ix]*fy[iy];
			}
		}
	}


	// get one dimensional Butterworth_Filter
	template<class T>
	void filter_Butterworth_1D(int nx, T dx, T Radius, int n, bool shift, T *fI)
	{	
		int nxh = nx/2;
		auto R02 = pow(Radius, 2);

		for(auto ix = 0; ix < nx; ix++)
		{
			auto Rx = FSC(ix, nxh, shift)*dx; 
			fI[ix] = 1.0/(1.0+pow(Rx*Rx/R02, n));
		}
	}

	// get two dimensional Butterworth_Filter
	template<class T>
	void filter_Butterworth_2D(int ny, int nx, T dx, T dy, T Radius, int n, bool shift, T *fI)
	{	
		int nxh = nx/2, nyh = ny/2;
		auto R02 = pow(Radius, 2);

		for(auto ix = 0; ix < nx; ix++)
		{
			for(auto iy = 0; iy < ny; iy++)
			{
				auto Rx = FSC(ix, nxh, shift)*dx; 
				auto Ry = FSC(iy, nyh, shift)*dy;
				auto R2 = Rx*Rx + Ry*Ry;
				fI[ix*ny+iy] = 1.0/(1.0+pow(R2/R02, n));
			}
		}
	}


	// get two dimensional radial distribution for regular grid
	void radial_distribution_2D(int nR, double *R, double *fR, int nRl, double *Rl, double *rl, double *frl, double *cfrl, bool reg, int typ)
	{	
		double Rlmin = Rl[0], Rlmax = Rl[nRl-1], dRl = Rl[1]-Rl[0];
 
		for(auto i = 0; i < nRl-1; i++)
		{
			rl[i] = 0.5*(Rl[i]+Rl[i+1]);
			frl[i] = 0.0;
			cfrl[i] = 0.0;
		}

		int j;
		for(auto i = 0; i < nR; i++)
		{
			if((Rlmin <= R[i])&&(R[i]<Rlmax))
			{
				j = (reg)?(int)floor((R[i]-Rlmin)/dRl):getIndex(0, nRl-1, Rl, 0, R[i]);
				frl[j] += fR[i];
				cfrl[j] += 1.0;
			}
		}

		if(typ == 0 )
			for(auto i = 0; i < nRl-1; i++)
			{
				if(cfrl[i]>0)
					frl[i] /= cfrl[i];
			}
	}

	template<class TGrid, class TVector>
	void radial_distribution_2D(TGrid &grid, TVector &Im, Value_type<TVector> x_i, 
	Value_type<TVector> y_i, Value_type<TVector> radius_i, TVector &rl, TVector &frl)
	{
		using T = Value_type<TVector>;

		T dR = grid.dRx;
		int nrl = static_cast<int>(floor(radius_i/dR+0.5));
		T R_max = dR*nrl;

		int ix0, ixe;
		get_bn(x_i, grid.nx, grid.dRx, R_max, grid.pbc_xy, ix0, ixe);
		ixe += ix0;

		int iy0, iye;
		get_bn(y_i, grid.ny, grid.dRy, R_max, grid.pbc_xy, iy0, iye);
		iye += iy0;

		// get radial distribution
		rl.resize(nrl);
		frl.resize(nrl);
		TVector cfrl(nrl);

		for(auto ir = 0; ir < rl.size(); ir++)
		{
			rl[ir] = grid.Rx(ir);
			frl[ir] = 0.0;
			cfrl[ir] = 0.0;
		}

		for(auto ix = ix0; ix < ixe; ix++)
		{
			for(auto iy = iy0; iy < iye; iy++)
			{
				T R = grid.R(ix, iy, x_i, y_i);
				if(R < R_max)
				{
					int ixy = grid.ind_col(ix, iy);
					auto ir = static_cast<int>(floor(R/dR));
					frl[ir] += Im[ixy];
					cfrl[ir] += 1.0;
				}
			}
		}

		for(auto ir = 0; ir < rl.size(); ir++)
		{
			frl[ir] /= ::fmax(1.0, cfrl[ir]);
		}
	}

	// get two dimensional radial distribution for regular grid
	void getCumRadDist_2D(int ny, int nx, bool shift, double *fI, int nr, double *r, double *fIr)
	{	
		int idx, nxh = nx/2, nyh = ny/2;
		Vector<double, e_host> cfIr(nr, 0.0);
		std::fill(fIr, fIr + nr, 0.0);

		double Rx, Ry, R;
		for(auto ix = 0; ix < nx; ix++)
		{
			Rx = FSC(ix, nxh, shift);
			for(auto iy = 0; iy < ny; iy++)
			{
				Ry = FSC(iy, nyh, shift);
				R = sqrt(Rx*Rx + Ry*Ry);
				if((0 <= R)&&(R<nr))
				{
					idx = (int)floor(R);
					fIr[idx] += fI[ix*ny+iy];
					cfIr[idx] += 1.0;
				}
			}
		}

		r[0] = 0; 
		fIr[0] /= cfIr[0];
		for(auto i = 0; i < nr; i++)
		{
			r[i] = i;
			if(cfIr[i]>0) 
				fIr[i] /= cfIr[i];
			fIr[i] += fIr[i-1];
		}
	}

	// mean smooth
	template<class TVector>
	void smooth(Stream<e_host> &stream, TVector &f_i, Value_type<TVector> nkr, TVector &f_o)
	{
		using T = Value_type<TVector>;

		int nf_i = f_i.size();
		int nk0 = -nkr;
		int nke = nkr+1;

		Vector<T, e_host> fv(nf_i);
		auto krn_mean = [&](const int &ix_i, TVector &f_i, Vector<T, e_host> &f_o)
		{
			int ix0 = max(ix_i+nk0, 0);
			int ixe = min(ix_i+nke, nf_i);

			T f_mean = 0;
			for (auto ix = ix0; ix < ixe; ix++)
			{
				f_mean += f_i[ix];
			}

			f_o[ix_i] = f_mean/(ixe-ix0);
		};

		auto thr_mean = [&](const Range &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				krn_mean(ixy, f_i, fv);
			}
		};

		stream.set_n_act_stream(nf_i);
		stream.set_grid(1, nf_i);
		stream.exec(thr_mean);

		thrust::copy(fv.begin(), fv.end(), f_o.begin());
	}

	// get weighted position
	template<class TGrid, class TVector>
	void Rx_Ry_weight(TGrid &grid, TVector &Im, Value_type<TVector> x_i, Value_type<TVector> y_i, 
	Value_type<TVector> radius_i, Value_type<TVector> &x_o, Value_type<TVector> &y_o)
	{
		using T = Value_type<TVector>;

		T dR = grid.dRx;
		int nrl = static_cast<int>(floor(radius_i/dR+0.5));
		T R_max = dR*nrl;

		int ix0, ixe;
		get_bn(x_i, grid.nx, grid.dRx, R_max, grid.pbc_xy, ix0, ixe);
		ixe += ix0;

		int iy0, iye;
		get_bn(y_i, grid.ny, grid.dRy, R_max, grid.pbc_xy, iy0, iye);
		iye += iy0;

		T R2_max = pow(R_max, 2);

		x_o = 0;
		y_o = 0;
		T alpha = 0.5/pow(R_max, 2);
		T wx = 0;
		T wy = 0;
		for(auto ix = ix0; ix < ixe; ix++)
		{
			for(auto iy = iy0; iy < iye; iy++)
			{
				if(grid.R2(ix, iy, x_i, y_i) < R2_max)
				{
					int ixy = grid.ind_col(ix, iy);
					auto imv = Im[ixy];
					T Rx = grid.Rx(ix, x_i);
					T Ry = grid.Ry(iy, y_i);
					auto wxx = exp(-alpha*Rx*Rx)*imv;
					auto wyy = exp(-alpha*Ry*Ry)*imv;
					x_o += wxx*(Rx+x_i); 
					y_o += wyy*(Ry+y_i); 
					wx += wxx;
					wy += wyy;
				}
			}
		}
		x_o /= wx; 
		y_o /= wy; 
	}
	// get weighted position
	template<class TGrid, class TVector>	
	void Rx_Ry_fit(TGrid &grid, TVector &Im, Value_type<TVector> x_i, Value_type<TVector> y_i, 
	Value_type<TVector> radius_i, Value_type<TVector> &x_o, Value_type<TVector> &y_o)
	{
		using T = Value_type<TVector>;

		T dR = grid.dRx;
		int nrl = static_cast<int>(floor(radius_i/dR+0.5));
		T R_max = dR*nrl;

		int ix0, ixe;
		get_bn(x_i, grid.nx, grid.dRx, R_max, grid.pbc_xy, ix0, ixe);
		ixe += ix0;

		int iy0, iye;
		get_bn(y_i, grid.ny, grid.dRy, R_max, grid.pbc_xy, iy0, iye);
		iye += iy0;

		T R2_max = pow(R_max, 2);
		//get number of elements
		int m = 0;
		int n = 6;
		for(auto ix = ix0; ix < ixe; ix++)
		{
			for(auto iy = iy0; iy < iye; iy++)
			{
				m += (grid.R2(ix, iy, x_i, y_i) < R2_max)?1:0;
			}
		}

		Vector<T, e_host> M(m*n);
		Vector<T, e_host> b(m);
		int ic = 0;
		for(auto ix = ix0; ix < ixe; ix++)
		{
			for(auto iy = iy0; iy < iye; iy++)
			{
				if(grid.R2(ix, iy, x_i, y_i) < R2_max)
				{
					T Rx = grid.Rx(ix);
					T Ry = grid.Ry(iy);

					//set M values
					M[0*m+ic] = Rx*Rx;
					M[1*m+ic] = Ry*Ry;
					M[2*m+ic] = Rx;
					M[3*m+ic] = Ry;
					M[4*m+ic] = Rx*Ry;
					M[5*m+ic] = 1.0;

					//set b values
					int ixy = grid.ind_col(ix, iy);
					b[ic] = log(Im[ixy]+1e-04);

					ic++;
				}
			}
		}

		Vector<T, e_host> x(n);
		lapack::GELS<T> gels;
		gels(m, n, M.data(), 1, b.data(), x.data());
		x_o = -(2*x[1]*x[2]-x[3]*x[4])/(4*x[0]*x[1]-x[4]*x[4]);
		y_o = -(2*x[0]*x[3]-x[2]*x[4])/(4*x[0]*x[1]-x[4]*x[4]);
		if(fabs(x_o-x_i)>radius_i)
		{
			x_o = x_i;
		}
		if(fabs(y_o-y_i)>radius_i)
		{
			y_o = y_i;
		}
	}

	template<class TGrid, class TVector>
	void max_pos(Stream<e_host> &stream, TGrid &grid, TVector &Im, Value_type<TVector> x_i, 
	Value_type<TVector> y_i, Value_type<TVector> radius_i, Value_type<TVector> &x_o, Value_type<TVector> &y_o)
	{
		using T = Value_type<TVector>;
	
		auto get_rmin = [&](TVector &Im, T x_i, T y_i, T radius_i)->T
		{
			// get radial distribution
			Vector<T, e_host> rl;
			Vector<T, e_host> frl;
			radial_distribution_2D(grid, Im, x_i, y_i, radius_i, rl, frl);

			// log
			thrust::transform(frl.begin(), frl.end(), frl.begin(), [](const T &x){return log(x+1e-04); });

			// smooth
			int nkr = 1;
			smooth(stream, frl, nkr, frl);

			// derivative and minimun
			int idr = 0;
			T fr_min = frl[1]-frl[0];
			for(auto ir=0; ir<frl.size()-1; ir++)
			{
				auto fr = frl[ir+1]-frl[ir];
				if(fr_min>fr)
				{
					fr_min = fr;
					idr = ir;
				}
			}
			return ::fmax(2.0, rl[idr]);
		};

		//// weighted average position
		//auto rmin = get_rmin(Im, x_i, y_i, radius_i);
		//Rx_Ry_weight(grid, Im, x_i, y_i, rmin, x_o, y_o);

		//// fit position
		//rmin = get_rmin(Im, x_o, y_o, radius_i);
		//Rx_Ry_fit(grid, Im, x_o, y_o, rmin, x_o, y_o);

		// weighted average position
		auto rmin = get_rmin(Im, x_i, y_i, radius_i);
		Rx_Ry_weight(grid, Im, x_i, y_i, rmin, x_o, y_o);
	}

	// phase correlation function
	template <class TGrid, class TVector>
	void PCF(Stream<e_host> &stream, FFT2<Value_type<TGrid>, e_host> &fft2, TGrid &grid, 
		TVector &M_1i, TVector &M_2i, Value_type<TGrid> k, Value_type<TGrid> sigma, TVector &M_o)
	{
		using T = Value_type<TVector>;
		using TVector_c = Vector<complex<T>, e_host>;

		// get mean values
		T mean_M1 = 0;
		T mean_M2 = 0;
		auto thr_mean = [&](const Range &range)
		{
			auto mean_M1_partial = thrust::reduce(M_1i.begin()+range.ixy_0, M_1i.begin()+range.ixy_e);
			auto mean_M2_partial = thrust::reduce(M_2i.begin()+range.ixy_0, M_2i.begin()+range.ixy_e);

			stream.stream_mutex.lock();
			mean_M1 += mean_M1_partial;
			mean_M2 += mean_M2_partial;
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_mean);

		mean_M1 /= grid.nxy();
		mean_M2 /= grid.nxy();

		// apply Hanning filter and copy data to complex matrix
		TVector_c M_1(grid.nxy());
		TVector_c M_2(grid.nxy());

		T cx = c_2Pi/grid.nx;
		T cy = c_2Pi/grid.ny;

		auto krn_Hanning = [&](const int &ix, const int &iy, const TGrid &grid, TVector &M_1i, TVector &M_2i, TVector_c &M_1, TVector_c &M_2)
		{
			int ixy = grid.ind_col(ix, iy);
			T fxy = 0.25*(1.0-cos(cx*grid.Rx(ix)))*(1.0-cos(cy*grid.Ry(iy)));
			fxy = (fxy>k)?1.0:fxy/k;
			M_1[ixy] = complex<T>((M_1i[ixy]-mean_M1)*fxy);
			M_2[ixy] = complex<T>((M_2i[ixy]-mean_M2)*fxy);
		};

		auto thr_Hanning = [&](const Range &range)
		{
			host_detail::matrix_iter(range, krn_Hanning, grid, M_1i, M_2i, M_1, M_2);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_Hanning);

		// shift matrix
		fft2_shift(stream, grid, M_1);
		fft2_shift(stream, grid, M_2);

		fft2.forward(M_1);
		fft2.forward(M_2);

		T alpha = 0.5/(sigma*sigma);
		T g2d = sigma*sigma;
		auto krn_PCF = [&](const int &ix, const int &iy, const TGrid &grid, TVector_c &M_1, TVector_c &M_2)
		{
			int ixy = grid.ind_col(ix, iy);
			T gx = grid.igx_shift(ix);
			T gy = grid.igy_shift(iy);
			complex<T> z = M_1[ixy]*conj(M_2[ixy]);
			T mz = abs(z);
			M_2[ixy] = (isZero(mz))?0:z/((1.0+pow((gx*gx+gy*gy)/g2d, 8))*mz);
			//M_2[ixy] = (isZero(mz))?0:z*exp(-alpha*(gx*gx+gy*gy))/mz;
		};

		auto thr_PCF = [&](const Range &range)
		{
			host_detail::matrix_iter(range, krn_PCF, grid, M_1, M_2);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_PCF);

		fft2.inverse(M_2);

		// get minimun value and set real part
		T m_min = M_2[0].real();
		for(auto ixy = 0; ixy < M_o.size(); ixy++)
		{
			T x = M_o[ixy] = M_2[ixy].real();
			if(x<m_min)
			{
				m_min = x;
			}
		}

		// subtract minimum value
		auto thr_sub = [&](const Range &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				M_o[ixy] -= m_min;
			}
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_sub);

		fft2_shift(stream, grid, M_o);
	}

	// Gaussian convolution
	template <class TGrid, class TVector_c>
	void gaussian_convolution(Stream<e_host> &stream, FFT2<Value_type<TGrid>, e_host> &fft2, TGrid &grid, Value_type<TGrid> sigma, TVector_c &image)
	{	
		using T = Value_type<TGrid>;

		fft2.forward(image);

		T alpha = 2.0*c_Pi2*sigma*sigma;

		auto krn_gaussian = [&](const int &ix, const int &iy, const TGrid &grid, TVector_c &image)
		{
			int ixy = grid.ind_col(ix, iy);
			T g2 = grid.g2_shift(ix, iy);
			image[ixy] *= Value_type<TVector_c>(exp(-alpha*g2)/grid.nxy());
		};

		auto thr_gaussian_convolution = [&](const Range &range)
		{
			host_detail::matrix_iter(range, krn_gaussian, grid, image);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_gaussian_convolution);

		fft2.inverse(image);
	}

	// add Poisson noise
	template <class TVector>
	Value_type<TVector> add_poisson_noise(Value_type<TVector> SNR, TVector &image)
	{	
		using value_type = Value_type<TVector>;

		std::mt19937_64 gen;
		std::poisson_distribution<int> randp;

		auto get_std = [](TVector &image)->value_type
		{
			value_type x_mean = 0;
			value_type x_std = 0;
			for(auto ixy = 0; ixy < image.size(); ixy++)
			{
				auto x = image[ixy];
				x_mean += x;
				x_std += x*x;
			}
			x_mean /= image.size();
			return sqrt(x_std/image.size()-x_mean*x_mean);
		};

		auto get_SNR = [&](TVector &image, value_type image_std, value_type k)->value_type
		{
			value_type x_mean_n = 0;
			value_type x_std_n = 0;
			for(auto ixy = 0; ixy < image.size(); ixy++)
			{
				auto y = k*image[ixy];
				randp.param(std::poisson_distribution<int>::param_type(y));
				auto x = y-randp(gen);
				x_mean_n += x;
				x_std_n += x*x;
			}
			x_mean_n /= image.size();
			x_std_n = sqrt(x_std_n/image.size()-x_mean_n*x_mean_n);

			return k*image_std/x_std_n;
		};

		auto image_std = get_std(image);

		value_type SNR_k = get_SNR(image, image_std, 1);
	
		value_type k_0, k_e;
		k_0 = k_e = 1;

		if(SNR_k<SNR)
		{
			do
			{
				k_e *= 2;
				SNR_k = get_SNR(image, image_std, k_e);
			} while (SNR_k < SNR);
			k_0 = k_e/2;	
		}
		else
		{
			do
			{
				k_0 /= 2;
				SNR_k = get_SNR(image, image_std, k_0);
			} while (SNR_k >= SNR);
			k_e = 2*k_0;
		}


		// bisection method
		int c = 0;
		value_type k_m;
		do
		{
			k_m = 0.5*(k_0 + k_e);
			auto SNR_k = get_SNR(image, image_std, k_m);

			if(SNR_k < SNR)
				k_0 = k_m;
			else
				k_e = k_m;
			c++;
		} while ((fabs(SNR-SNR_k)>0.05) && (c<10));


		// add Poisson noise
		for(auto ixy = 0; ixy < image.size(); ixy++)
		{
			auto y = k_m*image[ixy];
			randp.param(std::poisson_distribution<int>::param_type(y));
			image[ixy] = randp(gen);
		}
		return k_m;
	}

	// get index to maximun distance
	int get_dmaxIndex_Point2Line(double x1, double y1, double x2, double y2, int ix1, int ix2, double *x, double *y)
	{
		int imax;
		double x0, y0, dmax;
		double d, c1, c2, c3;

		c1 = y2-y1;
		c2 = -(x2-x1); 
		c3 = x2*y1-y2*x1;
		d = sqrt(c1*c1*+c2*c2);
		c1 /= d; c2 /= d; c3 /= d;

		imax = 0; dmax = 0;
		for(auto i = ix1; i < ix2+1; i++)
		{
			x0 = x[i]; y0 = y[i];
			d = fabs(c1*x0+c2*y0+c3);
			if(d>dmax)
			{
				dmax = d;
				imax = i;
			}
		}

		return imax;
	}

	double getLengthCurve(int ix1, int ix2, double *x, double *y)
	{
		double x1, y1, x2, y2;
		double d, dx, dy;
		d = 0;
		for(auto i = ix1; i < ix2-1; i++)
		{
			x1 = x[i]; y1 = y[i];
			x2 = x[i+1]; y2 = y[i+1];
			dx = x2-x1; dy = y2-y1;
			d = d + sqrt(dx*dx + dy*dy);
		}
		return d;
	}

	double getLengthCurve(int ix1, int ix2, double *x, double *y, double lmax, int &il)
	{
		double x1, y1, x2, y2;
		double l, dx, dy;
		if(ix1<ix2)
		{
			l = 0; il = ix2;
			for(int i =ix1; i<ix2-1; i++)
			{
				x1 = x[i]; y1 = y[i];
				x2 = x[i+1]; y2 = y[i+1];
				dx = x2-x1; dy = y2-y1;
				l = l + sqrt(dx*dx + dy*dy);
				if((lmax>0)&&(l>=lmax))
				{
					il = i;
					break;
				}
			}
		}
		else
		{
			l = 0; il = ix2;
			for(int i =ix1; i>ix2; i--)
			{
				x1 = x[i-1]; y1 = y[i-1];
				x2 = x[i]; y2 = y[i];
				dx = x2-x1; dy = y2-y1;
				l = l + sqrt(dx*dx + dy*dy);
				if((lmax>0)&&(l>=lmax))
				{
					il = i;
					break;
				}
			}
		}
		return l;
	}

	// get information limit for regular gridBT
	void get_sigma_by_row(int ny, int nx, double *g, double *fg, double *sigma)
	{
		int nr = nx/2-1;
		vector<double> rl, frl, cfrl;
		double g_min = 0, g_max = nr, dg = fabs(g[1]-g[0]);

		rl.resize(nr);
		frl.resize(nr);
		cfrl.resize(nr);

		for(int i = 0; i < nr; i++)
		{
			rl[i] = 0.5+i;
		}

		// Shift and Normalize
		auto shift_normalize = [&](const double &x, const double &x_min, const double &x_max)->double
		{
			return (x-x_min)/(x_max-x_min);
		};

		for(int iy =0; iy < ny; iy++)
		{
			std::fill(frl.begin(), frl.end(), 0.0);
			std::fill(cfrl.begin(), cfrl.end(), 0.0);
			for(int ix =0; ix < nx; ix++)
			{
				if((g_min <= g[ix])&&(g[ix]<g_max))
				{
					auto j = (int)floor((g[ix]-g_min)/dg);
					frl[j] += fg[ix];
					cfrl[j] += 1.0;
				}
			}

			for(auto i = 0; i < nr; i++)
			{
				if(cfrl[i]>0)
					frl[i] /= cfrl[i];
			}

			frl[0] = 1.01*frl[1];

			auto rl_min = rl.front();
			auto rl_max = rl.back();
			auto frl_min = *std::min_element(frl.begin(), frl.end());
			auto frl_max = *std::max_element(frl.begin(), frl.end());

			int k0 = 0;
			auto x = shift_normalize(rl[k0], rl_min, rl_max);
			auto y = shift_normalize(frl[k0], frl_min, frl_max);
			auto d2_min = x*x + y*y;
			for(int i =1; i < nr; i++)
			{
				auto x = shift_normalize(rl[i], rl_min, rl_max);
				auto y = shift_normalize(frl[i], frl_min, frl_max);
				auto d2 = x*x + y*y;
				if(d2 < d2_min)
				{
					d2_min = d2;
					k0 = i;
				}
			}

			int ke = static_cast<int>(std::max_element(frl.begin() + k0, frl.end()) - frl.begin());

			sigma[iy] = 0.5*(rl[ke]+rl[k0]);
		}
	}

	// get information limit for regular gridBT
	double sigma(int ny, int nx, double *g, double *fg)
	{
		auto dgx = fabs(g[1]-g[0]);
		auto dgy = fabs(g[ny]-g[0]);
		int nr = (nx*dgx < ny*dgy)?(nx/2-1):(ny/2-1);
		vector<double> gl, rl, frl, cfrl;

		gl.resize(nr+1);	
		rl.resize(nr);
		frl.resize(nr);
		cfrl.resize(nr);

		std::iota(gl.begin(), gl.end(), 0);

		// radial distribution
		auto nxy = nx*ny;
		radial_distribution_2D(nxy, g, fg, nr+1, gl.data(), rl.data(), frl.data(), cfrl.data(), true, 0);
		frl[0] = 1.01*frl[1];
		auto rl_min = rl.front();
		auto rl_max = rl.back();
		auto frl_min = *std::min_element(frl.begin(), frl.end());
		auto frl_max = *std::max_element(frl.begin(), frl.end());

		// Shift and Normalize
		auto shift_normalize = [&](const double &x, const double &x_min, const double &x_max)->double
		{
			return (x-x_min)/(x_max-x_min);
		};

		int k0 = 0;
		auto x = shift_normalize(rl[k0], rl_min, rl_max);
		auto y = shift_normalize(frl[k0], frl_min, frl_max);
		auto d2_min = x*x + y*y;
		for(int i =1; i < nr; i++)
		{
			auto x = shift_normalize(rl[i], rl_min, rl_max);
			auto y = shift_normalize(frl[i], frl_min, frl_max);
			auto d2 = x*x + y*y;
			if(d2 < d2_min)
			{
				d2_min = d2;
				k0 = i;
			}
		}

		int ke = static_cast<int>(std::max_element(frl.begin() + k0, frl.end()) - frl.begin());

		return rl[ke];
	}

	// get information limit for regular gridBT
	double FFT_information_limit_2D(int ny, int nx, bool shift, double *fI)
	{
		int nr = min(nx/2, ny/2)-1;
		double *r, *fIr;

		r = new double[nr];
		fIr = new double[nr]; 

		// Cumulative radial integration
		getCumRadDist_2D(ny, nx, shift, fI, nr, r, fIr);

		// Shift and Normalize
		double r0 = r[0], fIr0 = fIr[0];
		for(int i = 0; i < nr; i++)
		{
			r[i] = (r[i]-r0)/(r[nr-1]-r0);
			fIr[i] = (fIr[i]-fIr0)/(fIr[nr-1]-fIr0);
		}

		int ir1, ir2, irm;
		double x1, y1, x2, y2;

		ir1 = 0; ir2 = nr-1;
		x1 = r[ir1]; y1 = fIr[ir1];
		x2 = r[ir2]; y2 = fIr[ir2];
	
		irm = get_dmaxIndex_Point2Line(x1, y1, x2, y2, ir1, ir2, r, fIr);
		double fIr_lim = 0.45*fIr[irm];

		for(int i = 0; i < nr; i++)
			if(fIr[i]>fIr_lim)
			{
				irm = i-1;
				break;
			}

		delete [] r;
		delete [] fIr;

		return irm;
	}

	template<class InputIterator, class TVector>
	void match_vectors(InputIterator A_base_first, InputIterator A_base_last, TVector &A_search)
	{
		using T = Value_type<TVector>;

		std::size_t A_base_size = static_cast<int>(A_base_last - A_base_first);
		Vector<bool, e_host> A_base_c(A_base_size, true);
		int m_size = min(A_base_size, A_search.size());
		TVector A_d;
		A_d.reserve(m_size);

		for(auto i = 0; i<A_search.size(); i++)
		{
			T val = A_search[i];
			auto it = std::min_element(A_base_first, A_base_last, [&val](const T &a, const T &b)->bool{ return fabs(a-val+Epsilon<float>::rel)<fabs(b-val); });
			auto imin = static_cast<int>(it - A_base_first);

			if(A_base_c[imin])
			{
				A_base_c[imin] = false;
				A_d.push_back(*it);
			}
		}
		A_search = A_d;
		A_search.shrink_to_fit();
	}

	template<class InputIterator, class TVector>
	void match_vectors(InputIterator A_base_first, InputIterator A_base_last, TVector &A_search, Vector<int, e_host> &A_idx)
	{
		using T = Value_type<TVector>;

		std::size_t A_base_size = static_cast<int>(A_base_last - A_base_first);
		Vector<bool, e_host> A_base_c(A_base_size, true);
		int m_size = min(A_base_size, A_search.size());
		A_idx.clear();
		A_idx.reserve(m_size);
		TVector A_d;
		A_d.reserve(m_size);

		for(auto i = 0; i<A_search.size(); i++)
		{
			T val = A_search[i];
			auto it = std::min_element(A_base_first, A_base_last, [&val](const T &a, const T &b)->bool{ return fabs(a-val+Epsilon<float>::rel)<fabs(b-val); });
			auto imin = static_cast<int>(it - A_base_first);

			if(A_base_c[imin])
			{
				A_base_c[imin] = false;
				A_idx.push_back(imin);
				A_d.push_back(*it);
			}
		}
		A_search = A_d;
		A_search.shrink_to_fit();
		A_idx.shrink_to_fit();
	}

	int get_close_Prime(int64_t n)
	{
		std::vector<int64_t> p2(18), p3(5), p5(3), p7(2);

		for(auto i = 0; i<p2.size(); i++)
		{
			p2[i] = std::pow<int64_t>(2, i+1);

			if(i < p3.size())
			{
				p3[i] = std::pow<int64_t>(3, i);
			}
			if(i < p5.size())
			{
				p5[i] = std::pow<int64_t>(5, i);
			}
			if(i < p7.size())
			{
				p7[i] = std::pow<int64_t>(7, i);
			}
		}

		int64_t prime, prime_0 = 128, prime_e = std::pow<int64_t>(2, 18);
		int64_t n_c = 64, dn = prime_e, idn;

		for(auto ip7:p7)
		{
			for(auto ip5:p5)
			{
				for(auto ip3:p3)
				{
					for(auto ip2:p2)
					{
						prime = ip7*ip5*ip3*ip2;
						idn = std::abs(n-prime);
						if((prime_0<=prime)&&(prime<=prime_e)&&(idn<dn))
						{
							dn = idn;
							n_c = prime;
						}
					}
				}
			}
		}
		return static_cast<int>(n_c);
	}

	void get_prime_sampling(int &nx, int &ny)
	{
		nx = get_close_Prime(nx);
		ny = get_close_Prime(ny);
	}

	/***************************************************************************/
	/***************************************************************************/

	template<class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	assign(Stream<e_host> &stream, TVector_1 &M_i, TVector_2 &M_o, Vector<Value_type<TVector_2>, e_host> *M_i_h =nullptr)
	{
		M_o.resize(M_i.size());
		auto thr_assign = [&](const Range &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				M_o[ixy] = M_i[ixy];
			}
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_assign);
	}

	template<class TVector_1, class TVector_2>
	typename std::enable_if<is_host_vector_and_device_vector<TVector_1, TVector_2>::value 
	&& is_complex<Value_type<TVector_2>>::value && !std::is_same<Value_type<TVector_1>, Value_type<TVector_2>>::value, void>::type
	assign(Stream<e_host> &stream, TVector_1 &M_i, TVector_2 &M_o, Vector<Value_type<TVector_2>, e_host> *M_i_h =nullptr)
	{
		Vector<Value_type<TVector_2>, e_host> M_h;
		M_i_h = (M_i_h == nullptr)?&M_h:M_i_h;

		// copy data to the same output type
		assign(stream, M_i, *M_i_h);

		// data transfer from CPU to GPU
		M_o.assign(M_i_h->begin(), M_i_h->end());
	}

	template<class TVector_1, class TVector_2>
	typename std::enable_if<is_host_vector_and_device_vector<TVector_1, TVector_2>::value 
	&& (!is_complex<Value_type<TVector_2>>::value || std::is_same<Value_type<TVector_1>, Value_type<TVector_2>>::value), void>::type
	assign(Stream<e_host> &stream, TVector_1 &M_i, TVector_2 &M_o, Vector<Value_type<TVector_2>, e_host> *M_i_h =nullptr)
	{
		M_o.assign(M_i.begin(), M_i.end());
	}

	template<class TVector>
	enable_if_host_vector<TVector, void>
	fill(Stream<e_host> &stream, TVector &M_io, Value_type<TVector> value_i)
	{
		auto thr_fill = [&](const Range &range)
		{
			thrust::fill(M_io.begin()+range.ixy_0, M_io.begin()+range.ixy_e, value_i);
		};

		stream.set_n_act_stream(M_io.size());
		stream.set_grid(1, M_io.size());
		stream.exec(thr_fill);
	}

	template<class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	scale(Stream<e_host> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_scale = [&](const Range &range)
		{
			thrust::transform(M_i.begin()+range.ixy_0, M_i.begin()+range.ixy_e, 
				M_o.begin()+range.ixy_0, functor::scale<value_type>(w_i));
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_scale);
	}

	template<class TVector>
	enable_if_host_vector<TVector, void>
	scale(Stream<e_host> &stream, Value_type<TVector> w_i, TVector &M_io)
	{
		scale(stream, w_i, M_io, M_io);
	}

	template<class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	square(Stream<e_host> &stream, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_square = [&](const Range &range)
		{
			thrust::transform(M_i.begin()+range.ixy_0, M_i.begin()+range.ixy_e, 
				M_o.begin()+range.ixy_0, functor::square<value_type>());
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_square);
	}

	template<class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	square_scale(Stream<e_host> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_square_scale = [&](const Range &range)
		{
			thrust::transform(M_i.begin()+range.ixy_0, M_i.begin()+range.ixy_e, 
				M_o.begin()+range.ixy_0, functor::square_scale<value_type>(w_i));
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_square_scale);
	}

	template<class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	add(Stream<e_host> &stream, TVector_1 &M1_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_add = [&](const Range &range)
		{
			thrust::transform(M1_i.begin()+range.ixy_0, M1_i.begin()+range.ixy_e, 
				M2_i.begin()+range.ixy_0, M_o.begin()+range.ixy_0, functor::add<value_type>());
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_add);
	}

	template<class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	add(Stream<e_host> &stream, TVector_1 &M_i, TVector_2 &M_io)
	{
		add(stream, M_i, M_io, M_io);
	}

	template<class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	add_scale(Stream<e_host> &stream, Value_type<TVector_2> w1_i, TVector_1 &M1_i, 
	Value_type<TVector_2> w2_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_add_scale = [&](const Range &range)
		{
			thrust::transform(M1_i.begin()+range.ixy_0, M1_i.begin()+range.ixy_e, 
				M2_i.begin()+range.ixy_0, M_o.begin()+range.ixy_0, functor::add_scale_i<value_type>(w1_i, w2_i));
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_add_scale);
	}

	template<class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	add_scale(Stream<e_host> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_io)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_add_scale = [&](const Range &range)
		{
			thrust::transform(M_i.begin()+range.ixy_0, M_i.begin()+range.ixy_e, 
				M_io.begin()+range.ixy_0, M_io.begin()+range.ixy_0, functor::add_scale<value_type>(w_i));
		};

		stream.set_n_act_stream(M_io.size());
		stream.set_grid(1, M_io.size());
		stream.exec(thr_add_scale);
	}

	template<class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	add_square(Stream<e_host> &stream, TVector_1 &M1_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_add_square = [&](const Range &range)
		{
			thrust::transform(M1_i.begin()+range.ixy_0, M1_i.begin()+range.ixy_e, 
				M2_i.begin()+range.ixy_0, M_o.begin()+range.ixy_0, functor::add_square_i<value_type>());
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_add_square);
	}

	template<class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	add_square(Stream<e_host> &stream, TVector_1 &M_i, TVector_2 &M_io)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_add_square = [&](const Range &range)
		{
			thrust::transform(M_i.begin()+range.ixy_0, M_i.begin()+range.ixy_e, 
				M_io.begin()+range.ixy_0, M_io.begin()+range.ixy_0, functor::add_square<value_type>());
		};

		stream.set_n_act_stream(M_io.size());
		stream.set_grid(1, M_io.size());
		stream.exec(thr_add_square);
	}

	template<class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	add_square_scale(Stream<e_host> &stream, Value_type<TVector_2> w1_i, TVector_1 &M1_i, Value_type<TVector_2> w2_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_add_square_scale = [&](const Range &range)
		{
			thrust::transform(M1_i.begin()+range.ixy_0, M1_i.begin()+range.ixy_e, 
				M2_i.begin()+range.ixy_0, M_o.begin()+range.ixy_0, functor::add_square_scale_i<value_type>(w1_i, w2_i));
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_add_square_scale);
	}

	template<class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	add_square_scale(Stream<e_host> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_io)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_add_square_scale = [&](const Range &range)
		{
			thrust::transform(M_i.begin()+range.ixy_0, M_i.begin()+range.ixy_e, 
				M_io.begin()+range.ixy_0, M_io.begin()+range.ixy_0, functor::add_square_scale<value_type>(w_i));
		};

		stream.set_n_act_stream(M_io.size());
		stream.set_grid(1, M_io.size());
		stream.exec(thr_add_square_scale);
	}

	template<class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	multiply(Stream<e_host> &stream, TVector_1 &M1_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_multiply = [&](const Range &range)
		{
			thrust::transform(M1_i.begin()+range.ixy_0, M1_i.begin()+range.ixy_e, 
				M2_i.begin()+range.ixy_0, M_o.begin()+range.ixy_0, functor::multiply<value_type>());
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_multiply);
	}

	template<class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	multiply(Stream<e_host> &stream, TVector_1 &M_i, TVector_2 &M_io)
	{
		multiply(stream, M_i, M_io, M_io);
	}

	template<class TVector>
	enable_if_host_vector<TVector, Value_type<TVector>>
	sum(Stream<e_host> &stream, TVector &M_i)
	{
		using value_type = Value_type<TVector>;

		value_type sum_total = 0;
		auto thr_sum = [&](const Range &range)
		{
			auto sum_partial = thrust::reduce(M_i.begin()+range.ixy_0, M_i.begin()+range.ixy_e);

			stream.stream_mutex.lock();
			sum_total += sum_partial;
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(M_i.size());
		stream.set_grid(1, M_i.size());
		stream.exec(thr_sum);

		return sum_total;
	}

	template<class TVector>
	enable_if_host_vector<TVector, Value_type_r<TVector>>
	sum_square(Stream<e_host> &stream, TVector &M_i)
	{
		using value_type_r = Value_type_r<TVector>;

		value_type_r sum_total = 0;
		auto thr_sum_square = [&](const Range &range)
		{
			auto sum_partial = thrust::transform_reduce(M_i.begin()+range.ixy_0, M_i.begin()+range.ixy_e, 
			functor::square<value_type_r>(), value_type_r(0), functor::add<value_type_r>());

			stream.stream_mutex.lock();
			sum_total += sum_partial;
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(M_i.size());
		stream.set_grid(1, M_i.size());
		stream.exec(thr_sum_square);

		return sum_total;
	}

	template<class TVector>
	enable_if_host_vector<TVector, Value_type_r<TVector>>
	mean(Stream<e_host> &stream, TVector &M_i)
	{
		return sum(stream, M_i)/M_i.size();
	}

	template<class TVector>
	enable_if_host_vector<TVector, void>
	mean_std(Stream<e_host> &stream, TVector &M_i, Value_type_r<TVector> &x_mean, Value_type_r<TVector> &x_std)
	{
		using value_type_r = Value_type_r<TVector>;

		x_mean = mean(stream, M_i);

		x_std = 0;
		auto thr_std = [&](const Range &range)
		{
			auto x_std_partial = thrust::transform_reduce(M_i.begin()+range.ixy_0, M_i.begin()+range.ixy_e, 
			functor::square_dif<value_type_r>(x_mean), value_type_r(0), functor::add<value_type_r>());

			stream.stream_mutex.lock();
			x_std += x_std_partial;
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(M_i.size());
		stream.set_grid(1, M_i.size());
		stream.exec(thr_std);

		x_std = sqrt(x_std/M_i.size());
	}

	template<class TVector>
	enable_if_host_vector<TVector, Value_type_r<TVector>>
	std(Stream<e_host> &stream, TVector &M_i)
	{
		using value_type_r = Value_type_r<TVector>;

		value_type_r x_mean , x_std;
		mean_std(stream, M_i, x_mean, x_std);
		return x_std;
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	phase_factor_1D(Stream<e_host> &stream, TGrid &grid, Value_type<TGrid> x, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto thr_phase_factor_1D = [&](const Range &range)
		{
			thrust::counting_iterator<int> first(range.ixy_0);
			thrust::counting_iterator<int> last = first + range.ixy_e;
			thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(first, fPsi_i.begin()+range.ixy_0, fPsi_o.begin()+range.ixy_0)), 
							 thrust::make_zip_iterator(thrust::make_tuple(last, fPsi_i.begin()+range.ixy_e, fPsi_o.begin()+range.ixy_e)), 
							 functor::phase_factor<TGrid>(grid, x));
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(1, grid.nx);
		stream.exec(thr_phase_factor_1D);
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	phase_factor_2D(Stream<e_host> &stream, TGrid &grid, Value_type<TGrid> x, Value_type<TGrid> y, 
	TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto thr_phase_factor = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::phase_factor_2D<TGrid, TVector_c>, grid, x, y, fPsi_i, fPsi_o);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_phase_factor);
	}

	// shift 2d
	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	shift_2D(Stream<e_host> &stream, FFT2<Value_type<TGrid>, e_host> &fft2, TGrid &grid, 
	Value_type<TGrid> x, Value_type<TGrid> y, TVector_c &Im)
	{	
		fft2_shift(stream, grid, Im);
		fft2.forward(Im);
		phase_factor_2D(stream, grid, -c_2Pi*x, -c_2Pi*y, Im, Im);
		fft2.inverse(Im);
		fft2_shift(stream, grid, Im);
	}

	template<class TGrid, class TVector_r>
	enable_if_host_vector<TVector_r, Value_type<TVector_r>>
	atom_cost_function(TGrid &grid, const Atom_Ip<Value_type<TGrid>> &atom_Ip, TVector_r &M_i)
	{
		return host_detail::atom_cost_function<typename TVector_r::value_type>(grid, atom_Ip, M_i);
	}

	template<class TGrid, class TVector_r>
	enable_if_host_vector<TVector_r, void>
	subtract_atom(Stream<e_host> &stream, TGrid &grid, Vector<Atom_Ip<Value_type<TGrid>>, e_host> &atom_Ip, TVector_r &M_i)
	{
		if(stream.n_act_stream<=0)
		{
			return;
		}

		for(auto istream = 0; istream < stream.n_act_stream; istream++)
		{
			stream[istream] = std::thread(std::bind(host_detail::subtract_atom<typename TVector_r::value_type>, std::ref(stream), std::ref(grid), std::ref(atom_Ip[istream]), std::ref(M_i)));
		}
		stream.synchronize();
	}

	template<class TQ1>
	enable_if_host<TQ1, void>
	get_cubic_poly_coef_Vz(Stream<e_host> &stream, ePotential_Type potential_type, 
	TQ1 &qz, Vector<Atom_Vp<Value_type<TQ1>>, e_host> &atom_Vp)
	{
		using value_type_r = Value_type<TQ1>;

		if(stream.n_act_stream<=0)
		{
			return;
		}

		auto thr_cubic_poly_coef = [](const ePotential_Type &potential_type, TQ1 &qz, Atom_Vp<value_type_r> &atom_Vp)
		{
			if(atom_Vp.charge==0)
			{
				switch(potential_type)
				{
					case ePT_Doyle_0_4:
						host_detail::linear_Vz<ePT_Doyle_0_4, 0, typename TQ1::value_type>(qz, atom_Vp);
						break;
					case ePT_Peng_0_4:
						host_detail::linear_Vz<ePT_Peng_0_4, 0, typename TQ1::value_type>(qz, atom_Vp);
						break;
					case ePT_Peng_0_12:
						host_detail::linear_Vz<ePT_Peng_0_12, 0, typename TQ1::value_type>(qz, atom_Vp);
						break;
					case ePT_Kirkland_0_12:
						host_detail::linear_Vz<ePT_Kirkland_0_12, 0, typename TQ1::value_type>(qz, atom_Vp);
						break;
					case ePT_Weickenmeier_0_12:
						host_detail::linear_Vz<ePT_Weickenmeier_0_12, 0, typename TQ1::value_type>(qz, atom_Vp);
						break;
					case ePT_Lobato_0_12:
						host_detail::linear_Vz<ePT_Lobato_0_12, 0, typename TQ1::value_type>(qz, atom_Vp);
						break;
				}
			}
			else
			{
				switch(potential_type)
				{
					case ePT_Doyle_0_4:
						host_detail::linear_Vz<ePT_Doyle_0_4, 1, typename TQ1::value_type>(qz, atom_Vp);
						break;
					case ePT_Peng_0_4:
						host_detail::linear_Vz<ePT_Peng_0_4, 1, typename TQ1::value_type>(qz, atom_Vp);
						break;
					case ePT_Peng_0_12:
						host_detail::linear_Vz<ePT_Peng_0_12, 1, typename TQ1::value_type>(qz, atom_Vp);
						break;
					case ePT_Kirkland_0_12:
						host_detail::linear_Vz<ePT_Kirkland_0_12, 1, typename TQ1::value_type>(qz, atom_Vp);
						break;
					case ePT_Weickenmeier_0_12:
						host_detail::linear_Vz<ePT_Weickenmeier_0_12, 1, typename TQ1::value_type>(qz, atom_Vp);
						break;
					case ePT_Lobato_0_12:
						host_detail::linear_Vz<ePT_Lobato_0_12, 1, typename TQ1::value_type>(qz, atom_Vp);
						break;
				}
			}
			host_detail::cubic_poly_coef<typename TQ1::value_type>(atom_Vp); 
		};

		for(auto istream = 0; istream < stream.n_act_stream; istream++)
		{
			stream[istream] = std::thread(std::bind(thr_cubic_poly_coef, potential_type, std::ref(qz), std::ref(atom_Vp[istream])));
		}
		stream.synchronize();
	}

	template<class TGrid, class TVector_r>
	enable_if_host_vector<TVector_r, void>
	eval_cubic_poly(Stream<e_host> &stream, TGrid &grid, 
	Vector<Atom_Vp<Value_type<TGrid>>, e_host> &atom_Vp, TVector_r &V_0)
	{
		if(stream.n_act_stream<=0)
		{
			return;
		}

		for(auto istream = 0; istream < stream.n_act_stream; istream++)
		{
			stream[istream] = std::thread(std::bind(host_detail::eval_cubic_poly<typename TVector_r::value_type>, std::ref(stream), std::ref(grid), std::ref(atom_Vp[istream]), std::ref(V_0)));
		}
		stream.synchronize();
	}

	template<class TGrid, class TVector>
	enable_if_host_vector<TVector, void>
	fft2_shift(Stream<e_host> &stream, TGrid &grid, TVector &M_io)
	{
		auto thr_fft2_shift = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::fft2_shift<TGrid, TVector>, grid, M_io);
		};

		stream.set_n_act_stream(grid.nxh);
		stream.set_grid(grid.nxh, grid.nyh);
		stream.exec(thr_fft2_shift);
	}

	template<class TGrid, class TVector>
	enable_if_host_vector<TVector, Value_type<TVector>>
	sum_over_Det(Stream<e_host> &stream, TGrid &grid, Value_type<TGrid> g_min, Value_type<TGrid> g_max, TVector &M_i)
	{
		using value_type_r = Value_type<TGrid>;
		using value_type = Value_type<TVector>;

		value_type_r g2_min = pow(g_min, 2);
		value_type_r g2_max = pow(g_max, 2);
		value_type sum = 0;

		auto thr_sum_over_Det = [&](const Range &range)
		{
			value_type sum_partial = 0;
			host_detail::matrix_iter(range, host_device_detail::sum_over_Det<TGrid, TVector>, grid, g2_min, g2_max, M_i, sum_partial);

			stream.stream_mutex.lock();
			sum += sum_partial;
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_sum_over_Det);

		return sum;
	}

	template<class TGrid, class TVector>
	enable_if_host_vector<TVector, Value_type<TGrid>>
	sum_square_over_Det(Stream<e_host> &stream, TGrid &grid, Value_type<TGrid> g_min, Value_type<TGrid> g_max, TVector &M_i)
	{
		using value_type_r = Value_type<TGrid>;

		value_type_r g2_min = pow(g_min, 2);
		value_type_r g2_max = pow(g_max, 2);
		value_type_r sum = 0;

		auto thr_sum_square_over_Det = [&](const Range &range)
		{
			value_type_r sum_partial = 0;
			host_detail::matrix_iter(range, host_device_detail::sum_square_over_Det<TGrid, TVector>, grid, g2_min, g2_max, M_i, sum_partial);

			stream.stream_mutex.lock();
			sum += sum_partial;
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_sum_square_over_Det);

		return sum;
	}

	template<class TGrid, class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, Value_type<TGrid>>
	sum_square_over_Det(Stream<e_host> &stream, TGrid &grid, TVector_1 &S_i, TVector_2 &M_i)
	{
		using value_type_r = Value_type<TGrid>;

		value_type_r sum = 0;

		auto thr_sum_square_over_Det = [&](const Range &range)
		{
			value_type_r sum_partial = 0;
			host_detail::matrix_iter(range, host_device_detail::sum_square_over_Det<TGrid, TVector_1, TVector_2>, grid, S_i, M_i, sum_partial);

			stream.stream_mutex.lock();
			sum += sum_partial;
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_sum_square_over_Det);

		return sum;
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	bandwidth_limit(Stream<e_host> &stream, TGrid &grid, Value_type<TGrid> g_min, Value_type<TGrid> g_max, Value_type<TVector_c> w, TVector_c &M_io)
	{
		using value_type_r = Value_type<TGrid>;

		value_type_r g2_min = pow(g_min, 2);
		value_type_r g2_max = pow(g_max, 2);

		auto thr_bandwidth_limit = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::bandwidth_limit<TGrid, TVector_c>, grid, g2_min, g2_max, w, M_io);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_bandwidth_limit);
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	phase_components(Stream<e_host> &stream, TGrid &grid, Value_type<TGrid> gxu, Value_type<TGrid> gyu, TVector_c &V_x_o, TVector_c &V_y_o)
	{
		auto thr_phase_components = [&](const Range &range)
		{
			thrust::counting_iterator<int> first(range.ixy_0);
			thrust::counting_iterator<int> last = first + range.ixy_e;

			thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(first, V_x_o.begin()+range.ixy_0, V_y_o.begin()+range.ixy_0)), 
							 thrust::make_zip_iterator(thrust::make_tuple(last, V_x_o.begin()+range.ixy_e, V_y_o.begin()+range.ixy_e)), 
							 functor::phase_components<TGrid>(grid, c_2Pi, gxu, gyu));
		};

		stream.set_n_act_stream(grid.nx_ny_max());
		stream.set_grid(1, grid.nx_ny_max());
		stream.exec(thr_phase_components);
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	phase_multiplication(Stream<e_host> &stream, TGrid &grid, TVector_c &exp_x_i, TVector_c &exp_y_i, TVector_c &psi_i, TVector_c &psi_o)
	{
		auto thr_phase_multiplication = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::phase_multiplication<TGrid, TVector_c>, grid, exp_x_i, exp_y_i, psi_i, psi_o);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_phase_multiplication);
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	propagator_components(Stream<e_host> &stream, TGrid &grid, Value_type<TGrid> gxu, Value_type<TGrid> gyu, Value_type<TGrid> w, TVector_c &V_x_o, TVector_c &V_y_o)
	{
		auto thr_propagator_components = [&](const Range &range)
		{
			thrust::counting_iterator<int> first(range.ixy_0);
			thrust::counting_iterator<int> last = first + range.ixy_e;

			thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(first, V_x_o.begin()+range.ixy_0, V_y_o.begin()+range.ixy_0)), 
							 thrust::make_zip_iterator(thrust::make_tuple(last, V_x_o.begin()+range.ixy_e, V_y_o.begin()+range.ixy_e)), 
							 functor::propagator_components<TGrid>(grid, w, gxu, gyu));
		};

		stream.set_n_act_stream(grid.nx_ny_max());
		stream.set_grid(1, grid.nx_ny_max());
		stream.exec(thr_propagator_components);
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	propagator_multiplication(Stream<e_host> &stream, TGrid &grid, TVector_c &prop_x_i, TVector_c &prop_y_i, TVector_c &psi_i, TVector_c &psi_o)
	{
		auto thr_propagator_multiplication = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::propagator_multiplication<TGrid, TVector_c>, grid, prop_x_i, prop_y_i, psi_i, psi_o);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_propagator_multiplication);
	}

	template<class TGrid, class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	transmission_function(Stream<e_host> &stream, TGrid &grid, eElec_Spec_Int_Model elec_spec_int_model, Value_type<TGrid> w, TVector_1 &V0_i, TVector_2 &Trans_o)
	{	
		using value_type_r = Value_type<TGrid>;

		auto thr_transmission_funtion = [&](const Range &range)
		{
			thrust::transform(V0_i.begin()+range.ixy_0, V0_i.begin()+range.ixy_e, 
				Trans_o.begin()+range.ixy_0, functor::transmission_function<value_type_r>(w, elec_spec_int_model));
		};

		stream.set_n_act_stream(grid.nxy());
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_transmission_funtion);
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	probe(Stream<e_host> &stream, TGrid &grid, Lens<Value_type<TGrid>> &lens, Value_type<TGrid> x, Value_type<TGrid> y, TVector_c &fPsi_o)
	{
		auto thr_probe = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::probe<TGrid, TVector_c>, grid, lens, x, y, fPsi_o);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_probe);

		auto total = sum_square(stream, fPsi_o);
		multem::scale(stream, sqrt(1.0/total), fPsi_o);
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	apply_CTF(Stream<e_host> &stream, TGrid &grid, Lens<Value_type<TGrid>> &lens, Value_type<TGrid> gxu, Value_type<TGrid> gyu, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto thr_apply_CTF = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::apply_CTF<TGrid, TVector_c>, grid, lens, gxu, gyu, fPsi_i, fPsi_o);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_apply_CTF);
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	apply_PCTF(Stream<e_host> &stream, TGrid &grid, Lens<Value_type<TGrid>> &lens, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto thr_apply_PCTF = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::apply_PCTF<TGrid, TVector_c>, grid, lens, fPsi_i, fPsi_o);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_apply_PCTF);
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	kernel_xyz(Stream<e_host> &stream, TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_host> &fft2, TVector_c &k_x, TVector_c &k_y, TVector_c &k_z)
	{
		eels.factor = host_detail::Lorentz_factor(stream, grid, eels);

		auto thr_kernel_xyz = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::kernel_xyz<TGrid, TVector_c>, grid, eels, k_x, k_y, k_z);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_kernel_xyz);

		fft2.inverse(k_x);
		fft2.inverse(k_y);
		fft2.inverse(k_z);
	}
	
	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	kernel_x(Stream<e_host> &stream, TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_host> &fft2, TVector_c &k_x)
	{
		eels.factor = host_detail::Lorentz_factor(stream, grid, eels);

		auto thr_kernel_x = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::kernel_x<TGrid, TVector_c>, grid, eels, k_x);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_kernel_x);

		fft2.inverse(k_x);
	}
	
	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	kernel_y(Stream<e_host> &stream, TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_host> &fft2, TVector_c &k_y)
	{
		eels.factor = host_detail::Lorentz_factor(stream, grid, eels);

		auto thr_kernel_y = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::kernel_y<TGrid, TVector_c>, grid, eels, k_y);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_kernel_y);

		fft2.inverse(k_y);
	}
	
	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	kernel_z(Stream<e_host> &stream, TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_host> &fft2, TVector_c &k_z)
	{
		eels.factor = host_detail::Lorentz_factor(stream, grid, eels);

		auto thr_kernel_z = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::kernel_z<TGrid, TVector_c>, grid, eels, k_z);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_kernel_z);

		fft2.inverse(k_z);
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	kernel_mn1(Stream<e_host> &stream, TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_host> &fft2, TVector_c &k_mn1)
	{
		eels.factor = host_detail::Lorentz_factor(stream, grid, eels);

		auto thr_kernel_mn1 = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::kernel_mn1<TGrid, TVector_c>, grid, eels, k_mn1);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_kernel_mn1);

		fft2.inverse(k_mn1);
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	kernel_mp1(Stream<e_host> &stream, TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_host> &fft2, TVector_c &k_mp1)
	{
		eels.factor = host_detail::Lorentz_factor(stream, grid, eels);

		auto thr_kernel_mp1 = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::kernel_mp1<TGrid, TVector_c>, grid, eels, k_mp1);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_kernel_mp1);

		fft2.inverse(k_mp1);
	}

	/***************************************************************************/
	/***************************************************************************/
	template<class TVector_i, class TVector_o>
	enable_if_host_vector_and_host_vector<TVector_i, TVector_o, void>
	copy_to_host(Stream<e_host> &stream, TVector_i &M_i, TVector_o &M_o, 
	Vector<Value_type<TVector_i>, e_host> *M_i_h =nullptr)
	{
		multem::assign(stream, M_i, M_o, M_i_h);
	}

	template<class TVector_i, class TVector_o>
	enable_if_host_vector_and_host_vector<TVector_i, TVector_o, void>
	add_scale_to_host(Stream<e_host> &stream, Value_type<TVector_i> w_i, 
	TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_host> *M_i_h =nullptr)
	{
		multem::add_scale(stream, w_i, M_i, M_o);
	}

	template<class TVector_i, class TVector_o>
	enable_if_host_vector_and_host_vector<TVector_i, TVector_o, void>
	add_square_scale_to_host(Stream<e_host> &stream, Value_type<TVector_o> w_i, 
	TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_host> *M_i_h =nullptr)
	{
		multem::add_square_scale(stream, w_i, M_i, M_o);
	}

	template<class TVector_c_i, class TVector_r_o, class TVector_c_o>
	enable_if_host_vector_and_host_vector<TVector_c_i, TVector_c_o, void>
	add_scale_m2psi_psi_to_host(Stream<e_host> &stream, Value_type<TVector_r_o> w_i, 
	TVector_c_i &psi_i, TVector_r_o &m2psi_o, TVector_c_o &psi_o, Vector<Value_type<TVector_c_i>, e_host> *psi_i_h =nullptr)
	{
		multem::add_scale(stream, w_i, psi_i, psi_o);
		multem::add_square_scale(stream, w_i, psi_i, m2psi_o);
	}
} // namespace multem

#endif