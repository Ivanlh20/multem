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
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HOST_FUNCTIONS_H
#define HOST_FUNCTIONS_H

#include <thread>
#include <type_traits>
#include <algorithm>
#include <algorithm>
#include <iterator>

#include <fftw3.h>
#include "math.cuh"
#include "types.hpp"
#include "traits.cuh"
#include "fft2.cuh"
#include "stream.cuh"

#include "host_device_functions.cuh"

namespace multem
{
	// get index (with typ=0: bottom index for equal values and typ=1: upper index for equal values)
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

	//symmetric coordinates(Fourier space Coordinates)
	inline
	int FSC(int i, int nh, int shift=2)
	{
		int j;
		if(shift == 1)
			j = (i < nh)?i:i-2*nh;
		else
			j = i-nh;
		return j;
	}

	// get two dimensional Hanning_Filter
	void filter_Hanning_2D(int ny, int nx, double dx, double dy, double k, int shift, double *fI)
	{	
		int nxh = nx/2, nyh = ny/2;
		double Rx, Ry, f;
		double cx = c_2Pi/(dx*(nx-1));
		double cy = c_2Pi/(dy*(ny-1));

		Vector<double, e_Host> fx(nx);
		Vector<double, e_Host> fy(ny);

		for(auto ix = 0; ix < fx.size(); ix++)
		{
			Rx = FSC(ix, nxh, shift)*dx; 
			fx[ix] = 0.5*(1.0+cos(cx*Rx));
		}
		for(auto iy = 0; iy < fy.size(); iy++)
		{
			Ry = FSC(iy, nyh, shift)*dy; 
			fy[iy] = 0.5*(1.0+cos(cy*Ry));
		}

		for(auto ix = 0; ix < fx.size(); ix++)
		{
			for(auto iy = 0; iy < fy.size(); iy++)
			{
				f = fx[ix]*fy[iy];
				if(f>k)
					fI[ix*ny+iy] = 1;
				else
					fI[ix*ny+iy] = f/k;
			}
		}
	}

	// get two dimensional Gaussian_Filter
	void filter_Gaussian_2D(int ny, int nx, double dx, double dy, double Sigma, int shift, double *fI)
	{	
		int nxh = nx/2, nyh = ny/2;
		double Rx, Ry;
		double c = 0.5/(Sigma*Sigma);

		Vector<double, e_Host> fx(nx);
		Vector<double, e_Host> fy(ny);

		for(auto ix = 0; ix < fx.size(); ix++)
		{
			Rx = FSC(ix, nxh, shift)*dx; 
			fx[ix] = exp(-c*Rx*Rx);
		}

		for(auto iy = 0; iy < fy.size(); iy++)
		{
			Ry = FSC(iy, nyh, shift)*dy; 
			fy[iy] = exp(-c*Ry*Ry);
		}

		for(auto ix = 0; ix < fx.size(); ix++)
		{
			for(auto iy = 0; iy < fy.size(); iy++)
			{
				fI[ix*ny+iy] = fx[ix]*fy[iy];
			}
		}
	}

	// get two dimensional Butterworth_Filter
	void filter_Butterworth_2D(int ny, int nx, double dx, double dy, double Radius, int n, int lpf, int shift, double *fI)
	{	
		int nxh = nx/2, nyh = ny/2;
		double R02 = pow(Radius, 2);

		for(auto ix = 0; ix < nx; ix++)
		{
			for(auto iy = 0; iy < ny; iy++)
			{
				double Rx = FSC(ix, nxh, shift)*dx; 
				double Ry = FSC(iy, nyh, shift)*dy;
				double R2 = Rx*Rx + Ry*Ry;
				fI[ix*ny+iy] = (lpf == 1)?1.0/(1.0+pow(R2/R02, n)):(isZero(R2))?0.0:1.0/(1.0+pow(R02/R2, n));
			}
		}
	}

	// get two dimensional radial distribution for regular gridBT
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

	// get two dimensional radial distribution for regular gridBT
	void getCumRadDist_2D(int ny, int nx, int shift, double *fI, int nr, double *r, double *fIr)
	{	
		int idx, nxh = nx/2, nyh = ny/2;
		Vector<double, e_Host> cfIr(nr, 0.0);
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
			for(int i=ix1; i<ix2-1; i++)
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
			for(int i=ix1; i>ix2; i--)
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
	double FFT_information_limit_2D(int ny, int nx, int shift, double *fI)
	{
		int nr = min(nx/2,ny/2)-1;
		double *r, *fIr;

		r = new double[nr];
		fIr = new double[nr]; 

		// Cumulative radial integration
		getCumRadDist_2D(ny, nx, shift, fI, nr, r, fIr);

		// Shift and Normalize
		double r0 = r[0], fIr0 = fIr[0];
		for(int i=0; i < nr; i++)
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

		for(int i=0; i < nr; i++)
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
		Vector<bool, e_Host> A_base_c(A_base_size, true);
		int size = min(A_base_size, A_search.size());
		TVector A_d;
		A_d.reserve(size);

		for(auto i=0; i<A_search.size(); i++)
		{
			T val = A_search[i]+Epsilon<float>::rel;
			auto it = std::min_element(A_base_first, A_base_last, [&val](const T &a, const T &b)->bool{ return fabs(a-val)<fabs(b-val); });
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
	void match_vectors(InputIterator A_base_first, InputIterator A_base_last, TVector &A_search, Vector<int, e_Host> &A_idx)
	{
		using T = Value_type<TVector>;

		std::size_t A_base_size = static_cast<int>(A_base_last - A_base_first);
		Vector<bool, e_Host> A_base_c(A_base_size, true);
		int size = min(A_base_size, A_search.size());
		A_idx.clear();
		A_idx.reserve(size);
		TVector A_d;
		A_d.reserve(size);

		for(auto i=0; i<A_search.size(); i++)
		{
			T val = A_search[i];
			auto it = std::min_element(A_base_first, A_base_last, [&val](const T &a, const T &b)->bool{ return fabs(a-val)<fabs(b-val); });
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

		for(auto i=0; i<p2.size(); i++)
		{
			p2[i] = std::pow((int64_t)2, i+1);

			if(i < p3.size())
			{
				p3[i] = std::pow((int64_t)3, i);
			}
			if(i < p5.size())
			{
				p5[i] = std::pow((int64_t)5, i);
			}
			if(i < p7.size())
			{
				p7[i] = std::pow((int64_t)7, i);
			}
		}

		int64_t prime, prime_0 = 128, prime_e = std::pow((int64_t)2, 18);
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
		return n_c;
	}

	void get_prime_sampling(int &nx, int &ny)
	{
		nx = get_close_Prime(nx);
		ny = get_close_Prime(ny);
	}

	/***************************************************************************/
	/***************************************************************************/

	namespace host_detail
	{
		// Linear projected potential: V and zV
		template<ePotential_Type potential_type, class T> 
		void linear_Vz(const Q1<T, e_Host> &qz, Atom_Vp<T> &atom_Vp)
		{	
			for(auto iR=0; iR < c_nR; iR++)
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
					Vr_dVrir<potential_type, T>(r, atom_Vp.cl, atom_Vp.cnl, a*qz.w[ix], V, dVir);
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
						Vr_dVrir<potential_type, T>(r, atom_Vp.cl, atom_Vp.cnl, a*qz.w[ix], V, dVir);
						V0s += V;
						dV0s += dVir;
					}
				}

				atom_Vp.c0[iR] = V0s; 		// V0
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
		void eval_cubic_poly(const Grid<T> &grid, const Atom_Vp<T> &atom_Vp, const rVector<T> &V0g)
		{	
			for(auto ix0 = 0; ix0 < atom_Vp.ixn; ix0++)
			{
				int iyc = 0;
				for(auto iy0 = 0; iy0 < atom_Vp.iyn; iy0++)
				{
					int ix = ix0 + atom_Vp.ix0;
					int iy = iy0 + atom_Vp.iy0;

					T R2 = grid.R2(ix, iy, atom_Vp.x, atom_Vp.y);
					if (R2 < atom_Vp.R_max2)
					{
						int ixy;
						T V = host_device_detail::eval_cubic_poly(ix, iy, R2, grid, atom_Vp, ixy);
						
						atom_Vp.iv[iyc] = ixy;
						atom_Vp.v[iyc] = V;
						iyc++;
					}
				}

				multem_mutex.lock();
				for(auto iy0 = 0; iy0 < iyc; iy0++)
				{
					V0g.V[atom_Vp.iv[iy0]] += atom_Vp.v[iy0];
				}
				multem_mutex.unlock();
			}
		}

		template<class TGrid>
		Value_type<TGrid> Lorentz_factor(TGrid &grid, EELS<Value_type<TGrid>> &eels)
		{
			using value_type_r = Value_type<TGrid>;

			value_type_r sum = 0.0;
			for(auto ix = 0; ix < grid.nx; ix++)
			{
				for(auto iy = 0; iy < grid.ny; iy++)
				{
					value_type_r g2 = grid.g2_shift(ix, iy);
					if(g2 < eels.gc2)
					{
						sum += 1.0/(g2 + eels.ge2);
					}
				}
			}
			return sqrt(eels.occ)/sum;
		}

	} // host_detail

	template<class TQ1>
	enable_if_Host<TQ1, void>
	get_cubic_poly_coef_Vz(ePotential_Type potential_type, TQ1 &qz, 
	Stream<Value_type<TQ1>, e_Host> &stream, Vector<Atom_Vp<Value_type<TQ1>>, e_Host> &atom_Vp)
	{
		using value_type_r = Value_type<TQ1>;

		if(stream.n_act_stream<=0)
		{
			return;
		}

		auto cubic_poly_coef = [](const ePotential_Type &potential_type, TQ1 &qz, Atom_Vp<value_type_r> &atom_Vp)
		{
			switch(potential_type)
			{
				case ePT_Doyle_0_4:
					host_detail::linear_Vz<ePT_Doyle_0_4, typename TQ1::value_type>(qz, atom_Vp);
					break;
				case ePT_Peng_0_4:
					host_detail::linear_Vz<ePT_Peng_0_4, typename TQ1::value_type>(qz, atom_Vp);
					break;
				case ePT_Peng_0_12:
					host_detail::linear_Vz<ePT_Peng_0_12, typename TQ1::value_type>(qz, atom_Vp);
					break;
				case ePT_Kirkland_0_12:
					host_detail::linear_Vz<ePT_Kirkland_0_12, typename TQ1::value_type>(qz, atom_Vp);
					break;
				case ePT_Weickenmeier_0_12:
					host_detail::linear_Vz<ePT_Weickenmeier_0_12, typename TQ1::value_type>(qz, atom_Vp);
					break;
				case ePT_Lobato_0_12:
					host_detail::linear_Vz<ePT_Lobato_0_12, typename TQ1::value_type>(qz, atom_Vp);
					break;
			}
			host_detail::cubic_poly_coef<typename TQ1::value_type>(atom_Vp); 
		};

		for(auto istream = 0; istream < stream.n_act_stream; istream++)
		{
			stream[istream] = std::thread(std::bind(cubic_poly_coef, potential_type, std::ref(qz), std::ref(atom_Vp[istream])));
		}
		stream.synchronize();
	}

	template<class TGrid, class TVector_r>
	enable_if_host_vector<TVector_r, void>
	eval_cubic_poly(TGrid &grid, Stream<Value_type<TGrid>, e_Host> &stream, 
	Vector<Atom_Vp<Value_type<TGrid>>, e_Host> &atom_Vp, TVector_r &V0)
	{
		if(stream.n_act_stream<=0)
		{
			return;
		}

		for(auto istream = 0; istream < stream.n_act_stream; istream++)
		{
			stream[istream] = std::thread(std::bind(host_detail::eval_cubic_poly<typename TVector_r::value_type>, std::ref(grid), std::ref(atom_Vp[istream]), std::ref(V0)));
		}
		stream.synchronize();
	}

	//typename std::enable_if<is_Host<TVector>::value && !is_rmatrix_c<TVector>::value, void>::type
	template<class TGrid, class TVector>
	enable_if_Host<TVector, void>
	fft2_shift(const TGrid &grid, TVector &M_io)
	{
		for(auto ix = 0; ix <grid.nxh; ix++)
		{
			for(auto iy = 0; iy <grid.nyh; iy++)			
			{
				host_device_detail::fft2_shift(ix, iy, grid, M_io);
			}
		}
	}

	template<class TGrid, class TVector>
	enable_if_host_vector<TVector, Value_type<TVector>>
	sum_over_Det(const TGrid &grid, const Value_type<TGrid> &g_min, const Value_type<TGrid> &g_max, TVector &M_i)
	{
		using value_type_r = Value_type<TGrid>;
		using value_type = Value_type<TVector>;

		value_type_r g2_min = pow(g_min, 2);
		value_type_r g2_max = pow(g_max, 2);

		value_type sum = 0.0;
		for(auto ix = 0; ix < grid.nx; ix++)
		{
			for(auto iy = 0; iy < grid.ny; iy++)
			{				
				value_type_r g2 = grid.g2_shift(ix, iy);
				if((g2_min <= g2)&&(g2 <= g2_max))
				{
					int ixy = grid.ind_col(ix, iy); 
					sum += M_i[ixy];
				}
			}
		}

		return sum;
	}

	template<class TGrid, class TVector_r>
	enable_if_host_vector<TVector_r, Value_type<TGrid>>
	sum_square_over_Det(const TGrid &grid, const Value_type<TGrid> &g_min, const Value_type<TGrid> &g_max, TVector_r &M_i)
	{
		using value_type_r = Value_type<TGrid>;

		value_type_r g2_min = pow(g_min, 2);
		value_type_r g2_max = pow(g_max, 2);

		value_type_r sum = 0;
		for(auto ix = 0; ix < grid.nx; ix++)
		{
			for(auto iy = 0; iy < grid.ny; iy++)
			{	
				value_type_r g2 = grid.g2_shift(ix, iy);
				if((g2_min <= g2)&&(g2 <= g2_max))
				{
					int ixy = grid.ind_col(ix, iy);
					sum += thrust::norm(M_i[ixy]);
				}
			}
		}

		return sum;
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	bandwidth_limit(const TGrid &grid, const Value_type<TGrid> &g_min, const Value_type<TGrid> &g_max, Value_type<TVector_c> w, TVector_c &M_io)
	{
		using value_type_r = Value_type<TGrid>;

		value_type_r g2_min = pow(g_min, 2);
		value_type_r g2_max = pow(g_max, 2);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			for(auto iy = 0; iy < grid.ny; iy++)
			{
				host_device_detail::bandwidth_limit(ix, iy, grid, g2_min, g2_max, w, M_io);
			}
		}
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	phase_multiplication(const TGrid &grid, TVector_c &exp_x_i, TVector_c &exp_y_i, TVector_c &psi_i, TVector_c &psi_o)
	{
		for(auto ix = 0; ix < grid.nx; ix++)
		{
			for(auto iy = 0; iy < grid.ny; iy++)
			{
				host_device_detail::phase_multiplication(ix, iy, grid, exp_x_i, exp_y_i, psi_i, psi_o);
			}
		}
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	propagator_multiplication(const TGrid &grid, TVector_c &prop_x_i, TVector_c &prop_y_i, TVector_c &psi_i, TVector_c &psi_o)
	{
		for(auto ix = 0; ix < grid.nx; ix++)
		{
			for(auto iy = 0; iy < grid.ny; iy++)
			{
				host_device_detail::propagator_multiplication(ix, iy, grid, prop_x_i, prop_y_i, psi_i, psi_o);
			}
		}
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	probe(const TGrid &grid, const Lens<Value_type<TGrid>> &lens, Value_type<TGrid> x, Value_type<TGrid> y, TVector_c &fPsi_o)
	{
		using value_type_r = Value_type<TGrid>;

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			for(auto iy = 0; iy < grid.ny; iy++)
			{
				host_device_detail::probe(ix, iy, grid, lens, x, y, fPsi_o);
			}
		}

		value_type_r total = sum_square(grid, fPsi_o);
		multem::scale(fPsi_o, sqrt(value_type_r(grid.nxy())/total));
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	apply_CTF(const TGrid &grid, const Lens<Value_type<TGrid>> &lens, Value_type<TGrid> gxu, Value_type<TGrid> gyu, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		for(auto ix = 0; ix < grid.nx; ix++)
		{
			for(auto iy = 0; iy < grid.ny; iy++)
			{
				host_device_detail::apply_CTF(ix, iy, grid, lens, gxu, gyu, fPsi_i, fPsi_o);
			}
		}
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	apply_PCTF(const TGrid &grid, const Lens<Value_type<TGrid>> &lens, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		for(auto ix = 0; ix < grid.nx; ix++)
		{
			for(auto iy = 0; iy < grid.ny; iy++)
			{
				host_device_detail::apply_PCTF(ix, iy, grid, lens, fPsi_i, fPsi_o);
			}
		}
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	kernel_xyz(const TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Host> &fft2, TVector_c &k_x, TVector_c &k_y, TVector_c &k_z)
	{
		eels.factor = host_detail::Lorentz_factor(grid, eels);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			for(auto iy = 0; iy < grid.ny; iy++)
			{
				host_device_detail::kernel_xyz(ix, iy, grid, eels, k_x, k_y, k_z);
			}
		}

		fft2.inverse(k_x);
		fft2.inverse(k_y);
		fft2.inverse(k_z);
	}
	
	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	kernel_x(const TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Host> &fft2, TVector_c &k_x)
	{
		eels.factor = host_detail::Lorentz_factor(grid, eels);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			for(auto iy = 0; iy < grid.ny; iy++)
			{
				host_device_detail::kernel_x(ix, iy, grid, eels, k_x);
			}
		}

		fft2.inverse(k_x);
	}
	
	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	kernel_y(const TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Host> &fft2, TVector_c &k_y)
	{
		eels.factor = host_detail::Lorentz_factor(grid, eels);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			for(auto iy = 0; iy < grid.ny; iy++)
			{
				host_device_detail::kernel_y(ix, iy, grid, eels, k_y);
			}
		}

		fft2.inverse(k_y);
	}
	
	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	kernel_z(const TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Host> &fft2, TVector_c &k_z)
	{
		eels.factor = host_detail::Lorentz_factor(grid, eels);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			for(auto iy = 0; iy < grid.ny; iy++)
			{
				host_device_detail::kernel_z(ix, iy, grid, eels, k_z);
			}
		}

		fft2.inverse(k_z);
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	kernel_mn1(const TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Host> &fft2, TVector_c &k_mn1)
	{
		eels.factor = host_detail::Lorentz_factor(grid, eels);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			for(auto iy = 0; iy < grid.ny; iy++)
			{
				host_device_detail::kernel_mn1(ix, iy, grid, eels, k_mn1);
			}
		}

		fft2.inverse(k_mn1);
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	kernel_mp1(const TGrid &grid, EELS<Value_type<TGrid>> &eels, FFT2<Value_type<TGrid>, e_Host> &fft2, TVector_c &k_mp1)
	{
		eels.factor = host_detail::Lorentz_factor(grid, eels);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			for(auto iy = 0; iy < grid.ny; iy++)
			{
				host_device_detail::kernel_mp1(ix, iy, grid, eels, k_mp1);
			}
		}

		fft2.inverse(k_mp1);
	}

	template<class TVector_c>
	enable_if_host_vector<TVector_c, void>
	rmatrix_c_to_complex(rmatrix_c &M_i, TVector_c &M_o)
	{
		for(auto ixy = 0; ixy < M_i.size; ixy++)
		{
			M_o[ixy] = M_i[ixy];
		}
	}

	/***************************************************************************/
	/***************************************************************************/

	template<class TGrid, class TVector_i, class TVector_o>
	enable_if_host_vector<TVector_i, void>
	copy_to_host(const TGrid &grid, TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_Host> *M_i_h=nullptr)
	{
		for(auto ixy = 0; ixy < M_i.size(); ixy++)
		{
			M_o[ixy] = M_i[ixy];
		}
	}

	template<class TGrid, class TVector_i, class TVector_o>
	enable_if_host_vector<TVector_i, void>
	add_scale_to_host(TGrid &grid, Value_type<TVector_i> w_i, TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_Host> *M_i_h=nullptr)
	{
		for(auto i = 0; i < M_i_h->size(); i++)
		{
			M_o[i] += w_i*M_i[i];
		}
	}

	template<class TGrid, class TVector_c_i, class TVector_r_o, class TVector_c_o>
	enable_if_host_vector<TVector_c_i, void>
	add_scale_m2psi_psi(TGrid &grid, Value_type<TGrid> w_i, TVector_c_i &psi_i, TVector_r_o &m2psi_o, TVector_c_o &psi_o, Vector<Value_type<TVector_c_i>, e_Host> *psi_i_h=nullptr)
	{
		Value_type<TVector_c_i> w_i_c = w_i;
		for(auto i = 0; i < psi_i.size(); i++)
		{
			m2psi_o[i] += w_i*thrust::norm(psi_i[i]);
			psi_o[i] += w_i_c*psi_i[i];
		}
	}
} //namespace multem

#endif