/*
 * This file is part of MULTEM.
 * Copyright 2016 Ivan Lobato <Ivanlh20@gmail.com>
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

namespace mt
{
	// median filter 1d
	template<class TVector>
	TVector filter_median_1d(Stream<e_host> &stream, TVector &Im_i, int nkr);

	// median filter 2d
	template<class TGrid, class TVector>
	TVector filter_median_2d(Stream<e_host> &stream, TGrid &grid, TVector &Im_i, int nkr);

	// wiener filter 1d
	template<class TVector>
	TVector filter_wiener_1d(Stream<e_host> &stream, TVector &Im_i, int nkr);

	// wiener filter 2d
	template<class TGrid, class TVector>
	TVector filter_wiener_2d(Stream<e_host> &stream, TGrid &grid, TVector &Im_i, int nkr);

	// denoising poisson 
	template<class TVector>
	TVector filter_denoising_poisson_1d(Stream<e_host> &stream, TVector &Im_i, int nkr_w, int nkr_m);

	// get peak signal to noise ratio PSNR 
	template<class TVector>
	Value_type<TVector> get_PSNR(Stream<e_host> &stream, TVector &Im_i, TVector &Im_d);

	// gray opening
	template<class TVector>
	TVector morp_g_open(Stream<e_host> &stream, int ny_i, int nx_i, TVector &Im_i, int nkr);

	// thresholding
	template<class TVector>
	TVector thresholding(Stream<e_host> &stream, TVector &v_i, Value_type<TVector> threshold);

	template<class T>
	class Neigh_2d;

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
		void matrix_iter_yx(const Range &range, TFn &fn, TArg &...arg)
		{
			for(auto iy = range.iy_0; iy < range.iy_e; iy++)
			{
				for(auto ix = range.ix_0; ix < range.ix_e; ix++)
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
		T atom_cost_function(const Grid<T> &grid, const Atom_Sa<T> &atom_Ip, rVector<T> M_i)
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
		void subtract_atom(Stream<e_host> &stream, const Grid<T> &grid, const Atom_Sa<T> &atom_Ip, rVector<T> M_i)
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
					M_i.V[atom_Ip.iv[iy0]] -= atom_Ip.v[iy0];
				}
				stream.stream_mutex.unlock();
			}
		}

		// Linear projected potential: V and zV
		template<ePotential_Type potential_type, int charge, class TAtom>
		void linear_Vz(const Q1<Value_type<TAtom>, e_host> &qz, TAtom &atom)
		{	
			using T = Value_type<TAtom>;

			for(auto iR = 0; iR < c_nR; iR++)
			{
				T R2 = atom.R2[iR];
				T V = 0;
				T dVir = 0;

				T a = (atom.split)?(-atom.z0h):(atom.zeh-atom.z0h);
				T b = (atom.split)?(atom.z0h):(atom.zeh+atom.z0h);
				for(auto ix = 0; ix < qz.size(); ix++)
				{
					T z = a*qz.x[ix] + b;
					T r = sqrt(z*z + R2);
					T V0s, dV0s;
					Vr_dVrir<potential_type, charge, T>(r, atom.cl, atom.cnl, a*qz.w[ix], V0s, dV0s);
					V += V0s;
					dVir += dV0s;
				}

				if (atom.split)
				{
					T a = atom.zeh;
					T b = atom.zeh;
					for(auto ix = 0; ix < qz.size(); ix++)
					{
						T z = a*qz.x[ix] + b;
						T r = sqrt(z*z + R2);
						T V0s, dV0s;
						Vr_dVrir<potential_type, charge, T>(r, atom.cl, atom.cnl, a*qz.w[ix], V0s, dV0s);
						V += V0s;
						dVir += dV0s;
					}
				}

				dVir = 0.5*dVir;

				auto R2_tap = atom.R2_tap;
				auto tap_cf = atom.tap_cf;
				host_device_detail::apply_tapering(R2_tap, tap_cf, R2, V, dVir);
				atom.c0[iR] = V;			// V_0
				atom.c1[iR] = dVir; 		// dR2V0
			}
		}

		// Get Local interpolation coefficients
		template<class TAtom> 
		void cubic_poly_coef(TAtom &atom)
		{
			for(auto iR = 0; iR < c_nR-1; iR++)
			{
				host_device_detail::cubic_poly_coef(iR, atom);
			}
		}

		// Cubic polynomial evaluation
		template<bool TShift, class TAtom> 
		void eval_cubic_poly(Stream<e_host> &stream, const Grid<Value_type<TAtom>> &grid, const TAtom &atom, const rVector<Value_type<TAtom>> &M_o)
		{	
			using T = Value_type<TAtom>;

			for(auto ix0 = 0; ix0 < atom.ixn; ix0++)
			{
				int iyc = 0;
				for(auto iy0 = 0; iy0 < atom.iyn; iy0++)
				{
					int ix = ix0 + atom.ix0;
					int iy = iy0 + atom.iy0;

					T R2 = grid.R2(ix, iy, atom.x, atom.y);
					if (R2 < atom.R2_max)
					{
						int ixy;
						T V = host_device_detail::eval_cubic_poly<TShift>(ix, iy, grid, R2, atom, ixy);
						
						atom.iv[iyc] = ixy;
						atom.v[iyc] = V;
						iyc++;
					}
				}

				stream.stream_mutex.lock();
				for(auto iy0 = 0; iy0 < iyc; iy0++)
				{
					M_o.V[atom.iv[iy0]] += atom.v[iy0];
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

		// Linear Gaussian
		template<class TAtom>
		void linear_Gaussian(TAtom &atom)
		{	
			for(auto iR = 0; iR < c_nR; iR++)
			{
				host_device_detail::linear_Gaussian(iR, atom);
			}
		}

		template<class TVector>
		void fit_log_gaussian_a_sigma(const TVector &x2_i, const TVector &y_i, Value_type<TVector> &a1, Value_type<TVector> &a0)
		{
			using T = Value_type<TVector>;

			// get a, sigma
			T sx1x1 = 0;
			T sx2x2 = 0;
			T sx1x2 = 0;
			T sx1y = 0;
			T sx2y = 0;
			for (auto ix = 0; ix < x2_i.size(); ix++)
			{
				T f = x2_i[ix];
				T x1 = f;
				T x2 = 1;
				T y = y_i[ix];

				sx1x1 += x1*x1;
				sx2x2 += x2*x2;
				sx1x2 += x1*x2;
				sx1y += x1*y;
				sx2y += x2*y;
			}

			T det = sx1x1*sx2x2-sx1x2*sx1x2;
			a1 = (sx2x2*sx1y-sx1x2*sx2y)/det;
			a0 = (sx1x1*sx2y-sx1x2*sx1y)/det;
		}

		template<class TVector>
		void fit_gaussian_a_c(const TVector &x2_i, const TVector &y_i, Value_type<TVector> lambdai, 
		Value_type<TVector> sigma, Value_type<TVector> &a, Value_type<TVector> &c)
		{
			using T = Value_type<TVector>;

			// get a, c
			T c_0 = 0.5/pow(sigma, 2);
			T sx1x1 = 0;
			T sx2x2 = 0;
			T sx1x2 = 0;
			T sx1y = 0;
			T sx2y = 0;
			for (auto ix = 0; ix < x2_i.size(); ix++)
			{
				T f = exp(-c_0*x2_i[ix]);
				T x1 = f;
				T x2 = 1;
				T y = y_i[ix];

				sx1x1 += x1*x1;
				sx2x2 += x2*x2;
				sx1x2 += x1*x2;
				sx1y += x1*y;
				sx2y += x2*y;
			}
			T lambda = lambdai*(sx1x1+sx2x2);
			sx1x1 += lambda;
			sx2x2 += lambda;

			T det = sx1x1*sx2x2-sx1x2*sx1x2;
			a = (sx2x2*sx1y-sx1x2*sx2y)/det;
			c = (sx1x1*sx2y-sx1x2*sx1y)/det;
		}

		template<class TVector>
		Value_type<TVector> fit_gaussian_a(const TVector &x2_i, const TVector &y_i, Value_type<TVector> sigma)
		{
			using T = Value_type<TVector>;

			// get a = sum xy/sum x^2
			T c_0 = 0.5/pow(sigma, 2);

			T sx1x1 = 0;
			T sx1y = 0;
			for (auto ix = 0; ix < x2_i.size(); ix++)
			{
				T x1 = exp(-c_0*x2_i[ix]);
				T y = y_i[ix];

				sx1x1 += x1*x1;
				sx1y += x1*y;
			}

			return sx1y/sx1x1;
		}

		template<class TVector>
		Value_type<TVector> fit_gaussian_sigma(const TVector &x2_i, const TVector &y_i, Value_type<TVector> a, 
		Value_type<TVector> &sigma, Value_type<TVector> sigma_min, Value_type<TVector> sigma_max,
		int nit, Value_type<TVector> d_sigma_error)
		{
			using T = Value_type<TVector>;

			T sigma_o = sigma;
			// get b = sum xy/sum x^2
			for (auto it = 0; it < nit; it++)
			{
				T c_0 = 0.5/pow(sigma, 2);
				T c_1 = a/pow(sigma, 3);

				T sx1x1 = 0;
				T sx1y = 0;
				for (auto ix = 0; ix < x2_i.size(); ix++)
				{
					T f = exp(-c_0*x2_i[ix]);
					T x1 = c_1*x2_i[ix]*f;
					T y = y_i[ix]-a*f;

					sx1x1 += x1*x1;
					sx1y += x1*y;
				}
				T d_sigma = sx1y/sx1x1;
				sigma += d_sigma;
				//sigma = min(max(sigma, sigma_min), sigma_max);

				if(sigma <= sigma_min)
				{
					sigma = sigma_min;
					break;
				}

				if(sigma >= sigma_max)
				{
					sigma = sigma_max;
					break;
				}

				if (abs(d_sigma/sigma)<d_sigma_error)
				{
					break;
				}
			}

			return sigma-sigma_o;
		}

		template<class TVector>
		Value_type<TVector> fit_gaussian_w_a(const TVector &x2_i, const TVector &y_i, Value_type<TVector> sigma)
		{
			using T = Value_type<TVector>;

			// get a = sum xy/sum x^2
			T c_0 = 0.5/pow(sigma, 2);

			T sx1x1 = 0;
			T sx1y = 0;
			for (auto ix = 0; ix < x2_i.size(); ix++)
			{
				//T w = 1/::fmax(x2_i[ix], 0.001);
				T w = x2_i[ix];
				T x1 = exp(-c_0*x2_i[ix]);
				T y = w*y_i[ix];

				sx1x1 += w*x1*x1;
				sx1y += x1*y;
			}

			return sx1y/sx1x1;
		}

		template<class TVector>
		Value_type<TVector> fit_gaussian_w_sigma(const TVector &x2_i, const TVector &y_i, Value_type<TVector> a, 
		Value_type<TVector> &sigma, Value_type<TVector> sigma_min, Value_type<TVector> sigma_max,
		int nit, Value_type<TVector> d_sigma_error)
		{
			using T = Value_type<TVector>;

			T sigma_o = sigma;
			// get b = sum xy/sum x^2
			for (auto it = 0; it < nit; it++)
			{
				T c_0 = 0.5/pow(sigma, 2);
				T c_1 = a/pow(sigma, 3);

				T sx1x1 = 0;
				T sx1y = 0;
				for (auto ix = 0; ix < x2_i.size(); ix++)
				{
					//T w = 1/::fmax(x2_i[ix], 0.001);
					T w = x2_i[ix];
					T f = exp(-c_0*x2_i[ix]);
					T x1 = c_1*x2_i[ix]*f;
					T y = y_i[ix]-a*f;

					sx1x1 += w*x1*x1;
					sx1y += w*x1*y;
				}
				T d_sigma = sx1y/sx1x1;
				sigma += d_sigma;

				if(sigma <= sigma_min)
				{
					sigma = sigma_min;
					break;
				}

				if(sigma >= sigma_max)
				{
					sigma = sigma_max;
					break;
				}

				if (abs(d_sigma/sigma)<d_sigma_error)
				{
					break;
				}
			}

			return sigma-sigma_o;
		}

	} // host_detail

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

	template<class TVector_1, class TVector_2>
	typename std::enable_if<is_host_vector_and_host_vector<TVector_1, TVector_2>::value && is_complex<Value_type<TVector_1>>::value, void>::type
	assign_real(Stream<e_host> &stream, TVector_1 &M_i, TVector_2 &M_o, Vector<Value_type<TVector_2>, e_host> *M_i_h =nullptr)
	{
		M_o.resize(M_i.size());
		auto thr_assign_real = [&](const Range &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				M_o[ixy] = M_i[ixy].real();
			}
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_assign_real);
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
	enable_if_host_vector<TVector, Value_type<TVector>>
	min_element(Stream<e_host> &stream, TVector &M_i)
	{
		using T = Value_type<TVector>;

		vector<T> v_min;
		v_min.reserve(stream.size());

		auto thr_min = [&](const Range &range)
		{
			auto m_min = *std::min_element(M_i.begin()+range.ixy_0, M_i.begin()+range.ixy_e);
			stream.stream_mutex.lock();
			v_min.push_back(m_min);
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(M_i.size());
		stream.set_grid(1, M_i.size());
		stream.exec(thr_min);

		return *std::min_element(v_min.begin(), v_min.end());
	}

	template<class TVector>
	enable_if_host_vector<TVector, Value_type<TVector>>
	max_element(Stream<e_host> &stream, TVector &M_i)
	{
		using T = Value_type<TVector>;

		vector<T> v_max;
		v_max.reserve(stream.size());

		auto thr_max = [&](const Range &range)
		{
			T m_max = *std::max_element(M_i.begin()+range.ixy_0, M_i.begin()+range.ixy_e);
			stream.stream_mutex.lock();
			v_max.push_back(m_max);
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(M_i.size());
		stream.set_grid(1, M_i.size());
		stream.exec(thr_max);

		return *std::max_element(v_max.begin(), v_max.end());
	}

	template<class TVector>
	enable_if_host_vector<TVector, Value_type_r<TVector>>
	mean(TVector &M_i)
	{
		return thrust::reduce(M_i.begin(), M_i.end())/M_i.size();
	}

	template<class TVector>
	enable_if_host_vector<TVector, Value_type_r<TVector>>
	mean(Stream<e_host> &stream, TVector &M_i)
	{
		return sum(stream, M_i)/M_i.size();
	}

	template<class TVector>
	enable_if_host_vector<TVector, void>
	mean_var(TVector &M_i, Value_type_r<TVector> &x_mean, Value_type_r<TVector> &x_var)
	{
		using value_type_r = Value_type_r<TVector>;

		x_mean = mean(M_i);

		x_var = thrust::transform_reduce(M_i.begin(), M_i.end(), 
		functor::square_dif<value_type_r>(x_mean), value_type_r(0), functor::add<value_type_r>());

		x_var = x_var/M_i.size();
	}

	template<class TVector>
	enable_if_host_vector<TVector, void>
	mean_var(Stream<e_host> &stream, TVector &M_i, Value_type_r<TVector> &x_mean, Value_type_r<TVector> &x_var)
	{
		using value_type_r = Value_type_r<TVector>;

		x_mean = mean(stream, M_i);

		x_var = 0;
		auto thr_var = [&](const Range &range)
		{
			auto x_var_partial = thrust::transform_reduce(M_i.begin()+range.ixy_0, M_i.begin()+range.ixy_e, 
			functor::square_dif<value_type_r>(x_mean), value_type_r(0), functor::add<value_type_r>());

			stream.stream_mutex.lock();
			x_var += x_var_partial;
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(M_i.size());
		stream.set_grid(1, M_i.size());
		stream.exec(thr_var);

		x_var = x_var/M_i.size();
	}

	template<class TVector>
	enable_if_host_vector<TVector, Value_type_r<TVector>>
	variance(TVector &M_i)
	{
		using value_type_r = Value_type_r<TVector>;

		value_type_r x_mean , x_var;
		mean_var(M_i, x_mean, x_var);
		return x_var;
	}

	template<class TVector>
	enable_if_host_vector<TVector, Value_type_r<TVector>>
	variance(Stream<e_host> &stream, TVector &M_i)
	{
		using value_type_r = Value_type_r<TVector>;

		value_type_r x_mean , x_var;
		mean_var(stream, M_i, x_mean, x_var);
		return x_var;
	}

	template<class TVector>
	void rescale_data(const Value_type<TVector> &x_mean, const Value_type<TVector> &x_std, TVector &x)
	{
		using T = Value_type<TVector>;
		std::for_each(x.begin(), x.end(), [=](T &v){ v = (v-x_mean)/x_std; });
	}

	template<class TVector, class ...TArg>
	void rescale_data(const Value_type<TVector> &x_mean, const Value_type<TVector> &x_std, TVector &x, TArg &...arg)
	{
		rescale_data(x_mean, x_std, x);
		rescale_data(x_mean, x_std, arg...);
	}

	template<class TVector>
	Value_type<TVector> min_data(TVector &x)
	{
		return *std::min_element(x.begin(), x.end());
	}

	template<class TVector, class ...TArg>
	Value_type<TVector> min_data(TVector &x, TArg &...arg)
	{
		return ::fmin(min_data(x), min_data(arg...));
	}

	template<class TVector>
	Value_type<TVector> max_data(TVector &x)
	{
		return *std::max_element(x.begin(), x.end());
	}

	template<class TVector, class ...TArg>
	Value_type<TVector> max_data(TVector &x, TArg &...arg)
	{
		return ::fmax(max_data(x), max_data(arg...));
	}

	template<class TVector>
	Value_type<TVector> max_std_data(TVector &x)
	{
		return sqrt(variance(x));
	}

	template<class TVector, class ...TArg>
	Value_type<TVector> max_std_data(TVector &x, TArg &...arg)
	{
		return ::fmax(max_std_data(x), max_std_data(arg...));
	}
	
	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	phase_factor_1d(Stream<e_host> &stream, TGrid &grid, Value_type<TGrid> x, TVector_c &fPsi_i, TVector_c &fPsi_o)
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
	phase_factor_2d(Stream<e_host> &stream, TGrid &grid, Value_type<TGrid> x, Value_type<TGrid> y, 
	TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		auto thr_phase_factor = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::phase_factor_2d<TGrid, TVector_c>, grid, x, y, fPsi_i, fPsi_o);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_phase_factor);
	}

	template<class TGrid, class TVector_r>
	enable_if_host_vector<TVector_r, Value_type<TVector_r>>
	atom_cost_function(TGrid &grid, const Atom_Sa<Value_type<TGrid>> &atom_Ip, TVector_r &M_i)
	{
		return host_detail::atom_cost_function<typename TVector_r::value_type>(grid, atom_Ip, M_i);
	}

	template<class TGrid, class TVector_r>
	enable_if_host_vector<TVector_r, void>
	subtract_atom(Stream<e_host> &stream, TGrid &grid, Vector<Atom_Sa<Value_type<TGrid>>, e_host> &atom_Ip, TVector_r &M_i)
	{
		if(stream.n_act_stream<= 0)
		{
			return;
		}

		for(auto istream = 0; istream < stream.n_act_stream; istream++)
		{
			stream[istream] = std::thread(std::bind(host_detail::subtract_atom<typename TVector_r::value_type>, std::ref(stream), std::ref(grid), std::ref(atom_Ip[istream]), std::ref(M_i)));
		}
		stream.synchronize();
	}

	// Linear projected potential: V and zV
	template<class TQ1, class TVAtom>
	enable_if_host<TQ1, void>
	linear_Vz(Stream<e_host> &stream, ePotential_Type potential_type, TQ1 &qz, TVAtom &vatom)
	{
		using TAtom = Value_type<TVAtom>;

		if(stream.n_act_stream<= 0)
		{
			return;
		}

		auto thr_linear_Vz = [](const ePotential_Type &potential_type, TQ1 &qz, TAtom &atom)
		{
			if(atom.charge== 0)
			{
				switch(potential_type)
				{
					case ePT_Doyle_0_4:
						host_detail::linear_Vz<ePT_Doyle_0_4, 0, TAtom>(qz, atom);
						break;
					case ePT_Peng_0_4:
						host_detail::linear_Vz<ePT_Peng_0_4, 0, TAtom>(qz, atom);
						break;
					case ePT_Peng_0_12:
						host_detail::linear_Vz<ePT_Peng_0_12, 0, TAtom>(qz, atom);
						break;
					case ePT_Kirkland_0_12:
						host_detail::linear_Vz<ePT_Kirkland_0_12, 0, TAtom>(qz, atom);
						break;
					case ePT_Weickenmeier_0_12:
						host_detail::linear_Vz<ePT_Weickenmeier_0_12, 0, TAtom>(qz, atom);
						break;
					case ePT_Lobato_0_12:
						host_detail::linear_Vz<ePT_Lobato_0_12, 0, TAtom>(qz, atom);
						break;
				}
			}
			else
			{
				switch(potential_type)
				{
					case ePT_Doyle_0_4:
						host_detail::linear_Vz<ePT_Doyle_0_4, 1, TAtom>(qz, atom);
						break;
					case ePT_Peng_0_4:
						host_detail::linear_Vz<ePT_Peng_0_4, 1, TAtom>(qz, atom);
						break;
					case ePT_Peng_0_12:
						host_detail::linear_Vz<ePT_Peng_0_12, 1, TAtom>(qz, atom);
						break;
					case ePT_Kirkland_0_12:
						host_detail::linear_Vz<ePT_Kirkland_0_12, 1, TAtom>(qz, atom);
						break;
					case ePT_Weickenmeier_0_12:
						host_detail::linear_Vz<ePT_Weickenmeier_0_12, 1, TAtom>(qz, atom);
						break;
					case ePT_Lobato_0_12:
						host_detail::linear_Vz<ePT_Lobato_0_12, 1, TAtom>(qz, atom);
						break;
				}
			}
		};

		for(auto istream = 0; istream < stream.n_act_stream-1; istream++)
		{
			stream[istream] = std::thread(std::bind(thr_linear_Vz, potential_type, std::ref(qz), std::ref(vatom[istream])));
		}

		thr_linear_Vz(potential_type, qz, vatom[stream.n_act_stream-1]);

		stream.synchronize();
	}

	// Get Local interpolation coefficients
	template<class TVAtom> 
	enable_if_host<typename TVAtom::value_type, void>
	cubic_poly_coef(Stream<e_host> &stream, TVAtom &vatom)
	{
		using TAtom = Value_type<TVAtom>;

		if(stream.n_act_stream<= 0)
		{
			return;
		}

		for(auto istream = 0; istream < stream.n_act_stream-1; istream++)
		{
			stream[istream] = std::thread(std::bind(host_detail::cubic_poly_coef<TAtom>, std::ref(vatom[istream])));
		}

		host_detail::cubic_poly_coef<TAtom>(vatom[stream.n_act_stream-1]);

		stream.synchronize();
	}

	// Cubic polynomial evaluation
	template<bool TShift, class TGrid, class TVAtom, class TVector_r>
	enable_if_host_vector<TVector_r, void>
	eval_cubic_poly(Stream<e_host> &stream, TGrid &grid, TVAtom &vatom, TVector_r &M_o)
	{
		using TAtom = Value_type<TVAtom>;

		if(stream.n_act_stream<= 0)
		{
			return;
		}

		for(auto istream = 0; istream < stream.n_act_stream-1; istream++)
		{
			stream[istream] = std::thread(std::bind(host_detail::eval_cubic_poly<TShift, TAtom>, std::ref(stream), std::ref(grid), std::ref(vatom[istream]), std::ref(M_o)));
		}

		host_detail::eval_cubic_poly<TShift, TAtom>(stream, grid, vatom[stream.n_act_stream-1], M_o);

		stream.synchronize();
	}

	// Get linear Gaussian
	template<class TVAtom> 
	enable_if_host<typename TVAtom::value_type, void>
	linear_Gaussian(Stream<e_host> &stream, TVAtom &vatom)
	{
		using TAtom = Value_type<TVAtom>;

		if(stream.n_act_stream<= 0)
		{
			return;
		}

		for(auto istream = 0; istream < stream.n_act_stream-1; istream++)
		{
			stream[istream] = std::thread(std::bind(host_detail::linear_Gaussian<TAtom>, std::ref(vatom[istream])));
		}

		host_detail::linear_Gaussian<TAtom>(vatom[stream.n_act_stream-1]);

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
	bandwidth_limit(Stream<e_host> &stream, TGrid &grid, TVector_c &M_io)
	{
		using value_type_r = Value_type<TGrid>;

		auto thr_bandwidth_limit = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::bandwidth_limit<TGrid, TVector_c>, grid, M_io);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_bandwidth_limit);
	}

	template<class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
	hard_aperture(Stream<e_host> &stream, TGrid &grid, Value_type<TGrid> g_max, Value_type<TGrid> w, TVector_c &M_io)
	{
		using value_type_r = Value_type<TGrid>;

		value_type_r g2_max = pow(g_max, 2);

		auto thr_hard_aperture = [&](const Range &range)
		{
			host_detail::matrix_iter(range, host_device_detail::hard_aperture<TGrid, TVector_c>, grid, g2_max, w, M_io);
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_hard_aperture);
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
		mt::scale(stream, sqrt(1.0/total), fPsi_o);
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
		mt::assign(stream, M_i, M_o, M_i_h);
	}

	template<class TVector_i, class TVector_o>
	enable_if_host_vector_and_host_vector<TVector_i, TVector_o, void>
	add_scale_to_host(Stream<e_host> &stream, Value_type<TVector_i> w_i, 
	TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_host> *M_i_h =nullptr)
	{
		mt::add_scale(stream, w_i, M_i, M_o);
	}

	template<class TVector_i, class TVector_o>
	enable_if_host_vector_and_host_vector<TVector_i, TVector_o, void>
	add_square_scale_to_host(Stream<e_host> &stream, Value_type<TVector_o> w_i, 
	TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_host> *M_i_h =nullptr)
	{
		mt::add_square_scale(stream, w_i, M_i, M_o);
	}

	template<class TVector_c_i, class TVector_r_o, class TVector_c_o>
	enable_if_host_vector_and_host_vector<TVector_c_i, TVector_c_o, void>
	add_scale_m2psi_psi_to_host(Stream<e_host> &stream, Value_type<TVector_r_o> w_i, 
	TVector_c_i &psi_i, TVector_r_o &m2psi_o, TVector_c_o &psi_o, Vector<Value_type<TVector_c_i>, e_host> *psi_i_h =nullptr)
	{
		mt::add_scale(stream, w_i, psi_i, psi_o);
		mt::add_square_scale(stream, w_i, psi_i, m2psi_o);
	}

	/***************************************************************************/
	/***************************************************************************/

	template <class TVector>
	TVector transpose(Stream<e_host> &stream, const int &nrows, const int &ncols, TVector &M)
	{
		TVector M_t(nrows*ncols);
		auto thr_transpose = [nrows, ncols](const Range &range, TVector &M, TVector &M_t)
		{
			for(auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for(auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					int ixy = ix*nrows+iy;
					int ixy_t = iy*ncols+ix;
					M_t[ixy_t] = M[ixy];
				}
			}
		};

		stream.set_n_act_stream(ncols);
		stream.set_grid(ncols, nrows);
		stream.exec(thr_transpose, M, M_t);

		return M_t;
	}

	/***************************************************************************/
	/***************************************************************************/

	// get index (with typ = 0: bottom index for equal values and typ = 1: upper index for equal values)
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

	template <class TVector>
	TVector lsf_poly_n(TVector &x, TVector &y, int n)
	{
		using T = Value_type<TVector>;

		int m = x.size();
		n++; // add bias term

		TVector M(m*n);
		for(auto in = 0; in < n; in++)
		{
			for(auto im = 0; im<m; im++)
			{
				M[in*m+im] = (in==0)?1:x[im]*M[(in-1)*m+im];
			}
		}

		TVector coef(n);
		lapack::GELS<T> gels;
		gels(m, n, M.data(), 1, y.data(), coef.data());

		return coef;
	}

	/********************************************************************/
	// get one dimensional Hanning Filter
	template<class TVector>
	TVector func_hanning_1d(Stream<e_host> &stream, int nx, Value_type<TVector> dx, Value_type<TVector> k, bool shift,
	Value_type<TVector> x_0=0, Value_type<TVector> x_e=0)
	{	
		using T = Value_type<TVector>;

		int nxh = nx/2;

		const T lx = dx*nx - (x_0 + x_e);
		const T cx = c_2Pi/lx;
		x_e = lx + x_0;

		TVector f(nx);

		auto thr_hanning = [&](const Range &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				T Rx = RSC(ixy, nxh, shift)*dx; 
				T v = ((Rx<x_0)||(Rx>x_e))?0:0.5*(1.0-cos(cx*(Rx-x_0)));
				f[ixy] = (v>k)?1.0:v/k;
			}
		};

		stream.set_n_act_stream(nx);
		stream.set_grid(nx, 1);
		stream.exec(thr_hanning);

		return f;
	}

	// get two dimensional Hanning by row
	template<class TVector, class TGrid>
	TVector func_hanning_2d_by_row(Stream<e_host> &stream, TGrid &grid, Value_type<TVector> k, bool shift, 
	Value_type<TVector> x_0=0, Value_type<TVector> x_e=0)
	{	
		using T = Value_type<TVector>;

		const T lx = grid.lx - (x_0 + x_e);
		const T cx = c_2Pi/lx;
		x_e = lx + x_0;

		TVector fx;
		fx.reserve(grid.nx);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			T Rx = (shift)?grid.Rx_shift(ix):grid.Rx(ix);
			T v = ((Rx<x_0)||(Rx>x_e))?0:0.5*(1.0-cos(cx*(Rx-x_0)));
			fx.push_back(v);
		}

		TVector f(grid.nxy());

		auto thr_hanning = [&](const Range &range)
		{
			for(auto iy = range.ixy_0; iy < range.ixy_e; iy++)
			{
				for(auto ix = 0; ix < grid.nx; ix++)
				{
					T v = fx[ix];
					f[grid.ind_col(ix, iy)] = (v>k)?1.0:v/k;
				}
			}
		};

		stream.set_n_act_stream(grid.ny);
		stream.set_grid(1, grid.ny);
		stream.exec(thr_hanning);

		return f;
	}

	// get two dimensional Hanning Filter
	template<class TVector, class TGrid>
	TVector func_hanning_2d(Stream<e_host> &stream, TGrid &grid, Value_type<TVector> k, bool shift, 
	Value_type<TVector> x_0=0, Value_type<TVector> x_e=0, Value_type<TVector> y_0=0, Value_type<TVector> y_e=0)
	{	
		using T = Value_type<TVector>;

		const T lx = grid.lx - (x_0+x_e);
		const T ly = grid.ly - (y_0+y_e);
		const T cx = c_2Pi/lx;
		const T cy = c_2Pi/ly;
		x_e = lx + x_0;
		y_e = ly + y_0;

		TVector fx;
		fx.reserve(grid.nx);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			auto Rx = (shift)?grid.Rx_shift(ix):grid.Rx(ix); 
			T v = ((Rx<x_0)||(Rx>x_e))?0:0.5*(1.0-cos(cx*(Rx-x_0)));
			fx.push_back(v);
		}

		TVector fy;
		fy.reserve(grid.ny);

		for(auto iy = 0; iy < grid.ny; iy++)
		{
			auto Ry = (shift)?grid.Ry_shift(iy):grid.Ry(iy); 
			T v = ((Ry<y_0)||(Ry>y_e))?0:0.5*(1.0-cos(cy*(Ry-y_0)));
			fy.push_back(v);
		}
		TVector f(grid.nxy());

		auto thr_hanning = [&](const Range &range)
		{
			for(auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for(auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					auto v = fx[ix]*fy[iy];
					f[grid.ind_col(ix, iy)] = (v>k)?1.0:v/k;
				}
			}
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_hanning);

		return f;
	}

	/********************************************************************/
	// get two dimensional Gaussian by row
	template<class TVector>
	TVector func_gaussian_1d(Stream<e_host> &stream, int nx, Value_type<TVector> dx, Value_type<TVector> sigma, bool shift,
	Value_type<TVector> x_0=0, Value_type<TVector> x_e=0)
	{	
		using T = Value_type<TVector>;
		int nxh = nx/2;

		const T lx = dx*nx - (x_0 + x_e);
		auto alpha_x = 0.5/(sigma*sigma);
		x_e = x_0 + lx;
		const T x_c = x_0 + 0.5*lx;

		TVector f(nx);
		auto thr_gaussian = [&](const Range &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				auto Rx = RSC(ixy, nxh, shift)*dx; 
				T v = ((Rx<x_0)||(Rx>x_e))?0:exp(-alpha_x*pow(Rx-x_c, 2));
				f[ixy] = v;
			}
		};

		stream.set_n_act_stream(nx);
		stream.set_grid(nx, 1);
		stream.exec(thr_gaussian);

		return f;
	}
	
	// get two dimensional Gaussian Filter
	template<class TVector, class TGrid>
	TVector func_gaussian_2d_by_row(Stream<e_host> &stream, TGrid &grid, Value_type<TVector> sigma, bool shift, 
	Value_type<TVector> x_0=0, Value_type<TVector> x_e=0)
	{	
		using T = Value_type<TVector>;

		const T lx = grid.lx - (x_0 + x_e);
		auto alpha_x = 0.5/(sigma*sigma);
		x_e = x_0 + lx;
		const T x_c = x_0 + 0.5*lx;

		TVector fx;
		fx.reserve(grid.nx);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			auto Rx = (shift)?grid.Rx_shift(ix):grid.Rx(ix); 
			T v = ((Rx<x_0)||(Rx>x_e))?0:exp(-alpha_x*pow(Rx-x_c, 2));
			fx.push_back(v);
		}

		TVector f(grid.nxy());

		auto thr_gaussian = [&](const Range &range)
		{
			for(auto iy = range.ixy_0; iy < range.ixy_e; iy++)
			{
				for(auto ix = 0; ix < grid.nx; ix++)
				{
					f[grid.ind_col(ix, iy)] = fx[ix];
				}
			}
		};

		stream.set_n_act_stream(grid.ny);
		stream.set_grid(1, grid.ny);
		stream.exec(thr_gaussian);

		return f;
	}

	// get two dimensional Gaussian Filter
	template<class TVector, class TGrid>
	TVector func_gaussian_2d(Stream<e_host> &stream, TGrid &grid, Value_type<TVector> sigma, bool shift, 
	Value_type<TVector> x_0=0, Value_type<TVector> x_e=0, Value_type<TVector> y_0=0, Value_type<TVector> y_e=0)
	{	
		using T = Value_type<TVector>;

		const T lx = grid.lx - (x_0 + x_e);
		const T ly = grid.ly - (y_0 + y_e);
		const T alpha_x = 0.5/(sigma*sigma);
		const T alpha_y = 0.5/(sigma*sigma);
		x_e = x_0 + lx;
		y_e = y_0 + ly;
		const T x_c = x_0 + 0.5*lx;
		const T y_c = y_0 + 0.5*ly;

		TVector fx;
		fx.reserve(grid.nx);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			auto Rx = (shift)?grid.Rx_shift(ix):grid.Rx(ix); 
			T v = ((Rx<x_0)||(Rx>x_e))?0:exp(-alpha_x*pow(Rx-x_c, 2));
			fx.push_back(v);
		}

		TVector fy;
		fy.reserve(grid.ny);

		for(auto iy = 0; iy < grid.ny; iy++)
		{
			auto Ry = (shift)?grid.Ry_shift(iy):grid.Ry(iy); 
			T v = ((Ry<y_0)||(Ry>y_e))?0:exp(-alpha_y*pow(Ry-y_c, 2));
			fy.push_back(v);
		}

		TVector f(grid.nxy());

		auto thr_gaussian = [&](const Range &range)
		{
			for(auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for(auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					f[grid.ind_col(ix, iy)] = fx[ix]*fy[iy];
				}
			}
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_gaussian);

		return f;
	}

	/********************************************************************/
	// get one dimensional Butterworth Filter
	template<class TVector>
	TVector func_butterworth_1d(Stream<e_host> &stream, int nx, Value_type<TVector> dx, Value_type<TVector> Radius, 
	int n, bool shift, Value_type<TVector> x_0=0, Value_type<TVector> x_e=0)
	{
		using T = Value_type<TVector>;
		
		int nxh = nx/2;

		const T lx = dx*nx - (x_0 + x_e);
		auto R02 = pow(Radius, 2);
		x_e = x_0 + lx;
		const T x_c = x_0 + 0.5*lx;

		TVector f(nx);
		auto thr_butterworth = [&](const Range &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				auto Rx = RSC(ixy, nxh, shift)*dx;
				T v = ((Rx<x_0)||(Rx>x_e))?0:1.0/(1.0+pow((Rx-x_c)*(Rx-x_c)/R02, n));
				f[ixy] = v;
			}
		};

		stream.set_n_act_stream(nx);
		stream.set_grid(nx, 1);
		stream.exec(thr_butterworth);

		return f;
	}

	// get two dimensional Butterworth Filter
	template<class TVector, class TGrid>
	TVector func_butterworth_2d_by_row(Stream<e_host> &stream, TGrid &grid, Value_type<TVector> Radius, 
	int n, bool shift, Value_type<TVector> x_0=0, Value_type<TVector> x_e=0)
	{
		using T = Value_type<TVector>;

		const T lx = grid.lx - (x_0 + x_e);
		const T R02 = pow(Radius, 2);
		x_e = x_0 + lx;
		const T x_c = x_0 + 0.5*lx;

		TVector fx;
		fx.reserve(grid.nx);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			auto Rx = (shift)?grid.Rx_shift(ix):grid.Rx(ix); 
			T v = ((Rx<x_0)||(Rx>x_e))?0:1.0/(1.0+pow((Rx-x_c)*(Rx-x_c)/R02, n));
			fx.push_back(v);
		}

		TVector f(grid.nxy());

		auto thr_butterworth = [&](const Range &range)
		{
			for(auto iy = range.ixy_0; iy < range.ixy_e; iy++)
			{
				for(auto ix = 0; ix < grid.nx; ix++)
				{
					f[grid.ind_col(ix, iy)] = fx[ix];
				}
			}
		};

		stream.set_n_act_stream(grid.ny);
		stream.set_grid(1, grid.ny);
		stream.exec(thr_butterworth);

		return f;
	}

	// get two dimensional Butterworth Filter
	template<class TVector, class TGrid>
	TVector func_butterworth_2d(Stream<e_host> &stream, TGrid &grid, Value_type<TVector> Radius, 
	int n, bool shift, Value_type<TVector> x_0=0, Value_type<TVector> x_e=0, Value_type<TVector> y_0=0, 
	Value_type<TVector> y_e=0)
	{
		using T = Value_type<TVector>;

		const T lx = grid.lx - (x_0 + x_e);
		const T ly = grid.ly - (y_0 + y_e);
		const T R02 = pow(Radius, 2);
		x_e = x_0 + lx;
		y_e = y_0 + ly;
		const T x_c = x_0 + 0.5*lx;
		const T y_c = y_0 + 0.5*ly;

		TVector f(grid.nxy());

		auto thr_butterworth = [&](const Range &range)
		{
			for(auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for(auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					auto Rx = (shift)?grid.Rx_shift(ix):grid.Rx(ix); 
					auto Ry = (shift)?grid.Ry_shift(iy):grid.Ry(iy); 
					T v = ((Rx<x_0)||(Rx>x_e)||(Ry<y_0)||(Ry>y_e))?0:1.0/(1.0+pow((pow(Rx-x_c, 2)+pow(Ry-y_c, 2))/R02, n));
					f[grid.ind_col(ix, iy)] = v;
				}
			}
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_butterworth);

		return f;
	}

	/********************************************************************/
	// get two dimensional radial distribution for regular grid
	void radial_distribution_2d(int nR, double *R, double *fR, int nRl, double *Rl, double *rl, double *frl, double *cfrl, bool reg, int typ)
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
				frl[i] /= ::fmax(1.0, cfrl[i]);
			}
	}

	template<class TGrid, class TVector>
	void radial_distribution_2d(TGrid &grid, TVector &Im, Value_type<TVector> x_i, 
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
					auto ir = static_cast<int>(floor(R/dR));
					frl[ir] += Im[grid.ind_col(ix, iy)];
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
	void getCumRadDist_2d(int ny, int nx, bool shift, double *fI, int nr, double *r, double *fIr)
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
	TVector smooth(TVector &f_i, Value_type<TVector> nkr)
	{
		using T = Value_type<TVector>;

		int nf_i = f_i.size();
		int nk0 = -nkr;
		int nke = nkr+1;

		TVector f_o(nf_i);

		for(auto ix = 0; ix < nf_i; ix++)
		{
			int it0 = max(ix+nk0, 0);
			int ite = min(ix+nke, nf_i);

			T f_mean = 0;
			for (auto it = it0; it < ite; it++)
			{
				f_mean += f_i[it];
			}

			f_o[ix] = f_mean/(ite-it0);
		}

		return f_o;
	}

	// mean smooth
	template<class TVector>
	TVector smooth(Stream<e_host> &stream, TVector &f_i, Value_type<TVector> nkr)
	{
		using T = Value_type<TVector>;

		int nf_i = f_i.size();
		int nk0 = -nkr;
		int nke = nkr+1;

		auto krn_mean = [&](const int &ix_i, TVector &f_i, TVector &f_o)
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

		TVector fv(nf_i);

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

		return fv;
	}

	/*******************************************************************/
	// add Poisson noise
	template <class TVector>
	TVector add_poisson_noise(Stream<e_host> &stream, TVector &M_i, 
	Value_type<TVector> scf)
	{	
		using T = Value_type<TVector>;

		TVector M(M_i.size());
		auto thr_poisson_noise = [&](const Range &range)
		{
			std::mt19937_64 gen;
			std::poisson_distribution<int> randp;

			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				auto x0 = scf*M_i[ixy];
				randp.param(std::poisson_distribution<int>::param_type(x0));
				M[ixy] = randp(gen);
			}
		};

		stream.set_n_act_stream(M_i.size());
		stream.set_grid(M_i.size(), 1);
		stream.exec(thr_poisson_noise);

		return M;
	};

	// add Poisson noise
	template <class TVector>
	TVector add_poisson_noise_by_SNR(Stream<e_host> &stream, TVector &Im_i, 
		Value_type<TVector> SNR_i, Value_type<TVector> &scl_o)
	{	
		using T = Value_type<TVector>;

		auto get_SNR = [](Stream<e_host> &stream, TVector &Im, T Im_std, T scf)->T
		{
			T x_mean = 0;
			T x_var = 0;
			auto thr_SNR = [&](const Range &range)
			{
				std::mt19937_64 gen;
				std::poisson_distribution<int> rand;

				T x_mean_partial = 0;
				T x_var_partial = 0;
				for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
				{
					auto x0 = scf*Im[ixy];
					rand.param(std::poisson_distribution<int>::param_type(x0));
					auto xn = rand(gen)-x0;
					x_mean_partial += xn;
					x_var_partial += xn*xn;
				}

				stream.stream_mutex.lock();
				x_mean += x_mean_partial;
				x_var += x_var_partial;
				stream.stream_mutex.unlock();
			};

			stream.set_n_act_stream(Im.size());
			stream.set_grid(1, Im.size());
			stream.exec(thr_SNR);

			x_mean /= Im.size();
			x_var = x_var/Im.size()-x_mean*x_mean;

			return scf*Im_std/sqrt(x_var);
		};

		auto Im_i_std = sqrt(variance(stream, Im_i));

		T SNR_k = get_SNR(stream, Im_i, Im_i_std, 1);
	
		T k_0 = 1;
		T k_e = 1;

		if(SNR_k<SNR_i)
		{
			do
			{
				k_e *= 2;
				SNR_k = get_SNR(stream, Im_i, Im_i_std, k_e);
			}
			while (SNR_k < SNR_i);
			k_0 = k_e/2;	
		}
		else
		{
			do
			{
				k_0 /= 2;
				SNR_k = get_SNR(stream, Im_i, Im_i_std, k_0);
			}
			while (SNR_k >= SNR_i);
			k_e = 2*k_0;
		}


		// bisection method
		int ic = 0;
		do
		{
			scl_o = 0.5*(k_0 + k_e);
			auto SNR_k = get_SNR(stream, Im_i, Im_i_std, scl_o);

			if(SNR_k < SNR_i)
			{
				k_0 = scl_o;
			}
			else
			{
				k_e = scl_o;
			}
			ic++;
		} while ((fabs(SNR_i-SNR_k)>0.05) && (ic<10));

		// add Poisson noise
		return add_poisson_noise(stream, Im_i, scl_o);
	}

	/*******************************************************************/
	// add periodic boundary border to the image
	template <class TGrid, class TVector_i, class TVector_o>
	vector<int> add_PB_border(TGrid &grid_i, TVector_i &image_i, int border_x, int border_y, TGrid &grid_o, TVector_o &image_o)
	{
		Prime_Num CP;

		auto ix0 = border_x;
		auto nx = CP(grid_i.nx+2*border_x, eDST_Greater_Than);
		auto ixe = ix0 + grid_i.nx;
		auto bx_l = border_x;
		auto bx_r = nx - ixe;

		auto iy0 = border_y;
		auto ny = CP(grid_i.ny+2*border_y, eDST_Greater_Than);
		auto iye = iy0 + grid_i.ny;
		auto by_t = border_y;
		auto by_b = ny - iye;

		grid_o.set_input_data(nx, ny, nx*grid_i.dRx, ny*grid_i.dRy, grid_i.dz, grid_i.bwl, grid_i.pbc_xy);
		
		image_o.resize(grid_o.nxy());

		// copy central image
		for(auto ix_o=ix0; ix_o<ixe; ix_o++)
		{
			for(auto iy_o=iy0; iy_o<iye; iy_o++)
			{
				auto ix_i = ix_o-ix0;
				auto iy_i = iy_o-iy0;
				auto ixy_i = grid_i.ind_col(ix_i, iy_i);
				auto ixy_o = grid_o.ind_col(ix_o, iy_o);
				image_o[ixy_o] = image_i[ixy_i] ;
			}
		}

		// left
		for(auto ix_o=0; ix_o<bx_l; ix_o++)
		{
			for(auto iy_o=iy0; iy_o<iye; iy_o++)
			{
				auto ix_i = 2*bx_l-ix_o;
				auto iy_i = iy_o;
				image_o[grid_o.ind_col(ix_o, iy_o)] = image_o[grid_o.ind_col(ix_i, iy_i)] ;
			}
		}

		// right
		for(auto ix_o=ixe; ix_o<grid_o.nx; ix_o++)
		{
			for(auto iy_o=iy0; iy_o<iye; iy_o++)
			{
				auto ix_i = ixe-2-(ix_o-ixe);
				auto iy_i = iy_o;
				image_o[grid_o.ind_col(ix_o, iy_o)] = image_o[grid_o.ind_col(ix_i, iy_i)] ;
			}
		}

		// top
		for(auto ix_o=0; ix_o<grid_o.nx; ix_o++)
		{
			for(auto iy_o=0; iy_o<by_t; iy_o++)
			{
				auto ix_i = ix_o;
				auto iy_i = 2*by_t-iy_o;
				image_o[grid_o.ind_col(ix_o, iy_o)] = image_o[grid_o.ind_col(ix_i, iy_i)] ;
			}
		}

		// bottom
		for(auto ix_o=0; ix_o<grid_o.nx; ix_o++)
		{
			for(auto iy_o=iye; iy_o<grid_o.ny; iy_o++)
			{
				auto ix_i = ix_o;
				auto iy_i = iye-2-(iy_o-iye);
				image_o[grid_o.ind_col(ix_o, iy_o)] = image_o[grid_o.ind_col(ix_i, iy_i)] ;
			}
		}

		vector<int> points = {ix0, ixe, iy0, iye};

		return points;
	}

	/*******************************************************************/
	// add periodic boundary border to the image
	template <class TGrid, class TVector>
	void set_const_border(TGrid &grid, TVector &image, Value_type<TVector> xb_0, Value_type<TVector> xb_e, 
	Value_type<TVector> yb_0, Value_type<TVector> yb_e)
	{
		using T = Value_type<TVector>;

		int ix_0 = grid.ceil_dRx(xb_0);
		int ix_e = grid.nx-grid.ceil_dRx(xb_e);

		int iy_0 = grid.ceil_dRy(yb_0);
		int iy_e = grid.ny-grid.ceil_dRy(yb_e);

		//get average value
		int ix_i0 = ix_0 + max(10, (ix_e-ix_0)/10);
		int ix_ie = ix_e-ix_i0;
		int iy_i0 = iy_0 + max(10, (iy_e-iy_0)/10);
		int iy_ie = iy_e-iy_i0;

		T val = 0;
		T cval = 0;
		for(auto ix=ix_0; ix<ix_e; ix++)
		{
			for(auto iy=iy_0; iy<iy_e; iy++)
			{
				if((ix<=ix_i0)||(ix>=ix_ie)||(iy<=iy_i0)||(iy>=iy_ie))
				{
					val += image[grid.ind_col(ix, iy)];
					cval++;
				}
			}
		}
		val /= cval;

		// left
		for(auto ix=0; ix<ix_0; ix++)
		{
			for(auto iy=0; iy<grid.ny; iy++)
			{
				image[grid.ind_col(ix, iy)] = val;
			}
		}

		// right
		for(auto ix=ix_e; ix<grid.nx; ix++)
		{
			for(auto iy=0; iy<grid.ny; iy++)
			{
				image[grid.ind_col(ix, iy)] = val;
			}
		}

		// top
		for(auto ix=0; ix<grid.nx; ix++)
		{
			for(auto iy=0; iy<iy_0; iy++)
			{
				image[grid.ind_col(ix, iy)] = val;
			}
		}

		// bottom
		for(auto ix=0; ix<grid.nx; ix++)
		{
			for(auto iy=iy_e; iy<grid.ny; iy++)
			{
				image[grid.ind_col(ix, iy)] = val;
			}
		}
	}

	/*******************************************************************/	
	// extract region ix0<=x<ixe
	template<class TVector_o, class TVector_i>
	TVector_o extract_region_real_part(Stream<e_host> &stream, int nx_src, int ny_src, 
	TVector_i &Im_src, int ix0_src, int ixe_src, int iy0_src, int iye_src)
	{
		auto nx_dst = ixe_src - ix0_src;
		auto ny_dst = iye_src - iy0_src;

		TVector_o Im_dst(nx_dst*ny_dst);

		auto thr_extract_region = [&](const Range &range)
		{
			for (auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for (auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					Im_dst[ix*ny_dst+iy] = Im_src[(ix0_src+ix)*ny_src+(iy0_src+iy)].real();
				}
			}
		};

		stream.set_n_act_stream(nx_dst);
		stream.set_grid(nx_dst, ny_dst);
		stream.exec(thr_extract_region);

		return Im_dst;
	}

	// extract real part of vector
	template<class TVector_o, class TVector_i>
	TVector_o extract_real_part(Stream<e_host> &stream, TVector_i &Im_src)
	{
		TVector_o Im_dst(Im_src.size());

		auto thr_extract_real = [](const Range &range, TVector_i &Im_src, TVector_o &Im_dst)
		{
			for (auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				Im_dst[ixy] = Im_src[ixy].real();
			}
		};

		stream.set_n_act_stream(Im_src.size());
		stream.set_grid(Im_src.size(), 1);
		stream.exec(thr_extract_real, Im_src, Im_dst);

		return Im_dst;
	}

	// extract region ix0<=x<ixe
	template<class TVector_o, class TVector_i>
	TVector_o extract_region_abs(Stream<e_host> &stream, int nx_src, int ny_src, 
	TVector_i &Im_src, int ix0_src, int ixe_src, int iy0_src, int iye_src)
	{
		auto nx_dst = ixe_src - ix0_src;
		auto ny_dst = iye_src - iy0_src;

		TVector_o Im_dst(nx_dst*ny_dst);

		auto thr_extract_region = [&](const Range &range)
		{
			for (auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for (auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					Im_dst[ix*ny_dst+iy] = abs(Im_src[(ix0_src+ix)*ny_src+(iy0_src+iy)]);
				}
			}
		};

		stream.set_n_act_stream(nx_dst);
		stream.set_grid(nx_dst, ny_dst);
		stream.exec(thr_extract_region);

		return Im_dst;
	}

	// extract abs of vector
	template<class TVector_o, class TVector_i>
	TVector_o extract_abs(Stream<e_host> &stream, TVector_i &Im_src)
	{
		TVector_o Im_dst(Im_src.size());

		auto thr_extract_abs = [](const Range &range, TVector_i &Im_src, TVector_o &Im_dst)
		{
			for (auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				Im_dst[ixy] = abs(Im_src[ixy]);
			}
		};

		stream.set_n_act_stream(Im_src.size());
		stream.set_grid(Im_src.size(), 1);
		stream.exec(thr_extract_abs, Im_src, Im_dst);

		return Im_dst;
	}
	/*******************************************************************/
	// analytic convolution
	template <class TGrid, class TFn, class TVector>
	void analytic_conv(Stream<e_host> &stream, FFT2<Value_type<TGrid>, e_host> &fft2,
		TGrid &grid, TFn fn, TVector &image)
	{
		fft2_shift(stream, grid, image);
		fft2.forward(image);

		auto thr_analytic_convolution = [&](const Range &range)
		{
			for (auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for (auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					image[grid.ind_col(ix, iy)] *= fn(ix, iy, grid);
				}
			}
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_analytic_convolution);

		fft2.inverse(image);
		fft2_shift(stream, grid, image);
	}

	/*******************************************************************/
	// Gaussian convolution
	template <class TGrid, class TVector>
	TVector gaussian_conv(Stream<e_host> &stream, FFT2<Value_type<TGrid>, e_host> &fft2, 
	TGrid &grid_i, Value_type<TGrid> sigma_r, TVector &image_i)
	{	
		using T = Value_type<TVector>;
		using TVector_c = vector<complex<T>>;

		T alpha = 2.0*c_Pi2*sigma_r*sigma_r;
		auto sigma_r_lim = 4*sigma_r;
		auto border_x = grid_i.nx_dRx(sigma_r_lim);
		auto border_y = grid_i.ny_dRy(sigma_r_lim);

		TGrid grid;
		TVector_c image;
		// add borders
		auto pts = add_PB_border(grid_i, image_i, border_x, border_y, grid, image);

		auto krn_gaussian = [=](const int &ix, const int &iy, const TGrid &grid)
		{
			return exp(-alpha*grid.g2_shift(ix, iy))*grid.inxy;
		};

		// create 2d plan
		fft2.create_plan_2d(grid.ny, grid.nx, stream.size());

		// perform analytic convolution
		analytic_conv(stream, fft2, grid, krn_gaussian, image);

		int ix0 = pts[0], ixe = pts[1];
		int iy0 = pts[2], iye = pts[3];

		return extract_region_real_part<TVector>(stream, grid.nx, grid.ny, image, ix0, ixe, iy0, iye);
	}

	// Gaussian deconvolution
	template <class TGrid, class TVector>
	TVector gaussian_deconv(Stream<e_host> &stream, FFT2<Value_type<TGrid>, e_host> &fft2, 
	TGrid &grid_i, Value_type<TGrid> sigma_r, Value_type<TGrid> PSNR, TVector &image_i)
	{	
		using T = Value_type<TVector>;
		using TVector_c = vector<complex<T>>;

		// deconvolution
		T alpha = 2.0*c_Pi2*sigma_r*sigma_r;
		auto sigma_r_lim = 4.0*sigma_r;
		auto border_x = grid_i.nx_dRx(sigma_r_lim);
		auto border_y = grid_i.ny_dRy(sigma_r_lim);

		TGrid grid;
		TVector_c image;
		// add borders
		auto pts = add_PB_border(grid_i, image_i, border_x, border_y, grid, image);

		// sigma_r ---> sigma_f = 1/(2*pi*sigma_r)
		T sigma_f = 3.0/(c_2Pi*sigma_r);
		T g2_lim = pow(min(sigma_f, grid.g_max()), 2);
		auto krn_dc_gaussian = [=](const int &ix, const int &iy, const TGrid &grid)
		{
			auto g2 = grid.g2_shift(ix, iy);
			if(g2<g2_lim)
			{
				auto fg = exp(-alpha*g2);
				return fg/((fg*fg+PSNR)*grid.nxy());
			}
			else
			{
				return T(0);
			}
		};

		// create 2d plan
		fft2.create_plan_2d(grid.ny, grid.nx, stream.size());

		// perform analytic convolution
		analytic_conv(stream, fft2, grid, krn_dc_gaussian, image);

		int ix0 = pts[0], ixe = pts[1];
		int iy0 = pts[2], iye = pts[3];

		return extract_region_real_part<TVector>(stream, grid.nx, grid.ny, image, ix0, ixe, iy0, iye);
	}

	// Modified Gaussian deconvolution
	template <class TGrid, class TVector>
	TVector mod_gaussian_deconv(Stream<e_host> &stream, FFT2<Value_type<TGrid>, e_host> &fft2, 
	TGrid &grid_i, Value_type<TGrid> sigma_r, Value_type<TGrid> PSNR, TVector &image_i)
	{	
		using T = Value_type<TVector>;
		using TVector_c = vector<complex<T>>;

		// deconvolution
		T sigma_a = sigma_r;
		T sigma_b = 0.25*sigma_r;
		T alpha_a = 2.0*c_Pi2*pow(sigma_a, 2);
		T alpha_b = 2.0*c_Pi2*pow(sigma_b, 2);

		auto sigma_r_lim = 4.0*sigma_a;
		auto border_x = grid_i.nx_dRx(sigma_r_lim);
		auto border_y = grid_i.ny_dRy(sigma_r_lim);

		TGrid grid;
		TVector_c image;
		// add borders
		vector<int> pts;
		
		pts = add_PB_border(grid_i, image_i, border_x, border_y, grid, image);

		// sigma_r ---> sigma_f = 1/(2*pi*sigma_r)
		T sigma_f_lim = 3.0*(1.0/(c_2Pi*sigma_r));
		T g2_lim = pow(min(sigma_f_lim, grid.g_max()), 2);
		auto krn_dc_gaussian = [=](const int &ix, const int &iy, const TGrid &grid)
		{
			auto g2 = grid.g2_shift(ix, iy);
			if(g2<g2_lim)
			{
				auto fg_a = exp(-alpha_a*g2);
				auto fg_b = exp(-alpha_b*g2);
				return fg_a*fg_b/((fg_a*fg_a+PSNR)*grid.nxy());
			}
			else
			{
				return T(0);
			}
		};

		// create 2d plan
		fft2.create_plan_2d(grid.ny, grid.nx, stream.size());

		// perform analytic convolution
		analytic_conv(stream, fft2, grid, krn_dc_gaussian, image);

		int ix0 = pts[0], ixe = pts[1];
		int iy0 = pts[2], iye = pts[3];

		return extract_region_real_part<TVector>(stream, grid.nx, grid.ny, image, ix0, ixe, iy0, iye);

	}

	/*******************************************************************/
	// get weighted position
	template<class TGrid, class TVector>
	r2d<Value_type<TVector>> Rx_Ry_weight(TGrid &grid, TVector &Im, r2d<Value_type<TVector>> p_i, Value_type<TVector> radius_i, bool env)
	{
		using T = Value_type<TVector>;

		T dR = grid.dRx;
		int nrl = static_cast<int>(floor(radius_i/dR+0.5));
		T R_max = dR*nrl;

		auto ix_i = static_cast<int>(floor(p_i.x/dR));
		int ix0 = max(ix_i-nrl, 0);
		int ixe = min(ix_i+nrl+1, grid.nx);

		auto iy_i = static_cast<int>(floor(p_i.y/dR));
		int iy0 = max(iy_i-nrl, 0);
		int iye = min(iy_i+nrl+1, grid.ny);

		T R2_max = pow(R_max, 2);

		r2d<T> p(0, 0);

		T alpha = 0.5/pow(R_max, 2);
		T wt = 0;
		for(auto ix = ix0; ix < ixe; ix++)
		{
			for(auto iy = iy0; iy < iye; iy++)
			{
				auto R2 = grid.R2(ix, iy, p_i.x, p_i.y);
				if(R2 < R2_max)
				{
					int ixy = grid.ind_col(ix, iy);
					auto imv = ((env)?exp(-alpha*R2):1)*Im[ixy];
					p += imv*r2d<T>(grid.Rx(ix), grid.Ry(iy)); 
					wt += imv;
				}
			}
		}
		p /= wt; 

		if(module(p-p_i)>=radius_i)
		{
			p = p_i;
		}

		return p;
	}
	
	// get fit position
	template<class TGrid, class TVector>	
	TVector Rx_Ry_fit(TGrid &grid, TVector &Im_i, r2d<Value_type<TVector>> p_i, 
	Value_type<TVector> sigma_i, Value_type<TVector> radius_i)
	{
		using T = Value_type<TVector>;

		auto select_cir_reg = [](TGrid &grid, TVector &Im, r2d<T> p, 
		T radius, TVector &Rx, TVector &Ry, TVector &Ixy, T &Rx_sf, T &Ry_sf, T &Rxy_sc, T &Ixy_sc)
		{
			T R_max = radius;
			T R2_max = pow(R_max, 2);

			auto range = grid.index_range(p, R_max);

			Rx.clear();
			Rx.reserve(range.ixy_e);
			Ry.clear();
			Ry.reserve(range.ixy_e);
			Ixy.clear();
			Ixy.reserve(range.ixy_e);

			// select circular region
			for (auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for (auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					if (grid.R2(ix, iy, p.x, p.y) < R2_max)
					{
						Rx.push_back(grid.Rx(ix));
						Ry.push_back(grid.Ry(iy));
						Ixy.push_back(Im[grid.ind_col(ix, iy)]);
					}
				}
			}

			Rx_sf = p.x;
			Ry_sf = p.y;
			Rxy_sc = R_max;

			Ixy_sc = max_data(Ixy);
			int m = Ixy.size();
			for (auto ixy = 0; ixy < m; ixy++)
			{
				Rx[ixy] = (Rx[ixy]-Rx_sf)/Rxy_sc;
				Ry[ixy] = (Ry[ixy]-Ry_sf)/Rxy_sc;
				Ixy[ixy] = Ixy[ixy]/Ixy_sc;
			}

			Rx.shrink_to_fit();
			Ry.shrink_to_fit();
			Ixy.shrink_to_fit();
		};

		TVector Rx, Ry, Ixy;
		T Rx_sf, Ry_sf, Rxy_sc, Ixy_sc;

		select_cir_reg(grid, Im_i, p_i, radius_i, Rx, Ry, Ixy, Rx_sf, Ry_sf, Rxy_sc, Ixy_sc);

		T sigma = sigma_i/Rxy_sc;

		int m = Ixy.size();
		int n = 6;

		TVector J(m*n);
		TVector d_Ixy(m);

		auto get_chi2 = [](TVector &Rx, TVector &Ry, TVector &Ixy, TVector &coef)->T
		{
			T x_0 = coef[0];
			T y_0 = coef[1];

			T theta = coef[5];

			T A = coef[2];

			T cos_1t = cos(theta);
			T sin_1t = sin(theta);
			T cos_2t = cos(2*theta);
			T sin_2t = sin(2*theta);

			T sx2 = pow(coef[3], 2);
			T sy2 = pow(coef[4], 2);

			T a = 0.5*(cos_1t*cos_1t/sx2 + sin_1t*sin_1t/sy2);
			T b = 0.5*(sin_1t*sin_1t/sx2 + cos_1t*cos_1t/sy2);
			T c = 0.5*(-sin_2t/sx2 + sin_2t/sy2);

			T chi2 = 0;
			T chi2_ee = 0;
			const int m = Ixy.size();
			for (auto im = 0; im < m; im++)
			{
				T x = Rx[im]-x_0;
				T y = Ry[im]-y_0;
				T v = Ixy[im]-A*exp(-a*x*x-b*y*y-c*x*y);
				host_device_detail::kh_sum(chi2, v*v, chi2_ee);
			}
			return chi2;
		};

		auto get_dIxy_J = [](TVector &Rx, TVector &Ry, TVector &Ixy, TVector &coef, 
		TVector &dIxy, TVector &J)
		{
			T x_0 = coef[0];
			T y_0 = coef[1];

			T theta = coef[5];

			T A = coef[2];

			T cos_1t = cos(theta);
			T sin_1t = sin(theta);
			T cos_2t = cos(2*theta);
			T sin_2t = sin(2*theta);

			T sx2 = pow(coef[3],2);
			T sy2 = pow(coef[4],2);

			T sx3 = pow(coef[3],3);
			T sy3 = pow(coef[4],3);

			T a = 0.5*(cos_1t*cos_1t/sx2 + sin_1t*sin_1t/sy2);
			T b = 0.5*(sin_1t*sin_1t/sx2 + cos_1t*cos_1t/sy2);
			T c = 0.5*(-sin_2t/sx2 + sin_2t/sy2);

			const int m = Ixy.size();
			for (auto im = 0; im < m; im++)
			{
				T x = Rx[im]-x_0;
				T y = Ry[im]-y_0;

				T v = exp(-a*x*x-b*y*y-c*x*y);
				T f = A*v;

				J[0*m+im] = (2*a*x+c*y)*f;
				J[1*m+im] = (2*b*y+c*x)*f;
				J[2*m+im] = v;    
				J[3*m+im] = (pow(cos_1t*x,2)+pow(sin_1t*y,2)-sin_2t*x*y)*f/sx3;
				J[4*m+im] = (pow(sin_1t*x,2)+pow(cos_1t*y,2)+sin_2t*x*y)*f/sy3;
				J[5*m+im] = (sin_2t*x*x-sin_2t*y*y+2*cos_2t*x*y)*f*(0.5/sx2-0.5/sy2);

				dIxy[im] = Ixy[im] - f;
			}
		};

		TVector coef_0 = {0, 0, 1, sigma, sigma, 0};
		TVector coef_min = {-sigma/5, -sigma/5, 0.75, sigma/3, sigma/3, -T(c_Pi)};
		TVector coef_max = {sigma/5, sigma/5, 1.25, 3*sigma, 3*sigma, T(c_Pi)};
		TVector coef = coef_0;

		T chi2 = get_chi2(Rx, Ry, Ixy, coef);

		T lambda = 2;
		T lambda_f = 2;

		lapack::FLSF<T> flsf;

		const int niter = 100;
		for (auto iter = 0; iter<niter; iter++)
		{
			get_dIxy_J(Rx, Ry, Ixy, coef, d_Ixy, J);

			TVector d_coef(n);
			TVector D(n);
			TVector G(n);
			flsf(m, n, J.data(), d_Ixy.data(), d_coef.data(), lambda, D.data(), G.data());

			TVector coef_t = coef;
			T rho_f = 0;
			T G_max = 0;
			for(auto ic=0; ic<n; ic++)
			{
				coef_t[ic] += d_coef[ic];
				rho_f += coef_t[ic]*(D[ic]*coef_t[ic]+G[ic]);
				G_max = ::fmax(G_max, abs(G[ic]));
			}

			T chi2_t = get_chi2(Rx, Ry, Ixy, coef_t);
			T rho = (chi2-chi2_t)/rho_f;

			if ((G_max<1e-6)||abs(chi2-chi2_t)<1e-7)
			{
				break;
			}

			if(rho>0)
			{
				coef = coef_t;

				for(auto ic=0; ic<n; ic++)
				{
					coef[ic] = min(max(coef[ic], coef_min[ic]), coef_max[ic]);
				}

				chi2 = get_chi2(Rx, Ry, Ixy, coef);

				lambda = (rho>1e-6)?(::fmax(lambda/lambda_f, 1e-7)):lambda;
			}
			else
			{
				lambda = ::fmin(lambda*lambda_f, 1e+7);
			}
		}

		coef[0] = coef[0]*Rxy_sc + Rx_sf;
		coef[1] = coef[1]*Rxy_sc + Ry_sf;
		coef[2] = coef[2]*Ixy_sc;
		coef[3] = coef[3]*Rxy_sc;
		coef[4] = coef[4]*Rxy_sc;

		return coef;
	}

	template<class TGrid, class TVector>
	Value_type<TVector> neighbor_radius(TGrid &grid_i, TVector &Im_i, r2d<Value_type<TVector>> p_i, Value_type<TVector> radius_i)
	{
		using T = Value_type<TVector>;

		T R2_max = pow(radius_i, 2);

		int ix0 = grid_i.lb_index_x(p_i.x - radius_i);
		int ixe = grid_i.ub_index_x(p_i.x + radius_i);

		int iy0 = grid_i.lb_index_y(p_i.y - radius_i);
		int iye = grid_i.ub_index_y(p_i.y + radius_i);

		int nrl = grid_i.ceil_dR_min(radius_i);

		// get radial distribution
		TVector frl(nrl);
		TVector cfrl(nrl);

		T dR = grid_i.dR_min();

		for(auto ix = ix0; ix < ixe; ix++)
		{
			for(auto iy = iy0; iy < iye; iy++)
			{
				T R2_d = grid_i.R2(ix, iy, p_i.x, p_i.y);
				if(R2_d < R2_max)
				{
					auto ir = static_cast<int>(floor(sqrt(R2_d)/dR));
					frl[ir] += Im_i[grid_i.ind_col(ix, iy)];
					cfrl[ir] += 1.0;
				}
			}
		}

		for(auto ir=0; ir<frl.size(); ir++)
		{
			frl[ir] /= ::fmax(1.0, cfrl[ir]);
		}

		frl[0] = ::fmax(frl[0], 1.01*frl[1]);

		frl = smooth(frl, 1);

		// derivative and minimun
		int ir0 = frl.size()-1;
		for(auto ir=0; ir<frl.size()-1; ir++)
		{
			if(frl[ir]<frl[ir+1])
			{
				ir0 = ir+1;
				break;
			}
		}

		int irm = frl.size()-1;
		for(auto ir=ir0; ir<frl.size()-1; ir++)
		{
			if(frl[ir]>frl[ir+1])
			{
				irm = ir;
				break;
			}
		}

		return max(2, irm)*dR;
	}

	// select circular region
	template<class TGrid, class TVector>
	Value_type<TVector> mean_cir_reg(TGrid &grid_i, TVector &Im_i, r2d<Value_type<TVector>> p_i, 
	Value_type<TVector> radius_i, Value_type<TVector> bg_i)
	{
		using T = Value_type<TVector>;

		T R2_max = radius_i*radius_i;

		auto range = grid_i.index_range(p_i, radius_i);

		// select circular region
		T I_mean = 0;
		int Ic = 0;
		for (auto ix = range.ix_0; ix < range.ix_e; ix++)
		{
			for (auto iy = range.iy_0; iy < range.iy_e; iy++)
			{
				T R2_d = grid_i.R2(ix, iy, p_i.x, p_i.y);
				if (R2_d < R2_max)
				{
					I_mean += Im_i[grid_i.ind_col(ix, iy)]-bg_i;
					Ic++;
				}
			}
		}
		I_mean /= Ic;
		return I_mean;
	}

	/*******************************************************************/
	template<class TGrid, class TVector>
	TVector interp_profile(TGrid &grid, TVector &Im, const r2d<Value_type<TVector>> &p1, const r2d<Value_type<TVector>> &p2, int nr)
	{
		TVector v;

		if(!(grid.ckb_bound(p1) && grid.ckb_bound(p2)))
		{
			return v;
		}

		if(module(p1-p2)<grid.dR_min())
		{
			return v;
		}

		using T = Value_type<TGrid>;

		auto interp_2d = [](const r2d<T> &p, TGrid &grid, TVector &Im)->T
		{
			auto ix = grid.lb_index_x(p.x);
			auto iy = grid.lb_index_y(p.y);

			T f11 = Im[grid.ind_col(ix, iy)];
			T f12 = Im[grid.ind_col(ix, iy+1)];
			T f21 = Im[grid.ind_col(ix+1, iy)];
			T f22 = Im[grid.ind_col(ix+1, iy+1)];

			T x1 = grid.Rx(ix);
			T x2 = grid.Rx(ix+1);
			T y1 = grid.Ry(iy);
			T y2 = grid.Ry(iy+1);
					
			T dx1 = p.x-x1;
			T dx2 = x2-p.x;
			T dy1 = p.y-y1;
			T dy2 = y2-p.y;

			T f = (dx2*(f11*dy2 + f12*dy1)+dx1*(f21*dy2 + f22*dy1))/((x2-x1)*(y2-y1));
			return f;

		};

		r2d<T> p12 = p2-p1;
		T mp12 = p12.module();

		nr = (nr<=0)?static_cast<int>(ceil(mp12/grid.dR_min())):nr;
		nr = max(nr, 2);
		T dr = mp12/(nr-1);
		r2d<T> u = dr*normalized(p12);

		v.reserve(nr);
		for(auto ir=0; ir<nr; ir++)
		{
			r2d<T> p = p1 + T(ir)*u;
			v.push_back(interp_2d(p, grid, Im));
		}
		return v;
	}
	
	template<class TGrid, class TVector>
	void interp_profile(TGrid &grid, TVector &Im, Value_type<TVector> bg, r2d<Value_type<TVector>> p1, r2d<Value_type<TVector>> p2, 
	TVector &x, TVector &y)
	{
		x.clear();
		y.clear();

		if(!(grid.ckb_bound(p1) && grid.ckb_bound(p2)))
		{
			return;
		}

		if(module(p1-p2)<grid.dR_min())
		{
			return;
		}

		using T = Value_type<TGrid>;

		auto interp_2d = [](const r2d<T> &p, TGrid &grid, TVector &Im)->T
		{
			auto ix = grid.lb_index_x(p.x);
			auto iy = grid.lb_index_y(p.y);

			T f11 = Im[grid.ind_col(ix, iy)];
			T f12 = Im[grid.ind_col(ix, iy+1)];
			T f21 = Im[grid.ind_col(ix+1, iy)];
			T f22 = Im[grid.ind_col(ix+1, iy+1)];

			T x1 = grid.Rx(ix);
			T x2 = grid.Rx(ix+1);
			T y1 = grid.Ry(iy);
			T y2 = grid.Ry(iy+1);
					
			T dx1 = p.x-x1;
			T dx2 = x2-p.x;
			T dy1 = p.y-y1;
			T dy2 = y2-p.y;

			T f = (dx2*(f11*dy2 + f12*dy1)+dx1*(f21*dy2 + f22*dy1))/((x2-x1)*(y2-y1));
			return f;

		};

		int Ixy_1 = Im[grid.ixy(p1.x, p1.y)];
		int Ixy_2 = Im[grid.ixy(p2.x, p2.y)];

		r2d<T> p12 = p2-p1;
		T mp12 = p12.module();
		int nr = grid.ceil_dR_min(mp12);
		T dr = mp12/(nr-1);
		r2d<T> u = dr*normalized(p12);

		x.reserve(nr);
		y.reserve(nr);
		for(auto ir=0; ir<nr; ir++)
		{
			x.push_back(ir*dr);
			r2d<T> p = p1 + T(ir)*u;
			T v = interp_2d(p, grid, Im);
			y.push_back(v-bg);
		}
	}

	// get fit position
	template<class TVector>	
	TVector find_peaks_vector_typ_1(TVector &x, TVector &y, Value_type<TVector> y_thr)
	{
		using T = Value_type<TVector>;
		TVector x_peak(y.size());

		T x_max = 0;
		T y_max = y_thr;
		bool bb = y[0]>y_thr;
		int ipeak = 0;
		for(auto ix=0; ix<y.size(); ix++)
		{
			if(y[ix]>y_thr)
			{
				if(y_max < y[ix])
				{
					y_max = y[ix];
					x_max = x[ix];
				}
				bb = true;
			}
			else if(bb)
			{
				x_peak[ipeak++] = x_max;
				y_max = y_thr;	
				bb = false;
			}
		}

		if(bb)
		{
			x_peak[ipeak++] = x_max;
		}

		x_peak.resize(ipeak);
		x_peak.shrink_to_fit();

		return x_peak;
	}

	/*******************************************************************/
	// find and fit peak position
	template<class TGrid, class TVector>
	r2d<Value_type<TVector>> max_pos(TGrid &grid, TVector &Im, 
	r2d<Value_type<TVector>> p_i, Value_type<TVector> sigma_i, Value_type<TVector> radius_i)
	{
		using T = Value_type<TVector>;

		auto coef = Rx_Ry_fit(grid, Im, p_i, sigma_i, radius_i);
		return r2d<T>(coef[0], coef[1]);
	}

	// shift 2d
	template <class TGrid, class TVector>
	enable_if_complex_host_vector<TVector, void>
	shift_2d(Stream<e_host> &stream, FFT2<Value_type<TGrid>, e_host> &fft2, TGrid &grid, 
	r2d<Value_type<TGrid>> p, TVector &Im)
	{	
		fft2.create_plan_2d(grid.ny, grid.nx, stream.size());

		fft2_shift(stream, grid, Im);
		fft2.forward(Im);
		phase_factor_2d(stream, grid, -c_2Pi*p.x, -c_2Pi*p.y, Im, Im);
		fft2.inverse(Im);
		fft2_shift(stream, grid, Im);
	}

	// shift 2d
	template <class TGrid, class TVector>
	enable_if_real_host_vector<TVector, void>
	shift_2d(Stream<e_host> &stream, FFT2<Value_type<TGrid>, e_host> &fft2, TGrid &grid, 
	r2d<Value_type<TGrid>> p, TVector &Im)
	{	
		using T = Value_type<TVector>;
		vector<complex<T>> Im_c(Im.begin(), Im.end());

		T Im_min = min_data(Im);

		shift_2d(stream, fft2, grid, p, Im_c);
		assign_real(stream, Im_c, Im);

		if(Im_min>=0)
		{
			std::for_each(Im.begin(), Im.end(), [](T &v){ v = (v<0)?0:v; });
		}
	}

	// phase correlation function
	template <class TGrid, class TVector>
	TVector PCF(Stream<e_host> &stream, FFT2<Value_type<TGrid>, e_host> &fft2, 
	TGrid &grid, TVector &M1_i, TVector &M2_i, Value_type<TGrid> sigma_g, Border<Value_type<TGrid>> bd1, 
	Border<Value_type<TGrid>> bd2)
	{
		using T = Value_type<TVector>;
		using TVector_c = Vector<complex<T>, e_host>;

		// apply Hanning filter and copy data to complex matrix
		auto diff_x_Hanning = [&](Grid<T> &grid, TVector &M_i, Border<T> bd)->TVector_c
		{
			const T x_c = bd.x_c();
			const T y_c = bd.y_c();

			const T Radius = 0.9*::fmin(bd.lx_wb(), bd.ly_wb())/2;
			const T R02 = pow(Radius, 2);
			const int n = 4;

			TVector_c M(grid.nxy());

			auto thr_diff_x_Hanning = [=](const Range &range, Grid<T> &grid, TVector &M_i, TVector_c &M)
			{
				for (auto ix = range.ix_0; ix < range.ix_e; ix++)
				{
					for (auto iy = range.iy_0; iy < range.iy_e; iy++)
					{
						int ixy = grid.ind_col(ix, iy);
						int ix_n = (ix+1<grid.nx)?(ix+1):ix;
						int ixy_n = grid.ind_col(ix_n, iy);
						T Rx = grid.Rx(ix);
						T Ry = grid.Ry(iy);
						T fxy = (bd.chk_bound(Rx, Ry))?1.0/(1.0+pow((pow(Rx-x_c, 2)+pow(Ry-y_c, 2))/R02, n)):0;
						M[ixy] = complex<T>((M_i[ixy_n] - M_i[ixy])*fxy);
					}

				}
			};

			stream.set_n_act_stream(grid.nx);
			stream.set_grid(grid.nx, grid.ny);
			stream.exec(thr_diff_x_Hanning, grid, M_i, M);

			return M;
		};

		TVector_c M1 = diff_x_Hanning(grid, M1_i, bd1);
		TVector_c M2 = diff_x_Hanning(grid, M2_i, bd2);

		// create 2d plan
		fft2.create_plan_2d(grid.ny, grid.nx, stream.size());

		// shift matrix
		fft2_shift(stream, grid, M1);
		fft2_shift(stream, grid, M2);

		// fft2
		fft2.forward(M1);
		fft2.forward(M2);

		// apply symmetric bandwidth
		auto thr_PCF = [&](const Range &range)
		{
			T alpha = 0.5/(sigma_g*sigma_g);

			for(auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for(auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					int ixy = grid.ind_col(ix, iy);
					complex<T> z = conj(M1[ixy])*M2[ixy];
					auto mz = abs(z);
					M2[ixy] = (isZero(mz))?0:z*exp(-alpha*grid.g2_shift(ix, iy))/mz;
				}
			}
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_PCF);

		fft2.inverse(M2);

		TVector M_o(grid.nxy());
		assign_real(stream, M2, M_o);

		// shift M
		fft2_shift(stream, grid, M_o);

		std::for_each(M_o.begin(), M_o.end(), [&](T &v){ v = ::fmax(v, 0); });

		return M_o;
	}

	// get shift
	template <class TGrid, class TVector>
	r2d<Value_type<TGrid>> find_shift_2d(Stream<e_host> &stream, 
	FFT2<Value_type<TGrid>, e_host> &fft2, TGrid &grid, TVector &M1_i, 
	TVector &M2_i, Value_type<TGrid> sigma_g, Border<Value_type<TGrid>> bd1, 
	Border<Value_type<TGrid>> bd2, int nit=2)
	{
		using T = Value_type<TVector>;

		T sigma_r = ::fmax(3*grid.dR_min(), 1.0/(c_2Pi*sigma_g));
		T radius = 1.25*sigma_r;

		TVector M = M2_i;
		r2d<T> dr(0,0);
		for (auto it=0; it<nit; it++)
		{
			auto pcf = PCF(stream, fft2, grid, M1_i, M, sigma_g, bd1, bd2);

			// get maximum index position
			int ixy_max = std::max_element(pcf.begin(), pcf.end())-pcf.begin();

			int ix_max, iy_max;
			grid.col_row(ixy_max, ix_max, iy_max);

			r2d<T> p(grid.Rx(ix_max), grid.Ry(iy_max));

			p = max_pos(grid, pcf, p, sigma_r, radius);
			r2d<T> dr_t = p - r2d<T>(grid.lxh(), grid.lyh());

			if(it<nit-1)
			{
				shift_2d(stream, fft2, grid, -dr_t, M);
			}
			dr += dr_t;
			//calculated new borders
			bd1.shift(-dr_t);
			bd2.shift(-dr_t);
		}

		return dr;
	}

	// correct shift
	template <class TGrid, class TVector>
	r2d<Value_type<TVector>> correct_shift_2d(Stream<e_host> &stream, FFT2<Value_type<TGrid>, e_host> &fft2, 
	TGrid &grid, TVector &M1_i, TVector &M2_io, Value_type<TGrid> sigma_g, Border<Value_type<TGrid>> bd1, 
	Border<Value_type<TGrid>> bd2, int nit=3)
	{
		using T = Value_type<TVector>;

		auto dr = find_shift_2d(stream, fft2, grid, M1_i, M2_io, sigma_g, bd1, bd2, nit);
		shift_2d(stream, fft2, grid, -dr, M2_io);

		return dr;
	}

	/*******************************************************************/
	// get projective standard deviation
	template <class TGrid, class TVector>
	TVector projected_intensity(TGrid &grid, TVector &M_i, Value_type<TVector> np_min, Value_type<TVector> delta)
	{
		using T = Value_type<TVector>;

		vector<r2d<T>> pts_c = {r2d<T>(0, 0), r2d<T>(grid.nx-1, 0), r2d<T>(grid.nx-1, grid.ny-1), r2d<T>(0, grid.ny-1)};

		auto n_pp = static_cast<int>(ceil(module(r2d<T>(grid.nx, grid.ny))+2));
		TVector y_pt(n_pp);
		TVector c_pt(n_pp);

		T cos_d, sin_d;
		sincos(delta, &sin_d, &cos_d);

		// get projected points
		auto krn_proj_point = [&](const T &cos_d, const T &sin_d, const r2d<T> &p)->r2d<T>
		{
			return r2d<T>(cos_d*p.x+sin_d*p.y, -sin_d*p.x+cos_d*p.y);
		};

		// get reference corner point
		auto p_0 = krn_proj_point(cos_d, sin_d, pts_c[0]);
		for(auto ixy = 1; ixy < pts_c.size(); ixy++)
		{
			auto p_r = krn_proj_point(cos_d, sin_d, pts_c[ixy]);
			if(p_0.y > p_r.y)
			{
				p_0 = p_r;
			}
		}

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			for(auto iy = 0; iy < grid.ny; iy++)
			{
				auto ixy = grid.ind_col(ix, iy);
				auto p_r = krn_proj_point(cos_d, sin_d, r2d<T>(ix, iy))-p_0;
				auto j = static_cast<int>(floor(p_r.y));
				y_pt[j] += M_i[ixy];
				c_pt[j] += 1;
			}
		}

		TVector V_o;
		V_o.reserve(n_pp);
		for(auto j = 0; j < n_pp; j++)
		{
			if(c_pt[j]>np_min)
			{
				V_o.push_back(y_pt[j]/c_pt[j]);	
			}
		}
		V_o.shrink_to_fit();

		return V_o;
	}

	// get projective standard deviation
	template <class TGrid, class TVector>
	void PSD(Stream<e_host> &stream, TGrid &grid, TVector &M_i, Value_type<TVector> np_min, 
	Value_type<TVector> del_0, Value_type<TVector> del_e, Value_type<TVector> d_del, TVector &x_o, TVector &y_o)
	{
		// get number of angles
		auto n_delta = static_cast<int>(floor((del_e-del_0)/d_del + 0.5));
		x_o.resize(n_delta);
		y_o.resize(n_delta);

		auto thr_psd = [&](const Range &range)
		{
			for(auto idel=range.ixy_0; idel<range.ixy_e; idel++)
			{
				auto delta = (del_0 + idel*d_del)*c_Pi/180;
				x_o[idel] = delta*180/c_Pi;

				auto pIm = projected_intensity(grid, M_i, np_min, delta);

				y_o[idel] = variance(pIm);
			}
		};

		stream.set_n_act_stream(n_delta);
		stream.set_grid(1, n_delta);
		stream.exec(thr_psd);
	}

	// get projective standard deviation
	template <class TGrid, class TVector>
	TVector PSD_find_peaks(Stream<e_host> &stream, TGrid &grid, TVector &M_i, Value_type<TVector> np_min, 
	TVector x_i, TVector y_i)
	{
		using T = Value_type<TVector>;

		T d_del = abs(x_i[1]-x_i[0]);
		int nr = max(1, static_cast<int>(ceil(3.0/d_del)));
		auto y_f = filter_median_1d(stream, y_i, nr);

		T y_mean = 0;
		T y_max = y_i[0] - y_f[0];
		for(auto ix=0; ix<y_i.size(); ix++)
		{
			y_f[ix] = y_i[ix] - y_f[ix];
			y_mean += y_f[ix];
			y_max = max(y_max, y_f[ix]);
		}
		y_mean /= y_i.size();
		auto y_thr = y_mean + 0.2*(y_mean + y_max);

		auto x_peaks = find_peaks_vector_typ_1(x_i, y_f, y_thr);

		TVector x_ref, y_ref;
		for(auto ix=0; ix<x_peaks.size(); ix++)
		{
			auto delta = x_peaks[ix];
			auto delta_0 = delta - d_del;
			auto delta_e = delta + d_del;
			auto d_delta = 0.1*d_del;
			PSD(stream, grid, M_i, np_min, delta_0, delta_e, d_delta, x_ref, y_ref);
			std::for_each(y_ref.begin(), y_ref.end(), [](T &y){y = log(1+y);});
			auto coef = lsf_poly_n(x_ref, y_ref, 2);
			T px = -coef[1]/(2*coef[2]);
			if((px<=delta_0)||(delta_e<=px))
			{
				int idx = std::max_element(y_ref.begin(), y_ref.end())-y_ref.begin();
				px = x_ref[idx];
			}
			x_peaks[ix] = px;
		}

		return x_peaks;
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
			for(int i = ix1; i<ix2-1; i++)
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
			for(int i = ix1; i>ix2; i--)
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
	void get_sigma_2d_by_row(int ny, int nx, double *g, double *fg, double *sigma)
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

		for(int iy = 0; iy < ny; iy++)
		{
			std::fill(frl.begin(), frl.end(), 0.0);
			std::fill(cfrl.begin(), cfrl.end(), 0.0);
			for(int ix = 0; ix < nx; ix++)
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
			for(int i = 1; i < nr; i++)
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
		radial_distribution_2d(nxy, g, fg, nr+1, gl.data(), rl.data(), frl.data(), cfrl.data(), true, 0);
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
		for(int i = 1; i < nr; i++)
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
	double FFT_information_limit_2d(int ny, int nx, bool shift, double *fI)
	{
		int nr = min(nx/2, ny/2)-1;
		double *r, *fIr;

		r = new double[nr];
		fIr = new double[nr]; 

		// Cumulative radial integration
		getCumRadDist_2d(ny, nx, shift, fI, nr, r, fIr);

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

} // namespace mt

#endif