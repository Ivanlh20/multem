/*
 * This file is part of MULTEM.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef CPU_FCNS_H
#define CPU_FCNS_H

#include <thread>
#include <type_traits>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>

#include <fftw3.h>
#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "fft.cuh"
#include "stream.cuh"
#include "lapack.hpp"
#include "lin_alg_def.cuh"
#include "atomic_data_mt.hpp"

#include "cgpu_fcns.cuh"
#include "quadrature.hpp"

namespace mt
{
	// median filter 1d
	template <class TVector>
	TVector ftr_median_1d(Stream<e_host> &stream, TVector &Im_i, int nkr);

	// median filter 2d
	template <class TGrid, class TVector>
	TVector ftr_median_2d(Stream<e_host> &stream, TGrid &grid_2d, TVector &Im_i, int nkr);

	// wiener filter 1d
	template <class TVector>
	TVector ftr_wiener_1d(Stream<e_host> &stream, TVector &Im_i, int nkr);

	// wiener filter 2d
	template <class TGrid, class TVector>
	TVector ftr_wiener_2d(Stream<e_host> &stream, TGrid &grid_2d, TVector &Im_i, int nkr);

	// den poiss 
	template <class TVector>
	TVector ftr_poiss_dnois_1d(Stream<e_host> &stream, TVector &Im_i, int nkr_w, int nkr_m);

	// get peak signal to noise ratio PSNR 
	template <class TVector>
	Value_type<TVector> get_PSNR(Stream<e_host> &stream, TVector &Im_i, TVector &Im_d);

	// gray opening
	template <class TVector>
	TVector morp_g_open(Stream<e_host> &stream, int ny_i, int nx_i, TVector &Im_i, int nkr);

	// thresholding
	template <class TVector>
	TVector thresholding(Stream<e_host> &stream, TVector &v_i, Value_type<TVector> thr);

	template <class T>
	class Neigh_2d;

	namespace host_detail
	{
		template <class TFn, class ...TArg>
		void matrix_iter(const Range_2d &range, TFn &fn, TArg &...arg)
		{
			for (auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for (auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					fn(ix, iy, arg...);
				}
			}
		}

		template <class TFn, class ...TArg>
		void matrix_iter_yx(const Range_2d &range, TFn &fn, TArg &...arg)
		{
			for (auto iy = range.iy_0; iy < range.iy_e; iy++)
			{
				for (auto ix = range.ix_0; ix < range.ix_e; ix++)
				{
					fn(ix, iy, arg...);
				}
			}
		}

		template <class TFn, class ...TArg>
		void vector_iter(const Range_2d &range, TFn &fn, TArg &...arg)
		{
			for (auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				fn(ixy, arg...);
			}
		}

		template <class T>
		T atom_cost_function(const Grid_2d<T> &grid_2d, const Atom_Sa<T> &atom_Ip, rVector<T> M_i)
		{
			T sum = 0;

			for (auto ix_0 = 0; ix_0 < atom_Ip.ixn; ix_0++)
			{
				for (auto iy_0 = 0; iy_0 < atom_Ip.iyn; iy_0++)
				{
					int ix = ix_0 + atom_Ip.ix_0;
					int iy = iy_0 + atom_Ip.iy_0;

					T R2 = grid_2d.R2(ix, iy, atom_Ip.x, atom_Ip.y);
					if (R2 < atom_Ip.R2_max)
					{
						int ixy = grid_2d.ind_col(ix, iy);
						T M = M_i[ixy];
						const T V = host_device_detail::eval_cubic_poly(R2, atom_Ip);
						sum += (V - 2 * M)*V;
					}
				}
			}

			return sum;
		}

		template <class T>
		void subtract_atom(Stream<e_host> &stream, const Grid_2d<T> &grid_2d, const Atom_Sa<T> &atom_Ip, rVector<T> M_i)
		{
			for (auto ix_0 = 0; ix_0 < atom_Ip.ixn; ix_0++)
			{
				int iyc = 0;
				for (auto iy_0 = 0; iy_0 < atom_Ip.iyn; iy_0++)
				{
					int ix = ix_0 + atom_Ip.ix_0;
					int iy = iy_0 + atom_Ip.iy_0;

					T R2 = grid_2d.R2(ix, iy, atom_Ip.x, atom_Ip.y);
					if (R2 < atom_Ip.R2_max)
					{
						int ixy = grid_2d.ind_col(ix, iy);

						const T V = host_device_detail::eval_cubic_poly(R2, atom_Ip);

						atom_Ip.iv[iyc] = ixy;
						atom_Ip.v[iyc] = V;
						iyc++;
					}
				}

				stream.stream_mutex.lock();
				for (auto iy_0 = 0; iy_0 < iyc; iy_0++)
				{
					M_i.V[atom_Ip.iv[iy_0]] -= atom_Ip.v[iy_0];
				}
				stream.stream_mutex.unlock();
			}
		}

		// Linear projected potential: V and zV
		template <ePotential_Type potential_type, int charge, class TAtom>
		void linear_Vz(const Q1<Value_type<TAtom>, e_host> &qz, TAtom &atom)
		{
			using T = Value_type<TAtom>;

			for (auto iR = 0; iR < c_nR; iR++)
			{
				T R2 = atom.R2[iR];
				T V = 0;
				T dVir = 0;

				T a = (atom.split) ? (-atom.z0h) : (atom.zeh - atom.z0h);
				T b = (atom.split) ? (atom.z0h) : (atom.zeh + atom.z0h);
				for (auto ix = 0; ix < qz.size(); ix++)
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
					for (auto ix = 0; ix < qz.size(); ix++)
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
		template <class TAtom>
		void cubic_poly_coef(TAtom &atom)
		{
			for (auto iR = 0; iR < c_nR - 1; iR++)
			{
				host_device_detail::cubic_poly_coef(iR, atom);
			}
		}

		// Cubic polynomial evaluation
		template <class T>
		void eval_cubic_poly(Stream<e_host> &stream, Grid_2d<T> &grid_2d, Atom_Vp<T> &atom, rVector<T> M_o)
		{
			for (auto ix_0 = 0; ix_0 < atom.nx; ix_0++)
			{
				int iyc = 0;
				for (auto iy_0 = 0; iy_0 < atom.ny; iy_0++)
				{
					int ix = ix_0 + atom.ix_0;
					int iy = iy_0 + atom.iy_0;
					const auto R2 = grid_2d.R2(ix, iy, atom.x, atom.y);

					if (R2 < atom.R2_max)
					{
						const T V = atom.occ*host_device_detail::eval_cubic_poly(R2, atom);
						const int ixy = grid_2d.ind_col_pbc_shift(ix, iy);

						atom.iv[iyc] = ixy;
						atom.v[iyc] = V;
						iyc++;
					}
				}

				stream.stream_mutex.lock();
				for (auto iy_0 = 0; iy_0 < iyc; iy_0++)
				{
					M_o.V[atom.iv[iy_0]] += atom.v[iy_0];
				}
				stream.stream_mutex.unlock();
			}
		}

		// Gaussian evaluation
		template <class T>
		void gauss_eval(Stream<e_host> &stream, Grid_2d<T> &grid_2d, Gauss_Sp<T> &gauss, rVector<T> M_o)
		{
			for (auto ix_0 = 0; ix_0 < gauss.nx; ix_0++)
			{
				int iyc = 0;
				for (auto iy_0 = 0; iy_0 < gauss.ny; iy_0++)
				{
					int ix = ix_0 + gauss.ix_0;
					int iy = iy_0 + gauss.iy_0;

					T R2 = grid_2d.R2(ix, iy, gauss.x, gauss.y);
					if (R2 < gauss.R2_max)
					{
						gauss.iv[iyc] = grid_2d.ind_col_pbc(ix, iy);
						gauss.v[iyc] = gauss(R2);
						iyc++;
					}
				}

				stream.stream_mutex.lock();
				for (auto iy_0 = 0; iy_0 < iyc; iy_0++)
				{
					M_o.V[gauss.iv[iy_0]] += gauss.v[iy_0];
				}
				stream.stream_mutex.unlock();
			}
		}

		template <class TGrid>
		Value_type<TGrid> Lorentz_factor(Stream<e_host> &stream, TGrid &grid_2d, EELS<Value_type<TGrid>> &eels)
		{
			using T_r = Value_type<TGrid>;
			T_r sum = 0;

			auto thr_Lorentz_factor = [&](const Range_2d &range)
			{
				T_r sum_partial = 0;
				matrix_iter(range, host_device_detail::Lorentz_factor<TGrid>, grid_2d, eels.gc2, eels.ge2, sum_partial);

				stream.stream_mutex.lock();
				sum += sum_partial;
				stream.stream_mutex.unlock();
			};

			stream.set_n_act_stream(grid_2d.nx);
			stream.set_grid(grid_2d.nx, grid_2d.ny);
			stream.exec(thr_Lorentz_factor);

			return sqrt(eels.occ) / sum;
		}

		template <class TVector>
		void fit_log_gauss_a_sigma(const TVector &x2_i, const TVector &y_i, Value_type<TVector> &a1, Value_type<TVector> &a0)
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

			T det = sx1x1*sx2x2 - sx1x2*sx1x2;
			a1 = (sx2x2*sx1y - sx1x2*sx2y) / det;
			a0 = (sx1x1*sx2y - sx1x2*sx1y) / det;
		}

		template <class TVector>
		void fit_gauss_a_c(const TVector &x2_i, const TVector &y_i, Value_type<TVector> lambdai,
			Value_type<TVector> sigma, Value_type<TVector> &a, Value_type<TVector> &c)
		{
			using T = Value_type<TVector>;

			// get a, c
			T c_0 = 0.5 / pow(sigma, 2);
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
			T lambda = lambdai*(sx1x1 + sx2x2);
			sx1x1 += lambda;
			sx2x2 += lambda;

			T det = sx1x1*sx2x2 - sx1x2*sx1x2;
			a = (sx2x2*sx1y - sx1x2*sx2y) / det;
			c = (sx1x1*sx2y - sx1x2*sx1y) / det;
		}

		template <class TVector>
		Value_type<TVector> fit_gauss_a(const TVector &x2_i, const TVector &y_i, Value_type<TVector> sigma)
		{
			using T = Value_type<TVector>;

			// get a = sum xy/sum x^2
			T c_0 = 0.5 / pow(sigma, 2);

			T sx1x1 = 0;
			T sx1y = 0;
			for (auto ix = 0; ix < x2_i.size(); ix++)
			{
				T x1 = exp(-c_0*x2_i[ix]);
				T y = y_i[ix];

				sx1x1 += x1*x1;
				sx1y += x1*y;
			}

			return sx1y / sx1x1;
		}

		template <class TVector>
		Value_type<TVector> fit_gauss_sigma(const TVector &x2_i, const TVector &y_i, Value_type<TVector> a,
			Value_type<TVector> &sigma, Value_type<TVector> sigma_min, Value_type<TVector> sigma_max,
			int nit, Value_type<TVector> d_sigma_error)
		{
			using T = Value_type<TVector>;

			T sigma_o = sigma;
			// get b = sum xy/sum x^2
			for (auto it = 0; it < nit; it++)
			{
				T c_0 = 0.5 / pow(sigma, 2);
				T c_1 = a / pow(sigma, 3);

				T sx1x1 = 0;
				T sx1y = 0;
				for (auto ix = 0; ix < x2_i.size(); ix++)
				{
					T f = exp(-c_0*x2_i[ix]);
					T x1 = c_1*x2_i[ix] * f;
					T y = y_i[ix] - a*f;

					sx1x1 += x1*x1;
					sx1y += x1*y;
				}
				T d_sigma = sx1y / sx1x1;
				sigma += d_sigma;
				// sigma = min(max(sigma, sigma_min), sigma_max);

				if (sigma <= sigma_min)
				{
					sigma = sigma_min;
					break;
				}

				if (sigma >= sigma_max)
				{
					sigma = sigma_max;
					break;
				}

				if (fabs(d_sigma / sigma) < d_sigma_error)
				{
					break;
				}
			}

			return sigma - sigma_o;
		}

		template <class TVector>
		Value_type<TVector> fit_gauss_w_a(const TVector &x2_i, const TVector &y_i, Value_type<TVector> sigma)
		{
			using T = Value_type<TVector>;

			// get a = sum xy/sum x^2
			T c_0 = 0.5 / pow(sigma, 2);

			T sx1x1 = 0;
			T sx1y = 0;
			for (auto ix = 0; ix < x2_i.size(); ix++)
			{
				// T w = 1/::fmax(x2_i[ix], 0.001);
				T w = x2_i[ix];
				T x1 = exp(-c_0*x2_i[ix]);
				T y = w*y_i[ix];

				sx1x1 += w*x1*x1;
				sx1y += x1*y;
			}

			return sx1y / sx1x1;
		}

		template <class TVector>
		Value_type<TVector> fit_gauss_w_sigma(const TVector &x2_i, const TVector &y_i, Value_type<TVector> a,
			Value_type<TVector> &sigma, Value_type<TVector> sigma_min, Value_type<TVector> sigma_max,
			int nit, Value_type<TVector> d_sigma_error)
		{
			using T = Value_type<TVector>;

			T sigma_o = sigma;
			// get b = sum xy/sum x^2
			for (auto it = 0; it < nit; it++)
			{
				T c_0 = 0.5 / pow(sigma, 2);
				T c_1 = a / pow(sigma, 3);

				T sx1x1 = 0;
				T sx1y = 0;
				for (auto ix = 0; ix < x2_i.size(); ix++)
				{
					// T w = 1/::fmax(x2_i[ix], 0.001);
					T w = x2_i[ix];
					T f = exp(-c_0*x2_i[ix]);
					T x1 = c_1*x2_i[ix] * f;
					T y = y_i[ix] - a*f;

					sx1x1 += w*x1*x1;
					sx1y += w*x1*y;
				}
				T d_sigma = sx1y / sx1x1;
				sigma += d_sigma;

				if (sigma <= sigma_min)
				{
					sigma = sigma_min;
					break;
				}

				if (sigma >= sigma_max)
				{
					sigma = sigma_max;
					break;
				}

				if (fabs(d_sigma / sigma) < d_sigma_error)
				{
					break;
				}
			}

			return sigma - sigma_o;
		}

	} // host_detail

	/***************************************************************************/
	/***************************************************************************/

	template <class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
		assign(TVector_1 &M_i, TVector_2 &M_o, Vector<Value_type<TVector_2>, e_host> *M_i_h = nullptr)
	{
		M_o.assign(M_i.begin(), M_i.end());
	}

	template <class TVector_1, class TVector_2>
	typename std::enable_if<is_host_vector_and_device_vector<TVector_1, TVector_2>::value
		&& is_complex<Value_type<TVector_2>>::value && !std::is_same<Value_type<TVector_1>, Value_type<TVector_2>>::value, void>::type
		assign(TVector_1 &M_i, TVector_2 &M_o, Vector<Value_type<TVector_2>, e_host> *M_i_h = nullptr)
	{
		Vector<Value_type<TVector_2>, e_host> M_h;
		M_i_h = (M_i_h == nullptr) ? &M_h : M_i_h;

		// copy data to the same output type
		assign(M_i, *M_i_h);

		// data transfer from CPU to GPU
		M_o.assign(M_i_h->begin(), M_i_h->end());
	}

	template <class TVector_1, class TVector_2>
	typename std::enable_if<is_host_vector_and_device_vector<TVector_1, TVector_2>::value
		&& (!is_complex<Value_type<TVector_2>>::value || std::is_same<Value_type<TVector_1>, Value_type<TVector_2>>::value), void>::type
		assign(TVector_1 &M_i, TVector_2 &M_o, Vector<Value_type<TVector_2>, e_host> *M_i_h = nullptr)
	{
		M_o.assign(M_i.begin(), M_i.end());
	}

	template <class TVector>
	enable_if_host_vector<TVector, void>
	fill(Stream<e_host> &stream, TVector &M_io, Value_type<TVector> value_i)
	{
		auto thr_fill = [&](const Range_2d &range)
		{
			thrust::fill(M_io.begin() + range.ixy_0, M_io.begin() + range.ixy_e, value_i);
		};

		stream.set_n_act_stream(M_io.size());
		stream.set_grid(1, M_io.size());
		stream.exec(thr_fill);
	}

	template <class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	scale(Stream<e_host> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_scale = [&](const Range_2d &range)
		{
			thrust::transform(M_i.begin() + range.ixy_0, M_i.begin() + range.ixy_e,
				M_o.begin() + range.ixy_0, functor::scale<value_type>(w_i));
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_scale);
	}

	template <class TVector>
	enable_if_host_vector<TVector, void>
	scale(Stream<e_host> &stream, Value_type<TVector> w_i, TVector &M_io)
	{
		scale(stream, w_i, M_io, M_io);
	}

	template <class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	square(Stream<e_host> &stream, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_square = [&](const Range_2d &range)
		{
			thrust::transform(M_i.begin() + range.ixy_0, M_i.begin() + range.ixy_e,
				M_o.begin() + range.ixy_0, functor::square<value_type>());
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_square);
	}

	template <class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	square_scale(Stream<e_host> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_square_scale = [&](const Range_2d &range)
		{
			thrust::transform(M_i.begin() + range.ixy_0, M_i.begin() + range.ixy_e,
				M_o.begin() + range.ixy_0, functor::square_scale<value_type>(w_i));
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_square_scale);
	}

	template <class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	add(Stream<e_host> &stream, TVector_1 &M1_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_add = [&](const Range_2d &range)
		{
			thrust::transform(M1_i.begin() + range.ixy_0, M1_i.begin() + range.ixy_e,
				M2_i.begin() + range.ixy_0, M_o.begin() + range.ixy_0, functor::add<value_type>());
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_add);
	}

	template <class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	add(Stream<e_host> &stream, TVector_1 &M_i, TVector_2 &M_io)
	{
		add(stream, M_i, M_io, M_io);
	}

	template <class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	add_scale(Stream<e_host> &stream, Value_type<TVector_2> w1_i, TVector_1 &M1_i,
	Value_type<TVector_2> w2_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_add_scale = [&](const Range_2d &range)
		{
			thrust::transform(M1_i.begin() + range.ixy_0, M1_i.begin() + range.ixy_e,
				M2_i.begin() + range.ixy_0, M_o.begin() + range.ixy_0, functor::add_scale_i<value_type>(w1_i, w2_i));
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_add_scale);
	}

	template <class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	add_scale(Stream<e_host> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_io)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_add_scale = [&](const Range_2d &range)
		{
			thrust::transform(M_i.begin() + range.ixy_0, M_i.begin() + range.ixy_e,
				M_io.begin() + range.ixy_0, M_io.begin() + range.ixy_0, functor::add_scale<value_type>(w_i));
		};

		stream.set_n_act_stream(M_io.size());
		stream.set_grid(1, M_io.size());
		stream.exec(thr_add_scale);
	}

	template <class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	add_square(Stream<e_host> &stream, TVector_1 &M1_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_add_square = [&](const Range_2d &range)
		{
			thrust::transform(M1_i.begin() + range.ixy_0, M1_i.begin() + range.ixy_e,
				M2_i.begin() + range.ixy_0, M_o.begin() + range.ixy_0, functor::add_square_i<value_type>());
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_add_square);
	}

	template <class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
		add_square(Stream<e_host> &stream, TVector_1 &M_i, TVector_2 &M_io)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_add_square = [&](const Range_2d &range)
		{
			thrust::transform(M_i.begin() + range.ixy_0, M_i.begin() + range.ixy_e,
				M_io.begin() + range.ixy_0, M_io.begin() + range.ixy_0, functor::add_square<value_type>());
		};

		stream.set_n_act_stream(M_io.size());
		stream.set_grid(1, M_io.size());
		stream.exec(thr_add_square);
	}

	template <class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
		add_scale_square(Stream<e_host> &stream, Value_type<TVector_2> w1_i, TVector_1 &M1_i, Value_type<TVector_2> w2_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_add_scale_square = [&](const Range_2d &range)
		{
			thrust::transform(M1_i.begin() + range.ixy_0, M1_i.begin() + range.ixy_e,
				M2_i.begin() + range.ixy_0, M_o.begin() + range.ixy_0, functor::add_scale_square_i<value_type>(w1_i, w2_i));
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_add_scale_square);
	}

	template <class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
		add_scale_square(Stream<e_host> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_io)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_add_scale_square = [&](const Range_2d &range)
		{
			thrust::transform(M_i.begin() + range.ixy_0, M_i.begin() + range.ixy_e,
				M_io.begin() + range.ixy_0, M_io.begin() + range.ixy_0, functor::add_scale_square<value_type>(w_i));
		};

		stream.set_n_act_stream(M_io.size());
		stream.set_grid(1, M_io.size());
		stream.exec(thr_add_scale_square);
	}

	template <class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
		multiply(Stream<e_host> &stream, TVector_1 &M1_i, TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_multiply = [&](const Range_2d &range)
		{
			thrust::transform(M1_i.begin() + range.ixy_0, M1_i.begin() + range.ixy_e,
				M2_i.begin() + range.ixy_0, M_o.begin() + range.ixy_0, functor::multiply<value_type>());
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_multiply);
	}

	template <class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
		multiply(Stream<e_host> &stream, TVector_1 &M_i, TVector_2 &M_io)
	{
		multiply(stream, M_i, M_io, M_io);
	}

	template <class TVector>
	enable_if_host_vector<TVector, Value_type<TVector>>
		sum(Stream<e_host> &stream, TVector &M_i)
	{
		using value_type = Value_type<TVector>;

		value_type sum_total = 0;
		value_type sum_ee = 0;
		auto thr_sum = [&](const Range_2d &range)
		{
			auto sum_partial = thrust::reduce(M_i.begin() + range.ixy_0, M_i.begin() + range.ixy_e);

			stream.stream_mutex.lock();
			host_device_detail::kh_sum(sum_total, sum_partial, sum_ee);
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(M_i.size());
		stream.set_grid(1, M_i.size());
		stream.exec(thr_sum);

		return sum_total;
	}

	template <class TVector>
	enable_if_host_vector<TVector, Value_type_r<TVector>>
		sum_square(Stream<e_host> &stream, TVector &M_i)
	{
		using T_r = Value_type_r<TVector>;

		T_r sum_total = 0;
		T_r sum_ee = 0;
		auto thr_sum_square = [&](const Range_2d &range)
		{
			auto sum_partial = thrust::transform_reduce(M_i.begin() + range.ixy_0, M_i.begin() + range.ixy_e,
				functor::square<T_r>(), T_r(0), functor::add<T_r>());

			stream.stream_mutex.lock();
			host_device_detail::kh_sum(sum_total, sum_partial, sum_ee);
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(M_i.size());
		stream.set_grid(1, M_i.size());
		stream.exec(thr_sum_square);

		return sum_total;
	}

	template <class TVector>
	enable_if_host_vector<TVector, Value_type<TVector>>
		mean(Stream<e_host> &stream, TVector &M_i)
	{
		return sum(stream, M_i) / M_i.size();
	}

	template <class TVector>
	enable_if_host_vector<TVector, void>
		mean_var(Stream<e_host> &stream, TVector &M_i, Value_type<TVector> &x_mean, Value_type_r<TVector> &x_var)
	{
		using T = Value_type<TVector>;
		using T_r = Value_type_r<TVector>;

		x_mean = mean(stream, M_i);

		x_var = 0;
		auto thr_var = [&](const Range_2d &range)
		{
			auto x_var_partial = thrust::transform_reduce(M_i.begin() + range.ixy_0, M_i.begin() + range.ixy_e,
				functor::square_dif<T, T_r>(x_mean), T_r(0), functor::add<T_r>());

			stream.stream_mutex.lock();
			x_var += x_var_partial;
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(M_i.size());
		stream.set_grid(1, M_i.size());
		stream.exec(thr_var);

		x_var = x_var / M_i.size();
	}

	template <class TVector>
	enable_if_host_vector<TVector, Value_type_r<TVector>>
		variance(Stream<e_host> &stream, TVector &M_i)
	{
		using T = Value_type<TVector>;
		using T_r = Value_type_r<TVector>;

		T x_mean;
		T_r x_var;
		mean_var(stream, M_i, x_mean, x_var);
		return x_var;
	}

	template <class TVector>
	void rescale_data(const Value_type<TVector> &x_mean, const Value_type<TVector> &x_std, TVector &x)
	{
		using T = Value_type<TVector>;
		std::for_each(x.begin(), x.end(), [=](T &v) { v = (v - x_mean) / x_std; });
	}

	template <class TVector, class ...TArg>
	void rescale_data(const Value_type<TVector> &x_mean, const Value_type<TVector> &x_std, TVector &x, TArg &...arg)
	{
		rescale_data(x_mean, x_std, x);
		rescale_data(x_mean, x_std, arg...);
	}

	template <class TVector>
	Value_type<TVector> max_std_data(TVector &x)
	{
		return sqrt(variance(x));
	}

	template <class TVector, class ...TArg>
	Value_type<TVector> max_std_data(TVector &x, TArg &...arg)
	{
		return ::fmax(max_std_data(x), max_std_data(arg...));
	}

	/***********************************************************************/
	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		exp_r_factor_1d(TGrid &grid_1d, Value_type<TGrid> gx, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		for (auto ix = 0; ix < grid_1d.nx; ix++)
		{
			host_device_detail::exp_r_factor_1d(ix, grid_1d, gx, fPsi_i, fPsi_o);
		}
	}

	template <class TGrid, class TVector_r, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		exp_r_factor_2d_bc(Stream<e_host> &stream, TGrid &grid_2d, Value_type<TGrid> alpha,
			TVector_r &gy, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::exp_r_factor_2d_bc<TGrid, TVector_r, TVector_c>, grid_2d, alpha, gy, fPsi_i, fPsi_o);
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		exp_r_factor_2d(Stream<e_host> &stream, TGrid &grid_2d, Value_type<TGrid> gx, Value_type<TGrid> gy,
			TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::exp_r_factor_2d<TGrid, TVector_c>, grid_2d, gx, gy, fPsi_i, fPsi_o);
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		mul_exp_r_factor_2d(Stream<e_host> &stream, TGrid &grid_2d, Vector<Value_type<TGrid>, e_host> &gx,
			Vector<Value_type<TGrid>, e_host> &gy, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		using TVector_r = Vector<Value_type<TGrid>, e_host>;

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::mul_exp_r_factor_2d<TGrid, TVector_r, TVector_c>, grid_2d, gx, gy, fPsi_i, fPsi_o);
	}

	/***********************************************************************/
	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		exp_g_factor_1d(TGrid &grid_1d, Value_type<TGrid> x, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		for (auto ix = 0; ix < grid_1d.nx; ix++)
		{
			host_device_detail::exp_g_factor_1d(ix, grid_1d, x, fPsi_i, fPsi_o);
		}
	}

	template <class TGrid, class TVector_r, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		exp_g_factor_2d_bc(Stream<e_host> &stream, TGrid &grid_2d, Value_type<TGrid> alpha,
			TVector_r &y, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::exp_g_factor_2d_bc<TGrid, TVector_r, TVector_c>, grid_2d, alpha, y, fPsi_i, fPsi_o);
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		exp_g_factor_2d(Stream<e_host> &stream, TGrid &grid_2d, Value_type<TGrid> x, Value_type<TGrid> y,
			TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::exp_g_factor_2d<TGrid, TVector_c>, grid_2d, x, y, fPsi_i, fPsi_o);
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		mul_exp_g_factor_2d(Stream<e_host> &stream, TGrid &grid_2d, Vector<Value_type<TGrid>, e_host> &x,
			Vector<Value_type<TGrid>, e_host> &y, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		using TVector_r = Vector<Value_type<TGrid>, e_host>;

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::mul_exp_g_factor_2d<TGrid, TVector_r, TVector_c>, grid_2d, x, y, fPsi_i, fPsi_o);
	}

	/***********************************************************************/
	template <class TGrid, class TVector_r>
	enable_if_host_vector<TVector_r, Value_type<TVector_r>>
		atom_cost_function(TGrid &grid_2d, const Atom_Sa<Value_type<TGrid>> &atom_Ip, TVector_r &M_i)
	{
		return host_detail::atom_cost_function<typename TVector_r::value_type>(grid_2d, atom_Ip, M_i);
	}

	template <class TGrid, class TVector_r>
	enable_if_host_vector<TVector_r, void>
		subtract_atom(Stream<e_host> &stream, TGrid &grid_2d, Vector<Atom_Sa<Value_type<TGrid>>, e_host> &atom_Ip, TVector_r &M_i)
	{
		if (stream.n_act_stream <= 0)
		{
			return;
		}

		for (auto istream = 0; istream < stream.n_act_stream; istream++)
		{
			stream[istream] = std::thread(std::bind(host_detail::subtract_atom<typename TVector_r::value_type>, std::ref(stream), std::ref(grid_2d), std::ref(atom_Ip[istream]), std::ref(M_i)));
		}
		stream.synchronize();
	}

	// Linear projected potential: V and zV
	template <class TQ1, class TVAtom>
	enable_if_host<TQ1, void>
		linear_Vz(Stream<e_host> &stream, ePotential_Type potential_type, TQ1 &qz, TVAtom &vatom)
	{
		using TAtom = Value_type<TVAtom>;

		if (stream.n_act_stream <= 0)
		{
			return;
		}

		auto thr_linear_Vz = [](const ePotential_Type &potential_type, TQ1 &qz, TAtom &atom)
		{
			if (atom.charge == 0)
			{
				switch (potential_type)
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
				switch (potential_type)
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

		for (auto istream = 0; istream < stream.n_act_stream - 1; istream++)
		{
			stream[istream] = std::thread(std::bind(thr_linear_Vz, potential_type, std::ref(qz), std::ref(vatom[istream])));
		}

		thr_linear_Vz(potential_type, qz, vatom[stream.n_act_stream - 1]);

		stream.synchronize();
	}

	// Get Local interpolation coefficients
	template <class TVAtom>
	enable_if_host<typename TVAtom::value_type, void>
	cubic_poly_coef(Stream<e_host> &stream, TVAtom &vatom)
	{
		using TAtom = Value_type<TVAtom>;

		if (stream.n_act_stream <= 0)
		{
			return;
		}

		for (auto istream = 0; istream < stream.n_act_stream - 1; istream++)
		{
			stream[istream] = std::thread(std::bind(host_detail::cubic_poly_coef<TAtom>, std::ref(vatom[istream])));
		}

		host_detail::cubic_poly_coef<TAtom>(vatom[stream.n_act_stream - 1]);

		stream.synchronize();
	}

	template <class TGrid, class TVector>
	enable_if_host_vector<TVector, void>
	fft1_shift(TGrid &grid_1d, TVector &M_io)
	{
		for (auto ix = 0; ix < grid_1d.nxh; ix++)
		{
			host_device_detail::fft1_shift(ix, grid_1d, M_io);
		}
	}

	template <class TGrid, class TVector>
	enable_if_host_vector<TVector, void>
	fft2_sft_bc(Stream<e_host> &stream, TGrid &grid_2d, TVector &M_io)
	{
		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.nyh);
		stream.exec_matrix(host_device_detail::fft2_sft_bc<TGrid, TVector>, grid_2d, M_io);
	}

	template <class TGrid, class TVector>
	enable_if_host_vector<TVector, void>
	fft2_shift(Stream<e_host> &stream, TGrid &grid_2d, TVector &M_io)
	{
		stream.set_n_act_stream(grid_2d.nxh);
		stream.set_grid(grid_2d.nxh, grid_2d.nyh);
		stream.exec_matrix(host_device_detail::fft2_shift<TGrid, TVector>, grid_2d, M_io);
	}

	template <class TGrid, class TVector>
	enable_if_host_vector<TVector, Value_type<TVector>>
	sum_over_Det(Stream<e_host> &stream, TGrid &grid_2d, Value_type<TGrid> g_min, Value_type<TGrid> g_max, TVector &M_i)
	{
		using T_r = Value_type<TGrid>;
		using value_type = Value_type<TVector>;

		T_r g2_min = pow(g_min, 2);
		T_r g2_max = pow(g_max, 2);
		value_type sum = 0;

		auto thr_sum_over_Det = [&](const Range_2d &range)
		{
			value_type sum_partial = 0;
			host_detail::matrix_iter(range, host_device_detail::sum_over_Det<TGrid, TVector>, grid_2d, g2_min, g2_max, M_i, sum_partial);

			stream.stream_mutex.lock();
			sum += sum_partial;
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec(thr_sum_over_Det);

		return sum;
	}

	template <class TGrid, class TVector>
	enable_if_host_vector<TVector, Value_type<TGrid>>
		sum_square_over_Det(Stream<e_host> &stream, TGrid &grid_2d, Value_type<TGrid> g_min, Value_type<TGrid> g_max, TVector &M_i)
	{
		using T_r = Value_type<TGrid>;

		T_r g2_min = pow(g_min, 2);
		T_r g2_max = pow(g_max, 2);
		T_r sum = 0;

		auto thr_sum_square_over_Det = [&](const Range_2d &range)
		{
			T_r sum_partial = 0;
			host_detail::matrix_iter(range, host_device_detail::sum_square_over_Det<TGrid, TVector>, grid_2d, g2_min, g2_max, M_i, sum_partial);

			stream.stream_mutex.lock();
			sum += sum_partial;
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec(thr_sum_square_over_Det);

		return sum;
	}

	template <class TGrid, class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, Value_type<TGrid>>
		sum_square_over_Det(Stream<e_host> &stream, TGrid &grid_2d, TVector_1 &S_i, TVector_2 &M_i)
	{
		using T_r = Value_type<TGrid>;

		T_r sum = 0;

		auto thr_sum_square_over_Det = [&](const Range_2d &range)
		{
			T_r sum_partial = 0;
			host_detail::matrix_iter(range, host_device_detail::sum_square_over_Det<TGrid, TVector_1, TVector_2>, grid_2d, S_i, M_i, sum_partial);

			stream.stream_mutex.lock();
			sum += sum_partial;
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec(thr_sum_square_over_Det);

		return sum;
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		bandwidth_limit(Stream<e_host> &stream, TGrid &grid_2d, TVector_c &M_io)
	{
		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::bandwidth_limit<TGrid, TVector_c>, grid_2d, M_io);
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		hard_aperture(Stream<e_host> &stream, TGrid &grid_2d, Value_type<TGrid> g_max, Value_type<TGrid> w, TVector_c &M_io)
	{
		auto g2_max = pow(g_max, 2);

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::hard_aperture<TGrid, TVector_c>, grid_2d, g2_max, w, M_io);
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		propagate(Stream<e_host> &stream, TGrid &grid_2d, Value_type<TGrid> w,
			Value_type<TGrid> gxu, Value_type<TGrid> gyu, TVector_c &psi_i, TVector_c &psi_o)
	{
		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::propagate<TGrid, TVector_c>, grid_2d, w, gxu, gyu, psi_i, psi_o);
	}

	template <class TGrid, class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
		transmission_function(Stream<e_host> &stream, TGrid &grid_2d, eElec_Spec_Int_Model elec_spec_int_model,
			Value_type<TGrid> w, TVector_1 &V0_i, TVector_2 &Trans_o)
	{
		using T_r = Value_type<TGrid>;

		auto thr_transmission_funtion = [&](const Range_2d &range)
		{
			thrust::transform(V0_i.begin() + range.ixy_0, V0_i.begin() + range.ixy_e,
				Trans_o.begin() + range.ixy_0, functor::transmission_function<T_r>(w, elec_spec_int_model));
		};

		stream.set_n_act_stream(grid_2d.nxy());
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec(thr_transmission_funtion);
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		probe(Stream<e_host> &stream, TGrid &grid_2d, Lens<Value_type<TGrid>> &lens, Value_type<TGrid> x,
			Value_type<TGrid> y, Value_type<TGrid> gxu, Value_type<TGrid> gyu, TVector_c &fPsi_o)
	{
		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::probe<TGrid, TVector_c>, grid_2d, lens, x, y, gxu, gyu, fPsi_o);

		auto total = sum_square(stream, fPsi_o);
		mt::scale(stream, sqrt(1.0 / total), fPsi_o);
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		apply_CTF(Stream<e_host> &stream, TGrid &grid_2d, Lens<Value_type<TGrid>> &lens, Value_type<TGrid> gxu, Value_type<TGrid> gyu, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::apply_CTF<TGrid, TVector_c>, grid_2d, lens, gxu, gyu, fPsi_i, fPsi_o);
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		apply_PCTF(Stream<e_host> &stream, TGrid &grid_2d, Lens<Value_type<TGrid>> &lens, TVector_c &fPsi_i, TVector_c &fPsi_o)
	{
		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::apply_PCTF<TGrid, TVector_c>, grid_2d, lens, fPsi_i, fPsi_o);
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		kernel_xyz(Stream<e_host> &stream, TGrid &grid_2d, EELS<Value_type<TGrid>> &eels, FFT<Value_type<TGrid>, e_host> &fft_2d, TVector_c &k_x, TVector_c &k_y, TVector_c &k_z)
	{
		eels.factor = host_detail::Lorentz_factor(stream, grid_2d, eels);

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::kernel_xyz<TGrid, TVector_c>, grid_2d, eels, k_x, k_y, k_z);

		fft_2d.inverse(k_x);
		fft_2d.inverse(k_y);
		fft_2d.inverse(k_z);
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		kernel_x(Stream<e_host> &stream, TGrid &grid_2d, EELS<Value_type<TGrid>> &eels, FFT<Value_type<TGrid>, e_host> &fft_2d, TVector_c &k_x)
	{
		eels.factor = host_detail::Lorentz_factor(stream, grid_2d, eels);

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::kernel_x<TGrid, TVector_c>, grid_2d, eels, k_x);

		fft_2d.inverse(k_x);
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		kernel_y(Stream<e_host> &stream, TGrid &grid_2d, EELS<Value_type<TGrid>> &eels, FFT<Value_type<TGrid>, e_host> &fft_2d, TVector_c &k_y)
	{
		eels.factor = host_detail::Lorentz_factor(stream, grid_2d, eels);

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::kernel_y<TGrid, TVector_c>, grid_2d, eels, k_y);

		fft_2d.inverse(k_y);
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		kernel_z(Stream<e_host> &stream, TGrid &grid_2d, EELS<Value_type<TGrid>> &eels, FFT<Value_type<TGrid>, e_host> &fft_2d, TVector_c &k_z)
	{
		eels.factor = host_detail::Lorentz_factor(stream, grid_2d, eels);

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::kernel_z<TGrid, TVector_c>, grid_2d, eels, k_z);

		fft_2d.inverse(k_z);
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		kernel_mn1(Stream<e_host> &stream, TGrid &grid_2d, EELS<Value_type<TGrid>> &eels, FFT<Value_type<TGrid>, e_host> &fft_2d, TVector_c &k_mn1)
	{
		eels.factor = host_detail::Lorentz_factor(stream, grid_2d, eels);

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::kernel_mn1<TGrid, TVector_c>, grid_2d, eels, k_mn1);

		fft_2d.inverse(k_mn1);
	}

	template <class TGrid, class TVector_c>
	enable_if_host_vector<TVector_c, void>
		kernel_mp1(Stream<e_host> &stream, TGrid &grid_2d, EELS<Value_type<TGrid>> &eels, FFT<Value_type<TGrid>, e_host> &fft_2d, TVector_c &k_mp1)
	{
		eels.factor = host_detail::Lorentz_factor(stream, grid_2d, eels);

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::kernel_mp1<TGrid, TVector_c>, grid_2d, eels, k_mp1);

		fft_2d.inverse(k_mp1);
	}

	/***************************************************************************/
	/***************************************************************************/
	template <class TVector_i, class TVector_o>
	enable_if_host_vector_and_host_vector<TVector_i, TVector_o, void>
	copy_to_host(Stream<e_host> &stream, TVector_i &M_i, TVector_o &M_o,
	TVector_i *M_i_h = nullptr)
	{
		mt::assign(M_i, M_o, M_i_h);
	}

	template <class TVector_i, class TVector_o>
	enable_if_host_vector_and_host_vector<TVector_i, TVector_o, void>
	add_scale_to_host(Stream<e_host> &stream, Value_type<TVector_i> w_i,
	TVector_i &M_i, TVector_o &M_o, TVector_i *M_i_h = nullptr)
	{
		mt::add_scale(stream, w_i, M_i, M_o);
	}

	template <class TVector_i, class TVector_o>
	enable_if_host_vector_and_host_vector<TVector_i, TVector_o, void>
	add_scale_square_to_host(Stream<e_host> &stream, Value_type<TVector_o> w_i,
	TVector_i &M_i, TVector_o &M_o, TVector_i *M_i_h = nullptr)
	{
		mt::add_scale_square(stream, w_i, M_i, M_o);
	}

	template <class TVector_c_i, class TVector_r_o, class TVector_c_o>
	enable_if_host_vector_and_host_vector<TVector_c_i, TVector_c_o, void>
	add_scale_m2psi_psi_to_host(Stream<e_host> &stream, Value_type<TVector_r_o> w_i,
	TVector_c_i &psi_i, TVector_r_o &m2psi_o, TVector_c_o &psi_o, TVector_c_i *psi_i_h = nullptr)
	{
		mt::add_scale(stream, w_i, psi_i, psi_o);
		mt::add_scale_square(stream, w_i, psi_i, m2psi_o);
	}

 	/***************************************************************************/
	/****************************** Host to Host *******************************/
	/***************************************************************************/
  	template <class TGrid, class TVector>
	enable_if_host_vector<TVector, void>
	assign_shift_2d(TGrid &grid_2d, TVector &M_i, TVector &M_o, 
	Vector<Value_type<TVector>, e_host> *M_i_h = nullptr)
	{
		Stream<e_host> stream(1);
		stream.set_n_act_stream(grid_2d.nxh);
		stream.set_grid(grid_2d.nxh, grid_2d.nyh);
		stream.exec_matrix(host_device_detail::assign_shift_2d<TGrid, TVector>, grid_2d, M_i, M_o);	
	}

  	template <class TGrid, class TVector>
	enable_if_host_vector<TVector, void>
	add_scale_shift_2d(TGrid &grid_2d, Value_type<TVector> w, 
	TVector &M_i, TVector &M_o, Vector<Value_type<TVector>, e_host> *M_i_h = nullptr)
	{
		Stream<e_host> stream(1);
		stream.set_n_act_stream(grid_2d.nxh);
		stream.set_grid(grid_2d.nxh, grid_2d.nyh);
		stream.exec_matrix(host_device_detail::add_scale_shift_2d<TGrid, TVector>, grid_2d, w, M_i, M_o);
	}
 
  	template <class TGrid, class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	add_scale_square_shift_2d(TGrid &grid_2d, Value_type<TGrid> w, 
	TVector_1 &M_i, TVector_2 &M_o, Vector<Value_type<TVector_1>, e_host> *M_i_h = nullptr)
	{
		Stream<e_host> stream(1);
		stream.set_n_act_stream(grid_2d.nxh);
		stream.set_grid(grid_2d.nxh, grid_2d.nyh);
		stream.exec_matrix(host_device_detail::add_scale_square_shift_2d<TGrid, TVector_1, TVector_2>, grid_2d, w, M_i, M_o);
	}

	template <class TGrid, class TVector>
	enable_if_host_vector<TVector, void>
	assign_crop(TGrid &grid_2d, TVector &M_i, Range_2d &range, 
	TVector &M_o, Vector<Value_type<TVector>, e_host> *M_i_h = nullptr)
	{
		Stream<e_host> stream(1);
		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec_matrix(host_device_detail::assign_crop<TGrid, TVector>, grid_2d, M_i, range, M_o);	
	}

	template <class TGrid, class TVector>
	enable_if_host_vector<TVector, void>
	assign_crop_shift_2d(TGrid &grid_2d, TVector &M_i, Range_2d &range, 
	TVector &M_o, Vector<Value_type<TVector>, e_host> *M_i_h = nullptr)
	{
		Stream<e_host> stream(1);
		stream.set_n_act_stream(grid_2d.nxh);
		stream.set_grid(grid_2d.nxh, grid_2d.nyh);
		stream.exec_matrix(host_device_detail::assign_crop_shift_2d<TGrid, TVector>, grid_2d, M_i, range, M_o);	
	}

  	template <class TGrid, class TVector>
	enable_if_host_vector<TVector, void>
	add_scale_crop_shift_2d(TGrid &grid_2d, Value_type<TVector> w, 
	TVector &M_i, Range_2d &range, TVector &M_o, Vector<Value_type<TVector>, e_host> *M_i_h = nullptr)
	{
		Stream<e_host> stream(1);
		stream.set_n_act_stream(grid_2d.nxh);
		stream.set_grid(grid_2d.nxh, grid_2d.nyh);
		stream.exec_matrix(host_device_detail::add_scale_crop_shift_2d<TGrid, TVector>, grid_2d, w, M_i, range, M_o);
	}
 
  	template <class TGrid, class TVector_1, class TVector_2>
	enable_if_host_vector_and_host_vector<TVector_1, TVector_2, void>
	add_scale_square_crop_shift_2d(TGrid &grid_2d, Value_type<TGrid> w, 
	TVector_1 &M_i, Range_2d &range, TVector_2 &M_o, Vector<Value_type<TVector_1>, e_host> *M_i_h = nullptr)
	{
		Stream<e_host> stream(1);
		stream.set_n_act_stream(grid_2d.nxh);
		stream.set_grid(grid_2d.nxh, grid_2d.nyh);
		stream.exec_matrix(host_device_detail::add_scale_square_crop_shift_2d<TGrid, TVector_1, TVector_2>, grid_2d, w, M_i, range, M_o);
	}

	/***************************************************************************/
	/**********************************Transpose********************************/

	template <class TVector>
	enable_if_host_vector<TVector, void>
	trs(Stream<e_host> &stream, const int &nrows, const int &ncols, TVector &M)
	{
		TVector M_t(nrows*ncols);
		stream.set_n_act_stream(ncols);
		stream.set_grid(ncols, nrows);
		stream.exec_matrix(host_device_detail::trs<TVector>, ncols, nrows, M, M_t);
		M = M_t;
	}

	/***********************Temporal and Spatial quadratures*********************/
	template <class T>
	void obj_lens_temporal_spatial_quadratures(Lens<T> &lens, Q1<T, e_host> &qt, Q2<T, e_host> &qs)
	{
		/*********************Temporal quadrature**********************/
		Quadrature quadrature;
		quadrature(7, lens.ti_npts, qt); // 7: int_-infty^infty f(x) x^0 Exp[-x^2] dx
		std::for_each(qt.w.begin(), qt.w.end(), [](T &v) { v = v / c_Pii2; });

		/*********************Spatial quadrature**********************/
		qs.reserve((2 * lens.ngxs + 1)*(2 * lens.ngys + 1));
		int iqs = 0;
		T sum_w = 0;
		T sum_ee = 0;
		T alpha = 0.5 / pow(lens.si_sigma, 2);

		for (auto ix = -lens.ngxs; ix <= lens.ngxs; ix++)
		{
			for (auto iy = -lens.ngys; iy <= lens.ngys; iy++)
			{
				T gxs = lens.gxs(ix);
				T gys = lens.gys(iy);
				T g2s = gxs*gxs + gys*gys;
				if (g2s < lens.g2_maxs)
				{
					qs.x.push_back(gxs);
					qs.y.push_back(gys);
					T v = exp(-alpha*g2s);
					qs.w.push_back(v);
					host_device_detail::kh_sum(sum_w, v, sum_ee);
				}
			}
		}
		qs.resize(iqs);
		std::for_each(qs.w.begin(), qs.w.end(), [sum_w](T &v) { v = v / sum_w; });
	}
	
	template <class T>
	void cond_lens_temporal_spatial_quadratures(Lens<T> &lens, Q1<double, e_host> &qt, Q2<double, e_host> &qs)
	{
		/*********************Temporal quadrature**********************/
		bool bb_ti = (lens.ti_npts > 1) && nonZero(lens.ti_sigma);
		Quadrature quadrature;

		if (bb_ti)
		{
			double a = (isZero(lens.ti_sigma))?0:(lens.ti_a/(sqrt(c_2Pi)*lens.ti_sigma));
			double b = (isZero(lens.ti_sigma))?0:(0.5/lens.ti_sigma2());
			double c = (isZero(lens.ti_beta))?0:((1-lens.ti_a)/(2*lens.ti_beta));
			double d = (isZero(lens.ti_beta))?0:(1.0/lens.ti_beta);

			double sigma_t = sqrt(lens.ti_a*lens.ti_sigma2() + 2*(1-lens.ti_a)*lens.ti_beta2());
			double b_q = 0.5/(sigma_t*sigma_t);

			// 26, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
			quadrature(26, lens.ti_npts, qt, 0, 0, 0, b_q);	

			for (auto it = 0; it < lens.ti_npts; it++)
			{
				double r = abs(qt.x[it]); // positive
				double f = a*exp(-b*r*r + b_q*r*r) + c*exp(-d*r + b_q*r*r);
				qt.w[it] *= f;
			}
		}
		else
		{
			qt.resize(1);
			qt.x[0] = 0;
			qt.w[0] = 1;
		}

		/*********************Spatial quadrature**********************/
		bool bb_si = ((lens.si_rad_npts > 1) || (lens.si_azm_npts > 1)) && (nonZero(lens.si_sigma) || nonZero(lens.si_beta));

		if(bb_si)
		{
			double a = (isZero(lens.si_sigma))?0:(lens.si_a/(c_2Pi*lens.si_sigma2()));
			double b = (isZero(lens.si_sigma))?0:(0.5/lens.si_sigma2());
			double c = (isZero(lens.si_beta))?0:((1-lens.si_a)/(c_2Pi*lens.si_beta2()));
			double d = (isZero(lens.si_beta))?0:(1.0/lens.si_beta);

			double sigma_t = sqrt((2*lens.si_a*lens.si_sigma2()+ 6*(1-lens.si_a)*lens.si_beta2())/6);
			double b_q = 1.0/sigma_t;

			// radial part
			Q1<double, e_host> qr;
			// 25, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
			quadrature(25, lens.si_rad_npts, qr, 1, 0, 0, b_q);	

			for (auto ir = 0; ir < lens.si_rad_npts; ir++)
			{
				double r = abs(qr.x[ir]); // positive
				double f = a*exp(-b*r*r + b_q*r) + c*exp(-d*r + b_q*r);
				qr.w[ir] *= f;
			}

			// Azimuth part
			Q1<double, e_host> qa;
			qa.resize(lens.si_azm_npts);
			double h = c_2Pi/lens.si_azm_npts;
			for (auto ia = 0; ia < lens.si_azm_npts; ia++)
			{
				qa.x[ia] = ia*h;
				qa.w[ia] = h;
			}

			qs.reserve(lens.si_rad_npts*lens.si_azm_npts);
			for (auto ir = 0; ir < lens.si_rad_npts; ir++)
			{
				for (auto ia = 0; ia < lens.si_azm_npts; ia++)
				{
					double sin_theta, cos_theta;
					sincos(qa.x[ia], &sin_theta, &cos_theta);
					qs.x.push_back(qr.x[ir]*cos_theta);
					qs.y.push_back(qr.x[ir]*sin_theta);
					qs.w.push_back(qr.w[ir]*qa.w[ia]);
				}
			}
		}
		else
		{
			qs.resize(1);
			qs.x[0] = 0;
			qs.y[0] = 0;
			qs.w[0] = 1;
		}
	}

	/***************************************************************************/
	/***************************************************************************/
	// get index (with typ = 0: bottom index for equal values and typ = 1: upper index for equal values)
	inline
	int getIndex(int ixmin, int ixmax, double *x, int typ, double x0)
	{
		int ixmid;
		switch (typ)
		{
		case 0:
		{
			do {
				ixmid = (ixmin + ixmax) >> 1; 	// divide by 2
				if (x0 <= x[ixmid])
					ixmax = ixmid;
				else
					ixmin = ixmid;
			} while ((ixmax - ixmin) > 1);
		}
		break;
		case 1:
		{
			do {
				ixmid = (ixmin + ixmax) >> 1; 	// divide by 2
				if (x0 < x[ixmid])
					ixmax = ixmid;
				else
					ixmin = ixmid;
			} while ((ixmax - ixmin) > 1);
		}
		break;
		}

		if (x0 == x[ixmax])
			return ixmax;
		else
			return ixmin;
	}

	template <class T>
	void get_bn(const T &R, const int &nR, const T &dR, const T &R_max, const bool &pbc, int &iR0, int &iRn)
	{
		int iR_0 = static_cast<int>(floor((R - R_max) / dR));
		int iR_e = static_cast<int>(ceil((R + R_max) / dR));

		if (!pbc)
		{
			auto set_Bound = [](const int &i, const int &n)->int { return (i < 0) ? 0 : ((i >= n) ? n - 1 : i); };
			iR_0 = set_Bound(iR_0, nR);
			iR_e = set_Bound(iR_e, nR);
		}

		iR0 = iR_0;
		iRn = (iR_0 == iR_e) ? 0 : iR_e - iR_0 + 1;
	}

	template <class TVector>
	TVector lsf_poly_n(TVector &x, TVector &y, int n)
	{
		using T = Value_type<TVector>;

		int m = x.size();
		n++; // add bias term

		TVector M(m*n);
		for (auto in = 0; in < n; in++)
		{
			for (auto im = 0; im < m; im++)
			{
				M[in*m + im] = (in == 0) ? 1 : x[im] * M[(in - 1)*m + im];
			}
		}

		TVector coef(n);
		lapack::GELS<T> gels;
		gels(m, n, M.data(), 1, y.data(), coef.data());

		return coef;
	}

	template<class TVector>
	void rdf_3d(Atom_Data<Value_type<TVector>> &atoms, Value_type<TVector> r_max, int nr, TVector &r, TVector &rdf)
	{
		using T = Value_type<TVector>;

		const T dr = r_max / nr;
		for (auto ir = 0; ir < nr; ir++)
		{
			r[ir] = ir*dr;
		}
		thrust::fill(rdf.begin(), rdf.end(), T(0));

		const T r2_max = pow(r_max, 2);

		const int natoms = atoms.size();
		for (int iatoms = 0; iatoms < natoms; iatoms++)
		{
			const r3d<T> r_i = atoms.to_r3d(iatoms);
			for (int iatoms_s = 0; iatoms_s < natoms; iatoms_s++)
			{
				auto d2 = atoms.norm(iatoms_s, r_i);
				if ((iatoms != iatoms_s) && (d2 < r2_max))
				{
					auto ir = static_cast<int>(floor(sqrt(d2) / dr));
					rdf[ir] += 1;
				}
			}
		}

		rdf[0] = 0;
		for (auto ir = 1; ir < nr; ir++)
		{
			rdf[ir] = rdf[ir] / (4 * c_Pi*pow(r[ir], 2)*dr);
		}
	}

	/********************************************************************/
	// get two dimensional Gaussian by row
	template <class TVector, class TGrid>
	TVector func_gauss_1d(TGrid &grid_1d, Value_type<TVector> sigma, bool shift,
	Border_1d<Value_type<TVector>> bd, bool nb = true)
	{
		using T = Value_type<TVector>;

		Gauss_1d<T> gauss(bd, sigma, nb);
		TVector f(grid_1d.nx);

		for (auto ix = 0; ix < grid_1d.nx; ix++)
		{
			auto Rx2 = (shift) ? grid_1d.R2_shift(ix, gauss.x_c) : grid_1d.R2(ix, gauss.x_c);
			f[ix] = gauss.eval_norm(Rx2);
		}

		return f;
	}

	// get two dimensional Gaussian Filter
	template <class TVector, class TGrid>
	TVector func_gauss_2d_bc(Stream<e_host> &stream, TGrid &grid_2d, Value_type<TVector> sigma, bool shift,
	Border_1d<Value_type<TVector>> bd, bool nb = true)
	{
		using T = Value_type<TVector>;

		Gauss_1d<T> gauss(bd, sigma, nb);
		TVector fy;
		fy.reserve(grid_2d.ny);

		for (auto iy = 0; iy < grid_2d.ny; iy++)
		{
			auto Ry2 = (shift) ? grid_2d.Ry2_shift(iy, gauss.x_c) : grid_2d.Ry2(iy, gauss.x_c);
			fy.push_back(gauss.eval_norm(Ry2));
		}

		TVector f(grid_2d.nxy());

		auto thr_gaussian = [&](const Range_2d &range)
		{
			for (auto ix = range.ixy_0; ix < range.ixy_e; ix++)
			{
				for (auto iy = 0; iy < grid_2d.ny; iy++)
				{
					f[grid_2d.ind_col(ix, iy)] = fy[iy];
				}
			}
		};

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, 1);
		stream.exec(thr_gaussian);

		return f;
	}

	// get two dimensional Gaussian Filter
	template <class TVector, class TGrid>
	TVector func_gauss_2d(Stream<e_host> &stream, TGrid &grid_2d, Value_type<TVector> sigma, bool shift,
	Border_2d<Value_type<TVector>> bd, bool nb = true)
	{
		using T = Value_type<TVector>;

		Gauss_2d<T> gauss(bd, sigma, nb);
		TVector f(grid_2d.nxy());

		auto thr_gaussian = [&](const Range_2d &range)
		{
			for (auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for (auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					T R2 = (shift) ? grid_2d.R2_shift(ix, iy, gauss.x_c, gauss.y_c) : grid_2d.R2(ix, iy, gauss.x_c, gauss.y_c);
					f[grid_2d.ind_col(ix, iy)] = gauss.eval_norm(R2);
				}
			}
		};

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec(thr_gaussian);

		return f;
	}

	/********************************************************************/
	// get one dimensional Hanning Filter
	template <class TVector, class TGrid>
	TVector func_hanning_1d(TGrid &grid_1d, Value_type<TVector> k, bool shift,
	Border_1d<Value_type<TVector>> bd, bool nb = true)
	{
		using T = Value_type<TVector>;

		Hanning_1d<T> hann(bd, k, nb);
		TVector f(grid_1d.nx);

		for (auto ix = 0; ix < grid_1d.nx; ix++)
		{
			T Rx = (shift) ? grid_1d.Rx_shift(ix, hann.x_c) : grid_1d.Rx(ix, hann.x_c);
			f[ix] = hann.eval_norm(Rx);
		}

		return f;
	}

	// get two dimensional Hanning by row
	template <class TVector, class TGrid>
	TVector func_hanning_2d_bc(Stream<e_host> &stream, TGrid &grid_2d, Value_type<TVector> k, bool shift,
	Border_1d<Value_type<TVector>> bd, bool nb = true)
	{
		using T = Value_type<TVector>;

		Hanning_1d<T> hann(bd, k, nb);
		TVector fy;
		fy.reserve(grid_2d.ny);

		for (auto iy = 0; iy < grid_2d.ny; iy++)
		{
			T Ry = (shift) ? grid_2d.Ry_shift(iy, hann.x_c) : grid_2d.Ry(iy, hann.x_c);
			fy.push_back(hann.eval_norm(Ry));
		}

		TVector f(grid_2d.nxy());

		auto thr_hanning = [&](const Range_2d &range)
		{
			for (auto ix = range.ixy_0; ix < range.ixy_e; ix++)
			{
				for (auto iy = 0; iy < grid_2d.ny; iy++)
				{
					f[grid_2d.ind_col(ix, iy)] = fy[iy];
				}
			}
		};

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, 1);
		stream.exec(thr_hanning);

		return f;
	}

	// get two dimensional Hanning Filter
	template <class TVector, class TGrid>
	TVector func_hanning_2d(Stream<e_host> &stream, TGrid &grid_2d, Value_type<TVector> k, bool shift,
	Border_2d<Value_type<TVector>> bd, bool nb = true)
	{
		using T = Value_type<TVector>;

		Hanning_2d<T> hann(bd, k, nb);
		TVector f(grid_2d.nxy());

		auto thr_hanning = [&](const Range_2d &range)
		{
			for (auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for (auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					T R = (shift) ? grid_2d.R_shift(ix, iy, hann.x_c, hann.y_c) : grid_2d.R(ix, iy, hann.x_c, hann.y_c);
					f[grid_2d.ind_col(ix, iy)] = hann.eval_norm(R);
				}
			}
		};

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec(thr_hanning);

		return f;
	}

	/********************************************************************/
	// get one dimensional Butterworth Filter
	template <class TVector, class TGrid>
	TVector func_butterworth_1d(TGrid &grid_1d, Value_type<TVector> radius,
	int n, bool shift, Border_1d<Value_type<TVector>> bd, bool nb = true)
	{
		using T = Value_type<TVector>;

		Butterworth_1d<T> butwth(bd, radius, n, nb);
		TVector f(grid_1d.nx);

		for (auto ix = 0; ix < grid_1d.nx; ix++)
		{
			auto Rx2 = (shift) ? grid_1d.R2_shift(ix, butwth.x_c) : grid_1d.R2(ix, butwth.x_c);
			f[ix] = butwth.eval_norm(Rx2);
		}

		return f;
	}

	// get two dimensional Butterworth Filter
	template <class TVector, class TGrid>
	TVector func_butterworth_2d_bc(Stream<e_host> &stream, TGrid &grid_2d, Value_type<TVector> radius,
	int n, bool shift, Border_1d<Value_type<TVector>> bd, bool nb = true)
	{
		using T = Value_type<TVector>;

		Butterworth_1d<T> butwth(bd, radius, n, nb);
		TVector fy;

		fy.reserve(grid_2d.ny);

		for (auto iy = 0; iy < grid_2d.ny; iy++)
		{
			auto Ry2 = (shift) ? grid_2d.Ry2_shift(iy, butwth.x_c) : grid_2d.Ry2(iy, butwth.x_c);
			fy.push_back(butwth.eval_norm(Ry2));
		}

		TVector f(grid_2d.nxy());

		auto thr_butterworth = [&](const Range_2d &range)
		{
			for (auto ix = range.ixy_0; ix < range.ixy_e; ix++)
			{
				for (auto iy = 0; iy < grid_2d.ny; iy++)
				{
					f[grid_2d.ind_col(ix, iy)] = fy[iy];
				}
			}
		};

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, 1);
		stream.exec(thr_butterworth);

		return f;
	}

	// get two dimensional Butterworth Filter
	template <class TVector, class TGrid>
	TVector func_butterworth_2d(Stream<e_host> &stream, TGrid &grid_2d, Value_type<TVector> radius,
	int n, bool shift, Border_2d<Value_type<TVector>> bd, bool nb = true)
	{
		using T = Value_type<TVector>;

		Butterworth_2d<T> butwth(bd, radius, n, nb);

		TVector f(grid_2d.nxy());

		auto thr_butterworth = [&](const Range_2d &range)
		{
			for (auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for (auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					T R2 = (shift) ? grid_2d.R2_shift(ix, iy, butwth.x_c, butwth.y_c) : grid_2d.R2(ix, iy, butwth.x_c, butwth.y_c);
					f[grid_2d.ind_col(ix, iy)] = butwth.eval_norm(R2);
				}
			}
		};

		stream.set_n_act_stream(grid_2d.nx);
		stream.set_grid(grid_2d.nx, grid_2d.ny);
		stream.exec(thr_butterworth);

		return f;
	}

	/********************************************************************/
	// get two dimensional radial distribution for regular grid_2d
	inline
	void rad_dist_2d(int nR, double *R, double *fR, int nRl, double *Rl, double *rl, double *frl, double *cfrl, bool reg, int typ)
	{
		double Rlmin = Rl[0], Rlmax = Rl[nRl - 1], dRl = Rl[1] - Rl[0];

		for (auto i = 0; i < nRl - 1; i++)
		{
			rl[i] = 0.5*(Rl[i] + Rl[i + 1]);
			frl[i] = 0.0;
			cfrl[i] = 0.0;
		}

		int j;
		for (auto i = 0; i < nR; i++)
		{
			if ((Rlmin <= R[i]) && (R[i] < Rlmax))
			{
				j = (reg) ? (int)floor((R[i] - Rlmin) / dRl) : getIndex(0, nRl - 1, Rl, 0, R[i]);
				frl[j] += fR[i];
				cfrl[j] += 1.0;
			}
		}

		if (typ == 0)
			for (auto i = 0; i < nRl - 1; i++)
			{
				frl[i] /= ::fmax(1.0, cfrl[i]);
			}
	}

	template <class TGrid, class TVector>
	void rad_dist_2d(TGrid &grid_2d, TVector &Im, Value_type<TVector> x_i,
		Value_type<TVector> y_i, Value_type<TVector> radius_i, TVector &rl, TVector &frl)
	{
		using T = Value_type<TVector>;

		T dR = grid_2d.dRx;
		int nrl = static_cast<int>(floor(radius_i / dR + 0.5));
		T R_max = dR*nrl;

		int ix_0, ixe;
		get_bn(x_i, grid_2d.nx, grid_2d.dRx, R_max, grid_2d.pbc_xy, ix_0, ixe);
		ixe += ix_0;

		int iy_0, iye;
		get_bn(y_i, grid_2d.ny, grid_2d.dRy, R_max, grid_2d.pbc_xy, iy_0, iye);
		iye += iy_0;

		// get radial distribution
		rl.resize(nrl);
		frl.resize(nrl);
		TVector cfrl(nrl);

		for (auto ir = 0; ir < rl.size(); ir++)
		{
			rl[ir] = grid_2d.Rx(ir);
			frl[ir] = 0.0;
			cfrl[ir] = 0.0;
		}

		for (auto ix = ix_0; ix < ixe; ix++)
		{
			for (auto iy = iy_0; iy < iye; iy++)
			{
				T R = grid_2d.R(ix, iy, x_i, y_i);
				if (R < R_max)
				{
					auto ir = static_cast<int>(floor(R / dR));
					frl[ir] += Im[grid_2d.ind_col(ix, iy)];
					cfrl[ir] += 1.0;
				}
			}
		}

		for (auto ir = 0; ir < rl.size(); ir++)
		{
			frl[ir] /= ::fmax(1.0, cfrl[ir]);
		}
	}

	// cumulative radial distribution
	template <class TGrid, class TVector>
	void cum_rad_dist_2d(TGrid &grid_2d, TVector &Im, Value_type<TVector> x_i,
		Value_type<TVector> y_i, Value_type<TVector> radius_i, TVector &rl, TVector &frl)
	{
		rad_dist_2d(grid_2d, Im, x_i, y_i, radius_i, rl, frl);
		for (auto ir = 1; ir < frl.size(); ir++)
		{
			frl[ir] += frl[ir - 1];
		}
	}

	// mean smooth
	template <class TVector>
	TVector smooth(TVector &f_i, Value_type<TVector> nkr)
	{
		using T = Value_type<TVector>;

		int nf_i = f_i.size();
		int nk0 = -nkr;
		int nke = nkr + 1;

		TVector f_o(nf_i);

		for (auto ix = 0; ix < nf_i; ix++)
		{
			int it0 = max(ix + nk0, 0);
			int ite = min(ix + nke, nf_i);

			T f_mean = 0;
			for (auto it = it0; it < ite; it++)
			{
				f_mean += f_i[it];
			}

			f_o[ix] = f_mean / (ite - it0);
		}

		return f_o;
	}

	// mean smooth
	template <class TVector>
	TVector smooth(Stream<e_host> &stream, TVector &f_i, Value_type<TVector> nkr)
	{
		using T = Value_type<TVector>;

		int nf_i = f_i.size();
		int nk0 = -nkr;
		int nke = nkr + 1;

		auto krn_mean = [&](const int &ix_i, TVector &f_i, TVector &f_o)
		{
			int ix_0 = max(ix_i + nk0, 0);
			int ixe = min(ix_i + nke, nf_i);

			T f_mean = 0;
			for (auto ix = ix_0; ix < ixe; ix++)
			{
				f_mean += f_i[ix];
			}

			f_o[ix_i] = f_mean / (ixe - ix_0);
		};

		TVector fv(nf_i);

		auto thr_mean = [&](const Range_2d &range)
		{
			for (auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
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
	// add periodic boundary border to the image
	template <class TGrid, class TVector_i, class TVector_o>
	vector<int> add_PB_border(TGrid &grid_i, TVector_i &image_i, int border_x, int border_y, TGrid &grid_o, TVector_o &image_o)
	{
		Prime_Num CP;

		auto ix_0 = border_x;
		auto nx = CP(grid_i.nx + 2 * border_x, eDST_Greater_Than);
		auto ixe = ix_0 + grid_i.nx;
		auto bx_l = border_x;
		auto bx_r = nx - ixe;

		auto iy_0 = border_y;
		auto ny = CP(grid_i.ny + 2 * border_y, eDST_Greater_Than);
		auto iye = iy_0 + grid_i.ny;
		auto by_t = border_y;
		auto by_b = ny - iye;

		grid_o.set_input_data(nx, ny, nx*grid_i.dRx, ny*grid_i.dRy, grid_i.dz, grid_i.bwl, grid_i.pbc_xy);

		image_o.resize(grid_o.nxy());

		// copy central image
		for (auto ix_o = ix_0; ix_o < ixe; ix_o++)
		{
			for (auto iy_o = iy_0; iy_o < iye; iy_o++)
			{
				auto ix_i = ix_o - ix_0;
				auto iy_i = iy_o - iy_0;
				auto ixy_i = grid_i.ind_col(ix_i, iy_i);
				auto ixy_o = grid_o.ind_col(ix_o, iy_o);
				image_o[ixy_o] = image_i[ixy_i];
			}
		}

		// left
		for (auto ix_o = 0; ix_o < bx_l; ix_o++)
		{
			for (auto iy_o = iy_0; iy_o < iye; iy_o++)
			{
				auto ix_i = 2 * bx_l - ix_o;
				auto iy_i = iy_o;
				image_o[grid_o.ind_col(ix_o, iy_o)] = image_o[grid_o.ind_col(ix_i, iy_i)];
			}
		}

		// right
		for (auto ix_o = ixe; ix_o < grid_o.nx; ix_o++)
		{
			for (auto iy_o = iy_0; iy_o < iye; iy_o++)
			{
				auto ix_i = ixe - 2 - (ix_o - ixe);
				auto iy_i = iy_o;
				image_o[grid_o.ind_col(ix_o, iy_o)] = image_o[grid_o.ind_col(ix_i, iy_i)];
			}
		}

		// top
		for (auto ix_o = 0; ix_o < grid_o.nx; ix_o++)
		{
			for (auto iy_o = 0; iy_o < by_t; iy_o++)
			{
				auto ix_i = ix_o;
				auto iy_i = 2 * by_t - iy_o;
				image_o[grid_o.ind_col(ix_o, iy_o)] = image_o[grid_o.ind_col(ix_i, iy_i)];
			}
		}

		// bottom
		for (auto ix_o = 0; ix_o < grid_o.nx; ix_o++)
		{
			for (auto iy_o = iye; iy_o < grid_o.ny; iy_o++)
			{
				auto ix_i = ix_o;
				auto iy_i = iye - 2 - (iy_o - iye);
				image_o[grid_o.ind_col(ix_o, iy_o)] = image_o[grid_o.ind_col(ix_i, iy_i)];
			}
		}

		vector<int> points = { ix_0, ixe, iy_0, iye };

		return points;
	}

	/*******************************************************************/
	// add periodic boundary border to the image
	template <class TGrid, class TVector>
	void set_const_border(TGrid &grid_2d, TVector &image, Value_type<TVector> xb_0, Value_type<TVector> xb_e,
		Value_type<TVector> yb_0, Value_type<TVector> yb_e)
	{
		using T = Value_type<TVector>;

		int ix_0 = grid_2d.ceil_dRx(xb_0);
		int ix_e = grid_2d.nx - grid_2d.ceil_dRx(xb_e);

		int iy_0 = grid_2d.ceil_dRy(yb_0);
		int iy_e = grid_2d.ny - grid_2d.ceil_dRy(yb_e);

		T val_s = 0;
		T val_ee = 0;
		int val_c = 0;
		for (auto ix = ix_0; ix < ix_e; ix++)
		{
			for (auto iy = iy_0; iy < iy_e; iy++)
			{
				auto v = image[grid_2d.ind_col(ix, iy)];
				host_device_detail::kh_sum(val_s, v, val_ee);
				val_c++;
			}
		}
		val_s = val_s / val_c;

		// left
		for (auto ix = 0; ix < ix_0; ix++)
		{
			for (auto iy = 0; iy < grid_2d.ny; iy++)
			{
				image[grid_2d.ind_col(ix, iy)] = val_s;
			}
		}

		// right
		for (auto ix = ix_e; ix < grid_2d.nx; ix++)
		{
			for (auto iy = 0; iy < grid_2d.ny; iy++)
			{
				image[grid_2d.ind_col(ix, iy)] = val_s;
			}
		}

		// top
		for (auto ix = 0; ix < grid_2d.nx; ix++)
		{
			for (auto iy = 0; iy < iy_0; iy++)
			{
				image[grid_2d.ind_col(ix, iy)] = val_s;
			}
		}

		// bottom
		for (auto ix = 0; ix < grid_2d.nx; ix++)
		{
			for (auto iy = iy_e; iy < grid_2d.ny; iy++)
			{
				image[grid_2d.ind_col(ix, iy)] = val_s;
			}
		}
	}

	/*******************************************************************/
	// extract region ix_0<=x<ixe
	template <class TVector_o, class TVector_i>
	TVector_o extract_region_real_part(Stream<e_host> &stream, int nx_src, int ny_src,
		TVector_i &Im_src, int ix0_src, int ixe_src, int iy0_src, int iye_src)
	{
		auto nx_dst = ixe_src - ix0_src;
		auto ny_dst = iye_src - iy0_src;

		TVector_o Im_dst(nx_dst*ny_dst);

		auto thr_extract_region = [&](const Range_2d &range)
		{
			for (auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for (auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					Im_dst[ix*ny_dst + iy] = Im_src[(ix0_src + ix)*ny_src + (iy0_src + iy)].real();
				}
			}
		};

		stream.set_n_act_stream(nx_dst);
		stream.set_grid(nx_dst, ny_dst);
		stream.exec(thr_extract_region);

		return Im_dst;
	}

	// extract real part of vector
	template <class TVector_o, class TVector_i>
	TVector_o extract_real_part(Stream<e_host> &stream, TVector_i &Im_src)
	{
		TVector_o Im_dst(Im_src.size());

		auto thr_extract_real = [](const Range_2d &range, TVector_i &Im_src, TVector_o &Im_dst)
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

	// extract region ix_0<=x<ixe
	template <class TVector_o, class TVector_i>
	TVector_o extract_region_abs(Stream<e_host> &stream, int nx_src, int ny_src,
		TVector_i &Im_src, int ix0_src, int ixe_src, int iy0_src, int iye_src)
	{
		auto nx_dst = ixe_src - ix0_src;
		auto ny_dst = iye_src - iy0_src;

		TVector_o Im_dst(nx_dst*ny_dst);

		auto thr_extract_region = [&](const Range_2d &range)
		{
			for (auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for (auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					Im_dst[ix*ny_dst + iy] = fabs(Im_src[(ix0_src + ix)*ny_src + (iy0_src + iy)]);
				}
			}
		};

		stream.set_n_act_stream(nx_dst);
		stream.set_grid(nx_dst, ny_dst);
		stream.exec(thr_extract_region);

		return Im_dst;
	}

	// extract abs of vector
	template <class TVector_o, class TVector_i>
	TVector_o extract_abs(Stream<e_host> &stream, TVector_i &Im_src)
	{
		TVector_o Im_dst(Im_src.size());

		auto thr_extract_abs = [](const Range_2d &range, TVector_i &Im_src, TVector_o &Im_dst)
		{
			for (auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				Im_dst[ixy] = fabs(Im_src[ixy]);
			}
		};

		stream.set_n_act_stream(Im_src.size());
		stream.set_grid(Im_src.size(), 1);
		stream.exec(thr_extract_abs, Im_src, Im_dst);

		return Im_dst;
	}

	/*******************************************************************/
	// get weighted position
	template <class TGrid, class TVector>
	r2d<Value_type<TVector>> Rx_Ry_weight(TGrid &grid_2d, TVector &Im, r2d<Value_type<TVector>> p_i, Value_type<TVector> radius_i, bool env)
	{
		using T = Value_type<TVector>;

		T dR = grid_2d.dRx;
		int nrl = static_cast<int>(floor(radius_i / dR + 0.5));
		T R_max = dR*nrl;

		auto ix_i = static_cast<int>(floor(p_i.x / dR));
		int ix_0 = max(ix_i - nrl, 0);
		int ixe = min(ix_i + nrl + 1, grid_2d.nx);

		auto iy_i = static_cast<int>(floor(p_i.y / dR));
		int iy_0 = max(iy_i - nrl, 0);
		int iye = min(iy_i + nrl + 1, grid_2d.ny);

		T R2_max = pow(R_max, 2);

		r2d<T> p(0, 0);

		T alpha = 0.5 / pow(R_max, 2);
		T wt = 0;
		for (auto ix = ix_0; ix < ixe; ix++)
		{
			for (auto iy = iy_0; iy < iye; iy++)
			{
				auto R2 = grid_2d.R2(ix, iy, p_i.x, p_i.y);
				if (R2 < R2_max)
				{
					int ixy = grid_2d.ind_col(ix, iy);
					auto imv = ((env) ? exp(-alpha*R2) : 1)*Im[ixy];
					p += imv*r2d<T>(grid_2d.Rx(ix), grid_2d.Ry(iy));
					wt += imv;
				}
			}
		}
		p /= wt;

		if (module(p - p_i) >= radius_i)
		{
			p = p_i;
		}

		return p;
	}

	// fit gaussian 1d
	template <class TGrid, class TVector>
	TVector fit_gauss_1d(TGrid &grid_1d, TVector &Im_i, Value_type<TVector> x_i,
		Value_type<TVector> sigma_i, Value_type<TVector> radius_i)
	{
		using T = Value_type<TVector>;

		auto select_cir_reg = [](TGrid &grid_1d, TVector &Im, T x0,
			T radius, TVector &Rx, TVector &Ix, T &Rx_sf, T &Rx_sc, T &Ix_sc)
		{
			T R_max = radius;
			T R2_max = pow(R_max, 2);

			auto range = grid_1d.index_range(x0, R_max);

			Rx.clear();
			Rx.reserve(range.ixy_e);
			Ix.clear();
			Ix.reserve(range.ixy_e);

			// select circular region
			for (auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				if (grid_1d.R2(ix, x0) < R2_max)
				{
					Rx.push_back(grid_1d.Rx(ix));
					Ix.push_back(Im[ix]);
				}
			}

			Rx_sf = x0;
			Rx_sc = R_max;

			Ix_sc = max_element(Ix);
			int m = Ix.size();
			for (auto ix = 0; ix < m; ix++)
			{
				Rx[ix] = (Rx[ix] - Rx_sf) / Rx_sc;
				Ix[ix] = Ix[ix] / Ix_sc;
			}

			Rx.shrink_to_fit();
			Ix.shrink_to_fit();
		};

		TVector Rx, Ix;
		T Rx_sf, Rx_sc, Ix_sc;

		select_cir_reg(grid_1d, Im_i, x_i, radius_i, Rx, Ix, Rx_sf, Rx_sc, Ix_sc);

		T sigma = sigma_i / Rx_sc;

		int m = Ix.size();
		int n = 3;

		TVector J(m*n);
		TVector d_Ix(m);

		auto get_chi2 = [](TVector &Rx, TVector &Ix, TVector &coef)->T
		{
			T x0 = coef[0];
			T A = coef[1];
			T alpha = 0.5 / pow(coef[2], 2);

			T chi2 = 0;
			T chi2_ee = 0;
			const int m = Ix.size();
			for (auto im = 0; im < m; im++)
			{
				T x = Rx[im] - x0;
				T v = Ix[im] - A*exp(-alpha*x*x);
				host_device_detail::kh_sum(chi2, v*v, chi2_ee);
			}
			return chi2;
		};

		auto get_dIx_J = [](TVector &Rx, TVector &Ix, TVector &coef,
			TVector &dIx, TVector &J)
		{
			T x0 = coef[0];
			T A = coef[1];
			T alpha = 0.5 / pow(coef[2], 2);

			T c_x0 = 1.0 / pow(coef[2], 2);
			T c_A = 1;
			T c_sx = 1.0 / pow(coef[2], 3);

			const int m = Ix.size();
			for (auto im = 0; im < m; im++)
			{
				T x = Rx[im] - x0;
				T v = exp(-alpha*x*x);
				T f = A*v;

				J[0 * m + im] = c_x0*x*f;
				J[1 * m + im] = c_A*v;
				J[2 * m + im] = c_sx*x*x*f;

				dIx[im] = Ix[im] - f;
			}
		};

		T dx_max = ::fmax(grid_1d.dRx / Rx_sc, 0.25*::fmin(sigma, radius_i / Rx_sc));

		vector<T> c_coef_0 = { 0, 1, sigma };
		vector<T> c_coef_min = { -dx_max, 0.5, sigma / 3 };
		vector<T> c_coef_max = { dx_max, 1.25, 3 * sigma };

		TVector coef_0 = c_coef_0;
		TVector coef_min = c_coef_min;
		TVector coef_max = c_coef_max;
		TVector coef = coef_0;

		T chi2 = get_chi2(Rx, Ix, coef);

		T lambda = 2;
		T lambda_f = 2;

		lapack::FLSF<T> flsf;

		const int niter = 100;
		for (auto iter = 0; iter < niter; iter++)
		{
			get_dIx_J(Rx, Ix, coef, d_Ix, J);

			TVector d_coef(n);
			TVector D(n);
			TVector G(n);
			flsf(m, n, J.data(), d_Ix.data(), d_coef.data(), lambda, D.data(), G.data());

			TVector coef_t = coef;
			T rho_f = 0;
			T G_max = 0;
			for (auto ic = 0; ic < n; ic++)
			{
				coef_t[ic] += d_coef[ic];
				rho_f += coef_t[ic] * (D[ic] * coef_t[ic] + G[ic]);
				G_max = ::fmax(G_max, fabs(G[ic]));
			}

			T chi2_t = get_chi2(Rx, Ix, coef_t);
			T rho = (chi2 - chi2_t) / rho_f;

			if ((G_max < 5e-7) || fabs(chi2 - chi2_t) < 5e-8)
			{
				break;
			}

			if (rho > 0)
			{
				coef = coef_t;

				for (auto ic = 0; ic < n; ic++)
				{
					coef[ic] = min(max(coef[ic], coef_min[ic]), coef_max[ic]);
				}

				chi2 = get_chi2(Rx, Ix, coef);

				lambda = (rho > 1e-6) ? (::fmax(lambda / lambda_f, 1e-7)) : lambda;
			}
			else
			{
				lambda = ::fmin(lambda*lambda_f, 1e+7);
			}
		}

		coef[0] = coef[0] * Rx_sc + Rx_sf;
		coef[1] = coef[1] * Ix_sc;
		coef[2] = coef[2] * Rx_sc;

		return coef;
	}

	// fit gaussian 2d
	template <class TGrid, class TVector>
	TVector fit_ellipt_gauss_2d(TGrid &grid_2d, TVector &Im_i, r2d<Value_type<TVector>> p_i,
		Value_type<TVector> sigma_i, Value_type<TVector> radius_i)
	{
		using T = Value_type<TVector>;
		using TRegion = Region<TVector>;

		auto ellipse_var = [](const T &sx, const T &sy, const T &theta, T &a, T &b, T &c)
		{
			T cos_1t = cos(theta);
			T sin_1t = sin(theta);
			T cos_2t = cos(2 * theta);
			T sin_2t = sin(2 * theta);

			T sx2 = pow(sx, 2);
			T sy2 = pow(sy, 2);

			a = 0.5*(cos_1t*cos_1t / sx2 + sin_1t*sin_1t / sy2);
			b = 0.5*(sin_1t*sin_1t / sx2 + cos_1t*cos_1t / sy2);
			c = 0.5*(-sin_2t / sx2 + sin_2t / sy2);
		};

		auto select_cir_reg = [](TGrid &grid_2d, TVector &Im, r2d<T> p,
			T radius, TRegion &region)
		{
			T R_max = radius;
			T R2_max = pow(R_max, 2);
			T Rl2_max = pow(2.5*grid_2d.dR_min(), 2);

			auto range = grid_2d.index_range(p, R_max);

			region.clear();
			region.reserve(range.ixy_e);

			T I_m = 0;
			T I_ee = 0;
			int I_c = 0;
			// select circular region
			for (auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for (auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					T R2 = grid_2d.R2(ix, iy, p.x, p.y);
					if (R2 < R2_max)
					{
						T v = Im[grid_2d.ind_col(ix, iy)];
						region.Rx.push_back(grid_2d.Rx(ix));
						region.Ry.push_back(grid_2d.Ry(iy));
						region.Ixy.push_back(v);
						if (R2 < Rl2_max)
						{
							host_device_detail::kh_sum(I_m, v, I_ee);
							I_c++;
						}
					}

				}
			}
			I_m = I_m / I_c;

			region.R_max = R_max;
			region.Rx_sf = p.x;
			region.Ry_sf = p.y;
			region.Rxy_sc = R_max;

			region.Ixy_sc = (I_m + max_element(region.Ixy)) / 2;
			int m = region.Ixy.size();
			for (auto ixy = 0; ixy < m; ixy++)
			{
				region.Rx[ixy] = (region.Rx[ixy] - region.Rx_sf) / region.Rxy_sc;
				region.Ry[ixy] = (region.Ry[ixy] - region.Ry_sf) / region.Rxy_sc;
				region.Ixy[ixy] = region.Ixy[ixy] / region.Ixy_sc;
			}

			region.Rx.shrink_to_fit();
			region.Ry.shrink_to_fit();
			region.Ixy.shrink_to_fit();
		};

		auto select_ellip_reg = [&](TGrid &grid_2d, TVector &Im, TVector &coef,
			T f0, TRegion &region)
		{
			T x_c = coef[0];
			T y_c = coef[1];

			T a, b, c;
			ellipse_var(coef[3], coef[4], coef[5], a, b, c);

			T d = log(f0);
			T dd = c*c - 4 * a*b;

			T radius_y = sqrt(fabs(4 * a*d / dd));
			T radius_x = sqrt(fabs(4 * b*d / dd));

			T x_0 = x_c - radius_x;
			T x_e = x_c + radius_x;
			T y_0 = y_c - radius_y;
			T y_e = y_c + radius_y;

			int ix_0 = grid_2d.ceil_dRx(x_0);
			int ix_e = grid_2d.floor_dRx(x_e) + 1;

			int iy_0 = grid_2d.ceil_dRy(y_0);
			int iy_e = grid_2d.floor_dRy(y_e) + 1;

			region.clear();
			region.reserve((ix_e - ix_0)*(iy_e - iy_0));

			// select circular region
			for (auto ix = ix_0; ix < ix_e; ix++)
			{
				T dx = grid_2d.Rx(ix) - x_c;
				T ddy = sqrt(fabs(dd*dx*dx - 4 * b*d));

				T y_0 = (-c*dx - ddy) / (2 * b) + y_c;
				T y_e = (-c*dx + ddy) / (2 * b) + y_c;

				int iy_0 = grid_2d.lb_index_y(y_0);
				int iy_e = grid_2d.ub_index_y(y_e) + 1;

				for (auto iy = iy_0; iy < iy_e; iy++)
				{
					region.Rx.push_back(grid_2d.Rx(ix));
					region.Ry.push_back(grid_2d.Ry(iy));
					region.Ixy.push_back(Im[grid_2d.ind_col(ix, iy)]);
				}
			}

			region.R_max = max(coef[3], coef[4]);
			region.Rx_sf = x_c;
			region.Ry_sf = y_c;
			region.Rxy_sc = region.R_max;

			region.Ixy_sc = max_element(region.Ixy);
			int m = region.Ixy.size();
			for (auto ixy = 0; ixy < m; ixy++)
			{
				region.Rx[ixy] = (region.Rx[ixy] - region.Rx_sf) / region.Rxy_sc;
				region.Ry[ixy] = (region.Ry[ixy] - region.Ry_sf) / region.Rxy_sc;
				region.Ixy[ixy] = region.Ixy[ixy] / region.Ixy_sc;
			}

			region.shrink_to_fit();
		};

		TRegion region;

		select_cir_reg(grid_2d, Im_i, p_i, radius_i, region);

		T sigma = sigma_i / region.Rxy_sc;

		auto get_chi2 = [&](TVector &Rx, TVector &Ry, TVector &Ixy, TVector &coef)->T
		{
			T x_0 = coef[0];
			T y_0 = coef[1];

			T theta = coef[5];

			T A = coef[2];

			T a, b, c;

			ellipse_var(coef[3], coef[4], coef[5], a, b, c);

			T chi2 = 0;
			T chi2_ee = 0;
			const int m = Ixy.size();
			for (auto im = 0; im < m; im++)
			{
				T x = Rx[im] - x_0;
				T y = Ry[im] - y_0;
				T v = Ixy[im] - A*exp(-a*x*x - b*y*y - c*x*y);
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
			T cos_2t = cos(2 * theta);
			T sin_2t = sin(2 * theta);

			T sx2 = pow(coef[3], 2);
			T sy2 = pow(coef[4], 2);

			T sx3 = pow(coef[3], 3);
			T sy3 = pow(coef[4], 3);

			T a = 0.5*(cos_1t*cos_1t / sx2 + sin_1t*sin_1t / sy2);
			T b = 0.5*(sin_1t*sin_1t / sx2 + cos_1t*cos_1t / sy2);
			T c = 0.5*(-sin_2t / sx2 + sin_2t / sy2);

			const int m = Ixy.size();
			for (auto im = 0; im < m; im++)
			{
				T x = Rx[im] - x_0;
				T y = Ry[im] - y_0;

				T v = exp(-a*x*x - b*y*y - c*x*y);
				T f = A*v;

				J[0 * m + im] = (2 * a*x + c*y)*f;
				J[1 * m + im] = (2 * b*y + c*x)*f;
				J[2 * m + im] = v;
				J[3 * m + im] = (pow(cos_1t*x, 2) + pow(sin_1t*y, 2) - sin_2t*x*y)*f / sx3;
				J[4 * m + im] = (pow(sin_1t*x, 2) + pow(cos_1t*y, 2) + sin_2t*x*y)*f / sy3;
				J[5 * m + im] = (sin_2t*x*x - sin_2t*y*y + 2 * cos_2t*x*y)*f*(0.5 / sx2 - 0.5 / sy2);

				dIxy[im] = Ixy[im] - f;
			}
		};

		T dRx = grid_2d.dRx / region.Rxy_sc;
		T dRy = grid_2d.dRy / region.Rxy_sc;

		// 0.8 = exp(-0.5*f^2)
		T dx = ::fmax(0.668*sigma, 1.5*dRx);
		T dy = ::fmax(0.668*sigma, 1.5*dRy);

		T sigma_x_min = ::fmax(sigma / 4, 0.5*dRx);
		T sigma_y_min = ::fmax(sigma / 4, 0.5*dRy);

		T sigma_x_max = ::fmin(20 * sigma, grid_2d.lxh() / region.Rxy_sc);
		T sigma_y_max = ::fmin(20 * sigma, grid_2d.lyh() / region.Rxy_sc);

		vector<T> c_coef_0 = { 0, 0, 1, sigma, sigma, 0 };
		vector<T> c_coef_min = { -dx, -dy, 0.5, sigma_x_min, sigma_y_min, -T(c_Pi) };
		vector<T> c_coef_max = { dx, dy, 1.5, sigma_x_max, sigma_y_max, T(c_Pi) };

		TVector coef_0 = c_coef_0;
		TVector coef_min = c_coef_min;
		TVector coef_max = c_coef_max;
		TVector coef = coef_0;

		lapack::FLSF<T> flsf;

		int nit = 4;
		for (auto it = 0; it < nit; it++)
		{
			int m = region.Ixy.size();
			int n = 6;

			TVector J(m*n);
			TVector d_Ixy(m);

			T chi2 = get_chi2(region.Rx, region.Ry, region.Ixy, coef);

			T lambda = 2;
			T lambda_f = 2;

			const int niter = (it == nit - 1) ? 50 : 12;
			for (auto iter = 0; iter < niter; iter++)
			{
				get_dIxy_J(region.Rx, region.Ry, region.Ixy, coef, d_Ixy, J);

				TVector d_coef(n);
				TVector D(n);
				TVector G(n);
				flsf(m, n, J.data(), d_Ixy.data(), d_coef.data(), lambda, D.data(), G.data());

				TVector coef_t = coef;
				T rho_f = 0;
				T G_max = 0;
				for (auto ic = 0; ic < n; ic++)
				{
					coef_t[ic] += d_coef[ic];
					rho_f += coef_t[ic] * (D[ic] * coef_t[ic] + G[ic]);
					G_max = ::fmax(G_max, fabs(G[ic]));
				}

				T chi2_t = get_chi2(region.Rx, region.Ry, region.Ixy, coef_t);
				T rho = (chi2 - chi2_t) / rho_f;

				if ((G_max < 1e-5) || fabs(chi2 - chi2_t) < 1e-7)
				{
					break;
				}

				if (rho > 0)
				{
					coef = coef_t;

					for (auto ic = 0; ic < n - 1; ic++)
					{
						coef[ic] = min(max(coef[ic], coef_min[ic]), coef_max[ic]);
					}
					coef[n - 1] = acos(cos(coef[n - 1]));

					chi2 = get_chi2(region.Rx, region.Ry, region.Ixy, coef);

					lambda = (rho > 1e-6) ? (::fmax(lambda / lambda_f, 1e-7)) : lambda;
				}
				else
				{
					lambda = ::fmin(lambda*lambda_f, 1e+7);
				}
			}

			coef[0] = coef[0] * region.Rxy_sc + region.Rx_sf;
			coef[1] = coef[1] * region.Rxy_sc + region.Ry_sf;
			coef[2] = coef[2] * region.Ixy_sc;
			coef[3] = coef[3] * region.Rxy_sc;
			coef[4] = coef[4] * region.Rxy_sc;

			if (it < nit - 1)
			{
				T radius = ::fmax(coef[3], coef[4]);
				radius = ::fmax(radius_i, radius);

				// select_cir_reg(grid_2d, Im_i, r2d<T>(coef[0], coef[1]), radius, region);

				select_ellip_reg(grid_2d, Im_i, coef, 0.6, region);

				// scale
				coef[0] = (coef[0] - region.Rx_sf) / region.Rxy_sc;
				coef[1] = (coef[1] - region.Ry_sf) / region.Rxy_sc;
				coef[2] = coef[2] / region.Ixy_sc;
				coef[3] = coef[3] / region.Rxy_sc;
				coef[4] = coef[4] / region.Rxy_sc;

				coef_min[0] = -coef[3] / 4;
				coef_min[1] = -coef[4] / 4;
				coef_min[2] = 0.75*coef[2];
				coef_min[3] = coef[3] / 2;
				coef_min[4] = coef[4] / 2;

				coef_max[0] = coef[3] / 4;
				coef_max[1] = coef[4] / 4;
				coef_max[2] = 1.25*coef[2];
				coef_max[3] = 2 * coef[3];
				coef_max[4] = 2 * coef[4];
			}
		}

		return coef;
	}

	template <class TGrid, class TVector>
	Value_type<TVector> neighbor_radius(TGrid &grid_i, TVector &Im_i, r2d<Value_type<TVector>> p_i, Value_type<TVector> radius_i)
	{
		using T = Value_type<TVector>;

		T R2_max = pow(radius_i, 2);

		int ix_0 = grid_i.lb_index_x(p_i.x - radius_i);
		int ixe = grid_i.ub_index_x(p_i.x + radius_i);

		int iy_0 = grid_i.lb_index_y(p_i.y - radius_i);
		int iye = grid_i.ub_index_y(p_i.y + radius_i);

		int nrl = grid_i.ceil_dR_min(radius_i);

		// get radial distribution
		TVector frl(nrl);
		TVector cfrl(nrl);

		T dR = grid_i.dR_min();

		for (auto ix = ix_0; ix < ixe; ix++)
		{
			for (auto iy = iy_0; iy < iye; iy++)
			{
				T R2_d = grid_i.R2(ix, iy, p_i.x, p_i.y);
				if (R2_d < R2_max)
				{
					auto ir = static_cast<int>(floor(sqrt(R2_d) / dR));
					frl[ir] += Im_i[grid_i.ind_col(ix, iy)];
					cfrl[ir] += 1.0;
				}
			}
		}

		for (auto ir = 0; ir < frl.size(); ir++)
		{
			frl[ir] /= ::fmax(1.0, cfrl[ir]);
		}

		frl[0] = ::fmax(frl[0], 1.01*frl[1]);

		frl = smooth(frl, 1);

		// derivative and minimun
		int ir0 = frl.size() - 1;
		for (auto ir = 0; ir < frl.size() - 1; ir++)
		{
			if (frl[ir] < frl[ir + 1])
			{
				ir0 = ir + 1;
				break;
			}
		}

		int irm = frl.size() - 1;
		for (auto ir = ir0; ir < frl.size() - 1; ir++)
		{
			if (frl[ir] > frl[ir + 1])
			{
				irm = ir;
				break;
			}
		}

		return max(2, irm)*dR;
	}

	// select circular region
	template <class TGrid, class TVector>
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
					I_mean += Im_i[grid_i.ind_col(ix, iy)] - bg_i;
					Ic++;
				}
			}
		}
		I_mean /= Ic;
		return I_mean;
	}

	/*******************************************************************/
	template <class TGrid, class TVector>
	TVector interp_profile(TGrid &grid_2d, TVector &Im, const r2d<Value_type<TVector>> &p1, const r2d<Value_type<TVector>> &p2, int nr)
	{
		TVector v;

		if (!(grid_2d.ckb_bound(p1) && grid_2d.ckb_bound(p2)))
		{
			return v;
		}

		if (module(p1 - p2) < grid_2d.dR_min())
		{
			return v;
		}

		using T = Value_type<TGrid>;

		auto interp_2d = [](const r2d<T> &p, TGrid &grid_2d, TVector &Im)->T
		{
			auto ix = grid_2d.lb_index_x(p.x);
			auto iy = grid_2d.lb_index_y(p.y);

			T f11 = Im[grid_2d.ind_col(ix, iy)];
			T f12 = Im[grid_2d.ind_col(ix, iy + 1)];
			T f21 = Im[grid_2d.ind_col(ix + 1, iy)];
			T f22 = Im[grid_2d.ind_col(ix + 1, iy + 1)];

			T x1 = grid_2d.Rx(ix);
			T x2 = grid_2d.Rx(ix + 1);
			T y1 = grid_2d.Ry(iy);
			T y2 = grid_2d.Ry(iy + 1);

			T dx1 = p.x - x1;
			T dx2 = x2 - p.x;
			T dy1 = p.y - y1;
			T dy2 = y2 - p.y;

			T f = (dx2*(f11*dy2 + f12*dy1) + dx1*(f21*dy2 + f22*dy1)) / ((x2 - x1)*(y2 - y1));
			return f;

		};

		r2d<T> p12 = p2 - p1;
		T mp12 = p12.module();

		nr = (nr <= 0) ? static_cast<int>(ceil(mp12 / grid_2d.dR_min())) : nr;
		nr = max(nr, 2);
		T dr = mp12 / (nr - 1);
		r2d<T> u = dr*normalized(p12);

		v.reserve(nr);
		for (auto ir = 0; ir < nr; ir++)
		{
			r2d<T> p = p1 + T(ir)*u;
			v.push_back(interp_2d(p, grid_2d, Im));
		}
		return v;
	}

	template <class TGrid, class TVector>
	void interp_profile(Stream<e_host> &stream, TGrid &grid_2d, TVector &Im, Value_type<TVector> bg, r2d<Value_type<TVector>> p1, r2d<Value_type<TVector>> p2,
		TVector &x, TVector &y)
	{
		x.clear();
		y.clear();

		if (!(grid_2d.ckb_bound(p1) && grid_2d.ckb_bound(p2)))
		{
			return;
		}

		if (module(p1 - p2) < grid_2d.dR_min())
		{
			return;
		}

		using T = Value_type<TGrid>;

		auto interp_2d = [](const r2d<T> &p, TGrid &grid_2d, TVector &Im)->T
		{
			auto ix = grid_2d.lb_index_x(p.x);
			auto iy = grid_2d.lb_index_y(p.y);

			T f11 = Im[grid_2d.ind_col(ix, iy)];
			T f12 = Im[grid_2d.ind_col(ix, iy + 1)];
			T f21 = Im[grid_2d.ind_col(ix + 1, iy)];
			T f22 = Im[grid_2d.ind_col(ix + 1, iy + 1)];

			T x1 = grid_2d.Rx(ix);
			T x2 = grid_2d.Rx(ix + 1);
			T y1 = grid_2d.Ry(iy);
			T y2 = grid_2d.Ry(iy + 1);

			T dx1 = p.x - x1;
			T dx2 = x2 - p.x;
			T dy1 = p.y - y1;
			T dy2 = y2 - p.y;

			T f = (dx2*(f11*dy2 + f12*dy1) + dx1*(f21*dy2 + f22*dy1)) / ((x2 - x1)*(y2 - y1));
			return f;

		};

		int Ixy_1 = Im[grid_2d.ixy(p1.x, p1.y)];
		int Ixy_2 = Im[grid_2d.ixy(p2.x, p2.y)];

		r2d<T> p12 = p2 - p1;
		T mp12 = p12.module();
		int nr = grid_2d.ceil_dR_min(mp12);
		T dr = mp12 / (nr - 1);
		r2d<T> u = dr*normalized(p12);

		x.reserve(nr);
		y.reserve(nr);
		for (auto ir = 0; ir < nr; ir++)
		{
			x.push_back(ir*dr);
			r2d<T> p = p1 + T(ir)*u;
			T v = interp_2d(p, grid_2d, Im);
			y.push_back(v - bg);
		}
	}

	// find one maximum in all areas above the thr
	template <class TVector>
	TVector fd_peaks_vector_typ_1(TVector &x, TVector &y, Value_type<TVector> y_thr)
	{
		using T = Value_type<TVector>;
		TVector x_peak(y.size());

		T x_max = 0;
		T y_max = y_thr;
		bool bb = y[0] > y_thr;
		int ipeak = 0;
		for (auto ix = 0; ix < y.size(); ix++)
		{
			if (y[ix] > y_thr)
			{
				if (y_max < y[ix])
				{
					y_max = y[ix];
					x_max = x[ix];
				}
				bb = true;
			}
			else if (bb)
			{
				x_peak[ipeak++] = x_max;
				y_max = y_thr;
				bb = false;
			}
		}

		if (bb)
		{
			x_peak[ipeak++] = x_max;
		}

		x_peak.resize(ipeak);
		x_peak.shrink_to_fit();

		return x_peak;
	}

	// find all maximum in all areas above the thr
	template <class TVector>
	TVector fd_peaks_vector_typ_2(TVector &x, TVector &y, Value_type<TVector> y_thr)
	{
		using T = Value_type<TVector>;
		TVector x_peak;
		x_peak.reserve(y.size());

		for (auto ix = 1; ix < y.size() - 1; ix++)
		{
			if (y[ix] > y_thr)
			{
				if ((y[ix - 1] < y[ix]) && (y[ix + 1] < y[ix]))
				{
					x_peak.push_back(x[ix]);
				}
			}
		}
		x_peak.shrink_to_fit();

		return x_peak;
	}

	/*******************************************************************/
	// find peak position 1d
	template <class TGrid, class TVector>
	enable_if_host_vector<TVector, Value_type<TVector>>
		fit_max_pos_1d(TGrid &grid_1d, TVector &Im, Value_type<TVector> p_i,
			Value_type<TVector> sigma_i, Value_type<TVector> radius_i)
	{
		using T = Value_type<TVector>;

		auto range = grid_1d.index_range(p_i, radius_i);
		int ix_c = grid_1d.floor_dRx(p_i);

		int ix_0 = range.ix_0;
		for (auto ix = ix_c - 1; ix >= range.ix_0; ix--)
		{
			if (Im[ix - 1] > Im[ix])
			{
				ix_0 = ix;
				break;
			}
		}

		int ix_e = range.ix_e;
		for (auto ix = ix_c + 1; ix <= range.ix_e; ix++)
		{
			if (Im[ix] < Im[ix + 1])
			{
				ix_e = ix;
				break;
			}
		}

		// if there are few points
		if (fabs(ix_e - ix_0) <= 3)
		{
			T x = 0;
			T sI = 0;
			for (auto ix = ix_0; ix <= ix_e; ix++)
			{
				x += Im[ix] * grid_1d.Rx(ix);
				sI += Im[ix];
			}
			x = x / sI;

			return x;
		}

		if ((ix_0 = !range.ix_0) || (ix_e = !range.ix_e))
		{
			T x_min = grid_1d.Rx(ix_0);
			T x_max = grid_1d.Rx(ix_e);
			radius_i = ::fmin(p_i - x_min, x_max - p_i);
		}

		auto coef = fit_gauss_1d(grid_1d, Im, p_i, sigma_i, radius_i);

		return coef[0];
	}

	// find peak position 2d
	template <class TGrid, class TVector>
	enable_if_host_vector<TVector, r2d<Value_type<TVector>>>
		fit_max_pos_2d(TGrid &grid_2d, TVector &Im, r2d<Value_type<TVector>> p_i,
			Value_type<TVector> sigma_i, Value_type<TVector> radius_i)
	{
		using T = Value_type<TVector>;

		auto coef = fit_ellipt_gauss_2d(grid_2d, Im, p_i, sigma_i, radius_i);
		return r2d<T>(coef[0], coef[1]);
	}

	/*******************************************************************/
	// get projective standard deviation
	template <class TGrid, class TVector>
	TVector projected_intensity(TGrid &grid_2d, TVector &M_i, Value_type<TVector> np_min, Value_type<TVector> delta)
	{
		using T = Value_type<TVector>;

		vector<r2d<T>> pts_c = { r2d<T>(0, 0), r2d<T>(grid_2d.nx - 1, 0), r2d<T>(grid_2d.nx - 1, grid_2d.ny - 1), r2d<T>(0, grid_2d.ny - 1) };

		auto n_pp = static_cast<int>(ceil(module(r2d<T>(grid_2d.nx, grid_2d.ny)) + 2));
		TVector y_pt(n_pp);
		TVector c_pt(n_pp);

		T cos_d, sin_d;
		sincos(delta, &sin_d, &cos_d);

		// get projected points
		auto krn_proj_point = [&](const T &cos_d, const T &sin_d, const r2d<T> &p)->r2d<T>
		{
			return r2d<T>(cos_d*p.x + sin_d*p.y, -sin_d*p.x + cos_d*p.y);
		};

		// get reference corner point
		auto p_0 = krn_proj_point(cos_d, sin_d, pts_c[0]);
		for (auto ixy = 1; ixy < pts_c.size(); ixy++)
		{
			auto p_r = krn_proj_point(cos_d, sin_d, pts_c[ixy]);
			if (p_0.y > p_r.y)
			{
				p_0 = p_r;
			}
		}

		for (auto ix = 0; ix < grid_2d.nx; ix++)
		{
			for (auto iy = 0; iy < grid_2d.ny; iy++)
			{
				auto ixy = grid_2d.ind_col(ix, iy);
				auto p_r = krn_proj_point(cos_d, sin_d, r2d<T>(ix, iy)) - p_0;
				auto j = static_cast<int>(floor(p_r.y));
				y_pt[j] += M_i[ixy];
				c_pt[j] += 1;
			}
		}

		TVector V_o;
		V_o.reserve(n_pp);
		for (auto j = 0; j < n_pp; j++)
		{
			if (c_pt[j] > np_min)
			{
				V_o.push_back(y_pt[j] / c_pt[j]);
			}
		}
		V_o.shrink_to_fit();

		return V_o;
	}

	// get projective standard deviation
	template <class TGrid, class TVector>
	void PSD(Stream<e_host> &stream, TGrid &grid_2d, TVector &M_i, Value_type<TVector> np_min,
		Value_type<TVector> del_0, Value_type<TVector> del_e, Value_type<TVector> d_del, TVector &x_o, TVector &y_o)
	{
		// get number of angles
		auto n_delta = static_cast<int>(floor((del_e - del_0) / d_del + 0.5));
		x_o.resize(n_delta);
		y_o.resize(n_delta);

		auto thr_psd = [&](const Range_2d &range)
		{
			for (auto idel = range.ixy_0; idel < range.ixy_e; idel++)
			{
				auto delta = (del_0 + idel*d_del)*c_Pi / 180;
				x_o[idel] = delta * 180 / c_Pi;

				auto pIm = projected_intensity(grid_2d, M_i, np_min, delta);

				y_o[idel] = variance(pIm);
			}
		};

		stream.set_n_act_stream(n_delta);
		stream.set_grid(1, n_delta);
		stream.exec(thr_psd);
	}

	// get projective standard deviation
	template <class TGrid, class TVector>
	TVector PSD_fd_peaks(Stream<e_host> &stream, TGrid &grid_2d, TVector &M_i, Value_type<TVector> np_min,
		TVector x_i, TVector y_i)
	{
		using T = Value_type<TVector>;

		T d_del = fabs(x_i[1] - x_i[0]);
		int nr = max(1, static_cast<int>(ceil(3.0 / d_del)));
		auto y_f = ftr_median_1d(stream, y_i, nr);

		T y_mean = 0;
		T y_max = y_i[0] - y_f[0];
		for (auto ix = 0; ix < y_i.size(); ix++)
		{
			y_f[ix] = y_i[ix] - y_f[ix];
			y_mean += y_f[ix];
			y_max = max(y_max, y_f[ix]);
		}
		y_mean /= y_i.size();
		auto y_thr = y_mean + 0.2*(y_mean + y_max);

		auto x_peaks = fd_peaks_vector_typ_1(x_i, y_f, y_thr);

		TVector x_ref, y_ref;
		for (auto ix = 0; ix < x_peaks.size(); ix++)
		{
			auto delta = x_peaks[ix];
			auto delta_0 = delta - d_del;
			auto delta_e = delta + d_del;
			auto d_delta = 0.1*d_del;
			PSD(stream, grid_2d, M_i, np_min, delta_0, delta_e, d_delta, x_ref, y_ref);
			std::for_each(y_ref.begin(), y_ref.end(), [](T &y) {y = log(1 + y); });
			auto coef = lsf_poly_n(x_ref, y_ref, 2);
			T px = -coef[1] / (2 * coef[2]);
			if ((px <= delta_0) || (delta_e <= px))
			{
				int idx = std::max_element(y_ref.begin(), y_ref.end()) - y_ref.begin();
				px = x_ref[idx];
			}
			x_peaks[ix] = px;
		}

		return x_peaks;
	}

	// get index to maximun distance
	template <class TVector>
	int get_dmaxIndex_Point2Line(r2d<Value_type<TVector>> p1, r2d<Value_type<TVector>> p2,
		int ix1, int ix2, TVector &x, TVector &y)
	{
		using T = Value_type<TVector>;

		T c1 = p2.y - p1.y;
		T c2 = -(p2.x - p1.x);
		T c3 = p2.x*p1.y - p2.y*p1.x;
		T d = sqrt(c1*c1*+c2*c2);

		c1 /= d;
		c2 /= d;
		c3 /= d;

		int ix_max = 0;
		T dmax = 0;
		for (auto ir = ix1; ir < ix2; ir++)
		{
			T x0 = x[ir];
			T y0 = y[ir];
			T d = fabs(c1*x0 + c2*y0 + c3);
			if (d > dmax)
			{
				dmax = d;
				ix_max = ir;
			}
		}

		return ix_max;
	}

	inline
		double getLengthCurve(int ix1, int ix2, double *x, double *y)
	{
		double x1, y1, x2, y2;
		double d, dx, dy;
		d = 0;
		for (auto i = ix1; i < ix2 - 1; i++)
		{
			x1 = x[i]; y1 = y[i];
			x2 = x[i + 1]; y2 = y[i + 1];
			dx = x2 - x1; dy = y2 - y1;
			d = d + sqrt(dx*dx + dy*dy);
		}
		return d;
	}

	inline
		double getLengthCurve(int ix1, int ix2, double *x, double *y, double lmax, int &il)
	{
		double x1, y1, x2, y2;
		double l, dx, dy;
		if (ix1 < ix2)
		{
			l = 0; il = ix2;
			for (int i = ix1; i < ix2 - 1; i++)
			{
				x1 = x[i]; y1 = y[i];
				x2 = x[i + 1]; y2 = y[i + 1];
				dx = x2 - x1; dy = y2 - y1;
				l = l + sqrt(dx*dx + dy*dy);
				if ((lmax > 0) && (l >= lmax))
				{
					il = i;
					break;
				}
			}
		}
		else
		{
			l = 0; il = ix2;
			for (int i = ix1; i > ix2; i--)
			{
				x1 = x[i - 1]; y1 = y[i - 1];
				x2 = x[i]; y2 = y[i];
				dx = x2 - x1; dy = y2 - y1;
				l = l + sqrt(dx*dx + dy*dy);
				if ((lmax > 0) && (l >= lmax))
				{
					il = i;
					break;
				}
			}
		}
		return l;
	}

	// get information limit for regular grid
	template <class TGrid, class TVector>
	Value_type<TVector> info_limit_2d(TGrid &grid_2d, TVector &fM_i)
	{
		using T = Value_type<TVector>;

		TVector r, fr;

		// cumulative radial integration
		cum_rad_dist_2d(grid_2d, fM_i, grid_2d.lxh(), grid_2d.lyh(), grid_2d.lxh_lyh_min(), r, fr);

		// shift and normalize
		T r0 = r[0], fIr0 = fr[0];
		for (auto ir = 0; ir < fr.size(); ir++)
		{
			r[ir] = (r[ir] - r0) / (r.back() - r0);
			fr[ir] = (fr[ir] - fIr0) / (fr.back() - fIr0);
		}

		// smooth data
		int nkr = max<int>(5, fr.size() / 100);
		fr = smooth(fr, nkr);

		r2d<T> p1(r.front(), fr.front());
		r2d<T> p2(r.back(), fr.back());

		// get maximum curvature point
		auto ir_m = get_dmaxIndex_Point2Line(p1, p2, 0, r.size(), r, fr);

		return grid_2d.Rx(ir_m);
	}

	template <class TVector>
	void unique_vector(TVector &v)
	{
		using T = Value_type<TVector>;

		std::sort(v.begin(), v.end());
		v.erase(std::unique(v.begin(), v.end(), mt::isEqual<T>), v.end());
		v.shrink_to_fit();
	}

	template <class InputIterator, class TVector>
	void match_vectors(InputIterator vr_first, InputIterator vr_last, TVector &vs, Value_type<TVector> dv)
	{
		using T = Value_type<TVector>;

		dv = ::fmax(dv, T(0.1));

		TVector vr(vr_first, vr_last);

		// find min and max values of vr
		T vr_min, vr_max;
		minmax_element(vr, vr_min, vr_max);
		vr_min -= dv;
		vr_max += dv;

		// remove elements of vs outside the range
		auto remove = [vr_min, vr_max](T a)->bool
		{
			return ((a < vr_min) || (a > vr_max));
		};
		vs.erase(std::remove_if(vs.begin(), vs.end(), remove), vs.end());

		Vector<bool, e_host> vr_c(vr.size(), true);
		int m_size = min(vr.size(), vs.size());
		TVector A_d;
		A_d.reserve(m_size);

		for (auto i = 0; i < vs.size(); i++)
		{
			T val = vs[i];
			auto it = std::min_element(vr.begin(), vr.end(), functor::closest_element<T>(val));
			auto imin = static_cast<int>(it - vr.begin());

			if (vr_c[imin])
			{
				vr_c[imin] = false;
				A_d.push_back(*it);
			}
		}
		vs = A_d;
		vs.shrink_to_fit();
	}

	/************************ read matlab binary file ***********************/
	template<class T_i, class TVector>
	enable_if_real<T_i, void>
		read_mat_matrix(std::ifstream &bin_file, int nxy, TVector &matrix)
	{
		using T_o = Value_type<TVector>;

		vector<T_i> vect_r(nxy);
		bin_file.read(reinterpret_cast<char*>(vect_r.data()), nxy * sizeof(T_i));

		for (auto ik = 0; ik < nxy; ik++)
		{
			matrix[ik] = T_o(vect_r[ik]);
		}
	}

	template<class T_i, class TVector>
	enable_if_complex<T_i, void>
		read_mat_matrix(std::ifstream &bin_file, int nxy, TVector &matrix)
	{
		using T = Value_type_r<T_i>;
		using T_o = Value_type_r<TVector>;

		vector<T> vect_r(nxy);
		bin_file.read(reinterpret_cast<char*>(vect_r.data()), nxy * sizeof(T));
		vector<T> vect_i(nxy);
		bin_file.read(reinterpret_cast<char*>(vect_i.data()), nxy * sizeof(T));

		auto matrix_ptr = reinterpret_cast<complex<T_o>*>(matrix.data());

		for (auto ik = 0; ik < nxy; ik++)
		{
			matrix_ptr[ik].real(vect_r[ik]);
			matrix_ptr[ik].imag(vect_i[ik]);
		}
	}

	template<class TVector>
	void read_mat_binary_matrix(const char *filename, Grid_2d<Value_type_r<TVector>> &grid_2d, TVector &matrix)
	{
		int nx, ny;
		double dx, dy;
		int type;

		std::ifstream bin_file(filename, std::ofstream::binary);
		bin_file.read(reinterpret_cast<char*>(&nx), sizeof(int));
		bin_file.read(reinterpret_cast<char*>(&ny), sizeof(int));
		bin_file.read(reinterpret_cast<char*>(&dx), sizeof(double));
		bin_file.read(reinterpret_cast<char*>(&dy), sizeof(double));
		bin_file.read(reinterpret_cast<char*>(&type), sizeof(int));

		grid_2d.set_input_data(nx, ny, nx*dx, ny*dy);

		auto nxy = nx*ny;
		matrix.resize(nxy);

		switch (type)
		{
		case 1:
		{
			read_mat_matrix<float>(bin_file, nxy, matrix);
		}
		break;
		case 2:
		{
			read_mat_matrix<double>(bin_file, nxy, matrix);
		}
		break;
		case 3:
		{
			read_mat_matrix<complex<float>>(bin_file, nxy, matrix);
		}
		break;
		case 4:
		{
			read_mat_matrix<complex<double>>(bin_file, nxy, matrix);
		}
		break;
		}
		bin_file.close();
	}
	/************************ write matlab binary file ***********************/
	template<class TVector>
	enable_if_real_host_vector<TVector, void>
		write_mat_matrix(std::ofstream &bin_file, TVector &matrix)
	{
		//using T = Value_type<TVector>;

		//bin_file.write(reinterpret_cast<char*>(matrix.data()), matrix.size()*sizeof(T));
	}

	template<class TVector>
	enable_if_complex_host_vector<TVector, void>
		write_mat_matrix(std::ofstream &bin_file, TVector &matrix)
	{
		//using T = Value_type_r<TVector>;

		//vector<T> vect_r(nxy);
		//vector<T> vect_i(nxy);
		//for (auto ik = 0; ik<nxy; ik++)
		//{
		// vect_r[ik] = matrix[ik].real();
		// vect_i[ik] = matrix[ik].imag();
		//}

		//bin_file.write(reinterpret_cast<char*>(vect_r.data()), matrix.size()*sizeof(T));
		//bin_file.write(reinterpret_cast<char*>(vect_i.data()), matrix.size()*sizeof(T));
	}

	template<class TVector>
	void write_mat_binary_matrix(const char *filename, Grid_2d<Value_type_r<TVector>> &grid_2d, TVector &matrix)
	{
		//int type = matrix_type<TVector>;

		//std::ofstream bin_file(filename, std::ofstream::binary);
		//bin_file.write(reinterpret_cast<char*>(&(grid_2d.nx)), sizeof(int));
		//bin_file.write(reinterpret_cast<char*>(&(grid_2d.ny)), sizeof(int));
		//bin_file.write(reinterpret_cast<char*>(&(grid_2d.dRx)), sizeof(double));
		//bin_file.write(reinterpret_cast<char*>(&(grid_2d.dRy)), sizeof(double));
		//bin_file.write(reinterpret_cast<char*>(&type), sizeof(int));

		//switch (type)
		//{
		//case 1:
		//{
		// write_mat_matrix<float>(bin_file, matrix);
		//}
		//break;
		//case 2:
		//{
		// write_mat_matrix<double>(bin_file, matrix);
		//}
		//break;
		//case 3:
		//{
		// write_mat_matrix<complex<float>>(bin_file, matrix);
		//}
		//break;
		//case 4:
		//{
		// write_mat_matrix<complex<double>>(bin_file, matrix);
		//}
		//break;
		//}
		//bin_file.close();
	}

	/************************ extract real vector form complex vector**********************/
	template<class TVector_c, class TVector_r>
	void from_complex_to_real(eShow_CData show_data, TVector_c &cdata, TVector_r &data)
	{
		switch (show_data)
		{
		case eSCD_CReal:
		{
			for (auto ixy = 0; ixy < cdata.size(); ixy++)
				data[ixy] = cdata[ixy].real();
		}
		break;
		case eSCD_CImag:
		{
			for (auto ixy = 0; ixy < cdata.size(); ixy++)
				data[ixy] = cdata[ixy].imag();
		}
		break;
		case eSCD_CMod:
		{
			for (auto ixy = 0; ixy < cdata.size(); ixy++)
				data[ixy] = abs(cdata[ixy]);
		}
		break;
		case eSCD_CPhase:
		{
			for (auto ixy = 0; ixy < cdata.size(); ixy++)
				data[ixy] = arg(cdata[ixy]);
		}
		break;
		}
	}
} // namespace mt

#endif
