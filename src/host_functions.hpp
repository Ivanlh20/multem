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

namespace mt
{
	// median filter 1d
	template<class TVector>
	TVector filter_median_1d(Stream<e_host> &stream, TVector &Im_i, int nkr);

	// wiener filter 1d
	template<class TVector>
	TVector filter_wiener_1d(Stream<e_host> &stream, TVector &Im_i, int nkr);

	// denoising poisson 
	template<class TVector>
	TVector filter_denoising_poisson_1d(Stream<e_host> &stream, TVector &Im_i, int nkr_w, int nkr_m);

	// get peak signal to noise ratio PSNR 
	template<class TGrid, class TVector>
	Value_type<TVector> get_PSNR(Stream<e_host> &stream, TGrid &grid, TVector &Im_i, int nkr_w, int nkr_m);

	// thresholding
	template<class TVector>
	TVector thresholding(Stream<e_host> &stream, TVector &v_i, Value_type<TVector> threshold);

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
	TVector func_hanning_1d(Stream<e_host> &stream, int nx, Value_type<TVector> dx, Value_type<TVector> k, bool shift)
	{	
		TVector f(nx);

		int nxh = nx/2;
		auto cx = c_2Pi/(dx*nx);

		auto thr_hanning = [&](const Range &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				auto Rx = RSC(ixy, nxh, shift)*dx; 
				auto v = 0.5*(1.0-cos(cx*Rx));
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
	TVector func_hanning_2d_by_row(Stream<e_host> &stream, TGrid &grid, Value_type<TVector> k, bool shift)
	{	
		using T = Value_type<TVector>;

		const T cx = c_2Pi/grid.lx;

		TVector fx;
		fx.reserve(grid.nx);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			auto Rx = (shift)?grid.Rx_shift(ix):grid.Rx(ix); 
			fx.push_back(0.5*(1.0-cos(cx*Rx)));
		}

		TVector f(grid.nxy());

		auto thr_hanning = [&](const Range &range)
		{
			for(auto iy = range.ixy_0; iy < range.ixy_e; iy++)
			{
				for(auto ix = 0; ix < grid.nx; ix++)
				{
					auto v = fx[ix];
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
	TVector func_hanning_2d(Stream<e_host> &stream, TGrid &grid, Value_type<TVector> k, bool shift)
	{	
		using T = Value_type<TVector>;

		const T cx = c_2Pi/grid.lx;
		const T cy = c_2Pi/grid.ly;

		TVector fx;
		fx.reserve(grid.nx);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			auto Rx = (shift)?grid.Rx_shift(ix):grid.Rx(ix); 
			fx.push_back(0.5*(1.0-cos(cx*Rx)));
		}

		TVector fy;
		fy.reserve(grid.ny);

		for(auto iy = 0; iy < grid.ny; iy++)
		{
			auto Ry = (shift)?grid.Ry_shift(iy):grid.Ry(iy); 
			fy.push_back(0.5*(1.0-cos(cy*Ry)));
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
	TVector func_gaussian_1d(Stream<e_host> &stream, int nx, Value_type<TVector> dx, Value_type<TVector> sigma, bool shift)
	{	
		TVector f(nx);

		int nxh = nx/2;
		auto cx = 0.5/(sigma*sigma);

		auto thr_gaussian = [&](const Range &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				auto Rx = FSC(ixy, nxh, shift)*dx; 
				f[ixy] = exp(-cx*Rx*Rx);
			}
		};

		stream.set_n_act_stream(nx);
		stream.set_grid(nx, 1);
		stream.exec(thr_gaussian);

		return f;
	}
	
	// get two dimensional Gaussian Filter
	template<class TVector, class TGrid>
	TVector func_gaussian_2d_by_row(Stream<e_host> &stream, TGrid &grid, Value_type<TVector> sigma, bool shift)
	{	
		using T = Value_type<TVector>;

		const T alpha_x = 0.5/(sigma*sigma);

		const r2d<T> rc(grid.lx/2, grid.ly/2);

		TVector fx;
		fx.reserve(grid.nx);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			auto Rx = (shift)?grid.Rx_shift(ix, rc.x):grid.Rx(ix, rc.x); 
			fx.push_back(exp(-alpha_x*Rx*Rx));
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
	TVector func_gaussian_2d(Stream<e_host> &stream, TGrid &grid, Value_type<TVector> sigma, bool shift)
	{	
		using T = Value_type<TVector>;

		const T alpha_x = 0.5/(sigma*sigma);
		const T alpha_y = 0.5/(sigma*sigma);

		const r2d<T> rc(grid.lx/2, grid.ly/2);

		TVector fx;
		fx.reserve(grid.nx);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			auto Rx = (shift)?grid.Rx_shift(ix, rc.x):grid.Rx(ix, rc.x); 
			fx.push_back(exp(-alpha_x*Rx*Rx));
		}

		TVector fy;
		fy.reserve(grid.ny);

		for(auto iy = 0; iy < grid.ny; iy++)
		{
			auto Ry = (shift)?grid.Ry_shift(iy, rc.y):grid.Ry(iy, rc.y); 
			fy.push_back(exp(-alpha_y*Ry*Ry));
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
	TVector func_butterworth_1d(Stream<e_host> &stream, int nx, Value_type<TVector> dx, Value_type<TVector> Radius, int n, bool shift)
	{
		TVector f(nx);

		int nxh = nx/2;
		auto R02 = pow(Radius, 2);

		auto thr_butterworth = [&](const Range &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				auto Rx = FSC(ixy, nxh, shift)*dx; 
				f[ixy] = 1.0/(1.0+pow(Rx*Rx/R02, n));
			}
		};

		stream.set_n_act_stream(nx);
		stream.set_grid(nx, 1);
		stream.exec(thr_butterworth);

		return f;
	}

	// get two dimensional Butterworth Filter
	template<class TVector, class TGrid>
	TVector func_butterworth_2d_by_row(Stream<e_host> &stream, TGrid &grid, Value_type<TVector> Radius, int n, bool shift)
	{
		using T = Value_type<TVector>;

		const T R02 = pow(Radius, 2);

		const r2d<T> rc(grid.lx/2, grid.ly/2);

		TVector fx;
		fx.reserve(grid.nx);

		for(auto ix = 0; ix < grid.nx; ix++)
		{
			auto Rx = (shift)?grid.Rx_shift(ix, rc.x):grid.Rx(ix, rc.x);
			fx.push_back(1.0/(1.0+pow((Rx*Rx)/R02, n)));
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
	TVector func_butterworth_2d(Stream<e_host> &stream, TGrid &grid, Value_type<TVector> Radius, int n, bool shift)
	{
		using T = Value_type<TVector>;

		const T R02 = pow(Radius, 2);

		const r2d<T> rc(grid.lx/2, grid.ly/2);

		TVector f(grid.nxy());

		auto thr_butterworth = [&](const Range &range)
		{
			for(auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for(auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					auto R2 = (shift)?grid.R2_shift(ix, iy, rc.x, rc.y):grid.R2(ix, iy, rc.x, rc.y); 
					f[grid.ind_col(ix, iy)] = 1.0/(1.0+pow(R2/R02, n));
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
	TGrid &grid_i, Value_type<TGrid> sigma_r, TVector &image_i)
	{	
		using T = Value_type<TVector>;
		using TVector_c = vector<complex<T>>;

		// denoise image
		int nkr_w = max(1, static_cast<int>(floor(sigma_r/2+0.5)));
		int nkr_m = 0;

		// peak signal to noise ratio
		auto PSNR = get_PSNR(stream, grid_i, image_i, nkr_w, nkr_m);

		// deconvolution
		T alpha = 2.0*c_Pi2*sigma_r*sigma_r;
		auto sigma_r_lim = 4.0*sigma_r;
		auto border_x = grid_i.nx_dRx(sigma_r_lim);
		auto border_y = grid_i.ny_dRy(sigma_r_lim);

		TGrid grid;
		TVector_c image;
		// add borders
		auto pts = add_PB_border(grid_i, image_i, border_x, border_y, grid, image);

		//sigma_r ---> sigma_f = 1/(2*pi*sigma_r)
		T sigma_f = 4.0/(c_2Pi*sigma_r);
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
	TGrid &grid_i, Value_type<TGrid> sigma_r, TVector &image_i)
	{	
		using T = Value_type<TVector>;
		using TVector_c = vector<complex<T>>;

		// denoise image
		auto nkr_w = max(1, static_cast<int>(floor(sigma_r/2+0.5)));
		auto nkr_m = 0;

		// peak signal to noise ratio
		auto PSNR = get_PSNR(stream, grid_i, image_i, nkr_w, nkr_m);

		// deconvolution
		T sigma_a = sigma_r;
		T sigma_b = (0.5)*sigma_r;
		T alpha_a = 2.0*c_Pi2*pow(sigma_a, 2);
		T alpha_b = 2.0*c_Pi2*pow(sigma_b, 2);

		auto sigma_r_lim = 4.0*sigma_a;
		auto border_x = grid_i.nx_dRx(sigma_r_lim);
		auto border_y = grid_i.ny_dRy(sigma_r_lim);

		TGrid grid;
		TVector_c image;
		// add borders
		auto pts = add_PB_border(grid_i, image_i, border_x, border_y, grid, image);

		//sigma_r ---> sigma_f = 1/(2*pi*sigma_r)
		T sigma_f_lim = 4.0*(1.0/(c_2Pi*sigma_r));
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
	TVector Rx_Ry_fit(TGrid &grid, TVector &Im_i, r2d<Value_type<TVector>> p_i, Value_type<TVector> radius_i, bool env)
	{
		using T = Value_type<TVector>;
		TVector fit_par(6,1);

		auto get_gaussian_parameters = [&](const r2d<T> &p, TVector &Im_i)->T
		{
			T dR = grid.dRx;
			int nrl = 2;
			auto R2_max = pow(dR*nrl, 2);

			auto ix_i = grid.floor_dRx(p.x);
			int ix0 = max(ix_i-nrl, 0);
			int ixe = min(ix_i+nrl+1, grid.nx);

			auto iy_i = grid.floor_dRy(p.y);
			int iy0 = max(iy_i-nrl, 0);
			int iye = min(iy_i+nrl+1, grid.ny);

			T A = 0;
			int ic = 0;
			for(auto ix = ix0; ix < ixe; ix++)
			{
				for(auto iy = iy0; iy < iye; iy++)
				{
					auto R2 = grid.R2(ix, iy, p.x, p.y);
					if(R2 < R2_max)
					{
						A += Im_i[grid.ind_col(ix, iy)];
						ic++;
					}
				}
			}
			return A/ic;
		};

		T dR = grid.dRx;
		int nrl = static_cast<int>(floor(radius_i/dR+0.5));
		T R_max = dR*nrl;

		auto ix_i = grid.floor_dRy(p_i.x);
		int ix0 = max(ix_i-nrl, 0);
		int ixe = min(ix_i+nrl+1, grid.nx);

		auto iy_i = grid.floor_dRy(p_i.y);
		int iy0 = max(iy_i-nrl, 0);
		int iye = min(iy_i+nrl+1, grid.ny);

		T R2_max = pow(R_max, 2);
		//get number of elements
		int m = 0;
		int n = 6;
		for(auto ix = ix0; ix < ixe; ix++)
		{
			for(auto iy = iy0; iy < iye; iy++)
			{
				m += (grid.R2(ix, iy, p_i.x, p_i.y) < R2_max)?1:0;
			}
		}

		if(m<7)
		{
			auto p = Rx_Ry_weight(grid, Im_i, p_i, radius_i, env);
			fit_par[0] = p.x;
			fit_par[1] = p.y;
			fit_par[2] = get_gaussian_parameters(p, Im_i);
			fit_par[3] = 1/(2*pow(radius_i, 2));
			fit_par[4] = 1/(2*pow(radius_i, 2));
			fit_par[5] = 0;
			return fit_par;
		}

		T alpha = 0.5/pow(R_max, 2);
		TVector M(m*n);
		TVector b(m);
		int ic = 0;
		for(auto ix = ix0; ix < ixe; ix++)
		{
			for(auto iy = iy0; iy < iye; iy++)
			{
				auto R2 = grid.R2(ix, iy, p_i.x, p_i.y);
				if(R2 < R2_max)
				{
					auto Rx = grid.Rx(ix);
					auto Ry = grid.Ry(iy);
					//set M values
					M[0*m+ic] = Rx*Rx;
					M[1*m+ic] = Ry*Ry;
					M[2*m+ic] = Rx;
					M[3*m+ic] = Ry;
					M[4*m+ic] = Rx*Ry;
					M[5*m+ic] = 1.0;

					//set b values
					b[ic] = ((env)?exp(-alpha*R2):1)*Im_i[grid.ind_col(ix, iy)];

					ic++;
				}
			}
		}

		TVector x(n);
		lapack::GELS<T> gels;
		gels(m, n, M.data(), 1, b.data(), x.data());

		auto md = 4*x[0]*x[1]-x[4]*x[4];
		r2d<T> p(-(2*x[1]*x[2]-x[3]*x[4])/md, -(2*x[0]*x[3]-x[2]*x[4])/md);

		if(module(p-p_i)>=radius_i)
		{
			auto p = Rx_Ry_weight(grid, Im_i, p_i, radius_i, env);
			fit_par[0] = p.x;
			fit_par[1] = p.y;
			fit_par[2] = get_gaussian_parameters(p, Im_i);
			fit_par[3] = 1/(2*pow(max<T>(0.5*radius_i, 2),2));
			fit_par[4] = 1/(2*pow(max<T>(0.5*radius_i, 2),2));
			fit_par[5] = 0;
		}
		else
		{
			fit_par[0] = p.x;
			fit_par[1] = p.y;
			fit_par[2] = x[0]*p.x*p.x + x[1]*p.y*p.y + x[2]*p.x + x[3]*p.y+ x[4]*p.x*p.y + x[5];
			fit_par[3] = x[0];
			fit_par[4] = x[1];
			fit_par[5] = x[4];
		}

		return fit_par;
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

	// get fit position
	template<class TGrid, class TVector>
	TVector gaussian_fit(TGrid &grid_i, TVector &Im_i, r2d<Value_type<TVector>> p_i, 
	Value_type<TVector> sigma_i, int niter, Value_type<TVector> d_error)
	{
		using T = Value_type<TVector>;

		T radius_i = 2*sigma_i;
		T R2_max = pow(radius_i, 2);

		int ix0 = grid_i.lb_index_x(p_i.x - radius_i);
		int ixe = grid_i.ub_index_x(p_i.x + radius_i);

		int iy0 = grid_i.lb_index_y(p_i.y - radius_i);
		int iye = grid_i.ub_index_y(p_i.y + radius_i);

		int nxy = (ixe-ix0+1)*(iye-iy0+1);

		TVector R2;
		TVector z;

		R2.reserve(nxy);
		z.reserve(nxy);

		//get number of elements
		for (auto ix = ix0; ix < ixe; ix++)
		{
			for (auto iy = iy0; iy < iye; iy++)
			{
				T R2_d = grid_i.R2(ix, iy, p_i.x, p_i.y);
				if (R2_d < R2_max)
				{
					R2.push_back(R2_d);
					z.push_back(Im_i[grid_i.ind_col(ix, iy)]);
				}
			}
		}

		// set initial guess
		auto z_min_max = std::minmax_element(z.begin(), z.end());

		int m = R2.size();
		int n = 3;

		TVector M(m*2);
		TVector b(m);
		TVector v(m);
		TVector theta(n);
		TVector theta_min(n);
		TVector theta_max(n);
		TVector theta_0(n);
		lapack::GELS<T> gels;

		T f = 0.05;
		// set min values
		theta_min[0] = f*(*(z_min_max.second) - *(z_min_max.first));
		theta_min[1] = ::fmin(0, *(z_min_max.first));
		theta_min[2] = ::fmax(0.25*grid_i.dRx, f*sigma_i);

		// set max values
		theta_max[0] = (1.0 + f)*(*(z_min_max.second) - *(z_min_max.first));
		theta_max[1] = *(z_min_max.second);
		theta_max[2] = ::fmax(0.25*grid_i.dRx, 2.5*sigma_i);

		// set initial guess
		theta_0[0] = 0.95*(*(z_min_max.second) - *(z_min_max.first));
		theta_0[1] = *(z_min_max.first);
		theta_0[2] = ::fmax(0.25*grid_i.dRx, sigma_i);


		auto check_const = [=](TVector &theta)
		{
			for (auto in = 0; in<n; in++)
			{
				if ((theta[in]<theta_min[in]) || (theta_max[in]<theta[in]))
				{
					theta[in] = theta_0[in];
				}
			}
		};

		//a*exp(-r2/(s*b^2)+c)

		theta = theta_0;
		for (auto iter = 0; iter<niter; iter++)
		{
			// get a, c
			T c_0 = 0.5/pow(theta[2], 2);
			for (auto im = 0; im < m; im++)
			{
				v[im] = exp(-c_0*R2[im]);
				M[0*m + im] = v[im];
				M[1*m + im] = 1.0;
				b[im] = z[im];
			}
			gels(m, 2, M.data(), 1, b.data(), theta.data());

			//get b = sum xy/sum x^2
			T c_1 = theta[0]/pow(theta[2], 3);
			T sxy = 0;
			T sxx = 0;
			for (auto im = 0; im < m; im++)
			{
				T x = c_1*R2[im]*v[im];
				T y = z[im]-(theta[0]*v[im]+theta[1]);
				sxy += x*y;
				sxx += x*x;
			}
			T d_sigma = sxy/sxx;
			theta[2] += d_sigma;

			check_const(theta);

			if (abs(d_sigma)<d_error)
			{
				break;
			}
		}
		return theta;
	}

	/*******************************************************************/
	// find 2d peaks
	template<class TGrid, class TVector>
	void find_fast_peaks_2d(Stream<e_host> &stream, FFT2<Value_type<TVector>, e_host> &fft2, 
	TGrid &grid, TVector &image_i, Value_type<TVector> sigma_i, Value_type<TVector> thres, 
	TVector &x_o, TVector &y_o)
	{
		using T = Value_type<TVector>;

		auto Im_minmax = std::minmax_element(image_i.begin(), image_i.end());

		thres = *(Im_minmax.first) + thres*(*(Im_minmax.second)-*(Im_minmax.first));

		auto image = thresholding(stream, image_i, thres);

		// local maximum
		auto krn_maximum = [&](const int &ix_i, const int &iy_i, TVector &Im, r2d<T> &peak)->bool
		{
			auto val_max = Im[grid.ind_col(ix_i, iy_i)];
			peak = r2d<T>(grid.Rx(ix_i), grid.Ry(iy_i));

			if(val_max <= thres)
			{
				 return false;
			}

			const int ix0 = ix_i-1;
			const int ixe = ix_i+2;

			const int iy0 = iy_i-1;
			const int iye = iy_i+2;

			for (auto ix = ix0; ix < ixe; ix++)
			{
				for (auto iy = iy0; iy < iye; iy++)
				{
					if(val_max < Im[grid.ind_col(ix, iy)])
					{
						return false;
					}
				}
			}

			return true;
		};

		auto npeaks_m = static_cast<int>(ceil(grid.lx*grid.ly/(c_Pi*sigma_i*sigma_i)));

		x_o.reserve(2*npeaks_m);
		y_o.reserve(2*npeaks_m);

		// get local peaks
		auto thr_peaks = [&](const Range &range)
		{
			TVector x, y;

			x.reserve(npeaks_m);
			y.reserve(npeaks_m);

			auto ix_0 = 1 + range.ix_0;
			auto ix_e = 1 + range.ix_e;
			auto iy_0 = 1 + range.iy_0;
			auto iy_e = 1 + range.iy_e;

			for(auto ix = ix_0; ix < ix_e; ix++)
			{
				for(auto iy = iy_0; iy < iy_e; iy++)
				{
					r2d<T> peak;
					if(krn_maximum(ix, iy, image, peak))
					{
						x.push_back(peak.x);
						y.push_back(peak.y);
					}
				}
			}

			stream.stream_mutex.lock();
			x_o.insert(x_o.end(), x.begin(), x.end());
			y_o.insert(y_o.end(), y.begin(), y.end());
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(grid.nx-2);
		stream.set_grid(grid.nx-2, grid.ny-2);
		stream.exec(thr_peaks);

		x_o.shrink_to_fit();
		y_o.shrink_to_fit();
	}

	// find 2d peaks
	template<class TGrid, class TVector>
	void find_peaks_2d(Stream<e_host> &stream, FFT2<Value_type<TVector>, e_host> &fft2, 
	TGrid &grid, TVector &image_i, Value_type<TVector> sigma, Value_type<TVector> thres, 
	int niter, TVector &x_o, TVector &y_o, TVector &A_o, TVector &sigma_o)
	{
		using T = Value_type<TVector>;

		auto image = mod_gaussian_deconv(stream, fft2, grid, sigma, image_i);

		find_fast_peaks_2d(stream, fft2, grid, image, sigma, thres, x_o, y_o);

		//lapack::GELS<T> gels;

		T d_error = 1e-4*grid.dR_min();

		A_o.resize(x_o.size());
		sigma_o.resize(x_o.size());
		for(auto ipk=0; ipk<x_o.size(); ipk++)
		{
			r2d<T> p_i(x_o[ipk], y_o[ipk]);
			//T radius_n = neighbor_radius(grid, image_i, p_i, 5*sigma);
			auto theta = gaussian_fit(grid, image_i, p_i, sigma, niter, d_error);
			A_o[ipk] = theta[0];
			sigma_o[ipk] = theta[2];
			//sigma_o[ipk] = (theta[2]>radius_n)?radius_n:theta[2];
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
	r2d<Value_type<TVector>> max_pos(Stream<e_host> &stream, TGrid &grid, TVector &Im, 
	r2d<Value_type<TVector>> p_i, Value_type<TVector> radius_i)
	{
		using T = Value_type<TVector>;

		auto get_rmin = [&](TVector &Im, r2d<T> p_i, T radius_i)->T
		{
			// get radial distribution
			TVector rl;
			TVector frl;
			radial_distribution_2d(grid, Im, p_i.x, p_i.y, radius_i, rl, frl);

			// log
			thrust::transform(frl.begin(), frl.end(), frl.begin(), [](const T &x){return log(x); });

			// smooth
			int nkr = 1;
			frl = smooth(stream, frl, nkr);

			// derivative and minimun
			int idr = 0;
			T fr_min = frl[1]-frl[0];
			for(auto ir= 0; ir<frl.size()-1; ir++)
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
		auto rmin = get_rmin(Im, p_i, radius_i);
		auto p_o = Rx_Ry_weight(grid, Im, p_i, rmin, false);

		return p_o;
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
		vector<complex<Value_type<TVector>>> Im_c(Im.begin(), Im.end());

		shift_2d(stream, fft2, grid, p, Im_c);
		assign_real(stream, Im_c, Im);
	}

	// phase correlation function
	template <class TGrid, class TVector>
	TVector PCF(Stream<e_host> &stream, FFT2<Value_type<TGrid>, e_host> &fft2, TGrid &grid, 
	TVector &M1_i, TVector &M2_i, Value_type<TGrid> k, Value_type<TGrid> sigma)
	{
		using T = Value_type<TVector>;
		using TVector_c = Vector<complex<T>, e_host>;

		// apply Hanning filter and copy data to complex matrix
		TVector_c M_1(grid.nxy());
		TVector_c M_2(grid.nxy());

		auto thr_diff_x_Hanning = [&](const Range &range)
		{
			T cx = c_2Pi/grid.nx;
			T cy = c_2Pi/grid.ny;

			for(auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for(auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					int ixy = grid.ind_col(ix, iy);
					int ix_n = (ix+1<grid.nx)?(ix+1):ix;
					int ixy_n = grid.ind_col(ix_n, iy);
					T fxy = 0.25*(1.0-cos(cx*grid.Rx(ix)))*(1.0-cos(cy*grid.Ry(iy)));
					fxy = (fxy>k)?1.0:fxy/k;
					M_1[ixy] = complex<T>((M1_i[ixy_n]-M1_i[ixy])*fxy);
					M_2[ixy] = complex<T>((M2_i[ixy_n]-M2_i[ixy])*fxy);
				}

			}
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_diff_x_Hanning);

		// create 2d plan
		fft2.create_plan_2d(grid.ny, grid.nx, stream.size());

		// shift matrix
		fft2_shift(stream, grid, M_1);
		fft2_shift(stream, grid, M_2);

		// fft2
		fft2.forward(M_1);
		fft2.forward(M_2);

		// apply symmetric bandwidth
		sigma = sigma*grid.dg_min();

		auto thr_PCF = [&](const Range &range)
		{
			T alpha = 0.5/(sigma*sigma);

			for(auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for(auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					int ixy = grid.ind_col(ix, iy);
					complex<T> z = conj(M_1[ixy])*M_2[ixy];
					auto mz = abs(z);
					M_2[ixy] = (isZero(mz))?0:z*exp(-alpha*grid.g2_shift(ix, iy))/mz;
				}
			}
		};

		stream.set_n_act_stream(grid.nx);
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(thr_PCF);

		fft2.inverse(M_2);

		TVector M_o(grid.nxy());
		assign_real(stream, M_2, M_o);

		// shift M
		fft2_shift(stream, grid, M_o);

		T m_min = min_element(stream, M_o)-1;
		// subtract minimum value
		std::for_each(M_o.begin(), M_o.end(), [&](T &v){ v -= m_min; });

		return M_o;
	}

	// get shift
	template <class TGrid, class TVector>
	r2d<Value_type<TGrid>> find_shift_2d(Stream<e_host> &stream, FFT2<Value_type<TGrid>, e_host> &fft2, TGrid &grid, 
	TVector &M1_i, TVector &M2_i, Value_type<TGrid> k, Value_type<TGrid> sigma)
	{
		using T = Value_type<TVector>;

		auto M = PCF(stream, fft2, grid, M1_i, M2_i, k, sigma);

		//get maximum index position
		int ixy_max = std::max_element(M.begin(), M.end())-M.begin();

		int ix_max, iy_max;
		grid.row_col(ixy_max, ix_max, iy_max);

		T radius = std::min(grid.lxh(), grid.lyh());
		r2d<T> r(grid.Rx(ix_max), grid.Ry(iy_max));

		r = max_pos(stream, grid, M, r, radius);
		r -= r2d<T>(grid.lxh(), grid.lyh());

		return r;
	}

	// correct shift
	template <class TGrid, class TVector>
	r2d<Value_type<TVector>> correct_shift_2d(Stream<e_host> &stream, FFT2<Value_type<TGrid>, e_host> &fft2, TGrid &grid, 
	TVector &M1_i, TVector &M2_io, Value_type<TGrid> k, Value_type<TGrid> sigma)
	{
		using T = Value_type<TVector>;

		vector<complex<T>> M_c(M1_i.size());
		int n_it = 3;

		r2d<T> dr(0,0);

		for (auto it=1; it<n_it; it++)
		{
			auto dr_t = find_shift_2d(stream, fft2, grid, M1_i, M2_io, k, sigma);
			assign(stream, M2_io, M_c);
			shift_2d(stream, fft2, grid, -dr_t, M_c);
			assign_real(stream, M_c, M2_io);
			dr += dr_t;
		}

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