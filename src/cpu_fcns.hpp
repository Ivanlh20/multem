/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * Multem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version of the License, or
 * (at your option) any later version.
 *
 * Multem is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Multem. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef CPU_FCNS_H
	#define CPU_FCNS_H

	#include <thread>
	#include <type_traits>
	#include <tuple>
	#include <algorithm>
	#include <numeric>
	#include <fstream>

	#include "const_enum.cuh"
	#include "math.cuh"
	#include "type_traits_gen.cuh"
	#include "cgpu_vctr.cuh"
	#include "kahan_sum.cuh"
	#include "cgpu_stream.cuh"
	#include "cgpu_fft.cuh"
	#include "cgpu_fcns_gen.cuh"
	#include "grid.cuh"
	#include "cgpu_detail.cuh"
	#include "cpu_detail.hpp"
	#include "border.cuh"
	#include "region.cuh"
	#include "types.cuh"
	#include "fcns_wd.cuh"

	/* vector functions */
	namespace mt
	{
		/***************************************************************************************/
		template <class TVctr_c, class TVctr_r>
		enable_if_cvctr_cpu_and_rvctr_cpu<TVctr_c, TVctr_r, void>
		fcn_assign_real(TVctr_c& mx_i, TVctr_r& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_r>;
			auto thr_assign_real = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, mx_o.begin(), cgpu_fctr::assign_real<T>());
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_assign_real);
		}

		template <class TVctr_c, class TVctr_r>
		enable_if_cvctr_cpu_and_rvctr_cpu<TVctr_c, TVctr_r, void>
		fcn_assign_max_real(TVctr_c& mx_i, Value_type<TVctr_r> M_v, TVctr_r& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_r>;
			auto thr_assign_max_real = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, mx_o.begin() + ithread.ind_0, cgpu_fctr::assign_max_real<T>(M_v, M_v));
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_assign_max_real);
		}

		template <class TVctr_c, class TVctr_r>
		enable_if_cvctr_and_rvctr<TVctr_c, TVctr_r, void>
		fcn_assign_real(TVctr_c& mx_i, Value_type<TVctr_r> M_v, TVctr_r& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_r>;

			if (M_v>0)
			{
				fcn_assign_max_real(mx_i, T(0), mx_o, pstream);
			}
			else
			{
				fcn_assign_real(mx_i, mx_o, pstream);
			}
		}

		template <class TVctr_c, class TVctr_r>
		enable_if_cvctr_cpu_and_rvctr_cpu<TVctr_c, TVctr_r, void>
		fcn_assign_abs_real(TVctr_c& mx_i, TVctr_r& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_r>;
			auto thr_assign_abs_real = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, mx_o.begin() + ithread.ind_0, cgpu_fctr::assign_abs_real<T>());
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_assign_abs_real);
		}

		template <class TVctr_c, class TVctr_r>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_c, TVctr_r, void>
		fcn_assign_abs(TVctr_c& mx_i, TVctr_r& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_r>;
			auto thr_assign_abs = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, mx_o.begin() + ithread.ind_0, cgpu_fctr::assign_abs<T>());
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_assign_abs);
		}

		/***************************************************************************************/
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_fill(TVctr& mx_io, Value_type<TVctr> val, Stream_cpu* pstream = nullptr)
		{
			auto thr_fill = [&](const iThread_Rect_1d& ithread)
			{
				thrust::fill(mx_io.begin() + ithread.ind_0, mx_io.begin() + ithread.ind_e, val);
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_io.size_32(), thr_fill);
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_scale(Value_type<TVctr_2> w_i, TVctr_1& mx_i, TVctr_2& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			auto thr_scale = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, 
					mx_o.begin() + ithread.ind_0, cgpu_fctr::scale<T>(w_i));
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_scale);
		}

		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_scale(const Value_type<TVctr>& w_i, TVctr& mx_io, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			auto thr_scale = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_io.begin() + ithread.ind_0, mx_io.begin() + ithread.ind_e, 
					mx_io.begin() + ithread.ind_0, cgpu_fctr::scale<T>(w_i));
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_io.size_32(), thr_scale);
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_norm_2(TVctr_1& mx_i, TVctr_2& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			auto thr_norm_2 = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, 
					mx_o.begin() + ithread.ind_0, cgpu_fctr::norm_2<T>());
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_norm_2);
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_scale_norm_2(Value_type<TVctr_2> w_i, TVctr_1& mx_i, TVctr_2& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			auto thr_scale_norm_2 = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, 
					mx_o.begin() + ithread.ind_0, cgpu_fctr::scale_norm_2<T>(w_i));
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_scale_norm_2);
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_sub(TVctr_1&& mx_1i, TVctr_1&& mx_2i, TVctr_2&& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			auto thr_add = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_1i.begin() + ithread.ind_0, mx_1i.begin() + ithread.ind_e, 
					mx_2i.begin() + ithread.ind_0, mx_o.begin() + ithread.ind_0, cgpu_fctr::sub<T>());
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_add);
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_add(TVctr_1&& mx_1i, TVctr_1&& mx_2i, TVctr_2&& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			auto thr_sub = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_1i.begin() + ithread.ind_0, mx_1i.begin() + ithread.ind_e, 
					mx_2i.begin() + ithread.ind_0, mx_o.begin() + ithread.ind_0, cgpu_fctr::add<T>());
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_sub);
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_add(TVctr_1&& mx_i, TVctr_2&& mx_io, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			auto thr_add = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, 
					mx_io.begin() + ithread.ind_0, mx_io.begin() + ithread.ind_0, cgpu_fctr::add<T>());
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_io.size_32(), thr_add);
		}


		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_add_scale(Value_type<TVctr_1> w1_i, TVctr_1& mx_1i, 
		Value_type<TVctr_1> w2_i, TVctr_1& mx_2i, TVctr_2& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			auto thr_add_scale = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_1i.begin() + ithread.ind_0, mx_1i.begin() + ithread.ind_e, 
					mx_2i.begin() + ithread.ind_0, mx_o.begin() + ithread.ind_0, cgpu_fctr::add_scale_i<T>(w1_i, w2_i));
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_add_scale);
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_add_scale(Value_type<TVctr_1> w_i, TVctr_1& mx_i, TVctr_2& mx_io, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			auto thr_add_scale = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, 
					mx_io.begin() + ithread.ind_0, mx_io.begin() + ithread.ind_0, cgpu_fctr::add_scale<T>(w_i));
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_io.size_32(), thr_add_scale);
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_add_norm_2(TVctr_1& mx_1i, TVctr_1& mx_2i, TVctr_2& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			auto thr_add_norm_2 = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_1i.begin() + ithread.ind_0, mx_1i.begin() + ithread.ind_e, 
					mx_2i.begin() + ithread.ind_0, mx_o.begin() + ithread.ind_0, cgpu_fctr::add_norm_2_i<T>());
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_add_norm_2);
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_add_norm_2(TVctr_1& mx_i, TVctr_2& mx_io, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			auto thr_add_norm_2 = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, 
					mx_io.begin() + ithread.ind_0, mx_io.begin() + ithread.ind_0, cgpu_fctr::add_norm_2<T>());
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_io.size_32(), thr_add_norm_2);
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_add_scale_norm_2(Value_type<TVctr_2> w1_i, TVctr_1& mx_1i, Value_type<TVctr_2> w2_i, TVctr_1& mx_2i, TVctr_2& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			auto thr_add_scale_norm_2 = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_1i.begin() + ithread.ind_0, mx_1i.begin() + ithread.ind_e, 
					mx_2i.begin() + ithread.ind_0, mx_o.begin() + ithread.ind_0, cgpu_fctr::add_scale_norm_2_i<T>(w1_i, w2_i));
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_add_scale_norm_2);
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_add_scale_norm_2(Value_type<TVctr_2> w_i, TVctr_1& mx_i, TVctr_2& mx_io, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			auto thr_add_scale_norm_2 = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, 
					mx_io.begin() + ithread.ind_0, mx_io.begin() + ithread.ind_0, cgpu_fctr::add_scale_norm_2<T>(w_i));
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_io.size_32(), thr_add_scale_norm_2);
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_mult(TVctr_1& mx_1i, TVctr_1& mx_2i, TVctr_2& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			auto thr_mult = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_1i.begin() + ithread.ind_0, mx_1i.begin() + ithread.ind_e, 
					mx_2i.begin() + ithread.ind_0, mx_o.begin() + ithread.ind_0, cgpu_fctr::mult<T>());
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_mult);
		}

		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_mult(TVctr_1& mx_i, TVctr_2& mx_io, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_2>;
			auto thr_mult = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, 
					mx_io.begin() + ithread.ind_0, mx_io.begin() + ithread.ind_0, cgpu_fctr::mult<T>());
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_io.size_32(), thr_mult);
		}

		/***************************************************************************************/
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_div_sft(Value_type<TVctr> mx_sft, Value_type<TVctr> mx_div, TVctr& mx_io, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			auto thr_div_sft = [&](const iThread_Rect_1d& ithread)
			{
				thrust::transform(mx_io.begin() + ithread.ind_0, mx_io.begin() + ithread.ind_e, mx_io.begin() + ithread.ind_0, cgpu_fctr::div_sft<T>(mx_sft, mx_div));;
			};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_io.size_32(), thr_div_sft);
		}

		/***************************************************************************************/
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type<TVctr>>
		fcn_sum(TVctr& mx_i, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			KS<T> sum_total = T(0);

			if (pstream == nullptr)
			{
				sum_total = thrust::reduce(mx_i.begin(), mx_i.end());
			}
			else
			{
				auto thr_sum = [&](const iThread_Rect_1d& ithread)
				{
					auto sum_partial = thrust::reduce(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e);

					pstream->stream_mutex.lock();
					sum_total += sum_partial;
					pstream->stream_mutex.unlock();
				};

				pstream->exec_xd_fcn<edim_1>(mx_i.size_32(), thr_sum);
			}

			return sum_total;
		}

		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_sum_norm_2(TVctr& mx_i, Stream_cpu* pstream = nullptr)
		{
			using T_r = Value_type_r<TVctr>;

			KS<T_r> sum_total = T_r(0);

			if (pstream == nullptr)
			{
				sum_total = thrust::transform_reduce(mx_i.begin(), mx_i.end(), cgpu_fctr::norm_2<T_r>(), T_r(0), cgpu_fctr::add<T_r>());
			}
			else
			{
				auto thr_sum_norm_2 = [&](const iThread_Rect_1d& ithread)
				{
					auto sum_partial = thrust::transform_reduce(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, 
						cgpu_fctr::norm_2<T_r>(), T_r(0), cgpu_fctr::add<T_r>());

					pstream->stream_mutex.lock();
					sum_total += sum_partial;
					pstream->stream_mutex.unlock();
				};

				pstream->exec_xd_fcn<edim_1>(mx_i.size_32(), thr_sum_norm_2);
			}

			return sum_total;
		}

		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_sum_norm_2_sft(TVctr& mx_i, Value_type<TVctr> x_sft, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			using T_r = Value_type_r<TVctr>;

			KS<T_r> sum_total = T_r(0);

			if (pstream == nullptr)
			{
				sum_total = thrust::transform_reduce(mx_i.begin(), mx_i.end(), cgpu_fctr::norm_2_sft<T, T_r>(x_sft), T_r(0), cgpu_fctr::add<T_r>());
			}
			else
			{
				auto thr_sum_norm_2_sft = [&](const iThread_Rect_1d& ithread)
				{
					auto sum_partial = thrust::transform_reduce(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, 
						cgpu_fctr::norm_2_sft<T, T_r>(x_sft), T_r(0), cgpu_fctr::add<T_r>());

					pstream->stream_mutex.lock();
					sum_total += sum_partial;
					pstream->stream_mutex.unlock();
				};

				pstream->exec_xd_fcn<edim_1>(mx_i.size_32(), thr_sum_norm_2_sft);
			}

			return sum_total;
		}

 		template <class TVctr>
		enable_if_cvctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_sum_max_real(TVctr& mx_i, Value_type_r<TVctr> v_min, Stream_cpu* pstream = nullptr)
		{
			using T_r = Value_type_r<TVctr>;

			KS<T_r> sum_total = T_r(0);

			if (pstream == nullptr)
			{
				sum_total = thrust::transform_reduce(mx_i.begin(), mx_i.end(), cgpu_fctr::assign_max_real<T_r>(v_min, T_r(0)), T_r(0), cgpu_fctr::add<T_r>());
			}
			else
			{
				auto thr_sum_max_real = [&](const iThread_Rect_1d& ithread)
				{
					auto sum_partial = thrust::transform_reduce(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, 
						cgpu_fctr::assign_max_real<T_r>(v_min, T_r(0)), T_r(0), cgpu_fctr::add<T_r>());

					pstream->stream_mutex.lock();
					sum_total += sum_partial;
					pstream->stream_mutex.unlock();
				};

				pstream->exec_xd_fcn<edim_1>(mx_i.size_32(), thr_sum_max_real);
			}

			return sum_total;
		}

 		template <class TVctr>
		enable_if_cvctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_sum_abs_real(TVctr& mx_i, Stream_cpu* pstream = nullptr)
		{
			using T_r = Value_type_r<TVctr>;

			KS<T_r> sum_total = T_r(0);

			if (pstream == nullptr)
			{
				sum_total = thrust::transform_reduce(mx_i.begin(), mx_i.end(), cgpu_fctr::assign_abs_real<T_r>(), T_r(0), cgpu_fctr::add<T_r>());
			}
			else
			{
				auto thr_sum_abs_real = [&](const iThread_Rect_1d& ithread)
				{
					auto sum_partial = thrust::transform_reduce(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, 
						cgpu_fctr::assign_abs_real<T_r>(), T_r(0), cgpu_fctr::add<T_r>());

					pstream->stream_mutex.lock();
					sum_total += sum_partial;
					pstream->stream_mutex.unlock();
				};

				pstream->exec_xd_fcn<edim_1>(mx_i.size_32(), thr_sum_abs_real);
			}

			return sum_total;
		}

 		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_sum_abs(TVctr& mx_i, Stream_cpu* pstream = nullptr)
		{
			using T_r = Value_type_r<TVctr>;

			KS<T_r> sum_total = T_r(0);

			if (pstream == nullptr)
			{
				sum_total = thrust::transform_reduce(mx_i.begin(), mx_i.end(), cgpu_fctr::assign_abs<T_r>(), T_r(0), cgpu_fctr::add<T_r>());
			}
			else
			{
				auto thr_sum_abs = [&](const iThread_Rect_1d& ithread)
				{
					auto sum_partial = thrust::transform_reduce(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, 
						cgpu_fctr::assign_abs<T_r>(), T_r(0), cgpu_fctr::add<T_r>());

					pstream->stream_mutex.lock();
					sum_total += sum_partial;
					pstream->stream_mutex.unlock();
				};

				pstream->exec_xd_fcn<edim_1>(mx_i.size_32(), thr_sum_abs);
			}

			return sum_total;
		}

 		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_sum_abs_sft(TVctr& mx_i, Value_type<TVctr> x_sft, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			using T_r = Value_type_r<TVctr>;

			KS<T_r> sum_total = T_r(0);

			if (pstream == nullptr)
			{
				sum_total = thrust::transform_reduce(mx_i.begin(), mx_i.end(), cgpu_fctr::abs_sft<T, T_r>(), T_r(0), cgpu_fctr::add<T_r>());
			}
			else
			{
				auto thr_sum_abs_sft = [&](const iThread_Rect_1d& ithread)
				{
					auto sum_partial = thrust::transform_reduce(mx_i.begin() + ithread.ind_0, mx_i.begin() + ithread.ind_e, 
						cgpu_fctr::abs_sft<T, T_r>(), T_r(0), cgpu_fctr::add<T_r>());

					pstream->stream_mutex.lock();
					sum_total += sum_partial;
					pstream->stream_mutex.unlock();
				};

				pstream->exec_xd_fcn<edim_1>(mx_i.size_32(), thr_sum_abs_sft);
			}

			return sum_total;
		}

		/***************************************************************************************/
		/* calculate mean */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type<TVctr>>
		fcn_mean(TVctr& mx_i, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			return fcn_sum(mx_i, pstream)/T(mx_i.size());
		}

		/* calculate mean of the norm square */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_mean_norm_2(TVctr& mx_i, Stream_cpu* pstream = nullptr)
		{
			using T_r = Value_type_r<TVctr>;
			return fcn_sum_norm_2(mx_i, pstream)/T_r(mx_i.size());
		}

		/* calculate mean of the shifted norm square */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_mean_norm_2_sft(TVctr& mx_i, Value_type<TVctr> x_sft, Stream_cpu* pstream = nullptr)
		{
			using T_r = Value_type_r<TVctr>;
			return fcn_sum_norm_2_sft(mx_i, x_sft, pstream)/T_r(mx_i.size());
		}

 		template <class TVctr>
		enable_if_cvctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_mean_max_real(TVctr& mx_i, Value_type_r<TVctr> v_min, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			return fcn_sum_max_real(mx_i, v_min, pstream)/T(mx_i.size());
		}

 		template <class TVctr>
		enable_if_cvctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_mean_abs_real(TVctr& mx_i, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			return fcn_sum_abs_real(mx_i, pstream)/T(mx_i.size());
		}

 		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_mean_abs(TVctr& mx_i, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			return fcn_sum_abs(mx_i, pstream)/T(mx_i.size());
		}

 		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_mean_abs_sft(TVctr& mx_i, Value_type<TVctr> x_sft, Stream_cpu* pstream = nullptr)
		{
			using T_r = Value_type_r<TVctr>;
			return fcn_sum_abs_sft(mx_i, x_sft, pstream)/T_r(mx_i.size());
		}

		/***************************************************************************************/
		/* calculate mean and variance */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_mean_var(TVctr& mx_i, Value_type<TVctr>& x_mean, Value_type_r<TVctr>& x_var, Stream_cpu* pstream = nullptr)
		{
			x_mean = fcn_mean(mx_i, pstream);
			x_var = fcn_mean_norm_2_sft(mx_i, x_mean, pstream);
		}

		/* calculate mean and standard deviation */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_mean_std(TVctr& mx_i, Value_type<TVctr>& x_mean, Value_type_r<TVctr>& x_std, Stream_cpu* pstream = nullptr)
		{
			fcn_mean_var(mx_i, x_mean, x_std, pstream);
			x_std = ::sqrt(x_std);
		}

		// mean and first absolute moment
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_mean_moment_1a(TVctr& mx_i, Value_type<TVctr>& x_mean, Value_type_r<TVctr>& x_ma1, Stream_cpu* pstream = nullptr)
		{
			x_mean = fcn_mean(mx_i, pstream);
			x_ma1 = fcn_mean_abs_sft(mx_i, x_mean, pstream);
		}

		/* variance */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_variance(TVctr& mx_i, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			using T_r = Value_type_r<TVctr>;

			T x_mean;
			T_r x_var;
			fcn_mean_var(mx_i, x_mean, x_var, pstream);

			return x_var;
		}

		/* standard deviation */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_std(TVctr& mx_i, Stream_cpu* pstream = nullptr)
		{
			auto x_var = fcn_variance(mx_i, pstream);

			return ::sqrt(x_var);
		}

		// first absolute moment
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_moment_1a(TVctr& mx_i, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;
			using T_r = Value_type_r<TVctr>;

			T x_mean;
			T_r x_ma1;
			fcn_mean_moment_1a(mx_i, x_mean, x_ma1, pstream);

			return x_ma1;
		}

		/***************************************************************************************/
		/* minimum element of a vector */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type<TVctr>>
		fcn_min_element(const TVctr& x, Stream_cpu* pstream = nullptr)
		{
			auto x_min = x[0];			
			
			if (pstream == nullptr)
			{
				x_min = *thrust::min_element(x.begin(), x.end());
			}
			else
			{
				auto thr_min = [&](const iThread_Rect_1d& ithread)
				{
					auto x_min_partial = *thrust::min_element(x.begin() + ithread.ind_0, x.begin() + ithread.ind_e);

					pstream->stream_mutex.lock();
					x_min = fcn_min(x_min, x_min_partial);
					pstream->stream_mutex.unlock();
				};

				pstream->exec_xd_fcn<edim_1>(x.size_32(), thr_min);
			}

			return x_min;
		}

		/* maximum element of a vector */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type<TVctr>>
		fcn_max_element(const TVctr& x, Stream_cpu* pstream = nullptr)
		{
			auto x_max = x[0];

			if (pstream == nullptr)
			{
				x_max = *thrust::max_element(x.begin(), x.end());
			}
			else
			{
				auto thr_max = [&](const iThread_Rect_1d& ithread)
				{
					auto x_max_partial = *thrust::max_element(x.begin() + ithread.ind_0, x.begin() + ithread.ind_e);

					pstream->stream_mutex.lock();
					x_max = fcn_max(x_max, x_max_partial);
					pstream->stream_mutex.unlock();
				};

				pstream->exec_xd_fcn<edim_1>(x.size_32(), thr_max);
			}

			return x_max;
		}

		/* minimum and maximum element of a vector */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_minmax_element(const TVctr& x, Value_type<TVctr>& x_min, Value_type<TVctr>& x_max, Stream_cpu* pstream = nullptr)
		{
			x_min = x[0];
			x_max = x[0];

			if (pstream == nullptr)
			{
				auto x_min_max = thrust::minmax_element(x.begin(), x.end());
				x_min = *(x_min_max.first);
				x_max = *(x_min_max.second);
			}
			else
			{
				auto thr_min_max = [&](const iThread_Rect_1d& ithread)
				{
					auto x_min_max_partial = thrust::minmax_element(x.begin() + ithread.ind_0, x.begin() + ithread.ind_e);

					pstream->stream_mutex.lock();
					x_min = fcn_min(x_min, *(x_min_max_partial.first));
					x_max = fcn_max(x_max, *(x_min_max_partial.second));
					pstream->stream_mutex.unlock();
				};

				pstream->exec_xd_fcn<edim_1>(x.size_32(), thr_min_max);
			}
		}
	}

	/* add - assign - crop - norm_2 - fftsft */
	namespace mt
	{
		/* shift matrix respect to nx_h */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_fftsft_1d(TVctr& mx_io, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_io.igrid_1d();

			fcn_stream_exec_xd_krn<edim_1>(pstream, igrid.nx_h, cgpu_detail::fcn_fftsft_1d<T>, igrid, mx_io.m_data);
		}		
		
		/* shift matrix respect to ny_h */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_fftsft_bc_2d(TVctr& mx_io, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_io.igrid_2d();

			fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, igrid.ny_h, cgpu_detail::fcn_fftsft_bc_2d<T>, igrid, mx_io.m_data);
		}

		/* shift matrix respect to (nx_h, ny_h) */
		template <class TVctr>
		//enable_if_vctr_cpu<TVctr, void>
		void fcn_fftsft_2d(TVctr& mx_io, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_io.igrid_2d();

			// function fcn_fftsft_2d is overloading
 			void (*pfcn_fftsft_2d)(const dt_int32&, const dt_int32&, const iGrid_2d&, T*) = &cgpu_detail::fcn_fftsft_2d<T>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx_h, igrid.ny_h, pfcn_fftsft_2d, igrid, mx_io.m_data);
		}
		
		/* shift matrix respect to (nx_h, ny_h) */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_fftsft_2d(TVctr& mx_i, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			// function fcn_fftsft_2d is overloading
 			void (*pfcn_fftsft_2d)(const dt_int32&, const dt_int32&, const iGrid_2d&, T*, T*) = &cgpu_detail::fcn_fftsft_2d<T>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx_h, igrid.ny_h, pfcn_fftsft_2d, igrid, mx_i.m_data, mx_o.m_data);
		}

		/* shift matrix respect to (nx_h, ny_h, nz_h) */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_fftsft_3d(TVctr& mx_io, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_io.igrid_3d();

			// function fcn_fftsft_3d is overloading
 			void (*pfcn_fftsft_3d)(const dt_int32&, const dt_int32&, const iGrid_3d&, T*) = &cgpu_detail::fcn_fftsft_3d<T>;

			fcn_stream_exec_xd_krn<edim_3>(pstream, igrid.nx_h, igrid.ny_h, igrid.nz_h, pfcn_fftsft_3d, igrid, mx_io.m_data);
		}
		
		/* shift matrix respect to (nx_h, ny_h, nz_h) */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_fftsft_3d(TVctr& mx_i, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_3d();

			// function fcn_fftsft_3d is overloading
 			void (*pfcn_fftsft_3d)(const dt_int32&, const dt_int32&, const iGrid_3d&, T*, T*) = &cgpu_detail::fcn_fftsft_3d<T>;

			fcn_stream_exec_xd_krn<edim_3>(pstream, igrid.nx_h, igrid.ny_h, igrid.nz_h, pfcn_fftsft_3d, igrid, mx_i.m_data, mx_o.m_data);
		}

		/***************************************************************************************/
 		/* add, scale and shift */
 		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_add_sc_fftsft_2d(TVctr& mx_i, Value_type<TVctr> w, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx_h, igrid.ny_h, cgpu_detail::fcn_add_sc_fftsft_2d<T>, igrid, mx_i.m_data, w, mx_o.m_data);
		}
 
		/* add, scale, square and shift */
 		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_add_sc_norm_2_fftsft_2d(TVctr_1& mx_i, Value_type<TVctr_2> w, TVctr_2& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_1>;
			using U = Value_type<TVctr_2>;

			auto igrid = mx_i.igrid_2d();

			fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx_h, igrid.ny_h, cgpu_detail::fcn_add_sc_norm_2_fftsft_2d<T>, igrid, mx_i.m_data, w, mx_o.m_data);
		}

		/***************************************************************************************/
		/* assign and crop */
 		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_assign_crop_2d(TVctr& mx_i, TVctr& mx_o, iRegion_Rect_2d& iregion, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, igrid.ny, cgpu_detail::fcn_assign_crop_2d<T>, igrid, mx_i.m_data, iregion, mx_o.m_data);
		}

		/* assign, crop and shift */
 		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_assign_crop_fftsft_2d(TVctr& mx_i, iRegion_Rect_2d& iregion, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx_h, igrid.ny_h, cgpu_detail::fcn_assign_crop_fftsft_2d<T>, igrid, mx_i.m_data, iregion, mx_o.m_data);
		}

 		/* add, scale, crop and shift */
 		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_add_sc_crop_fftsft_2d(TVctr& mx_i, iRegion_Rect_2d& iregion, Value_type<TVctr> w, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx_h, igrid.ny_h,cgpu_detail::fcn_add_sc_crop_fftsft_2d<T>, igrid, mx_i.m_data, iregion, w, mx_o.m_data);
		}

		/* add, scale, square, crop and shift */
 		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_add_sc_norm_2_crop_fftsft_2d(TVctr_1& mx_i, iRegion_Rect_2d& iregion, Value_type<TVctr_2> w, TVctr_2& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr_1>;
			using U = Value_type<TVctr_2>;

			auto igrid = mx_i.igrid_2d();

			fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx_h, igrid.ny_h, cgpu_detail::fcn_add_sc_norm_2_crop_fftsft_2d<T, U>, igrid, mx_i.m_data, iregion, w, mx_o.m_data);
		}
	}

	/* transpose - element wise matrix op vector */
	namespace mt
	{
		/* transpose */
 		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_trs_2d(TVctr& mx_i, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			Vctr_cpu<T> mx_t(mx_i.shape_2d_trs());
			fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, igrid.ny, cgpu_detail::fcn_trs_2d<T>, igrid, mx_i.m_data, mx_t.m_data);

			return mx_t;
		}

		/* transpose */
 		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_trs_2d(TVctr& mx_i, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			if (&mx_i != &mx_o)
			{
				fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, igrid.ny, cgpu_detail::fcn_trs_2d<T>, igrid, mx_i.m_data, mx_o.m_data);
			}
			else
			{
				Vctr_cpu<T> mx_t(mx_i.shape_2d_trs());
				fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, igrid.ny, cgpu_detail::fcn_trs_2d<T>, igrid, mx_i.m_data, mx_t.m_data);
				mx_t.cpy_to_cpu_ptr(mx_o.m_data, mx_o.size());
			}
			
		}

		/* element wise addition: matrix + vector*/
 		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_ew_add_mx_vctr(TVctr& vctr, TVctr& mx_io, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_io.igrid_2d();

			if ((vctr.m_s0=mx_io.m_s0) && (vctr.m_s1==1))
			{
				fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, igrid.ny, cgpu_detail::fcn_ew_add_mx_vctr_col<T>, igrid, vctr.m_data, mx_io.m_data);
			}
			else if ((vctr.m_s0==1) && (vctr.m_s1==mx_io.m_s1) )
			{
				fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, igrid.ny, cgpu_detail::fcn_ew_add_mx_vctr_row<T>, igrid, vctr.m_data, mx_io.m_data);
			}
		}

		/* element wise subtraction: matrix - vector */
 		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_ew_sub_mx_vctr(TVctr& vctr, TVctr& mx_io, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_io.igrid_2d();

			if ((vctr.m_s0=mx_io.m_s0) && (vctr.m_s1==1))
			{
				fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, igrid.ny, cgpu_detail::fcn_ew_sub_mx_vctr_col<T>, igrid, vctr.m_data, mx_io.m_data);
			}
			else if ((vctr.m_s0==1) && (vctr.m_s1==mx_io.m_s1) )
			{
				fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, igrid.ny, cgpu_detail::fcn_ew_sub_mx_vctr_row<T>, igrid, vctr.m_data, mx_io.m_data);
			}
		}

		/* element wise multiplication matrix X vector */
 		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_ew_mult_mx_vctr(TVctr& vctr, TVctr& mx_io, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_io.igrid_2d();

			if ((vctr.m_s0=mx_io.m_s0) && (vctr.m_s1==1))
			{
				fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, igrid.ny, cgpu_detail::fcn_ew_mult_mx_vctr_col<T>, igrid, vctr.m_data, mx_io.m_data);
			}
			else if ((vctr.m_s0==1) && (vctr.m_s1==mx_io.m_s1) )
			{
				fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, igrid.ny, cgpu_detail::fcn_ew_mult_mx_vctr_row<T>, igrid, vctr.m_data, mx_io.m_data);
			}
		}
	}

	/* aperture functions */
	namespace mt
	{
		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_fermi_aperture(Grid_2d<T>& grid, T g2_cut, T alpha, T w, TVctr& mx_io, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail::fcn_fermi_aperture<T, U>, grid, g2_cut, alpha, w, mx_io.m_data);
		}

		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_hard_aperture(Grid_2d<T>& grid, T g2_cut, T w, TVctr& mx_io, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail::fcn_hard_aperture<T, U>, grid, g2_cut, w, mx_io.m_data);
		}
	}

	/* phase shifts real space*/
	namespace mt
	{
		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_rs_exp_factor_1d(Grid_1d<T>& grid, TVctr& psi_i, T gx, T w, TVctr& psi_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			fcn_stream_exec_xd_krn<edim_1>(pstream, grid.nx, cgpu_detail::fcn_rs_exp_factor_1d<T, U>, grid, psi_i.m_data, gx, w, psi_o.m_data);
		}

		template <class T, class TVctr_r, class TVctr_c>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_r, TVctr_c, void>
		fcn_rs_exp_factor_2d_bc(Grid_2d<T>& grid, TVctr_c& psi_i, T alpha, TVctr_r &gy, T w, TVctr_c& psi_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail::fcn_rs_exp_factor_2d_bc<T, U>, grid, alpha, psi_i.m_data, gy.m_data, w, psi_o.m_data);
		}

		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_rs_exp_factor_2d(Grid_2d<T>& grid, TVctr& psi_i, R_2d<T> g, T w, TVctr& psi_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail::fcn_rs_exp_factor_2d<T, U>, grid, psi_i.m_data, g, w, psi_o.m_data);
		}

		template <class T, class TVctr_r, class TVctr_c>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_r, TVctr_c, void>
		fcn_rs_mul_exp_factor_2d(Grid_2d<T>& grid, TVctr_c& psi_i, TVctr_r& g, T w, TVctr_c& psi_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail::fcn_rs_mul_exp_factor_2d<T, U>, grid, psi_i.m_data, g.m_data, w, psi_o.m_data);
		}
	}

	/* phase shifts fourier space*/
	namespace mt
	{
		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_fs_exp_factor_1d(Grid_1d<T>& grid, TVctr& psi_i, T rx, T w, TVctr& psi_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			fcn_stream_exec_xd_krn<edim_1>(pstream, grid.nx, cgpu_detail::fcn_fs_exp_factor_1d<T, U>, grid, psi_i.m_data, rx, w, psi_o.m_data);
		}

		template <class T, class TVctr_r, class TVctr_c>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_r, TVctr_c, void>
		fcn_fs_exp_factor_2d_bc(Grid_2d<T>& grid, TVctr_c& psi_i, T alpha, TVctr_r &ry, T w, TVctr_c& psi_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail::fcn_fs_exp_factor_2d_bc<T, U>, grid, alpha, psi_i.m_data, ry.m_data, w, psi_o.m_data);
		}

		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_fs_exp_factor_2d(Grid_2d<T>& grid, TVctr& psi_i, R_2d<T> r, T w, TVctr& psi_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail::fcn_fs_exp_factor_2d<T, U>, grid, psi_i.m_data, r, w, psi_o.m_data);
		}

		template <class T, class TVctr_r, class TVctr_c>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_r, TVctr_c, void>
		fcn_fs_mul_exp_factor_2d(Grid_2d<T>& grid, TVctr_c& psi_i, TVctr_r& r, T w, TVctr_c& psi_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail::fcn_fs_mul_exp_factor_2d<T, U>, grid, psi_i.m_data, r.m_data, w, psi_o.m_data);
		}
	}

	/* gradient */
	namespace mt
	{
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_grad_x(TVctr& mx_i, TVctr& dmx_x, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, igrid.ny, cgpu_detail::fcn_grad_x<T>, igrid, mx_i.m_data, dmx_x.m_data);
		}

		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_grad_y(TVctr& mx_i, TVctr& dmx_y, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, igrid.ny, cgpu_detail::fcn_grad_y<T>, igrid, mx_i.m_data, dmx_y.m_data);
		}

		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_grad(TVctr& mx_i, TVctr& dmx_x, TVctr& dmx_y, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();

			fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, igrid.ny, cgpu_detail::fcn_grad<T>, igrid, mx_i.m_data, dmx_x.m_data, dmx_y.m_data);
		}
	}

	/* function multiplication fourier space */
	namespace mt
	{
		#define FCN_MULT_FS_FCN_CPU(FN, FCN, POW, DIM)																													\
		template <class T, class TVctr_c>																																\
		enable_if_vctr_cpu<TVctr_c, void>																																\
		fcn_mult_fs_##FN##_##DIM##d(Grid_##DIM##d<T>& grid, FCN<T>& fcn, T w, TVctr_c& mx_io, Stream_cpu* pstream = nullptr)											\
		{																																								\
			using U = Value_type<TVctr_c>;																																\
																																										\
			fcn_stream_exec_xd_krn<EDIM(DIM)>(pstream, F_SHAPE_ND(grid, DIM), cgpu_detail::fcn_mult_fs_fcn_g##POW##_##DIM##d<T, U, FCN<T>>, grid, fcn, w, mx_io.m_data);		\
		}

		/* gaussian */
		FCN_MULT_FS_FCN_CPU(gauss, Fcn_Gauss, 2, 1);		// fcn_mult_fs_gauss_1d
		FCN_MULT_FS_FCN_CPU(gauss, Fcn_Gauss, 2, 2);		// fcn_mult_fs_gauss_2d
		FCN_MULT_FS_FCN_CPU(gauss, Fcn_Gauss, 2, 3);		// fcn_mult_fs_gauss_3d

		/* exponential */
		FCN_MULT_FS_FCN_CPU(exp, Fcn_Exp, , 1);				// fcn_mult_fs_exp_1d
		FCN_MULT_FS_FCN_CPU(exp, Fcn_Exp, , 2);				// fcn_mult_fs_exp_2d
		FCN_MULT_FS_FCN_CPU(exp, Fcn_Exp, , 3);				// fcn_mult_fs_exp_3d

		/* fermi */
		FCN_MULT_FS_FCN_CPU(fermi, Fcn_Fermi, 2, 1);		// fcn_mult_fs_fermi_1d
		FCN_MULT_FS_FCN_CPU(fermi, Fcn_Fermi, 2, 2);		// fcn_mult_fs_fermi_2d
		FCN_MULT_FS_FCN_CPU(fermi, Fcn_Fermi, 2, 3);		// fcn_mult_fs_fermi_3d

		/* butterworth */
		FCN_MULT_FS_FCN_CPU(butwth, Fcn_Butwth, 2, 1);		// fcn_mult_fs_butwth_1d
		FCN_MULT_FS_FCN_CPU(butwth, Fcn_Butwth, 2, 2);		// fcn_mult_fs_butwth_2d
		FCN_MULT_FS_FCN_CPU(butwth, Fcn_Butwth, 2, 3);		// fcn_mult_fs_butwth_3d

		/* hann */
		FCN_MULT_FS_FCN_CPU(hann, Fcn_Hann, , 1);			// fcn_mult_fs_hann_1d
		FCN_MULT_FS_FCN_CPU(hann, Fcn_Hann, , 2);			// fcn_mult_fs_hann_2d
		FCN_MULT_FS_FCN_CPU(hann, Fcn_Hann, , 3);			// fcn_mult_fs_hann_3d
	}	

	/* convolution */
	namespace mt
	{
		#define FCN_CV_FS_FCN_CPU(FN, FCN, POW, DIM)																													\
		template <class T, class TVctr_c>																																\
		enable_if_cvctr_cpu<TVctr_c, void>																																\
		fcn_cv_fs_##FN##_##DIM##d(Grid_##DIM##d<T>& grid, FFT_cpu<T>& fft, FCN<T>& fcn, TVctr_c& mx_io, Stream_cpu* pstream = nullptr)									\
		{																																								\
			using U = Value_type<TVctr_c>;																																\
																																										\
			fcn_fftsft_##DIM##d(mx_io, pstream);																														\
			fft.forward(mx_io);																																			\
																																										\
			const T w = grid.isize_r();																																	\
																																										\
			fcn_stream_exec_xd_krn<EDIM(DIM)>(pstream, F_SHAPE_ND(grid, DIM), cgpu_detail::fcn_mult_fs_fcn_g##POW##_##DIM##d<T, U, FCN<T>>, grid, fcn, w, mx_io.m_data);		\
																																										\
			fft.inverse(mx_io);																																			\
			fcn_fftsft_##DIM##d(mx_io, pstream);																														\
		}

		/* gaussian convolution */
		FCN_CV_FS_FCN_CPU(gauss, Fcn_Gauss, 2, 1);			// fcn_cv_fs_gauss_1d
		FCN_CV_FS_FCN_CPU(gauss, Fcn_Gauss, 2, 2);			// fcn_cv_fs_gauss_2d
		FCN_CV_FS_FCN_CPU(gauss, Fcn_Gauss, 2, 3);			// fcn_cv_fs_gauss_3d

		/* exponential convolution */
		FCN_CV_FS_FCN_CPU(exp, Fcn_Exp, , 1);				// fcn_cv_fs_exp_1d
		FCN_CV_FS_FCN_CPU(exp, Fcn_Exp, , 2);				// fcn_cv_fs_exp_2d
		FCN_CV_FS_FCN_CPU(exp, Fcn_Exp, , 3);				// fcn_cv_fs_exp_3d

		/* fermi convolution */
		FCN_CV_FS_FCN_CPU(fermi, Fcn_Fermi, 2, 1);			// fcn_cv_fs_fermi_1d
		FCN_CV_FS_FCN_CPU(fermi, Fcn_Fermi, 2, 2);			// fcn_cv_fs_fermi_2d
		FCN_CV_FS_FCN_CPU(fermi, Fcn_Fermi, 2, 3);			// fcn_cv_fs_fermi_3d

		/* butterworth convolution */
		FCN_CV_FS_FCN_CPU(butwth, Fcn_Butwth, 2, 1);		// fcn_cv_fs_butwth_1d
		FCN_CV_FS_FCN_CPU(butwth, Fcn_Butwth, 2, 2);		// fcn_cv_fs_butwth_2d
		FCN_CV_FS_FCN_CPU(butwth, Fcn_Butwth, 2, 3);		// fcn_cv_fs_butwth_3d

		/* hann convolution */
		FCN_CV_FS_FCN_CPU(hann, Fcn_Hann, , 1);				// fcn_cv_fs_hann_1d
		FCN_CV_FS_FCN_CPU(hann, Fcn_Hann, , 2);				// fcn_cv_fs_hann_2d
		FCN_CV_FS_FCN_CPU(hann, Fcn_Hann, , 3);				// fcn_cv_fs_hann_3d
	}
	
	/* deconvolution */
	namespace mt
	{
		#define FCN_DCV_FS_FCN_CPU(FN, FCN, POW, DIM)																														\
		template <class T, class TVctr_c>																																	\
		enable_if_cvctr_cpu<TVctr_c, void>																																	\
		fcn_dcv_fs_##FN##_##DIM##d(Grid_##DIM##d<T>& grid, FFT_cpu<T>& fft, FCN<T>& fcn, T psnr, TVctr_c& mx_io, Stream_cpu* pstream = nullptr)								\
		{																																									\
			using U = Value_type<TVctr_c>;																																	\
																																											\
			fcn_fftsft_##DIM##d(mx_io, pstream);																															\
			fft.forward(mx_io);																																				\
																																											\
			const T w = grid.isize_r();																																		\
																																											\
			fcn_stream_exec_xd_krn<EDIM(DIM)>(pstream, F_SHAPE_ND(grid, DIM), cgpu_detail::fcn_dcv_fs_fcn_g##POW##_##DIM##d<T, U, FCN<T>>, grid, fcn, psnr, w, mx_io.m_data);	\
																																											\
			fft.inverse(mx_io);																																				\
			fcn_fftsft_##DIM##d(mx_io, pstream);																															\
		}

		/* gaussian convolution */
		FCN_DCV_FS_FCN_CPU(gauss, Fcn_Gauss, 2, 1);			// fcn_dcv_fs_gauss_1d
		FCN_DCV_FS_FCN_CPU(gauss, Fcn_Gauss, 2, 2);			// fcn_dcv_fs_gauss_2d
		FCN_DCV_FS_FCN_CPU(gauss, Fcn_Gauss, 2, 3);			// fcn_dcv_fs_gauss_3d

		/* exponential convolution */
		FCN_DCV_FS_FCN_CPU(exp, Fcn_Exp, , 1);				// fcn_dcv_fs_exp_1d
		FCN_DCV_FS_FCN_CPU(exp, Fcn_Exp, , 2);				// fcn_dcv_fs_exp_2d
		FCN_DCV_FS_FCN_CPU(exp, Fcn_Exp, , 3);				// fcn_dcv_fs_exp_3d

		/* fermi convolution */
		FCN_DCV_FS_FCN_CPU(fermi, Fcn_Fermi, 2, 1);			// fcn_dcv_fs_fermi_1d
		FCN_DCV_FS_FCN_CPU(fermi, Fcn_Fermi, 2, 2);			// fcn_dcv_fs_fermi_2d
		FCN_DCV_FS_FCN_CPU(fermi, Fcn_Fermi, 2, 3);			// fcn_dcv_fs_fermi_3d

		/* butterworth convolution */
		FCN_DCV_FS_FCN_CPU(butwth, Fcn_Butwth, 2, 1);		// fcn_dcv_fs_butwth_1d
		FCN_DCV_FS_FCN_CPU(butwth, Fcn_Butwth, 2, 2);		// fcn_dcv_fs_butwth_2d
		FCN_DCV_FS_FCN_CPU(butwth, Fcn_Butwth, 2, 3);		// fcn_dcv_fs_butwth_3d

		/* hann convolution */
		FCN_DCV_FS_FCN_CPU(hann, Fcn_Hann, , 1);			// fcn_dcv_fs_hann_1d
		FCN_DCV_FS_FCN_CPU(hann, Fcn_Hann, , 2);			// fcn_dcv_fs_hann_2d
		FCN_DCV_FS_FCN_CPU(hann, Fcn_Hann, , 3);			// fcn_dcv_fs_hann_3d
	}

	/* window functions */
	namespace mt
	{
		//template <eDim Dim, eFcn_typ Fcn_typ, class T, class TVctr>
		//enable_if_vctr_cpu<TVctr, void>
		//fcn_wd_fcn_xd(Grid_xd<T, Dim>& grid, Wd_fcn_xd<T, Dim, Fcn_typ>& fcn, dt_bool sft, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		//{
		//	using U = Value_type<TVctr>;
		//	
		//	if (sft)
		//	{
		//		fcn_stream_exec_xd_krn<Dim>(pstream, F_SHAPE_ND(grid, DIM), cgpu_detail::fcn_wd_fcn_r2_sft_xd<T, U, Fcn_typ>, grid, fcn, mx_o.m_data);
		//	}
		//	else
		//	{
		//		fcn_stream_exec_xd_krn<Dim>(pstream, F_SHAPE_ND(grid, DIM), cgpu_detail::fcn_wd_fcn_r2_xd<T, U, Fcn_typ>, grid, fcn, mx_o.m_data)
		//	}
		//}

		// FCN_WD_FCN_CPU(gauss, Wd_Gauss_1d, 1);		// fcn_wd_gauss_1d
		// FCN_WD_FCN_CPU(gauss, Wd_Gauss_2d, 2);		// fcn_wd_gauss_2d
		// FCN_WD_FCN_CPU(gauss, Wd_Gauss_3d, 3);		// fcn_wd_gauss_3d

		// FCN_WD_FCN_CPU(exp, Wd_Exp_1d, 1);			// fcn_wd_exp_1d
		// FCN_WD_FCN_CPU(exp, Wd_Exp_2d, 2);			// fcn_wd_exp_2d
		// FCN_WD_FCN_CPU(exp, Wd_Exp_3d, 3);			// fcn_wd_exp_3d

		// FCN_WD_FCN_CPU(fermi, Wd_Fermi_1d, 1);		// fcn_wd_fermi_1d
		// FCN_WD_FCN_CPU(fermi, Wd_Fermi_2d, 2);		// fcn_wd_fermi_2d
		// FCN_WD_FCN_CPU(fermi, Wd_Fermi_3d, 3);		// fcn_wd_fermi_3d

		// FCN_WD_FCN_CPU(butwth, Wd_Butwth_1d, 1);	// fcn_wd_butwth_1d
		// FCN_WD_FCN_CPU(butwth, Wd_Butwth_2d, 2);	// fcn_wd_butwth_2d
		// FCN_WD_FCN_CPU(butwth, Wd_Butwth_3d, 3);	// fcn_wd_butwth_3d

		// FCN_WD_FCN_CPU(hann, Wd_Hann_1d, 1);		// fcn_wd_hann_1d
		// FCN_WD_FCN_CPU(hann, Wd_Hann_2d, 2);		// fcn_wd_hann_2d
		// FCN_WD_FCN_CPU(hann, Wd_Hann_3d, 3);		// fcn_wd_hann_3d

		#define FCN_WD_FCN_CPU(FN, FCN, DIM)																															\
		template <class T, class TVctr>																																	\
		enable_if_vctr_cpu<TVctr, void>																																	\
		fcn_wd_##FN##_##DIM##d(Grid_##DIM##d<T>& grid, FCN<T>& fcn, dt_bool sft, TVctr& mx_o, Stream_cpu* pstream = nullptr)											\
		{																																								\
			using U = Value_type<TVctr>;																																\
																																										\
			if (sft)																																					\
			{																																							\
				fcn_stream_exec_xd_krn<EDIM(DIM)>(pstream, F_SHAPE_ND(grid, DIM), cgpu_detail::fcn_wd_fcn_r2_sft_##DIM##d<T, U, FCN<T>>, grid, fcn, mx_o.m_data);					\
			}																																							\
			else																																						\
			{																																							\
				fcn_stream_exec_xd_krn<EDIM(DIM)>(pstream, F_SHAPE_ND(grid, DIM), cgpu_detail::fcn_wd_fcn_r2_##DIM##d<T, U, FCN<T>>, grid, fcn, mx_o.m_data);						\
			}																																							\
		}

		 FCN_WD_FCN_CPU(gauss, Wd_Gauss_1d, 1);			// fcn_wd_gauss_1d
		 FCN_WD_FCN_CPU(gauss, Wd_Gauss_2d, 2);			// fcn_wd_gauss_2d
		 FCN_WD_FCN_CPU(gauss, Wd_Gauss_3d, 3);			// fcn_wd_gauss_3d

		 FCN_WD_FCN_CPU(exp, Wd_Exp_1d, 1);				// fcn_wd_exp_1d
		 FCN_WD_FCN_CPU(exp, Wd_Exp_2d, 2);				// fcn_wd_exp_2d
		 FCN_WD_FCN_CPU(exp, Wd_Exp_3d, 3);				// fcn_wd_exp_3d

		 FCN_WD_FCN_CPU(fermi, Wd_Fermi_1d, 1);			// fcn_wd_fermi_1d
		 FCN_WD_FCN_CPU(fermi, Wd_Fermi_2d, 2);			// fcn_wd_fermi_2d
		 FCN_WD_FCN_CPU(fermi, Wd_Fermi_3d, 3);			// fcn_wd_fermi_3d

		 FCN_WD_FCN_CPU(butwth, Wd_Butwth_1d, 1);		// fcn_wd_butwth_1d
		 FCN_WD_FCN_CPU(butwth, Wd_Butwth_2d, 2);		// fcn_wd_butwth_2d
		 FCN_WD_FCN_CPU(butwth, Wd_Butwth_3d, 3);		// fcn_wd_butwth_3d

		 FCN_WD_FCN_CPU(hann, Wd_Hann_1d, 1);			// fcn_wd_hann_1d
		 FCN_WD_FCN_CPU(hann, Wd_Hann_2d, 2);			// fcn_wd_hann_2d
		 FCN_WD_FCN_CPU(hann, Wd_Hann_3d, 3);			// fcn_wd_hann_3d
	}

	/* phase correlation */
	namespace mt
	{
		/****************** pcf data processing real space *******************/
		template <class T, class TVctr_r, class TVctr_c>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_r, TVctr_c, void>
		fcn_rs_pcf_1d_dp(Grid_1d<T>& grid, TVctr_r& mx_i, Wd_Butwth_1d<T>& fcn, T w, TVctr_c& mx_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			fcn_stream_exec_xd_krn<edim_1>(pstream, grid.nx, cgpu_detail::fcn_rs_pcf_1d_dp<T, U>, grid, mx_i.m_data, fcn, w, mx_o.m_data);
		}

		template <class T, class TVctr_r, class TVctr_c>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_r, TVctr_c, void>
		fcn_rs_pcf_2d_dp(Grid_2d<T>& grid, TVctr_r& mx_i, Wd_Butwth_2d<T>& fcn, T w, TVctr_c& mx_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail::fcn_rs_pcf_2d_dp<T, U>, grid, mx_i.m_data, fcn, w, mx_o.m_data);
		}
		
		template <class T, class TVctr_r, class TVctr_c>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_r, TVctr_c, void>
		fcn_rs_pcf_3d_dp(Grid_3d<T>& grid, TVctr_r& mx_i, Wd_Butwth_3d<T>& fcn, T w, TVctr_c& mx_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			fcn_stream_exec_xd_krn<edim_3>(pstream, grid.nx, grid.ny, grid.nz, cgpu_detail::fcn_rs_pcf_3d_dp<T, U>, grid, mx_i.m_data, fcn, w, mx_o.m_data);
		}

		/***************** pcf data processing fourier space *****************/
		template <class T, class TVctr_r, class TVctr_c>
		enable_if_cvctr_gpu_and_vctr_rgpu<TVctr_c, TVctr_r, void>
		fcn_fs_pcf_1d_dp(Grid_1d<T>& grid, TVctr_r& mx_i, Wd_Gauss_1d<T>& fcn, T w, TVctr_c& mx_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			fcn_stream_exec_xd_krn<edim_1>(pstream, grid.nx, cgpu_detail::fcn_fs_pcf_1d_dp<T, U>, grid, mx_i.m_data, fcn, w, mx_o.m_data);
		}

		template <class T, class TVctr_r, class TVctr_c>
		enable_if_cvctr_gpu_and_vctr_rgpu<TVctr_c, TVctr_r, void>
		fcn_fs_pcf_2d_dp(Grid_2d<T>& grid, TVctr_r& mx_i, Wd_Gauss_2d<T>& fcn, T w, TVctr_c& mx_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail::fcn_fs_pcf_2d_dp<T, U>, grid, mx_i.m_data, fcn, w, mx_o.m_data);
		}
		
		template <class T, class TVctr_r, class TVctr_c>
		enable_if_cvctr_gpu_and_vctr_rgpu<TVctr_c, TVctr_r, void>
		fcn_fs_pcf_3d_dp(Grid_3d<T>& grid, TVctr_r& mx_i, Wd_Gauss_3d<T>& fcn, T w, TVctr_c& mx_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			fcn_stream_exec_xd_krn<edim_3>(pstream, grid.nx, grid.ny, grid.nz, cgpu_detail::fcn_fs_pcf_3d_dp<T, U>, grid, mx_i.m_data, fcn, w, mx_o.m_data);
		}
	}

	/* optical flow */
	namespace mt
	{
		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_opt_flow(TVctr& mx_s, TVctr& mx_m, T alpha, TVctr& dphi_x, TVctr& dphi_y, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			auto igrid = mx_s.igrid_2d();

			fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, igrid.ny, cgpu_detail::fcn_opt_flow<U>, igrid, mx_s.m_data, mx_m.m_data, alpha, dphi_x.m_data, dphi_y.m_data);
		}
	}

	/* bilinear interpolation */
	namespace mt
	{
		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_intrpl_bl_rg_2d(Grid_2d<T>& grid, TVctr& mx_i, TVctr& vx, TVctr& vy, T bg, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail::fcn_intrpl_bl_rg_2d<U>, grid, mx_i.m_data, vx.m_data, vy.m_data, bg, mx_o.m_data);
		}
	}

	/* histogram */
	namespace mt
	{
		template <class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		fcn_hist(const TVctr_1& vctr, const Value_type<TVctr_2>& n_bins, TVctr_2& y, TVctr_2& x)
		{
			using T1 = Value_type<TVctr_1>;
			using T2 = Value_type<TVctr_2>;
			using size_type = Size_type<TVctr_2>;

			T1 vctr_min = 0;
			T1 vctr_max = 0;
			fcn_minmax_element(vctr, vctr_min, vctr_max);
			auto dv = dt_float64(vctr_max-vctr_min)/dt_float64(n_bins);

			fcn_fill(y, T2(0));

			const auto vctr_size = vctr.size();

			for(size_type iv = 0; iv< vctr_size; iv++)
			{
				auto val = dt_float64(vctr[iv]);
				auto iy = fcn_bcfloor((val - dt_float64(vctr_min))/dv, size_type(0), size_type(n_bins-1));
				y[iy]++;
			}

			for(size_type ix = 0; ix < n_bins; ix++)
			{
				x[ix] = T2(vctr_min + dt_float64(ix)*dv);
			}
		}
	}
	
	/* anscombe transform */
	namespace mt
	{
		/* forward anscombe transform */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_anscombe(TVctr& mx_i, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

		 	auto thr_anscombe_fwd = [](const iThread_Rect_1d& ithread, pVctr_cpu_32<T>& mx_i, pVctr_cpu_32<T>& mx_o)
		 	{
		 		thrust::transform(mx_i.begin()+ithread.ind_0, mx_i.begin()+ithread.ind_e, 
		 		mx_o.begin()+ithread.ind_0, cgpu_fctr::anscombe_fwd<T>());
		 	};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_i.size_32(), thr_anscombe_fwd, mx_i.ptr_32(), mx_o.ptr_32());
		}

		/* forward anscombe transform */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_anscombe_inv(TVctr& mx_i, TVctr& mx_o, Stream_cpu* pstream = nullptr)
		{
			using T = Value_type<TVctr>;

		 	auto thr_anscombe_inv = [](const iThread_Rect_1d& ithread, pVctr_cpu_32<T>& mx_i, pVctr_cpu_32<T>& mx_o)
		 	{
		 		thrust::transform(mx_i.begin()+ithread.ind_0, mx_i.begin()+ithread.ind_e, 
		 		mx_o.begin()+ithread.ind_0, cgpu_fctr::anscombe_inv<T>());
		 	};

			fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_i.size_32(), thr_anscombe_inv, mx_i.ptr_32(), mx_o.ptr_32());
		}
	}

	/* border */
	namespace mt
	{
		 /* replace border */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_repl_bdr(TVctr& mx_i, iBorder_Rect_2d& iborder, eFil_Sel_Typ bg_opt=efst_min, Value_type<TVctr> bg=0, Stream_cpu* pstream = nullptr)
		{
		 	using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();
			iRegion_Rect_2d iregion({iborder.bx_0, igrid.nx - iborder.bx_e, iborder.by_0, igrid.ny - iborder.by_e});

		 	dt_int32 ic = 0;
		 	T v_min = mx_i(iregion.ry_0, iregion.rx_0);
		 	T v_max = v_min;
		 	KS<dt_float64> v_mean = 0;

			Vctr_cpu<T> mx_o(mx_i);

			if (bg_opt<6)
			{
				for(auto ix = iregion.rx_0; ix < iregion.rx_e; ix++)
				{
					for(auto iy = iregion.ry_0; iy < iregion.ry_e; iy++)
					{
						const auto v = mx_i(iy, ix);
						v_min = fcn_min(v_min, v);
						v_max = fcn_max(v_max, v);
						v_mean += dt_float64(v);
						ic++;
					}
				}
		 		v_mean = v_mean/T(ic);
			}

		 	auto val_f = fcn_select_bg(bg_opt, v_min, v_max, T(v_mean()), bg);

		 	// left
		 	for(auto ix = 0; ix < iregion.rx_0; ix++)
		 	{
		 		for(auto iy = 0; iy < igrid.ny; iy++)
		 		{
		 			mx_o(iy, ix) = val_f;
		 		}
		 	}

		 	// right
		 	for(auto ix = iregion.rx_e; ix < igrid.nx; ix++)
		 	{
		 		for(auto iy = 0; iy < igrid.ny; iy++)
		 		{
		 			mx_o(iy, ix) = val_f;
		 		}
		 	}

		 	// top
		 	for(auto ix = 0; ix < igrid.nx; ix++)
		 	{
		 		for(auto iy = 0; iy < iregion.ry_0; iy++)
		 		{
		 			mx_o(iy, ix) = val_f;
		 		}
		 	}

		 	// bottom
		 	for(auto ix = 0; ix < igrid.nx; ix++)
		 	{
		 		for(auto iy = iregion.ry_e; iy < igrid.ny; iy++)
		 		{
		 			mx_o(iy, ix) = val_f;
		 		}
		 	}

			return mx_o;
		 }

		 /* add pbc border */
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, Vctr_cpu<Value_type<TVctr>>>
		fcn_add_pbdr(TVctr& mx_i, iBorder_Rect_2d& iborder)
		 {
		 	using T = Value_type<TVctr>;

			auto igrid = mx_i.igrid_2d();
			iRegion_Rect_2d iregion({iborder.bx_0, igrid.nx + iborder.bx_0, iborder.by_0, igrid.ny + iborder.by_0});

			Vctr_cpu<T> mx_o(dt_shape(igrid.ny + iborder.by_sum(), igrid.nx + iborder.bx_sum()));
			igrid.set_size(mx_o.s1_32(), mx_o.s0_32());

		 	// copy central mx_i
		 	for(auto ix = iregion.rx_0; ix < iregion.rx_e; ix++)
		 	{
		 		for(auto iy = iregion.ry_0; iy < iregion.ry_e; iy++)
		 		{
		 			const auto ix_i = ix - iregion.rx_0;
		 			const auto iy_i = iy - iregion.ry_0;
					mx_o(iy, ix) = mx_i(iy_i, ix_i);
		 		}
		 	}

		 	// left
		 	for(auto ix = 0; ix < iregion.rx_0; ix++)
		 	{
		 		for(auto iy = iregion.ry_0; iy < iregion.ry_e; iy++)
		 		{
		 			const auto ix_i = 2*iregion.rx_0 - ix;
		 			mx_o(iy, ix) = mx_o(iy, ix_i);
		 		}
		 	}

		 	// right
		 	for(auto ix = iregion.rx_e; ix < igrid.nx; ix++)
		 	{
		 		for(auto iy = iregion.ry_0; iy < iregion.ry_e; iy++)
		 		{
		 			const auto ix_i = max(0, iregion.rx_e - 1 - (ix - iregion.rx_e));
		 			mx_o(iy, ix) = mx_o(iy, ix_i);
		 		}
		 	}

		 	// top
		 	for(auto ix = 0; ix < igrid.nx; ix++)
		 	{
		 		for(auto iy = 0; iy < iregion.ry_0; iy++)
		 		{
		 			const auto iy_i = 2*iregion.ry_0 - iy;
		 			mx_o(iy, ix) = mx_o(iy_i, ix);
		 		}
		 	}

		 	// bottom
		 	for(auto ix = 0; ix < igrid.nx; ix++)
		 	{
		 		for(auto iy = iregion.ry_e; iy < igrid.ny; iy++)
		 		{
		 			const auto iy_i = max(0, iregion.ry_e - 1 - (iy - iregion.ry_e));
		 			mx_o(iy, ix) = mx_o(iy_i, ix);
		 		}
		 	}

		 	return mx_o;
		 }
	}

	/* radial distribution */
	namespace mt
	{
		 // radial distribution: typ = 0; % 0: sum(v)/n, 1: sum(v), 2: cumsum(sum(v)/n), 3: cumsum(sum(v))
		 template <class TVctr, class TFcn>
		 enable_if_vctr_cpu<TVctr, void>
		 fcn_rad_dist(TVctr& r_i, TVctr& fr_i, Value_type<TVctr> r_0, Value_type<TVctr> r_e, 
		 dt_int32 n_r, TVctr& rl_o, TVctr& frl_o, TVctr& cfrl_o, dt_int32 typ, const TFcn& fcn)
		 {
			using T = Value_type<TVctr>;

		 	const T dr = (r_e-r_0)/T(n_r);

		 	for(auto ir = 0; ir < n_r; ir++)
		 	{
		 		rl_o[ir] = r_0 + T(ir)*dr;
		 		frl_o[ir] = T(0);
		 		cfrl_o[ir] = T(0);
		 	}

		 	for(auto ir = 0; ir < r_i.size(); ir++)
		 	{
				const auto r = r_i[ir];
		 		if (fcn_chk_bound(r, r_0, r_e))
		 		{
					const auto ik = fcn(r);
		 			frl_o[ik] += fr_i[ir];
		 			cfrl_o[ik] += T(1);
		 		}
		 	}

		 	if ((typ == 0) || (typ == 2))
		 	{
		 		for(auto ir = 0; ir < n_r; ir++)
		 		{
		 			frl_o[ir] /= fcn_max(T(1), cfrl_o[ir]);
		 		}
		 	}

		 	if (typ > 1)
		 	{
		 		for(auto ir = 1; ir < n_r; ir++)
		 		{
					frl_o[ir] += frl_o[ir - 1];
		 		}
		 	}
		 }

		 // radial distribution by division: typ = 0; % 0: sum(v)/n, 1: sum(v), 2: cumsum(sum(v)/n), 3: cumsum(sum(v))
		 template <class TVctr>
		 enable_if_vctr_cpu<TVctr, void>
		 fcn_rad_dist_by_div(TVctr& r_i, TVctr& fr_i, Value_type<TVctr> r_0, Value_type<TVctr> r_e, 
		 dt_int32 n_r, TVctr& rl_o, TVctr& frl_o, TVctr& cfrl_o, dt_int32 typ)
		 {
			using T = Value_type<TVctr>;
		 	const T dr = (r_e-r_0)/T(n_r);

			cgpu_fctr::rad_dist_ind_by_div<T> fcn(r_0, dr, 0, n_r-1);

			fcn_rad_dist(r_i, fr_i, r_0, r_e, n_r, rl_o, frl_o, cfrl_o, typ, fcn);
		 }		 
		 
		 // radial distribution by search: typ = 0; % 0: sum(v)/n, 1: sum(v), 2: cumsum(sum(v)/n), 3: cumsum(sum(v))
		 template <class TVctr>
		 enable_if_vctr_cpu<TVctr, void>
		 fcn_rad_dist_by_srch(TVctr& r_i, TVctr& fr_i, Value_type<TVctr> r_0, Value_type<TVctr> r_e, 
		 dt_int32 n_r, TVctr& rl_o, TVctr& frl_o, TVctr& cfrl_o, dt_int32 typ)
		 {
			using T = Value_type<TVctr>;

			cgpu_fctr::rad_dist_ind_by_srch<T> fcn(rl_o.data(), 0, n_r-1);

			fcn_rad_dist(r_i, fr_i, r_0, r_e, n_r, rl_o, frl_o, cfrl_o, typ, fcn);
		 }	
	}

	/* radial distribution 2d */
	namespace mt
	{		 
		 // radial distribution 2d: typ = 0; % 0: sum(v)/n, 1: sum(v), 2: cumsum(sum(v)/n), 3: cumsum(sum(v))
		 template <class TVctr, class TFcn>
		 enable_if_vctr_cpu<TVctr, void>
		 fcn_rad_dist_2d(Grid_2d<Value_type<TVctr>>& grid, TVctr& mx_i, R_2d<Value_type<TVctr>> p_c, 
		 Value_type<TVctr> radius_i, TVctr& rl_o, TVctr& frl_o, TVctr& cfrl_o, dt_int32 typ, const TFcn& fcn)
		 {
			using T = Value_type<TVctr>;

		 	dt_int32 ix_0, ix_e;
			grid.ix_0_ix_e(p_c.x, radius_i, ix_0, ix_e);

		 	dt_int32 iy_0, iy_e;
			grid.iy_0_iy_e(p_c.y, radius_i, iy_0, iy_e);

			const auto r2_max = ::square(radius_i);
			const auto dr = grid.drx;

		 	for(auto ir = 0; ir < rl_o.size(); ir++)
		 	{
		 		rl_o[ir] = grid.rx(ir);
		 		frl_o[ir] = T(0);
		 		cfrl_o[ir] = T(0);
		 	}

		 	for(auto ix = ix_0; ix < ix_e; ix++)
		 	{
		 		for(auto iy = iy_0; iy < iy_e; iy++)
		 		{
		 			const auto r2 = grid.r2(ix, iy, p_c);
		 			if (r2 < r2_max)
		 			{
		 				const auto ir = fcn(::sqrt(r2));
		 				frl_o[ir] += mx_i[grid.sub_2_ind(ix, iy)];
		 				cfrl_o[ir] += T(1);
		 			}
		 		}
		 	}

		 	if ((typ == 0) || (typ == 2))
		 	{
		 		for(auto ir = 0; ir < rl_o.size(); ir++)
		 		{
		 			frl_o[ir] /= fcn_max(T(1), cfrl_o[ir]);
		 		}
		 	}

		 	if (typ > 1)
		 	{
		 		for(auto ir = 1; ir < rl_o.size(); ir++)
		 		{
					frl_o[ir] += frl_o[ir - 1];
		 		}
		 	}
		 }


		 // radial distribution 2d by division
		 template <class TVctr>
		 enable_if_vctr_cpu<TVctr, void>
		 fcn_rad_dist_2d_by_div(Grid_2d<Value_type<TVctr>>& grid, TVctr& mx_i, R_2d<Value_type<TVctr>> p_c, 
		 Value_type<TVctr> radius_i, TVctr& rl_o, TVctr& frl_o, TVctr& cfrl_o, dt_int32 typ)
		 {
			using T = Value_type<TVctr>;

			cgpu_fctr::rad_dist_ind_by_div<T> fcn(0, grid.drx, 0, rl_o.size()-1);

			fcn_rad_dist_2d(grid, mx_i, p_c, radius_i, rl_o, frl_o, cfrl_o, typ, fcn);
		 }		 
		 
		 // radial distribution 2d by search
		 template <class TVctr>
		 enable_if_vctr_cpu<TVctr, void>
		 fcn_rad_dist_2d_by_srch(Grid_2d<Value_type<TVctr>>& grid, TVctr& mx_i, R_2d<Value_type<TVctr>> p_c, 
		 Value_type<TVctr> radius_i, TVctr& rl_o, TVctr& frl_o, TVctr& cfrl_o, dt_int32 typ)
		 {
			using T = Value_type<TVctr>;

			cgpu_fctr::rad_dist_ind_by_srch<T> fcn(rl_o.data(), 0, rl_o.size()-1);

			fcn_rad_dist_2d(grid, mx_i, p_c, radius_i, rl_o, frl_o, cfrl_o, typ, fcn);
		 }	
	}

	/* information limit for regular grid */
	namespace mt
	{	
		/* filter mean */
		namespace cpu_fcn_image
		{
			/* 1d */
			template <class TVctr>
			enable_if_vctr_cpu<TVctr, void>
			fltr_mean_1d(TVctr& mx_i, dt_int32 n_kr, TVctr& mx_o, Stream_cpu* pstream);
		}

		template <class TGrid, class TVctr>
		Value_type<TVctr> fcn_info_lim_2d(TGrid& grid, TVctr& mx_i)
		{
			using T = Value_type<TVctr>;

			const auto r_c = grid.bs_h();
			const auto radius = grid.bs_h_min();
			const auto n_r = fcn_min(grid.nx_h, grid.ny_h);

			Vctr_cpu<T> r(n_r);
			Vctr_cpu<T> fr(n_r);
			Vctr_cpu<T> cfr(n_r);

			// cumulative radial integration
			fcn_rad_dist_2d_by_div(grid, mx_i.ptr_64(), r_c, grid.bs_h_min(), r.ptr_64(), fr.ptr_64(), cfr.ptr_64(), 2);

			// shift and normalize
			const T r_0 = r.front();	 	
			const T fr_0 = fr.front();	
			const T dr = r.back() - r_0;
			const T dfr = fr.back() - fr_0;
			for(auto ir = 0; ir < fr.size(); ir++)
			{
				r[ir] = (r[ir] - r_0)/dr;
				fr[ir] = (fr[ir] - fr_0)/dfr;
			}

			// fcn_smooth_mean_1d data
			const dt_int32 nkr = max<dt_int32>(5, fr.size()/100);
			cpu_fcn_image::fltr_mean_1d(fr, nkr, cfr);
			fr = cfr;

			// get maximum curvature point
			dt_int32 ir_m = 0;
			auto d_max = fcn_get_max_crv_2_ln_dist(r.front(), fr.front(), r.back(), fr.back(), 0, r.size(), r.data(), fr.data(), ir_m);

			return grid.rx(ir_m);
		}
	}

	/* unique / match vector */
	namespace mt
	{
		template <class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_uniq_vctr(TVctr& vctr)
		{
			using T = Value_type<TVctr>;

			std::sort(vctr.begin(), vctr.end());
			const auto size = std::distance(vctr.begin(), std::unique(vctr.begin(), vctr.end(), mt::fcn_is_equal<T>));

			if (vctr.s0()==1)
			{
				vctr.resize({1, size});
			}
			else
			{
				vctr.resize({size, 1});
			}
		}

		template <class InputIterator, class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_match_vctr(InputIterator v_first, InputIterator v_last, TVctr& vs, Value_type<TVctr> dv=0.1)
		{
			using T = Value_type<TVctr>;

			TVctr vr(v_first, v_last);

			// find min and max values
			T v_min, v_max;
			fcn_minmax_element(vr, v_min, v_max);

			// add tolerance
			v_min -= dv;
			v_max += dv;

			// remove elements of vs outside the ithread
			auto fcn_rm = [v_min, v_max](T a)->dt_bool
			{
				return ((a < v_min) || (a > v_max));
			};

			const auto size = std::distance(vr.begin(), std::remove_if(vs.begin(), vs.end(), fcn_rm));
			vr.resize(size);

			vs.erase(std::remove_if(vs.begin(), vs.end(), fcn_rm), vs.end());

			Vctr<dt_bool, edev_cpu> v_c(vr.size(), true);
			dt_int32 m_size = min(vr.size(), vs.size());
			TVctr A_d;
			A_d.reserve(m_size);

			for(auto i = 0; i < vs.size(); i++)
			{
				T val = vs[i];
				auto it = std::min_element(vr.begin(), vr.end(), cgpu_fctr::closest_element<T>(val));
				auto imin = static_cast<dt_int32>(it - vr.begin());

				if (v_c[imin])
				{
					v_c[imin] = false;
					A_d.push_back(*it);
				}
			}
			vs = A_d;
			vs.shrink_to_fit();
		}
	}

	namespace mt
	{
		// /***************************************************************************************/
		// /***************************************************************************************/
		// template <class TVctr_i, class TVctr_o>
		// enable_if_vctr_cpu_and_vctr_cpu<TVctr_i, TVctr_o, void>
		// cpy_to_host(Stream_cpu& stream, TVctr_i &mx_i, TVctr_o &mx_o, 
		// TVctr_i *M_i_h = nullptr)
		// {
		// 	mt::assign(mx_i, mx_o, M_i_h);
		// }

		// template <class TVctr_i, class TVctr_o>
		// enable_if_vctr_cpu_and_vctr_cpu<TVctr_i, TVctr_o, void>
		// add_sc_to_host(Stream_cpu& stream, Value_type<TVctr_i> w_i, 
		// TVctr_i &mx_i, TVctr_o &mx_o, TVctr_i *M_i_h = nullptr)
		// {
		// 	mt::add_scale(stream, w_i, mx_i, mx_o);
		// }

		// template <class TVctr_i, class TVctr_o>
		// enable_if_vctr_cpu_and_vctr_cpu<TVctr_i, TVctr_o, void>
		// add_sc_square_to_host(Stream_cpu& stream, Value_type<TVctr_o> w_i, 
		// TVctr_i &mx_i, TVctr_o &mx_o, TVctr_i *M_i_h = nullptr)
		// {
		// 	mt::add_scale_norm_2(stream, w_i, mx_i, mx_o);
		// }

		// template <class TVctr_c_i, class TVctr_r_o, class TVctr_c_o>
		// enable_if_vctr_cpu_and_vctr_cpu<TVctr_c_i, TVctr_c_o, void>
		// add_sc_m2psi_psi_to_host(Stream_cpu& stream, Value_type<TVctr_r_o> w_i, 
		// TVctr_c_i &psi_i, TVctr_r_o &m2psi_o, TVctr_c_o &psi_o, TVctr_c_i *psi_i_h = nullptr)
		// {
		// 	mt::add_scale(stream, w_i, psi_i, psi_o);
		// 	mt::add_scale_norm_2(stream, w_i, psi_i, m2psi_o);
		// }
		// /***************************************************************************************/
		// /***************************************************************************************/
		// template <class TVctr>
		// TVctr lsf_poly_n(TVctr& x, TVctr& y, dt_int32 n)
		// {
		// 	using T = Value_type<TVctr>;

		// 	dt_int32 m = x.size();
		// 	n++;		// add bias term

		// 	TVctr M(m*n);
		// 	for(auto in = 0; in < n; in++)
		// 	{
		// 		for(auto im = 0; im < m; im++)
		// 		{
		// 			M[in*m + im] = (in == 0)?1:x[im] * M[(in - 1)*m + im];
		// 		}
		// 	}

		// 	TVctr coef(n);
		// 	lapack::GELS<T> gels;
		// 	gels(m, n, M.data(), 1, y.data(), coef.data());

		// 	return coef;
		// }

		// /***************************************************************************************/
		// // extract region ix_0<=x<ix_e
		// template <class TVctr_o, class TVctr_i>
		// TVctr_o extract_region_real_part(Stream_cpu& stream, dt_int32 nx_src, dt_int32 ny_src, 
		// 	TVctr_i &Im_src, dt_int32 ix0_src, dt_int32 ixe_src, dt_int32 iy0_src, dt_int32 iye_src)
		// {
		// 	auto nx_dst = ixe_src - ix0_src;
		// 	auto ny_dst = iye_src - iy0_src;

		// 	TVctr_o Im_dst(nx_dst*ny_dst);

		// 	auto thr_extract_region = [&](const iRegion_Rect_2d& ithread)
		// 	{
		// 		for(auto ix = ithread.ix_0; ix < ithread.ix_e; ix++)
		// 		{
		// 			for(auto iy = ithread.iy_0; iy < ithread.iy_e; iy++)
		// 			{
		// 				Im_dst[ix*ny_dst + iy] = Im_src[(ix0_src + ix)*ny_src + (iy0_src + iy)].real();
		// 			}
		// 		}
		// 	};

		// 	stream.set_n_stream_act(nx_dst);
		// 	stream.set_grid(nx_dst, ny_dst);
		// 	stream.exec_xd_fcn<edim_1>(thr_extract_region);

		// 	return Im_dst;
		// }

		// // extract real part of vector
		// template <class TVctr_o, class TVctr_i>
		// TVctr_o extract_real_part(Stream_cpu& stream, TVctr_i &Im_src)
		// {
		// 	TVctr_o Im_dst(Im_src.size());

		// 	auto thr_extract_real = [](const iThread_Rect_2d& ithread, TVctr_i &Im_src, TVctr_o &Im_dst)
		// 	{
		// 		for(auto ixy = ithread.ind_0; ixy < ithread.ind_e; ixy++)
		// 		{
		// 			Im_dst[ixy] = Im_src[ixy].real();
		// 		}
		// 	};

		// 	stream.set_n_stream_act(Im_src.size());
		// 	stream.set_grid(Im_src.size(), 1);
		// 	stream.exec_xd_fcn<edim_1>(thr_extract_real, Im_src, Im_dst);

		// 	return Im_dst;
		// }

		// // extract region ix_0<=x<ix_e
		// template <class TVctr_o, class TVctr_i>
		// TVctr_o extract_region_abs(Stream_cpu& stream, dt_int32 nx_src, dt_int32 ny_src, 
		// 	TVctr_i &Im_src, dt_int32 ix0_src, dt_int32 ixe_src, dt_int32 iy0_src, dt_int32 iye_src)
		// {
		// 	auto nx_dst = ixe_src - ix0_src;
		// 	auto ny_dst = iye_src - iy0_src;

		// 	TVctr_o Im_dst(nx_dst*ny_dst);

		// 	auto thr_extract_region = [&](const iThread_Rect_2d& ithread)
		// 	{
		// 		for(auto ix = ithread.ix_0; ix < ithread.ix_e; ix++)
		// 		{
		// 			for(auto iy = ithread.iy_0; iy < ithread.iy_e; iy++)
		// 			{
		// 				Im_dst[ix*ny_dst + iy] = fabs(Im_src[(ix0_src + ix)*ny_src + (iy0_src + iy)]);
		// 			}
		// 		}
		// 	};

		// 	stream.set_n_stream_act(nx_dst);
		// 	stream.set_grid(nx_dst, ny_dst);
		// 	stream.exec_xd_fcn<edim_1>(thr_extract_region);

		// 	return Im_dst;
		// }

		// // extract abs of vector
		// template <class TVctr_o, class TVctr_i>
		// TVctr_o extract_abs(Stream_cpu& stream, TVctr_i &Im_src)
		// {
		// 	TVctr_o Im_dst(Im_src.size());

		// 	auto thr_extract_abs = [](const iThread_Rect_2d& ithread, TVctr_i &Im_src, TVctr_o &Im_dst)
		// 	{
		// 		for(auto ixy = ithread.ind_0; ixy < ithread.ind_e; ixy++)
		// 		{
		// 			Im_dst[ixy] = fabs(Im_src[ixy]);
		// 		}
		// 	};

		// 	stream.set_n_stream_act(Im_src.size());
		// 	stream.set_grid(Im_src.size(), 1);
		// 	stream.exec_xd_fcn<edim_1>(thr_extract_abs, Im_src, Im_dst);

		// 	return Im_dst;
		// }

		// /***************************************************************************************/
		// // get weighted position
		// template <class TGrid, class TVctr>
		// R_2d<Value_type<TVctr>> Rx_Ry_weight(TGrid& grid_2d, TVctr& Im, R_2d<Value_type<TVctr>> p_i, Value_type<TVctr> radius_i, dt_bool env)
		// {
		// 	using T = Value_type<TVctr>;

		// 	T dR = grid_2d.drx;
		// 	dt_int32 nrl = static_cast<dt_int32>(::floor(radius_i/dR + 0.5));
		// 	T R_max = dR*nrl;

		// 	auto ix_i = static_cast<dt_int32>(::floor(p_i.x/dR));
		// 	dt_int32 ix_0 = max(ix_i - nrl, 0);
		// 	dt_int32 ix_e = min(ix_i + nrl + 1, grid_2d.nx);

		// 	auto iy_i = static_cast<dt_int32>(::floor(p_i.y/dR));
		// 	dt_int32 iy_0 = max(iy_i - nrl, 0);
		// 	dt_int32 iy_e = min(iy_i + nrl + 1, grid_2d.ny);

		// 	T R2_max = pow(R_max, 2);

		// 	R_2d<T> p(0, 0);

		// 	T alpha = 0.5/pow(R_max, 2);
		// 	T wt = 0;
		// 	for(auto ix = ix_0; ix < ix_e; ix++)
		// 	{
		// 		for(auto iy = iy_0; iy < iy_e; iy++)
		// 		{
		// 			auto R2 = grid_2d.r2(ix, iy, p_i.x, p_i.y);
		// 			if (R2 < R2_max)
		// 			{
		// 				dt_int32 ixy = grid_2d.sub_2_ind(ix, iy);
		// 				auto imv = ((env)?exp(-alpha*R2):1)*Im[ixy];
		// 				p += imv*R_2d<T>(grid_2d.rx(ix), grid_2d.ry(iy));
		// 				wt += imv;
		// 			}
		// 		}
		// 	}
		// 	p /= wt;

		// 	if (norm(p - p_i) >= radius_i)
		// 	{
		// 		p = p_i;
		// 	}

		// 	return p;
		// }

		// // fit gaussian 1d
		// template <class TGrid, class TVctr>
		// TVctr fit_gauss_1d(TGrid& grid_1d, TVctr& Im_i, Value_type<TVctr> x_i, 
		// 	Value_type<TVctr> sigma_i, Value_type<TVctr> radius_i)
		// {
		// 	using T = Value_type<TVctr>;

		// 	auto select_cir_reg = [](TGrid& grid_1d, TVctr& Im, T x0, 
		// 		T radius, TVctr& Rx, TVctr& Ix, T& Rx_sf, T& Rx_sc, T& Ix_sc)
		// 	{
		// 		T R_max = radius;
		// 		T R2_max = pow(R_max, 2);

		// 		auto ithread = grid_1d.region_ind(x0, R_max);

		// 		Rx.clear();
		// 		Rx.reserve(ithread.ind_e);
		// 		Ix.clear();
		// 		Ix.reserve(ithread.ind_e);

		// 		// select circular region
		// 		for(auto ix = ithread.ix_0; ix < ithread.ix_e; ix++)
		// 		{
		// 			if (grid_1d.r2(ix, x0) < R2_max)
		// 			{
		// 				Rx.push_back(grid_1d.rx(ix));
		// 				Ix.push_back(Im[ix]);
		// 			}
		// 		}

		// 		Rx_sf = x0;
		// 		Rx_sc = R_max;

		// 		Ix_sc = fcn_max_element(Ix);
		// 		dt_int32 m = Ix.size();
		// 		for(auto ix = 0; ix < m; ix++)
		// 		{
		// 			Rx[ix] = (Rx[ix] - Rx_sf)/Rx_sc;
		// 			Ix[ix] = Ix[ix]/Ix_sc;
		// 		}

		// 		Rx.shrink_to_fit();
		// 		Ix.shrink_to_fit();
		// 	};

		// 	TVctr Rx, Ix;
		// 	T Rx_sf, Rx_sc, Ix_sc;

		// 	select_cir_reg(grid_1d, Im_i, x_i, radius_i, Rx, Ix, Rx_sf, Rx_sc, Ix_sc);

		// 	T sigma = sigma_i/Rx_sc;

		// 	dt_int32 m = Ix.size();
		// 	dt_int32 n = 3;

		// 	TVctr J(m*n);
		// 	TVctr d_Ix(m);

		// 	auto get_chi2 = [](TVctr& Rx, TVctr& Ix, TVctr& coef)->T
		// 	{
		// 		T x0 = coef[0];
		// 		T A = coef[1];
		// 		T alpha = 0.5/pow(coef[2], 2);

		// 		T chi2 = 0;
		// 		T chi2_ee = 0;
		// 		const dt_int32 m = Ix.size();
		// 		for(auto im = 0; im < m; im++)
		// 		{
		// 			T x = Rx[im] - x0;
		// 			T v = Ix[im] - A*exp(-alpha*x*x);
		// 			fcn_kh_sum(chi2, v*v, chi2_ee);
		// 		}
		// 		return chi2;
		// 	};

		// 	auto get_dIx_J = [](TVctr& Rx, TVctr& Ix, TVctr& coef, 
		// 		TVctr& dIx, TVctr& J)
		// 	{
		// 		T x0 = coef[0];
		// 		T A = coef[1];
		// 		T alpha = 0.5/pow(coef[2], 2);

		// 		T c_x0 = 1.0/pow(coef[2], 2);
		// 		T c_A = 1;
		// 		T c_sx = 1.0/pow(coef[2], 3);

		// 		const dt_int32 m = Ix.size();
		// 		for(auto im = 0; im < m; im++)
		// 		{
		// 			T x = Rx[im] - x0;
		// 			T v = exp(-alpha*x*x);
		// 			T f = A*v;

		// 			J[0 * m + im] = c_x0*x*f;
		// 			J[1 * m + im] = c_A*v;
		// 			J[2 * m + im] = c_sx*x*x*f;

		// 			dIx[im] = Ix[im] - f;
		// 		}
		// 	};

		// 	T dx_max = ::fmax(grid_1d.drx/Rx_sc, 0.25*::fmin(sigma, radius_i/Rx_sc));

		// 	vector<T> c_coef_0 = { 0, 1, sigma };
		// 	vector<T> c_coef_min = { -dx_max, 0.5, sigma/3 };
		// 	vector<T> c_coef_max = { dx_max, 1.25, 3 * sigma };

		// 	TVctr coef_0 = c_coef_0;
		// 	TVctr coef_min = c_coef_min;
		// 	TVctr coef_max = c_coef_max;
		// 	TVctr coef = coef_0;

		// 	T chi2 = get_chi2(Rx, Ix, coef);

		// 	T lambda = 2;
		// 	T lambda_f = 2;

		// 	lapack::LSF_1<T> fls_1;

		// 	const dt_int32 niter = 100;
		// 	for(auto iter = 0; iter < niter; iter++)
		// 	{
		// 		get_dIx_J(Rx, Ix, coef, d_Ix, J);

		// 		TVctr d_coef(n);
		// 		TVctr D(n);
		// 		TVctr G(n);
		// 		fls_1(m, n, J.data(), d_Ix.data(), d_coef.data(), lambda, D.data(), G.data());

		// 		TVctr coef_t = coef;
		// 		T rho_f = 0;
		// 		T G_max = 0;
		// 		for(auto ic = 0; ic < n; ic++)
		// 		{
		// 			coef_t[ic] += d_coef[ic];
		// 			rho_f += coef_t[ic] * (D[ic] * coef_t[ic] + G[ic]);
		// 			G_max = ::fmax(G_max, fabs(G[ic]));
		// 		}

		// 		T chi2_t = get_chi2(Rx, Ix, coef_t);
		// 		T rho = (chi2 - chi2_t)/rho_f;

		// 		if ((G_max < 5e-7) || fabs(chi2 - chi2_t) < 5e-8)
		// 		{
		// 			break;
		// 		}

		// 		if (rho > 0)
		// 		{
		// 			coef = coef_t;

		// 			for(auto ic = 0; ic < n; ic++)
		// 			{
		// 				coef[ic] = min(max(coef[ic], coef_min[ic]), coef_max[ic]);
		// 			}

		// 			chi2 = get_chi2(Rx, Ix, coef);

		// 			lambda = (rho > 1e-6)?(::fmax(lambda/lambda_f, 1e-7)):lambda;
		// 		}
		// 		else
		// 		{
		// 			lambda = ::fmin(lambda*lambda_f, 1e+7);
		// 		}
		// 	}

		// 	coef[0] = coef[0] * Rx_sc + Rx_sf;
		// 	coef[1] = coef[1] * Ix_sc;
		// 	coef[2] = coef[2] * Rx_sc;

		// 	return coef;
		// }

		// template <class TGrid, class TVctr>
		// Value_type<TVctr> neighbor_radius(TGrid& grid_i, TVctr& Im_i, R_2d<Value_type<TVctr>> p_i, Value_type<TVctr> radius_i)
		// {
		// 	using T = Value_type<TVctr>;

		// 	T R2_max = pow(radius_i, 2);

		// 	dt_int32 ix_0 = grid_i.rx_2_irx_bfds(p_i.x - radius_i);
		// 	dt_int32 ix_e = grid_i.rx_2_irx_bcds(p_i.x + radius_i);

		// 	dt_int32 iy_0 = grid_i.ry_2_iry_bfds(p_i.y - radius_i);
		// 	dt_int32 iy_e = grid_i.ry_2_iry_bcds(p_i.y + radius_i);

		// 	dt_int32 nrl = grid_i.r_2_ir_cds_dr_min(radius_i);

		// 	// get radial distribution
		// 	TVctr frl(nrl);
		// 	TVctr cfrl(nrl);

		// 	T dR = grid_i.dR_min();

		// 	for(auto ix = ix_0; ix < ix_e; ix++)
		// 	{
		// 		for(auto iy = iy_0; iy < iy_e; iy++)
		// 		{
		// 			T R2_d = grid_i.r2(ix, iy, p_i.x, p_i.y);
		// 			if (R2_d < R2_max)
		// 			{
		// 				auto ir = static_cast<dt_int32>(::floor(::sqrt(R2_d)/dR));
		// 				frl[ir] += Im_i[grid_i.sub_2_ind(ix, iy)];
		// 				cfrl[ir] += 1.0;
		// 			}
		// 		}
		// 	}

		// 	for(auto ir = 0; ir < frl.size(); ir++)
		// 	{
		// 		frl[ir] /= ::fmax(1.0, cfrl[ir]);
		// 	}

		// 	frl[0] = ::fmax(frl[0], 1.01*frl[1]);

		// 	frl = fcn_smooth_mean_1d(frl, 1);

		// 	// derivative and minimum
		// 	dt_int32 ir0 = frl.size() - 1;
		// 	for(auto ir = 0; ir < frl.size() - 1; ir++)
		// 	{
		// 		if (frl[ir] < frl[ir + 1])
		// 		{
		// 			ir0 = ir + 1;
		// 			break;
		// 		}
		// 	}

		// 	dt_int32 irm = frl.size() - 1;
		// 	for(auto ir = ir0; ir < frl.size() - 1; ir++)
		// 	{
		// 		if (frl[ir] > frl[ir + 1])
		// 		{
		// 			irm = ir;
		// 			break;
		// 		}
		// 	}

		// 	return max(2, irm)*dR;
		// }

		// // select circular region
		// template <class TGrid, class TVctr>
		// Value_type<TVctr> mean_cir_reg(TGrid& grid_i, TVctr& Im_i, R_2d<Value_type<TVctr>> p_i, 
		// 	Value_type<TVctr> radius_i, Value_type<TVctr> bg_i)
		// {
		// 	using T = Value_type<TVctr>;

		// 	T R2_max = radius_i*radius_i;

		// 	auto ithread = grid_i.region_ind(p_i, radius_i);

		// 	// select circular region
		// 	T I_mean = 0;
		// 	dt_int32 Ic = 0;
		// 	for(auto ix = ithread.ix_0; ix < ithread.ix_e; ix++)
		// 	{
		// 		for(auto iy = ithread.iy_0; iy < ithread.iy_e; iy++)
		// 		{
		// 			T R2_d = grid_i.r2(ix, iy, p_i.x, p_i.y);
		// 			if (R2_d < R2_max)
		// 			{
		// 				I_mean += Im_i[grid_i.sub_2_ind(ix, iy)] - bg_i;
		// 				Ic++;
		// 			}
		// 		}
		// 	}
		// 	I_mean /= Ic;
		// 	return I_mean;
		// }

		// /***************************************************************************************/
		// template <class TGrid, class TVctr>
		// TVctr intrplprofile(TGrid& grid_2d, TVctr& Im, const R_2d<Value_type<TVctr>>& p1, const R_2d<Value_type<TVctr>>& p2, dt_int32 nr)
		// {
		// 	TVctr v;

		// 	if (!(grid_2d.chk_bound_eps(p1) && grid_2d.chk_bound_eps(p2)))
		// 	{
		// 		return v;
		// 	}

		// 	if (norm(p1 - p2) < grid_2d.dR_min())
		// 	{
		// 		return v;
		// 	}

		// 	using T = Value_type<TGrid>;

		// 	auto fcn_intrpl_bl_rg_2d = [](const R_2d<T>& p, TGrid& grid_2d, TVctr& Im)->T
		// 	{
		// 		auto ix = grid_2d.rx_2_irx_bfds(p.x);
		// 		auto iy = grid_2d.ry_2_iry_bfds(p.y);

		// 		T f11 = Im[grid_2d.sub_2_ind(ix, iy)];
		// 		T f12 = Im[grid_2d.sub_2_ind(ix, iy + 1)];
		// 		T f21 = Im[grid_2d.sub_2_ind(ix + 1, iy)];
		// 		T f22 = Im[grid_2d.sub_2_ind(ix + 1, iy + 1)];

		// 		T x1 = grid_2d.rx(ix);
		// 		T x2 = grid_2d.rx(ix + 1);
		// 		T y1 = grid_2d.ry(iy);
		// 		T y2 = grid_2d.ry(iy + 1);

		// 		T dx1 = p.x - x1;
		// 		T dx2 = x2 - p.x;
		// 		T dy1 = p.y - y1;
		// 		T dy2 = y2 - p.y;

		// 		T f = (dx2*(f11*dy2 + f12*dy1) + dx1*(f21*dy2 + f22*dy1))/((x2 - x1)*(y2 - y1));
		// 		return f;

		// 	};

		// 	R_2d<T> p12 = p2 - p1;
		// 	T mp12 = p12.norm();

		// 	nr = (nr <= 0)?static_cast<dt_int32>(::ceil(mp12/grid_2d.dR_min())):nr;
		// 	nr = max(nr, 2);
		// 	T dr = mp12/(nr - 1);
		// 	R_2d<T> u = dr*normalize(p12);

		// 	v.reserve(nr);
		// 	for(auto ir = 0; ir < nr; ir++)
		// 	{
		// 		R_2d<T> p = p1 + T(ir)*u;
		// 		v.push_back(fcn_intrpl_bl_rg_2d(p, grid_2d, Im));
		// 	}
		// 	return v;
		// }

		// template <class TGrid, class TVctr>
		// void intrplprofile(Stream_cpu& stream, TGrid& grid_2d, TVctr& Im, Value_type<TVctr> bg, R_2d<Value_type<TVctr>> p1, R_2d<Value_type<TVctr>> p2, 
		// 	TVctr& x, TVctr& y)
		// {
		// 	x.clear();
		// 	y.clear();

		// 	if (!(grid_2d.chk_bound_eps(p1) && grid_2d.chk_bound_eps(p2)))
		// 	{
		// 		return;
		// 	}

		// 	if (norm(p1 - p2) < grid_2d.dR_min())
		// 	{
		// 		return;
		// 	}

		// 	using T = Value_type<TGrid>;

		// 	auto fcn_intrpl_bl_rg_2d = [](const R_2d<T>& p, TGrid& grid_2d, TVctr& Im)->T
		// 	{
		// 		auto ix = grid_2d.rx_2_irx_bfds(p.x);
		// 		auto iy = grid_2d.ry_2_iry_bfds(p.y);

		// 		T f11 = Im[grid_2d.sub_2_ind(ix, iy)];
		// 		T f12 = Im[grid_2d.sub_2_ind(ix, iy + 1)];
		// 		T f21 = Im[grid_2d.sub_2_ind(ix + 1, iy)];
		// 		T f22 = Im[grid_2d.sub_2_ind(ix + 1, iy + 1)];

		// 		T x1 = grid_2d.rx(ix);
		// 		T x2 = grid_2d.rx(ix + 1);
		// 		T y1 = grid_2d.ry(iy);
		// 		T y2 = grid_2d.ry(iy + 1);

		// 		T dx1 = p.x - x1;
		// 		T dx2 = x2 - p.x;
		// 		T dy1 = p.y - y1;
		// 		T dy2 = y2 - p.y;

		// 		T f = (dx2*(f11*dy2 + f12*dy1) + dx1*(f21*dy2 + f22*dy1))/((x2 - x1)*(y2 - y1));
		// 		return f;

		// 	};

		// 	dt_int32 Ixy_1 = Im[grid_2d.rv_2_ir_bfds(p1.x, p1.y)];
		// 	dt_int32 Ixy_2 = Im[grid_2d.rv_2_ir_bfds(p2.x, p2.y)];

		// 	R_2d<T> p12 = p2 - p1;
		// 	T mp12 = p12.norm();
		// 	dt_int32 nr = grid_2d.r_2_ir_cds_dr_min(mp12);
		// 	T dr = mp12/(nr - 1);
		// 	R_2d<T> u = dr*normalize(p12);

		// 	x.reserve(nr);
		// 	y.reserve(nr);
		// 	for(auto ir = 0; ir < nr; ir++)
		// 	{
		// 		x.push_back(ir*dr);
		// 		R_2d<T> p = p1 + T(ir)*u;
		// 		T v = fcn_intrpl_bl_rg_2d(p, grid_2d, Im);
		// 		y.push_back(v - bg);
		// 	}
		// }

		// // find one maximum in all areas above the thr
		// template <class TVctr>
		// TVctr fd_peaks_vector_typ_1(TVctr& x, TVctr& y, Value_type<TVctr> y_thr)
		// {
		// 	using T = Value_type<TVctr>;
		// 	TVctr x_peak(y.size());

		// 	T x_max = 0;
		// 	T y_max = y_thr;
		// 	dt_bool bb = y[0] > y_thr;
		// 	dt_int32 ipeak = 0;
		// 	for(auto ix = 0; ix < y.size(); ix++)
		// 	{
		// 		if (y[ix] > y_thr)
		// 		{
		// 			if (y_max < y[ix])
		// 			{
		// 				y_max = y[ix];
		// 				x_max = x[ix];
		// 			}
		// 			bb = true;
		// 		}
		// 		else if (bb)
		// 		{
		// 			x_peak[ipeak++] = x_max;
		// 			y_max = y_thr;
		// 			bb = false;
		// 		}
		// 	}

		// 	if (bb)
		// 	{
		// 		x_peak[ipeak++] = x_max;
		// 	}

		// 	x_peak.resize(ipeak);
		// 	x_peak.shrink_to_fit();

		// 	return x_peak;
		// }

		// // find all maximum in all areas above the thr
		// template <class TVctr>
		// TVctr fd_peaks_vector_typ_2(TVctr& x, TVctr& y, Value_type<TVctr> y_thr)
		// {
		// 	using T = Value_type<TVctr>;
		// 	TVctr x_peak;
		// 	x_peak.reserve(y.size());

		// 	for(auto ix = 1; ix < y.size() - 1; ix++)
		// 	{
		// 		if (y[ix] > y_thr)
		// 		{
		// 			if ((y[ix - 1] < y[ix]) && (y[ix + 1] < y[ix]))
		// 			{
		// 				x_peak.push_back(x[ix]);
		// 			}
		// 		}
		// 	}
		// 	x_peak.shrink_to_fit();

		// 	return x_peak;
		// }

		// /***************************************************************************************/
		// // find peak position 1d
		// template <class TGrid, class TVctr>
		// enable_if_vctr_cpu<TVctr, Value_type<TVctr>>
		// fit_max_pos_1d(TGrid& grid_1d, TVctr& Im, Value_type<TVctr> p_i, 
		// Value_type<TVctr> sigma_i, Value_type<TVctr> radius_i)
		// {
		// 	using T = Value_type<TVctr>;

		// 	auto ithread = grid_1d.region_ind(p_i, radius_i);
		// 	dt_int32 ix_c = grid_1d.rx_2_irx_fds(p_i);

		// 	dt_int32 ix_0 = ithread.ix_0;
		// 	for(auto ix = ix_c - 1; ix >= ithread.ix_0; ix--)
		// 	{
		// 		if (Im[ix - 1] > Im[ix])
		// 		{
		// 			ix_0 = ix;
		// 			break;
		// 		}
		// 	}

		// 	dt_int32 ix_e = ithread.ix_e;
		// 	for(auto ix = ix_c + 1; ix <= ithread.ix_e; ix++)
		// 	{
		// 		if (Im[ix] < Im[ix + 1])
		// 		{
		// 			ix_e = ix;
		// 			break;
		// 		}
		// 	}

		// 	// if there are few points
		// 	if (fabs(ix_e - ix_0) <= 3)
		// 	{
		// 		T x = 0;
		// 		T sI = 0;
		// 		for(auto ix = ix_0; ix <= ix_e; ix++)
		// 		{
		// 			x += Im[ix] * grid_1d.rx(ix);
		// 			sI += Im[ix];
		// 		}
		// 		x = x/sI;

		// 		return x;
		// 	}

		// 	if ((ix_0 = !ithread.ix_0) || (ix_e = !ithread.ix_e))
		// 	{
		// 		T x_min = grid_1d.rx(ix_0);
		// 		T x_max = grid_1d.rx(ix_e);
		// 		radius_i = ::fmin(p_i - x_min, x_max - p_i);
		// 	}

		// 	auto coef = fit_gauss_1d(grid_1d, Im, p_i, sigma_i, radius_i);

		// 	return coef[0];
		// }

		// /***************************************************************************************/
		// // get projective standard deviation
		// template <class TGrid, class TVctr>
		// TVctr projected_intensity(TGrid& grid_2d, TVctr& mx_i, Value_type<TVctr> np_min, Value_type<TVctr> delta)
		// {
		// 	using T = Value_type<TVctr>;

		// 	vector<R_2d<T>> pts_c = { R_2d<T>(0, 0), R_2d<T>(grid_2d.nx - 1, 0), R_2d<T>(grid_2d.nx - 1, grid_2d.ny - 1), R_2d<T>(0, grid_2d.ny - 1) };

		// 	auto n_pp = static_cast<dt_int32>(::ceil(norm(R_2d<T>(grid_2d.nx, grid_2d.ny)) + 2));
		// 	TVctr y_pt(n_pp);
		// 	TVctr c_pt(n_pp);

		// 	T cos_d, sin_d;
		// 	sincos(delta, &sin_d, &cos_d);

		// 	// get projected points
		// 	auto krn_proj_point = [&](const T& cos_d, const T& sin_d, const R_2d<T>& p)->R_2d<T>
		// 	{
		// 		return R_2d<T>(cos_d*p.x + sin_d*p.y, -sin_d*p.x + cos_d*p.y);
		// 	};

		// 	// get reference corner point
		// 	auto p_0 = krn_proj_point(cos_d, sin_d, pts_c[0]);
		// 	for(auto ixy = 1; ixy < pts_c.size(); ixy++)
		// 	{
		// 		auto p_r = krn_proj_point(cos_d, sin_d, pts_c[ixy]);
		// 		if (p_0.y > p_r.y)
		// 		{
		// 			p_0 = p_r;
		// 		}
		// 	}

		// 	for(auto ix = 0; ix < grid_2d.nx; ix++)
		// 	{
		// 		for(auto iy = 0; iy < grid_2d.ny; iy++)
		// 		{
		// 			auto ixy = grid_2d.sub_2_ind(ix, iy);
		// 			auto p_r = krn_proj_point(cos_d, sin_d, R_2d<T>(ix, iy)) - p_0;
		// 			auto j = static_cast<dt_int32>(::floor(p_r.y));
		// 			y_pt[j] += mx_i[ixy];
		// 			c_pt[j] += 1;
		// 		}
		// 	}

		// 	TVctr V_o;
		// 	V_o.reserve(n_pp);
		// 	for(auto j = 0; j < n_pp; j++)
		// 	{
		// 		if (c_pt[j] > np_min)
		// 		{
		// 			V_o.push_back(y_pt[j]/c_pt[j]);
		// 		}
		// 	}
		// 	V_o.shrink_to_fit();

		// 	return V_o;
		// }

		// // get projective standard deviation
		// template <class TGrid, class TVctr>
		// void PSD(Stream_cpu& stream, TGrid& grid_2d, TVctr& mx_i, Value_type<TVctr> np_min, 
		// 	Value_type<TVctr> del_0, Value_type<TVctr> del_e, Value_type<TVctr> d_del, TVctr& x_o, TVctr& y_o)
		// {
		// 	// get number of angles
		// 	auto n_delta = static_cast<dt_int32>(::floor((del_e - del_0)/d_del + 0.5));
		// 	x_o.resize(n_delta);
		// 	y_o.resize(n_delta);

		// 	auto thr_psd = [&](const iThread_Rect_2d& ithread)
		// 	{
		// 		for(auto idel = ithread.ind_0; idel < ithread.ind_e; idel++)
		// 		{
		// 			auto delta = (del_0 + idel*d_del)*c_pi/180;
		// 			x_o[idel] = delta * 180/c_pi;

		// 			auto pIm = projected_intensity(grid_2d, mx_i, np_min, delta);

		// 			y_o[idel] = fcn_variance(pIm);
		// 		}
		// 	};

		// 	stream.set_n_stream_act(n_delta);
		// 	stream.set_grid(n_delta);
		// 	stream.exec_xd_fcn<edim_1>(thr_psd);
		// }

		// // get projective standard deviation
		// template <class TGrid, class TVctr>
		// TVctr PSD_fd_peaks(Stream_cpu& stream, TGrid& grid_2d, TVctr& mx_i, Value_type<TVctr> np_min, 
		// 	TVctr x_i, TVctr y_i)
		// {
		// 	using T = Value_type<TVctr>;

		// 	T d_del = fabs(x_i[1] - x_i[0]);
		// 	dt_int32 nr = max(1, static_cast<dt_int32>(::ceil(3.0/d_del)));
		// 	auto y_f = fltr_median_1d(stream, y_i, nr);

		// 	T y_mean = 0;
		// 	T y_max = y_i[0] - y_f[0];
		// 	for(auto ix = 0; ix < y_i.size(); ix++)
		// 	{
		// 		y_f[ix] = y_i[ix] - y_f[ix];
		// 		y_mean += y_f[ix];
		// 		y_max = max(y_max, y_f[ix]);
		// 	}
		// 	y_mean /= y_i.size();
		// 	auto y_thr = y_mean + 0.2*(y_mean + y_max);

		// 	auto x_peaks = fd_peaks_vector_typ_1(x_i, y_f, y_thr);

		// 	TVctr x_ref, y_ref;
		// 	for(auto ix = 0; ix < x_peaks.size(); ix++)
		// 	{
		// 		auto delta = x_peaks[ix];
		// 		auto delta_0 = delta - d_del;
		// 		auto delta_e = delta + d_del;
		// 		auto d_delta = 0.1*d_del;
		// 		PSD(stream, grid_2d, mx_i, np_min, delta_0, delta_e, d_delta, x_ref, y_ref);
		// 		std::for_each(y_ref.begin(), y_ref.end(), [](T& y) {y = log(1 + y); });
		// 		auto coef = lsf_poly_n(x_ref, y_ref, 2);
		// 		T px = -coef[1]/(2 * coef[2]);
		// 		if ((px <= delta_0) || (delta_e <= px))
		// 		{
		// 			dt_int32 idx = std::max_element(y_ref.begin(), y_ref.end()) - y_ref.begin();
		// 			px = x_ref[idx];
		// 		}
		// 		x_peaks[ix] = px;
		// 	}

		// 	return x_peaks;
		// }
		// /************************ read matlab binary file ***********************/
		// template <class T_i, class TVctr>
		// enable_if_float<T_i, void>
		// read_mat_matrix(std::ifstream &bin_file, dt_int32 nxy, TVctr& matrix)
		// {
		// 	using T_o = Value_type<TVctr>;

		// 	vector<T_i> vect_r(nxy);
		// 	bin_file.read(reinterpret_cast<char*>(vect_r.data()), nxy * sizeof(T_i));

		// 	for(auto ik = 0; ik < nxy; ik++)
		// 	{
		// 		matrix[ik] = T_o(vect_r[ik]);
		// 	}
		// }

		// template <class T_i, class TVctr>
		// enable_if_cfloat<T_i, void>
		// read_mat_matrix(std::ifstream &bin_file, dt_int32 nxy, TVctr& matrix)
		// {
		// 	using T = Value_type_r<T_i>;
		// 	using T_o = Value_type_r<TVctr>;

		// 	vector<T> vect_r(nxy);
		// 	bin_file.read(reinterpret_cast<char*>(vect_r.data()), nxy * sizeof(T));
		// 	vector<T> vect_i(nxy);
		// 	bin_file.read(reinterpret_cast<char*>(vect_i.data()), nxy * sizeof(T));

		// 	auto matrix_ptr = reinterpret_cast<complex<T_o>*>(matrix.data());

		// 	for(auto ik = 0; ik < nxy; ik++)
		// 	{
		// 		matrix_ptr[ik].real(vect_r[ik]);
		// 		matrix_ptr[ik].imag(vect_i[ik]);
		// 	}
		// }

		// template <class TVctr>
		// void read_mat_binary_matrix(const char *filename, Grid_2d<Value_type_r<TVctr>>& grid_2d, TVctr& matrix)
		// {
		// 	dt_int32 nx, ny;
		// 	dt_float64 dx, dy;
		// 	dt_int32 type;

		// 	std::ifstream bin_file(filename, std::ofstream::binary);
		// 	bin_file.read(reinterpret_cast<char*>(&nx), sizeof(dt_int32));
		// 	bin_file.read(reinterpret_cast<char*>(&ny), sizeof(dt_int32));
		// 	bin_file.read(reinterpret_cast<char*>(&dx), sizeof(dt_float64));
		// 	bin_file.read(reinterpret_cast<char*>(&dy), sizeof(dt_float64));
		// 	bin_file.read(reinterpret_cast<char*>(&type), sizeof(dt_int32));

		// 	grid_2d.set_in_data(nx, ny, nx*dx, ny*dy);

		// 	auto nxy = nx*ny;
		// 	matrix.resize(nxy);

		// 	switch (type)
		// 	{
		// 	case 1:
		// 	{
		// 		read_mat_matrix<dt_float32>(bin_file, nxy, matrix);
		// 	}
		// 	break;
		// 	case 2:
		// 	{
		// 		read_mat_matrix<dt_float64>(bin_file, nxy, matrix);
		// 	}
		// 	break;
		// 	case 3:
		// 	{
		// 		read_mat_matrix<complex<dt_float32>>(bin_file, nxy, matrix);
		// 	}
		// 	break;
		// 	case 4:
		// 	{
		// 		read_mat_matrix<complex<dt_float64>>(bin_file, nxy, matrix);
		// 	}
		// 	break;
		// 	}
		// 	bin_file.close();
		// }
		// /************************ write matlab binary file ***********************/
		// template <class TVctr>
		// enable_if_real_vctr_cpu<TVctr, void>
		// 	write_mat_matrix(std::ofstream &bin_file, TVctr& matrix)
		// {
		// 	// using T = Value_type<TVctr>;

		// 	// bin_file.write(reinterpret_cast<char*>(matrix.data()), matrix.size()*sizeof(T));
		// }

		// template <class TVctr>
		// enable_if_cfloat_vctr_cpu<TVctr, void>
		// 	write_mat_matrix(std::ofstream &bin_file, TVctr& matrix)
		// {
		// 	// using T = Value_type_r<TVctr>;

		// 	// vector<T> vect_r(nxy);
		// 	// vector<T> vect_i(nxy);
		// 	// for(auto ik = 0; ik<nxy; ik++)
		// 	// {
		// 	// vect_r[ik] = matrix[ik].real();
		// 	// vect_i[ik] = matrix[ik].imag();
		// 	// }

		// 	// bin_file.write(reinterpret_cast<char*>(vect_r.data()), matrix.size()*sizeof(T));
		// 	// bin_file.write(reinterpret_cast<char*>(vect_i.data()), matrix.size()*sizeof(T));
		// }

		// template <class TVctr>
		// void write_mat_binary_matrix(const char *filename, Grid_2d<Value_type_r<TVctr>>& grid_2d, TVctr& matrix)
		// {
		// 	// dt_int32 type = matrix_type<TVctr>;

		// 	// std::ofstream bin_file(filename, std::ofstream::binary);
		// 	// bin_file.write(reinterpret_cast<char*>(&(grid_2d.nx)), sizeof(dt_int32));
		// 	// bin_file.write(reinterpret_cast<char*>(&(grid_2d.ny)), sizeof(dt_int32));
		// 	// bin_file.write(reinterpret_cast<char*>(&(grid_2d.drx)), sizeof(dt_float64));
		// 	// bin_file.write(reinterpret_cast<char*>(&(grid_2d.dry)), sizeof(dt_float64));
		// 	// bin_file.write(reinterpret_cast<char*>(&type), sizeof(dt_int32));

		// 	// switch (type)
		// 	// {
		// 	// case 1:
		// 	// {
		// 	// write_mat_matrix<dt_float32>(bin_file, matrix);
		// 	// }
		// 	// break;
		// 	// case 2:
		// 	// {
		// 	// write_mat_matrix<dt_float64>(bin_file, matrix);
		// 	// }
		// 	// break;
		// 	// case 3:
		// 	// {
		// 	// write_mat_matrix<complex<dt_float32>>(bin_file, matrix);
		// 	// }
		// 	// break;
		// 	// case 4:
		// 	// {
		// 	// write_mat_matrix<complex<dt_float64>>(bin_file, matrix);
		// 	// }
		// 	// break;
		// 	// }
		// 	// bin_file.close();
		// }

		// /************************ extract real vector form complex vector**********************/
		// template <class TVctr_c, class TVctr_r>
		// void from_complex_to_real(eShow_CData show_data, TVctr_c& cdata, TVctr_r &data)
		// {
		// 	switch (show_data)
		// 	{
		// 	case escd_creal:
		// 	{
		// 		for(auto ixy = 0; ixy < cdata.size(); ixy++)
		// 			data[ixy] = cdata[ixy].real();
		// 	}
		// 	break;
		// 	case escs_cimag:
		// 	{
		// 		for(auto ixy = 0; ixy < cdata.size(); ixy++)
		// 			data[ixy] = cdata[ixy].imag();
		// 	}
		// 	break;
		// 	case escd_cmod:
		// 	{
		// 		for(auto ixy = 0; ixy < cdata.size(); ixy++)
		// 			data[ixy] = abs(cdata[ixy]);
		// 	}
		// 	break;
		// 	case escd_cphase:
		// 	{
		// 		for(auto ixy = 0; ixy < cdata.size(); ixy++)
		// 			data[ixy] = arg(cdata[ixy]);
		// 	}
		// 	break;
		// 	}
		// }
	}

	namespace mt
	{
		// /***************************************************************************************/
		// template <class TGrid, class TVctr_r>
		// enable_if_vctr_cpu<TVctr_r, Value_type<TVctr_r>>
		// 	atom_cost_function(TGrid& grid_2d, const Atom_Sa<Value_type<TGrid>>& atom_Ip, TVctr_r& mx_i)
		// {
		// 	return cpu_detail::atom_cost_function<typename TVctr_r::value_type>(grid_2d, atom_Ip, mx_i);
		// }

		// template <class TGrid, class TVctr_r>
		// enable_if_vctr_cpu<TVctr_r, void>
		// 	subtract_atom(Stream_cpu& stream, TGrid& grid_2d, Vctr<Atom_Sa<Value_type<TGrid>>, edev_cpu>& atom_Ip, TVctr_r& mx_i)
		// {
		// 	if (stream.n_stream_act <= 0)
		// 	{
		// 		return;
		// 	}

		// 	for(auto istm = 0; istm < stream.n_stream_act; istm++)
		// 	{
		// 		stream[istm] = std::thread(std::bind(cpu_detail::subtract_atom<typename TVctr_r::value_type>, std::ref(stream), std::ref(grid_2d), std::ref(atom_Ip[istm]), std::ref(mx_i)));
		// 	}
		// 	stream.synchronize();
		// }

		// // Linear projected potential: V and zV
		// template <class TQ1, class TVAtom>
		// enable_if_host<TQ1, void>
		// 	linear_Vz(Stream_cpu& stream, eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, TQ1 &qz, TVAtom &vatom)
		// {
		// 	using TAtom = Value_type<TVAtom>;

		// 	if (stream.n_stream_act <= 0)
		// 	{
		// 		return;
		// 	}

		// 	auto thr_linear_Vz = [](const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, TQ1 &qz, TAtom &atom)
		// 	{
		// 		if (atom.charge == 0)
		// 		{
		// 			switch (atomic_pot_parm_typ)
		// 			{
		// 			case eappt_doyle_0_4:
		// 				cpu_detail::linear_Vz<eappt_doyle_0_4, 0, TAtom>(qz, atom);
		// 				break;
		// 			case eappt_peng_0_4:
		// 				cpu_detail::linear_Vz<eappt_peng_0_4, 0, TAtom>(qz, atom);
		// 				break;
		// 			case eappt_peng_0_12:
		// 				cpu_detail::linear_Vz<eappt_peng_0_12, 0, TAtom>(qz, atom);
		// 				break;
		// 			case eappt_kirkland_0_12:
		// 				cpu_detail::linear_Vz<eappt_kirkland_0_12, 0, TAtom>(qz, atom);
		// 				break;
		// 			case eappt_weickenmeier_0_12:
		// 				cpu_detail::linear_Vz<eappt_weickenmeier_0_12, 0, TAtom>(qz, atom);
		// 				break;
		// 			case eappt_lobato_0_12:
		// 				cpu_detail::linear_Vz<eappt_lobato_0_12, 0, TAtom>(qz, atom);
		// 				break;
		// 			}
		// 		}
		// 		else
		// 		{
		// 			switch (atomic_pot_parm_typ)
		// 			{
		// 			case eappt_doyle_0_4:
		// 				cpu_detail::linear_Vz<eappt_doyle_0_4, 1, TAtom>(qz, atom);
		// 				break;
		// 			case eappt_peng_0_4:
		// 				cpu_detail::linear_Vz<eappt_peng_0_4, 1, TAtom>(qz, atom);
		// 				break;
		// 			case eappt_peng_0_12:
		// 				cpu_detail::linear_Vz<eappt_peng_0_12, 1, TAtom>(qz, atom);
		// 				break;
		// 			case eappt_kirkland_0_12:
		// 				cpu_detail::linear_Vz<eappt_kirkland_0_12, 1, TAtom>(qz, atom);
		// 				break;
		// 			case eappt_weickenmeier_0_12:
		// 				cpu_detail::linear_Vz<eappt_weickenmeier_0_12, 1, TAtom>(qz, atom);
		// 				break;
		// 			case eappt_lobato_0_12:
		// 				cpu_detail::linear_Vz<eappt_lobato_0_12, 1, TAtom>(qz, atom);
		// 				break;
		// 			}
		// 		}
		// 	};

		// 	for(auto istm = 0; istm < stream.n_stream_act - 1; istm++)
		// 	{
		// 		stream[istm] = std::thread(std::bind(thr_linear_Vz, atomic_pot_parm_typ, std::ref(qz), std::ref(vatom[istm])));
		// 	}

		// 	thr_linear_Vz(atomic_pot_parm_typ, qz, vatom[stream.n_stream_act - 1]);

		// 	stream.synchronize();
		// }

		// // Get local interpolation coefficients
		// template <class TVAtom>
		// enable_if_host<typename TVAtom::value_type, void>
		// fcn_vd_2_coef_poly3(Stream_cpu& stream, TVAtom &vatom)
		// {
		// 	using TAtom = Value_type<TVAtom>;

		// 	if (stream.n_stream_act <= 0)
		// 	{
		// 		return;
		// 	}

		// 	for(auto istm = 0; istm < stream.n_stream_act - 1; istm++)
		// 	{
		// 		stream[istm] = std::thread(std::bind(cpu_detail::fcn_vd_2_coef_poly3<TAtom>, std::ref(vatom[istm])));
		// 	}

		// 	cpu_detail::fcn_vd_2_coef_poly3<TAtom>(vatom[stream.n_stream_act - 1]);

		// 	stream.synchronize();
		// }
	}

#endif