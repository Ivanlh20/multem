/*
 * This file is part of Multem.
 * Copyright 2022 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * Multem is destroy software: you can redistribute it and/or modify
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

#ifndef MEMCPY_H
	#define MEMCPY_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include <typeinfo>

	#include "const_enum.h"
	#include "math_mt.h"
	#include "type_traits_gen.h"

	#ifdef __CUDACC__
		#include <cuda.h>
		#include <cuda_runtime.h>
	#endif

	/* gpu data cast */
	namespace mt
	{
		namespace gpu_detail
		{
		#ifdef __CUDACC__
			template <class Td, class Ts>
			__global__ void fcn_data_typ_cast(Td* data_dst, Ts* data_src, dt_int32 n_data)
			{
				FOR_IX_1DC(n_data)
				{
					data_dst[ix] = Td(data_src[ix]);
				}
			}

			template <class Td, class Ts>
			__global__ void fcn_real_gpu_gpu(Td* data_dst, Ts* data_src, dt_int32 n_data)
			{
				FOR_IX_1DC(n_data)
				{
					data_dst[ix] = Td(data_src[ix].real());
				}
			}

			template <class Td, class Ts>
			__global__ void fcn_imag_gpu_gpu(Td* data_dst, Ts* data_src, dt_int32 n_data)
			{
				FOR_IX_1DC(n_data)
				{
					data_dst[ix] = Td(data_src[ix].imag());
				}
			}
		#endif
		}
	}

	/* data copy */
	namespace mt
	{
		/******************************* (dst, src): cpu -> cpu ********************************/
		template <class Td, class Ts>
		void memcpy_cpu_cpu(Td* pcpu_dst, const Ts* pcpu_src, dt_uint64 n_size, Ts* pcpu_jk = nullptr)
		{
			for (dt_uint64 ik = 0; ik < n_size; ik++)
			{
				pcpu_dst[ik] = Td(pcpu_src[ik]);
			}
		}

	#ifdef __CUDACC__
		/******************************* (dst, src): gpu -> cpu ********************************/
		template <class Td, class Ts>
		enable_if_same_decay<Td, Ts, void>
		memcpy_gpu_cpu(Td* pcpu_dst, const Ts* pgpu_src, dt_uint64 n_size, Ts* pcpu_jk = nullptr)
		{
			const auto size_bytes = n_size*dt_uint64(sizeof(Ts));

			cudaMemcpy(pcpu_dst, pgpu_src, size_bytes, cudaMemcpyDeviceToHost);
		}		
		
		template <class Td, class Ts>
		enable_if_diff_decay<Td, Ts, void>
		memcpy_gpu_cpu(Td* pcpu_dst, const Ts* pgpu_src, dt_uint64 n_size, Ts* pcpu_jk = nullptr)
		{
			const auto size_bytes = n_size*dt_uint64(sizeof(Ts));

			if (pcpu_jk == nullptr)
			{
				Ts* pcpu_t = new Ts[n_size];
				cudaMemcpy(pcpu_t, pgpu_src, size_bytes, cudaMemcpyDeviceToHost);
				memcpy_cpu_cpu(pcpu_dst, pcpu_t, n_size);
				delete[] pcpu_t;
			}
			else
			{
				cudaMemcpy(pcpu_jk, pgpu_src, size_bytes, cudaMemcpyDeviceToHost);
				memcpy_cpu_cpu(pcpu_dst, pcpu_jk, n_size);
			}
		}

		/******************************* (dst, src): gpu -> gpu ********************************/
		// (dst, src): gpu -> gpu
		template <class Td, class Ts>
		enable_if_same_decay<Td, Ts, void>
		memcpy_gpu_gpu(Td* pgpu_dst, const Ts* pgpu_src, dt_uint64 n_size, Td* pcpu_jk = nullptr)
		{
			auto size_bytes = n_size*dt_uint64(sizeof(Td));
			cudaMemcpy(pgpu_dst, pgpu_src, size_bytes, cudaMemcpyDeviceToDevice);
		}		
		
		template <class Td, class Ts>
		enable_if_diff_decay<Td, Ts, void>
		memcpy_gpu_gpu(Td* pgpu_dst, const Ts* pgpu_src, dt_uint64 n_size, Td* pcpu_jk = nullptr)
		{
			// dt_int32 numSMs;
			// cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, devId);

			auto grid = fcn_cdg_1d(n_size);
			grid.x = min(128, grid.x);

			gpu_detail::fcn_data_typ_cast<Td, Ts><<<grid, c_thr_1d>>>(pgpu_dst, pgpu_src, dt_int32(n_size));
		}

		/******************************* (dst, src): cpu -> gpu ********************************/
		template <class Td, class Ts>
		enable_if_same_decay<Td, Ts, void>
		memcpy_cpu_gpu(Td* pgpu_dst, const Ts* pcpu_src, dt_uint64 n_size, Td* pcpu_jk = nullptr)
		{
			const auto size_bytes = n_size*dt_uint64(sizeof(Td));

			cudaMemcpy(pgpu_dst, pcpu_src, size_bytes, cudaMemcpyHostToDevice);
		}

		template <class Td, class Ts>
		enable_if_diff_decay<Td, Ts, void>
		memcpy_cpu_gpu(Td* pgpu_dst, const Ts* pcpu_src, dt_uint64 n_size, Td* pcpu_jk = nullptr)
		{
			const auto size_bytes = n_size*dt_uint64(sizeof(Td));

			if (pcpu_jk == nullptr)
			{
				Td* pcpu_t = new Td[n_size];
				memcpy_cpu_cpu(pcpu_t, pcpu_src, n_size);
				cudaMemcpy(pgpu_dst, pcpu_t, size_bytes, cudaMemcpyHostToDevice);
				delete[] pcpu_t;
			}
			else
			{
				memcpy_cpu_cpu(pcpu_jk, pcpu_src, n_size);
				cudaMemcpy(pgpu_dst, pcpu_jk, size_bytes, cudaMemcpyHostToDevice);
			}
		}

	#endif

	}

	/* complex real data copy */
	namespace mt
	{
		// (dst, src): cpu -> cpu
		template <class Td, class Ts>
		enable_if_cmplx<Ts, void>
		memcpy_real_cpu_cpu(Td* pcpu_dst, const Ts* pcpu_src, dt_uint64 n_size, Ts* pcpu_jk = nullptr)
		{
			if ((void*)pcpu_dst == (void*)pcpu_src) 
				return;

			for (dt_uint64 ik = 0; ik < n_size; ik++)
			{
				pcpu_dst[ik] = Td(pcpu_src[ik].real());
			}
		}

	#ifdef __CUDACC__
		// (dst, src): gpu -> cpu
		template <class Td, class Ts>
		enable_if_cmplx<Ts, void>
		memcpy_real_gpu_cpu(Td* pcpu_dst, const Ts* pgpu_src, dt_uint64 n_size, Ts* pcpu_jk = nullptr)
		{
			if ((void*)pcpu_dst == (void*)pgpu_src) 
				return;

			auto size_bytes = n_size*dt_uint64(sizeof(Ts));

			if (pcpu_jk == nullptr)
			{
				Ts* pcpu_t = new Ts[n_size];
				cudaMemcpy(pcpu_t, pgpu_src, size_bytes, cudaMemcpyDeviceToHost);
				memcpy_real_cpu_cpu(pcpu_dst, pcpu_t, n_size);
				delete[] pcpu_t;
			}
			else
			{
				cudaMemcpy(pcpu_jk, pgpu_src, size_bytes, cudaMemcpyDeviceToHost);
				memcpy_real_cpu_cpu(pcpu_dst, pcpu_jk, n_size);
			}
		}

		// (dst, src): gpu -> gpu
		template <class Td, class Ts>
		enable_if_cmplx<Ts, void>
		memcpy_real_gpu_gpu(Td* pgpu_dst, const Ts* pgpu_src, dt_uint64 n_size, Td* pcpu_jk = nullptr)
		{
			if ((void*)pgpu_dst == (void*)pgpu_src) 
				return;

			auto grid = fcn_cdg_1d(n_size);
			grid.x = min(128, grid.x);

			gpu_detail::fcn_real_gpu_gpu<<<grid, c_thr_1d>>>(pgpu_dst, pgpu_src, dt_int32(n_size));
		}

		// (dst, src): cpu -> gpu
		template <class Td, class Ts>
		enable_if_cmplx<Ts, void>
		memcpy_real_cpu_gpu(Td* pgpu_dst, const Ts* pcpu_src, dt_uint64 n_size, Td* pcpu_jk = nullptr)
		{
			if ((void*)pgpu_dst == (void*)pcpu_src) 
				return;

			auto size_bytes = n_size*dt_uint64(sizeof(Td));

			if (pcpu_jk == nullptr)
			{
				Td* pcpu_t = new Td[n_size];
				memcpy_real_cpu_cpu(pcpu_t, pcpu_src, n_size);
				cudaMemcpy(pgpu_dst, pcpu_t, size_bytes, cudaMemcpyHostToDevice);
				delete[] pcpu_t;
			}
			else
			{
				memcpy_real_cpu_cpu(pcpu_jk, pcpu_src, n_size);
				cudaMemcpy(pgpu_dst, pcpu_jk, size_bytes, cudaMemcpyHostToDevice);
			}
		}

	#endif

	}
	
	/* complex imag data copy */
	namespace mt
	{
		// (dst, src): cpu -> cpu
		template <class Td, class Ts>
		enable_if_cmplx<Ts, void>
		memcpy_imag_cpu_cpu(Td* pcpu_dst, const Ts* pcpu_src, dt_uint64 n_size, Ts* pcpu_jk = nullptr)
		{
			if ((void*)pcpu_dst == (void*)pcpu_src) 
				return;

			for (dt_uint64 ik = 0; ik < n_size; ik++)
			{
				pcpu_dst[ik] = Td(pcpu_src[ik].imag());
			}
		}

	#ifdef __CUDACC__
		// (dst, src): gpu -> cpu
		template <class Td, class Ts>
		enable_if_cmplx<Ts, void>
		memcpy_imag_gpu_cpu(Td* pcpu_dst, const Ts* pgpu_src, dt_uint64 n_size, Ts* pcpu_jk = nullptr)
		{
			if ((void*)pcpu_dst == (void*)pgpu_src) 
				return;

			auto size_bytes = n_size*dt_uint64(sizeof(Ts));

			if (pcpu_jk == nullptr)
			{
				Ts* pcpu_t = new Ts[n_size];
				cudaMemcpy(pcpu_t, pgpu_src, size_bytes, cudaMemcpyDeviceToHost);
				memcpy_imag_cpu_cpu(pcpu_dst, pcpu_t, n_size);
				delete[] pcpu_t;
			}
			else
			{
				cudaMemcpy(pcpu_jk, pgpu_src, size_bytes, cudaMemcpyDeviceToHost);
				memcpy_imag_cpu_cpu(pcpu_dst, pcpu_jk, n_size);
			}
		}

		// (dst, src): gpu -> gpu
		template <class Td, class Ts>
		enable_if_cmplx<Ts, void>
		memcpy_imag_gpu_gpu(Td* pgpu_dst, const Ts* pgpu_src, dt_uint64 n_size, Td* pcpu_jk = nullptr)
		{
			if ((void*)pgpu_dst == (void*)pgpu_src) 
				return;

			auto grid = fcn_cdg_1d(n_size);
			grid.x = min(128, grid.x);

			gpu_detail::fcn_imag_gpu_gpu<<<grid, c_thr_1d>>>(pgpu_dst, pgpu_src, dt_int32(n_size));
		}

		// (dst, src): cpu -> gpu
		template <class Td, class Ts>
		enable_if_cmplx<Ts, void>
		memcpy_imag_cpu_gpu(Td* pgpu_dst, const Ts* pcpu_src, dt_uint64 n_size, Td* pcpu_jk = nullptr)
		{
			if ((void*)pgpu_dst == (void*)pcpu_src) 
				return;

			auto size_bytes = n_size*dt_uint64(sizeof(Td));

			if (pcpu_jk == nullptr)
			{
				Td* pcpu_t = new Td[n_size];
				memcpy_imag_cpu_cpu(pcpu_t, pcpu_src, n_size);
				cudaMemcpy(pgpu_dst, pcpu_t, size_bytes, cudaMemcpyHostToDevice);
				delete[] pcpu_t;
			}
			else
			{
				memcpy_imag_cpu_cpu(pcpu_jk, pcpu_src, n_size);
				cudaMemcpy(pgpu_dst, pcpu_jk, size_bytes, cudaMemcpyHostToDevice);
			}
		}

	#endif

	}

#endif