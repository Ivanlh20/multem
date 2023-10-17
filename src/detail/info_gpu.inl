/*
* This file is part of Multem.
* Copyright 2022 Ivan Lobato <Ivanlh20@gmail.com>
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

#pragma once

#include <cstddef>
#include <string>
#include <thread>
#include <vector>
#include <algorithm>

#include "const_enum.h"
#include "info_gpu.h"

#ifdef __CUDACC__
	#include <cuda.h>
	#include <cuda_runtime.h>
#endif

namespace mt
{
#ifdef __CUDACC__
	Gpu_Info::Gpu_Info(): id(0), name(""), comp_cap(0), 
			tot_mem(0), free_mem(0) {}

	/* device info */
	namespace dev_info
	{
		void gpu_tot_free_mem(dt_float64 &tot, dt_float64 &free, dt_int32 idx)
		{
			dt_float64 s_bytes_2_mb(c_bytes_2_mb);

			tot = free = 0;
			dt_uint64 free_t, tot_t;
			if (cudaSuccess == cudaMemGetInfo(&free_t, &tot_t))
			{
				free = static_cast<dt_float64>(free_t)/s_bytes_2_mb;
				tot = static_cast<dt_float64>(tot_t)/s_bytes_2_mb;
			}
		}

		dt_float64 gpu_free_mem()
		{
			dt_float64 tot, free;
			gpu_tot_free_mem(tot, free);
			return tot;
		}

		dt_float64 gpu_tot_mem()
		{
			dt_float64 tot=0, free=0;
			gpu_tot_free_mem(tot, free);
			return free;
		}

		dt_bool is_gpu_avbl()
		{
			dt_bool is_available = false;
			try
			{
				dt_int32 device_count = 0;
				cudaError_t error_id = cudaGetDeviceCount(&device_count);

				is_available = !((error_id != cudaSuccess)||(device_count == 0));
			}
			catch(...)
			{
				is_available = false;
			}

			return is_available;
		}

		dt_int32 gpu_n_avbl()
		{
			dt_int32 device_count = 0;
			cudaError_t error_id = cudaGetDeviceCount(&device_count);

			return (error_id != cudaSuccess)?0:device_count;
		}

		std::vector<Gpu_Info> gpu_info()
		{ 
			std::vector<Gpu_Info> info;

			if (!is_gpu_avbl())
			{
				return info;
			}

			dt_int32 device_count = 0;
			cudaGetDeviceCount(&device_count);

			info.resize(device_count);
			for(auto idev = 0; idev < device_count; idev++)
			{
				cudaDeviceProp cuda_device_prop;
				cudaGetDeviceProperties(&cuda_device_prop, idev);

				info[idev].id = idev;
				info[idev].name = cuda_device_prop.name;
				info[idev].comp_cap = 10*cuda_device_prop.major+cuda_device_prop.minor;
				gpu_tot_free_mem(info[idev].tot_mem, info[idev].free_mem);
			}

			// auto compare_fn = [](const Gpu_Info &a, const Gpu_Info &b)->dt_bool
			// { 
			// 	return a.comp_cap > b.comp_cap;
			// };
			// std::sort(info.begin(), info.end(), compare_fn);

			return info;
		}
	}
	#endif
}