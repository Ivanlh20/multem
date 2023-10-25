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

#include "const_enum.h"
#include "vctr_cpu.h"
#include "info_cpu.h"
#include "info_gpu.h"

#include "system_config.h"

namespace mt
{
	System_Config::System_Config(): precision(eprc_float64), device(edev_cpu), cpu_n_proc(1), 
		cpu_n_thread(1), gpu_n_stream(1), gpu_n_avbl(0), n_stream(1), idx_0(0) {}

	System_Config::System_Config(const System_Config &system_config)
	{
		*this = system_config;
	}

	System_Config& System_Config::operator=(const System_Config &system_config)
	{
		precision = system_config.precision;
		device = system_config.device;
		cpu_n_proc = system_config.cpu_n_proc;
		cpu_n_thread = system_config.cpu_n_thread;
		gpu_device = system_config.gpu_device;
		gpu_n_stream = system_config.gpu_n_stream;
		n_stream = system_config.n_stream;
		idx_0 = system_config.idx_0;
		gpu_n_avbl = system_config.gpu_n_avbl;

		return *this;
	}

	void System_Config::set_dep_var()
	{
		// check precision
		if (!(is_float32() || is_float64()))
		{
			precision = eprc_float32;
		}

		// check cpu or gpu
		if (!(is_cpu() || is_gpu()))
		{
			device = edev_cpu;
		}
		if (is_gpu())
		{
			#ifdef __CUDACC__
				if (!dev_info::is_gpu_avbl())
				{
					device = edev_cpu;
					gpu_n_avbl = 0;
				}
				else
				{
					cpu_n_thread = 1;
					gpu_n_avbl = dev_info::gpu_n_avbl();
					for(auto gpu_ind=0; gpu_ind<gpu_device.size(); gpu_ind++)
					{
						gpu_device[gpu_ind] = min(max(0, gpu_device[gpu_ind]), gpu_n_avbl-1);
					}
				}
			#endif
		}

		cpu_n_thread = max(1, cpu_n_thread);
		gpu_n_stream = max(1, gpu_n_stream);
		n_stream = (is_cpu())?cpu_n_thread:gpu_n_stream;
	}

	void System_Config::set_gpu(dt_int32 gpu_ind)
	{
		if (is_gpu())
		{	
		#ifdef __CUDACC__
			gpu_ind = (gpu_ind<0)?gpu_device[0]:min(max(0, gpu_ind), gpu_n_avbl-1);
			cudaSetDevice(gpu_ind);
		#endif
		}
		else
		{
			device = edev_cpu;
		}
	}

	void System_Config::set_gpu_by_ind(dt_int32 gpu_ind)
	{
		if (is_gpu())
		{
		#ifdef __CUDACC__
			auto gpu_n_req = dt_int32(gpu_device.size());
			gpu_ind = min(max(0, gpu_ind), gpu_n_req-1);
			cudaSetDevice(gpu_device[gpu_ind]);
		#endif
		}
		else
		{
			device = edev_cpu;
		}
	}

	dt_int32 System_Config::get_n_sel_gpu()
	{
		return (dt_int32)gpu_device.size();
	}

 	dt_int32 System_Config::get_sel_gpu()
	{
		dt_int32 idx_dev = -1;
		if (is_gpu())
		{
		#ifdef __CUDACC__
			cudaGetDevice(&idx_dev);
		#endif
		}

		return idx_dev;
	}

	dt_bool System_Config::is_cpu() const
	{
		return device == edev_cpu;
	}

	dt_bool System_Config::is_gpu() const
	{
		return device == edev_gpu;
	}

	dt_bool System_Config::is_float32() const
	{
		return precision == eprc_float32;
	}

	dt_bool System_Config::is_float64() const
	{
		return precision == eprc_float64;
	}

	dt_bool System_Config::is_float32_cpu() const
	{
		return is_float32() && is_cpu();
	}

	dt_bool System_Config::is_float64_cpu() const
	{
		return is_float64() && is_cpu();
	}

	dt_bool System_Config::is_float32_gpu() const
	{
		return is_float32() && is_gpu();
	}

	dt_bool System_Config::is_float64_gpu() const
	{
		return is_float64() && is_gpu();
	}
}