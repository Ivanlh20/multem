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

namespace mt
{
	class System_Config
	{
	public:
		ePrecision precision;
		eDev device;							// eprc_float32 = 1, eprc_float64 = 2
		dt_int32 cpu_n_proc;					// number of Cores CPU
		dt_int32 cpu_n_thread;					// number of threads

		Vctr_cpu<dt_int32> gpu_device;			// gpu devices
		dt_int32 gpu_n_stream;					// number of streams

		dt_int32 gpu_n_avbl;					// number of selected gpus
		dt_int32 n_stream;						// number of streams
		dt_int32 idx_0;							// parameter position

		System_Config();

		System_Config(const System_Config &system_config);

		System_Config& operator=(const System_Config &system_config);

		void set_dep_var();

		void set_gpu(dt_int32 gpu_ind=-1);

		void set_gpu_by_ind(dt_int32 gpu_ind);

		dt_int32 get_n_sel_gpu();

 		dt_int32 get_sel_gpu();

		dt_bool is_cpu() const;

		dt_bool is_gpu() const;

		dt_bool is_float32() const;

		dt_bool is_float64() const;

		dt_bool is_float32_cpu() const;

		dt_bool is_float64_cpu() const;

		dt_bool is_float32_gpu() const;

		dt_bool is_float64_gpu() const;
	};
}

#include "../src/system_config.inl"