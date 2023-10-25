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
#include <vector>

namespace mt
{
#ifdef __CUDACC__
	class Gpu_Info
	{
	public:
		dt_int32 id;			// gpu id
		std::string name;		// gpu name
		dt_int32 comp_cap;		// compute capability
		dt_float64 tot_mem;		// size in Mb
		dt_float64 free_mem;	// size in Mb

		Gpu_Info();
	};

	/* device info */
	namespace dev_info
	{
		void gpu_tot_free_mem(dt_float64 &tot, dt_float64 &free, dt_int32 idx = 0);

		dt_float64 gpu_free_mem();

		dt_float64 gpu_tot_mem();

		dt_bool is_gpu_avbl();

		dt_int32 gpu_n_avbl();

		std::vector<Gpu_Info> gpu_info();
	}
#endif
}

#include "../src/info_gpu.inl"