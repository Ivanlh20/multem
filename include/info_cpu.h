/*
* This file is part of Multem.
* Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
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

namespace mt
{
	class Cpu_Info
	{
	public:
		dt_int32 n_proc;		// number of cores
		dt_int32 n_threads;		// number of threads
		dt_float64 tot_mem;		// size in Mb
		dt_float64 free_mem;	// size in Mb

		explicit Cpu_Info();
	};

	/* device info */
	namespace dev_info
	{
		void cpu_tot_free_mem(::dt_float64 &tot, ::dt_float64 &free);

		::dt_float64 cpu_free_mem();
		
		::dt_float64 cpu_tot_mem();

		mt::Cpu_Info cpu_info();
	}
}

#include "../src/info_cpu.inl"