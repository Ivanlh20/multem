/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef MEMORY_INFO_H
#define MEMORY_INFO_H

#include <Windows.h>
#include <cstddef>

#include "types.hpp"
#include <cuda.h>
#include <cuda_runtime.h>

namespace multem
{
	// the free and total amount of memory available (Mb)
	inline
	void memory_info(double &host_free, double &host_total, double &device_free, double &device_total)
	{
		MEMORYSTATUSEX status;
		status.dwLength = sizeof(status);
		GlobalMemoryStatusEx(&status);
		host_free = static_cast<double>(status.ullAvailPhys)/(1048576.0);
		host_total = static_cast<double>(status.ullTotalPhys)/(1048576.0);

		device_free = device_total = 0;
		size_t free, total;
		if(cudaSuccess==cudaMemGetInfo(&free, &total))
		{
			device_free = static_cast<double>(free)/(1048576.0);
			device_total = static_cast<double>(total)/(1048576.0);
		}
	}

	template<eDevice dev>
	void memory_info(double &free, double &total)
	{
		if(dev==Host)
		{
			MEMORYSTATUSEX status;
			status.dwLength = sizeof(status);
			GlobalMemoryStatusEx(&status);
			free = static_cast<double>(status.ullAvailPhys)/(1048576.0);
			total = static_cast<double>(status.ullTotalPhys)/(1048576.0);
		}
		else if(dev==Device)
		{
			free = total = 0;
			size_t free_t, total_t;
			if(cudaSuccess==cudaMemGetInfo(&free_t, &total_t))
			{
				free = static_cast<double>(free_t)/(1048576.0);
				total = static_cast<double>(total_t)/(1048576.0);
			}
		}
	}

} // namespace multem

#endif