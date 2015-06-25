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

#ifdef _WIN32
#include <Windows.h>
#else
#include <sys/types.h>
#ifdef __APPLE__
#include <sys/sysctl.h>
#else
#include <sys/sysinfo.h>
#include <unistd.h>
#endif
#endif
#include <cstddef>

#include <thread>
#include <vector>
#include <algorithm>

#include "types.hpp"
#include "device_functions.cuh"

#include <cuda.h>
#include <cuda_runtime.h>

namespace multem
{
	template<eDevice dev>
	void memory_info(double &total, double &free);

	template<>
	void memory_info<e_Host>(double &total, double &free)
	{
#ifdef _WIN32
		MEMORYSTATUSEX status;
		status.dwLength = sizeof(status);
		GlobalMemoryStatusEx(&status);
		free = static_cast<double>(status.ullAvailPhys)/(1048576.0);
		total = static_cast<double>(status.ullTotalPhys)/(1048576.0);
#else // unix
		struct sysinfo memInfo;
		sysinfo (&memInfo);
		long long totalPhysMem = memInfo.totalram;
		totalPhysMem *= memInfo.mem_unit;
		long long physMemFree = memInfo.freeram;
		physMemFree *= memInfo.mem_unit;
		free = static_cast<double>(physMemFree); // check if division by 1MB=1048576 is necessary
		total = static_cast<double>(totalPhysMem); // check if division by 1MB=1048576 is necessary
#endif
	}

	template<>
	void memory_info<e_Device>(double &total, double &free)
	{
		free = total = 0;
		size_t free_t, total_t;
		if(cudaSuccess == cudaMemGetInfo(&free_t, &total_t))
		{
			free = static_cast<double>(free_t)/(1048576.0);
			total = static_cast<double>(total_t)/(1048576.0);
		}
	}
    
    // the free and total amount of memory available (Mb)
	inline
	void memory_info(double &host_free, double &host_total, double &device_free, double &device_total)
	{
		memory_info<e_Host>(host_total, host_free);

		memory_info<e_Device>(device_total, device_free);
	}

	template<eDevice dev>
	double get_free_memory();

	template<>
	double get_free_memory<e_Host>()
	{
		double total, free;
		memory_info<e_Host>(total, free);
		return free;
	}

	template<>
	double get_free_memory<e_Device>()
	{
		double total, free;
		memory_info<e_Device>(total, free);
		return free;
	}

	inline
	void get_device_properties(std::vector<Device_Properties> &device_properties)
	{
		device_properties.clear();

		if (!is_gpu_available())
		{
			return;
		}

		int device_count = 0;
		cudaGetDeviceCount(&device_count);

		device_properties.resize(device_count);
		for (auto idev = 0; idev < device_count; idev++)
		{
			cudaSetDevice(idev);
			cudaDeviceProp cuda_device_prop;
			cudaGetDeviceProperties(&cuda_device_prop, idev);

			device_properties[idev].id = idev;
			device_properties[idev].name = cuda_device_prop.name;
			device_properties[idev].compute_capability = 10*cuda_device_prop.major+cuda_device_prop.minor;
			memory_info<e_Device>(device_properties[idev].total_memory_size, device_properties[idev].free_memory_size);
		}

		auto compare_fn = [](const Device_Properties &a, const Device_Properties &b)->bool{ return a.compute_capability > b.compute_capability; };
		std::sort(device_properties.begin(), device_properties.end(), compare_fn);
	}

	inline
	void get_host_properties(Host_Properties &host_properties)
	{
#ifdef _WIN32
		SYSTEM_INFO siSysInfo;
		GetSystemInfo(&siSysInfo);
		host_properties.nprocessors = siSysInfo.dwNumberOfProcessors;
#elif __APPLE__
		int count = 0;
		size_t count_len = sizeof(count);
		sysctlbyname("hw.logicalcpu", &count, &count_len, NULL, 0);
                host_properties.nprocessors = count;
#else
		host_properties.nprocessors = sysconf(_SC_NPROCESSORS_ONLN);
#endif
		host_properties.nthreads = std::thread::hardware_concurrency();
		memory_info<e_Host>(host_properties.total_memory_size, host_properties.free_memory_size);
	}

} // namespace multem

#endif