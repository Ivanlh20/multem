/*
 * This file is part of MULTEM.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef MEMORY_INFO_H
#define MEMORY_INFO_H

#ifdef _WIN32
	#include <Windows.h>
#else
	#include <sys/types.h>
	#ifdef __APPLE__
		#include <sys/sysctl.h>
		#include <mach/mach.h>
		#include <stdint.h>
		#include <unistd.h>
	#else
		#include <sys/sysinfo.h>
		#include <unistd.h>
	#endif
#endif
#include <cstddef>

#include <thread>
#include <vector>
#include <algorithm>

#include "types.cuh"

#include <cuda.h>
#include <cuda_runtime.h>

namespace mt
{
	bool is_gpu_available();

	template <eDevice dev>
	void memory_info(double &total, double &free);

	template <>
	inline
	void memory_info<e_host>(double &total, double &free)
	{
#if defined(_WIN32)
		MEMORYSTATUSEX status;
		status.dwLength = sizeof(status);
		GlobalMemoryStatusEx(&status);
		free = static_cast<double>(status.ullAvailPhys)/(1048576.0);
		total = static_cast<double>(status.ullTotalPhys)/(1048576.0);
#elif defined(__APPLE__)
		int mib[2];
		mib[0] = CTL_HW;
		mib[1] = HW_MEMSIZE;
				uint64_t physicalMem = 0;
		size_t returnSize = sizeof(physicalMem);
		sysctl(mib, 2, &physicalMem, &returnSize, NULL, 0);
				total = returnSize/(1048576.0);

		task_t targetTask = mach_task_self();
		struct task_basic_info ti;
		mach_msg_type_number_t count = TASK_BASIC_INFO_64_COUNT;
		task_info(targetTask, TASK_BASIC_INFO_64, (task_info_t) &ti, &count);
				free = ti.resident_size/(1048576.0);

#else // linux
		struct sysinfo memInfo;
		sysinfo (&memInfo);
		long long totalPhysMem = memInfo.totalram;
		totalPhysMem *= memInfo.mem_unit;
		long long physMemFree = memInfo.freeram;
		physMemFree *= memInfo.mem_unit;
		free = static_cast<double>(physMemFree)/(1048576.0);	// check if division by 1MB = 1048576 is necessary
		total = static_cast<double>(totalPhysMem)/(1048576.0);	// check if division by 1MB = 1048576 is necessary
#endif
	}

	template <>
	inline
	void memory_info<e_device>(double &total, double &free)
	{
		free = total = 0;
		size_t free_t, total_t;
		if(cudaSuccess == cudaMemGetInfo(&free_t, &total_t))
		{
			free = static_cast<double>(free_t)/(1048576.0);
			total = static_cast<double>(total_t)/(1048576.0);
		}
	}

	template <eDevice dev>
	double get_free_memory(){ return 0; }

	template <>
	inline
	double get_free_memory<e_host>()
	{
		double total, free;
		memory_info<e_host>(total, free);
		return free;
	}

	template <>
	inline
	double get_free_memory<e_device>()
	{
		double total, free;
		memory_info<e_device>(total, free);
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
			cudaDeviceProp cuda_device_prop;
			cudaGetDeviceProperties(&cuda_device_prop, idev);

			device_properties[idev].id = idev;
			device_properties[idev].name = cuda_device_prop.name;
			device_properties[idev].compute_capability = 10*cuda_device_prop.major+cuda_device_prop.minor;
			memory_info<e_device>(device_properties[idev].total_memory_size, device_properties[idev].free_memory_size);
		}

		//auto compare_fn = [](const Device_Properties &a, const Device_Properties &b)->bool
		//{ 
		//	return a.compute_capability > b.compute_capability; 
		//};
		//std::sort(device_properties.begin(), device_properties.end(), compare_fn);
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
		memory_info<e_host>(host_properties.total_memory_size, host_properties.free_memory_size);
	}

} // namespace mt

#endif