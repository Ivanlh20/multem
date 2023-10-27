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

#ifdef _WIN32
	#include <Windows.h>
#elif defined __APPLE__
	#include <sys/types.h>
	#include <sys/sysctl.h>
	#include <mach/mach.h>
	#include <stdint.h>
	#include <unistd.h>
#else
	#include <sys/types.h>
	#include <sys/sysinfo.h>
	#include <unistd.h>
#endif

#include <cstddef>
#include <string>
#include <thread>
#include <vector>
#include <algorithm>

#include "const_enum.h"
#include "info_cpu.h"

namespace mt
{
	Cpu_Info::Cpu_Info(): n_proc(0), n_threads(0), 
	tot_mem(0), free_mem(0) {}

	/* device info */
	namespace dev_info
	{
		void cpu_tot_free_mem(dt_float64 &tot, dt_float64 &free)
		{
			const dt_float64 s_bytes_2_mb(c_bytes_2_mb);

		#if defined(_WIN32)
			MEMORYSTATUSEX status;
			status.dwLength = sizeof(status);
			GlobalMemoryStatusEx(&status);
			free = static_cast<dt_float64>(status.ullAvailPhys)/s_bytes_2_mb;
			tot = static_cast<dt_float64>(status.ullTotalPhys)/s_bytes_2_mb;

		#elif defined(__APPLE__)
			dt_int32 mib[2];
			mib[0] = CTL_HW;
			mib[1] = HW_MEMSIZE;
					uint64_t physicalMem = 0;
			dt_uint64 returnSize = sizeof(physicalMem);
			sysctl(mib, 2, &physicalMem, &returnSize, NULL, 0);
			tot = returnSize/s_bytes_2_mb;

			task_t targetTask = mach_task_self();
			struct task_basic_info ti;
			mach_msg_type_number_t count = TASK_BASIC_INFO_64_COUNT;
			task_info(targetTask, TASK_BASIC_INFO_64, (task_info_t) &ti, &count);
			free = ti.resident_size/s_bytes_2_mb;

		#else // linux
			struct sysinfo memInfo;
			sysinfo (&memInfo);
			long long totalPhysMem = memInfo.totalram;
			totalPhysMem *= memInfo.mem_unit;
			long long physMemFree = memInfo.freeram;
			physMemFree *= memInfo.mem_unit;
			free = static_cast<dt_float64>(physMemFree)/s_bytes_2_mb;		// check if division by 1MB = 1048576 is necessary
			tot = static_cast<dt_float64>(totalPhysMem)/s_bytes_2_mb;		// check if division by 1MB = 1048576 is necessary
		#endif
		}

		dt_float64 cpu_free_mem()
		{
			dt_float64 tot, free;
			cpu_tot_free_mem(tot, free);
			return free;
		}
		
		dt_float64 cpu_tot_mem()
		{
			dt_float64 tot, free;
			cpu_tot_free_mem(tot, free);
			return tot;
		}

		mt::Cpu_Info cpu_info()
		{
			mt::Cpu_Info info;

		#ifdef _WIN32
			SYSTEM_INFO siSysInfo;
			GetSystemInfo(&siSysInfo);
			info.n_proc = siSysInfo.dwNumberOfProcessors;

		#elif defined(__APPLE__)
			dt_int32 count = 0;
			dt_uint64 count_len = sizeof(count);
			sysctlbyname("hw.logicalcpu", &count, &count_len, NULL, 0);
			info.n_proc = count;

		#else // linux
			info.n_proc = sysconf(_SC_NPROCESSORS_ONLN);

		#endif

			info.n_threads = std::thread::hardware_concurrency();
			cpu_tot_free_mem(info.tot_mem, info.free_mem);

			return info;
		}
	}
}