/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef CGPU_INFO_H
	#define CGPU_INFO_H

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

	#include "const_enum.cuh"
	#include "cgpu_vctr.cuh"

	#ifdef __CUDACC__
		#include <cuda.h>
		#include <cuda_runtime.h>
	#endif

	namespace mt
	{
		class Cpu_Info
		{
		public:
			dt_int32 n_proc;		// number of cores
			dt_int32 n_threads;		// number of threads
			dt_float64 tot_mem;		// size in Mb
			dt_float64 free_mem;	// size in Mb

			Cpu_Info(): n_proc(0), n_threads(0), 
				tot_mem(0), free_mem(0) {}
		};

		class Gpu_Info
		{
		public:
			dt_int32 id;			// gpu id
			std::string name;		// gpu name
			dt_int32 comp_cap;		// compute capability
			dt_float64 tot_mem;		// size in Mb
			dt_float64 free_mem;	// size in Mb

			Gpu_Info(): id(0), name(""), comp_cap(0), 
				tot_mem(0), free_mem(0) {}
		};

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

			Cpu_Info cpu_info()
			{
				Cpu_Info info;

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

		#ifdef __CUDACC__
			void gpu_tot_free_mem(dt_float64 &tot, dt_float64 &free, dt_int32 idx = 0)
			{
				dt_float64 s_bytes_2_mb(c_bytes_2_mb);

				tot = free = 0;
			#ifdef __CUDACC__
				dt_uint64 free_t, tot_t;
				if (cudaSuccess == cudaMemGetInfo(&free_t, &tot_t))
				{
					free = static_cast<dt_float64>(free_t)/s_bytes_2_mb;
					tot = static_cast<dt_float64>(tot_t)/s_bytes_2_mb;
				}
			#endif
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
		#endif
		}

		/******************************* device configuration **********************************/
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

				System_Config(): precision(eprc_float64), device(edev_cpu), cpu_n_proc(1), 
					cpu_n_thread(1), gpu_n_stream(1), gpu_n_avbl(0), n_stream(1), idx_0(0) {};

				System_Config(const System_Config &system_config)
				{
					*this = system_config;
				}

				System_Config& operator=(const System_Config &system_config)
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

				void set_dep_var()
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

				void set_gpu(dt_int32 gpu_ind=-1)
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

				void set_gpu_by_ind(dt_int32 gpu_ind)
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

				dt_int32 get_n_sel_gpu()
				{
					return (dt_int32)gpu_device.size();
				}

 				dt_int32 get_sel_gpu()
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

				dt_bool is_cpu() const
				{
					return device == edev_cpu;
				}

				dt_bool is_gpu() const
				{
					return device == edev_gpu;
				}

				dt_bool is_float32() const
				{
					return precision == eprc_float32;
				}

				dt_bool is_float64() const
				{
					return precision == eprc_float64;
				}

				dt_bool is_float32_cpu() const
				{
					return is_float32() && is_cpu();
				}

				dt_bool is_float64_cpu() const
				{
					return is_float64() && is_cpu();
				}

				dt_bool is_float32_gpu() const
				{
					return is_float32() && is_gpu();
				}

				dt_bool is_float64_gpu() const
				{
					return is_float64() && is_gpu();
				}
		};

	}

#endif