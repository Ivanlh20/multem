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

#ifndef TIMING_H
#define TIMING_H

#include "types.cuh"
#include <chrono>

#ifdef __CUDACC__
	#include <cuda.h>
	#include <cuda_runtime.h>
#endif

namespace mt
{
	template <eDev Dev>
	struct Timing;

	template <>
	struct Timing<edev_gpu>
	{
		public:

			Timing()
			{
				cudaEventCreate(&start);
				cudaEventCreate(&stop);
			}

			~Timing()
			{
				free();
			}

			void tic()
			{
				cudaEventRecord(start);
			}

			void toc()
			{
				cudaEventRecord(stop);
			}

			dt_float32 elapsed_ms()
			{
				// cudaDeviceSynchronize();
				cudaEventSynchronize(stop);

				dt_float32 milliseconds = 0;
				cudaEventElapsedTime(&milliseconds, start, stop);

				return milliseconds;
			}

			dt_float32 elapsed_s()
			{
				return elapsed_ms()/1000;
			}

		private:
			cudaEvent_t start;
			cudaEvent_t stop;

			void free()
			{
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
			}
	};


	template <>
	struct Timing<edev_cpu>
	{
		public:
			void tic()
			{
				start = std::chrono::high_resolution_clock::now();
			}

			void toc()
			{
				stop = std::chrono::high_resolution_clock::now();
			}

			dt_float32 elapsed_ms()
			{
				dt_float32 milliseconds = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count()*1e-3;

				return milliseconds;
			}

			dt_float32 elapsed_s()
			{
				return elapsed_ms()/1000;
			}

		private:
			std::chrono::high_resolution_clock::time_point start;
			std::chrono::high_resolution_clock::time_point stop;
	};
}

#endif