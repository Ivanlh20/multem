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
	template <eDevice dev>
	struct Timing;

	template<>
	struct Timing<e_device>
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

			float elapsed_ms()
			{
				//cudaDeviceSynchronize();
				cudaEventSynchronize(stop);

				float milliseconds = 0;
				cudaEventElapsedTime(&milliseconds, start, stop);

				return milliseconds;
			}

			float elapsed_s()
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


	template<>
	struct Timing<e_host>
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

			float elapsed_ms()
			{
				float milliseconds = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count()*1e-3;

				return milliseconds;
			}

			float elapsed_s()
			{
				return elapsed_ms()/1000;
			}

		private:
			std::chrono::high_resolution_clock::time_point start;
			std::chrono::high_resolution_clock::time_point stop;
	};
} // namespace mt

#endif