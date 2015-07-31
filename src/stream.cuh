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

#ifndef STREAM_H
#define STREAM_H

#include <thread>

#include <cuda.h>
#include <cuda_runtime.h>

namespace multem
{
	template<class T, eDevice dev>
	struct Stream;

	template<class T>
	struct Stream<T, e_Host>
	{
		public:
			using value_type = T;

			static const eDevice device = e_Host;

			Stream():nstream(0), n_act_stream(0), stream(nullptr){}

			~Stream(){ destroy(); nstream = 0; n_act_stream = 0; }

			int size() const
			{
				return nstream;
			}

			void resize(const int &new_size)
			{
				destroy();

				nstream = new_size;
				stream = new std::thread[nstream];
			}

			std::thread& operator[](const int i){ return stream[i]; }

			const std::thread& operator[](const int i) const { return stream[i]; }

			void synchronize()
			{
				destroy();

				stream = new std::thread[nstream];
			}

			int n_act_stream;
		private:
			int nstream;
			std::thread *stream;

			void destroy()
			{
				if(nstream == 0)
				{
					return;
				}

				for(auto i = 0; i < nstream; i++)
				{
					if(stream[i].joinable())
					{
						stream[i].join();
					}
				}

				delete [] stream;
			};
	};

	template<class T>
	struct Stream<T, e_Device>
	{
		public:
			using value_type = T;

			static const eDevice device = e_Device;

			Stream(): n_act_stream(0){}

			~Stream(){ destroy(); n_act_stream = 0; }

			int size() const
			{
				return stream.size();
			}

			void resize(const int &new_size)
			{
				destroy();

				stream.resize(new_size);

				for(auto i = 0; i < stream.size(); i++)
				{
					cudaStreamCreate(&(stream[i]));
				}
			}

			cudaStream_t& operator[](const int i){ return stream[i]; }

			const cudaStream_t& operator[](const int i) const { return stream[i]; }

			void synchronize()
			{
				cudaDeviceSynchronize();
			}

			int n_act_stream;
		private:
			std::vector<cudaStream_t> stream;

			void destroy()
			{
				if(stream.empty())
				{
					return;
				}

				cudaDeviceSynchronize();

				for(auto i = 0; i < stream.size(); i++)
				{
					cudaStreamDestroy(stream[i]);
				}
			}
	};
} // namespace multem

#endif