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

#ifndef STREAM_H
#define STREAM_H

#include <thread>

#include <cuda.h>
#include <cuda_runtime.h>

namespace mt
{
	template <eDevice dev>
	struct Stream;

	template <>
	struct Stream<e_host>
	{
		public:
			static const eDevice device = e_host;

			std::mutex stream_mutex;

			Stream(): nx(0), ny(0), nxy(0), nstream(0), n_act_stream(0), 
			stream(nullptr){}

			Stream(int new_nstream): nx(0), ny(0), nxy(0), nstream(0), 
			n_act_stream(0), stream(nullptr)
			{
				resize(new_nstream);
			}

			~Stream(){ destroy(); nstream = 1; n_act_stream = 1; }

			int size() const
			{
				return nstream;
			}

			void resize(int new_nstream)
			{
				new_nstream = max(1, new_nstream);

				destroy();

				nstream = new_nstream;
				if(nstream > 1)
				{
					stream = new std::thread[nstream-1];
				}

				set_n_act_stream(size());
			}

			std::thread& operator[](const int i){ return stream[i]; }

			const std::thread& operator[](const int i) const { return stream[i]; }

			void synchronize()
			{
				destroy();

				if(nstream > 1)
				{
					stream = new std::thread[nstream-1];
				}
			}

			void set_n_act_stream(const int &new_n_act_stream)
			{
				n_act_stream = (new_n_act_stream<0)?0:min(size(), new_n_act_stream);
			}

			void set_grid(const int &nx_i, const int &ny_i)
			{
				nx = nx_i;
				ny = ny_i;
				nxy = nx*ny;
			}

			Range_2d get_range(const int &istream)
			{
				Range_2d range;
				
				int qnxy = nxy/n_act_stream;
				range.ixy_0 = istream*qnxy;
				range.ixy_e = (istream+1)*qnxy;

				int qnx = nx/n_act_stream;
				range.ix_0 = istream*qnx;
				range.ix_e = (istream+1)*qnx;
				range.iy_0 = 0;
				range.iy_e = ny;

				if(istream == n_act_stream-1)
				{
					range.ix_e += (nx - qnx*n_act_stream);
					range.ixy_e += (nxy - qnxy*n_act_stream);
				}
				return range;
			}

			Range_2d get_range_yx(const int &istream)
			{
				Range_2d range;
				
				int qnxy = nxy/n_act_stream;
				range.ixy_0 = istream*qnxy;
				range.ixy_e = (istream+1)*qnxy;

				int qny = ny/n_act_stream;
				range.iy_0 = istream*qny;
				range.iy_e = (istream+1)*qny;
				range.ix_0 = 0;
				range.ix_e = nx;

				if(istream == n_act_stream-1)
				{
					range.iy_e += (ny - qny*n_act_stream);
					range.ixy_e += (nxy - qnxy*n_act_stream);
				}
				return range;
			}

			template <class TFn, class... TArgs>
			void exec(TFn &fn, TArgs &...arg)
			{
				if(n_act_stream < 1)
				{
					return;
				}

				for(auto istream = 0; istream < n_act_stream-1; istream++)
				{
					stream[istream] = std::thread(std::bind(fn, get_range(istream), std::ref<TArgs>(arg)...));
				}

				fn(get_range(n_act_stream-1), std::ref<TArgs>(arg)...);

				synchronize();
			}

			template <class TFn, class... TArgs>
			void exec_matrix(TFn &fn, TArgs &...arg)
			{
				if(n_act_stream < 1)
				{
					return;
				}

				auto thr_fn_m = [&](const Range_2d &range, TArgs &...arg)
				{
					for(auto ix = range.ix_0; ix < range.ix_e; ix++)
					{
						for(auto iy = range.iy_0; iy < range.iy_e; iy++)
						{
							fn(ix, iy, arg...);
						}
					}
				};

				for(auto istream = 0; istream < n_act_stream-1; istream++)
				{
					stream[istream] = std::thread(std::bind(thr_fn_m, get_range(istream), std::ref<TArgs>(arg)...));
				}

				thr_fn_m(get_range(n_act_stream-1), std::ref<TArgs>(arg)...);

				synchronize();
			}

			template <class TFn, class... TArgs>
			void exec_vector(TFn &fn, TArgs &...arg)
			{
				if(n_act_stream < 1)
				{
					return;
				}

				auto thr_fn_v = [&](const Range_2d &range, TArgs &...arg)
				{
					for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
					{
						fn(ixy, arg...);
					}
				};

				for(auto istream = 0; istream < n_act_stream-1; istream++)
				{
					stream[istream] = std::thread(std::bind(thr_fn_v, get_range(istream), std::ref<TArgs>(arg)...));
				}

				thr_fn_v(get_range(n_act_stream-1), std::ref<TArgs>(arg)...);

				synchronize();
			}

			int n_act_stream;
		private:
			int nx;
			int ny;
			int nxy;

			int nstream;
			std::thread *stream;

			void destroy()
			{
				if(nstream < 0)
				{
					return;
				}

				for(auto i = 0; i < nstream-1; i++)
				{
					if(stream[i].joinable())
					{
						stream[i].join();
					}
				}

				delete [] stream;
			};
	};

	template <>
	struct Stream<e_device>
	{
		public:
			static const eDevice device = e_device;

			Stream(): nx(0), ny(0), nxy(0), n_act_stream(0){}

			Stream(int new_nstream): nx(0), ny(0), nxy(0), n_act_stream(0)
			{
				resize(new_nstream);
			}

			~Stream(){ destroy(); n_act_stream = 0; }

			int size() const
			{
				return static_cast<int>(stream.size());
			}

			void resize(int new_nstream)
			{
				new_nstream = max(1, new_nstream);

				destroy();

				stream.resize(new_nstream);

				for(auto i = 0; i < stream.size(); i++)
				{
					cudaStreamCreate(&(stream[i]));
				}

				set_n_act_stream(size());
			}

			cudaStream_t& operator[](const int i){ return stream[i]; }

			const cudaStream_t& operator[](const int i) const { return stream[i]; }

			void synchronize()
			{
				cudaDeviceSynchronize();
			}

			void set_n_act_stream(const int &new_n_act_stream)
			{
				n_act_stream = (new_n_act_stream<0)?0:min(size(), new_n_act_stream);
			}

			void set_grid(const int &nx_i, const int &ny_i)
			{
				nx = nx_i;
				ny = ny_i;
				nxy = nx*ny;
			}

			Range_2d get_range(const int &istream)
			{
				Range_2d range;
				
				int qnxy = nxy/n_act_stream;
				range.ixy_0 = istream*qnxy;
				range.ixy_e = (istream+1)*qnxy;

				int qnx = nx/n_act_stream;
				range.ix_0 = istream*qnx;
				range.ix_e = (istream+1)*qnx;
				range.iy_0 = 0;
				range.iy_e = ny;

				if(istream == n_act_stream-1)
				{
					range.ix_e += (nx - qnx*n_act_stream);
					range.ixy_e += (nxy - qnxy*n_act_stream);
				}
				return range;
			}

			template <class TFn>
			void exec(TFn &fn)
			{
				// for(auto istream = 0; istream < n_act_stream; istream++)
				// {
				// 	stream[istream] = std::thread(fn, get_range(istream));
				// }
				// synchronize();
			}

			int n_act_stream;

			std::vector<cudaStream_t> stream;
		private:
			int nx;
			int ny;
			int nxy;

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
} // namespace mt

#endif