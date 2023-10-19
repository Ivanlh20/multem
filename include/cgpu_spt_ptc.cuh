/*
 * This file is part of Multem.
 * Copyright 2022 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef CGPU_SPT_PTC_H
	#define CGPU_SPT_PTC_H
#include <mex.h>
	#include "math_mt.h"
	#include "types.cuh"
	#include "type_traits_gen.h"
	#include "cgpu_stream.cuh"

	#ifdef __CUDACC__
		#include <cuda.h>
		#include <cuda_runtime.h>
		#include <cufft.h>
	#endif

	#include "in_classes.cuh"
	//#include "fcns_gpu.h"
	//#include "fcns_cpu.h"
	#include "particles.cuh"

	#ifdef __CUDACC__
		#include "fcns_gpu.h"
	#endif

	/* gpu detail */
	namespace mt
	{
		namespace gpu_detail
		{
			template <class T, eDim Dim, eFcn_typ Fcn_typ>
			__global__ void fcn_ptc_s_fcn_eval(Grid_2d<T> grid_2d, pVctr_gpu_32<Ptc_s_fcn_xd<T, Dim, Fcn_typ>> ptc, pVctr_gpu_32<T> mx_o)
			{
				const auto fcn_ptc = ptc.m_data[blockIdx.x];

				dt_int32 ix_t = threadIdx.y;
				while (ix_t < fcn_ptc.nx)
				{
					const dt_int32 ix = fcn_ptc.ix_0 + ix_t;

					dt_int32 iy_t = threadIdx.x;
					while (iy_t < fcn_ptc.ny)
					{
						const dt_int32 iy = fcn_ptc.iy_0 + iy_t;
						const auto r2 = grid_2d.r2(ix, iy, fcn_ptc.r);

						if (r2 < fcn_ptc.r_max_2)
						{
							const auto ixy = grid_2d.sub_2_ind_pbc(ix, iy);
							atomicAdd(&(mx_o.m_data[ixy]), fcn_ptc.eval_r2(r2));
						}

						iy_t += blockDim.x;
					}
					ix_t += blockDim.y;
				}
			}
		}
	}

	/* cpu detail */
	namespace mt
	{
		namespace cpu_detail
		{
			template <class T, eDim Dim, eFcn_typ Fcn_typ>
			void fcn_ptc_s_fcn_eval(Stream<edev_cpu>& stream, Grid_2d<T>& grid_2d, Ptc_s_fcn_xd<T, Dim, Fcn_typ>& fcn_ptc, Vctr_cpu<T>& mx_o)
			{
				for (auto ix_0 = 0; ix_0 < fcn_ptc.nx; ix_0++)
				{

					for (auto iy_0 = 0; iy_0 < fcn_ptc.ny; iy_0++)
					{
						const dt_int32 ix = ix_0 + fcn_ptc.ix_0;
						const dt_int32 iy = iy_0 + fcn_ptc.iy_0;

						const auto r2 = grid_2d.r2(ix, iy, fcn_ptc.r);
						if (r2 < fcn_ptc.r_max_2)
						{
							const auto ixy = grid_2d.sub_2_ind_pbc(ix, iy);
							mx_o[ixy] += fcn_ptc.eval_r2(r2);
						}
					}
				}
			}
		}
	}

	/***************************************************************************************/
	/********************************** superposition **************************************/
	/***************************************************************************************/
	namespace mt
	{	
		template <class T, eDim Dim, eFcn_typ Fcn_typ, eDev Dev> class Spt_Ptc_xd;

		template <class T, eDev Dev>
		using Spt_gauss_2d = Spt_Ptc_xd<T, edim_2, efcn_gauss, Dev>;
	}

	/* template specialization cpu */
	namespace mt
	{
		template <class T, eDim Dim, eFcn_typ Fcn_typ>
		class Spt_Ptc_xd<T, Dim, Fcn_typ, edev_cpu>
		{
			public:
				using T_r = T;
				using size_type = dt_uint64;
				static const eDev device = edev_cpu;

				using TVctr_r = Vctr<T, edev_cpu>;
				using Vctr_Ptc_dev = Vctr<Ptc_s_fcn_xd<T, Dim, Fcn_typ>, edev_cpu>;
				using Stream_dev = Stream<edev_cpu>;

				Spt_Ptc_xd(): in_spt_ptc(nullptr), stream(nullptr), n_ptc_blk(512) {}

				Spt_Ptc_xd(In_Spt_Ptc_xd<T, Dim, Fcn_typ> *in_spt_ptc, Stream_dev *stream): n_ptc_blk(512)
				{
					set_in_data(in_spt_ptc, stream);
				}

				void set_in_data(In_Spt_Ptc_xd<T, Dim, Fcn_typ> *in_spt_ptc, Stream_dev *stream)
				{	
					this->in_spt_ptc = in_spt_ptc;
					this->stream = stream;

					n_ptc_blk = stream->size();
					ptc_sp.resize(n_ptc_blk);
				}

				void operator()(TVctr_r &mx_io)
				{
					auto thr_gauss_eval = [](Stream_cpu& stream, Grid_2d<T>& grid, Vctr_Ptc_dev& ptc_sp, TVctr_r& mx_o)
					{
						if (stream.n_stream_act<= 0)
						{
							return;
						}

						for(auto istm = 0; istm < stream.n_stream_act-1; istm++)
						{
							stream[istm] = std::thread(cpu_detail::fcn_ptc_s_fcn_eval<T, Dim, Fcn_typ>, std::ref(stream), std::ref(grid), std::ref(ptc_sp[istm]), std::ref(mx_o));
						}

						cpu_detail::fcn_ptc_s_fcn_eval<T, Dim, Fcn_typ>(stream, grid, ptc_sp[stream.n_stream_act-1], mx_o);

						stream.synchronize();
					};

					mx_io.fill(0);

					stream->set_grid(in_spt_ptc->grid.shape());

					const dt_int32 iptc_0 = 0;
					const dt_int32 iptc_e = in_spt_ptc->ptc.size()-1;

					dt_int32 iptc = iptc_0;
					while (iptc <= iptc_e)
					{
						dt_int32 n_ptc = min(n_ptc_blk, iptc_e-iptc+1);
						stream->set_n_stream_act(n_ptc);
						set_ptc_s_fcn(iptc, n_ptc, ptc_sp);

						thr_gauss_eval(*stream, in_spt_ptc->grid, ptc_sp, mx_io);

						iptc += n_ptc;
					}
				}

				In_Spt_Ptc_xd<T, Dim, Fcn_typ> *in_spt_ptc;
				Stream_dev *stream;
			private:
				dt_int32 n_ptc_blk;

				void set_ptc_s_fcn(dt_int32 iptc, dt_int32 n_ptc_blk, Vctr_Ptc_dev& ptc_sp)
				{
					for(auto istm = 0; istm < n_ptc_blk; istm++)
					{
						ptc_sp[istm] = in_spt_ptc->ptc_s_fcn(iptc++, 0.85);
					}
				}

				Vctr_Ptc_dev ptc_sp;
		};
	}
	
	/* template specialization gpu */
	namespace mt
	{
		template <class T, eDim Dim, eFcn_typ Fcn_typ>
		class Spt_Ptc_xd<T, Dim, Fcn_typ, edev_gpu>
		{
			public:
				using T_r = T;
				using size_type = dt_uint64;
				static const eDev device = edev_cpu;

				using TVctr_r = Vctr<T, edev_gpu>;
				using Vctr_Ptc_cpu = Vctr<Ptc_s_fcn_xd<T, Dim, Fcn_typ>, edev_cpu>;
				using Vctr_Ptc_dev = Vctr<Ptc_s_fcn_xd<T, Dim, Fcn_typ>, edev_gpu>;
				using Stream_dev = Stream<edev_gpu>;

				Spt_Ptc_xd(): in_spt_ptc(nullptr), stream(nullptr), n_ptc_blk(512) {}

				Spt_Ptc_xd(In_Spt_Ptc_xd<T, Dim, Fcn_typ>*in_spt_ptc, Stream_dev *stream): n_ptc_blk(512)
				{
					set_in_data(in_spt_ptc, stream);
				}

				void set_in_data(In_Spt_Ptc_xd<T, Dim, Fcn_typ> *in_spt_ptc, Stream_dev *stream)
				{	
					this->in_spt_ptc = in_spt_ptc;
					this->stream = stream;

					n_ptc_blk = 512;
					ptc_sp_cpu.resize(n_ptc_blk);
					ptc_sp.resize(n_ptc_blk);
				}

				void operator()(TVctr_r &mx_io)
				{
					mx_io.fill(0);

					const dt_int32 iptc_0 = 0;
					const dt_int32 iptc_e = in_spt_ptc->ptc.size_32()-1;

					dt_int32 iptc = iptc_0;
					while (iptc <= iptc_e)
					{
						dt_int32 n_ptc = min(n_ptc_blk, iptc_e-iptc+1);
						set_ptc_s_fcn(iptc, n_ptc, ptc_sp);

						gpu_detail::fcn_ptc_s_fcn_eval<T, Dim, Fcn_typ><<<n_ptc, fcn_cdb_2d()>>>(in_spt_ptc->grid, ptc_sp.ptr_32(), mx_io.ptr_32());

						iptc += n_ptc;
					}
				}

				In_Spt_Ptc_xd<T, Dim, Fcn_typ> *in_spt_ptc;
				Stream_dev *stream;
			private:
				dt_int32 n_ptc_blk;

				void set_ptc_s_fcn(dt_int32 iptc, dt_int32 n_ptc_blk, Vctr_Ptc_dev& ptc_sp)
				{
					for(auto istm = 0; istm < n_ptc_blk; istm++)
					{
						ptc_sp_cpu[istm] = in_spt_ptc->ptc_s_fcn(iptc++, 0.85);
					}
					ptc_sp = ptc_sp_cpu;
				}

				Vctr_Ptc_cpu ptc_sp_cpu;
				Vctr_Ptc_dev ptc_sp;
		};
	}

#endif