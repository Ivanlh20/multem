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

//#define __CUDACC__

#pragma once

// https:// github.com/pfultz2/Cloak/wiki/C-Preprocessor-tricks, -tips, -and-idioms
#define COMMA ,

#define EVAL(...) __VA_ARGS__

#define PCAT2(a, ...) a ## __VA_ARGS__
#define CAT2(a, ...) PCAT2(a, __VA_ARGS__)
#define CAT2_2(a, ...) CAT2(a, __VA_ARGS__)

#define PCAT3(a, b, ...) a ## b ## __VA_ARGS__
#define CAT3(a, b, ...) PCAT3(a, b, __VA_ARGS__)
#define CAT3_2(a, b, ...) CAT3(a, b, __VA_ARGS__)

/***************************************************************************************/
#define IIF_1(t, ...) t
#define IIF_0(t, ...) __VA_ARGS__
#define IIF(c) PCAT2(IIF_, c)

/***************************************************************************************/
#define CHECK_N(x, n, ...) n
#define CHECK(...) CHECK_N(__VA_ARGS__, 0)
#define IS_1 x, 1
#define IS_ONE(x) CHECK(CAT2(IS_, x))

/***************************************************************************************/
#define REPV_1(N) N
#define REPV_2(N) REPV_1(N), N
#define REPV_3(N) REPV_2(N), N
#define REPV_4(N) REPV_3(N), N
#define REPV_5(N) REPV_4(N), N
#define REPV_6(N) REPV_5(N), N
#define REPV_7(N) REPV_6(N), N
#define REPV_8(N) REPV_7(N), N
#define REPV_9(N) REPV_8(N), N
#define REPV(V, N) CAT2(REPV_, N(V))

/***************************************************************************************/
#define IDX_1D ix
#define IDX_2D ix, iy
#define IDX_3D ix, iy, iz
#define IDX_ND(DIM) CAT3(IDX_, DIM, D)

/***************************************************************************************/
#define IDX_CDEF_1D const dt_int32& ix
#define IDX_CDEF_2D const dt_int32& ix, const dt_int32& iy
#define IDX_CDEF_3D const dt_int32& ix, const dt_int32& iy, const dt_int32& iz
#define IDX_CDEF_ND(DIM) CAT3(IDX_CDEF_, DIM, D)

/***************************************************************************************/
#define F_SHAPE_1D(G) G.nx
#define F_SHAPE_2D(G) G.nx, G.ny
#define F_SHAPE_3D(G) G.nx, G.ny, G.nz
#define F_SHAPE_ND(G, DIM) F_SHAPE_##DIM##D(G)

/***************************************************************************************/
#define EDIM(DIM) edim_##DIM

/***************************************************************************************/
// macro for int definition for blas and lapack
#ifdef MATLAB_BLAS_LAPACK
	#define MTL_BL_INT ptrdiff_t
#else
	#define MTL_BL_INT dt_int32
#endif

#define PASTER(pre, x, y) pre ## _ ## x ## _ ## y
#define EVAL_PASTE(pre, x, y) PASTER(pre, x, y)
#define EVAL_2_PASTE(pre, x, y) EVAL_PASTE(pre, x, y)

#define INC(X, K) EVAL_PASTE(INC, K, X)
#define DEC(X, K) EVAL_PASTE(DEC, K, X)

#define INC_1_0 1
#define INC_1_1 2
#define INC_1_2 3
#define INC_1_3 4
#define INC_1_4 5
#define INC_1_5 6
#define INC_1_6 7
#define INC_1_7 8
#define INC_1_8 9
#define INC_1_9 10
#define INC_1_10 11
#define INC_1_11 12
#define INC_1_12 13
#define INC_1_13 14
#define INC_1_14 15
#define INC_1_15 16

#define INC_2_0 2
#define INC_2_1 3
#define INC_2_2 4
#define INC_2_3 5
#define INC_2_4 6
#define INC_2_5 7
#define INC_2_6 8
#define INC_2_7 9
#define INC_2_8 10
#define INC_2_9 11
#define INC_2_10 12
#define INC_2_11 13
#define INC_2_12 14
#define INC_2_13 15
#define INC_2_14 16
#define INC_2_15 17

#define INC_3_0 3
#define INC_3_1 4
#define INC_3_2 5
#define INC_3_3 6
#define INC_3_4 7
#define INC_3_5 8
#define INC_3_6 9
#define INC_3_7 10
#define INC_3_8 11
#define INC_3_9 12
#define INC_3_10 13
#define INC_3_11 14
#define INC_3_12 15
#define INC_3_13 16
#define INC_3_14 17
#define INC_3_15 18

#define DEC_1_0 0
#define DEC_1_1 0
#define DEC_1_2 1
#define DEC_1_3 2
#define DEC_1_4 3
#define DEC_1_5 4
#define DEC_1_6 5
#define DEC_1_7 6
#define DEC_1_8 7
#define DEC_1_9 8
#define DEC_1_10 9
#define DEC_1_11 10
#define DEC_1_12 12
#define DEC_1_13 13
#define DEC_1_14 14
#define DEC_1_15 15

#define DEC_2_0 0
#define DEC_2_1 0
#define DEC_2_2 0
#define DEC_2_3 1
#define DEC_2_4 2
#define DEC_2_5 3
#define DEC_2_6 4
#define DEC_2_7 5
#define DEC_2_8 6
#define DEC_2_9 7
#define DEC_2_10 8
#define DEC_2_11 9
#define DEC_2_12 10
#define DEC_2_13 12
#define DEC_2_14 13
#define DEC_2_15 14

#define DEC_3_0 0
#define DEC_3_1 0
#define DEC_3_2 0
#define DEC_3_3 0
#define DEC_3_4 1
#define DEC_3_5 2
#define DEC_3_6 3
#define DEC_3_7 4
#define DEC_3_8 5
#define DEC_3_9 6
#define DEC_3_10 7
#define DEC_3_11 8
#define DEC_3_12 9
#define DEC_3_13 10
#define DEC_3_14 11
#define DEC_3_15 12

#ifdef __CUDACC__
	#include <cuda_runtime.h>

	#define CPU_EXEC __host__
	#define CPU_EXEC_INL __host__ inline

	#define GPU_EXEC __device__
	#define GPU_EXEC_INL __device__ inline

	#define CGPU_EXEC __host__ __device__
	#define CGPU_EXEC_INL __host__ __device__ inline

	#define FORCE_INLINE __forceinline__

	#define FOR_IX_1DC(nx) for(dt_int32 ix = threadIdx.x + blockIdx.x*blockDim.x; ix < nx; ix += blockDim.x*gridDim.x)

	#define FOR_IX_2DC(nx) for(dt_int32 ix = blockIdx.y*blockDim.y + threadIdx.y; ix < nx; ix += blockDim.y*gridDim.y)
	#define FOR_IY_2DC(ny) for(dt_int32 iy = blockIdx.x*blockDim.x + threadIdx.x; iy < ny; iy += blockDim.x*gridDim.x)

	#define FOR_LOOP_1DC(nx, fcn)																			\
	for(dt_int32 ix = threadIdx.x + blockIdx.x*blockDim.x; ix < nx; ix += blockDim.x*gridDim.x)				\
	{																										\
		fcn;																								\
	}

	#define FOR_LOOP_2DC(nx, ny, fcn)																		\
	for(dt_int32 ix = blockIdx.y*blockDim.y + threadIdx.y; ix < nx; ix += blockDim.y*gridDim.y)				\
	{																										\
		for(dt_int32 iy = blockIdx.x*blockDim.x + threadIdx.x; iy < ny; iy += blockDim.x*gridDim.x)			\
		{																									\
			fcn;																							\
		}																									\
	}

	//check out for(dt_int32 iz = blockIdx.z*blockDim.z + threadIdx.z; iz < nz; iz += blockDim.z*gridDim.z)
	#define FOR_LOOP_3DC(nx, ny, nz, fcn)																	\
	{																										\
		dt_int32 iz = threadIdx.z;																			\
		for(dt_int32 ix = blockIdx.y*blockDim.y + threadIdx.y; ix < nx; ix += blockDim.y*gridDim.y)			\
		{																									\
			for(dt_int32 iy = blockIdx.x*blockDim.x + threadIdx.x; iy < ny; iy += blockDim.x*gridDim.x)		\
			{																								\
				fcn;																						\
			}																								\
		}																									\
	}
		
	#define FOR_LOOP_NDC(DIM, ...) CAT3(FOR_LOOP_, DIM, DC)(__VA_ARGS__)

#else
	#define CPU_EXEC
	#define CPU_EXEC_INL inline

	#define GPU_EXEC
	#define GPU_EXEC_INL inline

	#define CGPU_EXEC
	#define CGPU_EXEC_INL inline

	#define FORCE_INLINE inline

	#ifdef _MSC_VER
		#define __restrict__ __restrict
	#endif
#endif