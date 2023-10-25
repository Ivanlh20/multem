/*
* This file is part of Multem.
* Copyright 2022 Ivan Lobato <Ivanlh20@gmail.com>
*
* Multem is destroy software: you can redistribute it and/or modify
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

#pragma once

#include <vector>

#include "macros.h"
#include "const_enum.h"
#include "type_traits_gen.h"
#include "igrid_1d.h"
#include "igrid_2d.h"
#include "igrid_3d.h"

#ifdef __CUDACC__
	#include <cuda.h>
	#include <cuda_runtime.h>
#endif

/* macro definition grid-block */
namespace mt
{
#define	FCNS_DEF_GPU_GRID_BLK_VCTR																\
 	dim3 d_blk_size();																			\
																								\
 	dim3 d_grid_size(const dim3 d_grid_max = dim3(128, 1, 1));									\
																								\
	/***************************************************************************************/	\
 	dim3 d_blk_1d();																			\
																								\
 	dim3 d_grid_1d(const dim3 d_grid_max = dim3(128, 1, 1));									\
																								\
 	D_Grid_Blk d_grid_blk_1d(const dim3 d_grid_max = dim3(128, 1, 1));							\
																								\
	/***************************************************************************************/	\
 	dim3 d_grid_1d_h(const dim3 d_grid_max = dim3(128, 1, 1));									\
																								\
	D_Grid_Blk d_grid_blk_1d_h(const dim3 d_grid_max = dim3(128, 1, 1));						\
																								\
	/***************************************************************************************/	\
 	dim3 d_blk_2d();																			\
																								\
	/***************************************************************************************/	\
 	dim3 d_grid_2d(const dim3 d_grid_max = dim3(64, 64, 0));									\
																								\
 	D_Grid_Blk d_grid_blk_2d(const dim3 d_grid_max = dim3(64, 64, 1));							\
																								\
	/***************************************************************************************/	\
 	dim3 d_grid_2d_h(const dim3 d_grid_max = dim3(64, 64, 1));									\
																								\
	D_Grid_Blk d_grid_blk_2d_h(const dim3 d_grid_max = dim3(64, 64, 1));						\
																								\
	/***************************************************************************************/	\
 	dim3 d_blk_3d();																			\
																								\
	/***************************************************************************************/	\
 	dim3 d_grid_3d(const dim3 d_grid_max = dim3(64, 64, 64));									\
																								\
 	D_Grid_Blk d_grid_blk_3d(const dim3 d_grid_max = dim3(64, 64, 64));							\
																								\
	/***************************************************************************************/	\
 	dim3 d_grid_3d_h(const dim3 d_grid_max = dim3(64, 64, 64));									\
																								\
	D_Grid_Blk d_grid_blk_h(const dim3 d_grid_max = dim3(64, 64, 64))
}

/* vector forward declaration */
namespace mt
{
#ifndef VCTR_DEC
	#define VCTR_DEC
	template <class T, eDev Dev> class Vctr;

	template <class T, eDev Dev, class ST> class pVctr;
#endif
}

/* derived class */
namespace mt
{
	template <class T, eDev Dev>
	using pVctr_32 = pVctr<T, Dev, dt_int32>;

	template <class T, eDev Dev>
	using pVctr_64 = pVctr<T, Dev, dt_int64>;

	template <class T>
	using pVctr_cpu_32 = pVctr<T, edev_cpu, dt_int32>;

	template <class T>
	using pVctr_cpu_64 = pVctr<T, edev_cpu, dt_int64>;

	template <class T>
	using pVctr_gpu_32 = pVctr<T, edev_gpu, dt_int32>;

	template <class T>
	using pVctr_gpu_64 = pVctr<T, edev_gpu, dt_int64>;
}

/* vector pointer */
namespace mt
{
	template <class T, eDev Dev, class ST>
	class pVctr
	{
	public:
		using value_type = T;
		using size_type = ST;
		static const eDev device = Dev;

		T* m_data;
		ST m_s0;
		ST m_s1;
		ST m_s2;
		ST m_s3;
		ST m_size;

		ST m_pitch_s1;
		ST m_pitch_s2;
		ST m_pitch_s3;

		/************************************* constructors ************************************/
		CGPU_EXEC
		pVctr();

		CPU_EXEC
		pVctr(T* data, dt_shape_64 shape);

		CGPU_EXEC
		pVctr(T* data, ST s0);

		CGPU_EXEC
		pVctr(T* data, ST s0, ST s1);

		CGPU_EXEC
		pVctr(T* data, ST s0, ST s1, ST s2);

		CGPU_EXEC
		pVctr(T* data, ST s0, ST s1, ST s2, ST s3);

		/************************** constructors *****************************/
		/* copy constructor */
		CGPU_EXEC
		pVctr(const pVctr<T, Dev, ST>& pvctr);

		/* Move constructor */
		CGPU_EXEC
		pVctr(pVctr<T, Dev, ST>&& pvctr);

		/* Converting constructor */
		template <class STU>
		CPU_EXEC
		explicit pVctr(const pVctr<T, Dev, STU>& pvctr);

		/* constructor from Vctr */
		CPU_EXEC
		explicit pVctr(const Vctr<T, Dev>& vctr);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		pVctr<T, Dev, ST>& operator=(const pVctr<T, Dev, ST>& pvctr);

		/* Move assignment operator */
		CGPU_EXEC
		pVctr<T, Dev, ST>& operator=(pVctr<T, Dev, ST>&& pvctr);

		/* Converting assignment operator */
		template <class STU>
		CPU_EXEC
		pVctr<T, Dev, ST>& operator=(const pVctr<T, Dev, STU>& pvctr);

		/* Assignment operator: Vctr -> pVctr */
		CPU_EXEC
		pVctr<T, Dev, ST>& operator=(const Vctr<T, Dev>& vctr);

		/**************** user define conversion operators *******************/
		pVctr<T, Dev, dt_int32> ptr_32() const;

		pVctr<T, Dev, dt_int64> ptr_64() const;

		operator pVctr<T, Dev, chg_btw_int32_int64<ST>>() const;

		CGPU_EXEC
		ST s0() const;

		CGPU_EXEC
		ST s1() const;

		CGPU_EXEC
		ST s2() const;

		CGPU_EXEC
		ST s3() const;

		CGPU_EXEC
		dt_int32 s0_32() const;

		CGPU_EXEC
		dt_int32 s1_32() const;

		CGPU_EXEC
		dt_int32 s2_32() const;

		CGPU_EXEC
		dt_int32 s3_32() const;
			
		CGPU_EXEC
		dt_int64 s0_64() const;

		CGPU_EXEC
		dt_int64 s1_64() const;

		CGPU_EXEC
		dt_int64 s2_64() const;

		CGPU_EXEC
		dt_int64 s3_64() const;

		CGPU_EXEC
		ST s0h() const;

		CGPU_EXEC
		ST s1h() const;

		CGPU_EXEC
		ST s2h() const;

		CGPU_EXEC
		ST s3h() const;

		dt_shape_st<ST> shape() const;

		dt_shape_st<ST> shape_2d_trs() const;

		CGPU_EXEC
		ST shape_size() const;

		CGPU_EXEC
		ST pitch_s1() const;

		CGPU_EXEC
		ST pitch_s2() const;

		CGPU_EXEC
		ST pitch_s3() const;

		CGPU_EXEC
		ST size() const;

		CGPU_EXEC
		dt_int32 size_32() const;

		CGPU_EXEC
		dt_int64 size_64() const;

		iGrid_1d igrid_1d() const;

		iGrid_2d igrid_2d() const;

		iGrid_3d igrid_3d() const;

		iGrid_1d_64 igrid_1d_64() const;

		iGrid_2d_64 igrid_2d_64() const;

		iGrid_3d_64 igrid_3d_64() const;

		CGPU_EXEC
		dt_bool empty() const;

		CGPU_EXEC
		dt_bool is_1d() const;

		CGPU_EXEC
		ST sub_2_ind(const ST& ix_0) const;

		CGPU_EXEC
		ST sub_2_ind(const ST& ix_0, const ST& ix_1) const;

		CGPU_EXEC
		ST sub_2_ind(const ST& ix_0, const ST& ix_1, const ST& ix_2) const;

		CGPU_EXEC
		ST sub_2_ind(const ST& ix_0, const ST& ix_1, const ST& ix_2, const ST& ix_3) const;

		CGPU_EXEC
		T& operator[](const ST& iy);

		CGPU_EXEC
		const T& operator[](const ST& iy) const;

		CGPU_EXEC
		T& operator()(const ST& iy);

		CGPU_EXEC
		const T& operator()(const ST& iy) const;

		CGPU_EXEC
		T& operator()(const ST& ix_0, const ST& ix_1);

		CGPU_EXEC
		const T& operator()(const ST& ix_0, const ST& ix_1) const ;

		CGPU_EXEC
		T& operator()(const ST& ix_0, const ST& ix_1, const ST& ix_2);

		CGPU_EXEC
		const T& operator()(const ST& ix_0, const ST& ix_1, const ST& ix_2) const;

		CGPU_EXEC
		T& operator()(const ST& ix_0, const ST& ix_1, const ST& ix_2, const ST& ix_3);

		CGPU_EXEC
		const T& operator()(const ST& ix_0, const ST& ix_1, const ST& ix_2, const ST& ix_3) const;

		CGPU_EXEC
		T* begin();

		CGPU_EXEC
		const T* begin() const;

		CGPU_EXEC
		T* end();

		CGPU_EXEC
		const T* end() const;

		CGPU_EXEC
		T* data();

		CGPU_EXEC
		const T* data() const;

		template <class U>
		U data_cast();

		template <class U>
		const U data_cast() const;

		CGPU_EXEC
		T front() const;

		CGPU_EXEC
		T back() const;

	#ifdef __CUDACC__
 		FCNS_DEF_GPU_GRID_BLK_VCTR;
	#endif

	private:
		CGPU_EXEC
		void set_picth();
	};
}

#include "../src/pvctr.inl"