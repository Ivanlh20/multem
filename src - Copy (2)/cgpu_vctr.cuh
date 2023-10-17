/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef CGPU_VCTR_H
	#define CGPU_VCTR_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include <vector>

	#include "macros.cuh"
	#include "const_enum.cuh"
	#include "type_traits_gen.cuh"
	#include "memcpy.cuh"
	#include "r_2d.cuh"
	#include "r_3d.cuh"
	#include "mx_2x2.cuh"
	#include "mx_3x3.cuh"
	#include "igrid.cuh"

	#ifdef __CUDACC__
		#include <cuda.h>
		#include <cuda_runtime.h>
		#include <thrust/detail/raw_pointer_cast.h>
		#include <thrust/device_vector.h>
		#include <thrust/host_vector.h>
	#endif

	/* vector copy */
	namespace mt
	{
		template <class T, eDev Dev> class Vctr;

		// (dst, src): cpu -> cpu
		template <class Td, class Ts>
		void vctr_cpy(Vctr<Td, edev_cpu>& vctr_cpu_dst, const Vctr<Ts, edev_cpu>& vctr_cpu_src, Ts* pvctr_cpu_jk = nullptr)
		{
			memcpy_cpu_cpu<Td, Ts>(vctr_cpu_dst.data(), vctr_cpu_src.data(), vctr_cpu_src.size(), pvctr_cpu_jk);
		}

	#ifdef __CUDACC__
		// (dst, src): gpu -> cpu
		template <class Td, class Ts>
		void vctr_cpy(Vctr<Td, edev_cpu>& vctr_cpu_dst, const Vctr<Ts, edev_gpu>& vctr_gpu_src, Ts* pvctr_cpu_jk = nullptr)
		{
			memcpy_gpu_cpu<Td, Ts>(vctr_cpu_dst.data(), vctr_gpu_src.data(), vctr_gpu_src.size(), pvctr_cpu_jk);
		}

		// (dst, src): cpu -> gpu
		template <class Td, class Ts>
		void vctr_cpy(Vctr<Td, edev_gpu>& vctr_gpu_dst, const Vctr<Ts, edev_cpu>& vctr_cpu_src, Td* pvctr_cpu_jk = nullptr)
		{
			memcpy_cpu_gpu<Td, Ts>(vctr_gpu_dst.data(), vctr_cpu_src.data(), vctr_cpu_src.size(), pvctr_cpu_jk);
		}

		// (dst, src): gpu -> gpu
		template <class Td, class Ts>
		void vctr_cpy(Vctr<Td, edev_gpu>& vctr_gpu_dst, const Vctr<Ts, edev_gpu>& vctr_gpu_src, Td* pvctr_cpu_jk = nullptr)
		{
			memcpy_gpu_gpu<Td, Ts>(vctr_gpu_dst.data(), vctr_gpu_src.data(), vctr_gpu_src.size(), pvctr_cpu_jk);
		}
	#endif
	}

	/* macro grid-block */
	namespace mt
	{
	#define	FCNS_GPU_GRID_BLK_VCTR																	\
 		dim3 d_blk_size()																			\
		{																							\
			return fcn_cdb_size();																	\
		}																							\
																									\
 		dim3 d_grid_size(const dim3 d_grid_max = dim3(128, 1, 1))									\
		{																							\
			auto grid = fcn_cdg_size(m_size);														\
																									\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;							\
																									\
			return grid;																			\
		}																							\
																									\
		/***************************************************************************************/	\
 		dim3 d_blk_1d()																				\
		{																							\
			return fcn_cdb_1d();																	\
		}																							\
																									\
 		dim3 d_grid_1d(const dim3 d_grid_max = dim3(128, 1, 1))										\
		{																							\
			auto grid = fcn_cdg_1d(m_s0);															\
																									\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;							\
																									\
			return grid;																			\
		}																							\
																									\
 		D_Grid_Blk d_grid_blk_1d(const dim3 d_grid_max = dim3(128, 1, 1))							\
		{																							\
			return D_Grid_Blk(d_grid_1d(d_grid_max), d_blk_1d());									\
		}																							\
																									\
		/***************************************************************************************/	\
 		dim3 d_grid_1d_h(const dim3 d_grid_max = dim3(128, 1, 1))									\
		{																							\
			auto grid = fcn_cdg_1d(m_s0/2);															\
																									\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;							\
																									\
			return grid;																			\
		}																							\
																									\
		D_Grid_Blk d_grid_blk_1d_h(const dim3 d_grid_max = dim3(128, 1, 1))							\
		{																							\
			return D_Grid_Blk(d_grid_1d_h(d_grid_max), d_blk_1d());									\
		}																							\
																									\
		/***************************************************************************************/	\
 		dim3 d_blk_2d()																				\
		{																							\
			return fcn_cdb_2d();																	\
		}																							\
																									\
		/***************************************************************************************/	\
 		dim3 d_grid_2d(const dim3 d_grid_max = dim3(64, 64, 0))										\
		{																							\
			auto grid = fcn_cdg_2d(m_s0, m_s1);														\
																									\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;							\
			grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;							\
																									\
			return grid;																			\
		}																							\
																									\
 		D_Grid_Blk d_grid_blk_2d(const dim3 d_grid_max = dim3(64, 64, 1))							\
		{																							\
			return D_Grid_Blk(d_grid_2d(d_grid_max), d_blk_2d());									\
		}																							\
																									\
		/***************************************************************************************/	\
 		dim3 d_grid_2d_h(const dim3 d_grid_max = dim3(64, 64, 1))									\
		{																							\
			auto grid = fcn_cdg_2d(m_s0/2, m_s1/2);													\
																									\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;							\
			grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;							\
																									\
			return grid;																			\
		}																							\
																									\
		D_Grid_Blk d_grid_blk_2d_h(const dim3 d_grid_max = dim3(64, 64, 1))							\
		{																							\
			return D_Grid_Blk(d_grid_2d_h(d_grid_max), d_blk_2d());									\
		}																							\
																									\
		/***************************************************************************************/	\
 		dim3 d_blk_3d()																				\
		{																							\
			return fcn_cdb_3d();																	\
		}																							\
																									\
		/***************************************************************************************/	\
 		dim3 d_grid_3d(const dim3 d_grid_max = dim3(64, 64, 64))									\
		{																							\
			auto grid = fcn_cdg_3d(m_s0, m_s1, m_s2);												\
																									\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;							\
			grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;							\
			grid.z = (d_grid_max.z > 0)?min(d_grid_max.z, grid.z):grid.z;							\
																									\
			return grid;																			\
		}																							\
																									\
 		D_Grid_Blk d_grid_blk_3d(const dim3 d_grid_max = dim3(64, 64, 64))							\
		{																							\
			return D_Grid_Blk(d_grid_3d(d_grid_max), d_blk_3d());									\
		}																							\
																									\
		/***************************************************************************************/	\
 		dim3 d_grid_3d_h(const dim3 d_grid_max = dim3(64, 64, 64))									\
		{																							\
			auto grid = fcn_cdg_3d(m_s0/2, m_s1/2, m_s2/2);											\
																									\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;							\
			grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;							\
			grid.z = (d_grid_max.z > 0)?min(d_grid_max.z, grid.z):grid.z;							\
																									\
			return grid;																			\
		}																							\
																									\
		D_Grid_Blk d_grid_blk_h(const dim3 d_grid_max = dim3(64, 64, 64))							\
		{																							\
			return D_Grid_Blk(d_grid_3d_h(d_grid_max), d_blk_3d());									\
		}
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
			pVctr(): m_data(nullptr), m_s0(0), m_s1(1), m_s2(1), m_s3(1), 
			m_size(0)
			{
				set_picth();
			}

			CPU_EXEC
			pVctr(T* data, dt_shape_64 shape): m_data(data), 
			m_s0(ST(shape[0])), m_s1(ST(shape[1])), 
			m_s2(ST(shape[2])), m_s3(ST(shape[3])), m_size(shape_size())
			{
				set_picth();
			}

			CGPU_EXEC
			pVctr(T* data, ST s0): m_data(data), 
			m_s0(s0), m_s1(1), m_s2(1), m_s3(1), m_size(shape_size())
			{
				set_picth();
			}

			CGPU_EXEC
			pVctr(T* data, ST s0, ST s1): m_data(data), 
			m_s0(s0), m_s1(s1), m_s2(1), m_s3(1), m_size(shape_size())
			{
				set_picth();
			}

			CGPU_EXEC
			pVctr(T* data, ST s0, ST s1, ST s2): m_data(data), 
			m_s0(s0), m_s1(s1), m_s2(s2), m_s3(1), m_size(shape_size())
			{
				set_picth();
			}

			CGPU_EXEC
			pVctr(T* data, ST s0, ST s1, ST s2, ST s3): m_data(data), 
			m_s0(s0), m_s1(s1), m_s2(s2), m_s3(s3), m_size(shape_size())
			{
				set_picth();
			}

			/************************** constructors *****************************/
			/* copy constructor */
			CGPU_EXEC
			pVctr(const pVctr<T, Dev, ST>& pvctr)
			{
				*this = pvctr;
			}

			/* Move constructor */
			CGPU_EXEC
			pVctr(pVctr<T, Dev, ST>&& pvctr)
			{
				*this = std::move(pvctr);
			}

			// ! Converting constructor
			template <class STU>
			CPU_EXEC
			explicit pVctr(const pVctr<T, Dev, STU>& pvctr)
			{
				*this = pvctr;
			}

			// ! constructor from Vctr
			CPU_EXEC
			explicit pVctr(const Vctr<T, Dev>& vctr)
			{
				*this = vctr;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			pVctr<T, Dev, ST>& operator=(const pVctr<T, Dev, ST>& pvctr)
			{
				if (this != &pvctr)
				{
					m_data = pvctr.m_data;
					m_s0 = pvctr.m_s0;
					m_s1 = pvctr.m_s1;
					m_s2 = pvctr.m_s2;
					m_s3 = pvctr.m_s3;
					m_size = pvctr.m_size;
					m_pitch_s1 = pvctr.m_pitch_s1;
					m_pitch_s2 = pvctr.m_pitch_s2;
					m_pitch_s3 = pvctr.m_pitch_s3;
				}

				return *this;
			}

			 /* Move assignment operator */
			CGPU_EXEC
			pVctr<T, Dev, ST>& operator=(pVctr<T, Dev, ST>&& pvctr)
			{
				if (this != &pvctr)
				{
					m_data = pvctr.m_data;
					m_s0 = pvctr.m_s0;
					m_s1 = pvctr.m_s1;
					m_s2 = pvctr.m_s2;
					m_s3 = pvctr.m_s3;
					m_size = pvctr.m_size;
					m_pitch_s1 = pvctr.m_pitch_s1;
					m_pitch_s2 = pvctr.m_pitch_s2;
					m_pitch_s3 = pvctr.m_pitch_s3;

					pvctr.m_data = nullptr;
					pvctr.m_s0 = 0;
					pvctr.m_s1 = 0;
					pvctr.m_s2 = 0;
					pvctr.m_s3 = 0;
					pvctr.m_size = 0;
					pvctr.m_pitch_s1 = 0;
					pvctr.m_pitch_s2 = 0;
					pvctr.m_pitch_s3 = 0;
				}

				return *this;
			}

			// ! Converting assignment operator
			template <class STU>
			CPU_EXEC
			pVctr<T, Dev, ST>& operator=(const pVctr<T, Dev, STU>& pvctr)
			{
				m_data = pvctr.m_data;
				m_s0 = ST(pvctr.m_s0);
				m_s1 = ST(pvctr.m_s1);
				m_s2 = ST(pvctr.m_s2);
				m_s3 = ST(pvctr.m_s3);
				m_size = ST(pvctr.m_size);
				m_pitch_s1 = ST(pvctr.m_pitch_s1);
				m_pitch_s2 = ST(pvctr.m_pitch_s2);
				m_pitch_s3 = ST(pvctr.m_pitch_s3);

				return *this;
			}

			 // ! Assignment operator: Vctr -> pVctr
			CPU_EXEC
			pVctr<T, Dev, ST>& operator=(const Vctr<T, Dev>& vctr)
			{
				m_data = vctr.m_data;
				m_s0 = ST(vctr.m_s0);
				m_s1 = ST(vctr.m_s1);
				m_s2 = ST(vctr.m_s2);
				m_s3 = ST(vctr.m_s3);
				m_size = ST(vctr.m_size);
				m_pitch_s1 = ST(vctr.m_pitch_s1);
				m_pitch_s2 = ST(vctr.m_pitch_s2);
				m_pitch_s3 = ST(vctr.m_pitch_s3);

				return *this;
			}

			/**************** user define conversion operators *******************/
			pVctr<T, Dev, dt_int32> ptr_32() const
			{
				return pVctr<T, Dev, dt_int32>(*this);
			}

			pVctr<T, Dev, dt_int64> ptr_64() const
			{
				return pVctr<T, Dev, dt_int64>(*this);
			}

			operator pVctr<T, Dev, chg_btw_int32_int64<ST>>() const
			{
				return pVctr<T, Dev, chg_btw_int32_int64<ST>>(*this);
			}

			CGPU_EXEC
			ST s0() const
			{
				return m_s0;
			}

			CGPU_EXEC
			ST s1() const
			{
				return m_s1;
			}

			CGPU_EXEC
			ST s2() const
			{
				return m_s2;
			}

			CGPU_EXEC
			ST s3() const
			{
				return m_s3;
			}

			CGPU_EXEC
			dt_int32 s0_32() const
			{
				return static_cast<dt_int32>(m_s0);
			}

			CGPU_EXEC
			dt_int32 s1_32() const
			{
				return static_cast<dt_int32>(m_s1);
			}

			CGPU_EXEC
			dt_int32 s2_32() const
			{
				return static_cast<dt_int32>(m_s2);
			}

			CGPU_EXEC
			dt_int32 s3_32() const
			{
				return static_cast<dt_int32>(m_s3);
			}
			
			CGPU_EXEC
			dt_int64 s0_64() const
			{
				return static_cast<dt_int64>(m_s0);
			}

			CGPU_EXEC
			dt_int64 s1_64() const
			{
				return static_cast<dt_int64>(m_s1);
			}

			CGPU_EXEC
			dt_int64 s2_64() const
			{
				return static_cast<dt_int64>(m_s2);
			}

			CGPU_EXEC
			dt_int64 s3_64() const
			{
				return static_cast<dt_int64>(m_s3);
			}

			CGPU_EXEC
			ST s0h() const
			{
				return m_s0/ST(2);
			}

			CGPU_EXEC
			ST s1h() const
			{
				return m_s1/ST(2);
			}

			CGPU_EXEC
			ST s2h() const
			{
				return m_s2/ST(2);
			}

			CGPU_EXEC
			ST s3h() const
			{
				return m_s3/ST(2);
			}

			dt_shape_st<ST> shape() const
			{
				return {m_s0, m_s1, m_s2, m_s3};
			}

			dt_shape_st<ST> shape_2d_trs() const
			{
				return {m_s1, m_s0, m_s2, m_s3};
			}

			CGPU_EXEC
			ST shape_size() const 
			{
				return m_s0*max(m_s1, ST(1))*max(m_s2, ST(1))*max(m_s3, ST(1));
			}

			CGPU_EXEC
			ST pitch_s1() const
			{
				return m_pitch_s1;
			}

			CGPU_EXEC
			ST pitch_s2() const
			{
				return m_pitch_s2;
			}

			CGPU_EXEC
			ST pitch_s3() const
			{
				return m_pitch_s3;
			}

			CGPU_EXEC
			ST size() const
			{
				return m_size;
			}

			CGPU_EXEC
			dt_int32 size_32() const
			{
				return static_cast<dt_int32>(m_size);
			}

			CGPU_EXEC
			dt_int64 size_64() const
			{
				return static_cast<dt_int64>(m_size);
			}

			iGrid_1d igrid_1d() const
			{
				return {size_32()};
			}

			iGrid_2d igrid_2d() const
			{
				return {s1_32(), s0_32()};
			}

			iGrid_2d igrid_3d() const
			{
				return {s1_32(), s0_32(), s2_32()};
			}

			iGrid_1d_64 igrid_1d_64() const
			{
				return {size_64()};
			}

			iGrid_2d_64 igrid_2d_64() const
			{
				return {s1_64(), s0_64()};
			}

			iGrid_2d_64 igrid_3d_64() const
			{
				return {s1_64(), s0_64(), s2_64()};
			}

			CGPU_EXEC
			dt_bool empty() const
			{
				return m_size == 0;
			}

			CGPU_EXEC
			dt_bool is_1d() const
			{
				return (m_s0 == 1) || (m_s1 == 1);
			}

			CGPU_EXEC
			ST sub_2_ind(const ST& ix_0) const 
			{ 
				return ix_0;
			}

			CGPU_EXEC
			ST sub_2_ind(const ST& ix_0, const ST& ix_1) const 
			{ 
				return ix_0 + m_s0*ix_1;
			}

			CGPU_EXEC
			ST sub_2_ind(const ST& ix_0, const ST& ix_1, const ST& ix_2) const 
			{ 
				return ix_0 + m_s0*(ix_1 + m_s1*ix_2);
			}

			CGPU_EXEC
			ST sub_2_ind(const ST& ix_0, const ST& ix_1, const ST& ix_2, const ST& ix_3) const 
			{ 
				return ix_0 + m_s0*(ix_1 + m_s1*(ix_2 + m_s2*ix_3));
			}

			CGPU_EXEC
			T& operator[](const ST& iy)
			{ 
				return m_data[iy];
			}

			CGPU_EXEC
			const T& operator[](const ST& iy) const 
			{ 
				return m_data[iy];
			}

			CGPU_EXEC
			T& operator()(const ST& iy)
			{ 
				return m_data[iy];
			}

			CGPU_EXEC
			const T& operator()(const ST& iy) const 
			{ 
				return m_data[iy];
			}

			CGPU_EXEC
			T& operator()(const ST& ix_0, const ST& ix_1)
			{ 
				return m_data[sub_2_ind(ix_0, ix_1)];
			}

			CGPU_EXEC
			const T& operator()(const ST& ix_0, const ST& ix_1) const 
			{ 
				return m_data[sub_2_ind(ix_0, ix_1)];
			}

			CGPU_EXEC
			T& operator()(const ST& ix_0, const ST& ix_1, const ST& ix_2)
			{ 
				return m_data[sub_2_ind(ix_0, ix_1, ix_2)];
			}

			CGPU_EXEC
			const T& operator()(const ST& ix_0, const ST& ix_1, const ST& ix_2) const 
			{ 
				return m_data[sub_2_ind(ix_0, ix_1, ix_2)];
			}

			CGPU_EXEC
			T& operator()(const ST& ix_0, const ST& ix_1, const ST& ix_2, const ST& ix_3)
			{ 
				return m_data[sub_2_ind(ix_0, ix_1, ix_2, ix_3)];
			}

			CGPU_EXEC
			const T& operator()(const ST& ix_0, const ST& ix_1, const ST& ix_2, const ST& ix_3) const 
			{ 
				return m_data[sub_2_ind(ix_0, ix_1, ix_2, ix_3)];
			}

			CGPU_EXEC
			T* begin()
			{ 
				return m_data;
			}

			CGPU_EXEC
			const T* begin() const 
			{ 
				return m_data;
			}

			CGPU_EXEC
			T* end()
			{ 
				return m_data + m_size;
			}

			CGPU_EXEC
			const T* end() const 
			{ 
				return m_data + m_size;
			}

			CGPU_EXEC
			T* data()
			{
				return m_data;
			}

			CGPU_EXEC
			const T* data() const
			{
				return m_data;
			}

			template <class U>
			U data_cast()
			{
				return reinterpret_cast<U>(m_data);
			}

			template <class U>
			const U data_cast() const
			{
				return reinterpret_cast<U>(m_data);
			}

			CGPU_EXEC
			T front() const 
			{ 
				return m_data[0];
			}

			CGPU_EXEC
			T back() const 
			{ 
				return m_data[m_size-1];
			}

	#ifdef __CUDACC__
			FCNS_GPU_GRID_BLK_VCTR;
	#endif

		private:

			CGPU_EXEC
			void set_picth()
			{
				m_pitch_s1 = m_s0*sizeof(T);
				m_pitch_s2 = m_pitch_s1*m_s1;
				m_pitch_s3 = m_pitch_s2*m_s2;
			}
		};

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

	/* cpu vector */
	namespace mt
	{
		template <class T>
		class Vctr<T, edev_cpu>
		{
		public:
			using value_type = T;
			using size_type = dt_int64;
			static const eDev device = edev_cpu;

			mutable T* m_data;
			size_type m_s0;
			size_type m_s1;
			size_type m_s2;
			size_type m_s3;
			size_type m_size;
			size_type m_capacity;

			size_type m_pitch_s1;
			size_type m_pitch_s2;
			size_type m_pitch_s3;

			/************************************* constructors ************************************/
			explicit Vctr(): m_data(nullptr), m_s0(0), m_s1(1), m_s2(1), m_s3(1), 
			m_size(shape_size()), m_capacity(0)
			{
				set_picth();
			}

			Vctr(const dt_init_list_f64& data): Vctr()
			{
				assign(data.begin(), data.end());
			}

			Vctr(size_type s0): Vctr()
			{
				resize({ s0, size_type(1), size_type(1), size_type(1) });
			}

			Vctr(size_type s0, const T& value): Vctr()
			{
				resize({ s0, size_type(1), size_type(1), size_type(1) }, value);
			}

			explicit Vctr(const dt_shape_st<size_type>& shape): Vctr()
			{
				resize(shape);
			}

			explicit Vctr(const dt_shape_st<size_type>& shape, const T& value): Vctr()
			{
				resize(shape, value);
			}

			/* copy constructor */
			Vctr(const Vctr<T, edev_cpu>& vctr): Vctr()
			{
				*this = vctr;
			}

			/* Move constructor */
			Vctr(Vctr<T, edev_cpu>&& vctr): Vctr()
			{
				*this = std::move(vctr);
			}

			/* converting constructor */
			template <class U>
			Vctr(const Vctr<U, edev_cpu>& vctr): Vctr()
			{
				assign(vctr);
			}

			template <class U>
			Vctr(U* first, U* last): Vctr()
			{
				assign(first, last);
			}

			template <class U, class V=T, typename = enable_if_r_nd<V>>
			Vctr(U *p, dt_int64 n_p, dt_int64 icol=0): Vctr()
			{
				set_pos_from_cpu_ptr(p, n_p, icol);
			}

			template <class U>
			Vctr(const std::vector<U>& vctr): Vctr()
			{
				assign(vctr);
			}

			// from cpu pVctr to Vctr
			template <class U, class STU>
			Vctr(const pVctr<U, edev_cpu, STU>& pvctr): Vctr()
			{
				assign(pvctr);
			}

	#ifdef __CUDACC__
			template <class U>
			Vctr(const Vctr<U, edev_gpu>& vctr): Vctr()
			{
				assign(vctr);
			}

			template <class U>
			Vctr(const thrust::host_vector<U>& vctr): Vctr()
			{
				assign(vctr);
			}

			template <class U>
			Vctr(const thrust::device_vector<U>& vctr): Vctr()
			{
				assign(vctr);
			}

			// from gpu pVctr to Vctr
			template <class U, class STU>
			Vctr(const pVctr<U, edev_gpu, STU>& pvctr): Vctr()
			{
				assign(pvctr);
			}
	#endif
			~Vctr()
			{
				destroy();
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			Vctr<T, edev_cpu>& operator=(const Vctr<T, edev_cpu>& vctr)
			{
				if (vctr.m_data == nullptr)
				{
					destroy();
				}
				else if (this != &vctr)
				{
					delete[] m_data;

					m_data = new T[vctr.m_capacity];
					memcpy_cpu_cpu(m_data, vctr.data(), vctr.size());

					m_s0 = vctr.m_s0;
					m_s1 = vctr.m_s1;
					m_s3 = vctr.m_s3;
					m_size = vctr.m_size;
					m_capacity = vctr.m_capacity;	
					m_pitch_s1 = vctr.m_pitch_s1;
					m_pitch_s2 = vctr.m_pitch_s2;
					m_pitch_s3 = vctr.m_pitch_s3;
				}

				return *this;
			}

			 /* Move assignment operator */
			Vctr<T, edev_cpu>& operator=(Vctr<T, edev_cpu>&& vctr)
			{
				if (vctr.m_data == nullptr)
				{
					destroy();
				}
				else if (this != &vctr)
				{
					delete[] m_data;

					m_data = vctr.m_data;
					m_s0 = vctr.m_s0;
					m_s1 = vctr.m_s1;
					m_s2 = vctr.m_s2;
					m_s3 = vctr.m_s3;
					m_size = vctr.m_size;
					m_capacity = vctr.m_capacity;
					m_pitch_s1 = vctr.m_pitch_s1;
					m_pitch_s2 = vctr.m_pitch_s2;
					m_pitch_s3 = vctr.m_pitch_s3;

					vctr.m_data = nullptr;
					vctr.m_s0 = 0;
					vctr.m_s1 = 0;
					vctr.m_s2 = 0;
					vctr.m_s3 = 0;
					vctr.m_size = 0;
					vctr.m_capacity = 0;
					vctr.m_pitch_s1 = 0;
					vctr.m_pitch_s2 = 0;
					vctr.m_pitch_s3 = 0;
				}

				return *this;
			}

			/* converting assignment operator */
			template <class U>
			Vctr<T, edev_cpu>& operator=(const Vctr<U, edev_cpu>& vctr) 
			{
				assign(vctr);
			
				return *this;
			}

			template <class U>
			Vctr<T, edev_cpu>& operator=(const std::vector<U>& vctr)
			{
				resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
				memcpy_cpu_cpu(m_data, vctr.data(), vctr.size());

				return *this;
			}

	#ifdef __CUDACC__
			template <class U>
			Vctr<T, edev_cpu>& operator=(const Vctr<U, edev_gpu>& vctr)
			{
				assign(vctr);
			
				return *this;
			}

			template <class U>
			Vctr<T, edev_cpu>& operator=(const thrust::host_vector<U>& vctr)
			{
				assign(vctr);
			
				return *this;
			}

			template <class U>
			Vctr<T, edev_cpu>& operator=(const thrust::device_vector<U>& vctr)
			{
				assign(vctr);
			
				return *this;
			}
	#endif

			template <class U>
			void assign(const Vctr<U, edev_cpu>& vctr, U* pvctr_cpu = nullptr)
			{
				if ((void*)this != (void*)&vctr)
				{
					this->allocate(vctr.shape());
					vctr_cpy(*this, vctr, pvctr_cpu);
				}
			}

			template <class U>
			void assign(U* first, U* last)
			{
				if ((void*)m_data != (void*)first)
				{
					auto m_size_i = std::distance(first, last);
					if (m_size_i>0)
					{
						auto m_size = size_type(m_size_i);
						resize({ m_size, size_type(1), size_type(1), size_type(1) });
						memcpy_cpu_cpu(m_data, first, m_size);
					}
				}
			}

			// from cpu pVctr to Vctr
			template <class U, class STU>
			void assign(const pVctr<U, edev_cpu, STU>& pvctr)
			{
				resize(pvctr.shape());
				memcpy_cpu_cpu(m_data, pvctr.data(), m_size);
			}

	#ifdef __CUDACC__
			template <class U>
			void assign(const thrust::device_ptr<U>& first, const thrust::device_ptr<U>& last, U* pvctr_cpu = nullptr)
			{
				auto m_size_i = thrust::distance(first, last);

				if (m_size_i>0)
				{
					auto m_size = size_type(m_size_i);
					resize({ m_size, size_type(1), size_type(1), size_type(1) });
					memcpy_gpu_cpu(m_data, (U*)first.get(), m_size, pvctr_cpu);
				}
			}

			// from gpu pVctr to Vctr
			template <class U, class STU>
			void assign(const pVctr<U, edev_gpu, STU>& pvctr)
			{
				resize(pvctr.shape());
				memcpy_gpu_cpu(m_data, pvctr.data(), m_size);
			}
	#endif

			template <class U>
			void assign(const std::vector<U>& vctr, U* pvctr_cpu = nullptr)
			{
				resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
				memcpy_cpu_cpu(m_data, vctr.data(), vctr.size(), pvctr_cpu);
			}

	#ifdef __CUDACC__
			template <class U>
			void assign(const Vctr<U, edev_gpu>& vctr, U* pvctr_cpu = nullptr)
			{
				this->allocate(vctr.shape());
				vctr_cpy(*this, vctr, pvctr_cpu);
			}

			template <class U>
			void assign(const thrust::host_vector<U>& vctr, U* pvctr_cpu = nullptr)
			{
				resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
				memcpy_cpu_cpu(m_data, vctr.data(), vctr.size(), pvctr_cpu);
			}

			template <class U>
			void assign(const thrust::device_vector<U>& vctr, U* pvctr_cpu = nullptr)
			{
				resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
				memcpy_gpu_cpu(m_data, (U*)vctr.data().get(), vctr.size(), pvctr_cpu);
			}
	#endif
			/**************** user define conversion operators *******************/
			pVctr_cpu_32<T> ptr_32() const
			{
				return pVctr_cpu_32<T>(*this);
			}

			pVctr_cpu_64<T> ptr_64() const
			{
				return pVctr_cpu_64<T>(*this);
			}

			operator pVctr_cpu_32<T>() const
			{
				return pVctr_cpu_32<T>(*this);
			}

			 operator pVctr_cpu_64<T>() const
			{
				return pVctr_cpu_64<T>(*this);
			}

			// ! user define conversion
			operator std::vector<T>() const
			{
				std::vector<T> vctr(m_size);
				memcpy_cpu_cpu(vctr.data(), m_data, m_size);

				return vctr;
			}

			// ! user define conversion in which T is the complemented precision
			template <class U = chg_2_compl_float_type<T>>
			operator std::vector<U>() const
			{
				std::vector<U> vctr(m_size);
				memcpy_cpu_cpu(vctr.data(), m_data, m_size);

				return vctr;
			}

	#ifdef __CUDACC__
			// ! user define conversion to output type std::vector<thrust::complex<dt_float32>>
			template <class U=T, typename = enable_if_std_cfloat<U>>
			operator std::vector<thrust::complex<dt_float32>>() const
			{
				std::vector<thrust::complex<dt_float32>> vctr(m_size);
				memcpy_cpu_cpu(vctr.data(), m_data, m_size);

				return vctr;
			}

			// ! user define conversion to output type std::vector<thrust::complex<dt_float64>>
			template <class U=T, typename = enable_if_std_cfloat<U>>
			operator std::vector<thrust::complex<dt_float64>>() const
			{
				std::vector<thrust::complex<dt_float64>>vctr(m_size);
				memcpy_cpu_cpu(vctr.data(), m_data, m_size);

				return vctr;
			}

			/***************************************************************************************/
			// ! user define conversion
			operator thrust::host_vector<T>() const
			{
				thrust::host_vector<T> vctr(m_size);
				memcpy_cpu_cpu(vctr.data(), m_data, m_size);

				return vctr;
			}

			// ! user define conversion in which T is the complemented precision
			template <class U = chg_2_compl_float_type<T>>
			operator thrust::host_vector<U>() const
			{
				thrust::host_vector<U> vctr(m_size);
				memcpy_cpu_cpu(vctr.data(), m_data, m_size);

				return vctr;
			}

			// ! user define conversion to output type thrust::host_vector<std::complex<dt_float32>>
			template <class U=T, typename = enable_if_thr_cfloat<U>>
			operator thrust::host_vector<std::complex<dt_float32>>() const
			{
				thrust::host_vector<std::complex<dt_float32>> vctr(m_size);
				memcpy_cpu_cpu(vctr.data(), m_data, m_size);

				return vctr;
			}

			// ! user define conversion to output type thrust::host_vector<std::complex<dt_float64>>
			template <class U=T, typename = enable_if_thr_cfloat<U>>
			operator thrust::host_vector<std::complex<dt_float64>>() const
			{
				thrust::host_vector<std::complex<dt_float64>> vctr(m_size);
				memcpy_cpu_cpu(vctr.data(), m_data, m_size);

				return vctr;
			}

			/***************************************************************************************/
			// ! user define conversion
			operator thrust::device_vector<T>() const
			{
				thrust::device_vector<T> vctr(m_size);
				memcpy_cpu_gpu((T*)vctr.data().get(), m_data, m_size);

				return vctr;
			}

			// ! user define conversion in which T is the complemented precision
			template <class U = chg_2_compl_float_type<T>>
			operator thrust::device_vector<U>() const
			{
				thrust::device_vector<U> vctr(m_size);
				memcpy_cpu_gpu((U*)vctr.data().get(), m_data, m_size);

				return vctr;
			}
	#endif
			/***************************************************************************************/
			template <class U>
			void cpy_to_cpu_ptr(U* pdata, size_type n_data, T* pvctr_cpu = nullptr)
			{
				n_data = std::min(m_size, n_data);
				memcpy_cpu_cpu(pdata, m_data, n_data, pvctr_cpu);
			}

			template <class U>
			void cpy_to_cpu_ptr(U* first, U* last, T* pvctr_cpu = nullptr)
			{
				auto n_data = std::distance(first, last);
				memcpy_cpu_cpu(first, m_data, n_data, pvctr_cpu);
			}

	#ifdef __CUDACC__
			template <class U>
			void cpy_to_gpu_ptr(U* pdata, size_type n_data, U* pvctr_cpu = nullptr)
			{
				n_data = std::min(m_size, n_data);
				memcpy_cpu_gpu(pdata, m_data, n_data, pvctr_cpu);
			}

			template <class U>
			void cpy_to_gpu_ptr(U* first, U* last, U* pvctr_cpu = nullptr)
			{
				auto n_data = std::distance(first, last);
				memcpy_cpu_gpu(first, m_data, n_data, pvctr_cpu);
			}
	#endif

			/***************************************************************************************/
			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_real_to_cpu_ptr(U* pdata, size_type n_data, T* pvctr_cpu = nullptr)
			{
				n_data = std::min(m_size, n_data);
				memcpy_real_cpu_cpu(pdata, m_data, n_data, pvctr_cpu);
			}

			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_real_to_cpu_ptr(U* first, U* last, T* pvctr_cpu = nullptr)
			{
				auto n_data = std::distance(first, last);
				memcpy_real_cpu_cpu(first, m_data, n_data, pvctr_cpu);
			}

	#ifdef __CUDACC__
			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_real_to_gpu_ptr(U* pdata, size_type n_data, U* pvctr_cpu = nullptr)
			{
				n_data = std::min(m_size, n_data);
				memcpy_real_cpu_gpu(pdata, m_data, n_data, pvctr_cpu);
			}

			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_real_to_gpu_ptr(U* first, U* last, U* pvctr_cpu = nullptr)
			{
				auto n_data = std::distance(first, last);
				memcpy_real_cpu_gpu(first, m_data, n_data, pvctr_cpu);
			}
	#endif

			/***************************************************************************************/
			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_imag_to_cpu_ptr(U* pdata, size_type n_data, T* pvctr_cpu = nullptr)
			{
				n_data = std::min(m_size, n_data);
				memcpy_imag_cpu_cpu(pdata, m_data, n_data, pvctr_cpu);
			}

			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_imag_to_cpu_ptr(U* first, U* last, T* pvctr_cpu = nullptr)
			{
				auto n_data = std::distance(first, last);
				memcpy_imag_cpu_cpu(first, m_data, n_data, pvctr_cpu);
			}

	#ifdef __CUDACC__
			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_imag_to_gpu_ptr(U* pdata, size_type n_data, U* pvctr_cpu = nullptr)
			{
				n_data = std::min(m_size, n_data);
				memcpy_imag_cpu_gpu(pdata, m_data, n_data, pvctr_cpu);
			}

			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_imag_to_gpu_ptr(U* first, U* last, U* pvctr_cpu = nullptr)
			{
				auto n_data = std::distance(first, last);
				memcpy_imag_cpu_gpu(first, m_data, n_data, pvctr_cpu);
			}
	#endif

			/***************************************************************************************/
			template <class U, class V=T, typename = enable_if_r_nd<V>>
			void set_pos_from_cpu_ptr(U *p, size_type n_p, dt_int64 icol=0)
			{
				resize({n_p});

				memcpy_pos_cpu_cpu(m_data, p, size(), icol);
			}

			template <class U, class V=T, typename = enable_if_r_nd<V>>
			void cpy_pos_to_cpu_ptr(U *p, size_type icol=0)
			{
				memcpy_pos_cpu_cpu(p, m_data, size(), icol);
			}

			/***************************************************************************************/
			void resize(const dt_shape_st<size_type>& shape)
			{
				this->allocate(shape);
			}

			void resize(const dt_shape_st<size_type>& shape, const T& value)
			{
				auto m_size_t = m_size;
				this->allocate(shape);
				std::fill(this->begin()+m_size_t, this->end(), value);
			}

			void reserve(const dt_shape_st<size_type>& shape)
			{
				this->allocate(shape, true);
			}

			void shrink_to_fit()
			{
				if (m_size < m_capacity)
				{
					if (m_size > 0)
					{
						// ! allocate memory and transfer data
						T* data_tmp = new T[m_size];
						memcpy_cpu_cpu(data_tmp, m_data, m_size);
						delete[] m_data;

						m_data = data_tmp;
						data_tmp = nullptr;
					}
					else
					{
						delete[] m_data;
						m_data = nullptr;
					}

					m_capacity = m_size;

					if (shape().dim()==2)
					{
						if (m_s1==1)
							m_s0 = m_size;

						if (m_s0==1)
							m_s1 = m_size;
					}
				}
			}

			template<class U>
			void push_back(const U& val)
			{
				if (m_size >= m_capacity)
				{
					dt_shape_st<size_type> shape{m_s0, m_s1, m_s2, m_s3};

					if (m_s3 == 1)
					{
						if (m_s2 == 1)
						{
							if (m_s1 == 1)
								shape[0] = max(size_type(1), 2*m_size);
							else
								shape[1] *= size_type(2);
						}
						else
						{
							shape[2] *= size_type(2);
						}
					}
					else
					{
						shape[3] *= size_type(2);
					}

					this->allocate(shape, true);
				}

				m_data[m_size++] = T(val);
			}

			template<class U>
			void push_back(const Vctr<U, edev_cpu>& vctr)
			{
				const size_type size_new = m_size + vctr.size();

				if (size_new >= m_capacity)
				{
					dt_shape_st<size_type> shape{m_s0, m_s1, m_s2, m_s3};

					if (m_s3 == 1)
					{
						if (m_s2 == 1)
						{
							if (m_s1 == 1)
								shape[0] = max(size_new, 2*shape[0]);
							else
								shape[1] = max(size_new/shape[0], 2*shape[1]);
						}
						else
						{
							shape[2] = max(size_new/(shape[0]*shape[1]), 2*shape[2]);
						}
					}
					else
					{
						shape[3] = max(size_new/(shape[0]*shape[1]*shape[2]), 2*shape[3]);
					}

					this->allocate(shape, true);
				}

				memcpy_cpu_cpu(m_data + m_size, vctr.data(), vctr.size());

				m_size = size_new;
			}

			void pop_back()
			{
				if (m_size >= 1)
				{
					m_size--;
				}
			}

			void fill(T val)
			{
				std::fill(this->begin(), this->end(), val);
			}

			/***************************************************************************************/
			size_type s0() const
			{
				return m_s0;
			}

			size_type s1() const
			{
				return m_s1;
			}

			size_type s2() const
			{
				return m_s2;
			}

			size_type s3() const
			{
				return m_s2;
			}			
			
			dt_int32 s0_32() const
			{
				return static_cast<dt_int32>(m_s0);
			}

			dt_int32 s1_32() const
			{
				return static_cast<dt_int32>(m_s1);
			}

			dt_int32 s2_32() const
			{
				return static_cast<dt_int32>(m_s2);
			}

			dt_int32 s3_32() const
			{
				return static_cast<dt_int32>(m_s3);
			}
			
			dt_int64 s0_64() const
			{
				return static_cast<dt_int64>(m_s0);
			}

			dt_int64 s1_64() const
			{
				return static_cast<dt_int64>(m_s1);
			}

			dt_int64 s2_64() const
			{
				return static_cast<dt_int64>(m_s2);
			}

			dt_int64 s3_64() const
			{
				return static_cast<dt_int64>(m_s3);
			}

			size_type s0h() const
			{
				return m_s0/size_type(2);
			}

			size_type s1h() const
			{
				return m_s1/size_type(2);
			}

			size_type s2h() const
			{
				return m_s2/size_type(2);
			}

			size_type s3h() const
			{
				return m_s3/size_type(2);
			}

			dt_shape_st<size_type> shape() const
			{
				return {m_s0, m_s1, m_s2, m_s3};
			}

			dt_shape_st<size_type> shape_2d_trs() const
			{
				return {m_s1, m_s0, m_s2, m_s3};
			}

			size_type shape_size() const 
			{
				return m_s0*max(m_s1, size_type(1))*max(m_s2, size_type(1))*max(m_s3, size_type(1));
			}

			size_type pitch_s1() const
			{
				return m_pitch_s1;
			}

			size_type pitch_s2() const
			{
				return m_pitch_s2;
			}

			size_type pitch_s3() const
			{
				return m_pitch_s3;
			}

			size_type size() const
			{
				return m_size;
			}

			dt_int32 size_32() const
			{
				return static_cast<dt_int32>(m_size);
			}

			dt_int64 size_64() const
			{
				return static_cast<dt_int64>(m_size);
			}

			iGrid_1d igrid_1d() const
			{
				return {size_32()};
			}

			iGrid_2d igrid_2d() const
			{
				return {s1_32(), s0_32()};
			}

			iGrid_2d igrid_3d() const
			{
				return {s1_32(), s0_32(), s2_32()};
			}

			iGrid_1d_64 igrid_1d_64() const
			{
				return {size_64()};
			}

			iGrid_2d_64 igrid_2d_64() const
			{
				return {s1_64(), s0_64()};
			}

			iGrid_2d_64 igrid_3d_64() const
			{
				return {s1_64(), s0_64(), s2_64()};
			}

			size_type capacity() const
			{
				return m_capacity;
			}

			dt_bool empty() const
			{
				return m_size == 0;
			}

			dt_bool is_1d() const
			{
				return (m_s0 == 1) || (m_s1 == 1);
			}

			void clear()
			{
				m_size = 0;
			}

			void clear_shrink_to_fit()
			{
				destroy();
			}

			size_type sub_2_ind(const size_type& ix_0) const 
			{ 
				return ix_0;
			}

			size_type sub_2_ind(const size_type& ix_0, const size_type& ix_1) const 
			{ 
				return ix_0 + m_s0*ix_1;
			}

			size_type sub_2_ind(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2) const 
			{ 
				return ix_0 + m_s0*(ix_1 + m_s1*ix_2);
			}

			size_type sub_2_ind(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2, const size_type& ix_3) const 
			{ 
				return ix_0 + m_s0*(ix_1 + m_s1*(ix_2 + m_s2*ix_3));
			}

			T& operator[](const size_type& iy)
			{ 
				return m_data[iy];
			}

			const T& operator[](const size_type& iy) const 
			{ 
				return m_data[iy];
			}

			T& operator()(const size_type& iy)
			{ 
				return m_data[iy];
			}

			const T& operator()(const size_type& iy) const 
			{ 
				return m_data[iy];
			}

			T& operator()(const size_type& ix_0, const size_type& ix_1)
			{ 
				return m_data[sub_2_ind(ix_0, ix_1)];
			}

			const T& operator()(const size_type& ix_0, const size_type& ix_1) const 
			{ 
				return m_data[sub_2_ind(ix_0, ix_1)];
			}

			T& operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2)
			{ 
				return m_data[sub_2_ind(ix_0, ix_1, ix_2)];
			}

			const T& operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2) const 
			{ 
				return m_data[sub_2_ind(ix_0, ix_1, ix_2)];
			}

			T& operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2, const size_type& ix_3)
			{ 
				return m_data[sub_2_ind(ix_0, ix_1, ix_2, ix_3)];
			}

			const T& operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2, const size_type& ix_3) const 
			{ 
				return m_data[sub_2_ind(ix_0, ix_1, ix_2, ix_3)];
			}

			T* begin()
			{ 
				return m_data;
			}

			const T* begin() const 
			{ 
				return m_data;
			}

			T* end()
			{ 
				return m_data + m_size;
			}

			const T* end() const
			{ 
				return m_data + m_size;
			}

			T* data()
			{
				return m_data;
			}

			const T* data() const
			{
				return m_data;
			}

			template <class U>
			U data_cast()
			{
				return reinterpret_cast<U>(m_data);
			}

			template <class U>
			const U data_cast() const
			{
				return reinterpret_cast<U>(m_data);
			}

			T& front()
			{ 
				return m_data[0];
			}

			const T& front() const 
			{ 
				return m_data[0];
			}

			T& back()
			{ 
				return m_data[m_size-1];
			}

			const T& back() const 
			{ 
				return m_data[m_size-1];
			}

			// set shape
			void set_shape(const dt_shape_st<size_type>& shape)
			{
				m_s0 = max(shape[0], size_type(0));
				m_s1 = max(shape[1], size_type(1));
				m_s2 = max(shape[2], size_type(1));
				m_s3 = max(shape[3], size_type(1));

				set_picth();
			}

			void trs_shape_2d()
			{
				set_shape({m_s1, m_s0, m_s2, m_s3});
			}

			template <class U, eDev Dev>
			void set_shape(const Vctr<U, Dev>& vctr, dt_bool bb_size = true)
			{
				this->set_shape(vctr.shape(), bb_size);
			}

	#ifdef __CUDACC__
			FCNS_GPU_GRID_BLK_VCTR;
	#endif

		private:

			void set_picth()
			{
				m_pitch_s1 = m_s0*sizeof(T);
				m_pitch_s2 = m_pitch_s1*m_s1;
				m_pitch_s3 = m_pitch_s2*m_s2;
			}

			void set_capacity(size_type size_r)
			{
				m_capacity = max(m_capacity, size_r);
			}

			void set_shape_cstr(dt_shape_st<size_type>& shape)
			{
				shape[0] = max(shape[0], size_type(0));
				shape[1] = max(shape[1], size_type(1));
				shape[2] = max(shape[2], size_type(1));
				shape[3] = max(shape[3], size_type(1));
			}

			// reallocate and copy memory
			void allocate(dt_shape_st<size_type> shape, dt_bool bb_reserve=false)
			{
				set_shape_cstr(shape);

				auto size_r = shape[0]*shape[1]*shape[2]*shape[3];
			
				if (size_r == 0)
				{
					return;
				}

				if (m_size == 0)
				{
					if (m_data != nullptr)
					{
						delete[] m_data;
					}

					m_data = new T[size_r];
				}
				else if (size_r>m_capacity)
				{
					// ! allocate	memory and transfer data
					T* data_tmp = new T[size_r];
					memcpy_cpu_cpu(data_tmp, m_data, m_size);

					if (m_data != nullptr)
					{
						delete[] m_data;
					}

					m_data = data_tmp;
					data_tmp = nullptr;
				}

				this->set_shape(shape);
				if (!bb_reserve)
				{
					m_size = size_r;
				}
				this->set_capacity(size_r);
			}

			// destroy memory on the device
			void init()
			{
				m_data = nullptr;
				m_s0 = 0; m_s1 = m_s2 = m_s3 = 1;
				m_pitch_s1 = m_pitch_s2 = m_pitch_s3 = 0;
				m_size = 0;
				m_capacity = 0;
			}

			// destroy memory on the device
			void destroy()
			{
				if (m_data != 0)
				{
					delete[] m_data;
				}

				init();
			}

		};
	}

	/* gpu vector */
	namespace mt
	{
	#ifdef __CUDACC__
		template <class T>
		class Vctr<T, edev_gpu>
		{
		public:
			using value_type = T;
			using size_type = dt_int64;
			static const eDev device = edev_gpu;

			mutable T* m_data;
			size_type m_s0;
			size_type m_s1;
			size_type m_s2;
			size_type m_s3;
			size_type m_size;
			size_type m_capacity;

			size_type m_pitch_s1;
			size_type m_pitch_s2;
			size_type m_pitch_s3;

			/************************************* constructors ************************************/
			explicit Vctr(): m_data(nullptr), m_s0(0), m_s1(1), m_s2(1), m_s3(1), 
			m_size(shape_size()), m_capacity(0)
			{
				set_picth();
			}

			Vctr(const dt_init_list_f64& data): Vctr()
			{
				assign(data.begin(), data.end());
			}

			Vctr(size_type s0): Vctr()
			{
				resize({ s0, size_type(1), size_type(1), size_type(1) });
			}

			Vctr(size_type s0, const T& value): Vctr()
			{
				resize({ s0, size_type(1), size_type(1), size_type(1)}, value);
			}

			explicit Vctr(const dt_shape_st<size_type>& shape): Vctr()
			{
				resize(shape);
			}

			explicit Vctr(const dt_shape_st<size_type>& shape, const T& value): Vctr()
			{
				resize(shape, value);
			}

			/* copy constructor */
			Vctr(const Vctr<T, edev_gpu>& vctr): Vctr()
			{
				*this = vctr;
			}

			/* Move constructor */
			Vctr(Vctr<T, edev_gpu>&& vctr): Vctr()
			{
				*this = std::move(vctr);
			}

			/* converting constructor */
			template <class U>
			Vctr(const Vctr<U, edev_gpu>& vctr): Vctr()
			{
				assign(vctr);
			}

			template <class U>
			Vctr(U* first, U* last): Vctr()
			{
				assign(first, last);
			}

			template <class U, class V=T, typename = enable_if_r_nd<V>>
			Vctr(U *p, dt_int64 n_p, size_type icol=0): Vctr()
			{
				set_pos_from_cpu_ptr(p, n_p, icol);
			}

			template <class U>
			Vctr(const std::vector<U>& vctr): Vctr()
			{
				assign(vctr);
			}

			// from gpu pVctr to Vctr
			template <class U, class STU>
			Vctr(const pVctr<U, edev_gpu, STU>& pvctr): Vctr()
			{
				assign(pvctr);
			}

			// from cpu pVctr to Vctr
			template <class U, class STU>
			Vctr(const pVctr<U, edev_cpu, STU>& pvctr): Vctr()
			{
				assign(pvctr);
			}

			template <class U>
			Vctr(const Vctr<U, edev_cpu>& vctr): Vctr()
			{
				assign(vctr);
			}

			template <class U>
			Vctr(const thrust::host_vector<U>& vctr): Vctr()
			{
				assign(vctr);
			}

			template <class U>
			Vctr(const thrust::device_vector<U>& vctr): Vctr()
			{
				assign(vctr);
			}

			~Vctr()
			{
				destroy();
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			Vctr<T, edev_gpu>& operator=(const Vctr<T, edev_gpu>& vctr)
			{
				if (vctr.m_data == nullptr)
				{
					destroy();
				}
				else if (this != &vctr)
				{
					fcn_cuda_free(m_data);

					fcn_cuda_malloc(m_data, vctr.m_capacity);

					memcpy_gpu_gpu(m_data, vctr.data(), vctr.size());

					m_s0 = vctr.m_s0;
					m_s1 = vctr.m_s1;
					m_s2 = vctr.m_s2;
					m_s3 = vctr.m_s3;
					m_size = vctr.m_size;
					m_capacity = vctr.m_capacity;	
					m_pitch_s1 = vctr.m_pitch_s1;
					m_pitch_s2 = vctr.m_pitch_s2;
					m_pitch_s3 = vctr.m_pitch_s3;
				}

				return *this;
			}

			/* Move assignment operator */
			Vctr<T, edev_gpu>& operator=(Vctr<T, edev_gpu>&& vctr)
			{
				if (vctr.m_data == nullptr)
				{
					destroy();
				}
				else if (this != &vctr)
				{
					fcn_cuda_free(m_data);

					m_data = vctr.m_data;
					m_s0 = vctr.m_s0;
					m_s1 = vctr.m_s1;
					m_s2 = vctr.m_s2;
					m_s3 = vctr.m_s3;
					m_size = vctr.m_size;
					m_capacity = vctr.m_capacity;
					m_pitch_s1 = vctr.m_pitch_s1;
					m_pitch_s2 = vctr.m_pitch_s2;
					m_pitch_s3 = vctr.m_pitch_s3;

					vctr.m_data = nullptr;
					vctr.m_s0 = 0;
					vctr.m_s1 = 0;
					vctr.m_s2 = 0;
					vctr.m_s3 = 0;
					vctr.m_size = 0;
					vctr.m_capacity = 0;
					vctr.m_pitch_s1 = 0;
					vctr.m_pitch_s2 = 0;
					vctr.m_pitch_s3 = 0;
				}

				return *this;
			}

			/* converting assignment operator */
			template <class U>
			Vctr<T, edev_gpu>& operator=(const Vctr<U, edev_gpu>& vctr)
			{
				assign(vctr);
			
				return *this;
			}

			template <class U>
			Vctr<T, edev_gpu>& operator=(const std::vector<U>& vctr)
			{
				resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
				memcpy_cpu_gpu(m_data, vctr.data(), vctr.size());

				return *this;
			}

			template <class U>
			Vctr<T, edev_gpu>& operator=(const Vctr<U, edev_cpu>& vctr)
			{
				assign(vctr);
			
				return *this;
			}

			template <class U>
			Vctr<T, edev_gpu>& operator=(const thrust::host_vector<U>& vctr)
			{
				assign(vctr);
			
				return *this;
			}

			template <class U>
			Vctr<T, edev_gpu>& operator=(const thrust::device_vector<U>& vctr)
			{
				assign(vctr);
			
				return *this;
			}

			template <class U>
			void assign(const Vctr<U, edev_gpu>& vctr, T* pvctr_cpu = nullptr)
			{
				if ((void*)this != (void*)&vctr)
				{
					this->allocate(vctr.shape());
					vctr_cpy(*this, vctr, pvctr_cpu);
				}
			}

			template <class U>
			void assign(U* first, U* last)
			{
				if ((void*)m_data != (void*)first) 
				{
					auto m_size_i = std::distance(first, last);
					if (m_size_i>0)
					{
						auto m_size = size_type(m_size_i);
						resize({ m_size, size_type(1), size_type(1), size_type(1) });
						memcpy_cpu_gpu(m_data, first, m_size);
					}
				}
			}
			
			// from gpu pVctr to Vctr
			template <class U, class STU>
			void assign(const pVctr<U, edev_gpu, STU>& pvctr)
			{
				resize(pvctr.shape());
				memcpy_gpu_gpu(m_data, pvctr.data(), m_size);
			}

			template <class U>
			void assign(const thrust::device_ptr<U>& first, const thrust::device_ptr<U>& last, T* pvctr_cpu = nullptr)
			{
				auto m_size_i = thrust::distance(first, last);
				if (m_size_i>0)
				{
					auto m_size = size_type(m_size_i);
					resize({ m_size, size_type(1), size_type(1), size_type(1) });
					memcpy_gpu_gpu(m_data, (U*)first.get(), m_size, pvctr_cpu);
				}
			}

			// from cpu pVctr to Vctr
			template <class U, class STU>
			void assign(const pVctr<U, edev_cpu, STU>& pvctr)
			{
				resize(pvctr.shape());
				memcpy_cpu_gpu(m_data, pvctr.data(), m_size);
			}

			template <class U>
			void assign(const std::vector<U>& vctr, T* pvctr_cpu = nullptr)
			{
				resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
				memcpy_cpu_gpu(m_data, vctr.data(), vctr.size(), pvctr_cpu);
			}

			template <class U>
			void assign(const Vctr<U, edev_cpu>& vctr, T* pvctr_cpu = nullptr)
			{
				this->allocate(vctr.shape());
				vctr_cpy(*this, vctr, pvctr_cpu);
			}

			template <class U>
			void assign(const thrust::host_vector<U>& vctr, T* pvctr_cpu = nullptr)
			{
				resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
				memcpy_cpu_gpu(m_data, vctr.data(), vctr.size(), pvctr_cpu);
			}

			template <class U>
			void assign(const thrust::device_vector<U>& vctr, T* pvctr_cpu = nullptr)
			{
				resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
				memcpy_gpu_gpu(m_data, (T*)vctr.data().get(), vctr.size(), pvctr_cpu);
			}

			/**************** user define conversion operators *******************/
			pVctr_gpu_32<T> ptr_32() const
			{
		
				return pVctr_gpu_32<T>(*this);
			}

			pVctr_gpu_64<T> ptr_64() const
			{
		
				return pVctr_gpu_64<T>(*this);
			}

			// ! user define conversion for pointer Vctr
			operator pVctr_gpu_32<T>() const
			{
				return pVctr_gpu_32<T>(*this);
			}

			operator pVctr_gpu_64<T>() const
			{
				return pVctr_gpu_64<T>(*this);
			}

			// ! user define conversion
			operator std::vector<T>() const
			{
				std::vector<T> vctr(m_size);
				memcpy_gpu_cpu(vctr.data(), m_data, m_size);

				return vctr;
			}

			// ! user define conversion in which T is the complemented precision
			template <class U = chg_2_compl_float_type<T>>
			operator std::vector<U>() const
			{
				std::vector<U> vctr(m_size);
				memcpy_gpu_cpu(vctr.data(), m_data, m_size);

				return vctr;
			}

			// ! user define conversion to output type std::vector<thrust::complex<dt_float32>>
			template <class U=T, typename = enable_if_std_cfloat<U>>
			operator std::vector<thrust::complex<dt_float32>>() const
			{
				std::vector<thrust::complex<dt_float32>> vctr(m_size);
				memcpy_gpu_cpu(vctr.data(), m_data, m_size);
				return vctr;
			}

			// ! user define conversion to output type std::vector<thrust::complex<dt_float64>>
			template <class U=T, typename = enable_if_std_cfloat<U>>
			operator std::vector<thrust::complex<dt_float64>>() const
			{
				std::vector<thrust::complex<dt_float64>>vctr(m_size);
				memcpy_gpu_cpu(vctr.data(), m_data, m_size);
				return vctr;
			}

			// ! user define conversion
			operator thrust::host_vector<T>() const
			{
				thrust::host_vector<T> vctr(m_size);
				memcpy_gpu_cpu(vctr.data(), m_data, m_size);

				return vctr;
			}

			// ! user define conversion in which T is the complemented precision
			template <class U = chg_2_compl_float_type<T>>
			operator thrust::host_vector<U>() const
			{
				thrust::host_vector<U> vctr(m_size);
				memcpy_gpu_cpu(vctr.data(), m_data, m_size);

				return vctr;
			}

			// ! user define conversion to output type thrust::host_vector<std::complex<dt_float32>>
			template <class U=T, typename = enable_if_thr_cfloat<U>>
			operator thrust::host_vector<std::complex<dt_float32>>() const
			{
				thrust::host_vector<std::complex<dt_float32>> vctr(m_size);
				memcpy_gpu_cpu(vctr.data(), m_data, m_size);

				return vctr;
			}

			// ! user define conversion to output type thrust::host_vector<std::complex<dt_float64>>
			template <class U=T, typename = enable_if_thr_cfloat<U>>
			operator thrust::host_vector<std::complex<dt_float64>>() const
			{
				thrust::host_vector<std::complex<dt_float64>> vctr(m_size);
				memcpy_gpu_cpu(vctr.data(), m_data, m_size);

				return vctr;
			}

			/***************************************************************************************/
			// ! user define conversion
			operator thrust::device_vector<T>() const
			{
				thrust::host_vector<T> vctr(m_size);
				memcpy_gpu_gpu((T*)vctr.data().get(), m_data, m_size);
				return vctr;
			}

			// ! user define conversion in which T is the complemented precision
			template <class U = chg_2_compl_float_type<T>>
			operator thrust::device_vector<U>() const
			{
				thrust::host_vector<U> vctr(m_size);
				memcpy_gpu_gpu((U*)vctr.data().get(), m_data, m_size);

				return vctr;
			}

			/***************************************************************************************/
			template <class U>
			void cpy_to_cpu_ptr(U* pdata, size_type n_data, T* pvctr_cpu = nullptr)
			{
				n_data = std::min(m_size, n_data);
				memcpy_gpu_cpu(pdata, m_data, n_data, pvctr_cpu);
			}

			template <class U>
			void cpy_to_cpu_ptr(U* first, U* last, T* pvctr_cpu = nullptr)
			{
				auto n_data = std::distance(first, last);
				memcpy_gpu_cpu(first, m_data, n_data, pvctr_cpu);
			}

			template <class U>
			void cpy_to_gpu_ptr(U* pdata, size_type n_data, U* pvctr_cpu = nullptr)
			{
				n_data = std::min(m_size, n_data);
				memcpy_gpu_gpu(pdata, m_data, n_data, pvctr_cpu);
			}

			template <class U>
			void cpy_to_gpu_ptr(U* first, U* last, U* pvctr_cpu = nullptr)
			{
				auto n_data = thrust::distance(first, last);
				memcpy_gpu_gpu(first, m_data, n_data, pvctr_cpu);
			}

			/***************************************************************************************/
			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_real_to_cpu_ptr(U* pdata, size_type n_data, T* pvctr_cpu = nullptr)
			{
				n_data = std::min(m_size, n_data);
				memcpy_real_gpu_cpu(pdata, m_data, n_data, pvctr_cpu);
			}

			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_real_to_cpu_ptr(U* first, U* last, T* pvctr_cpu = nullptr)
			{
				auto n_data = std::distance(first, last);
				memcpy_real_gpu_cpu(first, m_data, n_data, pvctr_cpu);
			}

			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_real_to_gpu_ptr(U* pdata, size_type n_data, U* pvctr_cpu = nullptr)
			{
				n_data = std::min(m_size, n_data);
				memcpy_real_gpu_gpu(pdata, m_data, n_data, pvctr_cpu);
			}

			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_real_to_gpu_ptr(U* first, U* last, U* pvctr_cpu = nullptr)
			{
				auto n_data = std::distance(first, last);
				memcpy_real_gpu_gpu(first, m_data, n_data, pvctr_cpu);
			}

			/***************************************************************************************/
			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_imag_to_cpu_ptr(U* pdata, size_type n_data, T* pvctr_cpu = nullptr)
			{
				n_data = std::min(m_size, n_data);
				memcpy_imag_gpu_cpu(pdata, m_data, n_data, pvctr_cpu);
			}

			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_imag_to_cpu_ptr(U* first, U* last, T* pvctr_cpu = nullptr)
			{
				auto n_data = std::distance(first, last);
				memcpy_imag_gpu_cpu(first, m_data, n_data, pvctr_cpu);
			}

			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_imag_to_gpu_ptr(U* pdata, size_type n_data, U* pvctr_cpu = nullptr)
			{
				n_data = std::min(m_size, n_data);
				memcpy_imag_gpu_gpu(pdata, m_data, n_data, pvctr_cpu);
			}

			template <class U, class V=T, typename = enable_if_cmplx<V>>
			void cpy_imag_to_gpu_ptr(U* first, U* last, U* pvctr_cpu = nullptr)
			{
				auto n_data = std::distance(first, last);
				memcpy_imag_gpu_gpu(first, m_data, n_data, pvctr_cpu);
			}

			/***************************************************************************************/
			template <class U, class V=T, typename = enable_if_r_nd<V>>
			void set_pos_from_cpu_ptr(U *p, size_type n_p, dt_int64 icol=0)
			{
				resize({n_p});

				memcpy_pos_cpu_gpu(m_data, p, size(), icol);
			}

			template <class U, class V=T, typename = enable_if_r_nd<V>>
			void cpy_pos_to_cpu_ptr(U *p, size_type icol=0)
			{
				memcpy_pos_gpu_cpu(p, m_data, size(), icol);
			}

			/***************************************************************************************/
			void resize(const dt_shape_st<size_type>& shape)
			{
				this->allocate(shape);
			}

			void resize(const dt_shape_st<size_type>& shape, const T& value)
			{
				auto m_size_t = m_size;
				this->allocate(shape);
				thrust::fill(this->begin()+m_size_t, this->end(), value);
			}

			void reserve(const dt_shape_st<size_type>& shape)
			{
				this->allocate(shape, true);
			}

			void shrink_to_fit()
			{
				if (m_size < m_capacity)
				{
					if (m_size > 0)
					{
						// ! allocate memory and transfer data
						T* data_tmp = nullptr;
						fcn_cuda_malloc(data_tmp, m_size);
						memcpy_gpu_gpu(data_tmp, m_data, m_size);
						fcn_cuda_free(m_data);

						m_data = data_tmp;
						data_tmp = nullptr;
					}
					else
					{
						fcn_cuda_free(m_data);
					}

					m_capacity = m_size;

					if (shape().dim()==2)
					{
						if (m_s1==1)
							m_s0 = m_size;

						if (m_s0==1)
							m_s1 = m_size;
					}
				}
			}

			template<class U>
			void push_back(const U& val)
			{
				if (m_size >= m_capacity)
				{
					dt_shape_st<size_type> shape{m_s0, m_s1, m_s2, m_s3};

					if (m_s3 == 1)
					{
						if (m_s2 == 1)
						{
							if (m_s1 == 1)
								shape[0] = max(size_type(1), 2*m_size);
							else
								shape[1] *= size_type(2);
						}
						else
						{
							shape[2] *= size_type(2);
						}
					}
					else
					{
						shape[3] *= size_type(2);
					}

					this->allocate(shape, true);
				}

				memcpy_cpu_gpu(m_data + m_size, &val, 1);

				m_size++;
			}

			template<class U>
			void push_back(const Vctr<U, edev_cpu>& vctr)
			{
				const size_type size_new = m_size + vctr.size();

				if (size_new >= m_capacity)
				{
					dt_shape_st<size_type> shape{m_s0, m_s1, m_s2, m_s3};

					if (m_s3 == 1)
					{
						if (m_s2 == 1)
						{
							if (m_s1 == 1)
								shape[0] = max(size_new, 2*shape[0]);
							else
								shape[1] = max(size_new/shape[0], 2*shape[1]);
						}
						else
						{
							shape[2] = max(size_new/(shape[0]*shape[1]), 2*shape[2]);
						}
					}
					else
					{
						shape[3] = max(size_new/(shape[0]*shape[1]*shape[2]), 2*shape[3]);
					}

					this->allocate(shape, true);
				}

				memcpy_cpu_gpu(m_data + m_size, vctr.data(), vctr.size());

				m_size = size_new;
			}

			void pop_back()
			{
				if (m_size >= 1)
				{
					m_size--;
				}
			}

			void fill(T val)
			{
				thrust::fill(this->begin(), this->end(), val);
			}

			/***************************************************************************************/
			size_type s0() const
			{
				return m_s0;
			}

			size_type s1() const
			{
				return m_s1;
			}

			size_type s2() const
			{
				return m_s2;
			}

			size_type s3() const
			{
				return m_s2;
			}

			dt_int32 s0_32() const
			{
				return static_cast<dt_int32>(m_s0);
			}

			dt_int32 s1_32() const
			{
				return static_cast<dt_int32>(m_s1);
			}

			dt_int32 s2_32() const
			{
				return static_cast<dt_int32>(m_s2);
			}

			dt_int32 s3_32() const
			{
				return static_cast<dt_int32>(m_s3);
			}
			
			dt_int64 s0_64() const
			{
				return static_cast<dt_int64>(m_s0);
			}

			dt_int64 s1_64() const
			{
				return static_cast<dt_int64>(m_s1);
			}

			dt_int64 s2_64() const
			{
				return static_cast<dt_int64>(m_s2);
			}

			dt_int64 s3_64() const
			{
				return static_cast<dt_int64>(m_s3);
			}

			size_type s0h() const
			{
				return m_s0/size_type(2);
			}

			size_type s1h() const
			{
				return m_s1/size_type(2);
			}

			size_type s2h() const
			{
				return m_s2/size_type(2);
			}

			size_type s3h() const
			{
				return m_s3/size_type(2);
			}

			dt_shape_st<size_type> shape() const
			{
				return {m_s0, m_s1, m_s2, m_s3};
			}

			dt_shape_st<size_type> shape_2d_trs() const
			{
				return {m_s1, m_s0, m_s2, m_s3};
			}

			size_type shape_size() const 
			{
				return m_s0*max(m_s1, size_type(1))*max(m_s2, size_type(1))*max(m_s3, size_type(1));
			}

			size_type pitch_s1() const
			{
				return m_pitch_s1;
			}

			size_type pitch_s2() const
			{
				return m_pitch_s2;
			}

			size_type pitch_s3() const
			{
				return m_pitch_s3;
			}

			size_type size() const
			{
				return m_size;
			}

			dt_int32 size_32() const
			{
				return static_cast<dt_int32>(m_size);
			}

			dt_int64 size_64() const
			{
				return static_cast<dt_int64>(m_size);
			}

			iGrid_1d igrid_1d() const
			{
				return {size_32()};
			}

			iGrid_2d igrid_2d() const
			{
				return {s1_32(), s0_32()};
			}

			iGrid_2d igrid_3d() const
			{
				return {s1_32(), s0_32(), s2_32()};
			}

			iGrid_1d_64 igrid_1d_64() const
			{
				return {size_64()};
			}

			iGrid_2d_64 igrid_2d_64() const
			{
				return {s1_64(), s0_64()};
			}

			iGrid_2d_64 igrid_3d_64() const
			{
				return {s1_64(), s0_64(), s2_64()};
			}

			size_type capacity() const
			{
				return m_capacity;
			}

			dt_bool empty() const
			{
				return m_size == 0;
			}

			dt_bool is_1d() const
			{
				return (m_s0 == 1) || (m_s1 == 1);
			}

			void clear()
			{
				m_size = 0;
			}

			void clear_shrink_to_fit()
			{
				destroy();
			}

			size_type sub_2_ind(const size_type& ix_0) const 
			{ 
				return ix_0;
			}

			size_type sub_2_ind(const size_type& ix_0, const size_type& ix_1) const 
			{ 
				return ix_0 + m_s0*ix_1;
			}

			size_type sub_2_ind(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2) const 
			{ 
				return ix_0 + m_s0*(ix_1 + m_s1*ix_2);
			}

			size_type sub_2_ind(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2, const size_type& ix_3) const 
			{ 
				return ix_0 + m_s0*(ix_1 + m_s1*(ix_2 + m_s2*ix_3));
			}

			// T& operator[](const size_type& iy)
			// { 
			// 	return m_data[iy];
			// }

			// const T& operator[](const size_type& iy) const 
			// { 
			// 	return m_data[iy];
			// }

			// T& operator()(const size_type& iy)
			// { 
			// 	return m_data[iy];
			// }

			// const T& operator()(const size_type& iy) const 
			// { 
			// 	return m_data[iy];
			// }

			// T& operator()(const size_type& ix_0, const size_type& ix_1)
			// { 
			// 	return m_data[sub_2_ind(ix_0, ix_1)];
			// }

			// const T& operator()(const size_type& ix_0, const size_type& ix_1) const 
			// { 
			// 	return m_data[sub_2_ind(ix_0, ix_1)];
			// }

			// T& operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2)
			// { 
			// 	return m_data[sub_2_ind(ix_0, ix_1, ix_2)];
			// }

			// const T& operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2) const 
			// { 
			// 	return m_data[sub_2_ind(ix_0, ix_1, ix_2)];
			// }

			// T& operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2, const size_type& ix_3)
			// { 
			// 	return m_data[sub_2_ind(ix_0, ix_1, ix_2, ix_3)];
			// }

			// const T& operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2, const size_type& ix_3) const 
			// { 
			// 	return m_data[sub_2_ind(ix_0, ix_1, ix_2, ix_3)];
			// }

			thrust::device_ptr<T> begin() noexcept 
			{ 
				return thrust::device_ptr<T>(m_data);
			}

			const thrust::device_ptr<T> begin() const
			{ 
				return thrust::device_ptr<T>(m_data);
			}

			thrust::device_ptr<T> end() noexcept 
			{ 
				return thrust::device_ptr<T>(m_data) + m_size;
			}

			const thrust::device_ptr<T> end() const
			{ 
				return thrust::device_ptr<T>(m_data) + m_size;
			}

			T* data()
			{
				return m_data;
			}

			const T* data() const
			{
				return m_data;
			}

			template <class U>
			U data_cast()
			{
				return reinterpret_cast<U>(m_data);
			}

			template <class U>
			const U data_cast() const
			{
				return reinterpret_cast<U>(m_data);
			}

			// T& front()
			// { 
			// 	return m_data[0];
			// }

			// const T& front() const 
			// { 
			// 	return m_data[0];
			// }

			// T& back(){ 
			// 	return m_data[m_size-1];
			// }

			// const T& back() const 
			// { 
			// 	return m_data[m_size-1];
			// }

			// set shape
			void set_shape(const dt_shape_st<size_type>& shape)
			{
				m_s0 = max(shape[0], size_type(0));
				m_s1 = max(shape[1], size_type(1));
				m_s2 = max(shape[2], size_type(1));
				m_s3 = max(shape[3], size_type(1));

				set_picth();
			}

			void trs_shape_2d()
			{
				set_shape({m_s1, m_s0, m_s2, m_s3});
			}

			template <class U, eDev Dev>
			void set_shape(const Vctr<U, Dev>& vctr, dt_bool bb_size = true)
			{
				this->set_shape(vctr.shape(), bb_size);
			}

			FCNS_GPU_GRID_BLK_VCTR;

		private:

			void set_picth()
			{
				m_pitch_s1 = m_s0*sizeof(T);
				m_pitch_s2 = m_pitch_s1*m_s1;
				m_pitch_s3 = m_pitch_s2*m_s2;
			}

			void set_capacity(size_type size_r)
			{
				m_capacity = max(m_capacity, size_r);
			}

			void set_shape_cstr(dt_shape_st<size_type>& shape)
			{
				shape[0] = max(shape[0], size_type(0));
				shape[1] = max(shape[1], size_type(1));
				shape[2] = max(shape[2], size_type(1));
				shape[3] = max(shape[3], size_type(1));
			}

			// reallocate and copy memory
			void allocate(dt_shape_st<size_type> shape, dt_bool bb_reserve=false)
			{
				set_shape_cstr(shape);

				auto size_r = shape[0]*shape[1]*shape[2]*shape[3];
			
				if (size_r == 0)
				{
					return;
				}

				if (m_size == 0)
				{
					fcn_cuda_free(m_data);

					fcn_cuda_malloc(m_data, size_r);
				}
				else if (size_r>m_capacity)
				{
					// allocate	memory and transfer data
					T* data_tmp = nullptr;
					fcn_cuda_malloc(data_tmp, size_r);
					memcpy_gpu_gpu(data_tmp, m_data, m_size);
					fcn_cuda_free(m_data);

					m_data = data_tmp;
					data_tmp = nullptr;
				}

				this->set_shape(shape);
				if (!bb_reserve)
				{
					m_size = size_r;
				}
				this->set_capacity(size_r);
			}

			// initialization
			void init()
			{
				m_data = nullptr;
				m_s0 = 0; m_s1 = m_s2 = m_s3 = 1;
				m_pitch_s1 = m_pitch_s2 = m_pitch_s3 = 0;
				m_size = 0;
				m_capacity = 0;
			}

			// destroy memory on the device
			void destroy()
			{
				if (m_data != nullptr)
				{
					fcn_cuda_free(m_data);
				}

				init();
			}
		};
	#endif
	}

	/* derive vectors */
	namespace mt
	{
		template <class T>
		using Vctr_std = std::vector<T>;

		template <class T>
		using Vctr_cpu = Vctr<T, edev_cpu>;

		template <class T>
		using Vctr_gpu = Vctr<T, edev_gpu>;

		/***************************************************************************************/
		/***************************************************************************************/
		using Vctr_uint32_cpu = Vctr_cpu<dt_uint32>;
		using Vctr_uint32_gpu = Vctr_gpu<dt_uint32>;

		using Vctr_int32_cpu = Vctr_cpu<dt_int32>;
		using Vctr_int32_gpu = Vctr_gpu<dt_int32>;

		using Vctr_uint64_cpu = Vctr_cpu<dt_uint64>;
		using Vctr_uint64_gpu = Vctr_gpu<dt_uint64>;

		using Vctr_int64_cpu = Vctr_cpu<dt_int64>;
		using Vctr_int64_gpu = Vctr_gpu<dt_int64>;

		/***************************************************************************************/
		/***************************************************************************************/
		template <class T, eDev Dev>
		using Vctr_r_2d = Vctr<R_2d<T>, Dev>;

		template <class T>
		using Vctr_r_2d_cpu = Vctr<R_2d<T>, edev_cpu>;

		template <class T>
		using Vctr_r_2d_gpu = Vctr<R_2d<T>, edev_gpu>;

		/***************************************************************************************/
		template <class T>
		struct is_vctr_cpu_r_2d: std::integral_constant<dt_bool, is_vctr_cpu<T>::value && is_r_2d<typename T::value_type>::value> {};

		template <class T, class U>
		struct is_vctr_cpu_r_2d_and_vctr_cpu: std::integral_constant<dt_bool, is_vctr_cpu_r_2d<T>::value && is_vctr_cpu<U>::value> {};

		template <class T>
		struct is_vctr_gpu_r_2d: std::integral_constant<dt_bool, is_vctr_gpu<T>::value && is_r_2d<typename T::value_type>::value> {};

		template <class T, class U>
		struct is_vctr_gpu_r_2d_and_vctr_gpu: std::integral_constant<dt_bool, is_vctr_gpu_r_2d<T>::value && is_vctr_gpu<U>::value> {};

		/***************************************************************************************/
		template <class T, class U=void>
		using enable_if_vctr_cpu_r_2d = typename std::enable_if<is_vctr_cpu_r_2d<T>::value, U>::type;

		template <class T, class U, class V=void>
		using enable_if_vctr_cpu_r_2d_and_vctr_cpu = typename std::enable_if<is_vctr_cpu_r_2d_and_vctr_cpu<T, U>::value, V>::type;
		
		template <class T, class U=void>
		using enable_if_vctr_gpu_r_2d = typename std::enable_if<is_vctr_gpu_r_2d<T>::value, U>::type;

		template <class T, class U, class V=void>
		using enable_if_vctr_gpu_r_2d_and_vctr_gpu = typename std::enable_if<is_vctr_gpu_r_2d_and_vctr_gpu<T, U>::value, V>::type;

		/***************************************************************************************/
		/***************************************************************************************/
		template <class T, eDev Dev>
		using Vctr_r_3d = Vctr<R_3d<T>, Dev>;
		
		template <class T>
		using Vctr_r_3d_cpu = Vctr<R_3d<T>, edev_cpu>;

		/***************************************************************************************/
		template <class T>
		struct is_vctr_cpu_r_3d: std::integral_constant<dt_bool, is_vctr_cpu<T>::value && is_r_3d<typename T::value_type>::value> {};

		template <class T, class U>
		struct is_vctr_cpu_r_3d_and_vctr_cpu: std::integral_constant<dt_bool, is_vctr_cpu_r_3d<T>::value && is_vctr_cpu<U>::value> {};

		template <class T>
		struct is_vctr_gpu_r_3d: std::integral_constant<dt_bool, is_vctr_gpu<T>::value && is_r_3d<typename T::value_type>::value> {};

		template <class T, class U>
		struct is_vctr_gpu_r_3d_and_vctr_gpu: std::integral_constant<dt_bool, is_vctr_gpu_r_3d<T>::value && is_vctr_gpu<U>::value> {};

		/***************************************************************************************/
		template <class T, class U=void>
		using enable_if_vctr_cpu_r_3d = typename std::enable_if<is_vctr_cpu_r_3d<T>::value, U>::type;

		template <class T, class U, class V=void>
		using enable_if_vctr_cpu_r_3d_and_vctr_cpu = typename std::enable_if<is_vctr_cpu_r_3d_and_vctr_cpu<T, U>::value, V>::type;
		
		template <class T, class U=void>
		using enable_if_vctr_gpu_r_3d = typename std::enable_if<is_vctr_gpu_r_3d<T>::value, U>::type;

		template <class T, class U, class V=void>
		using enable_if_vctr_gpu_r_3d_and_vctr_gpu = typename std::enable_if<is_vctr_gpu_r_3d_and_vctr_gpu<T, U>::value, V>::type;

		/***************************************************************************************/
		template <class T>
		struct is_vctr_cpu_r_nd: std::integral_constant<dt_bool, is_vctr_cpu<T>::value && is_r_nd<typename T::value_type>::value> {};

		template <class T, class U>
		struct is_vctr_cpu_r_nd_and_vctr_cpu: std::integral_constant<dt_bool, is_vctr_cpu_r_nd<T>::value && is_vctr_cpu<U>::value> {};

		template <class T>
		struct is_vctr_gpu_r_nd: std::integral_constant<dt_bool, is_vctr_gpu<T>::value && is_r_nd<typename T::value_type>::value> {};

		template <class T, class U>
		struct is_vctr_gpu_r_nd_and_vctr_gpu: std::integral_constant<dt_bool, is_vctr_gpu_r_nd<T>::value && is_vctr_gpu<U>::value> {};

		/***************************************************************************************/
		template <class T, class U=void>
		using enable_if_vctr_cpu_r_nd = typename std::enable_if<is_vctr_cpu_r_nd<T>::value, U>::type;

		template <class T, class U=void>
		using enable_if_vctr_gpu_r_nd = typename std::enable_if<is_vctr_gpu_r_nd<T>::value, U>::type;

		template <class T, class U, class V=void>
		using enable_if_vctr_cpu_r_nd_and_vctr_cpu = typename std::enable_if<is_vctr_cpu_r_nd_and_vctr_cpu<T, U>::value, V>::type;

		template <class T, class U, class V=void>
		using enable_if_vctr_gpu_r_nd_and_vctr_cpu = typename std::enable_if<is_vctr_gpu_r_nd_and_vctr_gpu<T, U>::value, V>::type;

		/***************************************************************************************/
		template <class T, eDev Dev>
		using Vctr_Mx_2x2 = Vctr<Mx_2x2<T>, Dev>;

		template <class T, eDev Dev>
		using Vctr_Mx_3x3 = Vctr<Mx_3x3<T>, Dev>;
	}
#endif