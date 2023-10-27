/*
 * This file is part of Multem.
 * Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
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

#include "pvctr.h"

/* macro definition grid-block */
namespace mt
{
#define	FCNS_IMP_GPU_GRID_BLK_VCTR(TDEF, TEVAL)													\
	TDEF																						\
 	dim3 TEVAL::d_blk_size()																	\
	{																							\
		return fcn_cdb_size();																	\
	}																							\
																								\
 	TDEF																						\
	dim3 TEVAL::d_grid_size(const dim3 d_grid_max)												\
	{																							\
		auto grid = fcn_cdg_size(m_size);														\
																								\
		grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;							\
																								\
		return grid;																			\
	}																							\
																								\
	/***************************************************************************************/	\
 	TDEF																						\
	dim3 TEVAL::d_blk_1d()																		\
	{																							\
		return fcn_cdb_1d();																	\
	}																							\
																								\
	TDEF																						\
 	dim3 TEVAL::d_grid_1d(const dim3 d_grid_max)												\
	{																							\
		auto grid = fcn_cdg_1d(m_s0);															\
																								\
		grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;							\
																								\
		return grid;																			\
	}																							\
																								\
	TDEF																						\
 	D_Grid_Blk TEVAL::d_grid_blk_1d(const dim3 d_grid_max)										\
	{																							\
		return D_Grid_Blk(d_grid_1d(d_grid_max), d_blk_1d());									\
	}																							\
																								\
	/***************************************************************************************/	\
	TDEF																						\
 	dim3 TEVAL::d_grid_1d_h(const dim3 d_grid_max)												\
	{																							\
		auto grid = fcn_cdg_1d(m_s0/2);															\
																								\
		grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;							\
																								\
		return grid;																			\
	}																							\
																								\
	TDEF																						\
	D_Grid_Blk TEVAL::d_grid_blk_1d_h(const dim3 d_grid_max)									\
	{																							\
		return D_Grid_Blk(d_grid_1d_h(d_grid_max), d_blk_1d());									\
	}																							\
																								\
	/***************************************************************************************/	\
	TDEF																						\
 	dim3 TEVAL::d_blk_2d()																		\
	{																							\
		return fcn_cdb_2d();																	\
	}																							\
																								\
	/***************************************************************************************/	\
	TDEF																						\
 	dim3 TEVAL::d_grid_2d(const dim3 d_grid_max)												\
	{																							\
		auto grid = fcn_cdg_2d(m_s0, m_s1);														\
																								\
		grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;							\
		grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;							\
																								\
		return grid;																			\
	}																							\
																								\
	TDEF																						\
 	D_Grid_Blk TEVAL::d_grid_blk_2d(const dim3 d_grid_max)										\
	{																							\
		return D_Grid_Blk(d_grid_2d(d_grid_max), d_blk_2d());									\
	}																							\
																								\
	/***************************************************************************************/	\
	TDEF																						\
 	dim3 TEVAL::d_grid_2d_h(const dim3 d_grid_max)												\
	{																							\
		auto grid = fcn_cdg_2d(m_s0/2, m_s1/2);													\
																								\
		grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;							\
		grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;							\
																								\
		return grid;																			\
	}																							\
																								\
	TDEF																						\
	D_Grid_Blk TEVAL::d_grid_blk_2d_h(const dim3 d_grid_max)									\
	{																							\
		return D_Grid_Blk(d_grid_2d_h(d_grid_max), d_blk_2d());									\
	}																							\
																								\
	/***************************************************************************************/	\
 	TDEF																						\
	dim3 TEVAL::d_blk_3d()																		\
	{																							\
		return fcn_cdb_3d();																	\
	}																							\
																								\
	/***************************************************************************************/	\
	TDEF																						\
 	dim3 TEVAL::d_grid_3d(const dim3 d_grid_max)												\
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
 	TDEF																						\
	D_Grid_Blk TEVAL::d_grid_blk_3d(const dim3 d_grid_max)										\
	{																							\
		return D_Grid_Blk(d_grid_3d(d_grid_max), d_blk_3d());									\
	}																							\
																								\
	/***************************************************************************************/	\
 	TDEF																						\
	dim3 TEVAL::d_grid_3d_h(const dim3 d_grid_max)												\
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
	TDEF																						\
	D_Grid_Blk TEVAL::d_grid_blk_h(const dim3 d_grid_max)										\
	{																							\
		return D_Grid_Blk(d_grid_3d_h(d_grid_max), d_blk_3d());									\
	}
}

/* vector pointer */
namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	pVctr<T, Dev, ST>::pVctr(): m_data(nullptr), m_s0(0), m_s1(1), m_s2(1), m_s3(1), 
	m_size(0)
	{
		set_picth();
	}

	template <class T, eDev Dev, class ST>
	CPU_EXEC
	pVctr<T, Dev, ST>::pVctr(T* data, dt_shape_64 shape): m_data(data), 
	m_s0(ST(shape[0])), m_s1(ST(shape[1])), 
	m_s2(ST(shape[2])), m_s3(ST(shape[3])), m_size(shape_size())
	{
		set_picth();
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	pVctr<T, Dev, ST>::pVctr(T* data, ST s0): m_data(data), 
	m_s0(s0), m_s1(1), m_s2(1), m_s3(1), m_size(shape_size())
	{
		set_picth();
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	pVctr<T, Dev, ST>::pVctr(T* data, ST s0, ST s1): m_data(data), 
	m_s0(s0), m_s1(s1), m_s2(1), m_s3(1), m_size(shape_size())
	{
		set_picth();
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	pVctr<T, Dev, ST>::pVctr(T* data, ST s0, ST s1, ST s2): m_data(data), 
	m_s0(s0), m_s1(s1), m_s2(s2), m_s3(1), m_size(shape_size())
	{
		set_picth();
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	pVctr<T, Dev, ST>::pVctr(T* data, ST s0, ST s1, ST s2, ST s3): m_data(data), 
	m_s0(s0), m_s1(s1), m_s2(s2), m_s3(s3), m_size(shape_size())
	{
		set_picth();
	}

	/************************** constructors *****************************/
	/* copy constructor */
	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	pVctr<T, Dev, ST>::pVctr(const pVctr<T, Dev, ST>& pvctr)
	{
		*this = pvctr;
	}

	/* Move constructor */
	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	pVctr<T, Dev, ST>::pVctr(pVctr<T, Dev, ST>&& pvctr)
	{
		*this = std::move(pvctr);
	}

	/* Converting constructor */
	template <class T, eDev Dev, class ST>
	template <class STU>
	CPU_EXEC
	pVctr<T, Dev, ST>::pVctr(const pVctr<T, Dev, STU>& pvctr)
	{
		*this = pvctr;
	}

	/* constructor from Vctr */
	template <class T, eDev Dev, class ST>
	CPU_EXEC
	pVctr<T, Dev, ST>:: pVctr(const Vctr<T, Dev>& vctr)
	{
		*this = vctr;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	pVctr<T, Dev, ST>& pVctr<T, Dev, ST>::operator=(const pVctr<T, Dev, ST>& pvctr)
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
	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	pVctr<T, Dev, ST>& pVctr<T, Dev, ST>::operator=(pVctr<T, Dev, ST>&& pvctr)
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

	/* Converting assignment operator */
	template <class T, eDev Dev, class ST>
	template <class STU>
	CPU_EXEC
	pVctr<T, Dev, ST>& pVctr<T, Dev, ST>::operator=(const pVctr<T, Dev, STU>& pvctr)
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

	/* Assignment operator: Vctr -> pVctr */
	template <class T, eDev Dev, class ST>
	CPU_EXEC
	pVctr<T, Dev, ST>& pVctr<T, Dev, ST>::operator=(const Vctr<T, Dev>& vctr)
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
	template <class T, eDev Dev, class ST>
	pVctr<T, Dev, dt_int32> pVctr<T, Dev, ST>::ptr_32() const
	{
		return pVctr<T, Dev, dt_int32>(*this);
	}

	template <class T, eDev Dev, class ST>
	pVctr<T, Dev, dt_int64> pVctr<T, Dev, ST>::ptr_64() const
	{
		return pVctr<T, Dev, dt_int64>(*this);
	}

	template <class T, eDev Dev, class ST>
	pVctr<T, Dev, ST>::operator pVctr<T, Dev, chg_btw_int32_int64<ST>>() const
	{
		return pVctr<T, Dev, chg_btw_int32_int64<ST>>(*this);
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::s0() const
	{
		return m_s0;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::s1() const
	{
		return m_s1;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::s2() const
	{
		return m_s2;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::s3() const
	{
		return m_s3;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	dt_int32 pVctr<T, Dev, ST>::s0_32() const
	{
		return static_cast<dt_int32>(m_s0);
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	dt_int32 pVctr<T, Dev, ST>::s1_32() const
	{
		return static_cast<dt_int32>(m_s1);
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	dt_int32 pVctr<T, Dev, ST>::s2_32() const
	{
		return static_cast<dt_int32>(m_s2);
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	dt_int32 pVctr<T, Dev, ST>::s3_32() const
	{
		return static_cast<dt_int32>(m_s3);
	}
		
	template <class T, eDev Dev, class ST>	
	CGPU_EXEC
	dt_int64 pVctr<T, Dev, ST>::s0_64() const
	{
		return static_cast<dt_int64>(m_s0);
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	dt_int64 pVctr<T, Dev, ST>::s1_64() const
	{
		return static_cast<dt_int64>(m_s1);
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	dt_int64 pVctr<T, Dev, ST>::s2_64() const
	{
		return static_cast<dt_int64>(m_s2);
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	dt_int64 pVctr<T, Dev, ST>::s3_64() const
	{
		return static_cast<dt_int64>(m_s3);
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::s0h() const
	{
		return m_s0/ST(2);
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::s1h() const
	{
		return m_s1/ST(2);
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::s2h() const
	{
		return m_s2/ST(2);
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::s3h() const
	{
		return m_s3/ST(2);
	}

	template <class T, eDev Dev, class ST>
	dt_shape_st<ST> pVctr<T, Dev, ST>::shape() const
	{
		return {m_s0, m_s1, m_s2, m_s3};
	}

	template <class T, eDev Dev, class ST>
	dt_shape_st<ST> pVctr<T, Dev, ST>::shape_2d_trs() const
	{
		return {m_s1, m_s0, m_s2, m_s3};
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::shape_size() const 
	{
		return m_s0*max(m_s1, ST(1))*max(m_s2, ST(1))*max(m_s3, ST(1));
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::pitch_s1() const
	{
		return m_pitch_s1;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::pitch_s2() const
	{
		return m_pitch_s2;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::pitch_s3() const
	{
		return m_pitch_s3;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::size() const
	{
		return m_size;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	dt_int32 pVctr<T, Dev, ST>::size_32() const
	{
		return static_cast<dt_int32>(m_size);
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	dt_int64 pVctr<T, Dev, ST>::size_64() const
	{
		return static_cast<dt_int64>(m_size);
	}

	template <class T, eDev Dev, class ST>
	iGrid_1d pVctr<T, Dev, ST>::igrid_1d() const
	{
		return {size_32()};
	}

	template <class T, eDev Dev, class ST>
	iGrid_2d pVctr<T, Dev, ST>::igrid_2d() const
	{
		return {s1_32(), s0_32()};
	}

	template <class T, eDev Dev, class ST>
	iGrid_3d pVctr<T, Dev, ST>::igrid_3d() const
	{
		return {s1_32(), s0_32(), s2_32()};
	}

	template <class T, eDev Dev, class ST>
	iGrid_1d_64 pVctr<T, Dev, ST>::igrid_1d_64() const
	{
		return {size_64()};
	}

	template <class T, eDev Dev, class ST>
	iGrid_2d_64 pVctr<T, Dev, ST>::igrid_2d_64() const
	{
		return {s1_64(), s0_64()};
	}

	template <class T, eDev Dev, class ST>
	iGrid_3d_64 pVctr<T, Dev, ST>::igrid_3d_64() const
	{
		return {s1_64(), s0_64(), s2_64()};
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	dt_bool pVctr<T, Dev, ST>::empty() const
	{
		return m_size == 0;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	dt_bool pVctr<T, Dev, ST>::is_1d() const
	{
		return (m_s0 == 1) || (m_s1 == 1);
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::sub_2_ind(const ST& ix_0) const 
	{ 
		return ix_0;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::sub_2_ind(const ST& ix_0, const ST& ix_1) const 
	{ 
		return ix_0 + m_s0*ix_1;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::sub_2_ind(const ST& ix_0, const ST& ix_1, const ST& ix_2) const 
	{ 
		return ix_0 + m_s0*(ix_1 + m_s1*ix_2);
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	ST pVctr<T, Dev, ST>::sub_2_ind(const ST& ix_0, const ST& ix_1, const ST& ix_2, const ST& ix_3) const 
	{ 
		return ix_0 + m_s0*(ix_1 + m_s1*(ix_2 + m_s2*ix_3));
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	T& pVctr<T, Dev, ST>::operator[](const ST& iy)
	{ 
		return m_data[iy];
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	const T& pVctr<T, Dev, ST>::operator[](const ST& iy) const 
	{ 
		return m_data[iy];
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	T& pVctr<T, Dev, ST>::operator()(const ST& iy)
	{ 
		return m_data[iy];
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	const T& pVctr<T, Dev, ST>::operator()(const ST& iy) const 
	{ 
		return m_data[iy];
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	T& pVctr<T, Dev, ST>::operator()(const ST& ix_0, const ST& ix_1)
	{ 
		return m_data[sub_2_ind(ix_0, ix_1)];
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	const T& pVctr<T, Dev, ST>::operator()(const ST& ix_0, const ST& ix_1) const 
	{ 
		return m_data[sub_2_ind(ix_0, ix_1)];
	}
	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	T& pVctr<T, Dev, ST>::operator()(const ST& ix_0, const ST& ix_1, const ST& ix_2)
	{ 
		return m_data[sub_2_ind(ix_0, ix_1, ix_2)];
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	const T& pVctr<T, Dev, ST>::operator()(const ST& ix_0, const ST& ix_1, const ST& ix_2) const 
	{ 
		return m_data[sub_2_ind(ix_0, ix_1, ix_2)];
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	T& pVctr<T, Dev, ST>::operator()(const ST& ix_0, const ST& ix_1, const ST& ix_2, const ST& ix_3)
	{ 
		return m_data[sub_2_ind(ix_0, ix_1, ix_2, ix_3)];
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	const T& pVctr<T, Dev, ST>::operator()(const ST& ix_0, const ST& ix_1, const ST& ix_2, const ST& ix_3) const 
	{ 
		return m_data[sub_2_ind(ix_0, ix_1, ix_2, ix_3)];
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	T* pVctr<T, Dev, ST>::begin()
	{ 
		return m_data;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	const T* pVctr<T, Dev, ST>::begin() const 
	{ 
		return m_data;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	T* pVctr<T, Dev, ST>::end()
	{ 
		return m_data + m_size;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	const T* pVctr<T, Dev, ST>::end() const 
	{ 
		return m_data + m_size;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	T* pVctr<T, Dev, ST>::data()
	{
		return m_data;
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	const T* pVctr<T, Dev, ST>::data() const
	{
		return m_data;
	}

	template <class T, eDev Dev, class ST>
	template <class U>
	U pVctr<T, Dev, ST>::data_cast()
	{
		return reinterpret_cast<U>(m_data);
	}

	template <class T, eDev Dev, class ST>
	template <class U>
	const U pVctr<T, Dev, ST>::data_cast() const
	{
		return reinterpret_cast<U>(m_data);
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	T pVctr<T, Dev, ST>::front() const 
	{ 
		return m_data[0];
	}

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	T pVctr<T, Dev, ST>::back() const 
	{ 
		return m_data[m_size-1];
	}

#ifdef __CUDACC__		
	FCNS_IMP_GPU_GRID_BLK_VCTR(template <class T COMMA eDev Dev COMMA class ST>, pVctr<T COMMA Dev COMMA ST>);
#endif

	template <class T, eDev Dev, class ST>
	CGPU_EXEC
	void pVctr<T, Dev, ST>::set_picth()
	{
		m_pitch_s1 = m_s0*sizeof(T);
		m_pitch_s2 = m_pitch_s1*m_s1;
		m_pitch_s3 = m_pitch_s2*m_s2;
	}
}