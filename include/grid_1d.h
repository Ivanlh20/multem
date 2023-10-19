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

#pragma once

#include "const_enum.h"
#include "math_mt.h"
#include "type_traits_gen.h"
#include "fcns_cgpu_gen.h"
#include "const_enum.h"
#include "igrid_1d.h"
#include "region.cuh"
#include "vctr_cpu.h"
#include "vctr_gpu.h"

/* template definition */
namespace mt
{
	template <class T, class ST, eDim Dim> class Grid_sxd;

	template <class T, eDim Dim>
	using Grid_xd = Grid_sxd<T, dt_int32, Dim>;
}

/* derived class */
namespace mt
{
	template <class T, class ST>
	using Grid_1d_st = Grid_sxd<T, ST, edim_1>;

	template <class T>
	using Grid_1d = Grid_sxd<T, dt_int32, edim_1>;

	template <class T>
	using Grid_1d_64 = Grid_sxd<T, dt_int64, edim_1>;
}

/* template specialization 1d */
namespace mt
{
	template <class T, class ST>
	class Grid_sxd<T, ST, edim_1>: public iGrid_sxd<ST, edim_1>
	{
	public:
		using value_type = T;
		using size_type = ST;

		T bs_x;				// simulation box size along x direction (Angstroms)

		T rx_0;				// reference coordinate system along x

		dt_bool pbc_x;		// peridic boundary condition along x

		dt_bool bwl;		// band-width limit

		T sli_thick;		// slice thicknes

		T drx;				// x-sampling in real space

		T dgx;				// x-sampling in reciprocal space

		/************************************* constructors ************************************/
		CGPU_EXEC
		Grid_sxd();

		Grid_sxd(const ST& nx);

		template <class U, class SU>
		Grid_sxd(const U& bs_x, const SU& nx);

		template <class U, class SU>
		Grid_sxd(const U& bs_x, const SU& nx, const U& rx_0, 
		dt_bool pbc_x = true, dt_bool bwl = false, U sli_thick = 0.5);

		/* copy constructor */
		CGPU_EXEC
		Grid_sxd(const Grid_sxd<T, ST, edim_1>& grid);

		/* converting constructor */
		template <class U, class SU>
		CGPU_EXEC
		Grid_sxd(const Grid_sxd<U, SU, edim_1>& grid);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		Grid_sxd<T, ST, edim_1>& operator=(const Grid_sxd<T, ST, edim_1>& grid);

		/* converting assignment operator */
		template <class U, class SU>
		CGPU_EXEC
		Grid_sxd<T, ST, edim_1>& operator=(const Grid_sxd<U, SU, edim_1>& grid);

		template <class U, class SU> 
		CGPU_EXEC
		void assign(const Grid_sxd<U, SU, edim_1>& grid);

		/************************ user define conversion operators ****************************/
		operator iRegion_Rect_xd<edim_1>() const;

		/***************************************************************************************/
		void set_in_data(const ST& nx);

		template <class U, class SU>
		void set_in_data(const U& bs_x, const SU& nx);

		template <class U, class SU>
		void set_in_data(const U& bs_x, const SU& nx, 
		const U& rx_0, dt_bool pbc_x = true, dt_bool bwl = false, U sli_thick = 0.5);

		void set_dep_var();

		CGPU_EXEC
		void set_r_0(const T& rx_0);

		/***************************************************************************************/
		CGPU_EXEC
		T nx_r() const;

		CGPU_EXEC
		T size_r() const;

		T isize_r() const;

		/***************************************************************************************/
		CGPU_EXEC
		T bs_x_h() const;

		CGPU_EXEC
		T bs_h() const;

		/***************************************************************************************/
		CGPU_EXEC
		T rx_c() const;

		CGPU_EXEC
		T rv_c() const;

		/***************************************************************************************/
		// maximum frequency
		CGPU_EXEC
		T g_max() const;

		// maximum square frequency
		CGPU_EXEC
		T g2_max() const;

		// maximum allowed frequency
		CGPU_EXEC
		T gl_max() const;

		// maximum square allowed frequency
		CGPU_EXEC
		T gl2_max() const;

		/***************************************************************************************/
		/****************************** Fourier space positions ********************************/
		/***************************************************************************************/
		CGPU_EXEC
		T gx(const ST& ix) const;

		CGPU_EXEC
		T gx2(const ST& ix) const;

		CGPU_EXEC
		T g2(const ST& ix) const;

		CGPU_EXEC
		T g(const ST& ix) const;

		/***************************************************************************************/
		CGPU_EXEC
		T gx(const ST& ix, const T& gx_0) const;

		CGPU_EXEC
		T gx2(const ST& ix, const T& gx_0) const;

		CGPU_EXEC
		T g2(const ST& ix, const T& gx_0) const;

		CGPU_EXEC
		T g(const ST& ix, const T& gx_0) const;

		/************************************* shift *******************************************/
		CGPU_EXEC
		T gx_sft(const ST& ix) const;

		CGPU_EXEC
		T gx2_sft(const ST& ix) const;

		CGPU_EXEC
		T g2_sft(const ST& ix) const;

		CGPU_EXEC
		T g_sft(const ST& ix) const;

		/***************************************************************************************/
		CGPU_EXEC
		T gx_sft(const ST& ix, const T& gx_0) const;

		CGPU_EXEC
		T gx2_sft(const ST& ix, const T& gx_0) const;

		CGPU_EXEC
		T g2_sft(const ST& ix, const T& gx_0) const;

		CGPU_EXEC
		T g_sft(const ST& ix, const T& gx_0) const;

		/***************************************************************************************/
		/******************************** real space positions *********************************/
		/***************************************************************************************/
		CGPU_EXEC
		T rx(const ST& ix) const;

		CGPU_EXEC
		T rx2(const ST& ix) const;

		CGPU_EXEC
		T r2(const ST& ix) const;

		CGPU_EXEC
		T r(const ST& ix) const;

		/***************************************************************************************/
		CGPU_EXEC
		T rx(const ST& ix, const T& x0) const;

		CGPU_EXEC
		T rx2(const ST& ix, const T& x0) const;

		CGPU_EXEC
		T r2(const ST& ix, const T& x0) const;

		CGPU_EXEC
		T r(const ST& ix, const T& x0) const;

		/***************************************************************************************/
		/***************************************************************************************/
		CGPU_EXEC
		T rx_sft(const ST& ix) const;

		CGPU_EXEC
		T rx2_sft(const ST& ix) const;

		CGPU_EXEC
		T r2_sft(const ST& ix) const;

		CGPU_EXEC
		T r_sft(const ST& ix) const;

		/***************************************************************************************/
		CGPU_EXEC
		T rx_sft(const ST& ix, const T& x0) const;

		CGPU_EXEC
		T rx2_sft(const ST& ix, const T& x0) const;

		CGPU_EXEC
		T r2_sft(const ST& ix, const T& x0) const;

		CGPU_EXEC
		T r_sft(const ST& ix, const T& x0) const;

		/***************************************************************************************/
		/******************************* from position to index ********************************/
		/***************************************************************************************/
		template <class SU>
		void ix_0_ix_n(const T& x, const T& x_max, SU& ix_0, SU& ix_n) const;
			
		template <class SU>
		void ix_0_ix_e(const T& x, const T& x_max, SU& ix_0, SU& ix_e) const;

		/************ fds = floor/division by pixel size ************/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_fd(const T& x) const;

		/********* bfds = bound/floor/division by pixel size ********/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_bfd(const T& x) const;

		/********* cds = ceil/division by pixel size/shift **********/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_cd(const T& x) const;

		/****** bcds = bound/ceil/division by pixel size/shift ******/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_bcd(const T& x) const;

		/********* fds = floor/division by pixel size/shift *********/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_fds(const T& x) const;

		/****** bfds = bound/floor/division by pixel size/shift ******/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_bfds(const T& x) const;

		/********* cds = ceil/division by pixel size/shift **********/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_cds(const T& x) const;

		/****** bcds = bound/ceil/division by pixel size/shift *******/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_bcds(const T& x) const;

		/************ From position to index by searching ***********/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_b(const T& x, ST ix_min=0, ST ix_max=0) const;

		/***************************************************************************************/
		/********************************** check bounds ***************************************/
		/***************************************************************************************/
		CGPU_EXEC
		dt_bool chk_bound_x(const T& x) const;

		CGPU_EXEC
		dt_bool chk_bound_x_eps(const T& x) const;

		/***************************************************************************************/
		/************************************ set bounds ***************************************/
		/***************************************************************************************/
		CGPU_EXEC
		T set_bound_x(const T& x) const;

		/***************************************************************************************/
		/*********************************** front/back ****************************************/
		/***************************************************************************************/
		CGPU_EXEC
		T gx_front() const;

		CGPU_EXEC
		T gx_back() const;

		/***************************************************************************************/
		CGPU_EXEC
		T rx_front() const;

		CGPU_EXEC
		T rx_back() const;

		/***************************************************************************************/
		/************************************** factors ****************************************/
		/***************************************************************************************/
		/* calculate fermi low-pass filter alpha parameter */
		CGPU_EXEC
		T fermi_lpf_alpha() const;

		CGPU_EXEC
		T factor_2pi_rx_ctr(const T& x) const;

		CPU_EXEC
		Vctr<T, edev_cpu> factor_2pi_rv_ctr(const Vctr<T, edev_cpu>& rv) const;

		/***************************************************************************************/
		iRegion_Rect_xd<edim_1> iregion_rect(const T& r, const T& radius) const;
	};
}


#include "../src/grid_1d.inl"