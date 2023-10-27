/*
* This file is part of Multem.
* Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
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

#include "grid_1d.h"
#include "igrid_2d.h"
#include "region_2d.h"

/* derived class */
namespace mt
{
	template <class T, class ST>
	using Grid_2d_st = Grid_sxd<T, ST, edim_2>;

	template <class T>
	using Grid_2d = Grid_sxd<T, dt_int32, edim_2>;

	template <class T>
	using Grid_2d_64 = Grid_sxd<T, dt_int64, edim_2>;
}

/* template specialization 2d */
namespace mt
{
	template <class T, class ST>
	class Grid_sxd<T, ST, edim_2>: public iGrid_sxd<ST, edim_2>
	{
	public:
		using value_type = T;
		using size_type = ST;

		T bs_x;				// simulation box size along x direction (Angstroms)
		T bs_y;				// simulation box size along y direction (Angstroms)

		T rx_0;				// reference coordinate system along x
		T ry_0;				// reference coordinate system along y

		dt_bool pbc_x;		// peridic boundary condition along x
		dt_bool pbc_y;		// peridic boundary condition along y

		dt_bool bwl;		// band-width limit

		T sli_thick;		// slice thicknes

		T drx;				// x-sampling in real space
		T dry;				// y-sampling in real space

		T dgx;				// x-sampling in reciprocal space
		T dgy;				// y-sampling in reciprocal space

		/************************************* constructors ************************************/
		Grid_sxd();

		Grid_sxd(const ST& nx, const ST& ny);

		template <class U, class SU>
		Grid_sxd(const U& bs_x, const U& bs_y, const SU& nx, const SU& ny);

		template <class U, class SU>
		Grid_sxd(const U& bs_x, const U& bs_y, const SU& nx, const SU& ny, const U& rx_0, const U& ry_0, 
		dt_bool pbc_x = true, dt_bool pbc_y = true, dt_bool bwl = false, U sli_thick = 0.5);

		/* copy constructor */
		CGPU_EXEC
		Grid_sxd(const Grid_sxd<T, ST, edim_2>& grid);

		/* converting constructor */
		template <class U, class SU>
		CGPU_EXEC
		Grid_sxd(const Grid_sxd<U, SU, edim_2>& grid);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		Grid_sxd<T, ST, edim_2>& operator=(const Grid_sxd<T, ST, edim_2>& grid);

		/* converting assignment operator */
		template <class U, class SU>
		CGPU_EXEC
		Grid_sxd<T, ST, edim_2>& operator=(const Grid_sxd<U, SU, edim_2>& grid);

		template <class U, class SU> 
		CGPU_EXEC
		void assign(const Grid_sxd<U, SU, edim_2>& grid);

		/************************** user define conversion operators ***************************/
		operator iRegion_Rect_xd<edim_2>() const;

		/***************************************************************************************/
		void set_in_data(const ST& nx, const ST& ny);

		template <class U, class SU>
		void set_in_data(const U& bs_x, const U& bs_y, const SU& nx, const SU& ny);

		template <class U, class SU>
		void set_in_data(const U& bs_x, const U& bs_y, const SU& nx, const SU& ny, 
		const U& rx_0, const U& ry_0, dt_bool pbc_x = true, dt_bool pbc_y = true, dt_bool bwl = false, U sli_thick = 0.5);

		void set_dep_var();

		CGPU_EXEC
		void set_r_0(const T& rx_0, const T& ry_0);

		CGPU_EXEC
		void set_r_0(const R_2d<T>& r_0);

		/***************************************************************************************/
		CGPU_EXEC
		T nx_r() const;

		CGPU_EXEC
		T ny_r() const;

		CGPU_EXEC
		T size_r() const;

		T isize_r() const;

		/***************************************************************************************/
		CGPU_EXEC
		T bs_x_h() const;

		CGPU_EXEC
		T bs_y_h() const;

		CGPU_EXEC
		R_2d<T> bs_h() const;

		/***************************************************************************************/
		CGPU_EXEC
		T bs_min() const;

		CGPU_EXEC
		T bs_max() const;

		CGPU_EXEC
		T bs_h_min() const;

		CGPU_EXEC
		T bs_h_max() const;

		CGPU_EXEC
		T rx_c() const;

		CGPU_EXEC
		T ry_c() const;

		CGPU_EXEC
		R_2d<T> rv_c() const;

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

		CGPU_EXEC
		T r_0_min() const;

		CGPU_EXEC
		T dr_min() const;

		CGPU_EXEC
		T dg_min() const;

		/***************************************************************************************/
		/***************************** Fourier space positions *********************************/
		/***************************************************************************************/
		CGPU_EXEC
		T gx(const ST& ix) const;

		CGPU_EXEC
		T gy(const ST& iy) const;

		CGPU_EXEC
		R_2d<T> gv(const ST& ix, const ST& iy) const;

		CGPU_EXEC
		T gx2(const ST& ix) const;

		CGPU_EXEC
		T gy2(const ST& iy) const;

		CGPU_EXEC
		T g2(const ST& ix, const ST& iy) const;

		CGPU_EXEC
		T g(const ST& ix, const ST& iy) const;

		/***************************************************************************************/
		CGPU_EXEC
		T gx(const ST& ix, const T& gx_0) const;

		CGPU_EXEC
		T gy(const ST& iy, const T& gy_0) const;

		CGPU_EXEC
		R_2d<T> gv(const ST& ix, const ST& iy, const T& gx_0, const T& gy_0) const;

		CGPU_EXEC
		R_2d<T> gv(const ST& ix, const ST& iy, const R_2d<T>& g_0) const;
			
		CGPU_EXEC
		T gx2(const ST& ix, const T& gx_0) const;

		CGPU_EXEC
		T gy2(const ST& iy, const T& gy_0) const;

		CGPU_EXEC
		T g2(const ST& ix, const ST& iy, const T& gx_0, const T& gy_0) const;

		CGPU_EXEC
		T g2(const ST& ix, const ST& iy, const R_2d<T>& g0) const;

		CGPU_EXEC
		T g(const ST& ix, const ST& iy, const T& gx_0, const T& gy_0) const;

		CGPU_EXEC
		T g(const ST& ix, const ST& iy, const R_2d<T>& g0) const;

		/************************************* shift *******************************************/
		CGPU_EXEC
		T gx_sft(const ST& ix) const;

		CGPU_EXEC
		T gy_sft(const ST& iy) const;

		CGPU_EXEC
		R_2d<T> gv_sft(const ST& ix, const ST& iy) const;

		CGPU_EXEC
		T gx2_sft(const ST& ix) const;

		CGPU_EXEC
		T gy2_sft(const ST& iy) const;

		CGPU_EXEC
		T g2_sft(const ST& ix, const ST& iy) const;

		CGPU_EXEC
		T g_sft(const ST& ix, const ST& iy) const;

		/***************************************************************************************/
		CGPU_EXEC
		T gx_sft(const ST& ix, const T& gx_0) const;

		CGPU_EXEC
		T gy_sft(const ST& iy, const T& gy_0) const;

		CGPU_EXEC
		R_2d<T> gv_sft(const ST& ix, const ST& iy, const T& gx_0, const T& gy_0) const;

		CGPU_EXEC
		R_2d<T> gv_sft(const ST& ix, const ST& iy, const R_2d<T>& g_0) const;
			
		CGPU_EXEC
		T gx2_sft(const ST& ix, const T& gx_0) const;

		CGPU_EXEC
		T gy2_sft(const ST& iy, const T& gy_0) const;

		CGPU_EXEC
		T g2_sft(const ST& ix, const ST& iy, const T& gx_0, const T& gy_0) const;

		CGPU_EXEC
		T g2_sft(const ST& ix, const ST& iy, const R_2d<T>& g0) const;

		CGPU_EXEC
		T g_sft(const ST& ix, const ST& iy, const T& gx_0, const T& gy_0) const;

		CGPU_EXEC
		T g_sft(const ST& ix, const ST& iy, const R_2d<T>& g0) const;

		/***************************************************************************************/
		/******************************* real space positions **********************************/
		/***************************************************************************************/
		CGPU_EXEC
		T rx(const ST& ix) const;

		CGPU_EXEC
		T ry(const ST& iy) const;

		CGPU_EXEC
		R_2d<T> rv(const ST& ix, const ST& iy) const;

		CGPU_EXEC
		T rx2(const ST& ix) const;

		CGPU_EXEC
		T ry2(const ST& iy) const;

		CGPU_EXEC
		T r2(const ST& ix, const ST& iy) const;

		CGPU_EXEC
		T r(const ST& ix, const ST& iy) const;

		/***************************************************************************************/
		CGPU_EXEC
		T rx(const ST& ix, const T& x0) const;

		CGPU_EXEC
		T ry(const ST& iy, const T& y0) const;

		CGPU_EXEC
		R_2d<T> rv(const ST& ix, const ST& iy, const T& x0, const T& y0) const;

		CGPU_EXEC
		R_2d<T> rv(const ST& ix, const ST& iy, const R_2d<T>& r_0) const;
			
		CGPU_EXEC
		T rx2(const ST& ix, const T& x0) const;

		CGPU_EXEC
		T ry2(const ST& iy, const T& y0) const;

		CGPU_EXEC
		T r2(const ST& ix, const ST& iy, const T& x0, const T& y0) const;

		CGPU_EXEC
		T r2(const ST& ix, const ST& iy, const R_2d<T>& r_0) const;

		CGPU_EXEC
		T r(const ST& ix, const ST& iy, const T& x0, const T& y0) const;

		CGPU_EXEC
		T r(const ST& ix, const ST& iy, const R_2d<T>& r_0) const;

		/************************************* shift *******************************************/
		CGPU_EXEC
		T rx_sft(const ST& ix) const;

		CGPU_EXEC
		T ry_sft(const ST& iy) const;

		CGPU_EXEC
		R_2d<T> rv_sft(const ST& ix, const ST& iy) const;

		CGPU_EXEC
		T rx2_sft(const ST& ix) const;

		CGPU_EXEC
		T ry2_sft(const ST& iy) const;

		CGPU_EXEC
		T r2_sft(const ST& ix, const ST& iy) const;

		CGPU_EXEC
		T r_sft(const ST& ix, const ST& iy) const;

		/***************************************************************************************/
		CGPU_EXEC
		T rx_sft(const ST& ix, const T& x0) const;

		CGPU_EXEC
		T ry_sft(const ST& iy, const T& y0) const;

		CGPU_EXEC
		R_2d<T> rv_sft(const ST& ix, const ST& iy, const T& x0, const T& y0) const;
			
		CGPU_EXEC
		R_2d<T> rv_sft(const ST& ix, const ST& iy, const R_2d<T>& r_0) const;

		CGPU_EXEC
		T rx2_sft(const ST& ix, const T& x0) const;

		CGPU_EXEC
		T ry2_sft(const ST& iy, const T& y0) const;

		CGPU_EXEC
		T r2_sft(const ST& ix, const ST& iy, const T& x0, const T& y0) const;

		CGPU_EXEC
		T r2_sft(const ST& ix, const ST& iy, const R_2d<T>& r_0) const;

		CGPU_EXEC
		T r_sft(const ST& ix, const ST& iy, const T& x0, const T& y0) const;

		CGPU_EXEC
		T r_sft(const ST& ix, const ST& iy, const R_2d<T>& r_0) const;

		/***************************************************************************************/
		/***************************** from position to index **********************************/
		/***************************************************************************************/
		template <class SU>
		void ix_0_ix_n(const T& x, const T& x_max, SU& ix_0, SU& ix_n) const;

		template <class SU>
		void iy_0_iy_n(const T& y, const T& y_max, SU& iy_0, SU& iy_n) const;

		template <class SU>
		void ix_0_ix_e(const T& x, const T& x_max, SU& ix_0, SU& ix_e) const;

		template <class SU>
		void iy_0_iy_e(const T& y, const T& y_max, SU& iy_0, SU& iy_e) const;

		/************ fds = floor/division by pixel size ************/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_fd(const T& x) const;

		// locate y -> iry
		CGPU_EXEC
		ST ry_2_iry_fd(const T& y) const;

		// locate x, y -> ind
		CGPU_EXEC
		ST rv_2_ir_fd(const T& x, const T& y) const;

		// locate r -> ir using dr_min
		CGPU_EXEC
		ST r_2_ir_fd_dr_min(const T& x) const;

		/********* bfds = bound/floor/division by pixel size ********/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_bfd(const T& x) const;

		// locate y -> iry
		CGPU_EXEC
		ST ry_2_iry_bfd(const T& y) const;

		// locate x, y -> ind
		CGPU_EXEC
		ST rv_2_ir_bfd(const T& x, const T& y) const;

		/********* cds = ceil/division by pixel size/shift **********/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_cd(const T& x) const;

		// locate y -> iry
		CGPU_EXEC
		ST ry_2_iry_cd(const T& y) const;

		// locate x, y -> ind
		CGPU_EXEC
		ST rv_2_ir_cd(const T& x, const T& y) const;

		// locate r -> ir using dr_min
		CGPU_EXEC
		ST r_2_ir_cd_dr_min(const T& x) const;

		/****** bcds = bound/ceil/division by pixel size/shift ******/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_bcd(const T& x) const;

		// locate y -> iry
		CGPU_EXEC
		ST ry_2_iry_bcd(const T& y) const;

		// locate x, y -> ind
		CGPU_EXEC
		ST rv_2_ir_bcd(const T& x, const T& y) const;

		/********* fds = floor/division by pixel size/shift *********/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_fds(const T& x) const;

		// locate y -> iry
		CGPU_EXEC
		ST ry_2_iry_fds(const T& y) const;

		// locate x, y -> ind
		CGPU_EXEC
		ST rv_2_ir_fds(const T& x, const T& y) const;

		// locate r -> ir using dr_min
		CGPU_EXEC
		ST r_2_ir_fds_dr_min(const T& x) const;

		/****** bfds = bound/floor/division by pixel size/shift ******/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_bfds(const T& x) const;

		// locate y -> iry
		CGPU_EXEC
		ST ry_2_iry_bfds(const T& y) const;

		// locate x, y -> ind
		CGPU_EXEC
		ST rv_2_ir_bfds(const T& x, const T& y) const;

		/********* cds = ceil/division by pixel size/shift **********/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_cds(const T& x) const;

		// locate y -> iry
		CGPU_EXEC
		ST ry_2_iry_cds(const T& y) const;

		// locate x, y -> ind
		CGPU_EXEC
		ST rv_2_ir_cds(const T& x, const T& y) const;

		// locate r -> ir using dr_min
		CGPU_EXEC
		ST r_2_ir_cds_dr_min(const T& x) const;

		/****** bcds = bound/ceil/division by pixel size/shift *******/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_bcds(const T& x) const;

		// locate y -> iry
		CGPU_EXEC
		ST ry_2_iry_bcds(const T& y) const;

		// locate x, y -> ind
		CGPU_EXEC
		ST rv_2_ir_bcds(const T& x, const T& y) const;

		/************ From position to index by searching ***********/
		// locate x -> irx
		CGPU_EXEC
		ST rx_2_irx_b(const T& x, ST ix_min=0, ST ix_max=0) const;

		// locate y -> iry
		CGPU_EXEC
		ST ry_2_iry_b(const T& y, ST iy_min=0, ST iy_max=0) const;

		/***************************************************************************************/
		/********************************** check bounds ***************************************/
		/***************************************************************************************/
		CGPU_EXEC
		dt_bool chk_bound_x(const T& x) const;

		CGPU_EXEC
		dt_bool chk_bound_y(const T& y) const;

		CGPU_EXEC
		dt_bool chk_bound(const R_2d<T>& r) const;

		CGPU_EXEC
		dt_bool chk_bound_x_eps(const T& x) const;

		CGPU_EXEC
		dt_bool chk_bound_y_eps(const T& y) const;

		CGPU_EXEC
		dt_bool chk_bound_eps(const R_2d<T>& r) const;

		/***************************************************************************************/
		/*********************************** set bounds ****************************************/
		/***************************************************************************************/
		CGPU_EXEC
		T set_bound_x(const T& x) const;

		CGPU_EXEC
		T set_bound_y(const T& y) const;

		CGPU_EXEC
		R_2d<T> set_bound(const R_2d<T>& r) const;

		/***************************************************************************************/
		/************************************ front/back ***************************************/
		/***************************************************************************************/
		CGPU_EXEC
		T gx_front() const;

		CGPU_EXEC
		T gx_back() const;

		CGPU_EXEC
		T gy_front() const;

		CGPU_EXEC
		T gy_back() const;

		/***************************************************************************************/
		CGPU_EXEC
		T rx_front() const;

		CGPU_EXEC
		T rx_back() const;

		CGPU_EXEC
		T ry_front() const;

		CGPU_EXEC
		T ry_back() const;

		/***************************************************************************************/
		/*********************************** factors *******************************************/
		/* calculate fermi low-pass filter alpha parameter */
		CGPU_EXEC
		T fermi_lpf_alpha() const;

		CGPU_EXEC
		T factor_2pi_rx_ctr(const T& x) const;

		CGPU_EXEC
		T factor_2pi_ry_ctr(const T& y) const;

		CGPU_EXEC
		R_2d<T> factor_2pi_rv_ctr(const R_2d<T>& r) const;

		CPU_EXEC
		Vctr<R_2d<T>, edev_cpu> factor_2pi_rv_ctr(const Vctr<R_2d<T>, edev_cpu>& rv) const;

		/***************************************************************************************/
		iRegion_Rect_xd<edim_2> iregion_rect(const R_2d<T>& r, const T& radius) const;

		iRegion_Rect_xd<edim_2> iregion_rect(const R_2d<T>& r, const T& f0, const T& a, const T& b, const T& c);
	};
}

#include "../src/grid_2d.inl"