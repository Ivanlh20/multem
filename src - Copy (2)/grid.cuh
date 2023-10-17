/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef GRID_H
	#define GRID_H

	#ifdef _MSC_VER
		#pragma once
	#endif

	#include "const_enum.cuh"
	#include "math.cuh"
	#include "type_traits_gen.cuh"
	#include "cgpu_fcns_gen.cuh"
	#include "igrid.cuh"
	#include "r_2d.cuh"
	#include "r_3d.cuh"
	#include "region.cuh"
	#include "cgpu_vctr.cuh"

	/***************************************************************************************/
	/*************************************** grid ******************************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T, class ST, eDim Dim> class Grid_sxd;

		template<class T, eDim Dim>
		using Grid_xd = Grid_sxd<T, dt_int32, Dim>;
		
		/* 1d */
		template<class T, class ST>
		using Grid_1d_st = Grid_sxd<T, ST, edim_1>;

		template<class T>
		using Grid_1d = Grid_sxd<T, dt_int32, edim_1>;

		template<class T>
		using Grid_1d_64 = Grid_sxd<T, dt_int64, edim_1>;

		/* 2d */
		template<class T, class ST>
		using Grid_2d_st = Grid_sxd<T, ST, edim_2>;

		template<class T>
		using Grid_2d = Grid_sxd<T, dt_int32, edim_2>;

		template<class T>
		using Grid_2d_64 = Grid_sxd<T, dt_int64, edim_2>;

		/* 3d */	
		template<class T, class ST>
		using Grid_3d_st = Grid_sxd<T, ST, edim_3>;

		template<class T>
		using Grid_3d = Grid_sxd<T, dt_int32, edim_3>;

		template<class T>
		using Grid_3d_64 = Grid_sxd<T, dt_int64, edim_3>;
	}

	/*template specialization 1d */
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

			T sli_thick;			// slice thicknes

			T drx;				// x-sampling in real space

			T dgx;				// x-sampling in reciprocal space

			/************************************* constructors ************************************/
			CGPU_EXEC
			Grid_sxd(): iGrid_sxd<ST, edim_1>(), bs_x(0), rx_0(0), pbc_x(true), 
				bwl(false), sli_thick(0), drx(0), dgx(0) {}

			Grid_sxd(const ST& nx)
			{
				set_in_data(nx);
			}

			template <class U, class SU>
			Grid_sxd(const U& bs_x, const SU& nx)
			{
				set_in_data(bs_x, nx);
			}

			template <class U, class SU>
			Grid_sxd(const U& bs_x, const SU& nx, const U& rx_0, 
			dt_bool pbc_x = true, dt_bool bwl = false, U sli_thick = 0.5)
			{
				set_in_data(bs_x, nx, rx_0, pbc_x, bwl, sli_thick);
			}

			/* copy constructor */
			CGPU_EXEC
			Grid_sxd(const Grid_sxd<T, ST, edim_1>& grid)
			{
				*this = grid;
			}

			/* converting constructor */
			template <class U, class SU>
			CGPU_EXEC
			Grid_sxd(const Grid_sxd<U, SU, edim_1>& grid)
			{
				*this = grid;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			Grid_sxd<T, ST, edim_1>& operator=(const Grid_sxd<T, ST, edim_1>& grid)
			{
				if (this != &grid)
				{
					iGrid_sxd<ST, edim_1>::operator=(grid);

					bs_x = grid.bs_x;
					rx_0 = grid.rx_0;
					pbc_x = grid.pbc_x;
					bwl = grid.bwl;
					sli_thick = grid.sli_thick;

					drx = grid.drx;
					dgx = grid.dgx;
				}

				return *this;
			}

			/* converting assignment operator */
			template <class U, class SU>
			CGPU_EXEC
			Grid_sxd<T, ST, edim_1>& operator=(const Grid_sxd<U, SU, edim_1>& grid)
			{
				iGrid_sxd<ST, edim_1>::operator=(grid);

				bs_x = T(grid.bs_x);
				rx_0 = T(grid.rx_0);
				pbc_x = grid.pbc_x;
				bwl = grid.bwl;
				sli_thick = T(grid.sli_thick);

				drx = T(grid.drx);
				dgx = T(grid.dgx);

				return *this;
			}

			template <class U, class SU> 
			CGPU_EXEC
			void assign(const Grid_sxd<U, SU, edim_1>& grid)
			{
				*this = grid;
			}

			/************************ user define conversion operators ****************************/
			operator iRegion_Rect_xd<edim_1>() const
			{
				return {ST(0), this->nx};
			}

			/***************************************************************************************/
			void set_in_data(const ST& nx)
			{
				set_in_data(T(nx), nx, T(0));
			}

			template <class U, class SU>
			void set_in_data(const U& bs_x, const SU& nx)
			{
				set_in_data(T(bs_x), ST(nx), T(0));
			}

			template <class U, class SU>
			void set_in_data(const U& bs_x, const SU& nx, 
			const U& rx_0, dt_bool pbc_x = true, dt_bool bwl = false, U sli_thick = 0.5)
			{
				this->set_size(nx);

				this->bs_x = T(bs_x);
				this->rx_0 = T(rx_0);
				this->pbc_x = pbc_x;
				this->bwl = bwl;
				this->sli_thick = T(sli_thick);

				set_dep_var();
			}

			void set_dep_var()
			{
				drx = mt::fcn_div(bs_x, this->nx);
				dgx = mt::fcn_div(T(1), bs_x);
			}

			CGPU_EXEC
			void set_r_0(const T& rx_0) 
			{	
				this->rx_0 = rx_0;
			}

			/***************************************************************************************/
			CGPU_EXEC
			T nx_r() const
			{ 
				return T(this->nx);
			}

			CGPU_EXEC
			T size_r() const 
			{ 
				return T(this->size());
			}

			T isize_r() const
			{ 
				return T(1)/size_r();
			}

			/***************************************************************************************/
			CGPU_EXEC
			T bs_x_h() const 
			{ 
				return T(0.5)*bs_x;
			}

			CGPU_EXEC
			T bs_h() const 
			{ 
				return bs_x_h();
			}

			/***************************************************************************************/
			CGPU_EXEC
			T rx_c() const 
			{ 
				return rx_0 + bs_x_h();
			}

			CGPU_EXEC
			T rv_c() const 
			{ 
				return rx_c();
			}

			/***************************************************************************************/
			// maximum frequency
			CGPU_EXEC
			T g_max() const 
			{ 
				return gx_back();
			}

			// maximum square frequency
			CGPU_EXEC
			T g2_max() const 
			{ 
				return ::square(g_max());
			}

			// maximum allowed frequency
			CGPU_EXEC
			T gl_max() const
			{
				return g_max()*T(2.0/3.0);
			}

			// maximum square allowed frequency
			CGPU_EXEC
			T gl2_max() const
			{
				return ::square(gl_max());
			}

			/***************************************************************************************/
			/****************************** Fourier space positions ********************************/
			/***************************************************************************************/
			CGPU_EXEC
			T gx(const ST& ix) const 
			{ 
				return T(this->igx(ix))*dgx;
			}

			CGPU_EXEC
			T gx2(const ST& ix) const 
			{ 
				return ::square(gx(ix));
			}

			CGPU_EXEC
			T g2(const ST& ix) const 
			{ 
				return gx2(ix);
			}

			CGPU_EXEC
			T g(const ST& ix) const 
			{ 
				return ::fabs(gx(ix));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T gx(const ST& ix, const T& gx_0) const 
			{ 
				return gx(ix) - gx_0;
			}

			CGPU_EXEC
			T gx2(const ST& ix, const T& gx_0) const 
			{ 
				return ::square(gx(ix, gx_0));
			}

			CGPU_EXEC
			T g2(const ST& ix, const T& gx_0) const 
			{ 
				return gx2(ix, gx_0);
			}

			CGPU_EXEC
			T g(const ST& ix, const T& gx_0) const 
			{ 
				return ::fabs(gx(ix, gx_0));
			}

			/************************************* shift *******************************************/
			CGPU_EXEC
			T gx_sft(const ST& ix) const 
			{ 
				return T(this->igx_sft(ix))*dgx;
			}

			CGPU_EXEC
			T gx2_sft(const ST& ix) const 
			{ 
				return ::square(gx_sft(ix));
			}

			CGPU_EXEC
			T g2_sft(const ST& ix) const 
			{ 
				return gx2_sft(ix);
			}

			CGPU_EXEC
			T g_sft(const ST& ix) const 
			{ 
				return ::fabs(gx_sft(ix));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T gx_sft(const ST& ix, const T& gx_0) const 
			{ 
				return gx_sft(ix) - gx_0;
			}

			CGPU_EXEC
			T gx2_sft(const ST& ix, const T& gx_0) const 
			{ 
				return ::square(gx_sft(ix, gx_0));
			}

			CGPU_EXEC
			T g2_sft(const ST& ix, const T& gx_0) const 
			{ 
				return gx2_sft(ix, gx_0);
			}

			CGPU_EXEC
			T g_sft(const ST& ix, const T& gx_0) const 
			{ 
				return ::fabs(gx_sft(ix, gx_0));
			}

			/***************************************************************************************/
			/******************************** real space positions *********************************/
			/***************************************************************************************/
			CGPU_EXEC
			T rx(const ST& ix) const 
			{ 
				return T(this->irx(ix))*drx + rx_0;
			}

			CGPU_EXEC
			T rx2(const ST& ix) const 
			{ 
				return ::square(rx(ix));
			}

			CGPU_EXEC
			T r2(const ST& ix) const 
			{ 
				return rx2(ix);
			}

			CGPU_EXEC
			T r(const ST& ix) const 
			{ 
				return ::fabs(rx(ix));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T rx(const ST& ix, const T& x0) const 
			{ 
				return rx(ix) - x0;
			}

			CGPU_EXEC
			T rx2(const ST& ix, const T& x0) const 
			{ 
				return ::square(rx(ix, x0));
			}

			CGPU_EXEC
			T r2(const ST& ix, const T& x0) const 
			{ 
				return rx2(ix, x0);
			}

			CGPU_EXEC
			T r(const ST& ix, const T& x0) const 
			{ 
				return ::fabs(rx(ix, x0));
			}

			/***************************************************************************************/
			/***************************************************************************************/
			CGPU_EXEC
			T rx_sft(const ST& ix) const 
			{ 
				return T(this->irx_sft(ix))*drx + rx_0;
			}

			CGPU_EXEC
			T rx2_sft(const ST& ix) const 
			{ 
				return ::square(rx_sft(ix));
			}

			CGPU_EXEC
			T r2_sft(const ST& ix) const 
			{ 
				return rx2_sft(ix);
			}

			CGPU_EXEC
			T r_sft(const ST& ix) const 
			{ 
				return ::fabs(rx_sft(ix));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T rx_sft(const ST& ix, const T& x0) const 
			{ 
				return rx_sft(ix) - x0;
			}

			CGPU_EXEC
			T rx2_sft(const ST& ix, const T& x0) const 
			{ 
				return ::square(rx_sft(ix, x0));
			}

			CGPU_EXEC
			T r2_sft(const ST& ix, const T& x0) const 
			{ 
				return rx2_sft(ix, x0);
			}

			CGPU_EXEC
			T r_sft(const ST& ix, const T& x0) const 
			{ 
				return ::fabs(rx_sft(ix, x0));
			}

			/***************************************************************************************/
			/******************************* from position to index ********************************/
			/***************************************************************************************/
			template <class SU>
			void ix_0_ix_n(const T& x, const T& x_max, SU& ix_0, SU& ix_n) const 
			{
				fcn_get_idx_0_idx_n(x, x_max, drx, pbc_x, this->nx-1, ix_0, ix_n);
			}	
			
			template <class SU>
			void ix_0_ix_e(const T& x, const T& x_max, SU& ix_0, SU& ix_e) const 
			{
				fcn_get_idx_0_idx_n(x, x_max, drx, pbc_x, this->nx-1, ix_0, ix_e);
				ix_e += ix_0;
			}

			/************ fds = floor/division by pixel size ************/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_fd(const T& x) const 
			{ 
				return fcn_cfloor<ST>(x/drx);
			}

			/********* bfds = bound/floor/division by pixel size ********/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_bfd(const T& x) const 
			{ 
				return fcn_set_bound(rx_2_irx_fd(x), ST(0), this->nx-1);
			}

			/********* cds = ceil/division by pixel size/shift **********/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_cd(const T& x) const 
			{ 
				return fcn_cceil<ST>(x/drx);
			}

			/****** bcds = bound/ceil/division by pixel size/shift ******/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_bcd(const T& x) const 
			{ 
				return fcn_set_bound(rx_2_irx_cd(x), ST(0), this->nx-1);
			}

			/********* fds = floor/division by pixel size/shift *********/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_fds(const T& x) const 
			{ 
				return fcn_cfloor<ST>((x - rx_0)/drx);
			}

			/****** bfds = bound/floor/division by pixel size/shift ******/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_bfds(const T& x) const 
			{ 
				return fcn_set_bound(rx_2_irx_fds(x), ST(0), this->nx-1);
			}

			/********* cds = ceil/division by pixel size/shift **********/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_cds(const T& x) const 
			{ 
				return fcn_cceil<ST>((x - rx_0)/drx);
			}

			/****** bcds = bound/ceil/division by pixel size/shift *******/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_bcds(const T& x) const 
			{ 
				return fcn_set_bound(rx_2_irx_cds(x), ST(0), this->nx-1);
			}

			/************ From position to index by searching ***********/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_b(const T& x, ST ix_min=0, ST ix_max=0) const 
			{
				if (ix_min == ix_max)
				{
					ix_min = ST(0);
					ix_max = this->nx-1;
				}

				return fcn_r_2_ir_b_by_fcn(rx, ix_min, ix_max);
			}

			/***************************************************************************************/
			/********************************** check bounds ***************************************/
			/***************************************************************************************/
			CGPU_EXEC
			dt_bool chk_bound_x(const T& x) const 
			{ 
				return fcn_chk_bound(x, rx_front(), rx_back());
			}

			CGPU_EXEC
			dt_bool chk_bound_x_eps(const T& x) const 
			{ 
				return fcn_chk_bound_eps(x, rx_front(), rx_back());
			}

			/***************************************************************************************/
			/************************************ set bounds ***************************************/
			/***************************************************************************************/
			CGPU_EXEC
			T set_bound_x(const T& x) const 
			{ 
				return fcn_set_bound(x, rx_front(), rx_back());
			}

			/***************************************************************************************/
			/*********************************** front/back ****************************************/
			/***************************************************************************************/
			CGPU_EXEC
			T gx_front() const 
			{ 
				return gx(ST(0));
			}

			CGPU_EXEC
			T gx_back() const 
			{ 
				return gx(this->nx-1);
			}

			/***************************************************************************************/
			CGPU_EXEC
			T rx_front() const
			{ 
				return rx(ST(0));
			}

			CGPU_EXEC
			T rx_back() const
			{ 
				return rx(this->nx-1);
			}

			/***************************************************************************************/
			/************************************** factors ****************************************/
			/***************************************************************************************/
			// ! calculate fermi low-pass filter alpha parameter
			CGPU_EXEC
			T fermi_lpf_alpha() const
			{ 
				return fcn_fermi_lpf_alpha(gl_max(), T(0.25), T(1e-02));
			}	

			CGPU_EXEC
			T factor_2pi_rx_ctr(const T& x) const 
			{ 
				return fcn_n_2pi_sft(x, bs_x_h());
			}

			CPU_EXEC
			Vctr<T, edev_cpu> factor_2pi_rv_ctr(const Vctr<T, edev_cpu>& rv) const
			{
				Vctr<T, edev_cpu> rv_o(rv.size());

				for(auto ik = 0; ik<rv.size(); ik++)
				{
					rv_o[ik] = factor_2pi_rx_ctr(rv[ik]);
				}

				return rv_o;
			}

			/***************************************************************************************/
			iRegion_Rect_xd<edim_1> iregion_rect(const T& r, const T& radius) const 
			{ 
				return {rx_2_irx_bfds(r - radius), rx_2_irx_bcds(r + radius)};
			}
		};
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

			T sli_thick;			// slice thicknes

			T drx;				// x-sampling in real space
			T dry;				// y-sampling in real space

			T dgx;				// x-sampling in reciprocal space
			T dgy;				// y-sampling in reciprocal space

			/************************************* constructors ************************************/
			Grid_sxd(): iGrid_sxd<ST, edim_2>(), bs_x(0), bs_y(0), 
				rx_0(0), ry_0(0), pbc_x(true), pbc_y(true), bwl(false), sli_thick(0), 
				drx(0), dry(0), dgx(0), dgy(0) {}

			Grid_sxd(const ST& nx, const ST& ny)
			{
				set_in_data(nx, ny);
			}

			template <class U, class SU>
			Grid_sxd(const U& bs_x, const U& bs_y, const SU& nx, const SU& ny)
			{
				set_in_data(bs_x, bs_y, nx, ny);
			}

			template <class U, class SU>
			Grid_sxd(const U& bs_x, const U& bs_y, const SU& nx, const SU& ny, const U& rx_0, const U& ry_0, 
			dt_bool pbc_x = true, dt_bool pbc_y = true, dt_bool bwl = false, U sli_thick = 0.5)
			{
				set_in_data(bs_x, bs_y, nx, ny, rx_0, ry_0, pbc_x, pbc_y, bwl, sli_thick);
			}

			/* copy constructor */
			CGPU_EXEC
			Grid_sxd(const Grid_sxd<T, ST, edim_2>& grid)
			{
				*this = grid;
			}

			/* converting constructor */
			template <class U, class SU>
			CGPU_EXEC
			Grid_sxd(const Grid_sxd<U, SU, edim_2>& grid)
			{
				*this = grid;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			Grid_sxd<T, ST, edim_2>& operator=(const Grid_sxd<T, ST, edim_2>& grid)
			{
				if (this != &grid)
				{
					iGrid_sxd<ST, edim_2>::operator=(grid);

					bs_x = grid.bs_x;
					bs_y = grid.bs_y;
					rx_0 = grid.rx_0;
					ry_0 = grid.ry_0;
					pbc_x = grid.pbc_x;
					pbc_y = grid.pbc_y;
					bwl = grid.bwl;
					sli_thick = grid.sli_thick;

					drx = grid.drx;
					dry = grid.dry;
					dgx = grid.dgx;
					dgy = grid.dgy;
				}

				return *this;
			}

			/* converting assignment operator */
			template <class U, class SU>
			CGPU_EXEC
			Grid_sxd<T, ST, edim_2>& operator=(const Grid_sxd<U, SU, edim_2>& grid)
			{
				iGrid_sxd<ST, edim_2>::operator=(grid);

				bs_x = T(grid.bs_x);
				bs_y = T(grid.bs_y);
				rx_0 = T(grid.rx_0);
				ry_0 = T(grid.ry_0);
				pbc_x = grid.pbc_x;
				pbc_y = grid.pbc_y;
				bwl = grid.bwl;
				sli_thick = T(grid.sli_thick);

				drx = T(grid.drx);
				dry = T(grid.dry);
				dgx = T(grid.dgx);
				dgy = T(grid.dgy);

				return *this;
			}

			template <class U, class SU> 
			CGPU_EXEC
			void assign(const Grid_sxd<U, SU, edim_2>& grid)
			{
				*this = grid;
			}

			/************************** user define conversion operators ***************************/
			operator iRegion_Rect_xd<edim_2>() const
			{
				return {ST(0), this->nx, ST(0), this->ny};
			}

			/***************************************************************************************/
			void set_in_data(const ST& nx, const ST& ny)
			{
				set_in_data(T(nx), T(ny), nx, ny, T(0), T(0));
			}

			template <class U, class SU>
			void set_in_data(const U& bs_x, const U& bs_y, const SU& nx, const SU& ny)
			{
				set_in_data(T(bs_x), T(bs_y), ST(nx), ST(ny), T(0), T(0));
			}

			template <class U, class SU>
			void set_in_data(const U& bs_x, const U& bs_y, const SU& nx, const SU& ny, 
			const U& rx_0, const U& ry_0, dt_bool pbc_x = true, dt_bool pbc_y = true, dt_bool bwl = false, U sli_thick = 0.5)
			{
				this->set_size(nx, ny);

				this->bs_x = T(bs_x);
				this->bs_y = T(bs_y);
				this->rx_0 = T(rx_0);
				this->ry_0 = T(ry_0);
				this->pbc_x = pbc_x;
				this->pbc_y = pbc_y;
				this->bwl = bwl;
				this->sli_thick = T(sli_thick);

				set_dep_var();
			}

			void set_dep_var()
			{
				drx = mt::fcn_div(bs_x, this->nx);
				dry = mt::fcn_div(bs_y, this->ny);
				dgx = mt::fcn_div(T(1), bs_x);
				dgy = mt::fcn_div(T(1), bs_y);
			}

			CGPU_EXEC
			void set_r_0(const T& rx_0, const T& ry_0) 
			{	
				this->rx_0 = rx_0;
				this->ry_0 = ry_0;
			}

			CGPU_EXEC
			void set_r_0(const R_2d<T>& r_0) 
			{	
				set_r_0(r_0.x, r_0.y);
			}

			/***************************************************************************************/
			CGPU_EXEC
			T nx_r() const
			{ 
				return T(this->nx);
			}

			CGPU_EXEC
			T ny_r() const 
			{ 
				return T(this->ny);
			}

			CGPU_EXEC
			T size_r() const 
			{ 
				return T(this->size());
			}

			T isize_r() const
			{ 
				return T(1)/size_r();
			}

			/***************************************************************************************/
			CGPU_EXEC
			T bs_x_h() const 
			{ 
				return T(0.5)*bs_x;
			}

			CGPU_EXEC
			T bs_y_h() const 
			{ 
				return T(0.5)*bs_y;
			}

			CGPU_EXEC
			R_2d<T> bs_h() const 
			{ 
				return {bs_x_h(), bs_y_h()};
			}

			/***************************************************************************************/
			CGPU_EXEC
			T bs_min() const 
			{ 
				return ::fmin(bs_x, bs_y);
			}

			CGPU_EXEC
			T bs_max() const 
			{ 
				return ::fmax(bs_x, bs_y);
			}

			CGPU_EXEC
			T bs_h_min() const 
			{ 
				return T(0.5)*bs_min();
			}

			CGPU_EXEC
			T bs_h_max() const 
			{ 
				return T(0.5)*bs_max();
			}

			CGPU_EXEC
			T rx_c() const 
			{ 
				return rx_0 + bs_x_h();
			}

			CGPU_EXEC
			T ry_c() const 
			{ 
				return ry_0 + bs_y_h();
			}

			CGPU_EXEC
			R_2d<T> rv_c() const 
			{ 
				return {rx_c(), ry_c()};
			}

			/***************************************************************************************/
			// maximum frequency
			CGPU_EXEC
			T g_max() const 
			{ 
				return ::fmin(gx_back(), gy_back());
			}

			// maximum square frequency
			CGPU_EXEC
			T g2_max() const 
			{ 
				return ::square(g_max());
			}

			// maximum allowed frequency
			CGPU_EXEC
			T gl_max() const
			{
				return g_max()*T(2.0/3.0);
			}

			// maximum square allowed frequency
			CGPU_EXEC
			T gl2_max() const
			{
				return ::square(gl_max());
			}

			CGPU_EXEC
			T r_0_min() const 
			{ 
				return ::fmin(rx_0, ry_0);
			}

			CGPU_EXEC
			T dr_min() const 
			{ 
				return ::fmin(drx, dry);
			}

			CGPU_EXEC
			T dg_min() const 
			{ 
				return ::fmin(dgx, dgy);
			}

			/***************************************************************************************/
			/***************************** Fourier space positions *********************************/
			/***************************************************************************************/
			CGPU_EXEC
			T gx(const ST& ix) const 
			{ 
				return T(this->igx(ix))*dgx;
			}

			CGPU_EXEC
			T gy(const ST& iy) const 
			{ 
				return T(this->igy(iy))*dgy;
			}

			CGPU_EXEC
			R_2d<T> gv(const ST& ix, const ST& iy) const 
			{ 
				return {gx(ix), gy(iy)};
			}

			CGPU_EXEC
			T gx2(const ST& ix) const 
			{ 
				return ::square(gx(ix));
			}

			CGPU_EXEC
			T gy2(const ST& iy) const 
			{ 
				return ::square(gy(iy));
			}

			CGPU_EXEC
			T g2(const ST& ix, const ST& iy) const 
			{ 
				return gx2(ix) + gy2(iy);
			}

			CGPU_EXEC
			T g(const ST& ix, const ST& iy) const 
			{ 
				return ::sqrt(g2(ix, iy));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T gx(const ST& ix, const T& gx_0) const 
			{ 
				return gx(ix) - gx_0;
			}

			CGPU_EXEC
			T gy(const ST& iy, const T& gy_0) const 
			{ 
				return gy(iy) - gy_0;
			}

			CGPU_EXEC
			R_2d<T> gv(const ST& ix, const ST& iy, const T& gx_0, const T& gy_0) const 
			{ 
				return {gx(ix, gx_0), gy(iy, gy_0)};
			}

			CGPU_EXEC
			R_2d<T> gv(const ST& ix, const ST& iy, const R_2d<T>& g_0) const 
			{ 
				return gv(ix, iy, g_0.x, g_0.y);
			}
			
			CGPU_EXEC
			T gx2(const ST& ix, const T& gx_0) const 
			{ 
				return ::square(gx(ix, gx_0));
			}

			CGPU_EXEC
			T gy2(const ST& iy, const T& gy_0) const 
			{ 
				return ::square(gy(iy, gy_0));
			}

			CGPU_EXEC
			T g2(const ST& ix, const ST& iy, const T& gx_0, const T& gy_0) const 
			{ 
				return gx2(ix, gx_0) + gy2(iy, gy_0);
			}

			CGPU_EXEC
			T g2(const ST& ix, const ST& iy, const R_2d<T>& g0) const 
			{ 
				return gx2(ix, iy, g0.x, g0.y);
			}

			CGPU_EXEC
			T g(const ST& ix, const ST& iy, const T& gx_0, const T& gy_0) const 
			{ 
				return ::sqrt(g2(ix, iy, gx_0, gy_0));
			}

			CGPU_EXEC
			T g(const ST& ix, const ST& iy, const R_2d<T>& g0) const 
			{ 
				return ::sqrt(g2(ix, iy, g0));
			}

			/************************************* shift *******************************************/
			CGPU_EXEC
			T gx_sft(const ST& ix) const 
			{ 
				return T(this->igx_sft(ix))*dgx;
			}

			CGPU_EXEC
			T gy_sft(const ST& iy) const 
			{ 
				return T(this->igy_sft(iy))*dgy;
			}

			CGPU_EXEC
			R_2d<T> gv_sft(const ST& ix, const ST& iy) const 
			{ 
				return {gx_sft(ix), gy_sft(iy)};
			}

			CGPU_EXEC
			T gx2_sft(const ST& ix) const 
			{ 
				return ::square(gx_sft(ix));
			}

			CGPU_EXEC
			T gy2_sft(const ST& iy) const 
			{ 
				return ::square(gy_sft(iy));
			}

			CGPU_EXEC
			T g2_sft(const ST& ix, const ST& iy) const 
			{ 
				return gx2_sft(ix) + gy2_sft(iy);
			}

			CGPU_EXEC
			T g_sft(const ST& ix, const ST& iy) const 
			{ 
				return ::sqrt(g2_sft(ix, iy));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T gx_sft(const ST& ix, const T& gx_0) const 
			{ 
				return gx_sft(ix) - gx_0;
			}

			CGPU_EXEC
			T gy_sft(const ST& iy, const T& gy_0) const 
			{ 
				return gy_sft(iy) - gy_0;
			}

			CGPU_EXEC
			R_2d<T> gv_sft(const ST& ix, const ST& iy, const T& gx_0, const T& gy_0) const 
			{ 
				return {gx_sft(ix, gx_0), gy_sft(iy, gy_0)};
			}

			CGPU_EXEC
			R_2d<T> gv_sft(const ST& ix, const ST& iy, const R_2d<T>& g_0) const 
			{ 
				return gv_sft(ix, iy, g_0.x, g_0.y);
			}
			
			CGPU_EXEC
			T gx2_sft(const ST& ix, const T& gx_0) const 
			{ 
				return ::square(gx_sft(ix, gx_0));
			}

			CGPU_EXEC
			T gy2_sft(const ST& iy, const T& gy_0) const 
			{ 
				return ::square(gy_sft(iy, gy_0));
			}

			CGPU_EXEC
			T g2_sft(const ST& ix, const ST& iy, const T& gx_0, const T& gy_0) const 
			{ 
				return gx2_sft(ix, gx_0) + gy2_sft(iy, gy_0);
			}

			CGPU_EXEC
			T g2_sft(const ST& ix, const ST& iy, const R_2d<T>& g0) const 
			{ 
				return g2_sft(ix, iy, g0.z, g0.y);
			}

			CGPU_EXEC
			T g_sft(const ST& ix, const ST& iy, const T& gx_0, const T& gy_0) const 
			{ 
				return ::sqrt(g2_sft(ix, iy, gx_0, gy_0));
			}

			CGPU_EXEC
			T g_sft(const ST& ix, const ST& iy, const R_2d<T>& g0) const 
			{ 
				return ::sqrt(g2_sft(ix, iy, g0));
			}

			/***************************************************************************************/
			/******************************* real space positions **********************************/
			/***************************************************************************************/
			CGPU_EXEC
			T rx(const ST& ix) const 
			{ 
				return T(this->irx(ix))*drx + rx_0;
			}

			CGPU_EXEC
			T ry(const ST& iy) const 
			{ 
				return T(this->iry(iy))*dry + ry_0;
			}

			CGPU_EXEC
			R_2d<T> rv(const ST& ix, const ST& iy) const 
			{ 
				return {rx(ix), ry(iy)};
			}

			CGPU_EXEC
			T rx2(const ST& ix) const 
			{ 
				return ::square(rx(ix));
			}

			CGPU_EXEC
			T ry2(const ST& iy) const 
			{ 
				return ::square(ry(iy));
			}

			CGPU_EXEC
			T r2(const ST& ix, const ST& iy) const 
			{ 
				return rx2(ix) + ry2(iy);
			}

			CGPU_EXEC
			T r(const ST& ix, const ST& iy) const 
			{ 
				return ::sqrt(r2(ix, iy));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T rx(const ST& ix, const T& x0) const 
			{ 
				return rx(ix) - x0;
			}

			CGPU_EXEC
			T ry(const ST& iy, const T& y0) const 
			{ 
				return ry(iy) - y0;
			}

			CGPU_EXEC
			R_2d<T> rv(const ST& ix, const ST& iy, const T& x0, const T& y0) const 
			{ 
				return {rx(ix, x0), ry(iy, y0)};
			}

			CGPU_EXEC
			R_2d<T> rv(const ST& ix, const ST& iy, const R_2d<T>& r_0) const 
			{ 
				return rv(ix, iy, r_0.x, r_0.y);
			}
			
			CGPU_EXEC
			T rx2(const ST& ix, const T& x0) const 
			{ 
				return ::square(rx(ix, x0));
			}

			CGPU_EXEC
			T ry2(const ST& iy, const T& y0) const 
			{ 
				return ::square(ry(iy, y0));
			}

			CGPU_EXEC
			T r2(const ST& ix, const ST& iy, const T& x0, const T& y0) const 
			{ 
				return rx2(ix, x0) + ry2(iy, y0);
			}

			CGPU_EXEC
			T r2(const ST& ix, const ST& iy, const R_2d<T>& r_0) const 
			{ 
				return r2(ix, iy, r_0.x, r_0.y);
			}

			CGPU_EXEC
			T r(const ST& ix, const ST& iy, const T& x0, const T& y0) const 
			{ 
				return ::sqrt(r2(ix, iy, x0, y0));
			}

			CGPU_EXEC
			T r(const ST& ix, const ST& iy, const R_2d<T>& r_0) const 
			{ 
				return ::sqrt(r2(ix, iy, r_0));
			}

			/************************************* shift *******************************************/
			CGPU_EXEC
			T rx_sft(const ST& ix) const 
			{ 
				return T(this->irx_sft(ix))*drx + rx_0;
			}

			CGPU_EXEC
			T ry_sft(const ST& iy) const 
			{ 
				return T(this->iry_sft(iy))*dry + ry_0;
			}

			CGPU_EXEC
			R_2d<T> rv_sft(const ST& ix, const ST& iy) const 
			{ 
				return {rx_sft(ix), ry_sft(iy)};
			}

			CGPU_EXEC
			T rx2_sft(const ST& ix) const 
			{ 
				return ::square(rx_sft(ix));
			}

			CGPU_EXEC
			T ry2_sft(const ST& iy) const 
			{ 
				return ::square(ry_sft(iy));
			}

			CGPU_EXEC
			T r2_sft(const ST& ix, const ST& iy) const 
			{ 
				return rx2_sft(ix) + ry2_sft(iy);
			}

			CGPU_EXEC
			T r_sft(const ST& ix, const ST& iy) const 
			{ 
				return ::sqrt(r2_sft(ix, iy));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T rx_sft(const ST& ix, const T& x0) const 
			{ 
				return rx_sft(ix) - x0;
			}

			CGPU_EXEC
			T ry_sft(const ST& iy, const T& y0) const 
			{ 
				return ry_sft(iy) - y0;
			}

			CGPU_EXEC
			R_2d<T> rv_sft(const ST& ix, const ST& iy, const T& x0, const T& y0) const 
			{ 
				return {rx_sft(ix, x0), ry_sft(iy, y0)};
			}
			
			CGPU_EXEC
			R_2d<T> rv_sft(const ST& ix, const ST& iy, const R_2d<T>& r_0) const 
			{ 
				return rv_sft(ix, iy, r_0.x, r_0.y);
			}

			CGPU_EXEC
			T rx2_sft(const ST& ix, const T& x0) const 
			{ 
				return ::square(rx_sft(ix, x0));
			}

			CGPU_EXEC
			T ry2_sft(const ST& iy, const T& y0) const 
			{ 
				return ::square(ry_sft(iy, y0));
			}

			CGPU_EXEC
			T r2_sft(const ST& ix, const ST& iy, const T& x0, const T& y0) const 
			{ 
				return rx2_sft(ix, x0) + ry2_sft(iy, y0);
			}

			CGPU_EXEC
			T r2_sft(const ST& ix, const ST& iy, const R_2d<T>& r_0) const 
			{ 
				return r2_sft(ix, iy, r_0.x, r_0.y);
			}

			CGPU_EXEC
			T r_sft(const ST& ix, const ST& iy, const T& x0, const T& y0) const 
			{ 
				return ::sqrt(r2_sft(ix, iy, x0, y0));
			}

			CGPU_EXEC
			T r_sft(const ST& ix, const ST& iy, const R_2d<T>& r_0) const 
			{ 
				return ::sqrt(r2_sft(ix, iy, r_0));
			}

			/***************************************************************************************/
			/***************************** from position to index **********************************/
			/***************************************************************************************/
			template <class SU>
			void ix_0_ix_n(const T& x, const T& x_max, SU& ix_0, SU& ix_n) const 
			{
				fcn_get_idx_0_idx_n(x, x_max, drx, pbc_x, this->nx-1, ix_0, ix_n);
			}

			template <class SU>
			void iy_0_iy_n(const T& y, const T& y_max, SU& iy_0, SU& iy_n) const 
			{
				fcn_get_idx_0_idx_n(y, y_max, dry, pbc_y, this->ny-1, iy_0, iy_n);
			}

			template <class SU>
			void ix_0_ix_e(const T& x, const T& x_max, SU& ix_0, SU& ix_e) const 
			{
				fcn_get_idx_0_idx_n(x, x_max, drx, pbc_x, this->nx-1, ix_0, ix_e);
				ix_e += ix_0;
			}

			template <class SU>
			void iy_0_iy_e(const T& y, const T& y_max, SU& iy_0, SU& iy_e) const 
			{
				fcn_get_idx_0_idx_n(y, y_max, dry, pbc_y, this->ny-1, iy_0, iy_e);
				iy_e += iy_0;
			}

			/************ fds = floor/division by pixel size ************/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_fd(const T& x) const 
			{ 
				return fcn_cfloor<ST>(x/drx);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_fd(const T& y) const 
			{ 
				return fcn_cfloor<ST>(y/dry);
			}

			// locate x, y -> ind
			CGPU_EXEC
			ST rv_2_ir_fd(const T& x, const T& y) const 
			{ 
				return this->sub_2_ind(rx_2_irx_fd(x), ry_2_iry_fd(y));
			}

			// locate r -> ir using dr_min
			CGPU_EXEC
			ST r_2_ir_fd_dr_min(const T& x) const 
			{ 
				return fcn_cfloor<ST>(x/dr_min());
			}

			/********* bfds = bound/floor/division by pixel size ********/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_bfd(const T& x) const 
			{ 
				return fcn_set_bound(rx_2_irx_fd(x), ST(0), this->nx-1);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_bfd(const T& y) const 
			{ 
				return fcn_set_bound(ry_2_iry_fd(y), ST(0), this->ny-1);
			}

			// locate x, y -> ind
			CGPU_EXEC
			ST rv_2_ir_bfd(const T& x, const T& y) const 
			{ 
				return this->sub_2_ind(rx_2_irx_bfd(x), ry_2_iry_bfd(y));
			}

			/********* cds = ceil/division by pixel size/shift **********/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_cd(const T& x) const 
			{ 
				return fcn_cceil<ST>(x/drx);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_cd(const T& y) const 
			{ 
				return fcn_cceil<ST>(y/dry);
			}

			// locate x, y -> ind
			CGPU_EXEC
			ST rv_2_ir_cd(const T& x, const T& y) const 
			{ 
				return this->sub_2_ind(rx_2_irx_cd(x), ry_2_iry_cd(y));
			}

			// locate r -> ir using dr_min
			CGPU_EXEC
			ST r_2_ir_cd_dr_min(const T& x) const 
			{ 
				return static_cast<ST>(::ceil(x/dr_min()));
			}

			/****** bcds = bound/ceil/division by pixel size/shift ******/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_bcd(const T& x) const 
			{ 
				return fcn_set_bound(rx_2_irx_cd(x), ST(0), this->nx-1);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_bcd(const T& y) const 
			{ 
				return fcn_set_bound(ry_2_iry_cd(y), ST(0), this->ny-1);
			}

			// locate x, y -> ind
			CGPU_EXEC
			ST rv_2_ir_bcd(const T& x, const T& y) const 
			{ 
				return this->sub_2_ind(rx_2_irx_bcd(x), ry_2_iry_bcd(y));
			}

			/********* fds = floor/division by pixel size/shift *********/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_fds(const T& x) const 
			{ 
				return fcn_cfloor<ST>((x - rx_0)/drx);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_fds(const T& y) const 
			{ 
				return fcn_cfloor<ST>((y - ry_0)/dry);
			}

			// locate x, y -> ind
			CGPU_EXEC
			ST rv_2_ir_fds(const T& x, const T& y) const 
			{ 
				return this->sub_2_ind(rx_2_irx_fds(x), ry_2_iry_fds(y));
			}

			// locate r -> ir using dr_min
			CGPU_EXEC
			ST r_2_ir_fds_dr_min(const T& x) const 
			{ 
				return fcn_cfloor<ST>((x - r_0_min())/dr_min());
			}

			/****** bfds = bound/floor/division by pixel size/shift ******/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_bfds(const T& x) const 
			{ 
				return fcn_set_bound(rx_2_irx_fds(x), ST(0), this->nx-1);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_bfds(const T& y) const 
			{ 
				return fcn_set_bound(ry_2_iry_fds(y), ST(0), this->ny-1);
			}

			// locate x, y -> ind
			CGPU_EXEC
			ST rv_2_ir_bfds(const T& x, const T& y) const 
			{ 
				return this->sub_2_ind(rx_2_irx_bfds(x), ry_2_iry_bfds(y));
			}

			/********* cds = ceil/division by pixel size/shift **********/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_cds(const T& x) const 
			{ 
				return fcn_cceil<ST>((x - rx_0)/drx);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_cds(const T& y) const 
			{ 
				return fcn_cceil<ST>((y - ry_0)/dry);
			}

			// locate x, y -> ind
			CGPU_EXEC
			ST rv_2_ir_cds(const T& x, const T& y) const 
			{ 
				return this->sub_2_ind(rx_2_irx_cds(x), ry_2_iry_cds(y));
			}

			// locate r -> ir using dr_min
			CGPU_EXEC
			ST r_2_ir_cds_dr_min(const T& x) const 
			{ 
				return fcn_cceil<ST>((x - r_0_min())/dr_min());
			}

			/****** bcds = bound/ceil/division by pixel size/shift *******/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_bcds(const T& x) const 
			{ 
				return fcn_set_bound(rx_2_irx_cds(x), ST(0), this->nx-1);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_bcds(const T& y) const 
			{ 
				return fcn_set_bound(ry_2_iry_cds(y), ST(0), this->ny-1);
			}

			// locate x, y -> ind
			CGPU_EXEC
			ST rv_2_ir_bcds(const T& x, const T& y) const 
			{ 
				return this->sub_2_ind(rx_2_irx_bcds(x), ry_2_iry_bcds(y));
			}

			/************ From position to index by searching ***********/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_b(const T& x, ST ix_min=0, ST ix_max=0) const 
			{
				if (ix_min == ix_max)
				{
					ix_min = ST(0);
					ix_max = this->nx-1;
				}

				return fcn_r_2_ir_b_by_fcn(rx, ix_min, ix_max);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_b(const T& y, ST iy_min=0, ST iy_max=0) const 
			{
				if (iy_min == iy_max)
				{
					iy_min = ST(0);
					iy_max = this->ny-1;
				}

				return fcn_r_2_ir_b_by_fcn(ry, iy_min, iy_max);
			}

			/***************************************************************************************/
			/********************************** check bounds ***************************************/
			/***************************************************************************************/
			CGPU_EXEC
			dt_bool chk_bound_x(const T& x) const 
			{ 
				return fcn_chk_bound(x, rx_front(), rx_back());
			}

			CGPU_EXEC
			dt_bool chk_bound_y(const T& y) const 
			{ 
				return fcn_chk_bound(y, ry_front(), ry_back());
			}

			CGPU_EXEC
			dt_bool chk_bound(const R_2d<T>& r) const 
			{ 
				return chk_bound_x(r.x) && chk_bound_y(r.y);
			}

			CGPU_EXEC
			dt_bool chk_bound_x_eps(const T& x) const 
			{ 
				return fcn_chk_bound_eps(x, rx_front(), rx_back());
			}

			CGPU_EXEC
			dt_bool chk_bound_y_eps(const T& y) const 
			{ 
				return fcn_chk_bound_eps(y, ry_front(), ry_back());
			}

			CGPU_EXEC
			dt_bool chk_bound_eps(const R_2d<T>& r) const 
			{ 
				return chk_bound_x_eps(r.x) && chk_bound_y_eps(r.y);
			}

			/***************************************************************************************/
			/*********************************** set bounds ****************************************/
			/***************************************************************************************/
			CGPU_EXEC
			T set_bound_x(const T& x) const 
			{ 
				return fcn_set_bound(x, rx_front(), rx_back());
			}

			CGPU_EXEC
			T set_bound_y(const T& y) const 
			{ 
				return fcn_set_bound(y, ry_front(), ry_back());
			}

			CGPU_EXEC
			R_2d<T> set_bound(const R_2d<T>& r) const 
			{ 
				return {set_bound_x(r.x), set_bound_y(r.y)};
			}

			/***************************************************************************************/
			/************************************ front/back ***************************************/
			/***************************************************************************************/
			CGPU_EXEC
			T gx_front() const 
			{ 
				return gx(ST(0));
			}

			CGPU_EXEC
			T gx_back() const 
			{ 
				return gx(this->nx-1);
			}

			CGPU_EXEC
			T gy_front() const 
			{ 
				return gy(ST(0));
			}

			CGPU_EXEC
			T gy_back() const 
			{ 
				return gy(this->ny-1);
			}

			/***************************************************************************************/
			CGPU_EXEC
			T rx_front() const 
			{ 
				return rx(ST(0));
			}

			CGPU_EXEC
			T rx_back() const 
			{ 
				return rx(this->nx-1);
			}

			CGPU_EXEC
			T ry_front() const 
			{ 
				return ry(ST(0));
			}

			CGPU_EXEC
			T ry_back() const 
			{ 
				return ry(this->ny-1);
			}

			/***************************************************************************************/
			/*********************************** factors *******************************************/
			// ! calculate fermi low-pass filter alpha parameter
			CGPU_EXEC
			T fermi_lpf_alpha() const
			{ 
				return fcn_fermi_lpf_alpha(gl_max(), T(0.25), T(1e-02));
			}

			CGPU_EXEC
			T factor_2pi_rx_ctr(const T& x) const 
			{ 
				return fcn_n_2pi_sft(x, bs_x_h());
			}

			CGPU_EXEC
			T factor_2pi_ry_ctr(const T& y) const 
			{ 
				return fcn_n_2pi_sft(y, bs_y_h());
			}

			CGPU_EXEC
			R_2d<T> factor_2pi_rv_ctr(const R_2d<T>& r) const
			{
				return {factor_2pi_rx_ctr(r.x), factor_2pi_ry_ctr(r.y)};
			}

			CPU_EXEC
			Vctr<R_2d<T>, edev_cpu> factor_2pi_rv_ctr(const Vctr<R_2d<T>, edev_cpu>& rv) const
			{
				Vctr<R_2d<T>, edev_cpu> rv_o(rv.size());

				for(auto ik = 0; ik<rv.size(); ik++)
				{
					rv_o[ik] = factor_2pi_rv_ctr(rv[ik]);
				}

				return rv_o;
			}

			/***************************************************************************************/
			iRegion_Rect_xd<edim_2> iregion_rect(const R_2d<T>& r, const T& radius) const 
			{ 
				return {rx_2_irx_bfds(r.x - radius), rx_2_irx_bcds(r.x + radius), ry_2_iry_bfds(r.y - radius), ry_2_iry_bcds(r.y + radius)};
			}

			iRegion_Rect_xd<edim_2> iregion_rect(const R_2d<T>& r, const T& f0, const T& a, const T& b, const T& c)
			{
				const T d = log(f0);
				const T dd = c*c-T(4)*a*b;

				const T radius_x = ::sqrt(T(4)*b*d/dd);
				const T radius_y = ::sqrt(T(4)*a*d/dd);

				return {rx_2_irx_bfds(r.x - radius_x), rx_2_irx_bcds(r.x + radius_x), ry_2_iry_bfds(r.y - radius_y), ry_2_iry_bcds(r.y + radius_y)};
			}
		};
	}

	/* template specialization n3d */
	namespace mt
	{
		template <class T, class ST>
		class Grid_sxd<T, ST, edim_3>: public iGrid_sxd<ST, edim_3>
		{
		public:
			using value_type = T;
			using size_type = ST;

			T bs_x;				// simulation box size along x direction (Angstroms)
			T bs_y;				// simulation box size along y direction (Angstroms)
			T bs_z;				// simulation box size along z direction (Angstroms)

			T rx_0;				// reference coordinate system along x
			T ry_0;				// reference coordinate system along y
			T rz_0;				// reference coordinate system along z

			dt_bool pbc_x;		// peridic boundary condition along x
			dt_bool pbc_y;		// peridic boundary condition along y
			dt_bool pbc_z;		// peridic boundary condition along z

			dt_bool bwl;		// band-width limit

			T sli_thick;			// slice thicknes

			T drx;				// x-sampling in real space
			T dry;				// y-sampling in real space
			T drz;				// z-sampling in real space

			T dgx;				// x-sampling in reciprocal space
			T dgy;				// y-sampling in reciprocal space
			T dgz;				// z-sampling in reciprocal space

			/************************************* constructors ************************************/
			Grid_sxd(): iGrid_sxd<ST, edim_3>(), bs_x(0), bs_y(0), bs_z(0), 
				rx_0(0), ry_0(0), rz_0(0), pbc_x(true), pbc_y(true), pbc_z(true), bwl(false), sli_thick(0), 
				drx(0), dry(0), drz(0), dgx(0), dgy(0), dgz(0){}

			Grid_sxd(const ST& nx, const ST& ny, const ST& nz)
			{
				set_in_data(nx, ny, nz);
			}

			template <class U, class SU>
			Grid_sxd(const U& bs_x, const U& bs_y, const U& bs_z, const SU& nx, const SU& ny, const ST& nz)
			{
				set_in_data(bs_x, bs_y, bs_z, nx, ny, nz);
			}

			template <class U, class SU>
			Grid_sxd(const U& bs_x, const U& bs_y, const U& bs_z, const SU& nx, const SU& ny, const SU& nz, 
			const U& rx_0, const U& ry_0, const U& rz_0, dt_bool pbc_x = true, dt_bool pbc_y = true, dt_bool pbc_z = true, dt_bool bwl = false, U sli_thick = 0.5)
			{
				set_in_data(bs_x, bs_y, bs_z, nx, ny, nz, rx_0, ry_0, rz_0, pbc_x, pbc_y, pbc_z, bwl, sli_thick);
			}

			/* copy constructor */
			CGPU_EXEC
			Grid_sxd(const Grid_sxd<T, ST, edim_3>& grid)
			{
				*this = grid;
			}

			/* converting constructor */
			template <class U, class SU>
			CGPU_EXEC
			Grid_sxd(const Grid_sxd<U, SU, edim_3>& grid)
			{
				*this = grid;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			Grid_sxd<T, ST, edim_3>& operator=(const Grid_sxd<T, ST, edim_3>& grid)
			{
				if (this != &grid)
				{
					iGrid_sxd<ST, edim_3>::operator=(grid);

					bs_x = grid.bs_x;
					bs_y = grid.bs_y;
					bs_z = grid.bs_z;
					rx_0 = grid.rx_0;
					ry_0 = grid.ry_0;
					rz_0 = grid.rz_0;
					pbc_x = grid.pbc_x;
					pbc_y = grid.pbc_y;
					pbc_z = grid.pbc_z;
					bwl = grid.bwl;
					sli_thick = grid.sli_thick;

					drx = grid.drx;
					dry = grid.dry;
					drz = grid.drz;
					dgx = grid.dgx;
					dgy = grid.dgy;
					dgz = grid.dgz;
				}

				return *this;
			}

			/* converting assignment operator */
			template <class U, class SU>
			CGPU_EXEC
			Grid_sxd<T, ST, edim_3>& operator=(const Grid_sxd<U, SU, edim_3>& grid)
			{
				iGrid_sxd<ST, edim_3>::operator=(grid);

				bs_x = T(grid.bs_x);
				bs_y = T(grid.bs_y);
				bs_z = T(grid.bs_z);
				rx_0 = T(grid.rx_0);
				ry_0 = T(grid.ry_0);
				rz_0 = T(grid.rz_0);
				pbc_x = grid.pbc_x;
				pbc_y = grid.pbc_y;
				pbc_z = grid.pbc_z;
				bwl = grid.bwl;
				sli_thick = T(grid.sli_thick);

				drx = T(grid.drx);
				dry = T(grid.dry);
				drz = T(grid.drz);
				dgx = T(grid.dgx);
				dgy = T(grid.dgy);
				dgz = T(grid.dgz);

				return *this;
			}

			template <class U, class SU> 
			CGPU_EXEC
			void assign(const Grid_sxd<U, SU, edim_3>& grid)
			{
				*this = grid;
			}

			/**************** user define conversion operators *******************/
			operator iRegion_Rect_xd<edim_3>() const
			{
				return {ST(0), this->nx, ST(0), this->ny, ST(0), this->nz};
			}

			/***************************************************************************************/
			void set_in_data(const ST& nx, const ST& ny, const ST& nz)
			{
				set_in_data(T(nx), T(ny), T(nz), nx, ny, nz, T(0), T(0), T(0));
			}

			template <class U, class SU>
			void set_in_data(const U& bs_x, const U& bs_y, const U& bs_z, const SU& nx, const SU& ny, const SU& nz)
			{
				set_in_data(T(bs_x), T(bs_y), T(bs_z), ST(nx), ST(ny), ST(nz), T(0), T(0), T(0));
			}

			template <class U, class SU>
			void set_in_data(const U& bs_x, const U& bs_y, const U& bs_z, const SU& nx, const SU& ny, const SU& nz, 
			const U& rx_0, const U& ry_0, const U& rz_0, dt_bool pbc_x = true, dt_bool pbc_y = true, dt_bool pbc_z = true, dt_bool bwl = false, U sli_thick = 0.5)
			{
				this->set_size(nx, ny, nz);

				this->bs_x = T(bs_x);
				this->bs_y = T(bs_y);
				this->bs_z = T(bs_z);
				this->rx_0 = T(rx_0);
				this->ry_0 = T(ry_0);
				this->rz_0 = T(rz_0);
				this->pbc_x = pbc_x;
				this->pbc_y = pbc_y;
				this->pbc_z = pbc_z;
				this->bwl = bwl;
				this->sli_thick = T(sli_thick);

				set_dep_var();
			}

			void set_dep_var()
			{
				drx = mt::fcn_div(bs_x, this->nx);
				dry = mt::fcn_div(bs_y, this->ny);
				drz = mt::fcn_div(bs_z, this->nz);
				dgx = mt::fcn_div(T(1), bs_x);
				dgy = mt::fcn_div(T(1), bs_y);
				dgz = mt::fcn_div(T(1), bs_z);
			}

			CGPU_EXEC
			void set_r_0(const T& rx_0, const T& ry_0, const T& rz_0) 
			{	
				this->rx_0 = rx_0;
				this->ry_0 = ry_0;
				this->rz_0 = rz_0;
			}

			CGPU_EXEC
			void set_r_0(const R_3d<T>& r_0) 
			{	
				set_r_0(r_0.x, r_0.y, r_0.z);
			}

			/***************************************************************************************/
			CGPU_EXEC
			T nx_r() const
			{ 
				return T(this->nx);
			}

			CGPU_EXEC
			T ny_r() const 
			{ 
				return T(this->ny);
			}

			CGPU_EXEC
			T nz_r() const 
			{ 
				return T(this->nz);
			}

			CGPU_EXEC
			T size_r() const 
			{ 
				return T(this->size());
			}

			T isize_r() const
			{ 
				return T(1)/size_r();
			}

			/***************************************************************************************/
			CGPU_EXEC
			T bs_x_h() const 
			{ 
				return T(0.5)*bs_x;
			}

			CGPU_EXEC
			T bs_y_h() const 
			{ 
				return T(0.5)*bs_y;
			}

			CGPU_EXEC
			T bs_z_h() const 
			{ 
				return T(0.5)*bs_z;
			}

			CGPU_EXEC
			R_3d<T> bs_h() const 
			{ 
				return {bs_x_h(), bs_y_h(), bs_z_h()};
			}

			/***************************************************************************************/
			CGPU_EXEC
			T bs_min() const 
			{ 
				return ::fmin(bs_x, ::fmin(bs_y, bs_z));
			}

			CGPU_EXEC
			T bs_max() const 
			{ 
				return ::fmax(bs_x, ::fmax(bs_y, bs_z));
			}

			CGPU_EXEC
			T bs_h_min() const 
			{ 
				return T(0.5)*bs_min();
			}

			CGPU_EXEC
			T bs_h_max() const 
			{ 
				return T(0.5)*bs_max();
			}

			CGPU_EXEC
			T rx_c() const 
			{ 
				return rx_0 + bs_x_h();
			}

			CGPU_EXEC
			T ry_c() const 
			{ 
				return ry_0 + bs_y_h();
			}

			CGPU_EXEC
			T rz_c() const 
			{ 
				return rz_0 + bs_z_h();
			}

			CGPU_EXEC
			R_3d<T> rv_c() const 
			{ 
				return {rx_c(), ry_c(), rz_c()};
			}

			/***************************************************************************************/
			// maximum frequency
			CGPU_EXEC
			T g_max() const 
			{ 
				return ::fmin(gx_back(), ::fmin(gy_back(), gz_back()));
			}

			// maximum square frequency
			CGPU_EXEC
			T g2_max() const 
			{ 
				return ::square(g_max());
			}

			// maximum allowed frequency
			CGPU_EXEC
			T gl_max() const
			{
				return g_max()*T(2.0/3.0);
			}

			// maximum square allowed frequency
			CGPU_EXEC
			T gl2_max() const
			{
				return ::square(gl_max());
			}

			CGPU_EXEC
			T r_0_min() const 
			{ 
				return ::fmin(rx_0, ::fmin(ry_0, rz_0));
			}

			CGPU_EXEC
			T dr_min() const 
			{ 
				return ::fmin(drx, ::fmin(dry, drz));
			}

			CGPU_EXEC
			T dg_min() const 
			{ 
				return ::fmin(dgx, ::fmin(dgy, dgz));
			}

			/***************************************************************************************/
			/********************** Fourier space positions **********************/
			/***************************************************************************************/
			CGPU_EXEC
			T gx(const ST& ix) const 
			{ 
				return T(this->igx(ix))*dgx;
			}

			CGPU_EXEC
			T gy(const ST& iy) const 
			{ 
				return T(this->igy(iy))*dgy;
			}

			CGPU_EXEC
			T gz(const ST& iz) const 
			{ 
				return T(this->igz(iz))*dgz;
			}

			CGPU_EXEC
			R_3d<T> gv(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return {gx(ix), gy(iy), gz(iz)};
			}

			CGPU_EXEC
			T gx2(const ST& ix) const 
			{ 
				return ::square(gx(ix));
			}

			CGPU_EXEC
			T gy2(const ST& iy) const 
			{ 
				return ::square(gy(iy));
			}

			CGPU_EXEC
			T gz2(const ST& iz) const 
			{ 
				return ::square(gz(iz));
			}

			CGPU_EXEC
			T g2(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return gx2(ix) + gy2(iy) + gz2(iz);
			}

			CGPU_EXEC
			T g(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return ::sqrt(g2(ix, iy, iz));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T gx(const ST& ix, const T& gx_0) const 
			{ 
				return gx(ix) - gx_0;
			}

			CGPU_EXEC
			T gy(const ST& iy, const T& gy_0) const 
			{ 
				return gy(iy) - gy_0;
			}

			CGPU_EXEC
			T gz(const ST& iz, const T& gz_0) const 
			{ 
				return gz(iz) - gz_0;
			}

			CGPU_EXEC
			R_3d<T> gv(const ST& ix, const ST& iy, const ST& iz, const T& gx_0, const T& gy_0, const T& gz_0) const 
			{ 
				return {gx(ix, gx_0), gy(iy, gy_0), gz(iz, gz_0)};
			}

			CGPU_EXEC
			R_3d<T> gv(const ST& ix, const ST& iy, const ST& iz, const R_3d<T>& g_0) const 
			{ 
				return gv(ix, iy, iz, g_0.x, g_0.y, g_0.z);
			}

			CGPU_EXEC
			T gx2(const ST& ix, const T& gx_0) const 
			{ 
				return ::square(gx(ix, gx_0));
			}

			CGPU_EXEC
			T gy2(const ST& iy, const T& gy_0) const 
			{ 
				return ::square(gy(iy, gy_0));
			}

			CGPU_EXEC
			T gz2(const ST& iz, const T& gz_0) const 
			{ 
				return ::square(gz(iz, gz_0));
			}

			CGPU_EXEC
			T g2(const ST& ix, const ST& iy, const ST& iz, const T& gx_0, const T& gy_0, const T& gz_0) const 
			{ 
				return gx2(ix, gx_0) + gy2(iy, gy_0) + gz2(iz, gz_0);
			}

			CGPU_EXEC
			T g2(const ST& ix, const ST& iy, const ST& iz, const R_3d<T>& g0) const 
			{ 
				return g2(ix, iy, iz, g0.x, g0.y, g0.z);
			}

			CGPU_EXEC
			T g(const ST& ix, const ST& iy, const ST& iz, const T& gx_0, const T& gy_0, const T& gz_0) const 
			{ 
				return ::sqrt(g2(ix, iy, iz, gx_0, gy_0, gz_0));
			}

			CGPU_EXEC
			T g(const ST& ix, const ST& iy, const ST& iz, const R_3d<T>& g0) const 
			{ 
				return ::sqrt(g2(ix, iy, iz, g0));
			}

			/************************************* shift *******************************************/
			CGPU_EXEC
			T gx_sft(const ST& ix) const 
			{ 
				return T(this->igx_sft(ix))*dgx;
			}

			CGPU_EXEC
			T gy_sft(const ST& iy) const 
			{ 
				return T(this->igy_sft(iy))*dgy;
			}

			CGPU_EXEC
			T gz_sft(const ST& iz) const 
			{ 
				return T(this->igz_sft(iz))*dgz;
			}

			CGPU_EXEC
			R_3d<T> gv_sft(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return {gx_sft(ix), gy_sft(iy), gz_sft(iz)};
			}

			CGPU_EXEC
			T gx2_sft(const ST& ix) const 
			{ 
				return ::square(gx_sft(ix));
			}

			CGPU_EXEC
			T gy2_sft(const ST& iy) const 
			{ 
				return ::square(gy_sft(iy));
			}

			CGPU_EXEC
			T gz2_sft(const ST& iz) const 
			{ 
				return ::square(gz_sft(iz));
			}

			CGPU_EXEC
			T g2_sft(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return gx2_sft(ix) + gy2_sft(iy) + gz2_sft(iz);
			}

			CGPU_EXEC
			T g_sft(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return ::sqrt(g2_sft(ix, iy, iz));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T gx_sft(const ST& ix, const T& gx_0) const 
			{ 
				return gx_sft(ix) - gx_0;
			}

			CGPU_EXEC
			T gy_sft(const ST& iy, const T& gy_0) const 
			{ 
				return gy_sft(iy) - gy_0;
			}

			CGPU_EXEC
			T gz_sft(const ST& iz, const T& gz_0) const 
			{ 
				return gz_sft(iz) - gz_0;
			}

			CGPU_EXEC
			R_3d<T> gv_sft(const ST& ix, const ST& iy, const ST& iz, const T& gx_0, const T& gy_0, const T& gz_0) const 
			{ 
				return {gx_sft(ix, gx_0), gy_sft(iy, gy_0), gz_sft(iz, gz_0)};
			}

			CGPU_EXEC
			R_3d<T> gv_sft(const ST& ix, const ST& iy, const ST& iz, const R_3d<T>& g_0) const 
			{ 
				return gv_sft(ix, iy, iz, g_0.x, g_0.y, g_0.z);
			}

			CGPU_EXEC
			T gx2_sft(const ST& ix, const T& gx_0) const 
			{ 
				return ::square(gx_sft(ix, gx_0));
			}

			CGPU_EXEC
			T gy2_sft(const ST& iy, const T& gy_0) const 
			{ 
				return ::square(gy_sft(iy, gy_0));
			}

			CGPU_EXEC
			T gz2_sft(const ST& iz, const T& gz_0) const 
			{ 
				return ::square(gz_sft(iz, gz_0));
			}

			CGPU_EXEC
			T g2_sft(const ST& ix, const ST& iy, const ST& iz, const T& gx_0, const T& gy_0, const T& gz_0) const 
			{ 
				return gx2_sft(ix, gx_0) + gy2_sft(iy, gy_0) + gz2_sft(iz, gz_0);
			}

			CGPU_EXEC
			T g2_sft(const ST& ix, const ST& iy, const ST& iz, const R_3d<T>& g0) const 
			{ 
				return g2_sft(ix, iy, iz, g0.x, g0.y, g0.z);
			}

			CGPU_EXEC
			T g_sft(const ST& ix, const ST& iy, const ST& iz, const T& gx_0, const T& gy_0, const T& gz_0) const 
			{ 
				return ::sqrt(g2_sft(ix, iy, iz, gx_0, gy_0, gz_0));
			}

			CGPU_EXEC
			T g_sft(const ST& ix, const ST& iy, const ST& iz, const R_3d<T>& g0) const 
			{ 
				return ::sqrt(g2_sft(ix, iy, iz, g0));
			}

			/***************************************************************************************/
			/************************ real space positions ***********************/
			/***************************************************************************************/
			CGPU_EXEC
			T rx(const ST& ix) const 
			{ 
				return T(this->irx(ix))*drx + rx_0;
			}

			CGPU_EXEC
			T ry(const ST& iy) const 
			{ 
				return T(this->iry(iy))*dry + ry_0;
			}

			CGPU_EXEC
			T rz(const ST& iz) const 
			{ 
				return T(this->irz(iz))*drz + rz_0;
			}

			CGPU_EXEC
			R_3d<T> rv(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return {rx(ix), ry(iy), rz(iz)};
			}

			CGPU_EXEC
			T rx2(const ST& ix) const 
			{ 
				return ::square(rx(ix));
			}

			CGPU_EXEC
			T ry2(const ST& iy) const 
			{ 
				return ::square(ry(iy));
			}

			CGPU_EXEC
			T rz2(const ST& iz) const 
			{ 
				return ::square(rz(iz));
			}

			CGPU_EXEC
			T r2(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return rx2(ix) + ry2(iy) + rz2(iz);
			}

			CGPU_EXEC
			T r(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return ::sqrt(r2(ix, iy, iz));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T rx(const ST& ix, const T& x0) const 
			{ 
				return rx(ix) - x0;
			}

			CGPU_EXEC
			T ry(const ST& iy, const T& y0) const 
			{ 
				return ry(iy) - y0;
			}

			CGPU_EXEC
			T rz(const ST& iz, const T& z0) const 
			{ 
				return rz(iz) - z0;
			}

			CGPU_EXEC
			R_3d<T> rv(const ST& ix, const ST& iy, const ST& iz, const T& x0, const T& y0, const T& z0) const 
			{ 
				return {rx(ix, x0), ry(iy, y0), rz(iz, z0)};
			}

			CGPU_EXEC
			R_3d<T> rv(const ST& ix, const ST& iy, const ST& iz, const R_3d<T>& r_0) const 
			{ 
				return rv(ix, iy, iz, r_0.x, r_0.y, r_0.z);
			}

			CGPU_EXEC
			T rx2(const ST& ix, const T& x0) const 
			{ 
				return ::square(rx(ix, x0));
			}

			CGPU_EXEC
			T ry2(const ST& iy, const T& y0) const 
			{ 
				return ::square(ry(iy, y0));
			}

			CGPU_EXEC
			T rz2(const ST& iz, const T& z0) const 
			{ 
				return ::square(rz(iz, z0));
			}

			CGPU_EXEC
			T r2(const ST& ix, const ST& iy, const ST& iz, const T& x0, const T& y0, const T& z0) const 
			{ 
				return rx2(ix, x0) + ry2(iy, y0) + rz2(iz, z0);
			}

			CGPU_EXEC
			T r2(const ST& ix, const ST& iy, const ST& iz, const R_3d<T>& r_0) const 
			{ 
				return r2(ix, iy, iz, r_0.x, r_0.y, r_0.z);
			}

			CGPU_EXEC
			T r(const ST& ix, const ST& iy, const ST& iz, const T& x0, const T& y0, const T& z0) const 
			{ 
				return ::sqrt(r2(ix, iy, iz, x0, y0, z0));
			}

			CGPU_EXEC
			T r(const ST& ix, const ST& iy, const ST& iz, const R_3d<T>& r_0) const 
			{ 
				return ::sqrt(r2(ix, iy, iz, r_0));
			}

			/************************************* shift *******************************************/
			CGPU_EXEC
			T rx_sft(const ST& ix) const 
			{ 
				return T(this->irx_sft(ix))*drx + rx_0;
			}

			CGPU_EXEC
			T ry_sft(const ST& iy) const 
			{ 
				return T(this->iry_sft(iy))*dry + ry_0;
			}

			CGPU_EXEC
			T rz_sft(const ST& iz) const 
			{ 
				return T(this->irz_sft(iz))*drz + rz_0;
			}

			CGPU_EXEC
			R_3d<T> rv_sft(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return {rx_sft(ix), ry_sft(iy), rz_sft(iz)};
			}

			CGPU_EXEC
			T rx2_sft(const ST& ix) const 
			{ 
				return ::square(rx_sft(ix));
			}

			CGPU_EXEC
			T ry2_sft(const ST& iy) const 
			{ 
				return ::square(ry_sft(iy));
			}

			CGPU_EXEC
			T rz2_sft(const ST& iz) const 
			{ 
				return ::square(rz_sft(iz));
			}

			CGPU_EXEC
			T r2_sft(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return rx2_sft(ix) + ry2_sft(iy) + rz2_sft(iz);
			}

			CGPU_EXEC
			T r_sft(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return ::sqrt(r2_sft(ix, iy, iz));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T rx_sft(const ST& ix, const T& x0) const 
			{ 
				return rx_sft(ix) - x0;
			}

			CGPU_EXEC
			T ry_sft(const ST& iy, const T& y0) const 
			{ 
				return ry_sft(iy) - y0;
			}

			CGPU_EXEC
			T rz_sft(const ST& iz, const T& z0) const 
			{ 
				return rz_sft(iz) - z0;
			}

			CGPU_EXEC
			R_3d<T> rv_sft(const ST& ix, const ST& iy, const ST& iz, const T& x0, const T& y0, const T& z0) const 
			{ 
				return {rx_sft(ix, x0), ry_sft(iy, y0), rz_sft(iz, z0)};
			}

			CGPU_EXEC
			R_3d<T> rv_sft(const ST& ix, const ST& iy, const ST& iz, const R_3d<T>& r_0) const 
			{ 
				return rv_sft(ix, iy, iz, r_0.x, r_0.y, r_0.z);
			}

			CGPU_EXEC
			T rx2_sft(const ST& ix, const T& x0) const 
			{ 
				return ::square(rx_sft(ix, x0));
			}

			CGPU_EXEC
			T ry2_sft(const ST& iy, const T& y0) const 
			{ 
				return ::square(ry_sft(iy, y0));
			}

			CGPU_EXEC
			T rz2_sft(const ST& iz, const T& z0) const 
			{ 
				return ::square(rz_sft(iz, z0));
			}

			CGPU_EXEC
			T r2_sft(const ST& ix, const ST& iy, const ST& iz, const T& x0, const T& y0, const T& z0) const 
			{ 
				return rx2_sft(ix, x0) + ry2_sft(iy, y0) + rz2_sft(iz, z0);
			}

			CGPU_EXEC
			T r2_sft(const ST& ix, const ST& iy, const ST& iz, const R_3d<T>& r_0) const 
			{ 
				return r2_sft(ix, iy, iz, r_0.x, r_0.y, r_0.z);
			}

			CGPU_EXEC
			T r_sft(const ST& ix, const ST& iy, const ST& iz, const T& x0, const T& y0, const T& z0) const 
			{ 
				return ::sqrt(r2_sft(ix, iy, iz, x0, y0, z0));
			}

			CGPU_EXEC
			T r_sft(const ST& ix, const ST& iy, const ST& iz, const R_3d<T>& r_0) const 
			{ 
				return ::sqrt(r2_sft(ix, iy, iz, r_0));
			}

			/***************************************************************************************/
			/******************** from position to index *************************/
			/***************************************************************************************/
			template <class SU>
			void ix_0_ix_n(const T& x, const T& x_max, SU& ix_0, SU& ix_n) const 
			{
				fcn_get_idx_0_idx_n(x, x_max, drx, pbc_x, this->nx-1, ix_0, ix_n);
			}

			template <class SU>
			void iy_0_iy_n(const T& y, const T& y_max, SU& iy_0, SU& iy_n) const 
			{
				fcn_get_idx_0_idx_n(y, y_max, dry, pbc_y, this->ny-1, iy_0, iy_n);
			}

			template <class SU>
			void iz_0_iz_n(const T& z, const T& z_max, SU& iz_0, SU& iz_n) const 
			{
				fcn_get_idx_0_idx_n(z, z_max, drz, pbc_z, this->nz-1, iz_0, iz_n);
			}

			template <class SU>
			void ix_0_ix_e(const T& x, const T& x_max, SU& ix_0, SU& ix_e) const 
			{
				fcn_get_idx_0_idx_n(x, x_max, drx, pbc_x, this->nx-1, ix_0, ix_e);
				ix_e += ix_0;
			}

			template <class SU>
			void iy_0_iy_e(const T& y, const T& y_max, SU& iy_0, SU& iy_e) const 
			{
				fcn_get_idx_0_idx_n(y, y_max, dry, pbc_y, this->ny-1, iy_0, iy_e);
				iy_e += iy_0;
			}

			template <class SU>
			void iz_0_iz_e(const T& z, const T& z_max, SU& iz_0, SU& iz_e) const 
			{
				fcn_get_idx_0_idx_n(z, z_max, drz, pbc_z, this->nz-1, iz_0, iz_e);
				iz_e += iz_0;
			}

			/*************** fds = floor/division by pixel size ******************/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_fd(const T& x) const 
			{ 
				return fcn_cfloor<ST>(x/drx);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_fd(const T& y) const 
			{ 
				return fcn_cfloor<ST>(y/dry);
			}

			// locate z -> irz
			CGPU_EXEC
			ST rz_2_irz_fd(const T& z) const 
			{ 
				return fcn_cfloor<ST>(z/drz);
			}

			// locate x, y, z -> ind
			CGPU_EXEC
			ST rv_2_ir_fd(const T& x, const T& y, const T& z) const 
			{ 
				return this->sub_2_ind(rx_2_irx_fd(x), ry_2_iry_fd(y), rz_2_irz_fd(z));
			}

			// locate r -> ir using dr_min
			CGPU_EXEC
			ST r_2_ir_fd_dr_min(const T& x) const 
			{ 
				return fcn_cfloor<ST>(x/dr_min());
			}

			/********* bfds = bound/floor/division by pixel size ********/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_bfd(const T& x) const 
			{ 
				return fcn_set_bound(rx_2_irx_fd(x), ST(0), this->nx-1);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_bfd(const T& y) const 
			{ 
				return fcn_set_bound(ry_2_iry_fd(y), ST(0), this->ny-1);
			}

			// locate z-> irz
			CGPU_EXEC
			ST rz_2_irz_bfd(const T& z) const 
			{ 
				return fcn_set_bound(rz_2_irz_fd(z), ST(0), this->nz-1);
			}

			// locate x, y, z -> ind
			CGPU_EXEC
			ST rv_2_ir_bfd(const T& x, const T& y, const T& z) const 
			{ 
				return this->sub_2_ind(rx_2_irx_bfd(x), ry_2_iry_bfd(y), rz_2_irz_bfd(z));
			}

			/********* cds = ceil/division by pixel size/shift **********/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_cd(const T& x) const 
			{ 
				return fcn_cceil<ST>(x/drx);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_cd(const T& y) const 
			{ 
				return fcn_cceil<ST>(y/dry);
			}

			// locate z -> irz
			CGPU_EXEC
			ST rz_2_irz_cd(const T& z) const 
			{ 
				return fcn_cceil<ST>(z/drz);
			}

			// locate x, y, z -> ind
			CGPU_EXEC
			ST rv_2_ir_cd(const T& x, const T& y, const T& z) const 
			{ 
				return this->sub_2_ind(rx_2_irx_cd(x), ry_2_iry_cd(y), rz_2_irz_cd(z));
			}

			// locate r -> ir using dr_min
			CGPU_EXEC
			ST r_2_ir_cd_dr_min(const T& x) const 
			{ 
				return static_cast<ST>(::ceil(x/dr_min()));
			}

			/****** bcds = bound/ceil/division by pixel size/shift ******/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_bcd(const T& x) const 
			{ 
				return fcn_set_bound(rx_2_irx_cd(x), ST(0), this->nx-1);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_bcd(const T& y) const 
			{ 
				return fcn_set_bound(ry_2_iry_cd(y), ST(0), this->ny-1);
			}

			// locate z -> iry
			CGPU_EXEC
			ST rz_2_irz_bcd(const T& z) const 
			{ 
				return fcn_set_bound(rz_2_irz_cd(z), ST(0), this->nz-1);
			}

			// locate x, y, z -> ind
			CGPU_EXEC
			ST rv_2_ir_bcd(const T& x, const T& y, const T& z) const 
			{ 
				return this->sub_2_ind(rx_2_irx_bcd(x), ry_2_iry_bcd(y), rz_2_irz_bcd(z));
			}

			/********* fds = floor/division by pixel size/shift *********/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_fds(const T& x) const 
			{ 
				return fcn_cfloor<ST>((x - rx_0)/drx);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_fds(const T& y) const 
			{ 
				return fcn_cfloor<ST>((y - ry_0)/dry);
			}

			// locate z -> irz
			CGPU_EXEC
			ST rz_2_irz_fds(const T& z) const 
			{ 
				return fcn_cfloor<ST>((z - rz_0)/drz);
			}

			// locate x, y, z -> ind
			CGPU_EXEC
			ST rv_2_ir_fds(const T& x, const T& y, const T& z) const 
			{ 
				return this->sub_2_ind(rx_2_irx_fds(x), ry_2_iry_fds(y), rz_2_irz_fds(z));
			}

			// locate r -> ir using dr_min
			CGPU_EXEC
			ST r_2_ir_fds_dr_min(const T& x) const 
			{ 
				return fcn_cfloor<ST>((x - r_0_min())/dr_min());
			}

			/****** bfds = bound/floor/division by pixel size/shift ******/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_bfds(const T& x) const 
			{ 
				return fcn_set_bound(rx_2_irx_fds(x), ST(0), this->nx-1);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_bfds(const T& y) const 
			{ 
				return fcn_set_bound(ry_2_iry_fds(y), ST(0), this->ny-1);
			}

			// locate z -> irz
			CGPU_EXEC
			ST rz_2_irz_bfds(const T& z) const 
			{ 
				return fcn_set_bound(rz_2_irz_fds(z), ST(0), this->nz-1);
			}

			// locate x, y, z -> ind
			CGPU_EXEC
			ST rv_2_ir_bfds(const T& x, const T& y, const T& z) const 
			{ 
				return this->sub_2_ind(rx_2_irx_bfds(x), ry_2_iry_bfds(y), rz_2_irz_bfds(z));
			}

			/********* cds = ceil/division by pixel size/shift **********/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_cds(const T& x) const 
			{ 
				return fcn_cceil<ST>((x - rx_0)/drx);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_cds(const T& y) const 
			{ 
				return fcn_cceil<ST>((y - ry_0)/dry);
			}

			// locate z -> irz
			CGPU_EXEC
			ST rz_2_irz_cds(const T& z) const 
			{ 
				return fcn_cceil<ST>((z - rz_0)/drz);
			}

			// locate x, y, z -> ind
			CGPU_EXEC
			ST rv_2_ir_cds(const T& x, const T& y, const T& z) const 
			{ 
				return this->sub_2_ind(rx_2_irx_cds(x), ry_2_iry_cds(y), rz_2_irz_cds(z));
			}

			// locate r -> ir using dr_min
			CGPU_EXEC
			ST r_2_ir_cds_dr_min(const T& x) const 
			{ 
				return fcn_cceil<ST>((x - r_0_min())/dr_min());
			}

			/****** bcds = bound/ceil/division by pixel size/shift *******/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_bcds(const T& x) const 
			{ 
				return fcn_set_bound(rx_2_irx_cds(x), ST(0), this->nx-1);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_bcds(const T& y) const 
			{ 
				return fcn_set_bound(ry_2_iry_cds(y), ST(0), this->ny-1);
			}

			// locate z -> irz
			CGPU_EXEC
			ST rz_2_irz_bcds(const T& z) const 
			{ 
				return fcn_set_bound(rz_2_irz_cds(z), ST(0), this->nz-1);
			}

			// locate x, y, z -> ind
			CGPU_EXEC
			ST rv_2_ir_bcds(const T& x, const T& y, const T& z) const 
			{ 
				return this->sub_2_ind(rx_2_irx_bcds(x), ry_2_iry_bcds(y), rz_2_irz_bcds(z));
			}

			/************ From position to index by searching ***********/
			// locate x -> irx
			CGPU_EXEC
			ST rx_2_irx_b(const T& x, ST ix_min=0, ST ix_max=0) const 
			{
				if (ix_min == ix_max)
				{
					ix_min = ST(0);
					ix_max = this->nx-1;
				}

				return fcn_r_2_ir_b_by_fcn(rx, ix_min, ix_max);
			}

			// locate y -> iry
			CGPU_EXEC
			ST ry_2_iry_b(const T& y, ST iy_min=0, ST iy_max=0) const 
			{
				if (iy_min == iy_max)
				{
					iy_min = ST(0);
					iy_max = this->ny-1;
				}

				return fcn_r_2_ir_b_by_fcn(ry, iy_min, iy_max);
			}

			// locate z-> irz
			CGPU_EXEC
			ST rz_2_irz_b(const T& z, ST iz_min=0, ST iz_max=0) const 
			{
				if (iz_min == iz_max)
				{
					iz_min = ST(0);
					iz_max = this->nz-1;
				}

				return fcn_r_2_ir_b_by_fcn(rz, iz_min, iz_max);
			}

			/***************************************************************************************/
			/********************************* check bounds ****************************************/
			/***************************************************************************************/
			CGPU_EXEC
			dt_bool chk_bound_x(const T& x) const 
			{ 
				return fcn_chk_bound(x, rx_front(), rx_back());
			}

			CGPU_EXEC
			dt_bool chk_bound_y(const T& y) const 
			{ 
				return fcn_chk_bound(y, ry_front(), ry_back());
			}

			CGPU_EXEC
			dt_bool chk_bound_z(const T& z) const 
			{ 
				return fcn_chk_bound(z, rz_front(), rz_back());
			}

			CGPU_EXEC
			dt_bool chk_bound(const R_3d<T>& r) const 
			{ 
				return chk_bound_x(r.x) && chk_bound_y(r.y) && chk_bound_z(r.z);
			}

			CGPU_EXEC
			dt_bool chk_bound_x_eps(const T& x) const 
			{ 
				return fcn_chk_bound_eps(x, rx_front(), rx_back());
			}

			CGPU_EXEC
			dt_bool chk_bound_y_eps(const T& y) const 
			{ 
				return fcn_chk_bound_eps(y, ry_front(), ry_back());
			}

			CGPU_EXEC
			dt_bool chk_bound_z_eps(const T& z) const 
			{ 
				return fcn_chk_bound_eps(z, rz_front(), rz_back());
			}

			CGPU_EXEC
			dt_bool chk_bound_eps(const R_3d<T>& r) const 
			{ 
				return chk_bound_x_eps(r.x) && chk_bound_y_eps(r.y) && chk_bound_z_eps(r.z);
			}

			/***************************************************************************************/
			/*********************************** set bounds ****************************************/
			/***************************************************************************************/
			CGPU_EXEC
			T set_bound_x(const T& x) const 
			{ 
				return fcn_set_bound(x, rx_front(), rx_back());
			}

			CGPU_EXEC
			T set_bound_y(const T& y) const 
			{ 
				return fcn_set_bound(y, ry_front(), ry_back());
			}

			CGPU_EXEC
			T set_bound_z(const T& z) const 
			{ 
				return fcn_set_bound(z, rz_front(), rz_back());
			}

			CGPU_EXEC
			R_3d<T> set_bound(const R_3d<T>& r) const 
			{ 
				return {set_bound_x(r.x), set_bound_y(r.y), set_bound_z(r.z)};
			}

			/***************************************************************************************/
			/*********************************** front/back ****************************************/
			/***************************************************************************************/
			CGPU_EXEC
			T gx_front() const 
			{ 
				return gx(ST(0));
			}

			CGPU_EXEC
			T gx_back() const 
			{ 
				return gx(this->nx-1);
			}

			CGPU_EXEC
			T gy_front() const 
			{ 
				return gy(ST(0));
			}

			CGPU_EXEC
			T gy_back() const 
			{ 
				return gy(this->ny-1);
			}

			CGPU_EXEC
			T gz_front() const 
			{ 
				return gz(ST(0));
			}

			CGPU_EXEC
			T gz_back() const 
			{ 
				return gz(this->nz-1);
			}

			/***************************************************************************************/
			CGPU_EXEC
			T rx_front() const 
			{ 
				return rx(ST(0));
			}

			CGPU_EXEC
			T rx_back() const 
			{ 
				return rx(this->nx-1);
			}

			CGPU_EXEC
			T ry_front() const 
			{ 
				return ry(ST(0));
			}

			CGPU_EXEC
			T ry_back() const 
			{ 
				return ry(this->ny-1);
			}

			CGPU_EXEC
			T rz_front() const 
			{ 
				return rz(ST(0));
			}

			CGPU_EXEC
			T rz_back() const 
			{ 
				return rz(this->nz-1);
			}

			/***************************************************************************************/
			/************************************* factors *****************************************/
			/***************************************************************************************/
			// ! calculate fermi low-pass filter alpha parameter
			CGPU_EXEC
			T fermi_lpf_alpha() const
			{ 
				return fcn_fermi_lpf_alpha(gl_max(), T(0.25), T(1e-02));
			}

			CGPU_EXEC
			T factor_2pi_rx_ctr(const T& x) const 
			{ 
				return fcn_n_2pi_sft(x, bs_x_h());
			}

			CGPU_EXEC
			T factor_2pi_ry_ctr(const T& y) const 
			{ 
				return fcn_n_2pi_sft(y, bs_y_h());
			}

			CGPU_EXEC
			T factor_2pi_rz_ctr(const T& z) const 
			{ 
				return fcn_n_2pi_sft(z, bs_z_h());
			}

			CGPU_EXEC
			R_3d<T> factor_2pi_rv_ctr(const R_3d<T>& r) const
			{
				return {factor_2pi_rx_ctr(r.x), factor_2pi_ry_ctr(r.y), factor_2pi_rz_ctr(r.z)};
			}

			CPU_EXEC
			Vctr<R_3d<T>, edev_cpu> factor_2pi_rv_ctr(const Vctr<R_3d<T>, edev_cpu>& rv) const
			{
				Vctr<R_3d<T>, edev_cpu> rv_o(rv.size());

				for(auto ik = 0; ik<rv.size(); ik++)
				{
					rv_o[ik] = factor_2pi_rv_ctr(rv[ik]);
				}

				return rv_o;
			}

			/***************************************************************************************/
			iRegion_Rect_xd<edim_3> iregion_rect(const R_3d<T>& r, const T& radius) const 
			{ 
				return {rx_2_irx_bfds(r.x - radius), rx_2_irx_bcds(r.x + radius), ry_2_iry_bfds(r.y - radius), 
					ry_2_iry_bcds(r.y + radius), rz_2_irz_bfds(r.z - radius), rz_2_irz_bcds(r.z + radius)};
			}

			iRegion_Rect_xd<edim_3> iregion_rect(const R_3d<T>& r, const T& f0, const T& a, const T& b, const T& c)
			{
				const T d = log(f0);
				const T dd = c*c-T(4)*a*b;

				const T radius_x = ::sqrt(T(4)*b*d/dd);
				const T radius_y = ::sqrt(T(4)*a*d/dd);
				const T radius_z = ::sqrt(T(4)*a*d/dd);

				return {nrx_2_irx_bfds(r.x - radius_x), rx_2_irx_bcds(r.x + radius_x), ry_2_iry_bfds(r.y - radius_y), 
					ry_2_iry_bcds(r.y + radius_y), rz_2_irz_bfds(r.z - radius_z), rz_2_irz_bcds(r.z + radius_z)};
			}
		};
	}

	/* traits */
	namespace mt
	{
		template <class T>
		struct is_grid_1d: std::integral_constant<dt_bool, std::is_same<T, Grid_1d_st<typename T::value_type, typename T::size_type>>::value> {};	

		template <class T>
		struct is_grid_2d: std::integral_constant<dt_bool, std::is_same<T, Grid_2d_st<typename T::value_type, typename T::size_type>>::value> {};		

		template <class T>
		struct is_grid_3d: std::integral_constant<dt_bool, std::is_same<T, Grid_3d_st<typename T::value_type, typename T::size_type>>::value> {};		
	
		template <class T>
		struct is_grid: std::integral_constant<dt_bool, is_grid_1d<T>::value && is_grid_2d<T>::value && is_grid_3d<T>::value> {};	

		/***************************************************************************************/
		template <class T, class U>
		struct is_grid_1d_and_vctr_cpu: std::integral_constant<dt_bool, is_grid_1d<T>::value && is_vctr_cpu<U>::value> {};

		template <class T, class U>
		struct is_grid_1d_and_vctr_gpu: std::integral_constant<dt_bool, is_grid_1d<T>::value && is_vctr_gpu<U>::value> {};

		template <class T, class U>
		struct is_grid_2d_and_vctr_cpu: std::integral_constant<dt_bool, is_grid_2d<T>::value && is_vctr_cpu<U>::value> {};

		template <class T, class U>
		struct is_grid_2d_and_vctr_gpu: std::integral_constant<dt_bool, is_grid_2d<T>::value && is_vctr_gpu<U>::value> {};

		template <class T, class U>
		struct is_grid_3d_and_vctr_cpu: std::integral_constant<dt_bool, is_grid_3d<T>::value && is_vctr_cpu<U>::value> {};

		template <class T, class U>
		struct is_grid_3d_and_vctr_gpu: std::integral_constant<dt_bool, is_grid_3d<T>::value && is_vctr_gpu<U>::value> {};

		template <class T, class U>
		struct is_grid_and_vctr_cpu: std::integral_constant<dt_bool, is_grid<T>::value && is_vctr_cpu<U>::value> {};

		template <class T, class U>
		struct is_grid_and_vctr_gpu: std::integral_constant<dt_bool, is_grid<T>::value && is_vctr_gpu<U>::value> {};

		/***************************************************************************************/
		template <class T, class U>
		struct is_grid_1d_and_cvctr_cpu: std::integral_constant<dt_bool, is_grid_1d<T>::value && is_cvctr_cpu<U>::value> {};

		template <class T, class U>
		struct is_grid_1d_and_cvctr_gpu: std::integral_constant<dt_bool, is_grid_1d<T>::value && is_cvctr_gpu<U>::value> {};

		template <class T, class U>
		struct is_grid_2d_and_cvctr_cpu: std::integral_constant<dt_bool, is_grid_2d<T>::value && is_cvctr_cpu<U>::value> {};

		template <class T, class U>
		struct is_grid_2d_and_cvctr_gpu: std::integral_constant<dt_bool, is_grid_2d<T>::value && is_cvctr_gpu<U>::value> {};

		template <class T, class U>
		struct is_grid_3d_and_cvctr_cpu: std::integral_constant<dt_bool, is_grid_3d<T>::value && is_cvctr_cpu<U>::value> {};

		template <class T, class U>
		struct is_grid_3d_and_cvctr_gpu: std::integral_constant<dt_bool, is_grid_3d<T>::value && is_cvctr_gpu<U>::value> {};

		template <class T, class U>
		struct is_grid_and_cvctr_cpu: std::integral_constant<dt_bool, is_grid<T>::value && is_cvctr_cpu<U>::value> {};

		template <class T, class U>
		struct is_grid_and_cvctr_gpu: std::integral_constant<dt_bool, is_grid<T>::value && is_cvctr_gpu<U>::value> {};

		/***************************************************************************************/
		template <class T, class U, class V=void>
		using enable_if_grid_1d_and_vctr_cpu = typename std::enable_if<is_grid_1d_and_vctr_cpu<T, U>::value, V>::type;

		template <class T, class U, class V=void>
		using enable_if_grid_1d_and_vctr_gpu = typename std::enable_if<is_grid_1d_and_vctr_gpu<T, U>::value, V>::type;

		template <class T, class U, class V=void>
		using enable_if_grid_2d_and_vctr_cpu = typename std::enable_if<is_grid_2d_and_vctr_cpu<T, U>::value, V>::type;

		template <class T, class U, class V=void>
		using enable_if_grid_2d_and_vctr_gpu = typename std::enable_if<is_grid_2d_and_vctr_gpu<T, U>::value, V>::type;

		template <class T, class U, class V=void>
		using enable_if_grid_3d_and_vctr_cpu = typename std::enable_if<is_grid_3d_and_vctr_cpu<T, U>::value, V>::type;

		template <class T, class U, class V=void>
		using enable_if_grid_3d_and_vctr_gpu = typename std::enable_if<is_grid_3d_and_vctr_gpu<T, U>::value, V>::type;

		template <class T, class U, class V=void>
		using enable_if_grid_nd_and_vctr_cpu = typename std::enable_if<is_grid_and_vctr_cpu<T, U>::value, V>::type;

		template <class T, class U, class V=void>
		using enable_if_grid_nd_and_vctr_gpu = typename std::enable_if<is_grid_and_vctr_gpu<T, U>::value, V>::type;

		/***************************************************************************************/
		template <class T, class U, class V=void>
		using enable_if_grid_1d_and_cvctr_cpu = typename std::enable_if<is_grid_1d_and_cvctr_cpu<T, U>::value, V>::type;

		template <class T, class U, class V=void>
		using enable_if_grid_1d_and_cvctr_gpu = typename std::enable_if<is_grid_1d_and_cvctr_gpu<T, U>::value, V>::type;

		template <class T, class U, class V=void>
		using enable_if_grid_2d_and_cvctr_cpu = typename std::enable_if<is_grid_2d_and_cvctr_cpu<T, U>::value, V>::type;

		template <class T, class U, class V=void>
		using enable_if_grid_2d_and_cvctr_gpu = typename std::enable_if<is_grid_2d_and_cvctr_gpu<T, U>::value, V>::type;

		template <class T, class U, class V=void>
		using enable_if_grid_3d_and_cvctr_cpu = typename std::enable_if<is_grid_3d_and_cvctr_cpu<T, U>::value, V>::type;

		template <class T, class U, class V=void>
		using enable_if_grid_3d_and_cvctr_gpu = typename std::enable_if<is_grid_3d_and_cvctr_gpu<T, U>::value, V>::type;

		template <class T, class U, class V=void>
		using enable_if_grid_nd_and_cvctr_cpu = typename std::enable_if<is_grid_and_cvctr_cpu<T, U>::value, V>::type;

		template <class T, class U, class V=void>
		using enable_if_grid_nd_and_cvctr_gpu = typename std::enable_if<is_grid_and_cvctr_gpu<T, U>::value, V>::type;
	}

#endif