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

#include "igrid_2d.h"

/* derived class */
namespace mt
{
	template <class ST>
	using iGrid_3d_st = iGrid_sxd<ST, edim_3>;

	using iGrid_3d = iGrid_sxd<dt_int32, edim_3>;

	using iGrid_3d_64 = iGrid_sxd<dt_int64, edim_3>;	
}

/* template specialization 3d */
namespace mt
{
	template <class ST>
	class iGrid_sxd<ST, edim_3>: public iGrid_sxd<ST, edim_2>
	{
	public:
		ST nz;			// number of pixels in z direction
		ST nz_h;		// half number of pixels in z direction

		/************************************* constructors ************************************/
		CGPU_EXEC
		iGrid_sxd();

		iGrid_sxd(const ST& nx, const ST& ny, const ST& nz);


		/* copy constructor */
		CGPU_EXEC
		iGrid_sxd(const iGrid_sxd<ST, edim_3>& grid);

		/* converting constructor */
		template <class SU>
		CGPU_EXEC
		iGrid_sxd(const iGrid_sxd<SU, edim_3>& grid);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		iGrid_sxd<ST, edim_3>& operator=(const iGrid_sxd<ST, edim_3>& grid);

		/* converting assignment operator */
		template <class SU>
		CGPU_EXEC
		iGrid_sxd<ST, edim_3>& operator=(const iGrid_sxd<SU, edim_3>& grid);

		template <class SU> 
		CGPU_EXEC
		void assign(const iGrid_sxd<SU, edim_3>& grid);

		/***************************************************************************************/
		void set_size(const ST& nx, const ST& ny, const ST& nz);

		CGPU_EXEC
		void clear();

		CGPU_EXEC
		ST size() const;

		template <class SU>
		CGPU_EXEC
		SU size_cast() const;

		template <class SU>
		CGPU_EXEC
		SU isize_cast() const;

		template <class SU>
		CGPU_EXEC
		SU nz_cast() const;
			
		dt_shape_st<ST> shape() const;

		/***************************************************************************************/
		/************************* apply boundary conditions to indices ************************/
		/*            https:// en.wikipedia.org/wiki/Periodic_boundary_conditionsnn            */
		/***************************************************************************************/
		CGPU_EXEC
		ST ind_z_pbc(const ST& iz) const;

		CGPU_EXEC
		ST ind_z_bd(const ST& iz) const;

		CGPU_EXEC
		void ind_pbc(ST& ix, ST& iy, ST& iz) const;

		/***************************************************************************************/
		/**************************** subindices -> linear indices *****************************/
		/***************************************************************************************/
		CGPU_EXEC
		ST operator()(const ST& ix, const ST& iy, const ST& iz) const;

		CGPU_EXEC
		ST sub_2_ind(const ST& ix, const ST& iy, const ST& iz) const;

		CGPU_EXEC
		ST sub_2_ind_bd(const ST& ix, const ST& iy, const ST& iz) const;

		CGPU_EXEC
		ST sub_2_ind_pbc(const ST& ix, const ST& iy, const ST& iz) const;

		CGPU_EXEC
		ST sub_2_ind_irv_sft_pbc(ST ix, ST iy, ST iz) const;

		/***************************************************************************************/
		//*************************** linear indices-> subindices ******************************/
		/***************************************************************************************/
		CGPU_EXEC
		void ind_2_sub(const ST& ind, ST& ix, ST& iy, ST& iz) const;

		/***************************************************************************************/
		/******************* subindices to linear indices using 2nd dimension ******************/
		/***************************************************************************************/
		CGPU_EXEC
		ST sub_2_ind_by_d2(const ST& ix, const ST& iy, const ST& iz) const;

		CGPU_EXEC
		ST sub_2_ind_by_d2_pbc(const ST& ix, const ST& iy, const ST& iz) const;

		/***************************************************************************************/
		/****************** linear indices-> subindices using 2nd dimension ********************/
		/***************************************************************************************/
		CGPU_EXEC
		void ind_2_sub_by_d2(const ST& ind, ST& ix, ST& iy, ST& iz) const;

		/***************************************************************************************/
		/****************************** Fourier space indices **********************************/
		/***************************************************************************************/
		CGPU_EXEC
		ST igz(const ST& iz) const;

		CGPU_EXEC
		void igv(ST& ix, ST& iy, ST& iz) const;

		/************************************* shift *******************************************/
		CGPU_EXEC
		ST igz_sft(const ST& iz) const;

		CGPU_EXEC
		void igv_sft(ST& ix, ST& iy, ST& iz) const;

		/***************************************************************************************/
		/******************************* real space indices ************************************/
		/***************************************************************************************/
		CGPU_EXEC
		ST irz(const ST& iz) const;

		CGPU_EXEC
		void irv(ST& ix, ST& iy, ST& iz) const;

		/************************************* shift *******************************************/
		CGPU_EXEC
		ST irz_sft(const ST& iz) const;

		CGPU_EXEC
		void irv_sft(ST& ix, ST& iy, ST& iz) const;

	#ifdef __CUDACC__
 		dim3 d_blk();
	
		/***************************************************************************************/
 		dim3 d_grid(const dim3 d_grid_max = dim3(64, 64, 64));
	
 		D_Grid_Blk d_grid_blk(const dim3 d_grid_max = dim3(64, 64, 64));
	
		/***************************************************************************************/
 		dim3 d_grid_h(const dim3 d_grid_max = dim3(64, 64, 64));
	
		D_Grid_Blk d_grid_blk_h(const dim3 d_grid_max = dim3(64, 64, 64));
	
		/***************************************************************************************/
 		dim3 d_grid_d0_h(const dim3 d_grid_max = dim3(64, 64, 64));
	
		D_Grid_Blk d_grid_blk_d0_h(const dim3 d_grid_max = dim3(64, 64, 64));
	#endif
	};													
}

#include "../src/igrid_3d.inl"