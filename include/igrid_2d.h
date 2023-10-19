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

#include "igrid_1d.h"

/* derived class */
namespace mt
{
	template <class ST>
	using iGrid_2d_st = iGrid_sxd<ST, edim_2>;

	using iGrid_2d = iGrid_sxd<dt_int32, edim_2>;

	using iGrid_2d_64 = iGrid_sxd<dt_int64, edim_2>;
}

/* template specialization 2d */
namespace mt
{
	template <class ST>
	class iGrid_sxd<ST, edim_2>: public iGrid_sxd<ST, edim_1>
	{
	public:
		ST ny;			// number of pixels in y direction
		ST ny_h;		// half number of pixels in y direction

		/************************************* constructors ************************************/
		CGPU_EXEC
		iGrid_sxd();

		iGrid_sxd(const ST& nx, const ST& ny);

		/* copy constructor */
		CGPU_EXEC
		iGrid_sxd(const iGrid_sxd<ST, edim_2>& grid);

		/* converting constructor */
		template <class SU>
		CGPU_EXEC
		iGrid_sxd(const iGrid_sxd<SU, edim_2>& grid);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		iGrid_sxd<ST, edim_2>& operator=(const iGrid_sxd<ST, edim_2>& grid);

		/* converting assignment operator */
		template <class SU>
		CGPU_EXEC
		iGrid_sxd<ST, edim_2>& operator=(const iGrid_sxd<SU, edim_2>& grid);

		template <class SU> 
		CGPU_EXEC
		void assign(const iGrid_sxd<SU, edim_2>& grid);

		/***************************************************************************************/
		void set_size(const ST& nx, const ST& ny);

		CGPU_EXEC
		void clear();

		CGPU_EXEC
		virtual ST size() const;

		template <class SU>
		CGPU_EXEC
		SU ny_cast() const;

		dt_shape_st<ST> shape() const;

		/***************************************************************************************/
		/*********************** apply boundary conditions to indices **************************/
		/*            https:// en.wikipedia.org/wiki/Periodic_boundary_conditionsnn            */
		/***************************************************************************************/
		CGPU_EXEC
		ST ind_y_pbc(const ST& iy) const;

		CGPU_EXEC
		ST ind_y_bd(const ST& iy) const;

		CGPU_EXEC
		void ind_pbc(ST& ix, ST& iy) const;	

		/***************************************************************************************/
		/**************************** subindices -> linear indices *****************************/
		/***************************************************************************************/
		CGPU_EXEC
		ST operator()(const ST& ix, const ST& iy) const;

		CGPU_EXEC
		ST sub_2_ind(const ST& ix, const ST& iy) const;

		CGPU_EXEC
		ST sub_2_ind_bd(const ST& ix, const ST& iy) const;

		CGPU_EXEC
		ST sub_2_ind_pbc(const ST& ix, const ST& iy) const;

		CGPU_EXEC
		ST sub_2_ind_irv_sft_pbc(ST ix, ST iy) const;

		/***************************************************************************************/
		/**************************** linear indices-> subindices ******************************/
		/***************************************************************************************/
		CGPU_EXEC
		void ind_2_sub(const ST& ind, ST& ix, ST& iy) const;

		/***************************************************************************************/
		/******************* subindices to linear indices using 2nd dimension ******************/
		/***************************************************************************************/
		CGPU_EXEC
		ST sub_2_ind_by_d2(const ST& ix, const ST& iy) const;

		CGPU_EXEC
		ST sub_2_ind_by_d2_pbc(const ST& ix, const ST& iy) const;

		/***************************************************************************************/
		/****************** linear indices-> subindices using 2nd dimension ********************/
		/***************************************************************************************/
		CGPU_EXEC
		void ind_2_sub_by_d2(const ST& ind, ST& ix, ST& iy) const;

		/***************************************************************************************/
		/******************************* Fourier space indices *********************************/
		/***************************************************************************************/
		CGPU_EXEC
		ST igy(const ST& iy) const;

		CGPU_EXEC
		void igv(ST& ix, ST& iy) const;

		/************************************* shift *******************************************/
		CGPU_EXEC
		ST igy_sft(const ST& iy) const;

		CGPU_EXEC
		void igv_sft(ST& ix, ST& iy) const;

		/***************************************************************************************/
		/******************************* real space indices ************************************/
		/***************************************************************************************/
		CGPU_EXEC
		ST iry(const ST& iy) const;

		CGPU_EXEC
		void irv(ST& ix, ST& iy) const;

		/***************************************** shift ***************************************/
		CGPU_EXEC
		ST iry_sft(const ST& iy) const;

		CGPU_EXEC
		void irv_sft(ST& ix, ST& iy) const;

	#ifdef __CUDACC__
 		dim3 d_blk();	

		/***************************************************************************************/
 		dim3 d_grid(const dim3 d_grid_max = dim3(64, 64, 0));		

 		D_Grid_Blk d_grid_blk(const dim3 d_grid_max = dim3(64, 64, 1));		

		/***************************************************************************************/
 		dim3 d_grid_h(const dim3 d_grid_max = dim3(64, 64, 1));		

		D_Grid_Blk d_grid_blk_h(const dim3 d_grid_max = dim3(64, 64, 1));	

		/***************************************************************************************/
 		dim3 d_grid_d0_h(const dim3 d_grid_max = dim3(64, 64, 1));		

		D_Grid_Blk d_grid_blk_d0_h(const dim3 d_grid_max = dim3(64, 64, 1));
	#endif
	};
}
	
#include "../src/igrid_2d.inl"