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
#include "fcns_cgpu_gen.h"

/* template definition */
namespace mt
{
	template <class ST, eDim Dim> class iGrid_sxd;

	template <eDim Dim>
	using iGrid_xd = iGrid_sxd<dt_int32, Dim>;
}

/* derived class */
namespace mt
{
	template <class ST>
	using iGrid_1d_st = iGrid_sxd<ST, edim_1>;

	using iGrid_1d = iGrid_sxd<dt_int32, edim_1>;

	using iGrid_1d_64 = iGrid_sxd<dt_int64, edim_1>;
}	

/* template specialization 1d */
namespace mt
{
	template <class ST>
	class iGrid_sxd<ST, edim_1>
	{
	public:
		ST nx;			// number of pixels in x direction
		ST nx_h;		// half number of pixels in x direction

		/************************************* constructors ************************************/
		CGPU_EXEC
		iGrid_sxd();

		iGrid_sxd(const ST& nx);

		/* copy constructor */
		CGPU_EXEC
		iGrid_sxd(const iGrid_sxd<ST, edim_1>& grid);

		/* converting constructor */
		template <class SU>
		CGPU_EXEC
		iGrid_sxd(const iGrid_sxd<SU, edim_1>& grid);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		iGrid_sxd<ST, edim_1>& operator=(const iGrid_sxd<ST, edim_1>& grid);

		/* converting assignment operator */
		template <class SU>
		CGPU_EXEC
		iGrid_sxd<ST, edim_1>& operator=(const iGrid_sxd<SU, edim_1>& grid);

		template <class SU> 
		CGPU_EXEC
		void assign(const iGrid_sxd<SU, edim_1>& grid);

		/***************************************************************************************/
		void set_size(const ST& nx);

		CGPU_EXEC
		void clear();

		CGPU_EXEC
		virtual ST size() const;

		template <class SU>
		CGPU_EXEC
		SU nx_cast() const;

		template <class SU>
		CGPU_EXEC
		SU size_cast() const;

		template <class SU>
		CGPU_EXEC
		SU isize_cast() const;

		dt_shape_st<ST> shape() const;

		/***************************************************************************************/
		/********************** apply boundary conditions to indices ***************************/
		/*********** https:// en.wikipedia.org/wiki/Periodic_boundary_conditionsnn *************/
		/***************************************************************************************/
		CGPU_EXEC
		ST ind_x_pbc(const ST& ix) const;	

		CGPU_EXEC
		ST ind_x_bd(const ST& ix) const;

		CGPU_EXEC
		void ind_pbc(ST& ix) const;

		/***************************************************************************************/
		/***************************** subindices -> linear indices ****************************/
		/***************************************************************************************/
		CGPU_EXEC
		ST operator()(const ST& ix) const;

		CGPU_EXEC
		ST sub_2_ind(const ST& ix) const;

		CGPU_EXEC
		ST sub_2_ind_bd(const ST& ix) const;

		CGPU_EXEC
		ST sub_2_ind_pbc(const ST& ix) const;

		CGPU_EXEC
		ST sub_2_ind_irv_sft_pbc(ST ix) const;

		/***************************************************************************************/
		/**************************** linear indices-> subindices ******************************/
		/***************************************************************************************/
		CGPU_EXEC
		void ind_2_sub(const ST& ind, ST& ix) const;

		/***************************************************************************************/
		/****************************** Fourier space indices **********************************/
		/***************************************************************************************/
		CGPU_EXEC
		ST igx(const ST& ix) const;

		/**************************************** shift ****************************************/
		CGPU_EXEC
		ST igx_sft(const ST& ix) const;

		/***************************************************************************************/
		/******************************** real space indices ***********************************/
		/***************************************************************************************/
		CGPU_EXEC
		ST irx(const ST& ix) const;

		/*************************************** shift *****************************************/
		CGPU_EXEC
		ST irx_sft(const ST& ix) const;

	#ifdef __CUDACC__
 		dim3 d_blk_size();
			
 		dim3 d_grid_size(const dim3 d_grid_max = dim3(128, 1, 1));
			
 		D_Grid_Blk d_grid_blk_size(const dim3 d_grid_max = dim3(128, 1, 1));
			
		/***************************************************************************************/
 		dim3 d_blk();
			
		/***************************************************************************************/		
 		dim3 d_grid(const dim3 d_grid_max = dim3(128, 1, 1));
			
 		D_Grid_Blk d_grid_blk(const dim3 d_grid_max = dim3(128, 1, 1));
			
		/***************************************************************************************/
 		dim3 d_grid_h(const dim3 d_grid_max = dim3(128, 1, 1));
			
		D_Grid_Blk d_grid_blk_h(const dim3 d_grid_max = dim3(128, 1, 1));
	#endif
	};
}

#include "../src/igrid_1d.inl"