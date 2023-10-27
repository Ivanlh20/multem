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

#include "igrid_1d.h"

/* template specialization igrid 1d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class ST>
	CGPU_EXEC
	iGrid_sxd<ST, edim_1>::iGrid_sxd(): nx(0), nx_h(0) {}

	template <class ST>
	iGrid_sxd<ST, edim_1>::iGrid_sxd(const ST& nx): nx(nx), nx_h(nx/ST(2)) {}

	/* copy constructor */
	template <class ST>
	CGPU_EXEC
	iGrid_sxd<ST, edim_1>::iGrid_sxd(const iGrid_sxd<ST, edim_1>& grid): iGrid_sxd()
	{
		*this = grid;
	}

	/* converting constructor */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	iGrid_sxd<ST, edim_1>::iGrid_sxd(const iGrid_sxd<SU, edim_1>& grid): iGrid_sxd()
	{
		*this = grid;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class ST>
	CGPU_EXEC
	iGrid_sxd<ST, edim_1>& iGrid_sxd<ST, edim_1>::operator=(const iGrid_sxd<ST, edim_1>& grid)
	{
		if (this != &grid)
		{
			nx = grid.nx;
			nx_h = grid.nx_h;
		}

		return *this;
	}

	/* converting assignment operator */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	iGrid_sxd<ST, edim_1>& iGrid_sxd<ST, edim_1>::operator=(const iGrid_sxd<SU, edim_1>& grid)
	{
		nx = ST(grid.nx);
		nx_h = ST(grid.nx_h);

		return *this;
	}

	template <class ST>
	template <class SU> 
	CGPU_EXEC
	void iGrid_sxd<ST, edim_1>::assign(const iGrid_sxd<SU, edim_1>& grid)
	{
		*this = grid;
	}

	/***************************************************************************************/
	template <class ST>
	void iGrid_sxd<ST, edim_1>::set_size(const ST& nx)
	{
		this->nx = nx;
		this->nx_h = nx/ST(2);
	}

	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_1>::clear()
	{
		nx = 0;
		nx_h = 0;
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_1>::size() const 
	{ 
		return nx;
	}

	template <class ST>
	template <class SU>
	CGPU_EXEC
	SU iGrid_sxd<ST, edim_1>::size_cast() const 
	{ 
		return SU(size());
	}

	template <class ST>
	template <class SU>
	CGPU_EXEC
	SU iGrid_sxd<ST, edim_1>::isize_cast() const
	{ 
		return SU(1)/size_cast<SU>();
	}

	template <class ST>
	template <class SU>
	CGPU_EXEC
	SU iGrid_sxd<ST, edim_1>::nx_cast() const
	{ 
		return SU(nx);
	}
	
	template <class ST>
	dt_shape_st<ST> iGrid_sxd<ST, edim_1>::shape() const
	{
		return {nx, ST(1), ST(1), ST(1)};
	}

	/***************************************************************************************/
	/********************** apply boundary conditions to indices ***************************/
	/*********** https:// en.wikipedia.org/wiki/Periodic_boundary_conditionsnn *************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_1>::ind_x_pbc(const ST& ix) const
	{
		return fcn_ind_pbc<ST>(ix, nx);
	}	

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_1>::ind_x_bd(const ST& ix) const
	{
		return fcn_set_bound(ix, ST(0), nx-ST(1));
	}	

	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_1>::ind_pbc(ST& ix) const
	{
		ix = ind_x_pbc(ix);
	}	

	/***************************************************************************************/
	/***************************** subindices -> linear indices ****************************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_1>::operator()(const ST& ix) const 
	{ 
		return ix;
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_1>::sub_2_ind(const ST& ix) const 
	{ 
		return ix;
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_1>::sub_2_ind_bd(const ST& ix) const 
	{ 
		return sub_2_ind(ind_x_bd(ix));
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_1>::sub_2_ind_pbc(const ST& ix) const 
	{ 
		return ind_x_pbc(ix);
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_1>::sub_2_ind_irv_sft_pbc(ST ix) const
	{ 
		ind_pbc(ix);
		irx_sft(ix);

		return sub_2_ind(ix);
	}

	/***************************************************************************************/
	/**************************** linear indices-> subindices ******************************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_1>::ind_2_sub(const ST& ind, ST& ix) const 
	{ 
		ix = ind;
	}

	/***************************************************************************************/
	/****************************** Fourier space indices **********************************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_1>::igx(const ST& ix) const 
	{ 
		return ix-nx_h;
	}
	/**************************************** shift ****************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_1>::igx_sft(const ST& ix) const 
	{ 
		return (ix<nx_h)?ix:(ix-nx);
	}

	/***************************************************************************************/
	/******************************** real space indices ***********************************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_1>::irx(const ST& ix) const 
	{ 
		return ix;
	}

	/*************************************** shift *****************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_1>::irx_sft(const ST& ix) const 
	{ 
		return (ix<nx_h)?(ix+nx_h):(ix-nx_h);
	}

#ifdef __CUDACC__
	template <class ST>
 	dim3 iGrid_sxd<ST, edim_1>::d_blk_size()
	{
		return fcn_cdb_size();
	}
			
	template <class ST>
 	dim3 iGrid_sxd<ST, edim_1>::d_grid_size(const dim3 d_grid_max)
	{
		auto grid = fcn_cdg_size(this->size());		
			
		grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;
			
		return grid;
	}
	
	template <class ST>
 	D_Grid_Blk iGrid_sxd<ST, edim_1>::d_grid_blk_size(const dim3 d_grid_max)
	{
		return D_Grid_Blk(d_grid_size(d_grid_max), d_blk_size());
	}
			
	/***************************************************************************************/
	template <class ST>
 	dim3 iGrid_sxd<ST, edim_1>::d_blk()	
	{
		return fcn_cdb_1d();		
	}
			
	/***************************************************************************************/
	template <class ST>
 	dim3 iGrid_sxd<ST, edim_1>::d_grid(const dim3 d_grid_max)
	{
		auto grid = fcn_cdg_1d(this->nx);			
			
		grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;
			
		return grid;
	}
			
	template <class ST>
 	D_Grid_Blk iGrid_sxd<ST, edim_1>::d_grid_blk(const dim3 d_grid_max)	
	{
		return D_Grid_Blk(d_grid(d_grid_max), d_blk());
	}
			
	/***************************************************************************************/
	template <class ST>
 	dim3 iGrid_sxd<ST, edim_1>::d_grid_h(const dim3 d_grid_max)
	{
		auto grid = fcn_cdg_1d(this->nx_h);
			
		grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;
			
		return grid;
	}
			
	template <class ST>
	D_Grid_Blk iGrid_sxd<ST, edim_1>::d_grid_blk_h(const dim3 d_grid_max)
	{
		return D_Grid_Blk(d_grid_h(d_grid_max), d_blk());
	}
#endif
}