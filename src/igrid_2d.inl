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

#include "igrid_2d.h"

/* template specialization 2d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class ST>
	CGPU_EXEC
	iGrid_sxd<ST, edim_2>::iGrid_sxd(): iGrid_sxd<ST, edim_1>(), ny(0), ny_h(0) {}

	template <class ST>
	iGrid_sxd<ST, edim_2>::iGrid_sxd(const ST& nx, const ST& ny): iGrid_sxd<ST, edim_1>(nx), ny(ny), ny_h(ny/ST(2)) {}

	/* copy constructor */
	template <class ST>
	CGPU_EXEC
	iGrid_sxd<ST, edim_2>::iGrid_sxd(const iGrid_sxd<ST, edim_2>& grid): iGrid_sxd()
	{
		*this = grid;
	}

	/* converting constructor */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	iGrid_sxd<ST, edim_2>::iGrid_sxd(const iGrid_sxd<SU, edim_2>& grid): iGrid_sxd()
	{
		*this = grid;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class ST>
	CGPU_EXEC
	iGrid_sxd<ST, edim_2>& iGrid_sxd<ST, edim_2>::operator=(const iGrid_sxd<ST, edim_2>& grid)
	{
		if (this != &grid)
		{
			iGrid_sxd<ST, edim_1>::operator=(grid);

			ny = grid.ny;
			ny_h = grid.ny_h;
		}

		return *this;
	}

	/* converting assignment operator */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	iGrid_sxd<ST, edim_2>& iGrid_sxd<ST, edim_2>::operator=(const iGrid_sxd<SU, edim_2>& grid)
	{
		iGrid_sxd<ST, edim_1>::operator=(grid);

		ny = ST(grid.ny);
		ny_h = ST(grid.ny_h);

		return *this;
	}

	template <class ST>
	template <class SU> 
	CGPU_EXEC
	void iGrid_sxd<ST, edim_2>::assign(const iGrid_sxd<SU, edim_2>& grid)
	{
		*this = grid;
	}

	/***************************************************************************************/
	template <class ST>
	void iGrid_sxd<ST, edim_2>::set_size(const ST& nx, const ST& ny)
	{
		iGrid_sxd<ST, edim_1>::set_size(nx);

		this->ny = ny;
		this->ny_h = ny/ST(2);
	}

	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_2>::clear()
	{
		iGrid_sxd<ST, edim_1>::clear();
		ny = 0;
		ny_h = 0;
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_2>::size() const 
	{ 
		return this->nx*ny;
	}

	template <class ST>
	template <class SU>
	CGPU_EXEC
	SU iGrid_sxd<ST, edim_2>::size_cast() const 
	{ 
		return SU(size());
	}

	template <class ST>
	template <class SU>
	CGPU_EXEC
	SU iGrid_sxd<ST, edim_2>::isize_cast() const
	{ 
		return SU(1)/size_cast<SU>();
	}

	template <class ST>
	template <class SU>
	CGPU_EXEC
	SU iGrid_sxd<ST, edim_2>::ny_cast() const 
	{ 
		return SU(ny);
	}

	template <class ST>
	dt_shape_st<ST> iGrid_sxd<ST, edim_2>::shape() const
	{
		return {ny, this->nx, ST(1), ST(1)};
	}

	/***************************************************************************************/
	/*********************** apply boundary conditions to indices **************************/
	/*           https:// en.wikipedia.org/wiki/Periodic_boundary_conditionsnn             */
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_2>::ind_y_pbc(const ST& iy) const
	{
		return fcn_ind_pbc<ST>(iy, ny);
	}	

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_2>::ind_y_bd(const ST& iy) const
	{
		return fcn_set_bound(iy, ST(0), ny-ST(1));
	}	

	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_2>::ind_pbc(ST& ix, ST& iy) const
	{
		ix = this->ind_x_pbc(ix);
		iy = ind_y_pbc(iy);
	}	

	/***************************************************************************************/
	/**************************** subindices -> linear indices *****************************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_2>::operator()(const ST& ix, const ST& iy) const 
	{ 
		return iy + ix*ny;
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_2>::sub_2_ind(const ST& ix, const ST& iy) const 
	{ 
		return iy + ix*ny;
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_2>::sub_2_ind_bd(const ST& ix, const ST& iy) const 
	{ 
		return sub_2_ind(this->ind_x_bd(ix), ind_y_bd(iy));
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_2>::sub_2_ind_pbc(const ST& ix, const ST& iy) const 
	{ 
		return ind_y_pbc(iy) + this->ind_x_pbc(ix)*ny;
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_2>::sub_2_ind_irv_sft_pbc(ST ix, ST iy) const
	{ 
		ind_pbc(ix, iy);
		irv_sft(ix, iy);

		return sub_2_ind(ix, iy);
	}

	/***************************************************************************************/
	/**************************** linear indices-> subindices ******************************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_2>::ind_2_sub(const ST& ind, ST& ix, ST& iy) const 
	{ 
		ix = ind/ny;
		iy = ind - ix*ny;
	}

	/***************************************************************************************/
	/******************* subindices to linear indices using 2nd dimension ******************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_2>::sub_2_ind_by_d2(const ST& ix, const ST& iy) const 
	{ 
		return iy*this->nx + ix;
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_2>::sub_2_ind_by_d2_pbc(const ST& ix, const ST& iy) const 
	{ 
		return ind_y_pbc(iy)*this->nx + this->ind_x_pbc(ix);
	}

	/***************************************************************************************/
	/****************** linear indices-> subindices using 2nd dimension ********************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_2>::ind_2_sub_by_d2(const ST& ind, ST& ix, ST& iy) const 
	{ 
		iy = ind/this->nx;
		ix = ind - iy*this->nx;
	}

	/***************************************************************************************/
	/******************************* Fourier space indices *********************************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_2>::igy(const ST& iy) const 
	{ 
		return iy - ny_h;
	}

	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_2>::igv(ST& ix, ST& iy) const 
	{ 
		ix = this->igx(ix);
		iy = igy(iy);
	}

	/************************************* shift *******************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_2>::igy_sft(const ST& iy) const 
	{ 
		return (iy<ny_h)?iy:(iy-ny);
	}

	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_2>::igv_sft(ST& ix, ST& iy) const 
	{ 
		ix = this->igx_sft(ix);
		iy = igy_sft(iy);
	}

	/***************************************************************************************/
	/******************************* real space indices ************************************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_2>::iry(const ST& iy) const 
	{ 
		return iy;
	}

	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_2>::irv(ST& ix, ST& iy) const 
	{ 
		ix = this->irx(ix);
		iy = iry(iy);
	}

	/***************************************** shift ***************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_2>::iry_sft(const ST& iy) const 
	{ 
		return (iy<ny_h)?(iy+ny_h):(iy-ny_h);
	}

	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_2>::irv_sft(ST& ix, ST& iy) const
	{ 
		ix = this->irx_sft(ix);
		iy = iry_sft(iy);
	}

#ifdef __CUDACC__
	template <class ST>
 	dim3 iGrid_sxd<ST, edim_2>::d_blk_size()
	{
		return fcn_cdb_size();
	}
			
	template <class ST>
 	dim3 iGrid_sxd<ST, edim_2>::d_grid_size(const dim3 d_grid_max)
	{
		auto grid = fcn_cdg_size(size());		
			
		grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;
			
		return grid;
	}
	
	template <class ST>
 	D_Grid_Blk iGrid_sxd<ST, edim_2>::d_grid_blk_size(const dim3 d_grid_max)
	{
		return D_Grid_Blk(d_grid_size(d_grid_max), d_blk_size());
	}
			
	/***************************************************************************************/
	template <class ST>
 	dim3 iGrid_sxd<ST, edim_2>::d_blk()				
	{		
		return fcn_cdb_2d();	
	}		

	/***************************************************************************************/
	template <class ST>
 	dim3 iGrid_sxd<ST, edim_2>::d_grid(const dim3 d_grid_max)				
	{		
		auto grid = fcn_cdg_2d(this->ny, this->nx);	
				
		grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;	
		grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;	
				
		return grid;			
	}		

	template <class ST>
 	D_Grid_Blk iGrid_sxd<ST, edim_2>::d_grid_blk(const dim3 d_grid_max)		
	{		
		return D_Grid_Blk(d_grid(d_grid_max), d_blk());
	}		

	/***************************************************************************************/
 	template <class ST>
	dim3  iGrid_sxd<ST, edim_2>::d_grid_h(const dim3 d_grid_max)				
	{		
		auto grid = fcn_cdg_2d(this->ny_h, this->nx_h);
				
		grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;	
		grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;	
				
		return grid;			
	}		

	template <class ST>
	D_Grid_Blk iGrid_sxd<ST, edim_2>::d_grid_blk_h(const dim3 d_grid_max)	
	{		
		return D_Grid_Blk(d_grid_h(d_grid_max), d_blk());				
	}		

	/***************************************************************************************/
	template <class ST>
 	dim3 iGrid_sxd<ST, edim_2>::d_grid_d0_h(const dim3 d_grid_max)			
	{		
		auto grid = fcn_cdg_2d(this->ny_h, this->nx);
				
		grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;	
		grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;	
				
		return grid;			
	}		

	template <class ST>
	D_Grid_Blk iGrid_sxd<ST, edim_2>::d_grid_blk_d0_h(const dim3 d_grid_max)	
	{		
		return D_Grid_Blk(d_grid_d0_h(d_grid_max), d_blk());			
	}
#endif
}