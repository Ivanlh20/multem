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

#include "igrid_3d.h"

/* template specialization 3d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class ST>
	CGPU_EXEC
	iGrid_sxd<ST, edim_3>::iGrid_sxd(): iGrid_sxd<ST, edim_2>(), nz(0), nz_h(0) {}

	template <class ST>
	iGrid_sxd<ST, edim_3>::iGrid_sxd(const ST& nx, const ST& ny, const ST& nz): iGrid_sxd<ST, edim_2>(nx, ny), nz(nz), nz_h(nz/ST(2)) {}

	/* copy constructor */
	template <class ST>
	CGPU_EXEC
	iGrid_sxd<ST, edim_3>::iGrid_sxd(const iGrid_sxd<ST, edim_3>& grid): iGrid_sxd()
	{
		*this = grid;
	}

	/* converting constructor */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	iGrid_sxd<ST, edim_3>::iGrid_sxd(const iGrid_sxd<SU, edim_3>& grid): iGrid_sxd()
	{
		*this = grid;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class ST>
	CGPU_EXEC
	iGrid_sxd<ST, edim_3>& iGrid_sxd<ST, edim_3>::operator=(const iGrid_sxd<ST, edim_3>& grid)
	{
		if (this != &grid)
		{
			iGrid_sxd<ST, edim_2>::operator=(grid);
			nz = grid.nz;
			nz_h = grid.nz_h;
		}

		return *this;
	}

	/* converting assignment operator */
	template <class ST>
	template <class SU>
	CGPU_EXEC
	iGrid_sxd<ST, edim_3>& iGrid_sxd<ST, edim_3>::operator=(const iGrid_sxd<SU, edim_3>& grid)
	{
		iGrid_sxd<ST, edim_2>::operator=(grid);
		nz = ST(grid.nz);
		nz_h = ST(grid.nz_h);

		return *this;
	}

	template <class ST>
	template <class SU> 
	CGPU_EXEC
	void iGrid_sxd<ST, edim_3>::assign(const iGrid_sxd<SU, edim_3>& grid)
	{
		*this = grid;
	}

	/***************************************************************************************/
	template <class ST>
	void iGrid_sxd<ST, edim_3>::set_size(const ST& nx, const ST& ny, const ST& nz)
	{
		iGrid_sxd<ST, edim_2>::set_size(nx, ny);

		this->nz = nz;
		this->nz_h = nz/ST(2);
	}

	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_3>::clear()
	{
		iGrid_sxd<ST, edim_2>::clear();
		nz = 0;
		nz_h = 0;
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_3>::size() const 
	{ 
		return this->nx*this->ny*nz;
	}

	template <class ST>
	template <class SU>
	CGPU_EXEC
	SU iGrid_sxd<ST, edim_3>::nz_cast() const 
	{ 
		return SU(nz);
	}
			
	template <class ST>
	dt_shape_st<ST> iGrid_sxd<ST, edim_3>::shape() const
	{
		return {this->ny, this->nx, nz, ST(1)};
	}

	/***************************************************************************************/
	/************************* apply boundary conditions to indices ************************/
	/*           https:// en.wikipedia.org/wiki/Periodic_boundary_conditionsnn             */
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_3>::ind_z_pbc(const ST& iz) const
	{
		return fcn_ind_pbc<ST>(iz, nz);
	}	

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_3>::ind_z_bd(const ST& iz) const
	{
		return fcn_set_bound(iz, ST(0), nz-ST(1));
	}	

	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_3>::ind_pbc(ST& ix, ST& iy, ST& iz) const
	{
		ix = ind_x_pbc(ix);
		iy = ind_y_pbc(iy);
		iz = ind_z_pbc(iz);
	}	

	/***************************************************************************************/
	/**************************** subindices -> linear indices *****************************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_3>::operator()(const ST& ix, const ST& iy, const ST& iz) const 
	{ 
		return iy + ix*this->ny + iz*this->ny*this->nx;
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_3>::sub_2_ind(const ST& ix, const ST& iy, const ST& iz) const 
	{ 
		return iy + ix*this->ny + iz*this->ny*this->nx;
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_3>::sub_2_ind_bd(const ST& ix, const ST& iy, const ST& iz) const 
	{ 
		return sub_2_ind(this->ind_x_bd(ix), this->ind_y_bd(iy), ind_z_bd(iz));
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_3>::sub_2_ind_pbc(const ST& ix, const ST& iy, const ST& iz) const 
	{ 
		return this->ind_y_pbc(iy) + this->ind_x_pbc(ix)*this->ny + ind_z_pbc(iz)*this->ny*this->nx;
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_3>::sub_2_ind_irv_sft_pbc(ST ix, ST iy, ST iz) const
	{ 
		ind_pbc(ix, iy, iz);
		irv_sft(ix, iy, iz);

		return sub_2_ind(ix, iy, iz);
	}

	/***************************************************************************************/
	/**************************** linear indices-> subindices ******************************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_3>::ind_2_sub(const ST& ind, ST& ix, ST& iy, ST& iz) const 
	{ 
		iz = ind/(this->ny*this->nx);
		iy = ind - iz*this->ny*this->nx;
		ix = iy/this->ny;
		iy = iy - ix*this->ny;
	}

	/***************************************************************************************/
	/******************* subindices to linear indices using 2nd dimension ******************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_3>::sub_2_ind_by_d2(const ST& ix, const ST& iy, const ST& iz) const 
	{ 
		return iy*this->nx + ix + iz*this->ny*this->nx;
	}

	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_3>::sub_2_ind_by_d2_pbc(const ST& ix, const ST& iy, const ST& iz) const 
	{ 
		return this->ind_y_pbc(iy)*this->nx + this->ind_x_pbc(ix) + ind_z_pbc(iz)*this->ny*this->nx;
	}

	/***************************************************************************************/
	/****************** linear indices-> subindices using 2nd dimension ********************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_3>::ind_2_sub_by_d2(const ST& ind, ST& ix, ST& iy, ST& iz) const 
	{ 
		iz = ind/(this->ny*this->nx);
		ix = ind - iz*this->ny*this->nx;
		iy = ix/this->nx;
		ix = ix - iy*this->nx;
	}

	/***************************************************************************************/
	/****************************** Fourier space indices **********************************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_3>::igz(const ST& iz) const 
	{ 
		return iz - nz_h;
	}

	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_3>::igv(ST& ix, ST& iy, ST& iz) const 
	{ 
		ix = this->igx(ix);
		iy = this->igy(iy);
		iz = igz(iz);
	}

	/************************************* shift *******************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_3>::igz_sft(const ST& iz) const 
	{ 
		return (iz<nz_h)?iz:(iz-nz);
	}

	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_3>::igv_sft(ST& ix, ST& iy, ST& iz) const 
	{ 
		ix = this->igx_sft(ix);
		iy = this->igy_sft(iy);
		iz = igz_sft(iz);
	}

	/***************************************************************************************/
	/******************************* real space indices ************************************/
	/***************************************************************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_3>::irz(const ST& iz) const 
	{ 
		return iz;
	}

	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_3>::irv(ST& ix, ST& iy, ST& iz) const 
	{ 
		ix = this->irx(ix);
		iy = this->iry(iy);
		iz = irz(iz);
	}

	/************************************* shift *******************************************/
	template <class ST>
	CGPU_EXEC
	ST iGrid_sxd<ST, edim_3>::irz_sft(const ST& iz) const 
	{ 
		return (iz<nz_h)?(iz+nz_h):(iz-nz_h);
	}

	template <class ST>
	CGPU_EXEC
	void iGrid_sxd<ST, edim_3>::irv_sft(ST& ix, ST& iy, ST& iz) const
	{ 
		ix = this->irx_sft(ix);
		iy = this->iry_sft(iy);
		iz = irz_sft(iz);
	}

#ifdef __CUDACC__
template <class ST>
 	dim3 iGrid_sxd<ST, edim_3>::d_blk()
	{
		return fcn_cdb_3d();
	}
	
	/***************************************************************************************/
	template <class ST>
 	dim3 iGrid_sxd<ST, edim_3>::d_grid(const dim3 d_grid_max)
	{
		auto grid = fcn_cdg_3d(this->ny, this->nx, this->nz);
		grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;
		grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;
		grid.z = (d_grid_max.z > 0)?min(d_grid_max.z, grid.z):grid.z;
	
		return grid;
	}
	
	template <class ST>
 	D_Grid_Blk iGrid_sxd<ST, edim_3>::d_grid_blk(const dim3 d_grid_max)
	{
		return D_Grid_Blk(d_grid(d_grid_max), d_blk());
	}
	
	/***************************************************************************************/
	template <class ST>
 	dim3 iGrid_sxd<ST, edim_3>::d_grid_h(const dim3 d_grid_max)
	{
		auto grid = fcn_cdg_3d(this->ny_h, this->nx_h, this->nz_h);
	
		grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;
		grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;
		grid.z = (d_grid_max.z > 0)?min(d_grid_max.z, grid.z):grid.z;
	
		return grid;
	}
	
	template <class ST>
	D_Grid_Blk iGrid_sxd<ST, edim_3>::d_grid_blk_h(const dim3 d_grid_max)
	{
		return D_Grid_Blk(d_grid_h(d_grid_max), d_blk());
	}
	
	/***************************************************************************************/
	template <class ST>
 	dim3 iGrid_sxd<ST, edim_3>::d_grid_d0_h(const dim3 d_grid_max)
	{
		auto grid = fcn_cdg_3d(this->ny_h, this->nx, this->nz);
	
		grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;
		grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;
		grid.z = (d_grid_max.z > 0)?min(d_grid_max.z, grid.z):grid.z;
	
		return grid;
	}
	
	template <class ST>
	D_Grid_Blk iGrid_sxd<ST, edim_3>::d_grid_blk_d0_h(const dim3 d_grid_max)
	{
		return D_Grid_Blk(d_grid_d0_h(d_grid_max), d_blk());
	}
#endif												
}