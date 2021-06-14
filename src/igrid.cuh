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

#ifndef IGRID_H
	#define IGRID_H
 
	#ifdef _MSC_VER
		#pragma once
	#endif

	#include "const_enum.cuh"
	#include "cgpu_fcns_gen.cuh"

	/* macro gpu grid-block */
	namespace mt
	{
	#define	FCNS_GPU_GRID_BLK_IGRID_1D																		\
 		dim3 d_blk_size()																					\
		{																									\
			return fcn_cdb_size();																			\
		}																									\
																											\
 		dim3 d_grid_size(const dim3 d_grid_max = dim3(128, 1, 1))											\
		{																									\
			auto grid = fcn_cdg_size(this->m_size);															\
																											\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;									\
																											\
			return grid;																					\
		}																									\
																											\
 		D_Grid_Blk d_grid_blk_size(const dim3 d_grid_max = dim3(128, 1, 1))									\
		{																									\
			return D_Grid_Blk(d_grid_size(d_grid_max), d_blk_size());										\
		}																									\
																											\
		/***************************************************************************************/			\
 		dim3 d_blk()																						\
		{																									\
			return fcn_cdb_1d();																			\
		}																									\
																											\
		/***************************************************************************************/			\
 		dim3 d_grid(const dim3 d_grid_max = dim3(128, 1, 1))												\
		{																									\
			auto grid = fcn_cdg_1d(this->nx);																\
																											\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;									\
																											\
			return grid;																					\
		}																									\
																											\
 		D_Grid_Blk d_grid_blk(const dim3 d_grid_max = dim3(128, 1, 1))										\
		{																									\
			return D_Grid_Blk(d_grid(d_grid_max), d_blk());													\
		}																									\
																											\
		/***************************************************************************************/			\
 		dim3 d_grid_h(const dim3 d_grid_max = dim3(128, 1, 1))												\
		{																									\
			auto grid = fcn_cdg_1d(this->nx_h);																\
																											\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;									\
																											\
			return grid;																					\
		}																									\
																											\
		D_Grid_Blk d_grid_blk_h(const dim3 d_grid_max = dim3(128, 1, 1))									\
		{																									\
			return D_Grid_Blk(d_grid_h(d_grid_max), d_blk());												\
		}
		
	#define	FCNS_GPU_GRID_BLK_IGRID_2D																		\
 		dim3 d_blk()																						\
		{																									\
			return fcn_cdb_2d();																			\
		}																									\
																											\
		/***************************************************************************************/			\
 		dim3 d_grid(const dim3 d_grid_max = dim3(64, 64, 0))												\
		{																									\
			auto grid = fcn_cdg_2d(this->ny, this->nx);														\
																											\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;									\
			grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;									\
																											\
			return grid;																					\
		}																									\
																											\
 		D_Grid_Blk d_grid_blk(const dim3 d_grid_max = dim3(64, 64, 1))										\
		{																									\
			return D_Grid_Blk(d_grid(d_grid_max), d_blk());													\
		}																									\
																											\
		/***************************************************************************************/			\
 		dim3 d_grid_h(const dim3 d_grid_max = dim3(64, 64, 1))												\
		{																									\
			auto grid = fcn_cdg_2d(this->ny_h, this->nx_h);													\
																											\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;									\
			grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;									\
																											\
			return grid;																					\
		}																									\
																											\
		D_Grid_Blk d_grid_blk_h(const dim3 d_grid_max = dim3(64, 64, 1))									\
		{																									\
			return D_Grid_Blk(d_grid_h(d_grid_max), d_blk());												\
		}																									\
																											\
		/***************************************************************************************/			\
 		dim3 d_grid_d0_h(const dim3 d_grid_max = dim3(64, 64, 1))											\
		{																									\
			auto grid = fcn_cdg_2d(this->ny_h, this->nx);													\
																											\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;									\
			grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;									\
																											\
			return grid;																					\
		}																									\
																											\
		D_Grid_Blk d_grid_blk_d0_h(const dim3 d_grid_max = dim3(64, 64, 1))									\
		{																									\
			return D_Grid_Blk(d_grid_d0_h(d_grid_max), d_blk());											\
		}

	#define	FCNS_GPU_GRID_BLK_IGRID_3D																		\
 		dim3 d_blk()																						\
		{																									\
			return fcn_cdb_3d();																			\
		}																									\
																											\
		/***************************************************************************************/			\
 		dim3 d_grid(const dim3 d_grid_max = dim3(64, 64, 64))												\
		{																									\
			auto grid = fcn_cdg_3d(this->ny, this->nx, this->nz);											\
																											\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;									\
			grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;									\
			grid.z = (d_grid_max.z > 0)?min(d_grid_max.z, grid.z):grid.z;									\
																											\
			return grid;																					\
		}																									\
																											\
 		D_Grid_Blk d_grid_blk(const dim3 d_grid_max = dim3(64, 64, 64))										\
		{																									\
			return D_Grid_Blk(d_grid(d_grid_max), d_blk());													\
		}																									\
																											\
		/***************************************************************************************/			\
 		dim3 d_grid_h(const dim3 d_grid_max = dim3(64, 64, 64))												\
		{																									\
			auto grid = fcn_cdg_3d(this->ny_h, this->nx_h, this->nz_h);										\
																											\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;									\
			grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;									\
			grid.z = (d_grid_max.z > 0)?min(d_grid_max.z, grid.z):grid.z;									\
																											\
			return grid;																					\
		}																									\
																											\
		D_Grid_Blk d_grid_blk_h(const dim3 d_grid_max = dim3(64, 64, 64))									\
		{																									\
			return D_Grid_Blk(d_grid_h(d_grid_max), d_blk());												\
		}																									\
																											\
		/***************************************************************************************/			\
 		dim3 d_grid_d0_h(const dim3 d_grid_max = dim3(64, 64, 64))											\
		{																									\
			auto grid = fcn_cdg_3d(this->ny_h, this->nx, this->nz);											\
																											\
			grid.x = (d_grid_max.x > 0)?min(d_grid_max.x, grid.x):grid.x;									\
			grid.y = (d_grid_max.y > 0)?min(d_grid_max.y, grid.y):grid.y;									\
			grid.z = (d_grid_max.z > 0)?min(d_grid_max.z, grid.z):grid.z;									\
																											\
			return grid;																					\
		}																									\
																											\
		D_Grid_Blk d_grid_blk_d0_h(const dim3 d_grid_max = dim3(64, 64, 64))								\
		{																									\
			return D_Grid_Blk(d_grid_d0_h(d_grid_max), d_blk());											\
		}
	}
	
	/***************************************************************************************/
	/*************************************** igrid *****************************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class ST, eDim Dim> class iGrid_sxd;

		template <eDim Dim>
		using iGrid_xd = iGrid_sxd<dt_int32, Dim>;

		/* 1d */
		template<class ST>
		using iGrid_1d_st = iGrid_sxd<ST, edim_1>;

		using iGrid_1d = iGrid_sxd<dt_int32, edim_1>;

		using iGrid_1d_64 = iGrid_sxd<dt_int64, edim_1>;

		/* 2d */
		template<class ST>
		using iGrid_2d_st = iGrid_sxd<ST, edim_2>;

		using iGrid_2d = iGrid_sxd<dt_int32, edim_2>;

		using iGrid_2d_64 = iGrid_sxd<dt_int64, edim_2>;	

		/* 3d */
		template<class ST>
		using iGrid_3d_st = iGrid_sxd<ST, edim_3>;

		using iGrid_3d = iGrid_sxd<dt_int32, edim_3>;

		using iGrid_3d_64 = iGrid_sxd<dt_int64, edim_3>;	
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
			iGrid_sxd(): nx(0), nx_h(0) {}

			iGrid_sxd(const ST& nx): nx(nx), nx_h(nx/ST(2)) {}

			/* copy constructor */
			CGPU_EXEC
			iGrid_sxd(const iGrid_sxd<ST, edim_1>& grid): iGrid_sxd()
			{
				*this = grid;
			}

			/* converting constructor */
			template <class SU>
			CGPU_EXEC
			iGrid_sxd(const iGrid_sxd<SU, edim_1>& grid): iGrid_sxd()
			{
				*this = grid;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			iGrid_sxd<ST, edim_1>& operator=(const iGrid_sxd<ST, edim_1>& grid)
			{
				if (this != &grid)
				{
					nx = grid.nx;
					nx_h = grid.nx_h;
				}

				return *this;
			}

			/* converting assignment operator */
			template <class SU>
			CGPU_EXEC
			iGrid_sxd<ST, edim_1>& operator=(const iGrid_sxd<SU, edim_1>& grid)
			{
				nx = ST(grid.nx);
				nx_h = ST(grid.nx_h);

				return *this;
			}

			template <class SU> 
			CGPU_EXEC
			void assign(const iGrid_sxd<SU, edim_1>& grid)
			{
				*this = grid;
			}

			/***************************************************************************************/
			void set_size(const ST& nx)
			{
				this->nx = nx;
				this->nx_h = nx/ST(2);
			}

			CGPU_EXEC
			void clear()
			{
				nx = 0;
				nx_h = 0;
			}

			CGPU_EXEC
			virtual ST size() const 
			{ 
				return nx;
			}

			template<class SU>
			CGPU_EXEC
			SU nx_cast() const
			{ 
				return SU(nx);
			}

			template<class SU>
			CGPU_EXEC
			SU size_cast() const 
			{ 
				return SU(size());
			}

			template<class SU>
			CGPU_EXEC
			SU isize_cast() const
			{ 
				return SU(1)/size_cast<SU>();
			}

			dt_shape_st<ST> shape() const
			{
				return {nx, ST(1), ST(1), ST(1)};
			}

			/***************************************************************************************/
			/********************** apply boundary conditions to indices ***************************/
			/*********** https:// en.wikipedia.org/wiki/Periodic_boundary_conditionsnn *************/
			/***************************************************************************************/
			CGPU_EXEC
			ST ind_x_pbc(const ST& ix) const
			{
				return fcn_ind_pbc<ST>(ix, nx);
			}	

			CGPU_EXEC
			ST ind_x_bd(const ST& ix) const
			{
				return fcn_set_bound(ix, ST(0), nx-ST(1));
			}	

			CGPU_EXEC
			void ind_pbc(ST& ix) const
			{
				ix = ind_x_pbc(ix);
			}	

			/***************************************************************************************/
			/***************************** subindices -> linear indices ****************************/
			/***************************************************************************************/
			CGPU_EXEC
			ST operator()(const ST& ix) const 
			{ 
				return ix;
			}

			CGPU_EXEC
			ST sub_2_ind(const ST& ix) const 
			{ 
				return ix;
			}

			CGPU_EXEC
			ST sub_2_ind_bd(const ST& ix) const 
			{ 
				return sub_2_ind(ind_x_bd(ix));
			}

			CGPU_EXEC
			ST sub_2_ind_pbc(const ST& ix) const 
			{ 
				return ind_x_pbc(ix);
			}

			CGPU_EXEC
			ST sub_2_ind_irv_sft_pbc(ST ix) const
			{ 
				ind_pbc(ix);
				irx_sft(ix);

				return sub_2_ind(ix);
			}

			/***************************************************************************************/
			/**************************** linear indices-> subindices ******************************/
			/***************************************************************************************/
			CGPU_EXEC
			void ind_2_sub(const ST& ind, ST& ix) const 
			{ 
				ix = ind;
			}

			/***************************************************************************************/
			/****************************** Fourier space indices **********************************/
			/***************************************************************************************/
			CGPU_EXEC
			ST igx(const ST& ix) const 
			{ 
				return ix-nx_h;
			}
			/**************************************** shift ****************************************/
			CGPU_EXEC
			ST igx_sft(const ST& ix) const 
			{ 
				return (ix<nx_h)?ix:(ix-nx);
			}

			/***************************************************************************************/
			/******************************** real space indices ***********************************/
			/***************************************************************************************/
			CGPU_EXEC
			ST irx(const ST& ix) const 
			{ 
				return ix;
			}

			/*************************************** shift *****************************************/
			CGPU_EXEC
			ST irx_sft(const ST& ix) const 
			{ 
				return (ix<nx_h)?(ix+nx_h):(ix-nx_h);
			}

		#ifdef __CUDACC__
			FCNS_GPU_GRID_BLK_IGRID_1D;
		#endif
		};
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
			iGrid_sxd(): iGrid_sxd<ST, edim_1>(), ny(0), ny_h(0) {}

			iGrid_sxd(const ST& nx, const ST& ny): iGrid_sxd<ST, edim_1>(nx), ny(ny), ny_h(ny/ST(2)) {}

			/* copy constructor */
			CGPU_EXEC
			iGrid_sxd(const iGrid_sxd<ST, edim_2>& grid): iGrid_sxd()
			{
				*this = grid;
			}

			/* converting constructor */
			template <class SU>
			CGPU_EXEC
			iGrid_sxd(const iGrid_sxd<SU, edim_2>& grid): iGrid_sxd()
			{
				*this = grid;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			iGrid_sxd<ST, edim_2>& operator=(const iGrid_sxd<ST, edim_2>& grid)
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
			template <class SU>
			CGPU_EXEC
			iGrid_sxd<ST, edim_2>& operator=(const iGrid_sxd<SU, edim_2>& grid)
			{
				iGrid_sxd<ST, edim_1>::operator=(grid);

				ny = ST(grid.ny);
				ny_h = ST(grid.ny_h);

				return *this;
			}

			template <class SU> 
			CGPU_EXEC
			void assign(const iGrid_sxd<SU, edim_2>& grid)
			{
				*this = grid;
			}

			/***************************************************************************************/
			void set_size(const ST& nx, const ST& ny)
			{
				iGrid_sxd<ST, edim_1>::set_size(nx);

				this->ny = ny;
				this->ny_h = ny/ST(2);
			}

			CGPU_EXEC
			void clear()
			{
				iGrid_sxd<ST, edim_1>::clear();
				ny = 0;
				ny_h = 0;
			}

			CGPU_EXEC
			virtual ST size() const 
			{ 
				return this->nx*ny;
			}

			template<class SU>
			CGPU_EXEC
			SU ny_cast() const 
			{ 
				return SU(ny);
			}

			dt_shape_st<ST> shape() const
			{
				return {ny, this->nx, ST(1), ST(1)};
			}

			/***************************************************************************************/
			/*********************** apply boundary conditions to indices **************************/
			/*nnhttps:// en.wikipedia.org/wiki/Periodic_boundary_conditionsnn */
			/***************************************************************************************/
			CGPU_EXEC
			ST ind_y_pbc(const ST& iy) const
			{
				return fcn_ind_pbc<ST>(iy, ny);
			}	

			CGPU_EXEC
			ST ind_y_bd(const ST& iy) const
			{
				return fcn_set_bound(iy, ST(0), ny-ST(1));
			}	

			CGPU_EXEC
			void ind_pbc(ST& ix, ST& iy) const
			{
				ix = this->ind_x_pbc(ix);
				iy = ind_y_pbc(iy);
			}	

			/***************************************************************************************/
			/**************************** subindices -> linear indices *****************************/
			/***************************************************************************************/
			CGPU_EXEC
			ST operator()(const ST& ix, const ST& iy) const 
			{ 
				return iy + ix*ny;
			}

			CGPU_EXEC
			ST sub_2_ind(const ST& ix, const ST& iy) const 
			{ 
				return iy + ix*ny;
			}

			CGPU_EXEC
			ST sub_2_ind_bd(const ST& ix, const ST& iy) const 
			{ 
				return sub_2_ind(this->ind_x_bd(ix), ind_y_bd(iy));
			}

			CGPU_EXEC
			ST sub_2_ind_pbc(const ST& ix, const ST& iy) const 
			{ 
				return ind_y_pbc(iy) + this->ind_x_pbc(ix)*ny;
			}

			CGPU_EXEC
			ST sub_2_ind_irv_sft_pbc(ST ix, ST iy) const
			{ 
				ind_pbc(ix, iy);
				irv_sft(ix, iy);

				return sub_2_ind(ix, iy);
			}

			/***************************************************************************************/
			/**************************** linear indices-> subindices ******************************/
			/***************************************************************************************/
			CGPU_EXEC
			void ind_2_sub(const ST& ind, ST& ix, ST& iy) const 
			{ 
				ix = ind/ny;
				iy = ind - ix*ny;
			}

			/***************************************************************************************/
			/******************* subindices to linear indices using 2nd dimension ******************/
			/***************************************************************************************/
			CGPU_EXEC
			ST sub_2_ind_by_d2(const ST& ix, const ST& iy) const 
			{ 
				return iy*this->nx + ix;
			}

			CGPU_EXEC
			ST sub_2_ind_by_d2_pbc(const ST& ix, const ST& iy) const 
			{ 
				return ind_y_pbc(iy)*this->nx + this->ind_x_pbc(ix);
			}

			/***************************************************************************************/
			/****************** linear indices-> subindices using 2nd dimension ********************/
			/***************************************************************************************/
			CGPU_EXEC
			void ind_2_sub_by_d2(const ST& ind, ST& ix, ST& iy) const 
			{ 
				iy = ind/this->nx;
				ix = ind - iy*this->nx;
			}

			/***************************************************************************************/
			/******************************* Fourier space indices *********************************/
			/***************************************************************************************/
			CGPU_EXEC
			ST igy(const ST& iy) const 
			{ 
				return iy - ny_h;
			}

			CGPU_EXEC
			void igv(ST& ix, ST& iy) const 
			{ 
				ix = this->igx(ix);
				iy = igy(iy);
			}

			/************************************* shift *******************************************/
			CGPU_EXEC
			ST igy_sft(const ST& iy) const 
			{ 
				return (iy<ny_h)?iy:(iy-ny);
			}

			CGPU_EXEC
			void igv_sft(ST& ix, ST& iy) const 
			{ 
				ix = this->igx_sft(ix);
				iy = igy_sft(iy);
			}

			/***************************************************************************************/
			/******************************* real space indices ************************************/
			/***************************************************************************************/
			CGPU_EXEC
			ST iry(const ST& iy) const 
			{ 
				return iy;
			}

			CGPU_EXEC
			void irv(ST& ix, ST& iy) const 
			{ 
				ix = this->irx(ix);
				iy = iry(iy);
			}
			/***************************************** shift ***************************************/
			CGPU_EXEC
			ST iry_sft(const ST& iy) const 
			{ 
				return (iy<ny_h)?(iy+ny_h):(iy-ny_h);
			}

			CGPU_EXEC
			void irv_sft(ST& ix, ST& iy) const
			{ 
				ix = this->irx_sft(ix);
				iy = iry_sft(iy);
			}

		#ifdef __CUDACC__
			FCNS_GPU_GRID_BLK_IGRID_2D;
		#endif
		};
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
			iGrid_sxd(): iGrid_sxd<ST, edim_2>(), nz(0), nz_h(0) {}

			iGrid_sxd(const ST& nx, const ST& ny, const ST& nz): iGrid_sxd<ST, edim_2>(nx, ny), nz(nz), nz_h(nz/ST(2)) {}


			/* copy constructor */
			CGPU_EXEC
			iGrid_sxd(const iGrid_sxd<ST, edim_3>& grid): iGrid_sxd()
			{
				*this = grid;
			}

			/* converting constructor */
			template <class SU>
			CGPU_EXEC
			iGrid_sxd(const iGrid_sxd<SU, edim_3>& grid): iGrid_sxd()
			{
				*this = grid;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			iGrid_sxd<ST, edim_3>& operator=(const iGrid_sxd<ST, edim_3>& grid)
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
			template <class SU>
			CGPU_EXEC
			iGrid_sxd<ST, edim_3>& operator=(const iGrid_sxd<SU, edim_3>& grid)
			{
				iGrid_sxd<ST, edim_2>::operator=(grid);
				nz = ST(grid.nz);
				nz_h = ST(grid.nz_h);

				return *this;
			}

			template <class SU> 
			CGPU_EXEC
			void assign(const iGrid_sxd<SU, edim_3>& grid)
			{
				*this = grid;
			}

			/***************************************************************************************/
			void set_size(const ST& nx, const ST& ny, const ST& nz)
			{
				iGrid_sxd<ST, edim_2>::set_size(nx, ny);

				this->nz = nz;
				this->nz_h = nz/ST(2);
			}

			CGPU_EXEC
			void clear()
			{
				iGrid_sxd<ST, edim_2>::clear();
				nz = 0;
				nz_h = 0;
			}

			CGPU_EXEC
			virtual ST size() const 
			{ 
				return this->nx*this->ny*nz;
			}

			template<class SU>
			CGPU_EXEC
			SU nz_cast() const 
			{ 
				return SU(nz);
			}
			
			dt_shape_st<ST> shape() const
			{
				return {this->ny, this->nx, nz, ST(1)};
			}

			/***************************************************************************************/
			/************************* apply boundary conditions to indices ************************/
			/*nnhttps:// en.wikipedia.org/wiki/Periodic_boundary_conditionsnn */
			/***************************************************************************************/
			CGPU_EXEC
			ST ind_z_pbc(const ST& iz) const
			{
				return fcn_ind_pbc<ST>(iz, nz);
			}	

			CGPU_EXEC
			ST ind_z_bd(const ST& iz) const
			{
				return fcn_set_bound(iz, ST(0), nz-ST(1));
			}	

			CGPU_EXEC
			void ind_pbc(ST& ix, ST& iy, ST& iz) const
			{
				ix = ind_x_pbc(ix);
				iy = ind_y_pbc(iy);
				iz = ind_z_pbc(iz);
			}	

			/***************************************************************************************/
			/**************************** subindices -> linear indices *****************************/
			/***************************************************************************************/
			CGPU_EXEC
			ST operator()(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return iy + ix*this->ny + iz*this->ny*this->nx;
			}

			CGPU_EXEC
			ST sub_2_ind(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return iy + ix*this->ny + iz*this->ny*this->nx;
			}

			CGPU_EXEC
			ST sub_2_ind_bd(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return sub_2_ind(this->ind_x_bd(ix), this->ind_y_bd(iy), ind_z_bd(iz));
			}

			CGPU_EXEC
			ST sub_2_ind_pbc(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return this->ind_y_pbc(iy) + this->ind_x_pbc(ix)*this->ny + ind_z_pbc(iz)*this->ny*this->nx;
			}

			CGPU_EXEC
			ST sub_2_ind_irv_sft_pbc(ST ix, ST iy, ST iz) const
			{ 
				ind_pbc(ix, iy, iz);
				irv_sft(ix, iy, iz);

				return sub_2_ind(ix, iy, iz);
			}

			/***************************************************************************************/
			//**************************** linear indices-> subindices ******************************/
			/***************************************************************************************/
			CGPU_EXEC
			void ind_2_sub(const ST& ind, ST& ix, ST& iy, ST& iz) const 
			{ 
				iz = ind/(this->ny*this->nx);
				iy = ind - iz*this->ny*this->nx;
				ix = iy/this->ny;
				iy = iy - ix*this->ny;
			}

			/***************************************************************************************/
			/******************* subindices to linear indices using 2nd dimension ******************/
			/***************************************************************************************/
			CGPU_EXEC
			ST sub_2_ind_by_d2(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return iy*this->nx + ix + iz*this->ny*this->nx;
			}

			CGPU_EXEC
			ST sub_2_ind_by_d2_pbc(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return this->ind_y_pbc(iy)*this->nx + this->ind_x_pbc(ix) + ind_z_pbc(iz)*this->ny*this->nx;
			}

			/***************************************************************************************/
			/****************** linear indices-> subindices using 2nd dimension ********************/
			/***************************************************************************************/
			CGPU_EXEC
			void ind_2_sub_by_d2(const ST& ind, ST& ix, ST& iy, ST& iz) const 
			{ 
				iz = ind/(this->ny*this->nx);
				ix = ind - iz*this->ny*this->nx;
				iy = ix/this->nx;
				ix = ix - iy*this->nx;
			}

			/***************************************************************************************/
			/****************************** Fourier space indices **********************************/
			/***************************************************************************************/
			CGPU_EXEC
			ST igz(const ST& iz) const 
			{ 
				return iz - nz_h;
			}

			CGPU_EXEC
			void igv(ST& ix, ST& iy, ST& iz) const 
			{ 
				ix = this->igx(ix);
				iy = this->igy(iy);
				iz = igz(iz);
			}

			/************************************* shift *******************************************/
			CGPU_EXEC
			ST igz_sft(const ST& iz) const 
			{ 
				return (iz<nz_h)?iz:(iz-nz);
			}

			CGPU_EXEC
			void igv_sft(ST& ix, ST& iy, ST& iz) const 
			{ 
				ix = this->igx_sft(ix);
				iy = this->igy_sft(iy);
				iz = igz_sft(iz);
			}

			/***************************************************************************************/
			/******************************* real space indices ************************************/
			/***************************************************************************************/
			CGPU_EXEC
			ST irz(const ST& iz) const 
			{ 
				return iz;
			}

			CGPU_EXEC
			void irv(ST& ix, ST& iy, ST& iz) const 
			{ 
				ix = this->irx(ix);
				iy = this->iry(iy);
				iz = irz(iz);
			}

			/************************************* shift *******************************************/
			CGPU_EXEC
			ST irz_sft(const ST& iz) const 
			{ 
				return (iz<nz_h)?(iz+nz_h):(iz-nz_h);
			}

			CGPU_EXEC
			void irv_sft(ST& ix, ST& iy, ST& iz) const
			{ 
				ix = this->irx_sft(ix);
				iy = this->iry_sft(iy);
				iz = irz_sft(iz);
			}

		#ifdef __CUDACC__
			FCNS_GPU_GRID_BLK_IGRID_3D
		#endif
		};													
	}
#endif