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

#ifndef RANGE_H
	#define RANGE_H

	#ifdef _MSC_VER
		#pragma once
	#endif

	#include "const_enum.cuh"
	#include "r_2d.cuh"
	#include "r_3d.cuh"
	#include "grid.cuh"	

	/***************************************************************************************/
	/************************** Range template forward declaration *************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class ST, eDim Dim> class Range_xd;

		template <eDim Dim>
		using Range = Range_xd<dt_int32, Dim>;

		/* 1d */
		template<class ST>
		using Range_1d_st = Range_xd<ST, edim_1>;

		using Range_1d = Range_xd<dt_int32, edim_1>;

		using Range_1d_64 = Range_xd<dt_int64, edim_1>;

		/* 2d */
		template<class ST>
		using Range_2d_st = Range_xd<ST, edim_2>;

		using Range_2d = Range_xd<dt_int32, edim_2>;

		using Range_2d_64 = Range_xd<dt_int64, edim_2>;

		/* 3d */
		template<class ST>
		using Range_3d_st = Range_xd<ST, edim_3>;

		using Range_3d = Range_xd<dt_int32, edim_3>;

		using Range_3d_64 = Range_xd<dt_int64, edim_3>;	
	}
	
	/* template specialization 1d */
	namespace mt
	{
		template <class ST>
		class Range_xd<ST, edim_1>
		{
		public:
			using size_type = ST;

			ST ix_0;		// initial x index
			ST nx;			// x points

			/************************************* constructors ************************************/
			CGPU_EXEC
			Range_xd() : ix_0(0), nx(0) {}

			Range_xd(const ST& ix_0, const ST& nx) : ix_0(ix_0), nx(nx) {}

			/* copy constructor */
			CGPU_EXEC
			Range_xd(const Range_xd<ST, edim_1>& range)
			{
				*this = range;
			}

			/* converting constructor */
			template <class SU>
			CGPU_EXEC
			Range_xd(const Range_xd<SU, edim_1>& range)
			{
				*this = range;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			Range_xd<ST, edim_1>& operator=(const Range_xd<ST, edim_1>& range)
			{
				if (this != &range)
				{
					ix_0 = range.ix_0;
					nx = range.nx;
				}

				return *this;
			}

			/* converting assignment operator */
			template <class SU>
			CGPU_EXEC
			Range_xd<ST, edim_1>& operator=(const Range_xd<SU, edim_1>& range)
			{
				if (this != &range)
				{
					ix_0 = ST(range.ix_0);
					nx = ST(range.nx);
				}

				return *this;
			}

			template <class SU>
			CGPU_EXEC
				void assign(const Range_xd<SU, edim_1>& range)
			{
				*this = range;
			}

			/***************************************************************************************/
			template <class SU>
			void set_in_data(const SU& ix_0, const SU& nx)
			{
				this->ix_0 = ST(ix_0);
				this->nx = ST(nx);
			}

			template <class T>
			void set_in_data(const T& r, const T& r_max, const Grid_1d<T>& grid)
			{
				grid.ix_0_ix_n(r, r_max, ix_0, nx);
			}

			CGPU_EXEC
				void clear()
			{
				ix_0 = ST(0);
				nx = ST(0);
			}
		};
	}

	/* template specialization 2d */
	namespace mt
	{
		template <class ST>
		class Range_xd<ST, edim_2>: public Range_xd<ST, edim_1>
		{
		public:
			using size_type = ST;

			ST iy_0;	// initial y index
			ST ny;		// y points

			/************************************* constructors ************************************/
			CGPU_EXEC
			Range_xd(): Range_xd<ST, edim_1>(), iy_0(0), ny(0) {}

			Range_xd(const ST& ix_0, const ST& nx, const ST& iy_0, const ST& ny)
				: Range_xd<ST, edim_1>(ix_0, nx), iy_0(iy_0), ny(ny) {}

			/* copy constructor */
			CGPU_EXEC
			Range_xd(const Range_xd<ST, edim_2>& range)
			{
				*this = range;
			}

			/* converting constructor */
			template <class SU>
			CGPU_EXEC
			Range_xd(const Range_xd<SU, edim_2>& range)
			{
				*this = range;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			Range_xd<ST, edim_2>& operator=(const Range_xd<ST, edim_2>& range)
			{
				if (this != &range)
				{
					Range_xd<ST, edim_1>::operator=(range);

					iy_0 = range.iy_0;
					ny = range.ny;
				}

				return *this;
			}
			
			/* converting assignment operator */
			template <class SU>
			CGPU_EXEC
			Range_xd<ST, edim_2>& operator=(const Range_xd<SU, edim_2>& range)
			{
				if (this != &range)
				{
					Range_xd<ST, edim_1>::operator=(range);

					iy_0 = ST(range.iy_0);
					ny = ST(range.ny);
				}

				return *this;
			}

			template <class SU> 
			CGPU_EXEC
			void assign(const Range_xd<SU, edim_2>& range)
			{
				*this = range;
			}

			/***************************************************************************************/
			template <class SU> 
			void set_in_data(const SU& ix_0, const SU& nx, const SU& iy_0, const SU& ny)
			{
				Range_xd<ST, edim_1>::set_in_data(ix_0, nx);

				this->iy_0 = ST(iy_0);
				this->ny = ST(ny);
			}

			template <class T> 
			void set_in_data(const R_2d<T>& r, const T& r_max, const Grid_2d<T>& grid)
			{
				grid.ix_0_ix_n(r.x, r_max, this->ix_0, this->nx);
				grid.iy_0_iy_n(r.y, r_max, iy_0, ny);
			}

			CGPU_EXEC
			void clear()
			{
				Range_xd<ST, edim_1>::clear();

				iy_0 = ST(0);
				ny = ST(0);
			}
		};
	}

	/* template specialization 3d */
	namespace mt
	{
		template <class ST>
		class Range_xd<ST, edim_3>: public Range_xd<ST, edim_2>
		{
		public:
			using size_type = ST;

			ST iz_0;	// initial z index
			ST nz;		// z points

			/************************************* constructors ************************************/
			CGPU_EXEC
			Range_xd(): Range_xd<ST, edim_2>(), iz_0(0), nz(0) {}

			Range_xd(const ST& ix_0, const ST& nx, const ST& iy_0, const ST& ny, const ST& iz_0, const ST& nz)
				: Range_xd<ST, edim_2>(ix_0, nx, iy_0, ny), iz_0(iz_0), nz(nz) {}

			/* copy constructor */
			CGPU_EXEC
			Range_xd(const Range_xd<ST, edim_3>& range)
			{
				*this = range;
			}

			/* converting constructor */
			template <class SU>
			CGPU_EXEC
			Range_xd(const Range_xd<SU, edim_3>& range)
			{
				*this = range;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			Range_xd<ST, edim_3>& operator=(const Range_xd<ST, edim_3>& range)
			{
				if (this != &range)
				{
					Range_xd<ST, edim_2>::operator=(range);

					iz_0 = range.iz_0;
					nz = range.nz;
				}

				return *this;
			}
			
			/* converting assignment operator */
			template <class SU>
			CGPU_EXEC
			Range_xd<ST, edim_3>& operator=(const Range_xd<SU, edim_3>& range)
			{
				if (this != &range)
				{
					Range_xd<ST, edim_2>::operator=(range);

					iz_0 = ST(range.iz_0);
					nz = ST(range.nz);
				}

				return *this;
			}

			template <class SU> 
			CGPU_EXEC
			void assign(const Range_xd<SU, edim_3>& range)
			{
				*this = range;
			}

			/***************************************************************************************/
			template <class SU> 
			void set_in_data(const SU& ix_0, const SU& nx, const SU& iy_0, const SU& ny, const SU& iz_0, const SU& nz)
			{
				Range_xd<ST, edim_2>::set_in_data(ix_0, nx, iy_0, ny);

				this->iz_0 = ST(iz_0);
				this->nz = ST(nz);
			}

			template <class T> 
			void set_in_data(const R_3d<T>& r, const T& r_max, const Grid_3d<T>& grid)
			{
				grid.ix_0_ix_n(r.x, r_max, this->ix_0, this->nx);
				grid.iy_0_iy_n(r.y, r_max, this->iy_0, this->ny);
				grid.iz_0_iz_n(r.z, r_max, iz_0, nz);
			}

			CGPU_EXEC
			void clear()
			{
				Range_xd<ST, edim_2>::clear();

				iz_0 = ST(0);
				nz = ST(0);
			}
		};	
	}

#endif 