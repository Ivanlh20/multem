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

#ifndef ITHREAD_RECT_H
	#define ITHREAD_RECT_H

	#ifdef _MSC_VER
		#pragma once
	#endif

	#include "const_enum.cuh"

	/***************************************************************************************/
	/********************** iThread_Rect template forward declaration **********************/
	/***************************************************************************************/
	namespace mt
	{
		template <class ST, eDim Dim> class iThread_Rect_sxd;

		template <eDim Dim>
		using iThread_Rect_xd = iThread_Rect_sxd<dt_int32, Dim>;

		/* 1d */
		template<class ST>
		using iThread_Rect_1d_st = iThread_Rect_sxd<ST, edim_1>;

		using iThread_Rect_1d = iThread_Rect_sxd<dt_int32, edim_1>;

		using iThread_Rect_1d_64 = iThread_Rect_sxd<dt_int64, edim_1>;

		/* 2d */
		template<class ST>
		using iThread_Rect_2d_st = iThread_Rect_sxd<ST, edim_2>;

		using iThread_Rect_2d = iThread_Rect_sxd<dt_int32, edim_2>;

		using iThread_Rect_2d_64 = iThread_Rect_sxd<dt_int64, edim_2>;	

		/* 3d */
		template<class ST>
		using iThread_Rect_3d_st = iThread_Rect_sxd<ST, edim_3>;

		using iThread_Rect_3d = iThread_Rect_sxd<dt_int32, edim_3>;

		using iThread_Rect_3d_64 = iThread_Rect_sxd<dt_int64, edim_3>;	
	}

	/* template specialization 1d */
	namespace mt
	{
		template <class ST>
		class iThread_Rect_sxd<ST, edim_1>
		{
		public:
			using size_type = ST;

			ST ind_0;			// initial linear index
			ST ind_e;			// final linear index

			ST ix_0;			// x initial index
			ST ix_e;			// x final index

			dt_int32 istm;	// stream number

			/************************************* constructors ************************************/
			iThread_Rect_sxd(): ix_0(0), ix_e(0), ind_0(0), ind_e(0), istm(0) {}

			iThread_Rect_sxd(const dt_init_list_int32& list): istm(0)
			{	
				auto ptr = list.begin();

				ix_0 = ptr[0];
				ix_e = ptr[1];

				ind_0 = ST(0);
				ind_e = size();
			}

			iThread_Rect_sxd(const ST& ix_0, const ST& ix_e): ix_0(ix_0), ix_e(ix_e), istm(0)
			{
				ind_0 = ST(0);
				ind_e = size();
			}

			iThread_Rect_sxd(const ST& nx): ix_0(0), ix_e(nx), istm(0)
			{
				ind_0 = ST(0);
				ind_e = size();
			}

			iThread_Rect_sxd(const ST& ix_0, const ST& ix_e, const ST& ind_0, const ST& ind_e)
				: ix_0(ix_0), ix_e(ix_e), ind_0(ind_0), ind_e(ind_e), istm(0) {}

			/* copy constructor */
			CGPU_EXEC
			iThread_Rect_sxd(const iThread_Rect_sxd<ST, edim_1>& ithread)
			{
				*this = ithread;
			}

			/* converting constructor */
			template <class SU>
			CGPU_EXEC
			iThread_Rect_sxd(const iThread_Rect_sxd<SU, edim_1>& ithread)
			{
				*this = ithread;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			iThread_Rect_sxd<ST, edim_1>& operator=(const iThread_Rect_sxd<ST, edim_1>& ithread)
			{
				if (this != &ithread)
				{
					ix_0 = ithread.ix_0;
					ix_e = ithread.ix_e;

					ind_0 = ithread.ind_0;
					ind_e = ithread.ind_e;

					istm = ithread.istm;
				}

				return *this;
			}			
			
			/* converting assignment operator */
			template <class SU>
			CGPU_EXEC
			iThread_Rect_sxd<ST, edim_1>& operator=(const iThread_Rect_sxd<SU, edim_1>& ithread)
			{
				if (this != &ithread)
				{
					ix_0 = ST(ithread.ix_0);
					ix_e = ST(ithread.ix_e);

					ind_0 = ST(ithread.ind_0);
					ind_e = ST(ithread.ind_e);

					istm = ithread.istm;
				}

				return *this;
			}

			template <class SU> 
			CGPU_EXEC
			void assign(const iThread_Rect_sxd<SU, edim_1>& ithread)
			{
				*this = ithread;
			}

			/***************************************************************************************/
			CGPU_EXEC
			void clear()
			{
				ix_0 = ST(0);
				ix_e = ST(0);

				ind_0 = ST(0);
				ind_e = ST(0);

				istm = 0;
			}

			CGPU_EXEC
			ST nx() const 
			{ 
				return ix_e - ix_0;
			}

			CGPU_EXEC
			ST size() const 
			{ 
				return nx();
			}

			CGPU_EXEC
 			dt_bool chk_bound_ix(const ST& ix) const
			{
				return fcn_chk_bound(ix, ix_0, ix_e);
			}

			CGPU_EXEC
			dt_bool chk_bound(const ST& ix) const
			{
				return chk_bound_ix(ix);
			}

			CGPU_EXEC
			void apply_bound_ix(const ST& ix_0, const ST& ix_e)
			{
				this->ix_0 = fcn_set_bound(this->ix_0, ix_0, ix_e);
				this->ix_e = fcn_set_bound(this->ix_e, ix_0, ix_e);
			}

 			CGPU_EXEC
			ST sub_2_ind(const ST& ix) const 
			{ 
				return ix-ix_0;
			}

			CGPU_EXEC
			void set_ind(const ST& ind_0, const ST& ind_e)
			{ 
				this->ind_0 = ind_0;
				this->ind_e = ind_e;
			}

			void set_ind_asc_order()
			{ 
				if (ix_0>ix_e)
				{
					std::swap(ix_0, ix_e);
				}
			}

		};
	}

	/* template specialization 2d */
	namespace mt
	{
		template <class ST>
		class iThread_Rect_sxd<ST, edim_2>: public iThread_Rect_sxd<ST, edim_1>
		{
		public:
			using size_type = ST;

			ST iy_0;	// y initial index
			ST iy_e;	// y final index

			/************************************* constructors ************************************/
			iThread_Rect_sxd(): iThread_Rect_sxd<ST, edim_1>(), iy_0(0), iy_e(0) {}

			iThread_Rect_sxd(const dt_init_list_int32& list)
			{	
				auto ptr = list.begin();

				this->ix_0 = ptr[0];
				this->ix_e = ptr[1];
				this->iy_0 = ptr[2];
				this->iy_e = ptr[3];

				this->istm = 0;

				this->ind_0 = ST(0);
				this->ind_e = size();
			}

			iThread_Rect_sxd(const ST& ix_0, const ST& ix_e, const ST& iy_0, const ST& iy_e)
				: iThread_Rect_sxd<ST, edim_1>(ix_0, ix_e), iy_0(iy_0), iy_e(iy_e)
			{
				this->ind_0 = ST(0);
				this->ind_e = size();
			}
			
			iThread_Rect_sxd(const ST& nx, const ST& ny)
				: iThread_Rect_sxd<ST, edim_1>(nx), iy_0(0), iy_e(ny)
			{
				this->ind_0 = ST(0);
				this->ind_e = size();
			}

			iThread_Rect_sxd(const ST& ix_0, const ST& ix_e, const ST& iy_0, const ST& iy_e, const ST& ind_0, const ST& ind_e)
				: iThread_Rect_sxd<ST, edim_1>(ix_0, ix_e, ind_0, ind_e), iy_0(iy_0), iy_e(iy_e) {}

			/* copy constructor */
			CGPU_EXEC
			iThread_Rect_sxd(const iThread_Rect_sxd<ST, edim_2>& ithread)
			{
				*this = ithread;
			}

			/* converting constructor */
			template <class SU>
			CGPU_EXEC
			iThread_Rect_sxd(const iThread_Rect_sxd<SU, edim_2>& ithread)
			{
				*this = ithread;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			iThread_Rect_sxd<ST, edim_2>& operator=(const iThread_Rect_sxd<ST, edim_2>& ithread)
			{
				if (this != &ithread)
				{
					iThread_Rect_sxd<ST, edim_1>::operator=(ithread);

					iy_0 = ithread.iy_0;
					iy_e = ithread.iy_e;
				}

				return *this;
			}
			
			/* converting assignment operator */
			template <class SU>
			CGPU_EXEC
			iThread_Rect_sxd<ST, edim_2>& operator=(const iThread_Rect_sxd<SU, edim_2>& ithread)
			{
				if (this != &ithread)
				{
					iThread_Rect_sxd<ST, edim_1>::operator=(ithread);

					iy_0 = ST(ithread.iy_0);
					iy_e = ST(ithread.iy_e);
				}

				return *this;
			}

			template <class SU> 
			CGPU_EXEC
			void assign(const iThread_Rect_sxd<SU, edim_2>& ithread)
			{
				*this = ithread;
			}

			/***************************************************************************************/
			CGPU_EXEC
			void clear()
			{
				iThread_Rect_sxd<ST, edim_1>::clear();

				iy_0 = ST(0);
				iy_e = ST(0);
			}


			CGPU_EXEC
			ST ny() const 
			{ 
				return iy_e - iy_0;
			}

			CGPU_EXEC
			ST size() const 
			{ 
				return iThread_Rect_sxd<ST, edim_1>::size()*ny();
			}


			CGPU_EXEC
 			dt_bool chk_bound_iy(const ST& iy) const
			{
				return fcn_chk_bound(iy, iy_0, iy_e);
			}

			CGPU_EXEC
			dt_bool chk_bound(const ST& ix, const ST& iy) const
			{
				return this->chk_bound_ix(ix) && chk_bound_iy(iy);
			}

			CGPU_EXEC
			void apply_bound_iy(const ST& iy_0, const ST& iy_e)
			{
				this->iy_0 = fcn_set_bound(this->iy_0, iy_0, iy_e);
				this->iy_e = fcn_set_bound(this->iy_e, iy_0, iy_e);
			}

 			CGPU_EXEC
			ST sub_2_ind(const ST& ix, const ST& iy) const 
			{ 
				return (iy-iy_0) + (ix-this->ix_0)*ny();
			}

			void set_ind_asc_order()
			{ 
				iThread_Rect_sxd<ST, edim_1>::set_ind_asc_order();

				if (iy_0>iy_e)
				{
					std::swap(iy_0, iy_e);
				}
			}

		};
	}

	/* template specialization 3d */
	namespace mt
	{
		template <class ST>
		class iThread_Rect_sxd<ST, edim_3>: public iThread_Rect_sxd<ST, edim_2>
		{
		public:
			using size_type = ST;

			ST iz_0;	// z initial index
			ST iz_e;	// z final index

			/************************************* constructors ************************************/
			iThread_Rect_sxd(): iThread_Rect_sxd<ST, edim_2>(), iz_0(0), iz_e(0) {}

			iThread_Rect_sxd(const dt_init_list_int32& list)
			{	
				auto ptr = list.begin();

				this->ix_0 = ptr[0];
				this->ix_e = ptr[1];
				this->iy_0 = ptr[2];
				this->iy_e = ptr[3];
				this->iz_0 = ptr[4];
				this->iz_e = ptr[5];

				this->istm = 0;

				this->ind_0 = ST(0);
				this->ind_e = size();
			}
			iThread_Rect_sxd(const ST& ix_0, const ST& ix_e, const ST& iy_0, const ST& iy_e, const ST& iz_0, const ST& iz_e)
				: iThread_Rect_sxd<ST, edim_2>(ix_0, ix_e, iy_0, iy_e), iz_0(iz_0), iz_e(iz_e)
			{
				this->ind_0 = ST(0);
				this->ind_e = size();
			}

			iThread_Rect_sxd(const ST& nx, const ST& ny, const ST& nz)
				: iThread_Rect_sxd<ST, edim_2>(nx, ny), iz_0(0), iz_e(nz)
			{
				this->ind_0 = ST(0);
				this->ind_e = size();
			}

			iThread_Rect_sxd(const ST& ix_0, const ST& ix_e, const ST& iy_0, const ST& iy_e, const ST& iz_0, const ST& iz_e, const ST& ind_0, const ST& ind_e)
				: iThread_Rect_sxd<ST, edim_2>(ix_0, ix_e, iy_0, iy_e, ind_0, ind_e), iz_0(iz_0), iz_e(iz_e) {}

			/* copy constructor */
			CGPU_EXEC
			iThread_Rect_sxd(const iThread_Rect_sxd<ST, edim_3>& ithread)
			{
				*this = ithread;
			}

			/* converting constructor */
			template <class SU>
			CGPU_EXEC
			iThread_Rect_sxd(const iThread_Rect_sxd<SU, edim_3>& ithread)
			{
				*this = ithread;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			iThread_Rect_sxd<ST, edim_3>& operator=(const iThread_Rect_sxd<ST, edim_3>& ithread)
			{
				if (this != &ithread)
				{
					iThread_Rect_sxd<ST, edim_2>::operator=(ithread);

					iz_0 = ithread.iz_0;
					iz_e = ithread.iz_e;
				}

				return *this;
			}
			
			/* converting assignment operator */
			template <class SU>
			CGPU_EXEC
			iThread_Rect_sxd<ST, edim_3>& operator=(const iThread_Rect_sxd<SU, edim_3>& ithread)
			{
				if (this != &ithread)
				{
					iThread_Rect_sxd<ST, edim_2>::operator=(ithread);

					iz_0 = ST(ithread.iz_0);
					iz_e = ST(ithread.iz_e);
				}

				return *this;
			}

			template <class SU> 
			CGPU_EXEC
			void assign(const iThread_Rect_sxd<SU, edim_3>& ithread)
			{
				*this = ithread;
			}

			/***************************************************************************************/
			CGPU_EXEC
			void clear()
			{
				iThread_Rect_sxd<ST, edim_2>::clear();

				iz_0 = ST(0);
				iz_e = ST(0);
			}

			CGPU_EXEC
			ST nz() const 
			{ 
				return iz_e - iz_0;
			}

			CGPU_EXEC
			ST size() const 
			{ 
				return iThread_Rect_sxd<ST, edim_2>::size()*nz();
			}

			CGPU_EXEC
 			dt_bool chk_bound_iz(const ST& iz) const
			{
				return fcn_chk_bound(iz, iz_0, iz_e);
			}

			CGPU_EXEC
			dt_bool chk_bound(const ST& ix, const ST& iy, const ST& iz) const
			{
				return this->chk_bound_ix(ix) && this->chk_bound_iy(iy) && chk_bound_iz(iz);
			}

			CGPU_EXEC
			void apply_bound_iz(const ST& iz_0, const ST& iz_e)
			{
				this->iz_0 = fcn_set_bound(this->iz_0, iz_0, iz_e);
				this->iz_e = fcn_set_bound(this->iz_e, iz_0, iz_e);
			}

 			CGPU_EXEC
			ST sub_2_ind(const ST& ix, const ST& iy, const ST& iz) const 
			{ 
				return (iy-this->iy_0) + (ix-this->ix_0)*this->ny_s() + (iz-iz_0)*this->ny()*this->nx();
			}

			void set_ind_asc_order()
			{ 
				iThread_Rect_sxd<ST, edim_2>::set_ind_asc_order();

				if (iz_0>iz_e)
				{
					std::swap(iz_0, iz_e);
				}
			}

		};
	}

#endif