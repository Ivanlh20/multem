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

#ifndef BORDER_H
	#define BORDER_H

	#ifdef _MSC_VER
		#pragma once
	#endif

	#include "math.cuh"
	#include "const_enum.cuh"
	#include "cgpu_fcns_gen.cuh"
	#include "r_2d.cuh"
	#include "r_3d.cuh"
	#include "cgpu_vctr.cuh"

	/***************************************************************************************/
	/******************************** rectangular border ***********************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T, eDim Dim> class Border_Rect_xd;

		template <eDim Dim>
		using iBorder_Rect_xd = Border_Rect_xd<dt_int32, Dim>;

		/* 1d */
		template<class T>
		using Border_Rect_1d = Border_Rect_xd<T, edim_1>;

		using iBorder_Rect_1d = Border_Rect_xd<dt_int32, edim_1>;

		using iBorder_Rect_1d_64 = Border_Rect_xd<dt_int64, edim_1>;

		/* 2d */
		template<class T>
		using Border_Rect_2d = Border_Rect_xd<T, edim_2>;

		using iBorder_Rect_2d = Border_Rect_xd<dt_int32, edim_2>;

		using iBorder_Rect_2d_64 = Border_Rect_xd<dt_int64, edim_2>;

		/* 3d */
		template<class T>
		using Border_Rect_3d = Border_Rect_xd<T, edim_3>;

		using iBorder_Rect_3d = Border_Rect_xd<dt_int32, edim_3>;

		using iBorder_Rect_3d_64 = Border_Rect_xd<dt_int64, edim_3>;
	}

	/* template specialization 1d */
	namespace mt
	{
		template <class T>
		class Border_Rect_xd<T, edim_1>
		{
		public:
			using value_type = T;
			using size_type = dt_int32;

			T bx_0;		// initial x position
			T bx_e;		// final x position

			/************************************* constructors ************************************/
			CGPU_EXEC
			Border_Rect_xd(): bx_0(0), bx_e(0) {}

			Border_Rect_xd(const T& bx_0, const T& bx_e)
			{
				set_in_data(bx_0, bx_e);
			}

			template <class U>
			Border_Rect_xd(const dt_init_list<U>& list)
			{
				set_in_data(list);
			}

			template <class U>
			Border_Rect_xd(const Vctr_cpu<U>& vctr)
			{
				set_in_data(vctr);
			}

			/* copy constructor */
			CGPU_EXEC
			Border_Rect_xd(const Border_Rect_xd<T, edim_1>& border)
			{
				*this = border;
			}

			/* converting constructor */
			template <class U>
			CGPU_EXEC
			Border_Rect_xd(const Border_Rect_xd<U, edim_1>& border)
			{
				*this = border;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			Border_Rect_xd<T, edim_1>& operator=(const Border_Rect_xd<T, edim_1>& border)
			{
				if (this != &border)
				{
					bx_0 = border.bx_0;
					bx_e = border.bx_e;
				}

				return *this;
			}			
			
			/* converting assignment operator */
			template <class U>
			CGPU_EXEC
			Border_Rect_xd<T, edim_1>& operator=(const Border_Rect_xd<U, edim_1>& border)
			{
				if (this != &border)
				{
					bx_0 = T(border.bx_0);
					bx_e = T(border.bx_e);
				}

				return *this;
			}

			template <class U> 
			CGPU_EXEC
			void assign(const Border_Rect_xd<U, edim_1>& border)
			{
				*this = border;
			}

			/***************************************************************************************/
			template <class U>
			void set_in_data(const U& bx_0, const U& bx_e)
			{
				this->bx_0 = bx_0;
				this->bx_e = bx_e;
			}

			template <class U> 
			void set_in_data(const dt_init_list<U>& list)
			{
				auto ptr = list.begin();

				bx_0 = T(ptr[0]);
				bx_e = T(ptr[1]);
			}

			template <class U>
			void set_in_data(const Vctr_cpu<U>& vctr)
			{
				if (vctr.size_32()==1)
				{
					set_in_data({vctr[0], vctr[0]});
				}
				else
				{
					set_in_data({vctr[0], vctr[1]});
				}
			}

			CGPU_EXEC
			void clear()
			{
				bx_0 = T(0);
				bx_e = T(0);
			}
			
			CGPU_EXEC
			T bx_sum() const
			{
				return bx_0 + bx_e;
			}

			CGPU_EXEC
			T bx_min() const
			{
				return fcn_min(bx_0, bx_e);
			}

			CGPU_EXEC
			T bx_max() const
			{
				return fcn_max(bx_0, bx_e);
			}
			CGPU_EXEC
			T bs_x(const T& bs_x_i) const 
			{ 
				return fcn_max(bs_x_i - bx_sum(), T(0));
			}

			CGPU_EXEC
			T bs(const T& bs_i) const 
			{ 
				return bs_x(bs_i);
			}

			CGPU_EXEC
			T bs_min(const T& bs) const 
			{ 
				return bs(bs);
			}

			CGPU_EXEC
			T bs_max(const T& bs) const 
			{ 
				return bs(bs);
			}

			CGPU_EXEC
			T bs_x_h(const T& bs_x) const 
			{ 
				return this->bs_x(bs_x)/T(2);
			}

			CGPU_EXEC
			T bs_h(const T& bs) const 
			{ 
				return bs_x_h(bs);
			}

			CGPU_EXEC
			T rx_c(const T& bs_x) const 
			{ 
				return bx_0 + bs_x_h(bs_x);
			}

			CGPU_EXEC
			T r_c(const T& bs) const 
			{ 
				return rx_c(bs);
			}

			CGPU_EXEC
			T radius_x(const T& bs_x) const
			{
				return bs_x_h(bs_x);
			}

			CGPU_EXEC
			T radius(const T& bs) const
			{
				return radius_x(bs);
			}

			CGPU_EXEC
			T radius_x_p(const T& bs_x, const T& p) const
			{
				return (T(1)-p)*radius_x(bs_x);
			}

			CGPU_EXEC
			T radius_p(const T& bs, const T& p) const
			{
				return radius_x_p(bs, p);
			}

			void set_by_sft(const T& bs, const T& dr)
			{
				bx_0 = fcn_set_bound(dr, T(0), bs);
				bx_e = fcn_set_bound(bs - dr, T(0), bs);
			}

			void sft_bdr(const T& bs, const T& dr)
			{
				bx_0 = fcn_set_bound(bx_0 + dr, T(0), bs);
				bx_e = fcn_set_bound(bx_e - dr, T(0), bs);
			}
		};
	}

	/* template specialization 2d */
	namespace mt
	{
		template <class T>
		class Border_Rect_xd<T, edim_2>: public Border_Rect_xd<T, edim_1>
		{		
		public:			
			using value_type = T;
			using size_type = dt_int32;

			T by_0;		// initial y position
			T by_e;		// final y position

			/************************************* constructors ************************************/
			CGPU_EXEC
			Border_Rect_xd(): Border_Rect_xd<T, edim_1>(), by_0(0), by_e(0) {}

			Border_Rect_xd(const T& bx_0, const T& bx_e, const T& by_0, const T& by_e)
			{
				set_in_data(bx_0, bx_e, by_0, by_e);
			}

			template <class U> 
			Border_Rect_xd(const dt_init_list<U>& list)
			{
				set_in_data(list);
			}

			template <class U> 
			Border_Rect_xd(const Vctr_cpu<U>& vctr)
			{
				set_in_data(vctr);
			}

			/* copy constructor */
			CGPU_EXEC
			Border_Rect_xd(const Border_Rect_xd<T, edim_2>& border)
			{
				*this = border;
			}

			/* converting constructor */
			template <class U>
			CGPU_EXEC
			Border_Rect_xd(const Border_Rect_xd<U, edim_2>& border)
			{
				*this = border;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			Border_Rect_xd<T, edim_2>& operator=(const Border_Rect_xd<T, edim_2>& border)
			{
				if (this != &border)
				{
					Border_Rect_xd<T, edim_1>::operator=(border);

					by_0 = border.by_0;
					by_e = border.by_e;
				}

				return *this;
			}			
			
			/* converting assignment operator */
			template <class U>
			CGPU_EXEC
			Border_Rect_xd<T, edim_2>& operator=(const Border_Rect_xd<U, edim_2>& border)
			{
				if (this != &border)
				{
					Border_Rect_xd<T, edim_1>::operator=(border);

					by_0 = T(border.by_0);
					by_e = T(border.by_e);
				}

				return *this;
			}

			template <class U> 
			CGPU_EXEC
			void assign(const Border_Rect_xd<U, edim_2>& border)
			{
				*this = border;
			}

			/***************************************************************************************/
			template <class U>
			void set_in_data(const U& bx_0, const U& bx_e, const U& by_0, const U& by_e)
			{
				Border_Rect_xd<T, edim_1>::set_in_data(bx_0, bx_e);

				this->by_0 = by_0;
				this->by_e = by_e;
			}

			template <class U> 
			void set_in_data(const dt_init_list<U>& list)
			{
				auto ptr = list.begin();

				Border_Rect_xd<T, edim_1>::set_in_data({ptr[0], ptr[1]});

				by_0 = T(ptr[2]);
				by_e = T(ptr[3]);
			}

			template <class U>
			void set_in_data(const Vctr_cpu<U>& vctr)
			{
				const auto shape = vctr.shape();

				if (shape[0]==1)
				{
					if (shape[1]<2)
						set_in_data({vctr[0], vctr[0], vctr[0], vctr[0]});
					else if (shape[1]<3)
						set_in_data({vctr[0], vctr[0], vctr[1], vctr[1]});
					else if (shape[1]<4)
						set_in_data({vctr[0], vctr[1], vctr[2], T(0)});			
					else
						set_in_data({vctr[0], vctr[1], vctr[2], vctr[3]});	
				}
				else if (shape[0]==2)
				{
					if (shape[1]<2)
						set_in_data({vctr[0], vctr[1], T(0), T(0)});
					else
						set_in_data({vctr[0], vctr[1], vctr[2], vctr[3]});
				}
			}

			CGPU_EXEC
			void clear()
			{
				Border_Rect_xd<T, edim_1>::clear();

				by_0 = T(0);
				by_e = T(0);
			}
	
			CGPU_EXEC
			T by_sum() const
			{
				return by_0 + by_e;
			}

			CGPU_EXEC
			T by_min() const
			{
				return fcn_min(by_0, by_e);
			}

			CGPU_EXEC
			T by_max() const
			{
				return fcn_max(by_0, by_e);
			}
			CGPU_EXEC
			T bs_y(const T& bs_y_i) const 
			{ 
				return fcn_max(bs_y_i - by_sum(), T(0));
			}

			CGPU_EXEC
			R_2d<T> bs(const R_2d<T>& bs_i) const 
			{ 
				return {this->bs_x(bs_i.x), bs_y(bs_i.y)};
			}

			CGPU_EXEC
			T bs_min(const R_2d<T>& bs) const 
			{ 
				return fmin(bs(bs));
			}

			CGPU_EXEC
			T bs_max(const R_2d<T>& bs) const 
			{ 
				return fmax(bs(bs));
			}

			CGPU_EXEC
			T bs_y_h(const T& bs_y) const 
			{ 
				return this->bs_y(bs_y)/T(2);
			}

			CGPU_EXEC
			R_2d<T> bs_h(const R_2d<T>& bs) const 
			{ 
				return {this->bs_x_h(bs.x), bs_y_h(bs.y)};
			}

			CGPU_EXEC
			T ry_c(const T& bs_y) const 
			{ 
				return by_0 + bs_y_h(bs_y);
			}

			CGPU_EXEC
			R_2d<T> r_c(const R_2d<T>& bs) const 
			{ 
				return {this->rx_c(bs.x), ry_c(bs.y)};
			}

			CGPU_EXEC
			T radius_y(const T& bs_y) const
			{
				return bs_y_h(bs_y);
			}

			CGPU_EXEC
			R_2d<T> radius(const R_2d<T>& bs) const
			{
				return {this->radius_x(bs.x), radius_y(bs.y)};
			}

			CGPU_EXEC
			T radius_y_p(const T& bs_y, const T& p) const
			{
				return (T(1)-p)*radius_y(bs_y);
			}

			CGPU_EXEC
			R_2d<T> radius_p(const R_2d<T>& bs, const T& p) const
			{
				return {this->radius_x_p(bs.x, p), radius_y_p(bs.y, p)};
			}

			void set_by_sft(const R_2d<T>& bs, const R_2d<T>& dr)
			{
				Border_Rect_xd<T, edim_1>::set_by_sft(bs.x, dr.x);

				by_0 = fcn_set_bound(dr.y, T(0), bs.y);
				by_e = fcn_set_bound(bs.y - dr.y, T(0), bs.y);
			}

			void sft_bdr(const R_2d<T>& bs, const R_2d<T>& dr)
			{
				Border_Rect_xd<T, edim_1>::sft_bdr(bs.x, dr.x);

				by_0 = fcn_set_bound(by_0 + dr.y, T(0), bs.y);
				by_e = fcn_set_bound(by_e - dr.y, T(0), bs.y);
			}
		};
	}

	/* template specialization 3d */
	namespace mt
	{
		template <class T>
		class Border_Rect_xd<T, edim_3>: public Border_Rect_xd<T, edim_2>
		{		
		public:			
			using value_type = T;
			using size_type = dt_int32;

			T bz_0;		// initial z position
			T bz_e;		// final z position

			/************************************* constructors ************************************/
			CGPU_EXEC
			Border_Rect_xd(): Border_Rect_xd<T, edim_2>(), bz_0(0), bz_e(0) {}

			Border_Rect_xd(const T& bx_0, const T& bx_e, const T& by_0, const T& by_e, const T& bz_0, const T& bz_e)
			{
				set_in_data(bx_0, bx_e, by_0, by_e, bz_0, bz_e);
			}

			template <class U> 
			Border_Rect_xd(const dt_init_list<U>& list)
			{
				set_in_data(list);
			}

			template <class U> 
			Border_Rect_xd(const Vctr_cpu<U>& vctr)
			{
				set_in_data(vctr);
			}

			/* copy constructor */
			CGPU_EXEC
			Border_Rect_xd(const Border_Rect_xd<T, edim_3>& border)
			{
				*this = border;
			}

			/* converting constructor */
			template <class U>
			CGPU_EXEC
			Border_Rect_xd(const Border_Rect_xd<U, edim_3>& border)
			{
				*this = border;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			Border_Rect_xd<T, edim_3>& operator=(const Border_Rect_xd<T, edim_3>& border)
			{
				if (this != &border)
				{
					Border_Rect_xd<T, edim_2>::operator=(border);

					bz_0 = border.bz_0;
					bz_e = border.bz_e;
				}

				return *this;
			}			
			
			/* converting assignment operator */
			template <class U>
			CGPU_EXEC
			Border_Rect_xd<T, edim_3>& operator=(const Border_Rect_xd<U, edim_3>& border)
			{
				if (this != &border)
				{
					Border_Rect_xd<T, edim_2>::operator=(border);

					bz_0 = T(border.bz_0);
					bz_e = T(border.bz_e);
				}

				return *this;
			}

			template <class U> 
			CGPU_EXEC
			void assign(const Border_Rect_xd<U, edim_3>& border)
			{
				*this = border;
			}

			/***************************************************************************************/
			template <class U>
			void set_in_data(const U& bx_0, const U& bx_e, const U& by_0, const U& by_e, const U& bz_0, const U& bz_e)
			{
				Border_Rect_xd<T, edim_2>::set_in_data(bx_0, bx_e, by_0, by_e);

				this->bz_0 = bz_0;
				this->bz_e = bz_e;
			}

			template <class U> 
			void set_in_data(const dt_init_list<U>& list)
			{
				auto ptr = list.begin();

				Border_Rect_xd<T, edim_2>::set_in_data({ptr[0], ptr[1], ptr[2], ptr[3]});

				bz_0 = T(ptr[4]);
				bz_e = T(ptr[5]);
			}

			template <class U>
			void set_in_data(const Vctr_cpu<U>& vctr)
			{
				const auto shape = vctr.shape();

				if (shape[0]==1)
				{
					if (shape[1]<2)
						set_in_data({vctr[0], vctr[0], vctr[0], vctr[0], T(0), T(0)});
					else if (shape[1]<3)
						set_in_data({vctr[0], vctr[0], vctr[1], vctr[1], T(0), T(0)});
					else if (shape[1]<4)
						set_in_data({vctr[0], vctr[0], vctr[1], vctr[1], vctr[2], vctr[2]});
					else if (shape[1]<5)
						set_in_data({vctr[0], vctr[1], vctr[2], vctr[3], T(0), T(0)});		
					else if (shape[1]<6)
						set_in_data({vctr[0], vctr[1], vctr[2], vctr[3], vctr[4], T(0)});			
					else
						set_in_data({vctr[0], vctr[1], vctr[2], vctr[3], vctr[4], vctr[5]});
				}
				else if (shape[0]==2)
				{
					if (shape[1]<2)
						set_in_data({vctr[0], vctr[1], T(0), T(0), T(0), T(0)});
					else if (shape[1]<3)
						set_in_data({vctr[0], vctr[1], vctr[2], vctr[3], T(0), T(0)});
					else
						set_in_data({vctr[0], vctr[1], vctr[2], vctr[3], vctr[4], vctr[5]});
				}
			}

			CGPU_EXEC
			void clear()
			{
				Border_Rect_xd<T, edim_2>::clear();

				bz_0 = T(0);
				bz_e = T(0);
			}

			CGPU_EXEC
			T bz_sum() const
			{
				return bz_0 + bz_e;
			}

			CGPU_EXEC
			T bz_min() const
			{
				return fcn_min(bz_0, bz_e);
			}

			CGPU_EXEC
			T bz_max() const
			{
				return fcn_max(bz_0, bz_e);
			}
			CGPU_EXEC
			T bs_z(const T& bs_z_i) const 
			{ 
				return fcn_max(bs_z_i - bz_sum(), T(0));
			}

			CGPU_EXEC
			R_3d<T> bs(const R_3d<T>& bs_i) const 
			{ 
				return {this->bs_x(bs_i.x), this->bs_y(bs_i.y), bs_z(bs_i.z)};
			}

			CGPU_EXEC
			T bs_min(const R_3d<T>& bs) const 
			{ 
				return fmin(bs(bs));
			}

			CGPU_EXEC
			T bs_max(const R_3d<T>& bs) const 
			{ 
				return fmax(bs(bs));
			}

			CGPU_EXEC
			T bs_z_h(const T& bs_z) const 
			{ 
				return this->bs_z(bs_z)/T(2);
			}

			CGPU_EXEC
			R_3d<T> bs_h(const R_3d<T>& bs) const 
			{ 
				return {this->bs_x_h(bs.x), this->bs_y_h(bs.y), bs_z_h(bs.z)};
			}

			CGPU_EXEC
			T rz_c(const T& bs_z) const 
			{ 
				return bz_0 + bs_z_h(bs_z);
			}

			CGPU_EXEC
			R_3d<T> r_c(const R_3d<T>& bs) const 
			{ 
				return {this->rx_c(bs.x), this->ry_c(bs.y), rz_c(bs.z)};
			}

			CGPU_EXEC
			T radius_z(const T& bs_z) const
			{
				return bs_z_h(bs_z);
			}

			CGPU_EXEC
			R_3d<T> radius(const R_3d<T>& bs) const
			{
				return {this->radius_x(bs.x), this->radius_y(bs.y), radius_z(bs.z)};
			}

			CGPU_EXEC
			T radius_z_p(const T& bs_z, const T& p) const
			{
				return (T(1)-p)*radius_z(bs_z);
			}

			CGPU_EXEC
			R_3d<T> radius_p(const R_3d<T>& bs, const T& p) const
			{
				return {this->radius_x_p(bs.x, p), this->radius_y_p(bs.y, p), radius_z_p(bs.z, p)};
			}

			void set_by_sft(const R_3d<T>& bs, const R_3d<T>& dr)
			{
				Border_Rect_xd<T, edim_2>::set_bz_sft({bs.x, bs.y}, {dr.x, dr.y});

				bz_0 = fcn_set_bound(dr.z, T(0), bs.z);
				bz_e = fcn_set_bound(bs.z - dr.z, T(0), bs.z);
			}

			void sft_bdr(const R_3d<T>& bs, const R_3d<T>& dr)
			{
				Border_Rect_xd<T, edim_2>::sft_bdr({bs.x, bs.y}, {dr.x, dr.y});

				bz_0 = fcn_set_bound(bz_0 + dr.z, T(0), bs.z);
				bz_e = fcn_set_bound(bz_e - dr.z, T(0), bs.z);
			}
		};
	}

#endif