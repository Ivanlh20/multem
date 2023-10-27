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

#include "r_3d.h"
#include "region_2d.h"

/* derived class */
namespace mt
{
	template <class T>
	using Region_Rect_3d = Region_Rect_xd<T, edim_3>;

	using iRegion_Rect_3d = Region_Rect_xd<dt_int32, edim_3>;

	using iRegion_Rect_3d_64 = Region_Rect_xd<dt_int64, edim_3>;
}

/* template specialization 3d */
namespace mt
{
	template <class T>
	class Region_Rect_xd<T, edim_3>: public Region_Rect_xd<T, edim_2>
	{		
	public:			
		using value_type = T;
		using size_type = dt_int32;

		T rz_0;		// initial z position
		T rz_e;		// final z position

		/************************************* constructors ************************************/
		CGPU_EXEC
		Region_Rect_xd();

		Region_Rect_xd(const T& rx_0, const T& rx_e, const T& ry_0, const T& ry_e, const T& rz_0, const T& rz_e);

		template <class U>
		Region_Rect_xd(U *data, const size_type& n_data);

		/* copy constructor */
		CGPU_EXEC
		Region_Rect_xd(const Region_Rect_xd<T, edim_3>& region);

		/* converting constructor */
		template <class U>
		CGPU_EXEC
		Region_Rect_xd(const Region_Rect_xd<U, edim_3>& region);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		Region_Rect_xd<T, edim_3>& operator=(const Region_Rect_xd<T, edim_3>& region);		
		
		/* converting assignment operator */
		template <class U>
		CGPU_EXEC
		Region_Rect_xd<T, edim_3>& operator=(const Region_Rect_xd<U, edim_3>& region);

		template <class U> 
		CGPU_EXEC
		void assign(const Region_Rect_xd<U, edim_3>& region);

		/***************************************************************************************/
		template <class U>
		void set_in_data(const U& rx_0, const U& rx_e, const U& ry_0, const U& ry_e, const U& rz_0, const U& rz_e);

		CGPU_EXEC
		void clear();

		CGPU_EXEC
		T bs_z() const;

		CGPU_EXEC
		R_3d<T> bs() const;

		CGPU_EXEC
		T bs_z_h() const;

		CGPU_EXEC
		R_3d<T> bs_h() const;

		CGPU_EXEC
		T rz_c() const;
		
		CGPU_EXEC
		R_3d<T> r_c() const;

		CGPU_EXEC
		T radius_z() const;

		CGPU_EXEC
		R_3d<T> radius() const;

		CGPU_EXEC
		T radius_z_p(const T& p) const;

		CGPU_EXEC
		R_3d<T> radius_p(const T& p) const;

		void sft_region(const R_3d<T>& bs, const R_3d<T>& dr);

		/***************************************************************************************/
		CGPU_EXEC
		dt_bool chk_bound_z(const T& z) const;

		CGPU_EXEC
		dt_bool chk_bound(const R_3d<T>& r) const;

		CGPU_EXEC
		dt_bool chk_bound(const T& rx, const T& ry, const T& rz) const;

		CGPU_EXEC
		dt_bool chk_bound_z_eps(const T& z) const;

		CGPU_EXEC
		dt_bool chk_bound_eps(const R_3d<T>& r) const;

		CGPU_EXEC
		dt_bool chk_bound_eps(const T& rx, const T& ry, const T& rz) const;

		template <class U=T, class = enable_if_int<U>>
		CGPU_EXEC
		T sub_2_ind(const T& ix, const T& iy, const T& iz) const;
	};
}

/***************************************************************************************/
/* symmetric region */
namespace mt
{
	/*
	template <class TVctr>
		struct Region_Rad_1d
		{
		public:
			using T = class TVctr::value_type;
			using size_type = dt_uint64;

			TVctr Rx;
			TVctr R2;
			TVctr Ixy;

			T R_max;
			T Rx_sf;

			T Ixy_sf;
			T Ixy_sc;
			T Ixy_sum;

			Region_Rad_1d(): Rx_sf(0), Ixy_sf(0), Ixy_sc(1), Ixy_sum(0)
			{
			}

			size_type size() const
			{
				return Ixy.size();
			}

			void clear()
			{
				Rx.clear();
				R2.clear();
				Ixy.clear();
			}

			void reserve(const size_type& new_size)
			{
				Rx.reserve(new_size);
				R2.reserve(new_size);
				Ixy.reserve(new_size);
			}

			void shrink_to_fit()
			{
				Rx.shrink_to_fit();
				R2.shrink_to_fit();
				Ixy.shrink_to_fit();
			}

			TVctr sft_Ixy(T bg)
			{
				TVctr Ixy_s;
				Ixy_s.reserve(Ixy.size());

				for(auto ixy=0; ixy<Ixy.size(); ixy++)
				{
					Ixy_s.push_back(Ixy[ixy]-bg);
				}
				return Ixy_s;
			}

			T sft_Rx(T x) const 
			{ 
				return (x-Rx_sf)/Rxy_sc;
			}
		};

		template <class TVctr>
		struct Region_Rad_2d
		{
		public:
			using T = class TVctr::value_type;
			using size_type = dt_uint64;

			TVctr Rx;
			TVctr Ry;
			TVctr R2;
			TVctr Ixy;

			T R_max;
			T Rx_sf;
			T Ry_sf;
			T Rxy_sc;

			T Ixy_sf;
			T Ixy_sc;
			T Ixy_sum;

			Region_Rad_2d(): Rx_sf(0), Ry_sf(0), Rxy_sc(1), 
			Ixy_sf(0), Ixy_sc(1), Ixy_sum(0) {}

			size_type size() const
			{
				return Ixy.size();
			}

			void clear()
			{
				Rx.clear();
				Ry.clear();
				R2.clear();
				Ixy.clear();
			}

			void reserve(const size_type& new_size)
			{
				Rx.reserve(new_size);
				Ry.reserve(new_size);
				R2.reserve(new_size);
				Ixy.reserve(new_size);
			}

			void shrink_to_fit()
			{
				Rx.shrink_to_fit();
				Ry.shrink_to_fit();
				R2.shrink_to_fit();
				Ixy.shrink_to_fit();
			}

			TVctr sft_Ixy(T bg)
			{
				TVctr Ixy_s;
				Ixy_s.reserve(Ixy.size());

				for(auto ixy=0; ixy<Ixy.size(); ixy++)
				{
					Ixy_s.push_back(Ixy[ixy]-bg);
				}
				return Ixy_s;
			}

			TVctr sub_region_to_Ixy(Grid_2d<T>& grid_2d, TVctr& Im_s, T x, T y)
			{
				TVctr v = Ixy;

				T R2_max = ::square(R_max);

				R_2d<T> p(x, y);
				auto range = grid_2d.ithread(p, R_max);
				dt_int32 iv = 0;
				for(auto ix = range.ix_0; ix < range.ix_e; ix++)
				{
					for(auto iy = range.iy_0; iy < range.iy_e; iy++)
					{
						T r2 = grid_2d.r2(ix, iy, p.x, p.y);
						if (r2 < R2_max)
						{
							v[iv++] -= Im_s[grid_2d.sub_2_ind(ix, iy)]/Ixy_sc;
						}
					}
				}

				return v;
			}

			TVctr sft_Ixy(Grid_2d<T>& grid_2d, TVctr& Im_s, T x, T y, T a, T s)
			{
				TVctr v = sub_region_to_Ixy(grid_2d, Im_s, x*Rxy_sc+Rx_sf, y*Rxy_sc+Ry_sf);

				T alpha = T(0.5)/::square(s);
				T r2_l = ::square(T(4)*s);
				for(auto im = 0; im < v.size(); im++)
				{
					T rx = Rx[im]-x;
					T ry = Ry[im]-y;
					T r2 = rx*rx+ry*ry;
					if (r2<r2_l)
					{
						v[im] += a*exp(-alpha*r2);
					}
				}

				return v;
			}

			TVctr sft_Ixy(Grid_2d<T>& grid_2d, TVctr& Im_s, TVctr& x, TVctr& y, TVctr& A, TVctr& S)
			{
				TVctr v = sub_region_to_Ixy(grid_2d, Im_s, x[0]*Rxy_sc+Rx_sf, y[0]*Rxy_sc+Ry_sf);

				for(auto ip = 0; ip < x.size(); ip++)
				{
					T a = A[ip];
					T b = S[ip];
					T alpha = 0.5/pow(b, 2);
					T r2_l = pow(4.0*b, 2);
					for(auto im = 0; im < v.size(); im++)
					{
						T rx = Rx[im]-x[ip];
						T ry = Ry[im]-y[ip];
						T r2 = rx*rx+ry*ry;
						if (r2<r2_l)
						{
							v[im] += a*exp(-alpha*r2);
						}
					}
				}

				return v;
			}

			T sft_Rx(T x) const 
			{ 
				return (x-Rx_sf)/Rxy_sc;
			}

			T sft_Ry(T y) const 
			{ 
				return (y-Ry_sf)/Rxy_sc;
			}

			R_2d<T> sft_x_y(T x, T y) const 
			{ 
				x = sft_Rx(x);
				y = sft_Ry(y);
				return R_2d<T>(x, y);
			}

			// not including R2
			void sf_sc()
			{ 
				KS<T> Ixy_sum_t = 0;
				for(auto ixy = 0; ixy < Ixy.size(); ixy++)
				{
					Rx[ixy] = (Rx[ixy] - Rx_sf)/Rxy_sc;
					Ry[ixy] = (Ry[ixy] - Ry_sf)/Rxy_sc;
					T Iv = Ixy[ixy]/Ixy_sc;
					Ixy[ixy] = Iv;
					Ixy_sum_t += Iv;
				}
				Ixy_sum = Ixy_sum_t;
			}
		};
		*/

}

#include "../src/region_3d.inl"