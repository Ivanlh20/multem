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

#include "region_3d.h"

/* template specialization 3d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class T>
	CGPU_EXEC
	Region_Rect_xd<T, edim_3>::Region_Rect_xd(): Region_Rect_xd<T, edim_2>(), rz_0(0), rz_e(0) {}

	template <class T>
	Region_Rect_xd<T, edim_3>::Region_Rect_xd(const T& rx_0, const T& rx_e, const T& ry_0, const T& ry_e, const T& rz_0, const T& rz_e)
	{
		set_in_data(rx_0, rx_e, ry_0, ry_e, rz_0, rz_e);
	}

	template <class T>
	template <class U>
	Region_Rect_xd<T, edim_3>::Region_Rect_xd(U *data, const size_type& n_data)
	{
		Region_Rect_xd<T, edim_2>(data, n_data);

		rz_0 = (n_data>4)?data[4]:T(0);
		rz_e = (n_data>5)?data[5]:T(0);
	}

	/* copy constructor */
	template <class T>
	CGPU_EXEC
	Region_Rect_xd<T, edim_3>::Region_Rect_xd(const Region_Rect_xd<T, edim_3>& region)
	{
		*this = region;
	}

	/* converting constructor */
	template <class T>
	template <class U>
	CGPU_EXEC
	Region_Rect_xd<T, edim_3>::Region_Rect_xd(const Region_Rect_xd<U, edim_3>& region)
	{
		*this = region;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T>
	CGPU_EXEC
	Region_Rect_xd<T, edim_3>& Region_Rect_xd<T, edim_3>::operator=(const Region_Rect_xd<T, edim_3>& region)
	{
		if (this != &region)
		{
			Region_Rect_xd<T, edim_2>::operator=(region);

			rz_0 = region.rz_0;
			rz_e = region.rz_e;
		}

		return *this;
	}			
	
	/* converting assignment operator */
	template <class T>
	template <class U>
	CGPU_EXEC
	Region_Rect_xd<T, edim_3>& Region_Rect_xd<T, edim_3>::operator=(const Region_Rect_xd<U, edim_3>& region)
	{
		if (this != &region)
		{
			Region_Rect_xd<T, edim_2>::operator=(region);

			rz_0 = T(region.rz_0);
			rz_e = T(region.rz_e);
		}

		return *this;
	}

	template <class T>
	template <class U> 
	CGPU_EXEC
	void Region_Rect_xd<T, edim_3>::assign(const Region_Rect_xd<U, edim_3>& region)
	{
		*this = region;
	}

	/***************************************************************************************/
	template <class T>
	template <class U>
	void Region_Rect_xd<T, edim_3>::set_in_data(const U& rx_0, const U& rx_e, const U& ry_0, const U& ry_e, const U& rz_0, const U& rz_e)
	{
		Region_Rect_xd<T, edim_2>::set_in_data(rx_0, rx_e, ry_0, ry_e);

		this->rz_0 = rz_0;
		this->rz_e = rz_e;
	}

	template <class T>
	CGPU_EXEC
	void Region_Rect_xd<T, edim_3>::clear()
	{
		Region_Rect_xd<T, edim_2>::clear();

		rz_0 = T(0);
		rz_e = T(0);
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_3>::bs_z() const 
	{ 
		return fcn_max(rz_e-rz_0, T(0));
	}

	template <class T>
	CGPU_EXEC
	R_3d<T> Region_Rect_xd<T, edim_3>::bs() const 
	{ 
		return {this->bs_x(), this->bs_y(), bs_z()};
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_3>::bs_z_h() const 
	{ 
		return bs_z()/T(2);
	}

	template <class T>
	CGPU_EXEC
	R_3d<T> Region_Rect_xd<T, edim_3>::bs_h() const 
	{ 
		return {this->bs_x_h(), this->bs_y_h(), bs_z_h()};
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_3>::rz_c() const 
	{ 
		return (rz_e+rz_0)/T(2);
	}
	
	template <class T>
	CGPU_EXEC
	R_3d<T> Region_Rect_xd<T, edim_3>::r_c() const 
	{ 
		return {this->rx_c(), this->ry_c(), rz_c()};
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_3>::radius_z() const
	{
		return bs_z_h();
	}

	template <class T>
	CGPU_EXEC
	R_3d<T> Region_Rect_xd<T, edim_3>::radius() const
	{
		return {this->radius_x(), this->radius_y(), radius_z()};
	}

	template <class T>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_3>::radius_z_p(const T& p) const
	{
		return (T(1)-p)*radius_z();
	}

	template <class T>
	CGPU_EXEC
	R_3d<T> Region_Rect_xd<T, edim_3>::radius_p(const T& p) const
	{
		return {this->radius_x_p(p), this->radius_y_p(p), radius_z_p(p)};
	}

	template <class T>
	void Region_Rect_xd<T, edim_3>::sft_region(const R_3d<T>& bs, const R_3d<T>& dr)
	{
		Region_Rect_xd<T, edim_2>::sft_region({bs.x, bs.y}, {dr.x, dr.y});

		rz_0 = fcn_set_bound(rz_0 + dr.z, T(0), bs.z);
		rz_e = fcn_set_bound(rz_e + dr.z, T(0), bs.z);
	}

	/***************************************************************************************/
	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_3>::chk_bound_z(const T& z) const
	{
		return fcn_chk_bound(z, rz_0, rz_e);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_3>::chk_bound(const R_3d<T>& r) const
	{
		return this->chk_bound_x(r.x) && this->chk_bound_y(r.y) && chk_bound_z(r.z);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_3>::chk_bound(const T& rx, const T& ry, const T& rz)const
	{
		return this->chk_bound_x(rx) && this->chk_bound_y(ry) && chk_bound_z(rz);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_3>::chk_bound_z_eps(const T& z) const 
	{ 
		return fcn_chk_bound_eps(z, rz_0, rz_e);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_3>::chk_bound_eps(const R_3d<T>& r) const 
	{ 
		return this->chk_bound_x_eps(r.x) && this->chk_bound_y_eps(r.y) && chk_bound_z_eps(r.z);
	}

	template <class T>
	CGPU_EXEC
	dt_bool Region_Rect_xd<T, edim_3>::chk_bound_eps(const T& rx, const T& ry, const T& rz)const
	{
		return this->chk_bound_x_eps(rx) && this->chk_bound_y_eps(ry) && chk_bound_z_eps(rz);
	}

	template <class T>
	template <class U, class>
	CGPU_EXEC
	T Region_Rect_xd<T, edim_3>::sub_2_ind(const T& ix, const T& iy, const T& iz) const 
	{ 
		return (iy-this->ry_0) + (ix-this->rx_0)*this->bs_y() + (iz-rz_0)*this->bs_y()*this->bs_x();
	}
}

// /***************************************************************************************/
// /* symmetric region */
// namespace mt
// {
// /*
// template <class TVctr>
// 	struct Region_Rad_1d
// 	{
// 	public:
// 		using T = class TVctr::value_type;
// 		using size_type = dt_uint64;

// 		TVctr Rx;
// 		TVctr R2;
// 		TVctr Ixy;

// 		T R_max;
// 		T Rx_sf;

// 		T Ixy_sf;
// 		T Ixy_sc;
// 		T Ixy_sum;

// 		Region_Rad_1d(): Rx_sf(0), Ixy_sf(0), Ixy_sc(1), Ixy_sum(0)
// 		{
// 		}

// 		size_type size() const
// 		{
// 			return Ixy.size();
// 		}

// 		void clear()
// 		{
// 			Rx.clear();
// 			R2.clear();
// 			Ixy.clear();
// 		}

// 		void reserve(const size_type& new_size)
// 		{
// 			Rx.reserve(new_size);
// 			R2.reserve(new_size);
// 			Ixy.reserve(new_size);
// 		}

// 		void shrink_to_fit()
// 		{
// 			Rx.shrink_to_fit();
// 			R2.shrink_to_fit();
// 			Ixy.shrink_to_fit();
// 		}

// 		TVctr sft_Ixy(T bg)
// 		{
// 			TVctr Ixy_s;
// 			Ixy_s.reserve(Ixy.size());

// 			for(auto ixy=0; ixy<Ixy.size(); ixy++)
// 			{
// 				Ixy_s.push_back(Ixy[ixy]-bg);
// 			}
// 			return Ixy_s;
// 		}

// 		T sft_Rx(T x) const 
// 		{ 
// 			return (x-Rx_sf)/Rxy_sc;
// 		}
// 	};

// 	template <class TVctr>
// 	struct Region_Rad_2d
// 	{
// 	public:
// 		using T = class TVctr::value_type;
// 		using size_type = dt_uint64;

// 		TVctr Rx;
// 		TVctr Ry;
// 		TVctr R2;
// 		TVctr Ixy;

// 		T R_max;
// 		T Rx_sf;
// 		T Ry_sf;
// 		T Rxy_sc;

// 		T Ixy_sf;
// 		T Ixy_sc;
// 		T Ixy_sum;

// 		Region_Rad_2d(): Rx_sf(0), Ry_sf(0), Rxy_sc(1), 
// 		Ixy_sf(0), Ixy_sc(1), Ixy_sum(0) {}

// 		size_type size() const
// 		{
// 			return Ixy.size();
// 		}

// 		void clear()
// 		{
// 			Rx.clear();
// 			Ry.clear();
// 			R2.clear();
// 			Ixy.clear();
// 		}

// 		void reserve(const size_type& new_size)
// 		{
// 			Rx.reserve(new_size);
// 			Ry.reserve(new_size);
// 			R2.reserve(new_size);
// 			Ixy.reserve(new_size);
// 		}

// 		void shrink_to_fit()
// 		{
// 			Rx.shrink_to_fit();
// 			Ry.shrink_to_fit();
// 			R2.shrink_to_fit();
// 			Ixy.shrink_to_fit();
// 		}

// 		TVctr sft_Ixy(T bg)
// 		{
// 			TVctr Ixy_s;
// 			Ixy_s.reserve(Ixy.size());

// 			for(auto ixy=0; ixy<Ixy.size(); ixy++)
// 			{
// 				Ixy_s.push_back(Ixy[ixy]-bg);
// 			}
// 			return Ixy_s;
// 		}

// 		TVctr sub_region_to_Ixy(Grid_2d<T>& grid_2d, TVctr& Im_s, T x, T y)
// 		{
// 			TVctr v = Ixy;

// 			T R2_max = ::square(R_max);

// 			R_2d<T> p(x, y);
// 			auto range = grid_2d.ithread(p, R_max);
// 			dt_int32 iv = 0;
// 			for(auto ix = range.ix_0; ix < range.ix_e; ix++)
// 			{
// 				for(auto iy = range.iy_0; iy < range.iy_e; iy++)
// 				{
// 					T r2 = grid_2d.r2(ix, iy, p.x, p.y);
// 					if (r2 < R2_max)
// 					{
// 						v[iv++] -= Im_s[grid_2d.sub_2_ind(ix, iy)]/Ixy_sc;
// 					}
// 				}
// 			}

// 			return v;
// 		}

// 		TVctr sft_Ixy(Grid_2d<T>& grid_2d, TVctr& Im_s, T x, T y, T a, T s)
// 		{
// 			TVctr v = sub_region_to_Ixy(grid_2d, Im_s, x*Rxy_sc+Rx_sf, y*Rxy_sc+Ry_sf);

// 			T alpha = T(0.5)/::square(s);
// 			T r2_l = ::square(T(4)*s);
// 			for(auto im = 0; im < v.size(); im++)
// 			{
// 				T rx = Rx[im]-x;
// 				T ry = Ry[im]-y;
// 				T r2 = rx*rx+ry*ry;
// 				if (r2<r2_l)
// 				{
// 					v[im] += a*exp(-alpha*r2);
// 				}
// 			}

// 			return v;
// 		}

// 		TVctr sft_Ixy(Grid_2d<T>& grid_2d, TVctr& Im_s, TVctr& x, TVctr& y, TVctr& A, TVctr& S)
// 		{
// 			TVctr v = sub_region_to_Ixy(grid_2d, Im_s, x[0]*Rxy_sc+Rx_sf, y[0]*Rxy_sc+Ry_sf);

// 			for(auto ip = 0; ip < x.size(); ip++)
// 			{
// 				T a = A[ip];
// 				T b = S[ip];
// 				T alpha = 0.5/pow(b, 2);
// 				T r2_l = pow(4.0*b, 2);
// 				for(auto im = 0; im < v.size(); im++)
// 				{
// 					T rx = Rx[im]-x[ip];
// 					T ry = Ry[im]-y[ip];
// 					T r2 = rx*rx+ry*ry;
// 					if (r2<r2_l)
// 					{
// 						v[im] += a*exp(-alpha*r2);
// 					}
// 				}
// 			}

// 			return v;
// 		}

// 		T sft_Rx(T x) const 
// 		{ 
// 			return (x-Rx_sf)/Rxy_sc;
// 		}

// 		T sft_Ry(T y) const 
// 		{ 
// 			return (y-Ry_sf)/Rxy_sc;
// 		}

// 		R_2d<T> sft_x_y(T x, T y) const 
// 		{ 
// 			x = sft_Rx(x);
// 			y = sft_Ry(y);
// 			return R_2d<T>(x, y);
// 		}

// 		// not including R2
// 		void sf_sc()
// 		{ 
// 			KS<T> Ixy_sum_t = 0;
// 			for(auto ixy = 0; ixy < Ixy.size(); ixy++)
// 			{
// 				Rx[ixy] = (Rx[ixy] - Rx_sf)/Rxy_sc;
// 				Ry[ixy] = (Ry[ixy] - Ry_sf)/Rxy_sc;
// 				T Iv = Ixy[ixy]/Ixy_sc;
// 				Ixy[ixy] = Iv;
// 				Ixy_sum_t += Iv;
// 			}
// 			Ixy_sum = Ixy_sum_t;
// 		}
// 	};
// 	*/

// }