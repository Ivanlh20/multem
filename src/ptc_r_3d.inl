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

#include "ptc_r_3d.h"

/* template specialization 3d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class T>
	Ptc_R_xd<T, edim_3>::Ptc_R_xd(): bs(), x_lim(), y_lim(), z_lim(), r_mean(), r_std(), sz() {}

	template <class T>
	template <class U>
	Ptc_R_xd<T, edim_3>::Ptc_R_xd(const pVctr_cpu_64<U>& ptc, const dt_int64& icol, const R_3d<U>& bs, dt_bool pbc_xy, dt_bool b_statistic)
	{
		set_ptc(ptc, icol, bs, pbc_xy, b_statistic);
	}

	/* copy constructor */
	template <class T>
	Ptc_R_xd<T, edim_3>::Ptc_R_xd(const Ptc_R_xd<T, edim_3>& ptc)
	{
		*this = ptc;
	}

	/* converting constructor */
	template <class T>
	template <class U>
	Ptc_R_xd<T, edim_3>::Ptc_R_xd(const Ptc_R_xd<U, edim_3>& ptc)
	{
		*this = ptc;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T>
	Ptc_R_xd<T, edim_3>& Ptc_R_xd<T, edim_3>::operator=(const Ptc_R_xd<T, edim_3>& ptc)
	{
		assign(ptc);

		return *this;
	}

	/* converting assignment operator */
	template <class T>
	template <class U> 
	Ptc_R_xd<T, edim_3>& Ptc_R_xd<T, edim_3>::operator=(const Ptc_R_xd<U, edim_3>& ptc)
	{
		assign(ptc);

		return *this;
	}

	template <class T>
	template <class U>
	void Ptc_R_xd<T, edim_3>::assign(const Ptc_R_xd<U, edim_3>& ptc)
	{
		if ((void*)this != (void*)&ptc)
		{
			bs = ptc.bs;

			x = ptc.x;
			y = ptc.y;
			z = ptc.z;

			x_lim = ptc.x_lim;
			y_lim = ptc.y_lim;
			z_lim = ptc.z_lim;

			r_mean = ptc.r_mean;
			r_std = ptc.r_std;
			sz = ptc.sz;
		}
	}

	/***************************************************************************************/
	template <class T>
	dt_shape_st<typename Ptc_R_xd<T, edim_3>::size_type> Ptc_R_xd<T, edim_3>::shape() const
	{
		return {size(), cols(), 1, 1};
	}

	template <class T>
	typename Ptc_R_xd<T, edim_3>::size_type Ptc_R_xd<T, edim_3>::size() const
	{
		return x.size();
	}
	
	template <class T>
	dt_int32 Ptc_R_xd<T, edim_3>::size_32() const
	{
		return static_cast<dt_int32>(x.size());
	}

	template <class T>
	typename Ptc_R_xd<T, edim_3>::size_type Ptc_R_xd<T, edim_3>::cols() const
	{
		return 3;
	}

	template <class T>
	dt_bool Ptc_R_xd<T, edim_3>::empty() const
	{
		return size() == 0;
	}

	template <class T>
	void Ptc_R_xd<T, edim_3>::clear()
	{
		bs = 0;

		x.clear();
		y.clear();
		z.clear();

		x_lim = 0;
		y_lim = 0;
		z_lim = 0;

		r_mean = 0;
		r_std = 0;
		sz = 0;
	}

	template <class T>
	void Ptc_R_xd<T, edim_3>::resize(size_type new_size)
	{
		new_size = max(size_type(0), new_size);

		x.resize(new_size);
		y.resize(new_size);
		z.resize(new_size);
	}

	template <class T>
	void Ptc_R_xd<T, edim_3>::reserve(size_type new_size)
	{
		new_size = max(size_type(0), new_size);

		x.reserve(new_size);
		y.reserve(new_size);
		z.reserve(new_size);
	}

	template <class T>
	void Ptc_R_xd<T, edim_3>::shrink_to_fit()
	{
		x.shrink_to_fit();
		y.shrink_to_fit();
		z.shrink_to_fit();
	}

	template <class T>
	void Ptc_R_xd<T, edim_3>::push_back(const R_3d<T>& r)
	{
		x.push_back(r.x);
		y.push_back(r.y);
		z.push_back(r.z);
	}

	template <class T>
	template <class U>
	void Ptc_R_xd<T, edim_3>::set_bs(const R_3d<U>& bs)
	{
		this->bs = bs;
	}

	template <class T>
	R_3d<T> Ptc_R_xd<T, edim_3>::get(const size_type& ia) const
	{
		return {x[ia], y[ia], z[ia]};
	}

	template <class T>
	void Ptc_R_xd<T, edim_3>::set(const size_type& ia, const R_3d<T>& r)
	{
		x[ia] = r.x;
		y[ia] = r.y;
		z[ia] = r.z;
	}

	template <class T>
	R_3d<T> Ptc_R_xd<T, edim_3>::get_pos(const size_type& ia) const
	{
		return {x[ia], y[ia], z[ia]};
	}

	template <class T>
	void Ptc_R_xd<T, edim_3>::set_pos(const size_type& ia, const R_3d<T>& r)
	{
		x[ia] = r.x;
		y[ia] = r.y;
		z[ia] = r.z;
	}

	template <class T>
	template <class U>
	void Ptc_R_xd<T, edim_3>::set_ptc(const Ptc_R_xd<U, edim_3>& ptc, dt_bool pbc_xy, dt_bool b_statistic)
	{
		clear();
		reserve(ptc.size());

		mt::ptc_detail::set_ptc_pbc_xy(ptc, pbc_xy, b_statistic, *this);
	}

	template <class T>
	template <class U>
	void Ptc_R_xd<T, edim_3>::set_ptc(const pVctr_cpu_64<U>& ptc, const dt_int64& icol, const R_3d<U>& bs, dt_bool pbc_xy, dt_bool b_statistic)
	{
		clear();
		reserve(ptc.size());

		mt::ptc_detail::set_ptc_pbc_xy(ptc, icol, bs, pbc_xy, b_statistic, *this);
	}
		
	/* copy data to pointer */
	template <class T>
	template <class U>
	dt_int32 Ptc_R_xd<T, edim_3>::cpy_to_ptr(U *ptc, size_type n_ptc, dt_int32 is_0, dt_int32 is_e) const
	{
		if (is_0>is_e)
		{
			std::swap(is_0, is_e);
		}

		auto n_data = min(n_ptc, size());
		dt_int32 is = 0;

		if (fcn_chk_bound(0, is_0, is_e))
			memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), x.data(), n_data);				// x-position

		if (fcn_chk_bound(1, is_0, is_e))
			memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), y.data(), n_data);				// y-position

		if (fcn_chk_bound(2, is_0, is_e))
			memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), z.data(), n_data);				// z-position

		return is;
	}

	// sort by x
	template <class T>
	void Ptc_R_xd<T, edim_3>::sort_by_x()
	{
		sort_by_idx(0);
	}

	// sort by y
	template <class T>
	void Ptc_R_xd<T, edim_3>::sort_by_y()
	{
		sort_by_idx(1);
	}

	// sort by z
	template <class T>
	void Ptc_R_xd<T, edim_3>::sort_by_z()
	{
		sort_by_idx(2);
	}

	// sort by idx
	template <class T>
	void Ptc_R_xd<T, edim_3>::sort_by_idx(const dt_int32& idx)
	{
		auto first = fcn_mkzipiter_begin(this->x, this->y, this->z);
		auto last = fcn_mkzipiter_end(this->x, this->y, this->z);

		switch(idx)
		{
			case 0: 
			{ 
				thrust::sort(first, last, mt::cgpu_fctr::less_soa<0>()); 
			} 
			break;
			case 1:
			{ 
				thrust::sort(first, last, mt::cgpu_fctr::less_soa<1>()); 
			} 
			break;
			case 2:
			{ 
				thrust::sort(first, last, mt::cgpu_fctr::less_soa<2>()); 
			} 
			break;
		}
	}

	/***************************************************************************************/
	template <class T>
	T Ptc_R_xd<T, edim_3>::norm_2_pbc_xy(const size_type& ia, const R_3d<T>& r_0) const
	{
		auto r = get_pos(ia) - r_0;

		r.x = fabs(r.x);
		r.y = fabs(r.y);

		r.x = ::fmin(r.x, fabs(r.x-bs.x));
		r.y = ::fmin(r.y, fabs(r.y-bs.y));

		return mt::norm_2(r);
	}				
				
	template <class T>
	T Ptc_R_xd<T, edim_3>::norm_2_pbc(const size_type& ia, const R_3d<T>& r_0) const
	{
		auto r = get_pos(ia) - r_0;

		r.x = fabs(r.x);
		r.y = fabs(r.y);
		r.z = fabs(r.z);

		r.x = ::fmin(r.x, fabs(r.x-bs.x));
		r.y = ::fmin(r.y, fabs(r.y-bs.y));
		r.z = ::fmin(r.z, fabs(r.y-bs.z));

		return mt::norm_2(r);
	}

	template <class T>
	T Ptc_R_xd<T, edim_3>::norm_2(const size_type& ia, const R_3d<T>& r_0) const
	{
		return mt::norm_2(get_pos(ia) - r_0);
	}

	template <class T>
	T Ptc_R_xd<T, edim_3>::norm_2(const size_type& ia, const T& x, const T& y, const T& z) const
	{
		return mt::norm_2(get_pos(ia) - R_3d<T>(x, y, z));
	}

	template <class T>
	T Ptc_R_xd<T, edim_3>::norm_2(const size_type& ia_0, const size_type& ia_e) const
	{
		return mt::norm_2(get_pos(ia_0) - get_pos(ia_e));
	}

	/***************************************************************************************/
	template <class T>
	T Ptc_R_xd<T, edim_3>::norm_pbc_xy(const size_type& ia, const R_3d<T>& r_0) const
	{
		return ::sqrt(this->norm_2_pbc_xy(ia, r_0));
	}				
				
	template <class T>
	T Ptc_R_xd<T, edim_3>::norm_pbc(const size_type& ia, const R_3d<T>& r_0) const
	{
		return ::sqrt(this->norm_2_pbc(ia, r_0));
	}

	template <class T>
	T Ptc_R_xd<T, edim_3>::norm(const size_type& ia, const R_3d<T>& r_0) const
	{
		return ::sqrt(this->norm_2(ia, r_0));
	}

	template <class T>
	T Ptc_R_xd<T, edim_3>::norm(const size_type& ia, const T& x, const T& y, const T& z) const
	{
		return ::sqrt(this->norm_2(ia, x, y, z));
	}

	template <class T>
	T Ptc_R_xd<T, edim_3>::norm(const size_type& ia_0, const size_type& ia_e) const
	{
		return ::sqrt(this->norm_2(ia_0, ia_e));
	}

	/***************************************************************************************/
	template <class T>
	void Ptc_R_xd<T, edim_3>::get_statistic()
	{
		fcn_ptc_pos_statistic(*this);
	}

	template <class T>
	void Ptc_R_xd<T, edim_3>::shift(const R_3d<T>& r_sft)
	{
		mt::fcn_ptc_pos_shift(r_sft, *this);
	}

	template <class T>
	void Ptc_R_xd<T, edim_3>::recenter(const R_3d<T>& bs)
	{
		mt::fcn_ptc_pos_recenter(bs, *this);
	}

	template <class T>
	void Ptc_R_xd<T, edim_3>::recenter()
	{
		mt::fcn_ptc_pos_recenter(bs, *this);
	}

	template <class T>
	void Ptc_R_xd<T, edim_3>::recenter_xy(const R_2d<T>& bs)
	{
		mt::fcn_ptc_pos_recenter_xy(bs, *this);
	}

	template <class T>
	void Ptc_R_xd<T, edim_3>::recenter_xy()
	{
		mt::fcn_ptc_pos_recenter_xy({bs.x, bs.y}, *this);
	}

	template <class T>
	void Ptc_R_xd<T, edim_3>::apply_ltf(const Mx_3x3<T>& mx, const R_3d<T>& p)
	{
		mt::fcn_ptc_pos_apply_ltf(mx, p, *this); 
	}

	template <class T>
	void Ptc_R_xd<T, edim_3>::rotate(const T& theta, const R_3d<T>& u_0, const R_3d<T>& p)
	{
		mt::fcn_ptc_pos_rotate(theta, u_0, p, *this);
	}
}

/* fcns 3d */
namespace mt
{
	template <class T>
	void fcn_ptc_pos_statistic(Ptc_R_3d<T>& ptc)
	{
		if (ptc.empty())
		{
			return;
		}

		mt::fcn_minmax_element(ptc.x, ptc.x_lim.x, ptc.x_lim.y);
		mt::(ptc.y, ptc.y_lim.x, ptc.y_lim.y);
		mt::(ptc.z, ptc.z_lim.x, ptc.z_lim.y);
				
		mt::fcn_mean_std(ptc.x, ptc.r_mean.x, ptc.r_std.x);
		mt::fcn_mean_std(ptc.y, ptc.r_mean.y, ptc.r_std.y);
		mt::fcn_mean_std(ptc.z, ptc.r_mean.z, ptc.r_std.z);

		ptc.sz = R_3d<T>(ptc.x_lim.y - ptc.x_lim.x, ptc.y_lim.y - ptc.y_lim.x, ptc.z_lim.y - ptc.z_lim.x);
	}
		
	template <class T>
	void fcn_ptc_pos_shift(const R_3d<T>& r_sft, Ptc_R_3d<T>& ptc)
	{
		for(auto ia = 0; ia < ptc.size(); ia++)
		{
			ptc.x[ia] += r_sft.x;
			ptc.y[ia] += r_sft.y;
			ptc.z[ia] += r_sft.z;
		}

		fcn_ptc_pos_statistic(ptc);
	}

	template <class T>
	void fcn_ptc_pos_recenter(const R_3d<T>& bs, Ptc_R_3d<T>& ptc)
	{
		const R_3d<T> r_sft = (bs-ptc.sz)/T(2) - R_3d<T>(ptc.x_lim.x, ptc.y_lim.x, ptc.z_lim.x);

		fcn_ptc_pos_shift(r_sft, ptc);
	}

	template <class T>
	void fcn_ptc_pos_recenter_xy(const R_2d<T>& bs, Ptc_R_3d<T>& ptc)
	{
		const R_2d<T> r_sft = (bs-R_2d<T>(ptc.sz.x, ptc.sz.y))/T(2) - R_2d<T>(ptc.x_lim.x, ptc.y_lim.x);

		fcn_ptc_pos_shift({r_sft.x, r_sft.y, T(0)}, ptc);
	}

	template <class T>
	void fcn_ptc_pos_apply_ltf(const Mx_3x3<T>& mx, const R_3d<T>& p, Ptc_R_3d<T>& ptc)
	{
		for(dt_int64 ia = 0; ia < ptc.size(); ia++)
		{
			auto r = mx*ptc.get_pos(ia) + p;

			ptc.x[ia] = r.x;
			ptc.y[ia] = r.y;
		}

		fcn_ptc_pos_statistic(ptc);
	}
				
	template <class T>
	void fcn_ptc_pos_rotate(const T& theta, const R_3d<T>& u_0, const R_3d<T>& p, Ptc_R_3d<T>& ptc)
	{
		const auto Rm = fcn_rot_mx_3d(theta, u_0);
		const auto p_sft = p - Rm*p;

		fcn_ptc_pos_apply_ltf(Rm, p_sft, ptc);
	}
}

