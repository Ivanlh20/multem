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

#include "ptc_r_2d.h"

/* template specialization 2d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class T>
	Ptc_R_xd<T, edim_2>::Ptc_R_xd(): bs(), x_lim(), y_lim(), r_mean(), r_std(), sz() {}

	template <class T>
	template <class U>
	Ptc_R_xd<T, edim_2>::Ptc_R_xd(const pVctr_cpu_64<U>& ptc, const dt_int64& icol, const R_2d<U>& bs, dt_bool pbc_xy, dt_bool b_statistic)
	{
		set_ptc(ptc, icol, bs, pbc_xy, b_statistic);
	}

	/* copy constructor */
	template <class T>
	Ptc_R_xd<T, edim_2>::Ptc_R_xd(const Ptc_R_xd<T, edim_2>& ptc)
	{
		*this = ptc;
	}

	/* converting constructor */
	template <class T>
	template <class U>
	Ptc_R_xd<T, edim_2>::Ptc_R_xd(const Ptc_R_xd<U, edim_2>& ptc)
	{
		*this = ptc;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T>
	Ptc_R_xd<T, edim_2>& Ptc_R_xd<T, edim_2>::operator=(const Ptc_R_xd<T, edim_2>& ptc)
	{
		assign(ptc);

		return *this;
	}

	/* converting assignment operator */
	template <class T>
	template <class U> 
	Ptc_R_xd<T, edim_2>& Ptc_R_xd<T, edim_2>::operator=(const Ptc_R_xd<U, edim_2>& ptc)
	{
		assign(ptc);

		return *this;
	}

	template <class T>
	template <class U>
	void Ptc_R_xd<T, edim_2>::assign(const Ptc_R_xd<U, edim_2>& ptc)
	{
		if ((void*)this != (void*)&ptc)
		{
			bs = ptc.bs;

			x = ptc.x;
			y = ptc.y;

			x_lim = ptc.x_lim;
			y_lim = ptc.y_lim;

			r_mean = ptc.r_mean;
			r_std = ptc.r_std;
			sz = ptc.sz;
		}
	}

	/***************************************************************************************/
	template <class T>
	dt_shape_st<typename Ptc_R_xd<T, edim_2>:: size_type> Ptc_R_xd<T, edim_2>::shape() const
	{
		return {size(), cols(), 1, 1};
	}

	template <class T>
	typename Ptc_R_xd<T, edim_2>::size_type Ptc_R_xd<T, edim_2>::size() const
	{
		return x.size();
	}

	template <class T>
	dt_int32 Ptc_R_xd<T, edim_2>::size_32() const
	{
		return static_cast<dt_int32>(x.size());
	}	

	template <class T>
	typename Ptc_R_xd<T, edim_2>::size_type Ptc_R_xd<T, edim_2>::cols() const
	{
		return 2;
	}

	template <class T>
	dt_bool Ptc_R_xd<T, edim_2>::empty() const
	{
		return size() == 0;
	}

	template <class T>
	void Ptc_R_xd<T, edim_2>::clear()
	{
		bs = 0;

		x.clear();
		y.clear();

		x_lim = 0;
		y_lim = 0;

		r_mean = 0;
		r_std = 0;
		sz = 0;
	}

	template <class T>
	void Ptc_R_xd<T, edim_2>::resize(size_type new_size)
	{
		new_size = max(size_type(0), new_size);

		x.resize(new_size);
		y.resize(new_size);
	}

	template <class T>
	void Ptc_R_xd<T, edim_2>::reserve(size_type new_size)
	{
		new_size = max(size_type(0), new_size);

		x.reserve(new_size);
		y.reserve(new_size);
	}

	template <class T>
	void Ptc_R_xd<T, edim_2>::shrink_to_fit()
	{
		x.shrink_to_fit();
		y.shrink_to_fit();
	}

	template <class T>
	void Ptc_R_xd<T, edim_2>::push_back(const R_2d<T>& r)
	{
		x.push_back(r.x);
		y.push_back(r.y);
	}

	template <class T>
	template <class U>
	void Ptc_R_xd<T, edim_2>::set_bs(const R_2d<U>& bs)
	{
		this->bs = bs;
	}

	template <class T>
	R_2d<T> Ptc_R_xd<T, edim_2>::get(const size_type& ia) const
	{
		return {x[ia], y[ia]};
	}

	template <class T>
	void Ptc_R_xd<T, edim_2>::set(const size_type& ia, const R_2d<T>& r)
	{
		x[ia] = r.x;
		y[ia] = r.y;
	}

	template <class T>
	R_2d<T> Ptc_R_xd<T, edim_2>::get_pos(const size_type& ia) const
	{
		return {x[ia], y[ia]};
	}

	template <class T>
	void Ptc_R_xd<T, edim_2>::set_pos(const size_type& ia, const R_2d<T>& r)
	{
		x[ia] = r.x;
		y[ia] = r.y;
	}

	template <class T>
	template <class U>
	void Ptc_R_xd<T, edim_2>::set_ptc(const Ptc_R_xd<U, edim_2>& ptc, dt_bool pbc_xy, dt_bool b_statistic)
	{
		mt::ptc_detail::set_ptc_pbc_xy(ptc, pbc_xy, b_statistic, *this);
	}

	template <class T>
	template <class U>
	void Ptc_R_xd<T, edim_2>::set_ptc(const pVctr_cpu_64<U>& ptc, const dt_int64& icol, const R_2d<U>& bs, dt_bool pbc_xy, dt_bool b_statistic)
	{
		mt::ptc_detail::set_ptc_pbc_xy(ptc, icol, bs, pbc_xy, b_statistic, *this);
	}
		
	/* copy data to pointer */
	template <class T>
	template <class U>
	dt_int32 Ptc_R_xd<T, edim_2>::cpy_to_ptr(U *ptc, size_type n_ptc, dt_int32 is_0, dt_int32 is_e) const
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

		return is;
	}

	// sort by x
	template <class T>
	void Ptc_R_xd<T, edim_2>::sort_by_x()
	{
		sort_by_idx(0);
	}

	// sort by y
	template <class T>
	void Ptc_R_xd<T, edim_2>::sort_by_y()
	{
		sort_by_idx(1);
	}

	// sort by idx
	template <class T>
	void Ptc_R_xd<T, edim_2>::sort_by_idx(const dt_int32& idx)
	{
		auto first = fcn_mkzipiter_begin(this->x, this->y);
		auto last = fcn_mkzipiter_end(this->x, this->y);

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
		}
	}

	/***************************************************************************************/
	template <class T>
	T Ptc_R_xd<T, edim_2>::norm_2_pbc_xy(const size_type& ia, const R_2d<T>& r_0) const
	{
		auto r = get_pos(ia) - r_0;

		r.x = ::fabs(r.x);
		r.y = ::fabs(r.y);

		r.x = ::fmin(r.x, ::fabs(r.x-bs.x));
		r.y = ::fmin(r.y, ::fabs(r.y-bs.y));

		return mt::norm_2(r);
	}

	template <class T>
	T Ptc_R_xd<T, edim_2>::norm_2_pbc(const size_type& ia, const R_2d<T>& r_0) const
	{
		auto r = get_pos(ia) - r_0;

		r.x = fabs(r.x);
		r.y = fabs(r.y);

		r.x = ::fmin(r.x, fabs(r.x-bs.x));
		r.y = ::fmin(r.y, fabs(r.y-bs.y));

		return mt::norm_2(r);
	}

	template <class T>
	T Ptc_R_xd<T, edim_2>::norm_2(const size_type& ia, const R_2d<T>& r_0) const
	{
		return mt::norm_2(get_pos(ia) - r_0);
	}

	template <class T>
	T Ptc_R_xd<T, edim_2>::norm_2(const size_type& ia, const T& x, const T& y) const
	{
		return mt::norm_2(get_pos(ia) - R_2d<T>(x, y));
	}

	template <class T>
	T Ptc_R_xd<T, edim_2>::norm_2(const size_type& ia_0, const size_type& ia_e) const
	{
		return mt::norm_2(get_pos(ia_0) - get_pos(ia_e));
	}

	/***************************************************************************************/
	template <class T>
	T Ptc_R_xd<T, edim_2>::norm_pbc_xy(const size_type& ia, const R_2d<T>& r_0) const
	{
		return ::sqrt(this->norm_2_pbc_xy(ia, r_0));
	}				
				
	template <class T>
	T Ptc_R_xd<T, edim_2>::norm_pbc(const size_type& ia, const R_2d<T>& r_0) const
	{
		return ::sqrt(this->norm_2_pbc(ia, r_0));
	}

	template <class T>
	T Ptc_R_xd<T, edim_2>::norm(const size_type& ia, const R_2d<T>& r_0) const
	{
		return ::sqrt(this->norm_2(ia, r_0));
	}

	template <class T>
	T Ptc_R_xd<T, edim_2>::norm(const size_type& ia, const T& x, const T& y) const
	{
		return ::sqrt(this->norm_2(ia, x, y));
	}

	template <class T>
	T Ptc_R_xd<T, edim_2>::norm(const size_type& ia_0, const size_type& ia_e) const
	{
		return ::sqrt(this->norm_2(ia_0, ia_e));
	}

	/***************************************************************************************/
	template <class T>
	void Ptc_R_xd<T, edim_2>::get_statistic()
	{
		mt::fcn_ptc_pos_statistic(*this);
	}

	template <class T>
	void Ptc_R_xd<T, edim_2>::shift(const R_2d<T>& r_sft)
	{
		mt::fcn_ptc_pos_shift(r_sft, *this);
	}

	template <class T>
	void Ptc_R_xd<T, edim_2>::recenter(const R_2d<T>& bs)
	{
		mt::fcn_ptc_pos_recenter(bs, *this);
	}

	template <class T>
	void Ptc_R_xd<T, edim_2>::recenter()
	{
		mt::fcn_ptc_pos_recenter(bs, *this);
	}

	template <class T>
	void Ptc_R_xd<T, edim_2>::apply_ltf(const Mx_2x2<T>& mx, const R_2d<T>& p)
	{
		mt::fcn_ptc_pos_apply_ltf(mx, p, *this); 
	}

	template <class T>
	void Ptc_R_xd<T, edim_2>::rotate(const T& theta, const R_2d<T>& p)
	{
		mt::fcn_ptc_pos_rotate(theta, p, *this);
	}
}

/* fcns 2d */
namespace mt
{
	template <class T>
	void fcn_ptc_pos_statistic(Ptc_R_2d<T>& ptc)
	{
		if (ptc.empty())
		{
			return;
		}

		mt::fcn_minmax_element(ptc.x, ptc.x_lim.x, ptc.x_lim.y);
		mt::fcn_minmax_element(ptc.y, ptc.y_lim.x, ptc.y_lim.y);
				
		mt::fcn_mean_std(ptc.x, ptc.r_mean.x, ptc.r_std.x);
		mt::fcn_mean_std(ptc.y, ptc.r_mean.y, ptc.r_std.y);

		ptc.sz = R_2d<T>(ptc.x_lim.y - ptc.x_lim.x, ptc.y_lim.y - ptc.y_lim.x);
	}

	template <class T>
	void fcn_ptc_pos_shift(const R_2d<T>& r_sft, Ptc_R_2d<T>& ptc)
	{
		for(auto ia = 0; ia < ptc.size(); ia++)
		{
			ptc.x[ia] += r_sft.x;
			ptc.y[ia] += r_sft.y;
		}

			fcn_ptc_pos_statistic(ptc);
	}

	template <class T>
	void fcn_ptc_pos_recenter(const R_2d<T>& bs, Ptc_R_2d<T>& ptc)
	{
		const R_2d<T> r_sft = (bs-ptc.sz)/T(2) - R_2d<T>(ptc.x_lim.x, ptc.y_lim.x);

		fcn_ptc_pos_shift(r_sft, ptc);
	}

	template <class T>
	void fcn_ptc_pos_apply_ltf(const Mx_2x2<T>& mx, const R_2d<T>& p, Ptc_R_2d<T>& ptc)
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
	void fcn_ptc_pos_rotate(const T& theta, const R_2d<T>& p, Ptc_R_2d<T>& ptc)
	{
		const auto Rm = fcn_rot_mx_2d(theta);
		const auto p_sft = p - Rm*p;

		fcn_ptc_pos_apply_ltf(Rm, p_sft, ptc);
	}
}