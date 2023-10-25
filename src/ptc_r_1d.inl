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

#include "ptc_r_1d.h"

namespace mt
{
	/************************************* constructors ************************************/
	template <class T>
	Ptc_R_xd<T, edim_1>::Ptc_R_xd(): bs(), x_lim(), r_mean(), r_std(), sz() {}

	template <class T>
	template <class U>
	Ptc_R_xd<T, edim_1>::Ptc_R_xd(const pVctr_cpu_64<U>& ptc, const dt_int64& icol, const R_1d<U>& bs, dt_bool pbc_x, dt_bool b_statistic)
	{
		set_ptc(ptc, icol, bs, pbc_x, b_statistic);
	}

	/* copy constructor */
	template <class T>
	Ptc_R_xd<T, edim_1>::Ptc_R_xd(const Ptc_R_xd<T, edim_1>& ptc)
	{
		*this = ptc;
	}

	/* converting constructor */
	template <class T>
	template <class U>
	Ptc_R_xd<T, edim_1>::Ptc_R_xd(const Ptc_R_xd<U, edim_1>& ptc)
	{
		*this = ptc;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T>
	Ptc_R_xd<T, edim_1>& Ptc_R_xd<T, edim_1>::operator=(const Ptc_R_xd<T, edim_1>& ptc)
	{
		assign(ptc);

		return *this;
	}

	/* converting assignment operator */
	template <class T>
	template <class U> 
	Ptc_R_xd<T, edim_1>& Ptc_R_xd<T, edim_1>::operator=(const Ptc_R_xd<U, edim_1>& ptc)
	{
		assign(ptc);

		return *this;
	}

	template <class T>
	template <class U>
	void Ptc_R_xd<T, edim_1>::assign(const Ptc_R_xd<U, edim_1>& ptc)
	{
		if ((void*)this != (void*)&ptc)
		{
			bs = ptc.bs;

			x = ptc.x;

			x_lim = ptc.x_lim;

			r_mean = ptc.r_mean;
			r_std = ptc.r_std;
			sz = ptc.sz;
		}
	}

	/***************************************************************************************/
	template <class T>
	dt_shape_st<typename Ptc_R_xd<T, edim_1>::size_type> Ptc_R_xd<T, edim_1>::shape() const
	{
		return {size(), cols(), 1, 1};
	}

	template <class T>
	typename Ptc_R_xd<T, edim_1>::size_type Ptc_R_xd<T, edim_1>::size() const
	{
		return x.size();
	}

	template <class T>
	dt_int32 Ptc_R_xd<T, edim_1>::size_32() const
	{
		return static_cast<dt_int32>(x.size());
	}	

	template <class T>
	typename Ptc_R_xd<T, edim_1>::size_type Ptc_R_xd<T, edim_1>::cols() const
	{
		return 1;
	}

	template <class T>
	dt_bool Ptc_R_xd<T, edim_1>::empty() const
	{
		return size() == 0;
	}

	template <class T>
	void Ptc_R_xd<T, edim_1>::clear()
	{
		bs = 0;

		x.clear();

		x_lim = 0;

		r_mean = 0;
		r_std = 0;
		sz = 0;
	}

	template <class T>
	void Ptc_R_xd<T, edim_1>::resize(size_type new_size)
	{
		new_size = max(size_type(0), new_size);

		x.resize(new_size);
	}

	template <class T>
	void Ptc_R_xd<T, edim_1>::reserve(size_type new_size)
	{
		new_size = max(size_type(0), new_size);

		x.reserve(new_size);
	}

	template <class T>
	void Ptc_R_xd<T, edim_1>::shrink_to_fit()
	{
		x.shrink_to_fit();
	}

	template <class T>
	void Ptc_R_xd<T, edim_1>::push_back(const R_1d<T>& r)
	{
		x.push_back(r);
	}

	template <class T>
	template <class U>
	void Ptc_R_xd<T, edim_1>::set_bs(const R_2d<U>& bs)
	{
		this->bs = bs;
	}

	template <class T>
	R_1d<T> Ptc_R_xd<T, edim_1>::get(const size_type& ia) const
	{
		return x[ia];
	}

	template <class T>
	void Ptc_R_xd<T, edim_1>::set(const size_type& ia, const R_1d<T>& r)
	{
		x[ia] = r;
	}

	template <class T>
	R_1d<T> Ptc_R_xd<T, edim_1>::get_pos(const size_type& ia) const
	{
		return x[ia];
	}

	template <class T>
	void Ptc_R_xd<T, edim_1>::set_pos(const size_type& ia, const R_1d<T>& r)
	{
		x[ia] = r;
	}

	template <class T>
	template <class U>
	void Ptc_R_xd<T, edim_1>::set_ptc(const Ptc_R_xd<U, edim_1>& ptc, dt_bool pbc_x, dt_bool b_statistic)
	{
		mt::set_ptc_pbc_xy(ptc, pbc_x, b_statistic, *this);
	}

	template <class T>
	template <class U>
	void Ptc_R_xd<T, edim_1>::set_ptc(const pVctr_cpu_64<U>& ptc, const dt_int64& icol, const R_1d<U>& bs, dt_bool pbc_x, dt_bool b_statistic)
	{
		mt::set_ptc_pbc_xy(ptc, icol, bs, pbc_x, b_statistic, *this);
	}
		
	/* copy data to pointer */
	template <class T>
	template <class U>
	dt_int32 Ptc_R_xd<T, edim_1>::cpy_to_ptr(U *ptc, size_type n_ptc, dt_int32 is_0, dt_int32 is_e) const
	{
		if (is_0>is_e)
		{
			std::swap(is_0, is_e);
		}

		auto n_data = min(n_ptc, size());
		dt_int32 is = 0;

		if (fcn_chk_bound(0, is_0, is_e))
			memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), x.data(), n_data);				// x-position

		return is;
	}

	// sort by x
	template <class T>
	void Ptc_R_xd<T, edim_1>::sort_by_x()
	{
		sort_by_idx(0);
	}

	// sort by idx
	template <class T>
	void Ptc_R_xd<T, edim_1>::sort_by_idx(const dt_int32& idx)
	{
		auto first = std::begin(this->x);
		auto last = std::end(this->x);

		switch(idx)
		{
			case 0: 
			{
				thrust::sort(x, last, mt::cgpu_fctr::less<T>()); 
			}
			break;
		}
	}

	/***************************************************************************************/
	template <class T>
	T Ptc_R_xd<T, edim_1>::norm_2_pbc_x(const size_type& ia, const R_1d<T>& r_0) const
	{
		auto r = get_pos(ia) - r_0;

		r = fabs(r);

		r = ::fmin(r, fabs(r-bs));

		return mt::norm_2(r);
	}				
				
	template <class T>
	T Ptc_R_xd<T, edim_1>::norm_2_pbc(const size_type& ia, const R_1d<T>& r_0) const
	{
		auto r = get_pos(ia) - r_0;

		r = fabs(r);

		r = ::fmin(r, fabs(r-bs));

		return mt::norm_2(r);
	}

	template <class T>
	T Ptc_R_xd<T, edim_1>::norm_2(const size_type& ia, const R_1d<T>& r_0) const
	{
		return mt::norm_2(get_pos(ia) - r_0);
	}

	template <class T>
	T Ptc_R_xd<T, edim_1>::norm_2(const size_type& ia_0, const size_type& ia_e) const
	{
		return mt::norm_2(get_pos(ia_0) - get_pos(ia_e));
	}

	/***************************************************************************************/
	template <class T>
	T Ptc_R_xd<T, edim_1>::norm_pbc_x(const size_type& ia, const R_1d<T>& r_0) const
	{
		return ::sqrt(this->norm_2_pbc_x(ia, r_0));
	}				
				
	template <class T>
	T Ptc_R_xd<T, edim_1>::norm_pbc(const size_type& ia, const R_1d<T>& r_0) const
	{
		return ::sqrt(this->norm_2_pbc(ia, r_0));
	}

	template <class T>
	T Ptc_R_xd<T, edim_1>::norm(const size_type& ia, const R_2d<T>& r_0) const
	{
		return ::sqrt(this->norm_2(ia, r_0));
	}

	template <class T>
	T Ptc_R_xd<T, edim_1>::norm(const size_type& ia_0, const size_type& ia_e) const
	{
		return ::sqrt(this->norm_2(ia_0, ia_e));
	}

	/***************************************************************************************/
	template <class T>
	void Ptc_R_xd<T, edim_1>::get_statistic()
	{
		mt::fcn_ptc_pos_statistic(*this);
	}

	template <class T>
	void Ptc_R_xd<T, edim_1>::shift(const R_1d<T>& r_sft)
	{
		mt::fcn_ptc_pos_shift(r_sft, *this);
	}

	template <class T>
	void Ptc_R_xd<T, edim_1>::recenter(const R_1d<T>& bs)
	{
		mt::fcn_ptc_pos_recenter(bs, *this);
	}

	template <class T>
	void Ptc_R_xd<T, edim_1>::recenter()
	{
		mt::fcn_ptc_pos_recenter(bs, *this);
	}

	template <class T>
	void Ptc_R_xd<T, edim_1>::apply_ltf(const T mx, const R_1d<T>& p)
	{
		mt::fcn_ptc_pos_apply_ltf(mx, p, *this); 
	}
}

/* fcns */
namespace mt
{
	template <class T, class U>
	void set_ptc_pbc_x(const Ptc_R_1d<T>& ptc_i, const dt_bool& pbc_x, const dt_bool& b_statistic, Ptc_R_1d<U>& ptc_o)
	{
		ptc_o.bs = ptc_i.bs;

		if (ptc_i.size()==0) 
			return;

		if (!pbc_xy)
		{
			ptc_o = ptc_i;
		}
		else
		{
			const U bs_e = ptc_o.bs - c_dflt_pos_ee;

			for(dt_int64 ia = 0; ia < ptc_i.size(); ia++)
			{
				auto bb_x = fcn_chk_bound(ptc_i.x[ia], U(0), bs_e);

				if (bb_x)
				{
					ptc_o.push_back(ptc_i.get(ia));
				}
			}

			ptc_o.shrink_to_fit();
		}

		if (b_statistic)
		{
			ptc_o.get_statistic();
		}
	}

	template <class T, class U>
	void set_ptc_pbc_x(const pVctr_cpu_64<T>& pvctr, const dt_int64& icol, const R_1d<U>& bs, const dt_bool& pbc_x, const dt_bool& b_statistic, Ptc_R_1d<U>& ptc_o)
	{
		ptc_o.bs = bs;

		if (p_ptc.empty()) 
			return;

		const auto s_0 = p_ptc.s0();
		const auto s_1 = p_ptc.s1();

		if (!pbc_xy)
		{
			for(dt_int64 ia = 0; ia < s_0; ia++)
			{
				ptc_o.push_back(Ptc_s(p_ptc.data(), s_0, s_1, ia, icol));
			}
		}
		else
		{
			auto bs_e = ptc_o.bs - c_dflt_pos_ee;

			for(dt_int64 ia = 0; ia < s_0; ia++)
			{
				auto r = pvctr(ia, icol);

				if (bb_x)
				{
					ptc_o.push_back(r);
				}
			}
		}

		ptc_o.shrink_to_fit();

		if (b_statistic)
		{
			ptc_o.get_statistic();
		}
	}
		
		template <class T>
	void fcn_ptc_pos_statistic(Ptc_R_1d<T>& ptc)
	{
		if (ptc.empty())
		{
			return;
		}

		mt::fcn_minmax_element(ptc.x, ptc.x_lim.x, ptc.x_lim.y);
				
		mt::fcn_mean_std(ptc.x, ptc.r_mean, ptc.r_std);

		ptc.sz = ptc.x_lim.y - ptc.x_lim.x;
	}

	template <class T>
	void fcn_ptc_pos_shift(const R_1d<T>& r_sft, Ptc_R_1d<T>& ptc)
	{
		for(auto ia = 0; ia < ptc.size(); ia++)
		{
			ptc.x[ia] += r_sft.x;
		}

			fcn_ptc_pos_statistic(ptc);
	}

	template <class T>
	void fcn_ptc_pos_recenter(const R_1d<T>& bs, Ptc_R_1d<T>& ptc)
	{
		const R_1d<T> r_sft = (bs-ptc.sz)/T(2) - ptc.x_lim.x;

		fcn_ptc_pos_shift(r_sft, ptc);
	}

	template <class T>
	void fcn_ptc_pos_apply_ltf(const T& mx, const R_1d<T>& p, Ptc_R_1d<T>& ptc)
	{
		for(auto ia = 0; ia < ptc.size(); ia++)
		{
			auto r = mx*ptc.get_pos(ia) + p;

			ptc.x[ia] = r;
		}

		fcn_ptc_pos_statistic(ptc);
	}
}
