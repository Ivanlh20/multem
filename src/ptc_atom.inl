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

#include "ptc_atom.h"

/* template Ptc_s_Atom */
namespace mt
{
	template <class T>
	Ptc_s_Atom<T>::Ptc_s_Atom():Z(0), Ptc_s_3d_0<T>(), sigma(0), occ(0), tag(0), charge(0) {};

	template <class T>
	Ptc_s_Atom<T>::Ptc_s_Atom(const dt_int32& Z, const T& x, const T& y, const T& z, 
	dt_float32 sigma, dt_float32 occ, dt_int32 tag, dt_int32 charge):
	Z(Z), Ptc_s_3d_0<T>(x, y, z), sigma(sigma), occ(occ), tag(tag), charge(charge) {};
			
	template <class T>
	Ptc_s_Atom<T>::Ptc_s_Atom(const dt_int32& Z, const R_3d<T>& r, 
	dt_float32 sigma, dt_float32 occ, dt_int32 tag, dt_int32 charge):
	Z(Z), Ptc_s_3d_0<T>(r.x, r.y, r.z), sigma(sigma), occ(occ), tag(tag), charge(charge) {};

	/* constructor by pointer */
	template <class T>
	template <class U>
	Ptc_s_Atom<T>::Ptc_s_Atom(U* v, const dt_int64& n_r, const dt_int64& n_c, const dt_int64& idx, dt_int64 icol)
	{
		const auto ip = icol*n_r + idx;

		Z = dt_int32(v[ip + 0*n_r]);								// atomic number
		this->x = (n_c>1)?T(v[ip + 1*n_r]):T(0);					// x-position
		this->y = (n_c>2)?T(v[ip + 2*n_r]):T(0);					// y-position
		this->z = (n_c>3)?T(v[ip + 3*n_r]):T(0);					// z-position

		sigma = (n_c>4)?dt_float32(v[ip + 4*n_r]):c_dflt_rms3d;		// standard deviation
		occ = (n_c>5)?dt_float32(v[ip + 5*n_r]):c_dflt_occ;			// occupancy
		tag = (n_c>6)?dt_int32(v[ip + 6*n_r]):c_dflt_tag;			// tag
		charge = (n_c>7)?dt_int32(v[ip + 7*n_r]):c_dflt_charge;		// charge
	}

	/* copy constructor */
	template <class T>
	Ptc_s_Atom<T>::Ptc_s_Atom(const Ptc_s_Atom<T>& ptc_s)
	{
		*this = ptc_s;
	}

	/* converting constructor */
	template <class T>
	template <class U> 
	Ptc_s_Atom<T>::Ptc_s_Atom(const Ptc_s_Atom<U>& ptc_s)
	{
		*this = ptc_s;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T>
	Ptc_s_Atom<T>& Ptc_s_Atom<T>::operator=(const Ptc_s_Atom<T>& ptc_s)
	{
		if (this != &ptc_s)
		{
			Z = ptc_s.Z;
			this->x = ptc_s.x;
			this->y = ptc_s.y;
			this->z = ptc_s.z;
			sigma = ptc_s.sigma;
			occ = ptc_s.occ;
			tag = ptc_s.tag;
			charge = ptc_s.charge;
		}
			
		return *this;
	}

	/* converting assignment operator */
	template <class T>
	template <class U> 
	Ptc_s_Atom<T>& Ptc_s_Atom<T>::operator=(const Ptc_s_Atom<U>& ptc_s)
	{
		assign(ptc_s);
			
		return *this;
	}

	template <class T>
	template <class U> 
	void Ptc_s_Atom<T>::assign(const Ptc_s_Atom<U>& ptc_s)
	{ 
		if ((void*)this != (void*)&ptc_s)
		{
			Z = ptc_s.Z;
			this->x = T(ptc_s.x);
			this->y = T(ptc_s.y);
			this->z = T(ptc_s.z);
			sigma = ptc_s.sigma;
			occ = ptc_s.occ;
			tag = ptc_s.tag;
			charge = ptc_s.charge;
		}
	}
}

/* template Ptc_Atom */
namespace mt
{
	/************************************* constructors ************************************/
	template <class T>
	Ptc_Atom<T>::Ptc_Atom(): Ptc_R_3d<T>(), cols_used(4), Z_lim(), sigma_lim(), occ_lim(), tag_lim(), charge_lim() {}

	template <class T>
	template <class U>
	Ptc_Atom<T>::Ptc_Atom(const pVctr_cpu_64<U>& ptc, const R_3d<U>& bs, dt_bool pbc_xy, dt_bool b_statistic)
	{
		set_ptc(ptc, bs, pbc_xy, b_statistic);
	}

	/* copy constructor */
	template <class T>
	Ptc_Atom<T>::Ptc_Atom(const Ptc_Atom<T>& ptc)
	{
		*this = ptc;
	}

	/* converting constructor */
	template <class T>
	template <class U>
	Ptc_Atom<T>::Ptc_Atom(const Ptc_Atom<U>& ptc)
	{
		*this = ptc;
	}

	/******************************** assignment operators *********************************/
	/* assignment operator */
	template <class T>
	Ptc_Atom<T>& Ptc_Atom<T>::operator=(const Ptc_Atom<T>& ptc)
	{
		assign(ptc);

		return *this;
	}

	/* converting assignment operator */
	template <class T>
	template <class U> 
	Ptc_Atom<T>& Ptc_Atom<T>::operator=(const Ptc_Atom<U>& ptc)
	{
		assign(ptc);

		return *this;
	}

	template <class T>
	template <class U>
	void Ptc_Atom<T>::assign(const Ptc_Atom<U>& ptc)
	{
		if ((void*)this != (void*)&ptc)
		{
			cols_used = ptc.cols_used;

			Z = ptc.Z;

			Ptc_R_3d<T>::assign(ptc);

			sigma = ptc.sigma;
			occ = ptc.occ;
			tag = ptc.tag;
			charge = ptc.charge;

			Z_lim = ptc.Z_lim;
			sigma_lim = ptc.sigma_lim;
			occ_lim = ptc.occ_lim;
			tag_lim = ptc.tag_lim;
			charge_lim = ptc.charge_lim;
		}
	}

	/***************************************************************************************/
	template <class T>
	Ptc_Atom<T>::size_type Ptc_Atom<T>::cols() const
	{
		return 8;	
	}		

	template <class T>
	void Ptc_Atom<T>::clear()
	{
		cols_used = 4;

		Z.clear();

		Ptc_R_3d<T>::clear();

		sigma.clear();
		occ.clear();
		tag.clear();
		charge.clear();

		Z_lim = 0;
		sigma_lim = 0;
		occ_lim = 0;
		tag_lim = 0;
		charge_lim = 0;
	}

	template <class T>
	void Ptc_Atom<T>::resize(size_type new_size)
	{
		new_size = max(size_type(0), new_size);

		Z.resize(new_size);

		Ptc_R_3d<T>::resize(new_size);

		if (cols_used>4)
			sigma.resize(new_size);

		if (cols_used>5)
			occ.resize(new_size);

		if (cols_used>6)
			tag.resize(new_size);

		if (cols_used>7)
			charge.resize(new_size);
	}

	template <class T>
	void Ptc_Atom<T>::reserve(size_type new_size)
	{
		new_size = max(size_type(0), new_size);

		Z.reserve(new_size);

		Ptc_R_3d<T>::reserve(new_size);

		if (cols_used>4)
			sigma.reserve(new_size);

		if (cols_used>5)
			occ.reserve(new_size);

		if (cols_used>6)
			tag.reserve(new_size);

		if (cols_used>7)
			charge.reserve(new_size);
	}

	template <class T>
	void Ptc_Atom<T>::shrink_to_fit()
	{
		Z.shrink_to_fit();

		Ptc_R_3d<T>::shrink_to_fit();		

		sigma.shrink_to_fit();
		occ.shrink_to_fit();
		tag.shrink_to_fit();
		charge.shrink_to_fit();
	}

	template <class T>
	void Ptc_Atom<T>::push_back(const Ptc_s& ptc_s)
	{
		Z.push_back(ptc_s.Z);							// atomic number

		Ptc_R_3d<T>::push_back(ptc_s);					// xyz

		if (cols_used>4)
			sigma.push_back(ptc_s.sigma);				// standard deviation

		if (cols_used>5)
			occ.push_back(ptc_s.occ);					// occupancy

		if (cols_used>6)
			tag.push_back(abs(ptc_s.tag));				// tag

		if (cols_used>7)
			charge.push_back(ptc_s.charge);				// charge
	}

	template <class T>
	dt_float32 Ptc_Atom<T>::get_sigma(const size_type& ia) const
	{
		return (cols_used>4)?sigma[ia]:c_dflt_rms3d;		// standard deviation
	}

	template <class T>
	dt_float32 Ptc_Atom<T>::get_occ(const size_type& ia) const
	{
		return (cols_used>5)?occ[ia]:c_dflt_occ;			// occupancy
	}

	template <class T>
	dt_int32 Ptc_Atom<T>::get_tag(const size_type& ia) const
	{
		return (cols_used>6)?tag[ia]:c_dflt_tag;			// tag
	}

	template <class T>
	dt_int32 Ptc_Atom<T>::get_charge(const size_type& ia) const
	{
		return (cols_used>7)?charge[ia]:c_dflt_charge;		// charge
	}

	template <class T>
	Ptc_Atom<T>::Ptc_s Ptc_Atom<T>::get(const size_type& ia) const
	{
		return {Z[ia], this->x[ia], this->y[ia], this->z[ia], get_sigma(ia), get_occ(ia), get_tag(ia), get_charge(ia)};
	}

	template <class T>
	void Ptc_Atom<T>::set(const size_type& ia, const Ptc_s& ptc_s)
	{
		Z[ia] = ptc_s.Z;					// atomic number

		Ptc_R_3d<T>::set(ia, ptc_s);		// xyz

		if (cols_used>4)					// standard deviation
			sigma[ia] = ptc_s.sigma;

		if (cols_used>5)					// occupancy
			occ[ia] = ptc_s.occ;

		if (cols_used>6)					// tag
			tag[ia] = ptc_s.tag;

		if (cols_used>7)					// charge
			charge[ia] = ptc_s.charge;
	}

	template <class T>
	template <class U>
	void Ptc_Atom<T>::set_ptc(const Ptc_Atom<U>& ptc, dt_bool pbc_xy, dt_bool b_statistic)
	{
		clear();
		cols_used = ptc.cols_used;
		reserve(ptc.size());

		mt::ptc_detail::set_ptc_pbc_xy(ptc, pbc_xy, b_statistic, *this);
	}

	template <class T>
	template <class U>
	void Ptc_Atom<T>::set_ptc(const pVctr_cpu_64<U>& ptc, const R_3d<U>& bs, dt_bool pbc_xy, dt_bool b_statistic)
	{
		clear();
		cols_used = ptc.s1_32();
		reserve(ptc.size());

		mt::ptc_detail::set_ptc_pbc_xy(ptc, 0, bs, pbc_xy, b_statistic, *this);
	}

	/* copy data to pointer */
	template <class T>
	template <class U>
	dt_int32 Ptc_Atom<T>::cpy_to_ptr(U *ptc, size_type n_ptc, dt_int32 is_0, dt_int32 is_e)
	{
		is_e = min(is_e, cols_used);

		if (is_0>is_e)
		{
			std::swap(is_0, is_e);
		}

		auto n_data = min(n_ptc, this->size());
		dt_int32 is = 0;

		if (fcn_chk_bound(0, is_0, is_e))
			memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), Z.data(), n_data);				// atomic number

		if (fcn_chk_bound(1, is_0, is_e))
			memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), this->x.data(), n_data);		// x-position

		if (fcn_chk_bound(2, is_0, is_e))
			memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), this->y.data(), n_data);		// y-position

		if (fcn_chk_bound(3, is_0, is_e))
			memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), this->z.data(), n_data);		// z-position

		if (fcn_chk_bound(4, is_0, is_e))
			memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), sigma.data(), n_data);			// standard deviation

		if (fcn_chk_bound(5, is_0, is_e))
			memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), occ.data(), n_data);			// occupancy

		if (fcn_chk_bound(6, is_0, is_e))
			memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), tag.data(), n_data);			// tag

		if (fcn_chk_bound(7, is_0, is_e))
			memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), charge.data(), n_data);		// charge

		return is;
	}

	// sort by x
	template <class T>
	void Ptc_Atom<T>::sort_by_x()
	{
		sort_by_idx(1);
	}

	// sort by y
	template <class T>
	void Ptc_Atom<T>::sort_by_y()
	{
		sort_by_idx(2);
	}

	// sort by z
	template <class T>
	void Ptc_Atom<T>::sort_by_z()
	{
		sort_by_idx(3);
	}

	// sort by idx
	template <class T>
	void Ptc_Atom<T>::sort_by_idx(const dt_int32 idx)
	{
		if (cols_used==4)
		{
			auto first = fcn_mkzipiter_begin(Z, this->x, this->y, this->z);
			auto last = fcn_mkzipiter_end(Z, this->x, this->y, this->z);

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
				case 3:
				{ 
					thrust::sort(first, last, mt::cgpu_fctr::less_soa<3>()); 
				} 
				break;
			}
		}
		else if (cols_used==5)
		{
			auto first = fcn_mkzipiter_begin(Z, this->x, this->y, this->z, sigma);
			auto last = fcn_mkzipiter_end(Z, this->x, this->y, this->z, sigma);

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
				case 3:
				{ 
					thrust::sort(first, last, mt::cgpu_fctr::less_soa<3>()); 
				} 
				break;
				case 4:
				{ 
					thrust::sort(first, last, mt::cgpu_fctr::less_soa<4>()); 
				} 
				break;
			}
		}
		else if (cols_used==6)
		{
			auto first = fcn_mkzipiter_begin(Z, this->x, this->y, this->z, sigma, occ);
			auto last = fcn_mkzipiter_end(Z, this->x, this->y, this->z, sigma, occ);

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
				case 3:
				{ 
					thrust::sort(first, last, mt::cgpu_fctr::less_soa<3>()); 
				} 
				break;
				case 4:
				{ 
					thrust::sort(first, last, mt::cgpu_fctr::less_soa<4>()); 
				} 
				break;				
				case 5:
				{ 
					thrust::sort(first, last, mt::cgpu_fctr::less_soa<5>()); 
				} 
				break;
			}
		}
		else if (cols_used==7)
		{
			auto first = fcn_mkzipiter_begin(Z, this->x, this->y, this->z, sigma, occ, tag);
			auto last = fcn_mkzipiter_end(Z, this->x, this->y, this->z, sigma, occ, tag);

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
				case 3:
				{ 
					thrust::sort(first, last, mt::cgpu_fctr::less_soa<3>()); 
				} 
				break;
				case 4:
				{ 
					thrust::sort(first, last, mt::cgpu_fctr::less_soa<4>()); 
				} 
				break;				
				case 5:
				{ 
					thrust::sort(first, last, mt::cgpu_fctr::less_soa<5>()); 
				} 
				break;				
				case 6:
				{ 
					thrust::sort(first, last, mt::cgpu_fctr::less_soa<6>()); 
				} 
				break;
			}
		}
		else if (cols_used==8)
		{
			auto first = fcn_mkzipiter_begin(Z, this->x, this->y, this->z, sigma, occ, tag, charge);
			auto last = fcn_mkzipiter_end(Z, this->x, this->y, this->z, sigma, occ, tag, charge);

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
				case 3:
				{ 
					thrust::sort(first, last, mt::cgpu_fctr::less_soa<3>()); 
				} 
				break;
				case 4:
				{ 
					thrust::sort(first, last, mt::cgpu_fctr::less_soa<4>()); 
				} 
				break;				
				case 5:
				{ 
					thrust::sort(first, last, mt::cgpu_fctr::less_soa<5>()); 
				} 
				break;
				case 6:
				{ 
					thrust::sort(first, last, mt::cgpu_fctr::less_soa<6>()); 
				} 
				break;
				case 7:
				{ 
					thrust::sort(first, last, mt::cgpu_fctr::less_soa<7>()); 
				} 
				break;
			}
		}
	}

	/***************************************************************************************/
	template <class T>
	void Ptc_Atom<T>::get_statistic()
	{
		if (this->empty())
		{
			return;
		}

		mt::fcn_minmax_element(Z, Z_lim.x, Z_lim.y);

		Ptc_R_3d<T>::get_statistic();	

		sigma_lim.x = c_dflt_rms3d;
		sigma_lim.y = c_dflt_rms3d;

		if (cols_used>4)
		{
			mt::fcn_minmax_element(sigma, sigma_lim.x, sigma_lim.y);
		}

		occ_lim.x = c_dflt_occ;
		occ_lim.y = c_dflt_occ;

		if (cols_used>5)
		{
			mt::fcn_minmax_element(occ, occ_lim.x, occ_lim.y);
		}

		tag_lim.x = c_dflt_tag;
		tag_lim.y = c_dflt_tag;

		if (cols_used>6)
		{
			mt::fcn_minmax_element(tag, tag_lim.x, tag_lim.y);
		}

		charge_lim.x = c_dflt_charge;
		charge_lim.y = c_dflt_charge;

		if (cols_used>7)
		{
			mt::fcn_minmax_element(tag, charge_lim.x, charge_lim.y);
		}

		this->bs.z = ::fmax(this->sz.z, this->bs.z);
	}

	// max z value within a tag
	template <class T>
	void Ptc_Atom<T>::minmax_z_by_region(const T& tag_v, T& z_min, T& z_max)
	{
		z_min = 1e20;
		z_max = -1e20;
		for(auto iz = 0; iz < this->size(); iz++)
		{
			if (tag[iz]==tag_v)
			{
				z_max = ::fmax(z_max, this->z[iz]);
				z_min = ::fmin(z_min, this->z[iz]);
			}
		} 
	}

	template <class T>
	void Ptc_Atom<T>::remove_ptc_out_z_range(const T& z_0, const T& z_e)
	{
		mt::remove_ptc_out_z_range(z_0, z_e, *this);
	}
}

/* fcns */
namespace mt
{
	template <class T>
	void remove_ptc_out_z_range(const T& z_0, const T& z_e, Ptc_Atom<T>& ptc)
	{
		dt_int32 ia_z = 0;
		for(dt_int64 ia = 0; ia < ptc.size(); ia++)
		{
			if (mt::fcn_chk_bound(ptc.z[ia], z_0, z_e))
			{
				ptc.Z[ia_z] = ptc.Z[ia];
				ptc.x[ia_z] = ptc.x[ia];
				ptc.y[ia_z] = ptc.y[ia];
				ptc.z[ia_z] = ptc.z[ia];

				if (ptc.cols_used>4)
					ptc.sigma[ia_z] = ptc.sigma[ia];

				if (ptc.cols_used>5)
					ptc.occ[ia_z] = ptc.occ[ia];

				if (ptc.cols_used>6)
					ptc.tag[ia_z] = ptc.tag[ia];

				if (ptc.cols_used>7)
					ptc.charge[ia_z] = ptc.charge[ia];

				ia_z++;
			}
		}

		ptc.resize(ia_z);
		ptc.shrink_to_fit();
	}
}