/*
 * This file is part of MULTEM.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * MULTEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MULTEM is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef ATOMIC_DATA_MT_H
#define ATOMIC_DATA_MT_H

#ifdef _MSC_VER
#pragma once
#endif// _MSC_VER

#include <numeric>
#include <algorithm>

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "lin_alg_def.cuh"
#include <atom_data_api.h>
#include <thrust/sort.h>

namespace mt
{
	namespace host_device_detail
	{
		template <class T>
		DEVICE_CALLABLE FORCE_INLINE 
		void kh_sum(T &sum_v, T v, T &error);
	}

	template <class TVector>
	void minmax_element(TVector &x, Value_type<TVector> &x_min, Value_type<TVector> &x_max);

	template <class TVector>
	Value_type<TVector> mean(TVector &M_i);

	template <class TVector>
	void mean_std(TVector &M_i, Value_type<TVector> &x_mean, Value_type_r<TVector> &x_std);

	template<class TVector>
	void scale(Value_type<TVector> f, TVector &x);

	template <class T>
	class Atom_Data;
	
  template <class T>
	void rotate_atoms(Atom_Data<T> &atoms, T theta, r3d<T> u0, r3d<T> p0)
	{
		const auto Rm = get_rotation_matrix(theta, u0);

		for(int iatoms = 0; iatoms<atoms.size(); iatoms++)
		{

			auto r = atoms.to_r3d(iatoms).rotate(Rm, p0);

			atoms.x[iatoms] = r.x;
			atoms.y[iatoms] = r.y;
			atoms.z[iatoms] = r.z;
		}
	}

	template <class T>
	void remove_atoms_outside_z_range(Atom_Data<T> &atoms, T z_0, T z_e)
	{
		int iatoms_z = 0;
		for(int iatoms = 0; iatoms<atoms.size(); iatoms++)
		{
			const auto z = atoms.z[iatoms];
			if((z_0<z) && (z<z_e))
			{
				atoms.Z[iatoms_z] = atoms.Z[iatoms];
				atoms.x[iatoms_z] = atoms.x[iatoms];
				atoms.y[iatoms_z] = atoms.y[iatoms];
				atoms.z[iatoms_z] = atoms.z[iatoms];
				atoms.sigma[iatoms_z] = atoms.sigma[iatoms];
				atoms.occ[iatoms_z] = atoms.occ[iatoms];
				atoms.region[iatoms_z] = atoms.region[iatoms];
				atoms.charge[iatoms_z] = atoms.charge[iatoms];

				iatoms_z++;
			}
		}

		atoms.resize(iatoms_z);
		atoms.shrink_to_fit();
	}


  template <typename T>
  struct Atom_Data<T>::sort_atoms_by_z
  {
    template <class Ttuple1, class Ttuple2>
    DEVICE_CALLABLE
    bool operator()(const Ttuple1 &t1, const Ttuple2 &t2)
    {
      return thrust::get<3>(t1) < thrust::get<3>(t2);
    }
  };

  template <typename T>
  void Atom_Data<T>::get_statistic(std::vector<Atom_Type<T, e_host>> *atom_type_ptr)
  {
    if(empty())
    {
      return;
    }

    get_Z_unique();

    auto Z_min_max = std::minmax_element(Z.begin(), Z.end(), [](const int &a, const int &b){ return (a % 1000)<(b % 1000); });
    Z_min = *(Z_min_max.first);
    Z_max = *(Z_min_max.second);

    mt::minmax_element(x, x_min, x_max);
    mt::minmax_element(y, y_min, y_max);
    mt::minmax_element(z, z_min, z_max);
    mt::minmax_element(sigma, sigma_min, sigma_max);
    mt::minmax_element(occ, occ_min, occ_max);
    mt::minmax_element(region, region_min, region_max);

    mt::mean_std(x, x_mean, x_std);
    mt::mean_std(y, y_mean, y_std);
    mt::mean_std(z, z_mean, z_std);

    bool bAtomTypes = (atom_type_ptr == nullptr)?false:true;
    R_int_min = R_int_max = 2.5;
    if(bAtomTypes)
    {
      R_int_min = R_int_max = (*atom_type_ptr)[get_Z(Z[0])-1].coef[0].R_max;
      for(auto iatoms = 0; iatoms < size(); iatoms++)
      {
        R_int_min = min((*atom_type_ptr)[get_Z(Z[iatoms])-1].coef[0].R_max, R_int_min);
        R_int_max = max((*atom_type_ptr)[get_Z(Z[iatoms])-1].coef[0].R_max, R_int_max);
      }
    }

    s_x = x_max - x_min;
    s_y = y_max - y_min;
    s_z = z_max - z_min;

    x_int_min = x_min - R_int_max;
    x_int_max = x_max + R_int_max;

    y_int_min = y_min - R_int_max;
    y_int_max = y_max + R_int_max;

    z_int_min = z_min - R_int_max;
    z_int_max = z_max + R_int_max;

    s_x_int = x_int_max - x_int_min;
    s_y_int = y_int_max - y_int_min;
    s_z_int = z_int_max - z_int_min;

    if(isZero(l_x))
    {
      l_x = s_x;
    }

    if(isZero(l_y))
    {
      l_y = s_y;
    }

    l_z = ::fmax(s_z, l_z);

    l_x_int = l_x + 2.0*R_int_max;
    l_y_int = l_y + 2.0*R_int_max;
    l_z_int = l_z + 2.0*R_int_max;

    if(isZero(ct_a))
    {
      ct_a = l_x;
    }

    if(isZero(ct_b))
    {
      ct_b = l_y;
    }

    if(isZero(ct_c))
    {
      ct_c = l_z;
    }
  }
  
  template <typename T>
  void Atom_Data<T>::get_statistic() {
    get_statistic(nullptr);
  }

  // Sort atoms along z-axis.
  template <typename T>
  void Atom_Data<T>::sort_by_z()
  {
    // if(!thrust::is_sorted(z.begin(), z.end()));
    auto first = thrust::make_zip_iterator(thrust::make_tuple(Z.begin(), x.begin(), y.begin(), z.begin(), sigma.begin(), occ.begin(), region.begin(), charge.begin()));
    auto last = thrust::make_zip_iterator(thrust::make_tuple(Z.end(), x.end(), y.end(), z.end(), sigma.end(), occ.end(), region.end(), charge.end()));

    thrust::sort(first, last, sort_atoms_by_z());
  }

  // max z value within a region
  template <typename T>
  void Atom_Data<T>::minmax_z_by_region(T region_v, T &zr_min, T &zr_max)
  {
    zr_min = 1e20;
    zr_max = -1e20;
    for(auto iz=0; iz<size(); iz++)
    {
      if(region[iz]==region_v)
      {
        zr_max = ::fmax(zr_max, z[iz]);
        zr_min = ::fmin(zr_min, z[iz]);
      }
    } 
  }

  template <typename T>
  void Atom_Data<T>::validate_amorphous_parameters()
  {
    // delete not valid regions
    int il = 0;
    while(il<amorp_lay_info.size())
    {
      if(amorp_lay_info[il].lz()<1e-4)
      {
        amorp_lay_info.erase(amorp_lay_info.begin()+il);
      }
      else
      {
        amorp_lay_info[il].set_region(*this);
        il++;
      }
    }

    // sort amorphous layers
    auto sort_z_0 = [](Amorp_Lay_Info<T> &v1, Amorp_Lay_Info<T> &v2){ return v1.z_0<v2.z_0; };
    std::sort(amorp_lay_info.begin(), amorp_lay_info.end(), sort_z_0);

    // set amorphous layer type
    if(amorp_lay_info.size()>0)
    {
      T zc_min, zc_max;
      int region_c = 0;
      minmax_z_by_region(region_c, zc_min, zc_max);

      for(auto il=0; il<amorp_lay_info.size(); il++)
      {
        if(isZero<T>(amorp_lay_info[il].dz))
        {
          amorp_lay_info[il].dz = 2.0;
        }

        auto z_0 = amorp_lay_info[il].z_0;
        if(z_0<zc_min)
        {
          amorp_lay_info[il].type = eALT_Top;
        }
        else if(z_0<zc_max)
        {
          amorp_lay_info[il].type = eALT_Middle;
        }
        else
        {
          amorp_lay_info[il].type = eALT_Bottom;
        }
      }
    }
  }


} // namespace mt

#endif
