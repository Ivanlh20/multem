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

#ifndef SPEC_SLIC_H
	#define SPEC_SLIC_H

	#include "const_enum_mt.cuh"
	#include "math_mt.h"
	#include "kahan_sum.h"
	#include "r_2d.h"
	#include "types_mt.cuh"
	#include "fcns_cgpu_gen.h"
	#include "particles.cuh"
	#include "fcns_cpu.h"
	#include "cgpu_vctr.cuh"
	#include "spec_slic_in_parm.hpp"

	namespace mt
	{

		/*************************** identify planes *************************/
		class Ident_Pla
		{
			public:
				using T = dt_float64;
				
				Ident_Pla(): dv(0.1) {}

				// identify planes: require v to be sorted
				template <class U>
				Vctr_cpu<U> operator()(pVctr_cpu_32<U>& v, T v_min=0, T v_max=0)
				{
					if (v.size()==0)
					{
						return Vctr_cpu<U>();
					}

					if (v_min >= v_max)
					{
						v_min = T(v.front());
						v_max = T(v.back());
					}

					// calculate hist and correct it
					auto v_hist = hist(v, dv, v_min, v_max);

					if (v_hist.size()==1)
					{
						return Vctr_cpu<U>({fcn_mean(v)});
					}

					// calculate layer limits
					Vctr_cpu<T> v_lim;
					v_lim.reserve(v_hist.size()+1);

					for(auto iv = 0; iv < v_hist.size()-1; iv++)
					{
						if ((v_hist[iv]>0) && (v_hist[iv+1]==0))
						{
							v_lim.push_back(v_min + T(iv+1)*dv);
						}
					}
					v_hist.clear_shrink_to_fit(); // free v_hist

					v_lim.push_back(v_max+dv);
					v_lim.shrink_to_fit();

					// calculate planes
					Vctr_cpu<U> v_pln;
					v_pln.reserve(v_lim.size());

					KS<T> v_s = 0;
					dt_int32 v_c = 1;
					dt_int32 ip = 0;
					for(auto iv = 0; iv < v.size(); iv++)
					{
						const auto v_ih = v[iv];

						if (!fcn_chk_bound_eps(v_ih, v_min, v_max))
						{
							continue;
						}

						if (v_ih<v_lim[ip])
						{
							v_s += v_ih;
							v_c++;
						}
						else
						{
							v_pln.push_back(v_s/T(v_c));
							v_s = v_ih;
							v_c = 1;
							ip++;
						}
					} 

					v_pln.push_back(v_s/T(v_c));
					v_pln.shrink_to_fit();

					return v_pln;
				}
				
				// calculate planes
				template <class U>
				Vctr_cpu<U> operator()(U v_min, U v_max, U dv, eMatch_Bdr match_bdr=emb_minmax)
				{
					if (fcn_is_equal(v_min, v_max))
					{
						return Vctr_cpu<U>({U(0.5)*(v_min+v_max)});
					}

					if (v_min >= v_max)
					{
						return Vctr_cpu<U>();
					}

					const auto v_l = T(v_max-v_min);

					const dt_int32 nv = max(0, fcn_cceil(v_l/T(dv)))+1;

					// calculate planes
					Vctr_cpu<U> v_pln;
					v_pln.reserve(nv);

					switch (match_bdr)
					{
						case emb_min:
						{
							for(auto iv=0; iv<nv; iv++)
							{
								v_pln.push_back(v_min + U(iv)*dv);
							}
						}
						break;
						case emb_max:
						{
							for(auto iv=nv-1; iv>=0; iv--)
							{
								v_pln.push_back(v_max - U(iv)*dv);
							}
						}
						break;
						case emb_minmax:
						{
							const auto dv_b = dv + 0.5*(v_l-(nv-1)*dv);

							v_pln.push_back(v_min);
							for(auto iv=1; iv<nv; iv++)
							{
								const auto dv_t = ((iv==1)||(iv==nv-1))?dv_b:dv;
								v_pln.push_back(v_pln[iv-1] + dv_t);
							}
						}
						break;
					}

					return v_pln;
				}

			private:
				T dv;

				// calculate corrected histogram
				template <class U>
				Vctr_cpu<dt_int32> hist(pVctr_cpu_32<U>& v, const T& dv, const T& v_min, const T& v_max)
				{
					const auto v_l = ::fmax(v_max-v_min, dv);

					const auto n_bins = fcn_cceil<dt_int32>(v_l/dv);

					auto fcn_v_r = [v_min, dv](const dt_int32& ik, T& v_r_0, T& v_r_e)
					{
						 v_r_0 = v_min + T(ik-1)*dv;
						 v_r_e = v_r_0 + dv;
					};

					Vctr_cpu<dt_int32> v_hist(n_bins, 0);

					// get histogram
					for(auto iv = 0; iv< v.size(); iv++)
					{
						T v_r_0, v_r_e;
						const auto v_id = T(v[iv]);
						auto ih = fcn_cfloor<dt_int32>((v_id-v_min)/dv);

						fcn_v_r(ih-1, v_r_0, v_r_e);

						if (v_id<v_r_0)
						{
							for(auto ik = ih; ik >= 0; ik--)
							{
								fcn_v_r(ik-1, v_r_0, v_r_e);
								if (fcn_chk_bound(v_id, v_r_0, v_r_e))
								{
									ih = ik-1;
									break;
								}
							}
						}
						else if (v_id>v_r_e)
						{
							for(auto ik = ih; ik < n_bins; ik++)
							{
								fcn_v_r(ik, v_r_0, v_r_e);
								if (fcn_chk_bound(v_id, v_r_0, v_r_e))
								{
									ih = ik;
									break;
								}
							}
						}
						ih = fcn_set_bound(ih, 0, n_bins-1);

						v_hist[ih]++;
					}

					return v_hist;
				}
		};

		/************************* slice thickness ***************************/
		template <class T>
		class Out_Slic
		{
			public:
				Out_Slic(): z_pln(0), z_lim(0, 0), z_int_lim(0, 0), 
				iatom_lim(0, 0), ithk(-1) {}

				R_2d<T> z_pln;				// z-position
				R_2d<T> z_lim;				// z-position limit
				R_2d<T> z_int_lim;			// z-int-position limit
				R_2d<dt_int32> iatom_lim;	// index to z-position limit
				dt_int32 ithk;				// thick index

				T sli_thick() const { return ::fabs(z_lim.y-z_lim.x); }
		};

		template <class T>
		class Spec_Slic
		{
			public:
				using value_type = T;
				using size_type = dt_uint64;

				Vctr_cpu<T> z_plns;
				Vctr_cpu<Out_Slic<T>> slice;
				Vctr_cpu<Thick<T>> thick;

				Spec_Slic(): patoms_r(nullptr), patoms(nullptr), z_eps(1e-3) {}

				void set_in_data(const Vctr_cpu<Spec_Slic_In_Parm>& in_slic, Ptc_Atom<T>* patoms_r=nullptr, Ptc_Atom<T>* patoms=nullptr)
				{
					patoms_r = patoms_r;
					if (fcn_is_null_ptr(patoms))
					{
						patoms = patoms_r;
					}

					z_plns = get_z_plns(in_slic, patoms_r);
				}

				Vctr_cpu<T> get_z_plns(const Vctr_cpu<Spec_Slic_In_Parm>& in_slic, Ptc_Atom<T>* patoms)
				{
					// get planes
					Vctr_cpu<Vctr_cpu<T>> z_pln_s(in_slic.size());

					for(auto ir=0; ir<in_slic.size(); ir++)
					{
						z_pln_s[ir] = get_z_plns(in_slic[ir], patoms);
					}

					// planes fusion
					Vctr_cpu<Vctr_cpu<T>> z_pln(in_slic.size());

					for(auto ir=0; ir<z_pln_s.size(); ir++)
					{
						z_pln.push_back(get_z_plns(in_slic[ir], patoms));
					}

					return z_pln;
				}

				//void match_thickness(eSpec_Slic_Typ spec_slic_typ, Ptc_Atom<T>& patoms, 
				//eSim_Thick_Typ thick_type, Vctr_cpu<T> &thick)
				//{
				//	if (thick_type == estt_whole_spec)
				//	{
				//		thick.resize(1);
				//		thick[0] = patoms->z_max;
				//		return;
				//	}

				//	auto z_plns = get_z_plane(spec_slic_typ, patoms);

				//	if (thick_type == estt_through_slices)
				//	{
				//		z_plns = get_z_slice(spec_slic_typ, z_plns, patoms);
				//	}

				//	fcn_match_vctr(z_plns.begin(), z_plns.end(), thick, 0.25*patoms->sli_thick);
				//}

				//T sli_thick(dt_int32 islice_0, dt_int32 islice_e)
				//{
				//	return (islice_e<slice.size())?(slice[islice_e].z_e - slice[islice_0].z_0):0.0;
				//}

				//T sli_thick(dt_int32 islice)
				//{
				//	return sli_thick(islice, islice);
				//}

				//T z_m(dt_int32 islice)
				//{
				//	return (islice<slice.size())?slice[islice].z_m():0.0;
				//}

				//T dz_m(dt_int32 islice_0, dt_int32 islice_e)
				//{
				//	return fabs(z_m(islice_e) - z_m(islice_0));
				//}

				//void calculate()
				//{
				//	m_z_slice = get_z_slice(m_multem_in_parm->spec_slic_typ, z_plns, *patoms);

				//	thick = get_thick(m_multem_in_parm, m_z_slice, *patoms_r);

				//	slice = get_slicing(m_multem_in_parm, m_z_slice, thick, *patoms);
				//}

			private:
				const T z_eps;

				Ptc_Atom<T>* patoms_r;
				Ptc_Atom<T>* patoms;

				Vctr_cpu<T> m_z_slice;
				Ident_Pla ident_pln;

				// get z positions by tag
				Vctr_cpu<T> get_z_pos_by_tag(Ptc_Atom<T>* patoms, const dt_int32& tag)
				{
					Vctr_cpu<T> z;
					z.reserve(patoms->size());

					for(auto iz = 0; iz<patoms->size(); iz++)
					{
						if (patoms->get_tag(iz)==tag)
						{
							z.push_back(patoms->z[iz]);
						}
					}
					z.shrink_to_fit();
					std::sort(z.begin(), z.end());

					return z;
				}

				// get z positions by z range and atomic number
				Vctr_cpu<T> get_z_pos_by_rng(Ptc_Atom<T>* patoms, const dt_int32& Z, const R_2d<T>& z_lim)
				{
					Vctr_cpu<T> z;
					z.reserve(patoms->size());

					if Z<=0
					{
						for(auto iz = 0; iz<patoms->size(); iz++)
						{
							if (fcn_chk_bound_eps(patoms->z[iz], z_lim.x, z_lim.y))
							{
								z.push_back(patoms->z[iz]);
							}
						}
					}
					else
					{
						for(auto iz = 0; iz<patoms->size(); iz++)
						{
							if ((patoms->Z[iz]==Z) && fcn_chk_bound_eps(patoms->z[iz], z_lim.x, z_lim.y))
							{
								z.push_back(patoms->z[iz]);
							}
						}
					}

					z.shrink_to_fit();
					std::sort(z.begin(), z.end());

					return z;
				}

				// identify planes: atoms have to be sorted along z
				Vctr_cpu<T> get_z_plns(const Spec_Slic_In_Parm& in_slic, Ptc_Atom<T>* patoms)
				{
					if (patoms->size() == 0)
					{
						return Vctr_cpu<T>();
					}

					// get z planes
					if (in_slic.is_spec_slic_by_user_def())
					{
						return in_slic.z_plns;
					}
					else
					{
						Vctr_cpu<T> z;

						if (in_slic.is_spec_slic_sel_typ_by_tag())
						{
							z = get_z_pos_by_tag(patoms, in_slic.sel_tag);
						}
						else if (in_slic.is_spec_slic_sel_typ_by_z())
						{
							z = get_z_pos_by_rng(patoms, in_slic.sel_Z, in_slic.sel_z_lim);
						}

						if (in_slic.is_spec_slic_by_planes())
						{
							return ident_pln(z);
						}
						else
						{
							return ident_pln(z_ct_min, z_ct_max, in_slic.sli_thick);
						}
					}

					return Vctr_cpu<T>();
				}


				//// get spacing
				//T get_spacing(size_type ix, Vctr_cpu<T> &x, T sli_thick=0.5)
				//{
				//	if (x.size()==1)
				//	{
				//		return sli_thick;
				//	}

				//	ix = (ix <= 0)?1:min(ix, x.size()-1);
				//	return (x.size()>1)?x[ix]-x[ix-1]:0.0;
				//}

				//// get z slicing
				//Vctr_cpu<T> get_z_slice(eSpec_Slic_Typ spec_slic_typ, Vctr_cpu<T> &z_plns, Ptc_Atom<T>& patoms)
				//{
				//	Vctr_cpu<T> z_slice;

				//	if ((spec_slic_typ!=esst_dz_sub) && (z_plns.size() == 1))
				//	{
				//		z_slice.resize(2);
				//		z_slice[0] = patoms->z_int_min;
				//		z_slice[1] = patoms->z_int_max;
				//		z_slice.shrink_to_fit();

				//		return z_slice;
				//	}

				//	Vctr_cpu<T> z_plane_sli = z_plns;

				//	z_slice.resize(z_plane_sli.size()+1);
				//	z_slice[0] = z_plane_sli[0]-0.5*get_spacing(0, z_plane_sli, patoms->sli_thick);
				//	for(auto iz=1; iz<z_slice.size()-1; iz++)
				//	{
				//		z_slice[iz] = 0.5*(z_plane_sli[iz]+z_plane_sli[iz-1]);
				//	}
				//	z_slice[z_slice.size()-1] = z_slice[z_slice.size()-2]+get_spacing(z_plane_sli.size(), z_plane_sli, patoms->sli_thick);

				//	if (spec_slic_typ==esst_dz_sub)
				//	{
				//		T dz_b = get_spacing(1, z_plane_sli, patoms->sli_thick);
				//		if (patoms->z_int_min<z_slice.front()-dz_b)
				//		{
				//			T dz_s = get_spacing(2, z_plane_sli, patoms->sli_thick);
				//			auto z_slice_top = ident_pln(patoms->z_int_min, z_slice.front()-dz_b, dz_s);
				//			z_slice.insert(z_slice.begin(), z_slice_top.begin(), z_slice_top.end());
				//		}
				//		else
				//		{
				//			z_slice[0] = patoms->z_int_min;
				//		}

				//		dz_b = get_spacing(z_plane_sli.size()-1, z_plane_sli, patoms->sli_thick);
				//		if (z_slice.back()+dz_b<patoms->z_int_max)
				//		{
				//			T dz_s = get_spacing(z_plane_sli.size()-2, z_plane_sli, patoms->sli_thick);
				//			auto z_slice_bottom = ident_pln(z_slice.back()+dz_b, patoms->z_int_max, dz_s);
				//			z_slice.insert(z_slice.end(), z_slice_bottom.begin(), z_slice_bottom.end());
				//		}
				//		else
				//		{
				//			z_slice[z_slice.size()-1] = patoms->z_int_max;
				//		}
				//	}

				//	z_slice.shrink_to_fit();

				//	return z_slice;
				//}

				//// get thick
				//Vctr<Thick<T>, edev_cpu> get_thick(Multem_In_Parm<T> *multem_in_parm, 
				//Vctr_cpu<T> &z_slice, Ptc_Atom<T>& patoms)
				//{
				//	const auto thick_type = multem_in_parm->thick_type;

				//	auto get_islice = [thick_type](Vctr_cpu<T> &z_slice, const T& z)->dt_int32
				//	{
				//		if (thick_type==estt_through_slices)
				//		{
				//			for(auto i = 0; i<z_slice.size()-1; i++)
				//			{
				//				if (fabs(z_slice[i+1]-z)<Epsilon<dt_float32>::rel)
				//				{
				//					return i;
				//				}
				//			}
				//			return 0;
				//		}
				//		else
				//		{
				//			for(auto i = 0; i<z_slice.size()-1; i++)
				//			{
				//				if ((z_slice[i] < z) && (z <= z_slice[i+1]))
				//				{
				//					return i;
				//				}
				//			}
				//			return 0;
				//		}
				//	};

				//	auto b_sws = multem_in_parm->is_spec_slic_by_dz_sub_whole_spec();

				//	Vctr<Thick<T>, edev_cpu> thick(multem_in_parm->thick.size());
				//	for(auto ik = 0; ik<thick.size(); ik++)
				//	{
				//		thick[ik].z = multem_in_parm->thick[ik];
				//		auto islice = (b_sws)?(z_slice.size()-2):get_islice(z_slice, thick[ik].z);
				//		thick[ik].islice = islice;

				//		auto iatom_e = fd_by_z(patoms->z, z_slice[islice+1], false);
				//		thick[ik].iatom_e = iatom_e;

				//		thick[ik].z_zero_def_plane = multem_in_parm->obj_lens.get_zero_def_plane(patoms->z[0], patoms->z[iatom_e]);
				//		if (multem_in_parm->is_sim_through_slices())
				//		{
				//			thick[ik].z_back_prop = 0;
				//		}
				//		else
				//		{
				//			// I need to recheck this part, the average should be replace by the plane
				//			thick[ik].z_back_prop = thick[ik].z_zero_def_plane - 0.5*(z_slice[islice]+z_slice[islice+1]);
				//			if (fabs(thick[ik].z_back_prop)<Epsilon<dt_float32>::rel)
				//			{
				//				thick[ik].z_back_prop = 0;
				//			}
				//		}
				//	}

				//	if (!multem_in_parm->is_multislice())
				//	{
				//		for(auto ithick = 0; ithick<thick.size(); ithick++)
				//		{
				//			thick[ithick].islice = ithick;
				//			thick[ithick].z_back_prop = thick[ithick].z_zero_def_plane - thick[ithick].z;
				//		}
				//	}

				//	return thick;
				//}

				//// slicing
				//Vctr<Out_Slic<T>, edev_cpu> get_slicing(Multem_In_Parm<T> *multem_in_parm, 
				//Vctr_cpu<T> &z_slice, Vctr<Thick<T>, edev_cpu>& thick, Ptc_Atom<T>& patoms)
				//{	
				//	if (!multem_in_parm->is_multislice())
				//	{
				//		Vctr<Out_Slic<T>, edev_cpu> slice(thick.size());
				//		for(auto islice = 0; islice<slice.size(); islice++)
				//		{
				//			slice[islice].z_0 = (islice == 0)?patoms->z_int_min:thick[islice-1].z;
				//			slice[islice].z_e = (thick.size() == 1)?patoms->z_int_max:thick[islice].z;
				//			slice[islice].z_int_0 = slice[islice].z_0;
				//			slice[islice].z_int_e = slice[islice].z_e;
				//			slice[islice].iatom_0 = (islice == 0)?0:(thick[islice-1].iatom_e+1);
				//			slice[islice].iatom_e = thick[islice].iatom_e;
				//			slice[islice] .ithk = islice;
				//		}
				//		return slice;
				//	}

				//	Vctr<Out_Slic<T>, edev_cpu> slice(z_slice.size()-1);
				//	for(auto islice = 0; islice<slice.size(); islice++)
				//	{
				//		dt_bool Inc_Borders = false;
				//		slice[islice].z_0 = z_slice[islice];
				//		slice[islice].z_e = z_slice[islice + 1];
				//		switch(multem_in_parm->spec_slic_typ)
				//		{
				//			case  esst_plns_proj:
				//			{
				//				slice[islice].z_int_0 = slice[islice].z_0;
				//				slice[islice].z_int_e = slice[islice].z_e;
				//				Inc_Borders = false;
				//			}
				//			break;
				//			case esst_dz_proj:
				//			{
				//				slice[islice].z_int_0 = slice[islice].z_0;
				//				slice[islice].z_int_e = slice[islice].z_e;
				//				Inc_Borders = false;
				//			}
				//			break;
				//			case esst_dz_sub:
				//			{
				//				T z_m = slice[islice].z_m();

				//				slice[islice].z_int_0 = ::fmin(z_m - patoms->R_int_max, slice[islice].z_0);
				//				slice[islice].z_int_e = ::fmax(z_m + patoms->R_int_max, slice[islice].z_e);
				//				Inc_Borders = true;
				//			}
				//			break;
				//		}

				//		fd_by_z(patoms->z, slice[islice].z_int_0, slice[islice].z_int_e, slice[islice].iatom_0, slice[islice].iatom_e, Inc_Borders);
				//	
				//		slice[islice].ithk = -1;
				//	}

				//	// set thickness index
				//	for(auto ithk = 0; ithk<thick.size(); ithk++)
				//	{
				//		slice[thick[ithk].islice].ithk = ithk;
				//	}

				//	return slice;
				//}

				//// find patoms in slice
				//dt_int32 fd_by_z(Vctr_cpu<T> &z, T z_e, dt_bool Inc_Borders)
				//{
				//	dt_int32 iz_e =-1;
				//	z_e = (Inc_Borders)?z_e+Epsilon<T>::rel:z_e;

				//	if (z_e<z.front())
				//	{
				//		return iz_e;
				//	}

				//	iz_e = (z.back() <= z_e)?(z.size()-1):(std::lower_bound(z.begin(), z.end(), z_e)-z.begin()-1);

				//	if (z[iz_e] < z.front())
				//	{
				//		iz_e = -1;
				//	}

				//	return iz_e;
				//}

				//// find patoms in slice
				//void fd_by_z(Vctr_cpu<T> &z, T z_0, T z_e, dt_int32& iz_0, dt_int32& iz_e, dt_bool Inc_Borders)
				//{
				//	z_0 = (Inc_Borders)?(z_0-Epsilon<T>::rel):z_0;
				//	z_e = (Inc_Borders)?(z_e+Epsilon<T>::rel):z_e;

				//	if ((z_0>z_e)||(z.back()<z_0)||(z_e<z.front()))
				//	{
				//		iz_0 = 1;
				//		iz_e = 0;
				//		return;
				//	}

				//	iz_0 = (z_0 <= z.front())?0:(std::lower_bound(z.begin(), z.end(), z_0)-z.begin());
				//	iz_e = (z.back() <= z_e)?(z.size()-1):(std::lower_bound(z.begin(), z.end(), z_e)-z.begin()-1);

				//	if ((iz_0 > iz_e)||(z[iz_e] < z_0)||(z_e < z[iz_0]))
				//	{
				//		iz_0 = 1;
				//		iz_e = 0;
				//	}
				//}
		};

	}

#endif
