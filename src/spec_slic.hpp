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

#ifndef SPEC_SLIC_H
	#define SPEC_SLIC_H

	#include "const_enum.cuh"
	#include "math.cuh"
	#include "particles.cuh"
	#include "cpu_fcns.hpp"
	#include "cgpu_vctr.cuh"

	namespace mt
	{

		/************************* slice thickness ***************************/
		template <class T>
		struct Slice
		{
			Slice(): z_0(0), z_e(0), z_int_0(0), 
			z_int_e(0), iatom_0(1), iatom_e(0), ithk(-1) {}

			T z_0;					// Initial z-position
			T z_e;					// Final z-position
			T z_int_0;				// Initial z-position
			T z_int_e;				// Final z-position
			dt_int32 iatom_0;		// Index to initial z-position
			dt_int32 iatom_e;		// Index to final z-position
			dt_int32 ithk;			// thick index

			T sli_thk() const { return fabs(z_e-z_0); }

			T z_m() const { return 0.5*(z_e+z_0); }
		};

		template <class T>
		class Spec_Slic
		{
			public:
				using value_type = T;
				using size_type = dt_uint64;
				using TVctr_r = Vctr_cpu<T>;
				using TVctr_i = Vctr_cpu<dt_int32>;

				Spec_Slic(): patoms_r(nullptr), patoms(nullptr), z_eps(1e-3) {}

				void set_in_data(Ptc_Atom<T>* patoms_r, const eSpec_Slic_Typ& pot_slic_typ, Ptc_Atom<T>* patoms=nullptr)
				{
					patoms_r = patoms_r;
					patoms = (fcn_is_null_ptr(patoms))?patoms_r:patoms;

					z_plane = get_z_plane(m_in_multem->pot_slic_typ, *patoms_r);
				}

				void match_thickness(eSpec_Slic_Typ pot_slic_typ, Ptc_Atom<T>& patoms, 
				eSim_Thick_Typ thick_type, TVctr_r &thick)
				{
					if (thick_type == estt_whole_spec)
					{
						thick.resize(1);
						thick[0] = patoms.z_max;
						return;
					}

					auto z_plane = get_z_plane(pot_slic_typ, patoms);

					if (thick_type == estt_through_slices)
					{
						z_plane = get_z_slice(pot_slic_typ, z_plane, patoms);
					}

					fcn_match_vctr(z_plane.begin(), z_plane.end(), thick, 0.25*patoms.sli_thk);
				}

				T sli_thk(dt_int32 islice_0, dt_int32 islice_e)
				{
					return (islice_e<slice.size())?(slice[islice_e].z_e - slice[islice_0].z_0):0.0;
				}

				T sli_thk(dt_int32 islice)
				{
					return sli_thk(islice, islice);
				}

				T z_m(dt_int32 islice)
				{
					return (islice<slice.size())?slice[islice].z_m():0.0;
				}

				T dz_m(dt_int32 islice_0, dt_int32 islice_e)
				{
					return fabs(z_m(islice_e) - z_m(islice_0));
				}

				void calculate()
				{
					m_z_slice = get_z_slice(m_in_multem->pot_slic_typ, z_plane, *patoms);

					thick = get_thick(m_in_multem, m_z_slice, *patoms_r);

					slice = get_slicing(m_in_multem, m_z_slice, thick, *patoms);
				}

				TVctr_r z_plane;
				Vctr<Slice<T>, edev_cpu> slice;
				Vctr<Thick<T>, edev_cpu> thick;

			private:
				const T z_eps;

				In_Multem<T> *m_in_multem;

				Ptc_Atom<T>* patoms_r;
				Ptc_Atom<T>* patoms;

				TVctr_r m_z_slice;
				Identify_Planes<T> identify_planes;

				// get spacing
				T get_spacing(size_type ix, TVctr_r &x, T sli_thk=0.5)
				{
					if (x.size()==1)
					{
						return sli_thk;
					}

					ix = (ix <= 0)?1:min(ix, x.size()-1);
					return (x.size()>1)?x[ix]-x[ix-1]:0.0;
				}

				// Identify planes: Require that the patoms to be sorted along z
				TVctr_r get_z_plane(eSpec_Slic_Typ pot_slic_typ, Ptc_Atom<T>& patoms)
				{
					TVctr_r z_plane;

					if (patoms.size() == 0)
					{
						return z_plane;
					}

					patoms.validate_amorphous_parameters();

					const dt_int32 region_ct = 0;

					// select z values of the xtl region
					TVctr_r z_ct;
					z_ct.reserve(patoms.size());

					for(auto iz = 0; iz<patoms.size(); iz++)
					{
						if (patoms.region[iz]==region_ct)
						{
							z_ct.push_back(patoms.z[iz]);
						}
					}
					std::sort(z_ct.begin(), z_ct.end());

					// calculate z xtl limits
					T z_ct_min = z_ct.front();
					T z_ct_max = z_ct.back();

					if (pot_slic_typ==esst_planes_proj)
					{
						z_plane = identify_planes(z_ct);
					}
					else
					{
						z_plane = identify_planes(z_ct_min, z_ct_max, patoms.sli_thk);
					}

					std::vector<dt_float64> zt2(z_plane.begin(), z_plane.end());
					auto bb_ali = !patoms.spec_lay_info.empty();
					if (bb_ali && (fabs(patoms.z_min-z_plane.front())>z_eps))
					{
						T dz_b = get_spacing(1, z_plane);
						auto amorp = patoms.spec_lay_info.front();
						auto z_plane_top = identify_planes(amorp.z_0, amorp.z_e-dz_b, amorp.sli_thk);
						z_plane.insert(z_plane.begin(), z_plane_top.begin(), z_plane_top.end());
					}

					if (bb_ali && (fabs(z_plane.back()-patoms.z_max)>z_eps))
					{
						T dz_b = get_spacing(z_plane.size()-1, z_plane);
						auto amorp = patoms.spec_lay_info.back();
						auto z_plane_bottom = identify_planes(amorp.z_0+dz_b, amorp.z_e, amorp.sli_thk);
						z_plane.insert(z_plane.end(), z_plane_bottom.begin(), z_plane_bottom.end());
					}

					// unique vector
					fcn_uniq_vctr(z_plane);

					return z_plane;
				}

				// get z slicing
				TVctr_r get_z_slice(eSpec_Slic_Typ pot_slic_typ, TVctr_r &z_plane, Ptc_Atom<T>& patoms)
				{
					TVctr_r z_slice;

					if ((pot_slic_typ!=esst_dz_sub) && (z_plane.size() == 1))
					{
						z_slice.resize(2);
						z_slice[0] = patoms.z_int_min;
						z_slice[1] = patoms.z_int_max;
						z_slice.shrink_to_fit();

						return z_slice;
					}

					TVctr_r z_plane_sli = z_plane;

					z_slice.resize(z_plane_sli.size()+1);
					z_slice[0] = z_plane_sli[0]-0.5*get_spacing(0, z_plane_sli, patoms.sli_thk);
					for(auto iz=1; iz<z_slice.size()-1; iz++)
					{
						z_slice[iz] = 0.5*(z_plane_sli[iz]+z_plane_sli[iz-1]);
					}
					z_slice[z_slice.size()-1] = z_slice[z_slice.size()-2]+get_spacing(z_plane_sli.size(), z_plane_sli, patoms.sli_thk);

					if (pot_slic_typ==esst_dz_sub)
					{
						T dz_b = get_spacing(1, z_plane_sli, patoms.sli_thk);
						if (patoms.z_int_min<z_slice.front()-dz_b)
						{
							T dz_s = get_spacing(2, z_plane_sli, patoms.sli_thk);
							auto z_slice_top = identify_planes(patoms.z_int_min, z_slice.front()-dz_b, dz_s);
							z_slice.insert(z_slice.begin(), z_slice_top.begin(), z_slice_top.end());
						}
						else
						{
							z_slice[0] = patoms.z_int_min;
						}

						dz_b = get_spacing(z_plane_sli.size()-1, z_plane_sli, patoms.sli_thk);
						if (z_slice.back()+dz_b<patoms.z_int_max)
						{
							T dz_s = get_spacing(z_plane_sli.size()-2, z_plane_sli, patoms.sli_thk);
							auto z_slice_bottom = identify_planes(z_slice.back()+dz_b, patoms.z_int_max, dz_s);
							z_slice.insert(z_slice.end(), z_slice_bottom.begin(), z_slice_bottom.end());
						}
						else
						{
							z_slice[z_slice.size()-1] = patoms.z_int_max;
						}
					}

					z_slice.shrink_to_fit();

					return z_slice;
				}

				// get thick
				Vctr<Thick<T>, edev_cpu> get_thick(In_Multem<T> *in_multem, 
				TVctr_r &z_slice, Ptc_Atom<T>& patoms)
				{
					const auto thick_type = in_multem->thick_type;

					auto get_islice = [thick_type](TVctr_r &z_slice, const T& z)->dt_int32
					{
						if (thick_type==estt_through_slices)
						{
							for(auto i = 0; i<z_slice.size()-1; i++)
							{
								if (fabs(z_slice[i+1]-z)<Epsilon<dt_float32>::rel)
								{
									return i;
								}
							}
							return 0;
						}
						else
						{
							for(auto i = 0; i<z_slice.size()-1; i++)
							{
								if ((z_slice[i] < z) && (z <= z_slice[i+1]))
								{
									return i;
								}
							}
							return 0;
						}
					};

					auto b_sws = in_multem->is_spec_slic_dz_sub_whole_spec();

					Vctr<Thick<T>, edev_cpu> thick(in_multem->thick.size());
					for(auto ik = 0; ik<thick.size(); ik++)
					{
						thick[ik].z = in_multem->thick[ik];
						auto islice = (b_sws)?(z_slice.size()-2):get_islice(z_slice, thick[ik].z);
						thick[ik].islice = islice;

						auto iatom_e = fd_by_z(patoms.z, z_slice[islice+1], false);
						thick[ik].iatom_e = iatom_e;

						thick[ik].z_zero_def_plane = in_multem->obj_lens.get_zero_def_plane(patoms.z[0], patoms.z[iatom_e]);
						if (in_multem->is_sim_through_slices())
						{
							thick[ik].z_back_prop = 0;
						}
						else
						{
							// I need to recheck this part, the average should be replace by the plane
							thick[ik].z_back_prop = thick[ik].z_zero_def_plane - 0.5*(z_slice[islice]+z_slice[islice+1]);
							if (fabs(thick[ik].z_back_prop)<Epsilon<dt_float32>::rel)
							{
								thick[ik].z_back_prop = 0;
							}
						}
					}

					if (!in_multem->is_multislice())
					{
						for(auto ithick = 0; ithick<thick.size(); ithick++)
						{
							thick[ithick].islice = ithick;
							thick[ithick].z_back_prop = thick[ithick].z_zero_def_plane - thick[ithick].z;
						}
					}

					return thick;
				}

				// slicing
				Vctr<Slice<T>, edev_cpu> get_slicing(In_Multem<T> *in_multem, 
				TVctr_r &z_slice, Vctr<Thick<T>, edev_cpu>& thick, Ptc_Atom<T>& patoms)
				{	
					if (!in_multem->is_multislice())
					{
						Vctr<Slice<T>, edev_cpu> slice(thick.size());
						for(auto islice = 0; islice<slice.size(); islice++)
						{
							slice[islice].z_0 = (islice == 0)?patoms.z_int_min:thick[islice-1].z;
							slice[islice].z_e = (thick.size() == 1)?patoms.z_int_max:thick[islice].z;
							slice[islice].z_int_0 = slice[islice].z_0;
							slice[islice].z_int_e = slice[islice].z_e;
							slice[islice].iatom_0 = (islice == 0)?0:(thick[islice-1].iatom_e+1);
							slice[islice].iatom_e = thick[islice].iatom_e;
							slice[islice] .ithk = islice;
						}
						return slice;
					}

					Vctr<Slice<T>, edev_cpu> slice(z_slice.size()-1);
					for(auto islice = 0; islice<slice.size(); islice++)
					{
						dt_bool Inc_Borders = false;
						slice[islice].z_0 = z_slice[islice];
						slice[islice].z_e = z_slice[islice + 1];
						switch(in_multem->pot_slic_typ)
						{
							case esst_planes_proj:
							{
								slice[islice].z_int_0 = slice[islice].z_0;
								slice[islice].z_int_e = slice[islice].z_e;
								Inc_Borders = false;
							}
							break;
							case esst_dz_proj:
							{
								slice[islice].z_int_0 = slice[islice].z_0;
								slice[islice].z_int_e = slice[islice].z_e;
								Inc_Borders = false;
							}
							break;
							case esst_dz_sub:
							{
								T z_m = slice[islice].z_m();

								slice[islice].z_int_0 = ::fmin(z_m - patoms.R_int_max, slice[islice].z_0);
								slice[islice].z_int_e = ::fmax(z_m + patoms.R_int_max, slice[islice].z_e);
								Inc_Borders = true;
							}
							break;
						}

						fd_by_z(patoms.z, slice[islice].z_int_0, slice[islice].z_int_e, slice[islice].iatom_0, slice[islice].iatom_e, Inc_Borders);
					
						slice[islice].ithk = -1;
					}

					// set thickness index
					for(auto ithk = 0; ithk<thick.size(); ithk++)
					{
						slice[thick[ithk].islice].ithk = ithk;
					}

					return slice;
				}

				// find patoms in slice
				dt_int32 fd_by_z(TVctr_r &z, T z_e, dt_bool Inc_Borders)
				{
					dt_int32 iz_e =-1;
					z_e = (Inc_Borders)?z_e+Epsilon<T>::rel:z_e;

					if (z_e<z.front())
					{
						return iz_e;
					}

					iz_e = (z.back() <= z_e)?(z.size()-1):(std::lower_bound(z.begin(), z.end(), z_e)-z.begin()-1);

					if (z[iz_e] < z.front())
					{
						iz_e = -1;
					}

					return iz_e;
				}

				// find patoms in slice
				void fd_by_z(TVctr_r &z, T z_0, T z_e, dt_int32& iz_0, dt_int32& iz_e, dt_bool Inc_Borders)
				{
					z_0 = (Inc_Borders)?(z_0-Epsilon<T>::rel):z_0;
					z_e = (Inc_Borders)?(z_e+Epsilon<T>::rel):z_e;

					if ((z_0>z_e)||(z.back()<z_0)||(z_e<z.front()))
					{
						iz_0 = 1;
						iz_e = 0;
						return;
					}

					iz_0 = (z_0 <= z.front())?0:(std::lower_bound(z.begin(), z.end(), z_0)-z.begin());
					iz_e = (z.back() <= z_e)?(z.size()-1):(std::lower_bound(z.begin(), z.end(), z_e)-z.begin()-1);

					if ((iz_0 > iz_e)||(z[iz_e] < z_0)||(z_e < z[iz_0]))
					{
						iz_0 = 1;
						iz_e = 0;
					}
				}
		};

	}

#endif
