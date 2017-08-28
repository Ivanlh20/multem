/*
 * This file is part of MULTEM.
 * Copyright 2017 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef SLICING_H
#define SLICING_H

#include "math.cuh"
#include "types.cuh"
#include "atom_data.hpp"
#include "host_functions.hpp"

namespace mt
{
	template <class T>
	class Input_Multislice;

	template <class T>
	struct Slicing
	{
		public:
			using value_type = T;
			using size_type = std::size_t;
			using TVector_r = Vector<T, e_host>;
			using TVector_i = Vector<int, e_host>;

			Slicing(): m_input_multislice(nullptr), m_atoms_r(nullptr), m_atoms(nullptr), z_eps(1e-3){}

			void set_input_data(Input_Multislice<T> *input_multislice, Atom_Data<T> *atoms_r, Atom_Data<T> *atoms=nullptr)
			{
				m_input_multislice = input_multislice;

				m_atoms_r = atoms_r;
				m_atoms = (atoms==nullptr)?atoms_r:atoms;

				z_plane = get_z_plane(m_input_multislice->potential_slicing, *m_atoms_r);
			}

			void match_thickness(ePotential_Slicing pot_sli, Atom_Data<T> &atoms, 
			eThick_Type thick_type, TVector_r &thick)
			{
				if(thick_type == eTT_Whole_Spec)
				{
					thick.resize(1);
					thick[0] = atoms.z_max;
					return;
				}

				auto z_plane = get_z_plane(pot_sli, atoms);

				if(thick_type == eTT_Through_Slices)
				{
					z_plane = get_z_slice(pot_sli, z_plane, atoms);
				}

				match_vectors(z_plane.begin(), z_plane.end(), thick, 0.25*atoms.dz);
			}

			T dz(int islice_0, int islice_e)
			{
				return (islice_e<slice.size())?(slice[islice_e].z_e - slice[islice_0].z_0):0.0;
			}

			T dz(int islice)
			{
				return dz(islice, islice);
			}

			T z_m(int islice)
			{
				return (islice<slice.size())?slice[islice].z_m():0.0;
			}

			T dz_m(int islice_0, int islice_e)
			{
				return fabs(z_m(islice_e) - z_m(islice_0));
			}

			void calculate()
			{
				m_z_slice = get_z_slice(m_input_multislice->potential_slicing, z_plane, *m_atoms);

				thick = get_thick(m_input_multislice, m_z_slice, *m_atoms_r);

				slice = get_slicing(m_input_multislice, m_z_slice, thick, *m_atoms);
			}

			TVector_r z_plane;
			Vector<Slice<T>, e_host> slice;
			Vector<Thick<T>, e_host> thick;

		private:
			const T z_eps;

			Input_Multislice<T> *m_input_multislice; 	

			Atom_Data<T> *m_atoms_r;
			Atom_Data<T> *m_atoms;

			TVector_r m_z_slice;
			Identify_Planes<T> identify_planes;

			// get spacing
			T get_spacing(size_type ix, TVector_r &x, T dz=0.5)
			{
				if(x.size()==1)
				{
					return dz;
				}

				ix = (ix <= 0)?1:min(ix, x.size()-1);
				return (x.size()>1)?x[ix]-x[ix-1]:0.0;
			}

			// Identify planes: Require that the atoms to be sorted along z
			TVector_r get_z_plane(ePotential_Slicing pot_sli, Atom_Data<T> &atoms)
			{
				TVector_r z_plane;

				if(atoms.size() == 0)
				{
					return z_plane;
				}

				atoms.validate_amorphous_parameters();

				const int region_ct = 0;

				// select z values of the crystal region
				TVector_r z_ct;
				z_ct.reserve(atoms.size());

				for(auto iz=0; iz<atoms.size(); iz++)
				{
					if(atoms.region[iz]==region_ct)
					{
						z_ct.push_back(atoms.z[iz]);
					}
				}
				std::sort(z_ct.begin(), z_ct.end());

				// calculate z crystal limits
				T z_ct_min = z_ct.front();
				T z_ct_max = z_ct.back();

				if(pot_sli==ePS_Planes)
				{
					z_plane = identify_planes(z_ct);
				}
				else
				{
					z_plane = identify_planes(z_ct_min, z_ct_max, atoms.dz);
				}

				std::vector<double> zt2(z_plane.begin(), z_plane.end());
	auto bb_ali = !atoms.amorp_lay_info.empty();
	if(bb_ali && (fabs(atoms.z_min-z_plane.front())>z_eps))
				{
					T dz_b = get_spacing(1, z_plane);
					auto amorp = atoms.amorp_lay_info.front();
					auto z_plane_top = identify_planes(amorp.z_0, amorp.z_e-dz_b, amorp.dz);
					z_plane.insert(z_plane.begin(), z_plane_top.begin(), z_plane_top.end());
				}

	 if(bb_ali && (fabs(z_plane.back()-atoms.z_max)>z_eps))
				{
					T dz_b = get_spacing(z_plane.size()-1, z_plane);
					auto amorp = atoms.amorp_lay_info.back();
					auto z_plane_bottom = identify_planes(amorp.z_0+dz_b, amorp.z_e, amorp.dz);
					z_plane.insert(z_plane.end(), z_plane_bottom.begin(), z_plane_bottom.end());
				}

				// unique vector
				unique_vector(z_plane);

				return z_plane;
			}

			// get z slicing
			TVector_r get_z_slice(ePotential_Slicing pot_sli, TVector_r &z_plane, Atom_Data<T> &atoms)
			{
				TVector_r z_slice;

				if((pot_sli!=ePS_dz_Sub)&&(z_plane.size() == 1))
				{
					z_slice.resize(2);
					z_slice[0] = atoms.z_int_min;
					z_slice[1] = atoms.z_int_max;
					z_slice.shrink_to_fit();

					return z_slice;
				}

				TVector_r z_plane_sli = z_plane;

				z_slice.resize(z_plane_sli.size()+1);
				z_slice[0] = z_plane_sli[0]-0.5*get_spacing(0, z_plane_sli, atoms.dz);
				for(auto iz=1; iz<z_slice.size()-1; iz++)
				{
					z_slice[iz] = 0.5*(z_plane_sli[iz]+z_plane_sli[iz-1]);
				}
				z_slice[z_slice.size()-1] = z_slice[z_slice.size()-2]+get_spacing(z_plane_sli.size(), z_plane_sli, atoms.dz);

				if(pot_sli==ePS_dz_Sub)
				{
					T dz_b = get_spacing(1, z_plane_sli, atoms.dz);
					if(atoms.z_int_min<z_slice.front()-dz_b)
					{
						T dz_s = get_spacing(2, z_plane_sli, atoms.dz);
						auto z_slice_top = identify_planes(atoms.z_int_min, z_slice.front()-dz_b, dz_s);
						z_slice.insert(z_slice.begin(), z_slice_top.begin(), z_slice_top.end());
					}
					else
					{
						z_slice[0] = atoms.z_int_min;
					}

					dz_b = get_spacing(z_plane_sli.size()-1, z_plane_sli, atoms.dz);
					if(z_slice.back()+dz_b<atoms.z_int_max)
					{
						T dz_s = get_spacing(z_plane_sli.size()-2, z_plane_sli, atoms.dz);
						auto z_slice_bottom = identify_planes(z_slice.back()+dz_b, atoms.z_int_max, dz_s);
						z_slice.insert(z_slice.end(), z_slice_bottom.begin(), z_slice_bottom.end());
					}
					else
					{
						z_slice[z_slice.size()-1] = atoms.z_int_max;
					}
				}

				z_slice.shrink_to_fit();

				return z_slice;
			}

			// get thick
			Vector<Thick<T>, e_host> get_thick(Input_Multislice<T> *input_multislice, 
			TVector_r &z_slice, Atom_Data<T> &atoms)
			{
				const auto thick_type = input_multislice->thick_type;

				auto get_islice = [thick_type](TVector_r &z_slice, const T &z)->int
				{
					if(thick_type==eTT_Through_Slices)
					{
						for(auto i = 0; i<z_slice.size()-1; i++)
						{
							if(fabs(z_slice[i+1]-z)<Epsilon<float>::rel)
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
							if((z_slice[i] < z)&&(z <= z_slice[i+1]))
							{
								return i;
							}
						}
						return 0;
					}
				};

				auto b_sws = input_multislice->is_subslicing_whole_spec();

				Vector<Thick<T>, e_host> thick(input_multislice->thick.size());
				for(auto ik = 0; ik<thick.size(); ik++)
				{
					thick[ik].z = input_multislice->thick[ik];
					auto islice = (b_sws)?(z_slice.size()-2):get_islice(z_slice, thick[ik].z);
					thick[ik].islice = islice;

					auto iatom_e = fd_by_z(atoms.z, z_slice[islice+1], false);
					thick[ik].iatom_e = iatom_e;

					thick[ik].z_zero_def_plane = input_multislice->obj_lens.get_zero_defocus_plane(atoms.z[0], atoms.z[iatom_e]);
					if(input_multislice->is_through_slices())
					{
						thick[ik].z_back_prop = 0;
					}
					else
					{
						// I need to recheck this part, the average should be replace by the plane
						thick[ik].z_back_prop = thick[ik].z_zero_def_plane - 0.5*(z_slice[islice]+z_slice[islice+1]);
						if(fabs(thick[ik].z_back_prop)<Epsilon<float>::rel)
						{
							thick[ik].z_back_prop = 0;
						}
					}
				}

				if(!input_multislice->is_multislice())
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
			Vector<Slice<T>, e_host> get_slicing(Input_Multislice<T> *input_multislice, 
			TVector_r &z_slice, Vector<Thick<T>, e_host> &thick, Atom_Data<T> &atoms)
			{	
				if(!input_multislice->is_multislice())
				{
					Vector<Slice<T>, e_host> slice(thick.size());
					for(auto islice = 0; islice<slice.size(); islice++)
					{
						slice[islice].z_0 = (islice == 0)?atoms.z_int_min:thick[islice-1].z; 
						slice[islice].z_e = (thick.size() == 1)?atoms.z_int_max:thick[islice].z;
						slice[islice].z_int_0 = slice[islice].z_0;
						slice[islice].z_int_e = slice[islice].z_e;
						slice[islice].iatom_0 = (islice == 0)?0:(thick[islice-1].iatom_e+1);
						slice[islice].iatom_e = thick[islice].iatom_e;
						slice[islice] .ithk = islice;
					}
					return slice;
				}

				Vector<Slice<T>, e_host> slice(z_slice.size()-1);
				for(auto islice = 0; islice<slice.size(); islice++)
				{
					bool Inc_Borders = false;
					slice[islice].z_0 = z_slice[islice]; 
					slice[islice].z_e = z_slice[islice + 1];
					switch(input_multislice->potential_slicing)
					{
						case ePS_Planes:
						{
							slice[islice].z_int_0 = slice[islice].z_0;
							slice[islice].z_int_e = slice[islice].z_e;
							Inc_Borders = false;
						}
						break;
						case ePS_dz_Proj:
						{
							slice[islice].z_int_0 = slice[islice].z_0;
							slice[islice].z_int_e = slice[islice].z_e;
							Inc_Borders = false;
						}
						break;
						case ePS_dz_Sub:
						{
							T z_m = slice[islice].z_m();

							slice[islice].z_int_0 = ::fmin(z_m - atoms.R_int_max, slice[islice].z_0);
							slice[islice].z_int_e = ::fmax(z_m + atoms.R_int_max, slice[islice].z_e);
							Inc_Borders = true;
						}
						break;
					}

					fd_by_z(atoms.z, slice[islice].z_int_0, slice[islice].z_int_e, slice[islice].iatom_0, slice[islice].iatom_e, Inc_Borders);
					
					slice[islice].ithk = -1;
				}

				//set thickness index
				for(auto ithk = 0; ithk<thick.size(); ithk++)
				{
					slice[thick[ithk].islice].ithk = ithk;
				}

				return slice;
			}

			// find atoms in slice
			int fd_by_z(TVector_r &z, T z_e, bool Inc_Borders)
			{
				int iz_e =-1;
				z_e = (Inc_Borders)?z_e+Epsilon<T>::rel:z_e;

				if(z_e<z.front())
				{
					return iz_e;
				}

				iz_e = (z.back() <= z_e)?(z.size()-1):(std::lower_bound(z.begin(), z.end(), z_e)-z.begin()-1);

				if(z[iz_e] < z.front())
				{
					iz_e = -1;
				}

				return iz_e;
			}

			// find atoms in slice
			void fd_by_z(TVector_r &z, T z_0, T z_e, int &iz_0, int &iz_e, bool Inc_Borders)
			{
				z_0 = (Inc_Borders)?(z_0-Epsilon<T>::rel):z_0;
				z_e = (Inc_Borders)?(z_e+Epsilon<T>::rel):z_e;

				if((z_0>z_e)||(z.back()<z_0)||(z_e<z.front()))
				{
					iz_0 = 1;
					iz_e = 0;
					return;
				}

				iz_0 = (z_0 <= z.front())?0:(std::lower_bound(z.begin(), z.end(), z_0)-z.begin());
				iz_e = (z.back() <= z_e)?(z.size()-1):(std::lower_bound(z.begin(), z.end(), z_e)-z.begin()-1);

				if((iz_0 > iz_e)||(z[iz_e] < z_0)||(z_e < z[iz_0]))
				{
					iz_0 = 1; 
					iz_e = 0;
				}
			}
	};

} // namespace mt

#endif
