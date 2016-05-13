/*
 * This file is part of MULTEM.
 * Copyright 2016 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef SPECIMEN_H
#define SPECIMEN_H

#include <random>
#include <algorithm>

#include "math.cuh"
#include "types.cuh"
#include "lin_alg_def.cuh"
#include "atomic_data.hpp"
#include "atom_data.hpp"
#include "input_multislice.cuh"
#include "host_device_functions.cuh"
#include "host_functions.hpp"

namespace mt
{
	template<class T>
	class Specimen
	{
		public:
			using value_type_r = T;
			using size_type = std::size_t;

			Specimen(): input_multislice(nullptr){}

			void set_input_data(Input_Multislice<T> *input_multislice_i)
			{
				input_multislice = input_multislice_i;

				/***************************************************************************/
				Atomic_Data atomic_data(input_multislice->potential_type);

				atom_type.resize(c_nAtomsTypes); 
				for(auto i = 0; i<atom_type.size(); i++)
				{
					atomic_data.To_atom_type_CPU(i+1, input_multislice->Vrl, input_multislice->nR, input_multislice->grid.dR_min(), atom_type[i]);
				}

				/***************************************************************************/
				thickness.resize(input_multislice->thickness.size());
				thickness.z.assign(input_multislice->thickness.begin(), input_multislice->thickness.end());

				/***************************************************************************/
				atoms_u.set_Atoms(input_multislice->atoms, input_multislice->grid.pbc_xy , &atom_type);
				atoms_u.Sort_by_z();
				atoms_u.get_z_layer();

				//convergent beam
				input_multislice->cond_lens.zero_defocus_plane = input_multislice->cond_lens.get_zero_defocus_plane(atoms_u.z.front(), atoms_u.z.back());
				/***************************************************************************/
				if((atoms_u.s_z_Int < 2.0*input_multislice->grid.dz) || ((atoms_u.z_layer.size() == 1) && input_multislice->is_slicing_by_planes()))
				{
					input_multislice->grid.dz = atoms_u.s_z_Int;
					input_multislice->interaction_model = eESIM_Phase_Object;
					input_multislice->islice = 0;
					input_multislice->fp_dim.z = false;
					if(input_multislice->is_through_slices())
					{
						input_multislice->thickness_type = eTT_Through_Thickness;
					}
					input_multislice->slice_storage = input_multislice->slice_storage || !input_multislice->is_whole_specimen();
				}

				atoms.set_Atoms(atoms_u, false, &atom_type);
				// This is needed for memory preallocation in Transmission function
				get_slicing(slice);
				/***************************************************************************/
				rand.set_input_data(input_multislice->fp_seed, input_multislice->fp_dim);
			}

			/* Move atoms (random distribution will be included in the future) */
			void move_atoms(const int &fp_iconf)
			{
				// set configuration
				if(input_multislice->is_frozen_phonon())
				{
					rand.set_configuration(fp_iconf);
				}

				Vector<T, e_host> Rm;
				if(input_multislice->is_tomography())
				{
					Rm = get_rotation_matrix(input_multislice->tm_theta, input_multislice->tm_u0);
				}
				// Move atoms
				for(int iatoms = 0; iatoms<atoms_u.size(); iatoms++)
				{
					atoms.Z[iatoms] = atoms_u.Z[iatoms];
					r3d<T> r = atoms_u.to_r3d(iatoms);

					if(input_multislice->is_tomography())
					{
						r = r.rotate(Rm, input_multislice->tm_p0);
					}

					if(input_multislice->is_frozen_phonon())
					{
						r += rand.dr(atoms_u.sigma[iatoms], atoms_u.sigma[iatoms], atoms_u.sigma[iatoms]);
					}
					atoms.x[iatoms] = r.x;
					atoms.y[iatoms] = r.y;
					atoms.z[iatoms] = r.z;
					atoms.sigma[iatoms] = atoms_u.sigma[iatoms];
					atoms.occ[iatoms] = atoms_u.occ[iatoms];
					atoms.charge[iatoms] = atoms_u.charge[iatoms];
				}
				// get atom information
				atoms.get_Statistic(&atom_type);
				// Ascending sort by z
				atoms.Sort_by_z();
				// Slicing procedure
				get_slicing(slice);
			}

			T dz(const int &islice)
			{
				return slice.dz(islice)/cos(input_multislice->theta);
			}

			T dz_m(const int &islice_0, const int &islice_e)
			{
				return slice.dz_m(islice_0, islice_e)/cos(input_multislice->theta);
			}

			Vector<Atom_Type<T, e_host>, e_host>* ptr_atom_type()
			{
				return &atom_type;
			}

			// int IsInThickness(int islice);

			Input_Multislice<T> *input_multislice; 			

			Vector<Atom_Type<T, e_host>, e_host> atom_type;		// Atom types
			Atom_Data<T> atoms; 								// displaced atoms
			Slice<T, e_host> slice; 							// Slicing procedure
			Thickness<T, e_host> thickness; 					// Thicknesses

		private:
			// get thickness
			void get_thickness(const Vector<T, e_host> &z_slice, Thickness<T, e_host> &thickness)
			{
				auto get_islice = [&](const Vector<T, e_host> &z_slice, const T &z)->int
				{
					if(input_multislice->is_through_slices())
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

				auto b_sws = input_multislice->is_subslicing_whole_specimen();

				for(auto i = 0; i<thickness.size(); i++)
				{
					auto islice = (b_sws)?(z_slice.size()-2):get_islice(z_slice, thickness.z[i]);
					thickness.islice[i] = islice;

					auto iatom_e = atoms_u.find_by_z(z_slice[islice+1], false);
					thickness.iatom_e[i] = iatom_e;

					thickness.z_zero_def_plane[i] = input_multislice->obj_lens.get_zero_defocus_plane(atoms_u.z[0], atoms_u.z[iatom_e]);
					if(input_multislice->is_through_slices())
					{
						thickness.z_back_prop[i] = 0;
					}
					else
					{
						thickness.z_back_prop[i] = thickness.z_zero_def_plane[i] - z_slice[islice+1];
						if(fabs(thickness.z_back_prop[i])<Epsilon<float>::rel)
						{
							thickness.z_back_prop[i] = 0;
						}
					}
				}
			}

			// Slicing
			void get_slicing(Slice<T, e_host> &slice)
			{				
				atoms_u.get_z_slice(input_multislice->potential_slicing, input_multislice->grid.dz, atoms, z_slice);
				get_thickness(z_slice, thickness);

				if(input_multislice->interaction_model != eESIM_Multislice){
					slice.resize(thickness.size());
					std::fill(slice.ithk.begin(), slice.ithk.end(), -1);

					for(auto islice = 0; islice<slice.size(); islice++)
					{
						thickness.islice[islice] = islice;
						thickness.z_back_prop[islice] = thickness.z_zero_def_plane[islice] - thickness.z[islice];
						slice.z_0[islice] = (islice == 0)?atoms.z_Int_min:thickness.z[islice-1]; 
						slice.z_e[islice] = (thickness.size() == 1)?atoms.z_Int_max:thickness.z[islice];
						slice.z_int_0[islice] = slice.z_0[islice];
						slice.z_int_e[islice] = slice.z_e[islice];
						slice.iatom_0[islice] = (islice == 0)?0:thickness.iatom_e[islice-1]+1;
						slice.iatom_e[islice] = thickness.iatom_e[islice];
						slice.ithk[islice] = islice;
					}
					return;
				}

				slice.resize(z_slice.size()-1);
				std::fill(slice.ithk.begin(), slice.ithk.end(), -1);
				for(auto ithk = 0; ithk<thickness.size(); ithk++)
				{
					slice.ithk[thickness.islice[ithk]] = ithk;
				}
				for(auto islice = 0; islice<slice.size(); islice++)
				{
					slice.z_0[islice] = z_slice[islice]; 
					slice.z_e[islice] = z_slice[islice + 1];
					switch(input_multislice->potential_slicing)
					{
						case ePS_Planes:
						{
							slice.z_int_0[islice] = slice.z_0[islice];
							slice.z_int_e[islice] = slice.z_e[islice];

							atoms.find_by_z(slice.z_int_0[islice], slice.z_int_e[islice], slice.iatom_0[islice], slice.iatom_e[islice], false);
						}
						break;
						case ePS_dz_Proj:
						{
							slice.z_int_0[islice] = slice.z_0[islice];
							slice.z_int_e[islice] = slice.z_e[islice];

							atoms.find_by_z(slice.z_int_0[islice], slice.z_int_e[islice], slice.iatom_0[islice], slice.iatom_e[islice], false);
						}
						break;
						case ePS_dz_Sub:
						{
							T zm = 0.5*(slice.z_0[islice] + slice.z_e[islice]);

							slice.z_int_0[islice] = ::fmin(zm - atoms.R_Int_max, slice.z_0[islice]);
							slice.z_int_e[islice] = ::fmax(zm + atoms.R_Int_max, slice.z_e[islice]);

							atoms.find_by_z(slice.z_int_0[islice], slice.z_int_e[islice], slice.iatom_0[islice], slice.iatom_e[islice]);
						}
						break;
					}
				}
			}

			struct Random
			{
				public:
					void set_input_data(int &fp_seed_i, FP_Dim &fp_dim_i)
					{
						fp_seed = fp_seed_i;
						fp_dim = fp_dim_i;
					}

					void set_configuration(int fp_iconf)
					{	
						gen_u.seed(fp_seed);
						rand_u.reset();
						unsigned int seed_x, seed_y, seed_z;
						for(auto i = 0; i<fp_iconf; i++)
						{
							seed_x = rand_u(gen_u);
							seed_y = rand_u(gen_u);
							seed_z = rand_u(gen_u);
						}
						gen_x.seed(seed_x);
						randn_x.reset();

						gen_y.seed(seed_y);
						randn_y.reset();

						gen_z.seed(seed_z);
						randn_z.reset();
					}

					T dx(const T &sigma)
					{
						return (fp_dim.x)?sigma*randn_x(gen_x):0.0;
					}

					T dy(const T &sigma)
					{
						return (fp_dim.y)?sigma*randn_y(gen_y):0.0;
					}

					T dz(const T &sigma)
					{
						return (fp_dim.z)?sigma*randn_z(gen_z):0.0;
					}

					r3d<T> dr(const T &sigma_x, const T &sigma_y, const T &sigma_z)
					{
						return r3d<T>(dx(sigma_x), dy(sigma_y), dz(sigma_z));
					}

				private:
					int fp_seed;
					FP_Dim fp_dim;
					std::mt19937_64 gen_u;
					std::mt19937_64 gen_x;
					std::mt19937_64 gen_y;
					std::mt19937_64 gen_z;
					std::uniform_int_distribution<int> rand_u;
					std::normal_distribution<T> randn_x;
					std::normal_distribution<T> randn_y;
					std::normal_distribution<T> randn_z;
					// std::uniform_int_distribution<int> rand_u(0);
					// std::normal_distribution<T> randn_x(0.0, 1.0);
					// std::normal_distribution<T> randn_y(0.0, 1.0);
					// std::normal_distribution<T> randn_z(0.0, 1.0);
			};

			Random rand;
			// Undisplaced atoms
			Atom_Data<T> atoms_u;
			Vector<T, e_host> z_slice;
	};

} // namespace mt

#endif