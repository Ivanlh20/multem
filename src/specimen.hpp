/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
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
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SPECIMEN_H
#define SPECIMEN_H

#include <random>
#include <algorithm>

#include "math.cuh"
#include "types.hpp"
#include "atomic_data.hpp"
#include "atom_data.hpp"
#include "input_multislice.hpp"
#include "host_device_functions.cuh"
#include "host_functions.hpp"
#include "device_functions.cuh"

namespace multem
{
	template<class T, eDevice dev>
	class Specimen{
		public:
			using value_type_r = T;
			using size_type = std::size_t;

			static const eDevice device = dev;

			Specimen():input_multislice(nullptr){}

			void set_input_data(Input_Multislice<T, dev> *input_multislice_i)
			{
				input_multislice = input_multislice_i;

				/***************************************************************************/
				Atomic_Data atomic_data;
				atomic_data.Load_Data(input_multislice->potential_type);

				atom_type.resize(c_nAtomsTypes); 
				for(auto i=0; i<atom_type.size(); i++)
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

				/***************************************************************************/
				if(atoms_u.s_z_Int < 2.0*input_multislice->grid.dz)
				{
					input_multislice->grid.dz = atoms_u.s_z_Int;
					input_multislice->interaction_model = eESIM_Phase_Object;
				}

				atoms.set_Atoms(atoms_u, false , &atom_type);

				/***************************************************************************/
				get_slicing(slice);

				/***************************************************************************/
				rand.set_input_data(input_multislice->fp_seed, input_multislice->fp_dim);
			}

			/* Move atoms (ramdom distribution will be included in the future) */
			void move_atoms(const int &iconf)
			{
				if(!input_multislice->is_frozen_phonon()||(iconf <= 0))
				{
					return;
				}

				// set configuration
				rand.set_configuration(iconf);

				// Move atoms
				for(int iAtoms = 0; iAtoms<atoms_u.size(); iAtoms++)
				{
					atoms.Z[iAtoms] = atoms_u.Z[iAtoms];
					atoms.x[iAtoms] = atoms_u.x[iAtoms] + rand.dx(atoms_u.sigma[iAtoms]);
					atoms.y[iAtoms] = atoms_u.y[iAtoms] + rand.dy(atoms_u.sigma[iAtoms]);
					atoms.z[iAtoms] = atoms_u.z[iAtoms] + rand.dz(atoms_u.sigma[iAtoms]);
					atoms.sigma[iAtoms] = atoms_u.sigma[iAtoms];
					atoms.occ[iAtoms] = atoms_u.occ[iAtoms];
				}
				// get atom information
				atoms.get_Statistic(&atom_type);
				if(input_multislice->is_multislice())
				{
					// Ascending sort by z
					atoms.Sort_by_z();
					// Slicing procedure
					get_slicing(slice);
				}
			}

			T dz(const int &islice)
			{
				return slice.dz(islice)/cos(input_multislice->theta);
			}

			T dz_m(const int &islice_0, const int &islice_e)
			{
				return slice.dz_m(islice_0, islice_e)/cos(input_multislice->theta);
			}

			Vector<Atom_Type<T, e_Host>, e_Host>* ptr_atom_type()
			{
				return &atom_type;
			}

			//int IsInThickness(int islice);

			Input_Multislice<T, dev> *input_multislice; 			

			Vector<Atom_Type<T, e_Host>, e_Host> atom_type;		// Atom types
			Atom_Data<T> atoms; 								// displaced atoms
			Slice<T, e_Host> slice; 							// Slicing procedure
			Thickness<T, e_Host> thickness; 					// Thicknesses

		private:
			// get thickness
			void get_thickness(const Vector<T, e_Host> &z_slice, Thickness<T, e_Host> &thickness)
			{
				auto get_zero_defocus_plane = [&](const T &z_min, const T &z_max)->T
				{
					T zero_defocus_plane;
					switch(input_multislice->zero_defocus_type)
					{
						case eZDT_First:
						{
							zero_defocus_plane = z_min;
						}
						break;
						case eZDT_Middle:
						{
							zero_defocus_plane = 0.5*(z_min + z_max);
						}
						break;
						case eZDT_Last:
						{
							zero_defocus_plane = z_max;
						}
						break;
						default:
						{
							zero_defocus_plane = input_multislice->zero_defocus_plane;
						}
					}
					return zero_defocus_plane;
				};

				for(auto i=0; i<thickness.size(); i++)
				{
					int iatom_e, islice = 0;
					for(auto j=0; j<z_slice.size()-1; j++)
					{
						if((z_slice[j] < thickness.z[i])&&(thickness.z[i] <= z_slice[j+1]))
						{
							islice = j;
							break;
						}
					}

					thickness.islice[i] = islice;
					thickness.iatom_e[i] = iatom_e = atoms_u.find_by_z(z_slice[islice+1], false);
					thickness.z_zero_def_plane[i] = get_zero_defocus_plane(atoms_u.z[0], atoms_u.z[iatom_e]);
					if(input_multislice->is_through_slices())
					{
						thickness.z_back_prop[i] = 0;
					}
					else
					{
						thickness.z_back_prop[i] = thickness.z_zero_def_plane[i] - z_slice[islice+1];
						if(abs(thickness.z_back_prop[i])<Epsilon<float>::rel)
						{
							thickness.z_back_prop[i] = 0;
						}
					}
				}
			}

			// Slicing
			void get_slicing(Slice<T, e_Host> &slice)
			{				
				
				atoms_u.get_z_slice(input_multislice->potential_slicing, input_multislice->grid.dz, atoms, z_slice);
				get_thickness(z_slice, thickness);

				if(input_multislice->interaction_model != eESIM_Multislice){
					slice.resize(thickness.size());
					std::fill(slice.ithk.begin(), slice.ithk.end(), -1);

					for(auto islice=0; islice<slice.size(); islice++)
					{
						thickness.islice[islice] = islice;
						thickness.z_back_prop[islice] = thickness.z_zero_def_plane[islice] - thickness.z[islice];
						slice.z_0[islice] = (islice==0)?atoms.z_Int_min:thickness.z[islice-1]; 
						slice.z_e[islice] = thickness.z[islice]; 
						slice.z_int_0[islice] = slice.z_0[islice];
						slice.z_int_e[islice] = slice.z_e[islice];
						slice.iatom_0[islice] = (islice==0)?0:thickness.iatom_e[islice-1]+1;
						slice.iatom_e[islice] = thickness.iatom_e[islice];
						slice.ithk[islice] = islice;
					}
					return;
				}

				slice.resize(z_slice.size()-1);
				std::fill(slice.ithk.begin(), slice.ithk.end(), -1);
				for(auto ithk=0; ithk<thickness.size(); ithk++)
				{
					slice.ithk[thickness.islice[ithk]] = ithk;
				}

				for(auto islice=0; islice<slice.size(); islice++)
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
							slice.z_int_0[islice] = slice.z_0[islice] ;
							slice.z_int_e[islice] = slice.z_e[islice] ;

							atoms.find_by_z(slice.z_int_0[islice], slice.z_int_e[islice], slice.iatom_0[islice], slice.iatom_e[islice], false);
						}
						break;
						case ePS_dz_Sub:
						{
							T zm = 0.5*(slice.z_0[islice]+ slice.z_e[islice]);

							slice.z_int_0[islice] = min(zm - atoms.R_Int_max, slice.z_0[islice]);
							slice.z_int_e[islice] = max(zm + atoms.R_Int_max, slice.z_e[islice]);

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

					void set_configuration(int iconf)
					{	
						gen_u.seed(fp_seed);
						rand_u.reset();
						unsigned int seed_x, seed_y, seed_z;
						for(auto i=0; i<iconf; i++)
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
					//std::uniform_int_distribution<int> rand_u(0);
					//std::normal_distribution<T> randn_x(0.0, 1.0);
					//std::normal_distribution<T> randn_y(0.0, 1.0);
					//std::normal_distribution<T> randn_z(0.0, 1.0);
			};

			Random rand;
			// Undisplaced atoms
			Atom_Data<T> atoms_u;
			Vector<T, e_Host> z_slice;
	};

} // namespace multem

#endif