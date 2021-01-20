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

#ifndef SPECIMEN_H
#define SPECIMEN_H

#include "math.cuh"
#include "types.cuh"
#include "lin_alg_def.cuh"
#include "cgpu_rand.cuh"
#include "atomic_data.hpp"
#include "atomic_data_mt.hpp"
#include "input_multislice.cuh"
#include "cgpu_fcns.cuh"
#include "cpu_fcns.hpp"
#include "slicing.hpp"

namespace mt
{
	template <class T>
	class Spec
	{
		public:
			using T_r = T;
			using size_type = std::size_t;

			Spec(): input_multislice(nullptr){}

			void set_input_data(Input_Multislice<T> *input_multislice_i)
			{
				input_multislice = input_multislice_i;

				/***************************************************************************/
				Atomic_Data atomic_data(input_multislice->potential_type);

				atom_type.resize(c_nAtomsTypes); 
				for(auto i = 0; i<atom_type.size(); i++)
				{
					atomic_data.To_atom_type_CPU(i+1, input_multislice->Vrl, input_multislice->nR, input_multislice->grid_2d.dR_min(), atom_type[i]);
				}

				/***************************************************************************/
				atoms_u.set_atoms(input_multislice->atoms, input_multislice->grid_2d.pbc_xy, &atom_type);
				atoms_u.sort_by_z();

				/***************************************************************************/
				slicing.set_input_data(input_multislice, &atoms_u, &atoms);

				/***************************************************************************/
				if((input_multislice->is_phase_object()) || (atoms_u.s_z_int < 2.0*input_multislice->grid_2d.dz) || ((slicing.z_plane.size() == 1) && input_multislice->is_slicing_by_planes()))
				{
					input_multislice->grid_2d.dz = atoms_u.s_z_int;
					input_multislice->interaction_model = eESIM_Phase_Object;
					input_multislice->islice = 0;
					input_multislice->pn_dim.z = false;
					if(input_multislice->is_through_slices())
					{
						input_multislice->thick_type = eTT_Through_Thick;
					}
					input_multislice->slice_storage = input_multislice->slice_storage || !input_multislice->is_whole_spec();
					atoms_u.dz = input_multislice->grid_2d.dz;
				}

				atoms.set_atoms(atoms_u, false, &atom_type);
				// This is needed for memory preallocation in Transmission function
				slicing.calculate();
			}

			/* Move atoms (cgpu_rand distribution will be included in the future) */
			void move_atoms(const int &fp_iconf)
			{
				// set phonon configuration
				if(input_multislice->is_frozen_phonon())
				{
					rand.seed(input_multislice->pn_seed, fp_iconf);
					rand.set_activation(input_multislice->pn_dim.x, input_multislice->pn_dim.y, input_multislice->pn_dim.z);
				}

				// move atoms
				for(int iatoms = 0; iatoms<atoms_u.size(); iatoms++)
				{
					atoms.Z[iatoms] = atoms_u.Z[iatoms];
					auto r = atoms_u.to_r3d(iatoms);

					if(input_multislice->is_frozen_phonon())
					{
						auto sigma_x = atoms_u.sigma[iatoms];
						auto sigma_y = atoms_u.sigma[iatoms];
						auto sigma_z = atoms_u.sigma[iatoms];
						r += rand(sigma_x, sigma_y, sigma_z);
					}

					atoms.x[iatoms] = r.x;
					atoms.y[iatoms] = r.y;
					atoms.z[iatoms] = r.z;
					atoms.sigma[iatoms] = atoms_u.sigma[iatoms];
					atoms.occ[iatoms] = atoms_u.occ[iatoms];
					atoms.region[iatoms] = atoms_u.region[iatoms];
					atoms.charge[iatoms] = atoms_u.charge[iatoms];
				}

				if(input_multislice->pn_dim.z)
				{
					atoms.sort_by_z();
				}

				// get atom information
				atoms.get_statistic(&atom_type);

				// slicing procedure
				slicing.calculate();
			}

			T dz(const int &islice)
			{
				return slicing.dz(islice)/cos(input_multislice->theta);
			}

			T dz_m(const int &islice_0, const int &islice_e)
			{
				return slicing.dz_m(islice_0, islice_e)/cos(input_multislice->theta);
			}

			Input_Multislice<T> *input_multislice; 			

			Atom_Data<T> atoms; 								// displaced atoms
			Slicing<T> slicing; 								// slicing procedure
			Vector<Atom_Type<T, e_host>, e_host> atom_type;		// Atom types
		private:
			Randn_3d<T, e_host> rand;
			Atom_Data<T> atoms_u;
	};

} // namespace mt

#endif
