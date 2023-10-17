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

#ifndef SPEC_H
#define SPEC_H

#include "math.cuh"
#include "types.cuh"
#include "cgpu_rand.cuh"
#include "atomic_data_mt.cuh"
#include "particles.cuh"
#include "in_classes.cuh"
#include "cgpu_fcns.cuh"
#include "cpu_fcns.hpp"
#include "spec_slic.hpp"

namespace mt
{
	template <class T>
	class Spec
	{
		public:
			using T_r = T;
			using size_type = dt_uint64;

			Spec(): multem_in_parm(nullptr) {}

			void set_in_data(Multem_In_Parm<T> *multem_in_parm_i)
			{
				multem_in_parm = multem_in_parm_i;

				 /***************************************************************************************/
				Atomic_Data atomic_data_mt(multem_in_parm->atomic_pot_parm_typ);

				atom_type.resize(c_n_atom_typ);
				for(auto i = 0; i < atom_type.size(); i++)
				{
					atomic_data_mt.To_atom_type_CPU(i+1, multem_in_parm->Vrl, multem_in_parm->nR, multem_in_parm->grid_2d.dR_min(), atom_type[i]);
				}

				 /***************************************************************************************/
				atoms_u.set_ptc(multem_in_parm->atoms, multem_in_parm->grid_2d.pbc_xy, &atom_type);
				atoms_u.sort_by_z();

				 /***************************************************************************************/
				slicing.set_in_data(multem_in_parm, &atoms_u, &atoms);

				 /***************************************************************************************/
				if ((atoms_u.s_z_int < 2.0*multem_in_parm->grid_2d.sli_thick) || ((slicing.z_plane.size() == 1) && multem_in_parm->is_spec_slic_by_plns_proj()))
				{
					multem_in_parm->grid_2d.sli_thick = atoms_u.s_z_int;
					multem_in_parm->elec_spec_interact_mod = eesim_phase_object;
					multem_in_parm->islice = 0;
					multem_in_parm->atomic_vib.dim_z = false;
					if (multem_in_parm->is_sim_through_slices())
					{
						multem_in_parm->thick_type = estt_through_thick;
					}
					multem_in_parm->slice_storage = multem_in_parm->slice_storage || !multem_in_parm->is_sim_whole_spec();
					atoms_u.sli_thick = multem_in_parm->grid_2d.sli_thick;
				}

				atoms.set_ptc(atoms_u, false, &atom_type);
				// This is needed for memory preallocation in Transmission function
				slicing.calculate();
			}

			/* Move atoms (random distribution will be included in the future) */
			void move_atoms(const dt_int32& fp_iconf)
			{
				// set atomic_vib configuration
				if (multem_in_parm->is_avm_frozen_phonon())
				{
					rand.seed(multem_in_parm->atomic_vib.seed, fp_iconf);
					rand.set_act_dim(multem_in_parm->atomic_vib.dim_x, multem_in_parm->atomic_vib.dim_y, multem_in_parm->atomic_vib.dim_z);
				}

				// move atoms
				for(dt_int32 iatoms = 0; iatoms<atoms_u.size(); iatoms++)
				{
					atoms.Z[iatoms] = atoms_u.Z[iatoms];
					auto r = atoms_u.get_pos(iatoms);

					if (multem_in_parm->is_avm_frozen_phonon())
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
					atoms.tag[iatoms] = atoms_u.tag[iatoms];
					atoms.charge[iatoms] = atoms_u.charge[iatoms];
				}

				if (multem_in_parm->atomic_vib.dim_z)
				{
					atoms.sort_by_z();
				}

				// get atom information
				atoms.get_statistic(&atom_type);

				// slicing procedure
				slicing.calculate();
			}

			T sli_thick(const dt_int32& islice)
			{
				return slicing.sli_thick(islice)/cos(multem_in_parm->theta);
			}

			T dz_m(const dt_int32& islice_0, const dt_int32& islice_e)
			{
				return slicing.dz_m(islice_0, islice_e)/cos(multem_in_parm->theta);
			}

			Multem_In_Parm<T> *multem_in_parm;

			Ptc_Atom<T> atoms;		// displaced atoms
			Spec_Slic<T> slicing;		// slicing procedure
			Vctr<Atom_Typ_cpu<T>, edev_cpu> atom_type;		// Atom types
		private:
			Randn_3d<T, edev_cpu> rand;
			Ptc_Atom<T> atoms_u;
	};

}

#endif