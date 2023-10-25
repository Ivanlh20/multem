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

#ifndef AMORPHOUS_BUILD_H
	#define AMORPHOUS_BUILD_H

	#include "math_mt.h"
	#include "types.cuh"
	#include "types_mt.cuh"
	#include "particles.cuh"
	#include "atomic_data_mt.cuh"
	#include "cgpu_rand.cuh"
	#include "box_occ.hpp"

	namespace mt
	{
		template <class T>
		class Amorp_Build {
			public:
				Amorp_Build():n_trial(500), randu_3d() {}

				Ptc_Atom<T> operator()(dt_int32 Z, T rms_3d, T occ, T d_min, T rho, dt_int32 seed, Spec_Lay_Info<T> spec_lay_info)
				{
					// calculate number of atoms for the amorphous spec
					const auto n_atoms_amorp = n_amorp_atoms(Z, rho, spec_lay_info.bs);

					Ptc_Atom<T> atoms;
					atoms.bs = spec_lay_info.bs;
					atoms.cols_used = 7;
					atoms.reserve(n_atoms_amorp);

					// set box occupancy input data
					box_occ.set_in_data(d_min, spec_lay_info.bs, spec_lay_info.r_0);

					// set random input data
					randu_3d.set_in_data(seed, spec_lay_info.bs, spec_lay_info.r_0); 

					const dt_int32 tag = spec_lay_info.tag;

					dt_int32 iatom_c = 0;
					for(dt_int32 iatom = 0; iatom < n_atoms_amorp; iatom++)
					{
						R_3d<T> r;
						if (rnd_point(atoms, r))
						{
							const Ptc_s_Atom<T> atom_ik(Z, r, rms_3d, occ, tag, c_dflt_charge);

							atoms.push_back(atom_ik);

							box_occ.set(r, iatom_c);
							iatom_c++;
						}
					}

					atoms.shrink_to_fit();
					atoms.sort_by_z();

					return atoms;
				}

				void operator()(Ptc_Atom<T>& atoms, dt_int32 Z, T rms_3d, T occ, T d_min, T rho, dt_int32 seed, Spec_Lay_Info<T> spec_lay_info)
				{
					// calculate number of atoms for the amorphous spec
					const auto n_atoms_amorp = n_amorp_atoms(Z, rho, spec_lay_info.bs);

					Ptc_Atom<T> m_atoms;
					m_atoms.bs = fmax(atoms.bs, spec_lay_info.bs);
					m_atoms.cols_used = atoms.cols_used;

					// set box occupancy input data
					box_occ.set_in_data(d_min, spec_lay_info.bs, spec_lay_info.r_0);

					// select atoms within the depth
					const T depth = 0.1; // set depth
					const T z_0 = spec_lay_info.r_0.z - depth;
					const T z_e = spec_lay_info.z_e() + depth;

					m_atoms.reserve(atoms.size());

					dt_int32 iatom_0 = 0;
					for(auto iatom = 0; iatom < atoms.size(); iatom++)
					{
						if (fcn_chk_bound(atoms.z[iatom], z_0, z_e))
						{
							const auto atom_ik = atoms.get(iatom);
							m_atoms.push_back(atom_ik);
							// set occupancy
							box_occ.set(atom_ik.get_pos(), iatom_0);
							iatom_0++;
						}
					}
					const dt_int32 n_atoms = iatom_0 + n_atoms_amorp;
					m_atoms.reserve(n_atoms);
					atoms.reserve(atoms.size() + n_atoms_amorp);

					// set random input data
					randu_3d.set_in_data(seed, spec_lay_info.bs, spec_lay_info.r_0); 

					const dt_int32 tag = spec_lay_info.tag;

					dt_int32 iatom_c = iatom_0;
					for(dt_int32 iatom = iatom_0; iatom < n_atoms; iatom++)
					{
						R_3d<T> r;
						if (rnd_point(m_atoms, r))
						{
							const Ptc_s_Atom<T> atom_ik(Z, r, rms_3d, occ, tag, c_dflt_charge);

							m_atoms.push_back(atom_ik);
							atoms.push_back(atom_ik);

							box_occ.set(r, iatom_c);
							iatom_c++;
						}
					}
					m_atoms.clear();
					m_atoms.shrink_to_fit();

					atoms.shrink_to_fit();
					atoms.sort_by_z();
				}

			private:
				const dt_int32 n_trial;

				Rndu_3d_cpu<T> randu_3d;

				Box_Occ_3d<dt_float64> box_occ;
				
				T vol(const R_3d<T>& bs) const
				{ 
					return bs.x*bs.y*bs.z;
				}

				T rho_n_A3(dt_int32 Z, T rho) const
				{
					Atomic_Info_cpu<T> atomic_info = Atomic_Data(Z);
					const T m = atomic_info.atomic_mass();

					return rho*c_Na<T>/(m*c_cm3_A3<T>);
				}

				dt_int32 n_amorp_atoms(dt_int32 Z, const T& rho, const R_3d<T>& bs) const 
				{ 
					return fcn_cceil<dt_int32>(rho_n_A3(Z, rho)*vol(bs));	
				}					

				dt_bool rnd_point(const Ptc_Atom<T>& atoms, R_3d<T>& r_o)
				{
					for(auto itrial = 0; itrial < n_trial; itrial++)
					{
						const auto r = randu_3d();
						if (box_occ.check_r_min(atoms, r))
						{
							r_o = r;
							return true;
						}
					}

					return false;
				}
		};
	}
#endif