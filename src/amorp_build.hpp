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

#ifndef AMORPHOUS_BUILD_H
	#define AMORPHOUS_BUILD_H

	#include "math.cuh"
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
					Ptc_Atom<T> atoms;
					this->operator()(atoms, Z, rms_3d, occ, d_min, rho, seed, spec_lay_info);

					return atoms;
				}

				void operator()(Ptc_Atom<T>& atoms, dt_int32 Z, T rms_3d, T occ, T d_min, T rho, dt_int32 seed, Spec_Lay_Info<T> spec_lay_info)
				{
					// calculate density = #atoms/A^3
					const auto rho_n_A3 = density_n_A3(Z, rho);

					// calculate number of atoms for the amorphous spec
					const auto n_atoms_amorp = n_amorp_atoms(rho_n_A3, spec_lay_info.bs);

					m_atoms.bs = atoms.bs;
					m_atoms.cols_used = atoms.cols_used;

					box_occ.set_in_data(d_min, spec_lay_info.bs, spec_lay_info.r_0);

					// set atoms
					dt_int32 iatom_0 = 0;
					dt_int32 n_atoms = iatom_0 + n_atoms_amorp;

					if (atoms.empty())
					{
						m_atoms.reserve(n_atoms);
					}
					else
					{
						const dt_int32 n_atoms_0 = atoms.size();
						m_atoms.reserve(n_atoms_0);

						T depth = 0.1; // set depth
						// select atoms within the depth
						const T z_0 = spec_lay_info.r_0.z - depth;
						const T z_e = spec_lay_info.r_0.z + spec_lay_info.bs.z + depth;

						for(auto iatom = 0; iatom < n_atoms_0; iatom++)
						{
							if (fcn_chk_bound(atoms.z[iatom], z_0, z_e))
							{
								m_atoms.push_back(atoms.get(iatom));
								iatom_0++;
							}
						}
						n_atoms = iatom_0 + n_atoms_amorp;
						m_atoms.reserve(n_atoms);
	
						// set occupancy
						for (auto iatom = 0; iatom < iatom_0; iatom++)
						{
							box_occ.set(m_atoms.get_pos(iatom), iatom);
						}
					}

					// set random input data
					randu_3d.set_in_data(seed, spec_lay_info.bs, spec_lay_info.r_0); 

					const dt_int32 region = spec_lay_info.region;

					atoms.reserve(atoms.size() + n_atoms_amorp);

					dt_int32 iatom_c = iatom_0;
					for(dt_int32 iatom = iatom_0; iatom < n_atoms; iatom++)
					{
						R_3d<T> r;
						if (rnd_point(m_atoms, r))
						{
							const Ptc_s_Atom<T> atom_ik(Z, r, rms_3d, occ, region, c_dflt_charge);

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

				Ptc_Atom<T> m_atoms;

			private:
				const dt_int32 n_trial;

				Randu_3d_cpu<T> randu_3d;

				Box_Occ_3d<dt_float64> box_occ;

				dt_int32 n_amorp_atoms(const T& rho_n_A3, const R_3d<T>& bs) const 
				{ 
					return fcn_cceil<dt_int32>(rho_n_A3*vol(bs));	
				}					
				
				T vol(const R_3d<T>& bs) const
				{ 
					return bs.x*bs.y*bs.z;
				}

				T density_n_A3(dt_int32 Z, T m_rho)
				{
					Atomic_Info_cpu<T> atomic_info = Atomic_Data(Z);
					const T m = atomic_info.atomic_mass();

					return m_rho*c_Na<T>/(m*c_cm3_A3<T>);
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