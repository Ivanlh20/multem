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
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef AMORPHOUS_SPECIMEN_H
#define AMORPHOUS_SPECIMEN_H

#include "math.cuh"
#include "types.cuh"
#include "lin_alg_def.cuh"
#include "atomic_data_mt.hpp"
#include "atomic_data.hpp"
#include "cgpu_rand.cuh"
#include "box_occ.hpp"

namespace mt
{
	template <class T>
	class Amorp_Spec {
	public:
		Amorp_Spec() : m_ntrial(500), m_depth(0.01), 
		m_l_x(0), m_l_y(0), m_l_z(0), m_al_type(eALT_Top) {}

		void create(Atom_Data<T> &atoms, T d_min, int Z, 
		T rms_3d, T rho, int seed = 300183)
		{
			T depth = m_depth;

			const int iatoms_0 = set_init_values(atoms, d_min, Z, rms_3d, rho, seed, depth);
			const int region = atoms.amorp_lay_info[0].region;

			int iatoms_c = iatoms_0;
			for (int iatoms = iatoms_0; iatoms < m_atoms.size(); iatoms++)
			{
				r3d<T> r;
				if (rand_point(m_atoms, r))
				{
					box.set_occ(r, iatoms_c);
					m_atoms.Z[iatoms_c] = Z;
					m_atoms.x[iatoms_c] = r.x;
					m_atoms.y[iatoms_c] = r.y;
					m_atoms.z[iatoms_c] = r.z;
					m_atoms.sigma[iatoms_c] = rms_3d;
					m_atoms.occ[iatoms_c] = 1.0;
					m_atoms.region[iatoms_c] = region;
					m_atoms.charge[iatoms_c] = 0;
					iatoms_c++;
				}
			}
			m_atoms.resize(iatoms_c);
			m_atoms.shrink_to_fit();

			depth = box.l_z-atoms.amorp_lay_info[0].lz();
			const int natoms_am = m_atoms.size() - iatoms_0;
			iatoms_c = atoms.size();
			atoms.resize(iatoms_c+natoms_am);
			for (int iatoms = iatoms_0; iatoms < m_atoms.size(); iatoms++)
			{
				atoms.Z[iatoms_c] = m_atoms.Z[iatoms];
				atoms.x[iatoms_c] = m_atoms.x[iatoms];
				atoms.y[iatoms_c] = m_atoms.y[iatoms];
				atoms.z[iatoms_c] = (m_al_type==eALT_Top)?(atoms.z_min+depth-m_atoms.z[iatoms]):(atoms.z_max-depth+m_atoms.z[iatoms]);
				atoms.sigma[iatoms_c] = m_atoms.sigma[iatoms];
				atoms.occ[iatoms_c] = m_atoms.occ[iatoms];
				atoms.region[iatoms_c] = m_atoms.region[iatoms];
				atoms.charge[iatoms_c] = m_atoms.charge[iatoms];
				iatoms_c++;
			}
			atoms.sort_by_z();
		}

		Atom_Data<T> m_atoms;

	private:
		const int m_ntrial;
		const T m_depth;

		T m_l_x;
		T m_l_y;
		T m_l_z;
		T m_rho;
		eAmorp_Lay_Type m_al_type;

		Box_Occ<T> box;
		Rand_3d<T, e_host> rand_3d;

		int set_init_values(Atom_Data<T> &atoms, T d_min, int Z, 
		T rms_3d, T rho, int seed, T depth)
		{
			m_l_x = atoms.l_x;
			m_l_y = atoms.l_y;
			m_l_z = atoms.amorp_lay_info[0].lz();
			m_al_type = atoms.amorp_lay_info[0].type;
			m_rho = convert_density(Z, rho);

			// set cgpu_rand set
			rand_3d.seed(seed);

			// calculate number of atoms for the amorphous spec
			int natoms_am = static_cast<int>(ceil(m_rho*vol()));
			m_atoms.l_x = m_l_x;
			m_atoms.l_y = m_l_y;
			m_atoms.l_z = m_l_z;

			if(atoms.empty())
			{
				m_atoms.resize(natoms_am);

				// set amorphous box size
				box.set_input_data(d_min, m_l_x, m_l_y, m_l_z);
				
				// set cgpu_rand box size
				rand_3d.set_box_size(m_l_x, m_l_y, m_l_z);

				return 0;
			}

			const int natoms_ct = atoms.size();
			m_atoms.resize(natoms_ct);

			// select atoms within the depth
			T z_e = (m_al_type==eALT_Top)?(atoms.z_min+depth):(atoms.z_max-depth);

			int iatoms_c = 0;
			for(auto iatoms = 0; iatoms < natoms_ct; iatoms++)
			{
				auto bb = (m_al_type==eALT_Top)?(atoms.z[iatoms]<z_e):(atoms.z[iatoms]>z_e);
				if(bb)
				{
					m_atoms.Z[iatoms_c] = atoms.Z[iatoms];
					m_atoms.x[iatoms_c] = atoms.x[iatoms];
					m_atoms.y[iatoms_c] = atoms.y[iatoms];
					m_atoms.z[iatoms_c] = atoms.z[iatoms];
					m_atoms.sigma[iatoms_c] = atoms.sigma[iatoms];
					m_atoms.occ[iatoms_c] = atoms.occ[iatoms];
					m_atoms.region[iatoms_c] = atoms.region[iatoms];
					m_atoms.charge[iatoms_c] = atoms.charge[iatoms];
					iatoms_c++;
				}
			}
			m_atoms.resize(iatoms_c);

			T z_min, z_max;
			minmax_element(m_atoms.z, z_min, z_max);

			auto depth_ct = (m_al_type==eALT_Top)?(z_max-atoms.z_min):(atoms.z_max-z_min);

			// set amorphous box size
			box.set_input_data(d_min, m_l_x, m_l_y, m_l_z+depth_ct);

			// set cgpu_rand box size
			rand_3d.set_box_size(m_l_x, m_l_y, m_l_z+depth_ct);

			// set number of atoms
			const int natoms_sct = iatoms_c;
			m_atoms.resize(natoms_sct+natoms_am);
			m_atoms.shrink_to_fit();

			// change reference system
			for(auto iatoms = 0; iatoms < natoms_sct; iatoms++)
			{
				m_atoms.z[iatoms] = (m_al_type==eALT_Top)?(z_max-m_atoms.z[iatoms]):(m_atoms.z[iatoms]-z_min);
				box.set_occ(m_atoms.to_r3d(iatoms), iatoms);
			}

			 return natoms_sct;
		}

		T convert_density(int Z, T m_rho)
		{
			const double Na = 6.022140857e23;
			const double cm3_A3 = 1e24;

			Atomic_Data atomic_data;
			const T m = atomic_data.atomic_mass(Z);

			return m_rho*Na/(m*cm3_A3);
		}
		
		T vol() const 
		{ 
			return m_l_x*m_l_y*m_l_z; 
		}

		T density(int natoms) const
		{
			return natoms/vol();
		}

		bool rand_point(Atom_Data<T> &atoms, r3d<T> &r_o)
		{
			for(auto itrial = 0; itrial < m_ntrial; itrial++)
			{
				auto r = rand_3d();
				if(box.check_r_min(atoms, r))
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