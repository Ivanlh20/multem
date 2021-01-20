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

#ifndef MULTEM_AMORP_SPEC_H
#define MULTEM_AMORP_SPEC_H

#include <multem/math.h>
#include <multem/types.h>
#include <multem/lin_alg_def.h>
#include <multem/config.h>
#include <multem/atom_data.h>
#include <multem/atomic_data.h>
#include <multem/cgpu_rand.h>
#include <multem/box_occ.h>

namespace mt
{
	template <class T>
	class Amorp_Spec {
	public:
		
    DLL_PUBLIC 
    Amorp_Spec();

		DLL_PUBLIC 
    void create(Atom_Data<T> &atoms, T d_min, int Z, T rms_3d, T rho, int seed = 300183);

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

		DLL_PUBLIC int set_init_values(Atom_Data<T> &atoms, T d_min, int Z, T rms_3d, T rho, int seed, T depth);
		DLL_PUBLIC T convert_density(int Z, T m_rho);
		DLL_PUBLIC T vol() const;
		DLL_PUBLIC T density(int natoms) const;
		DLL_PUBLIC bool rand_point(Atom_Data<T> &atoms, r3d<T> &r_o);
	};
}
#endif

