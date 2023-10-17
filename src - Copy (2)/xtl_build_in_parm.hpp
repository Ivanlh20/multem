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

#ifndef XTL_BUILD_IN_PARM_H
	#define XTL_BUILD_IN_PARM_H

	#include "particles.cuh"

	/* xtl */
	namespace mt
	{
		/* Asymmetric Unit
		The asymmetric unit is the smallest fraction of the unit cell that
		can be rotated and translated using only the symmetry operators
		allowed by the xtllographic symmetry to generate one unit cell.
		*/
		template <class T>
		class Xtl_Build_In_Parm
		{
		public:
			using value_type = T;

			T a;					// lattice constant a
			T b;					// lattice constant b
			T c;					// lattice constant c

			T alpha;				// angle between b & c
			T beta;					// angle between c & a
			T gamma;				// angle between a & b

			dt_int32 n_a;
			dt_int32 n_b;
			dt_int32 n_c;

			dt_int32 sgn;			// space group number
			dt_bool pbc;

			Ptc_Atom<T> asym_uc;	// normalized positions in the asymmetric unit cell
			Ptc_Atom<T> base;		// normalized positions of atom in the unit cell

			Xtl_Build_In_Parm() :
				a(0), b(0), c(0), alpha(0), beta(0), gamma(0),
				n_a(0), n_b(0), n_c(0), sgn(1), pbc(false) {};
		};
	}
#endif