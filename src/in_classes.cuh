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

#ifndef IN_CLASSES_H
	#define IN_CLASSES_H

	#include <algorithm>

	#include "math.cuh"
	#include "types.cuh"
	#include "particles.cuh"
	#include "cgpu_info.cuh"
	#include "cgpu_vctr.cuh"
	#include "grid.cuh"

	/* xtl */
	namespace mt
	{
		/********************* Input xtlline specimen *********************/
		/* Asymmetric Unit
		The asymmetric unit is the smallest fraction of the unit cell that
		can be rotated and translated using only the symmetry operators
		allowed by the xtllographic symmetry to generate one unit cell.
		*/
		template <class T>
		class In_Xtl_Build
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

			In_Xtl_Build() :
				a(0), b(0), c(0), alpha(0), beta(0), gamma(0),
				n_a(0), n_b(0), n_c(0), sgn(1), pbc(false) {};
		};
	}

	/***************************************************************************************/
	/***************************** Input particle superposition ****************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T, eDim Dim, eFcn_typ Fcn_typ> class In_Spt_Ptc_xd;

		template <class T>
		using In_Spt_Ptc_2d = In_Spt_Ptc_xd<T, edim_2, efcn_gauss>;
	}

	namespace mt
	{
		 template <class T, eDim Dim>
		 class In_Spt_Ptc_xd<T, Dim, efcn_gauss>
		 {
	 		public:
	 			using value_type = T;

	 			T scf_radius;
	 			Grid_xd<T, Dim> grid;
	 			Ptc_2d_2<T> ptc;

	 			In_Spt_Ptc_xd(): scf_radius(3.5) {};

				Ptc_s_fcn_xd<T, Dim, efcn_gauss> ptc_s_fcn(const dt_int32& iptc, const T& tap_ftr)
				{
					const auto r_max = scf_radius*ptc.c_2[iptc];
					const auto r_tap = tap_ftr*r_max;

					return {ptc.get_pos(iptc), ptc.c_1[iptc], ptc.c_2[iptc], r_tap, r_max, grid};
				}
		 };
	}

#endif