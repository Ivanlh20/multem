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

#ifndef MULTEM_ATOM_CAL_H
#define MULTEM_ATOM_CAL_H

#include <functional>

#include <multem/config.h>
#include <multem/math.h>
#include <multem/types.h>

namespace mt
{
	template <class T>
	class Atom_Cal{
		public:
			using value_type = T;

			DLL_PUBLIC Atom_Cal();

			DLL_PUBLIC void Set_Atom_Type(const ePotential_Type &PotPar_i, const int &charge_i, Atom_Type<T, e_host> *atom_type_CPU_i);

			DLL_PUBLIC inline void feg(const T &g, T &y);
			DLL_PUBLIC void feg(const int &ng, T *g, T *y);

      DLL_PUBLIC inline void feg_dfeg(const T &g, T &y, T &dy);
			DLL_PUBLIC void feg_dfeg(const int &ng, T *g, T *y, T *dy);

      DLL_PUBLIC inline void fxg(const T &g, T &y);
			DLL_PUBLIC void fxg(const int &ng, T *g, T *y);

			DLL_PUBLIC inline void fxg_dfxg(const T &g, T &y, T &dy);
			DLL_PUBLIC void fxg_dfxg(const int &ng, T *g, T *y, T *dy);

			DLL_PUBLIC inline void Pr(const T &r, T &y);
			DLL_PUBLIC void Pr(const int &nr, T *r, T *y);

			DLL_PUBLIC inline void Pr_dPr(const T &r, T &y, T &dy);
			DLL_PUBLIC void Pr_dPr(const int &nr, T *r, T *y, T *dy);

			DLL_PUBLIC inline void Vr(const T &r, T &y);
			DLL_PUBLIC void Vr(const int &nr, T *r, T *y);

			DLL_PUBLIC inline void Vr_dVr(const T &r, T &y, T &dy);
			DLL_PUBLIC void Vr_dVr(const int &nr, T *r, T *y, T *dy);

			DLL_PUBLIC inline void VR(const T &R, T &y);
			DLL_PUBLIC void VR(const int &nR, T *R, T *y);

			DLL_PUBLIC inline void VR_dVR(const T &R, T &y, T &dy);
			DLL_PUBLIC void VR_dVR(const int &nR, T *R, T *y, T *dy);

			DLL_PUBLIC inline void Vz(const T &z0, const T &ze, const T &R, T &y);
			DLL_PUBLIC void Vz(const T &z0, const T &ze, const int &nR, T *R, T *y);

			DLL_PUBLIC inline void Vz_dVz(const T &z0, const T &ze, const T &R, T &y, T &dy);
			DLL_PUBLIC void Vz_dVz(const T &z0, const T &ze, const int &nR, T *R, T *y, T *dy);

			DLL_PUBLIC T AtomicRadius_rms(const int &Dim);
			DLL_PUBLIC T AtomicRadius_Cutoff(const int &Dim, const T &Vrl);

		private:
			ePotential_Type potential_type;
			int charge;

			Atom_Type<T, e_host> *atom_type;
			PP_Coef<T, e_host>	*pr_feg;
			PP_Coef<T, e_host>	*pr_fxg;
			PP_Coef<T, e_host>	*pr_Pr;
			PP_Coef<T, e_host>	*pr_Vr;
			PP_Coef<T, e_host>	*pr_VR;

			Q1<T, e_host> Qz_a_b;
			Q1<T, e_host> Qz_0_I;
	};

} // namespace mt

#endif

