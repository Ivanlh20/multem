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

#ifndef ATOM_CAL_H
#define ATOM_CAL_H

#include <functional>

#include "math.cuh"
#include "types.hpp"
#include "quadrature.hpp"
#include "host_device_functions.cuh"

namespace multem
{
	template<class T>
	class Atom_Cal{
		public:
			using value_type = typename T;

			Atom_Cal():potential_type(ePT_Lobato_0_12)
			{
				Quadrature quadrature;
				quadrature.get(0, c_nqz, Qz_a_b); 	// 0: int_-1^1 y(x) dx - TanhSinh quadrature
				quadrature.get(1, c_nqz, Qz_0_I); 	// 1: int_0^infty y(x) dx - ExpSinh quadrature
			}

			void Set_Atom_Type(const ePotential_Type &PotPar_i, Atom_Type<T, Host> *atom_type_CPU_i)
			{
				potential_type = PotPar_i;
				atom_type = atom_type_CPU_i;
			}

			// Electron scattering factors calculation (feg)
			void feg(const T &g, T &y)
			{
				multem::feg<T>(potential_type, g, atom_type->feg, y);
			}

			void feg(const int &ng, T *g, T *y)
			{
				for(auto i = 0; i < ng; i++)
				{
					feg(g[i], y[i]);
				}
			}

			// Electron scattering factor(feg, dfeg) where dfg is the first derivative along g
			void feg_dfeg(const T &g, T &y, T &dy)
			{
				multem::feg_dfeg<T>(potential_type, g, atom_type->feg, y, dy);
			}

			void feg_dfeg(const int &ng, T *g, T *y, T *dy)
			{
				for(auto i = 0; i < ng; i++)
				{
					feg_dfeg(g[i], y[i], dy[i]);
				}
			}

			// Electron scattering factor(fg)
			void fxg(const T &g, T &y)
			{
				multem::fxg<T>(potential_type, g, atom_type->Z, atom_type->fxg, y);
			}

			void fxg(const int &ng, T *g, T *y)
			{
				for(auto i = 0; i < ng; i++)
					fxg(g[i], y[i]);
			}

			// Electron scattering factor(fg, dfg) where dfg is the first derivative along g
			void fxg_dfxg(const T &g, T &y, T &dy)
			{
				multem::fxg_dfxg<T>(potential_type, g, atom_type->Z, atom_type->fxg, y, dy);
			}

			void fxg_dfxg(const int &ng, T *g, T *y, T *dy)
			{
				for(auto i = 0; i < ng; i++)
					fxg_dfxg(g[i], y[i], dy[i]);
			}

			// Electron density (Pr)
			void Pr(const T &r, T &y)
			{
				multem::Pr<T>(potential_type, r, atom_type->Pr, y);
			}

			void Pr(const int &nr, T *r, T *y)
			{
				for(auto i = 0; i < nr; i++)
					Pr(r[i], y[i]);
			}

			// Electron density (Pr, dPr) where dPr is the first derivative along r
			void Pr_dPr(const T &r, T &y, T &dy)
			{
				multem::Pr_dPr<T>(potential_type, r, atom_type->Pr, y, dy);
			}

			void Pr_dPr(const int &nr, T *r, T *y, T *dy)
			{
				for(auto i = 0; i < nr; i++)
					Pr_dPr(r[i], y[i], dy[i]);
			}

			// Potential calculation(Vr)
			void Vr(const T &r, T &y)
			{
				multem::Vr<T>(potential_type, r, atom_type->Vr, y);
			}

			void Vr(const int &nr, T *r, T *y)
			{
				for(auto i = 0; i < nr; i++)
					Vr(r[i], y[i]);
			}

			// Potential calculation (Vr, dVr) where dVr is the first derivative along r
			void Vr_dVr(const T &r, T &y, T &dy)
			{
				multem::Vr_dVr<T>(potential_type, r, atom_type->Vr, y, dy);
			}

			void Vr_dVr(const int &nr, T *r, T *y, T *dy)
			{
				for(auto i = 0; i < nr; i++)
					Vr_dVr(r[i], y[i], dy[i]);
			}

			// Projected potential (VR)
			void VR(const T &R, T &y)
			{
				multem::VR<T>(potential_type, R, atom_type->VR, Qz_0_I, y);
			}

			void VR(const int &nR, T *R, T *y)
			{
				for(auto i = 0; i < nR; i++)
					VR(R[i], y[i]);
			}

			// Projected potential (VR, dVR) where dVr is the first derivative along R
			void VR_dVR(const T &R, T &y, T &dy)
			{
				multem::VR_dVR<T>(potential_type, R, atom_type->VR, Qz_0_I, y, dy);
			}

			void VR_dVR(const int &nR, T *R, T *y, T *dy)
			{
				for(auto i = 0; i < nR; i++)
					VR_dVR(R[i], y[i], dy[i]);
			}

			// Projected potential (Vz)[z0, ze]
			void Vz(const T &z0, const T &ze, const T &R, T &y)
			{
				multem::Vz<T>(potential_type, z0, ze, R, atom_type->Vr, Qz_a_b, y);
			}

			void Vz(const T &z0, const T &ze, const int &nR, T *R, T *y)
			{
				for(auto i = 0; i < nR; i++)
					Vz(R[i], y[i]);
			}

			// Projected potential (Vz, dVz)[z0, ze] where dVr is the first derivative along R
			void Vz_dVz(const T &z0, const T &ze, const T &R, T &y, T &dy)
			{
				multem::Vz_dVz<T>(potential_type, z0, ze, R, atom_type->Vr, Qz_a_b, y, dy);
			}

			void Vz_dVz(const T &z0, const T &ze, const int &nR, T *R, T *y, T *dy)
			{
				for(auto i = 0; i < nR; i++)
					Vz_dVz(z0, ze, R[i], y[i], dy[i]);
			}

			T AtomicRadius_rms(const int &Dim)
			{
				if(isZero(atom_type->Vr.cl[0]))
				{
					return 0.0;
				}
				else
				{
					auto f = [&](const T &x)->T
					{	
						T V;
						if(Dim==3) Vr(x, V);
						else VR(x, V);
						return V;
					};

					T iVr, iVrn, Vrt;
					T rmin, ri, Vri;

					rmin = atom_type->rn_c;
					Vri = f(rmin);

					iVr = Vri*pow(rmin, Dim)/Dim;
					iVrn = Vri*pow(rmin, Dim+2)/(Dim+2);

					/************************************************************/
					for(auto i = 0; i< Qz_0_I.size(); i++)
					{
						ri = Qz_0_I.x[i] + rmin;
						Vri = f(ri);
						iVr += Vrt = Qz_0_I.w[i]*Vri*pow(ri, Dim-1);
						iVrn += Vrt*ri*ri;
					}
					return sqrt(iVrn/iVr);
				}
			}

			T AtomicRadius_Cutoff(const int &Dim, const T &Vrl)
			{
				if(!nonZero(Vrl, atom_type->Vr.cl[0]))
				{
					return 0;
				}
				else
				{
					auto f = [&](const T &x)->T
					{	
						T V;
						if(Dim==3) Vr(x, V);
						else VR(x, V);
						return V-Vrl;
					};

					T rmin = atom_type->rn_c; 
					T rmax = 25.0;

					return host_device_detail::Root_Finder(f, rmin, rmax);
				}
			}

		private:
			ePotential_Type potential_type;
			Atom_Type<T, Host> *atom_type;
			Q1<T, Host> Qz_a_b;
			Q1<T, Host> Qz_0_I;
	};

} // namespace multem

#endif