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

#ifndef ATOMIC_FCNS_MT_H
#define ATOMIC_FCNS_MT_H

#ifdef _MSC_VER
#pragma once
#endif// _MSC_VER

#include <functional>

#include "math.cuh"
#include "types.cuh"
#include "quadrature.hpp"
#include "cgpu_fcns.cuh"

namespace mt
{
	template <class T>
	class Atom_Cal{
		public:
			using value_type = T;

			Atom_Cal(): potential_type(ePT_Lobato_0_12), charge(0), atom_type(nullptr), 
			pr_feg(nullptr), pr_fxg(nullptr), pr_Pr(nullptr), pr_Vr(nullptr), pr_VR(nullptr)
			{
				Quadrature quadrature;
				quadrature(0, c_nqz, Qz_a_b); 	// 0: int_-1^1 y(x) dx - TanhSinh quadrature
				quadrature(1, c_nqz, Qz_0_I); 	// 1: int_0^infty y(x) dx - ExpSinh quadrature
			}

			void Set_Atom_Type(const ePotential_Type &PotPar_i, const int &charge_i, Atom_Type<T, e_host> *atom_type_CPU_i)
			{
				potential_type = PotPar_i;
				charge = atom_type_CPU_i->check_charge(charge_i);
				atom_type = atom_type_CPU_i;
				pr_feg = atom_type->feg(charge);
				pr_fxg = atom_type->fxg(charge);
				pr_Pr = atom_type->Pr(charge);
				pr_Vr = atom_type->Vr(charge);
				pr_VR = atom_type->VR(charge);
			}

			// Electron scattering factors calculation (feg)
			inline void feg(const T &g, T &y)
			{
				mt::feg<T>(potential_type, charge, g, *pr_feg, y);
			}

			void feg(const int &ng, T *g, T *y)
			{
				for(auto i = 0; i < ng; i++)
				{
					feg(g[i], y[i]);
				}
			}

			// Electron scattering factor(feg, dfeg) where dfg is the first derivative along g
			inline void feg_dfeg(const T &g, T &y, T &dy)
			{
				mt::feg_dfeg<T>(potential_type, charge, g, *pr_feg, y, dy);
			}

			void feg_dfeg(const int &ng, T *g, T *y, T *dy)
			{
				for(auto i = 0; i < ng; i++)
				{
					feg_dfeg(g[i], y[i], dy[i]);
				}
			}

			// Electron scattering factor(fg)
			inline void fxg(const T &g, T &y)
			{
				mt::fxg<T>(potential_type, charge, atom_type->Z, g, *pr_fxg, y);
			}

			void fxg(const int &ng, T *g, T *y)
			{
				for(auto i = 0; i < ng; i++)
				{
					fxg(g[i], y[i]);
				}
			}

			// Electron scattering factor(fg, dfg) where dfg is the first derivative along g
			inline void fxg_dfxg(const T &g, T &y, T &dy)
			{
				mt::fxg_dfxg<T>(potential_type, charge, atom_type->Z, g, *pr_fxg, y, dy);
			}

			void fxg_dfxg(const int &ng, T *g, T *y, T *dy)
			{
				for(auto i = 0; i < ng; i++)
				{
					fxg_dfxg(g[i], y[i], dy[i]);
				}
			}

			// Electron density (Pr)
			inline void Pr(const T &r, T &y)
			{
				mt::Pr<T>(potential_type, charge, r, *pr_Pr, y);
			}

			void Pr(const int &nr, T *r, T *y)
			{
				for(auto i = 0; i < nr; i++)
				{
					Pr(r[i], y[i]);
				}
			}

			// Electron density (Pr, dPr) where dPr is the first derivative along r
			inline void Pr_dPr(const T &r, T &y, T &dy)
			{
				mt::Pr_dPr<T>(potential_type, charge, r, *pr_Pr, y, dy);
			}

			void Pr_dPr(const int &nr, T *r, T *y, T *dy)
			{
				for(auto i = 0; i < nr; i++)
				{
					Pr_dPr(r[i], y[i], dy[i]);
				}
			}

			// Projected_Potential calculation(Vr)
			inline void Vr(const T &r, T &y)
			{
				mt::Vr<T>(potential_type, charge, r, *pr_Vr, y);
			}

			void Vr(const int &nr, T *r, T *y)
			{
				for(auto i = 0; i < nr; i++)
				{
					Vr(r[i], y[i]);
				}
			}

			// Projected_Potential calculation (Vr, dVr) where dVr is the first derivative along r
			inline void Vr_dVr(const T &r, T &y, T &dy)
			{
				mt::Vr_dVr<T>(potential_type, charge, r, *pr_Vr, y, dy);
			}

			void Vr_dVr(const int &nr, T *r, T *y, T *dy)
			{
				for(auto i = 0; i < nr; i++)
				{
					Vr_dVr(r[i], y[i], dy[i]);
				}
			}

			// Projected potential (VR)
			inline void VR(const T &R, T &y)
			{
				mt::VR<T>(potential_type, charge, R, *pr_VR, Qz_0_I, y);
			}

			void VR(const int &nR, T *R, T *y)
			{
				for(auto i = 0; i < nR; i++)
				{
					VR(R[i], y[i]);
				}
			}

			// Projected potential (VR, dVR) where dVr is the first derivative along R
			inline void VR_dVR(const T &R, T &y, T &dy)
			{
				mt::VR_dVR<T>(potential_type, charge, R, *pr_VR, Qz_0_I, y, dy);
			}

			void VR_dVR(const int &nR, T *R, T *y, T *dy)
			{
				for(auto i = 0; i < nR; i++)
				{
					VR_dVR(R[i], y[i], dy[i]);
				}
			}

			// Projected potential (Vz)[z0, ze]
			inline void Vz(const T &z0, const T &ze, const T &R, T &y)
			{
				mt::Vz<T>(potential_type, charge, z0, ze, R, *pr_Vr, Qz_a_b, y);
			}

			void Vz(const T &z0, const T &ze, const int &nR, T *R, T *y)
			{
				for(auto i = 0; i < nR; i++)
				{
					Vz(z0, ze, R[i], y[i]);
				}
			}

			// Projected potential (Vz, dVz)[z0, ze] where dVr is the first derivative along R
			inline void Vz_dVz(const T &z0, const T &ze, const T &R, T &y, T &dy)
			{
				mt::Vz_dVz<T>(potential_type, charge, z0, ze, R, *pr_Vr, Qz_a_b, y, dy);
			}

			void Vz_dVz(const T &z0, const T &ze, const int &nR, T *R, T *y, T *dy)
			{
				for(auto i = 0; i < nR; i++)
				{
					Vz_dVz(z0, ze, R[i], y[i], dy[i]);
				}
			}

			T AtomicRadius_rms(const int &Dim)
			{
		
				if(isZero(pr_Vr->cl[0]))
				{
					return 0.0;
				}
				else
				{
					auto fx = [=](const T &x)->T
					{	
						T V;
						if(Dim == 3)
						{
							this->Vr(x, V);
						}
						else 
						{
							this->VR(x, V);
						}
						return V;
					};
					
					auto rmin = atom_type->rn_c;
					auto Vri = fx(rmin);

					auto iVr = Vri*pow(rmin, Dim)/Dim;
					auto iVrn = Vri*pow(rmin, Dim+2)/(Dim+2);

					for(auto i = 0; i< Qz_0_I.size(); i++)
					{
						auto ri = Qz_0_I.x[i] + rmin;
						Vri = fx(ri);
						T Vrt;
						iVr += Vrt = Qz_0_I.w[i]*Vri*pow(ri, Dim-1);
						iVrn += Vrt*ri*ri;
					}
					return sqrt(iVrn/iVr);
				}
			}

			T AtomicRadius_Cutoff(const int &Dim, const T &Vrl)
			{
				if(!nonZero(Vrl, pr_Vr->cl[0]))
				{
					return 0;
				}
				else
				{
					auto f = [=](const T &x)->T
					{	
						T V;
						if(Dim == 3)
						{
							this->Vr(x, V);
						}
						else 
						{
							this->VR(x, V);
						}
						return V-Vrl;
					};

					T rmin = atom_type->rn_c; 
					T rmax = 25.0;

					return host_device_detail::Root_Finder(f, rmin, rmax);
				}
			}

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