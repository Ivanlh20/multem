/*
 * This file is part of MULTEM.
 * Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
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

#include <cstring>
#include "hmathCPU.h"
#include "hConstTypes.h"
#include "hMT_AtomTypes_CPU.h"
#include "hfxg_CPU.h"

// Constructor
cfxg_CPU::cfxg_CPU(){
	PotPar = 0;
	MT_AtomTypes_CPU = 0;
	cl = cnl = 0;
}

// Destructor
cfxg_CPU::~cfxg_CPU(){
	PotPar = 0;
	MT_AtomTypes_CPU = 0;
	cl = cnl = 0;
}

// Set Atom type
void cfxg_CPU::SetAtomT(int PotPar_i, cMT_AtomTypes_CPU *MT_AtomTypes_CPU_i){
	PotPar = PotPar_i;
	MT_AtomTypes_CPU = MT_AtomTypes_CPU_i;
}

// x-ray scattering factor(fx, dfx) where dfx is the first derivative along g
void cfxg_CPU::fxg(double g, double &f, double &df){
	int i;
	g2 = g*g;
	f = df = 0.0;
	switch (PotPar){
		case 1:
			// 1: Doyle and Tugner pagametegization - 4 Gaussians - [0, 4]
			for (i=0; i<4; i++){
				cl = MT_AtomTypes_CPU->cfxg.cl[i]; cnl = MT_AtomTypes_CPU->cfxg.cnl[i];
				f += ft = cl*exp(-cnl*g2);
				df += ft*(1.0-cnl*g2);
			}
			f = MT_AtomTypes_CPU->Z - g2*f;
			df = -2.0*g*df;
			break;
		case 2:
			// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
			for (i=0; i<5; i++){
				cl = MT_AtomTypes_CPU->cfxg.cl[i]; cnl = MT_AtomTypes_CPU->cfxg.cnl[i];
				f += ft = cl*exp(-cnl*g2);
				df += ft*(1.0-cnl*g2);
			}
			f = MT_AtomTypes_CPU->Z - g2*f;
			df = -2.0*g*df;
			break;
		case 3:
			// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
			for (i=0; i<5; i++){
				cl = MT_AtomTypes_CPU->cfxg.cl[i]; cnl = MT_AtomTypes_CPU->cfxg.cnl[i];
				f += ft = cl*exp(-cnl*g2);
				df += ft*(1.0-cnl*g2);
			}
			f = MT_AtomTypes_CPU->Z - g2*f;
			df = -2.0*g*df;
			break;
		case 4:
			// 4: Kirkland parameterization - 3 Yukawa + 3 Gaussians - [0, 12]					
			for (i=0; i<3;i++){
				cl = MT_AtomTypes_CPU->cfxg.cl[i]; cnl = MT_AtomTypes_CPU->cfxg.cnl[i];
				t = 1.0/(cnl + g2);
				f += ft = cl*t;
				df += ft*cnl*t;
			}

			for (i=3; i<6; i++){
				cl = MT_AtomTypes_CPU->cfxg.cl[i]; cnl = MT_AtomTypes_CPU->cfxg.cnl[i];
				f += ft = cl*exp(-cnl*g2);
				df += ft*(1.0-cnl*g2);
			}
			f = MT_AtomTypes_CPU->Z - g2*f;
			df = -2.0*g*df;
			break;
		case 5:
			// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
			for (i=0; i<6;i++){
				cl = MT_AtomTypes_CPU->cfxg.cl[i]; cnl = MT_AtomTypes_CPU->cfxg.cnl[i];
				f += ft = cl*exp(-cnl*g2);
				df += cnl*ft;
			}
			df = -2.0*g*df;
			break;
		case 6:
			// 6: Lobato parameterization - 5 Hydrogen fe - [0, 12]
			for (i=0; i<5; i++){
				cl = MT_AtomTypes_CPU->cfxg.cl[i]; cnl = MT_AtomTypes_CPU->cfxg.cnl[i];
				t = 1.0/(1.0+cnl*g2);
				f += ft = cl*t*t;
				df += cnl*ft*t;
			}
			df = -4.0*g*df;
			break;
	}
}

/***************************************************************************/
// 3D electron scattering factors calculation (feg, dfeg) where dfeg is the first derivative along g
void cfxg_CPU::fxg(int ng, double *g, double *f, double *df){
	for (int i=0; i<ng; i++)
		fxg(g[i], f[i], df[i]);
}