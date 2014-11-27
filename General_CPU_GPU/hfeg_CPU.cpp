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

#include "hmathCPU.h"
#include <cstring>
#include "hConstTypes.h"
#include "hfeg_CPU.h"

// Constructor
cfeg_CPU::cfeg_CPU(){
	PotPar = 0;
	cl = cnl = 0;
}

// Destructor
cfeg_CPU::~cfeg_CPU(){
	PotPar = 0;
	cl = cnl = 0;
}

// Set Atom type
void cfeg_CPU::SetAtomT(int PotPari, sAtomTypesCPU AtomTypesCPUi){
	PotPar = PotPari;
	memcpy(&AtomTypesCPU, &AtomTypesCPUi, sizeof(sAtomTypesCPU));
}

// Electron scattering factor(fg, dfg) where dfg is the first derivative along g
void cfeg_CPU::feg(double g, double &f, double &df){
	int i;
	g2 = g*g;
	f = df = 0.0;
	switch (PotPar){
		case 1:
			// 1: Doyle and Tugner pagametegization - 4 Gaussians - [0, 4]
			for (i=0; i<4; i++){
				cl = AtomTypesCPU.cfeg.cl[i]; cnl = AtomTypesCPU.cfeg.cnl[i];
				f += ft = cl*exp(-cnl*g2);
				df += cnl*ft;
			}
			df = -2.0*g*df;
			break;
		case 2:
			// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
			for (i=0; i<5; i++){
				cl = AtomTypesCPU.cfeg.cl[i]; cnl = AtomTypesCPU.cfeg.cnl[i];
				f += ft = cl*exp(-cnl*g2);
				df += cnl*ft;
			}
			df = -2.0*g*df;
			break;
		case 3:
			// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
			for (i=0; i<5; i++){
				cl = AtomTypesCPU.cfeg.cl[i]; cnl = AtomTypesCPU.cfeg.cnl[i];
				f += ft = cl*exp(-cnl*g2);
				df += cnl*ft;
			}
			df = -2.0*g*df;
			break;
		case 4:
			// 4: Kirkland parameterization - 3 Yukawa + 3 Gaussians - [0, 12]					
			for (i=0; i<3;i++){
				cl = AtomTypesCPU.cfeg.cl[i]; cnl = AtomTypesCPU.cfeg.cnl[i];
				t = 1.0/(cnl + g2);
				f += ft = cl*t;
				df += ft*t;
			}

			for (i=3; i<6; i++){
				cl = AtomTypesCPU.cfeg.cl[i]; cnl = AtomTypesCPU.cfeg.cnl[i];
				f += ft = cl*exp(-cnl*g2);
				df += cnl*ft;
			}
			df = -2.0*g*df;
			break;
		case 5:
			// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
			if (g != 0){
				for (i=0; i<6;i++){
					cl = AtomTypesCPU.cfeg.cl[i]; cnl = AtomTypesCPU.cfeg.cnl[i];
					t = exp(-cnl*g2);
					f += ft = cl*(1.0-t);
					df += cl*(1.0-(1.0+cnl*g2)*t);
				}
				f = f/g2;
				df = -2.0*df/(g*g2);
			}else
				for (i=0; i<6;i++){
					cl = AtomTypesCPU.cfeg.cl[i]; cnl = AtomTypesCPU.cfeg.cnl[i];
					f += cl*cnl;
				}
			break;
		case 6:
			// 6: Lobato parameterization - 5 Hydrogen fe - [0, 12]
			for (i=0; i<5; i++){
				cl = AtomTypesCPU.cfeg.cl[i]; cnl = AtomTypesCPU.cfeg.cnl[i];
				t = 1.0/(1.0 + cnl*g2);
				f += cl*t*(t + 1.0);
				df += cl*cnl*t*t*(2.0*t + 1.0);
			}
			df = -2.0*g*df;
			break;
	}
}

/***************************************************************************/
// 3D electron scattering factors calculation (feg, dfeg) where dfeg is the first derivative along g
void cfeg_CPU::feg(int ng, double *g, double *f, double *df){
	for (int i=0; i<ng; i++)
		feg(g[i], f[i], df[i]);
}