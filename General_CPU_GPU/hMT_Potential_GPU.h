/**
 *  This file is part of MULTEM.
 *  Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
 *
 *  MULTEM is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MULTEM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MULTEM.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef hMT_Potential_GPU_H
#define hMT_Potential_GPU_H

#include "hConstTypes.h"
#include "hMT_Specimen_CPU.h"
#include "hMT_AtomTypes_GPU.h"

class cMT_Potential_GPU{
	private:
		sQ1 Qz;

		int ncVph;
		scVp *cVph;

		sMGP MGP;						// Multislice general parameters
		sGP GP;							// xy-Grid properties

		int nAtomTypesGPU;				// Number of atom types
		cMT_AtomTypes_GPU *AtomTypesGPU;	// Atom types

		void SetPotPar(int PotParh);
		int CheckGridLimits(int i, int n);
		void getbn(sGP &GP, double x, double y, double Rmax, sbn &bnx, sbn &bny);
		void setcVp(int ApproxModel, cMT_Specimen_CPU &MT_Specimen_CPU, cMT_AtomTypes_GPU *&AtomTypesGPU, int iSlice, int iatom, int nsatom, dim3 &BPot, dim3 &TPot, dim3 &BCoef, dim3 &TCoef);
	public:
		cMT_Specimen_CPU MT_Specimen_CPU;	// MT_Specimen_CPU

		double *V0;				// Zero moment of the potential
		double *V1;				// first moment of the potential
		double *V2;				// Second moment of the potential

		void freeMemory();	
		cMT_Potential_GPU();
		~cMT_Potential_GPU();

		void SetInputData(sMGP &MGP_io, sGP &GP_i, int nAtomsM_i, double *AtomsM_i);
		void ProjectedPotential(int iSlice);
};

#endif