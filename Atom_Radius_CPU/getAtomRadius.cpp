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

#include "hConstTypes.h"
#include "hMT_AtomTypes_CPU.h"
#include "hMT_General_CPU.h"
#include "hPotential_CPU.h"
#include <mex.h>

void GetAtomicRadius(int PotPar, int Dim, double Vrl, double *r){
	int nMT_AtomTypes;
	cMT_AtomTypes_CPU *MT_AtomTypes_CPU;

	nMT_AtomTypes = NE;
	MT_AtomTypes_CPU = new cMT_AtomTypes_CPU[nMT_AtomTypes];
	for(int i=0; i<nMT_AtomTypes; i++)
		MT_AtomTypes_CPU[i].SetAtomTypes(i+1, PotPar, Vrl, stnR, 0);

	cPotential_CPU Potential_CPU;
	for (int i=0; i<NE; i++){	
		Potential_CPU.SetAtomTypes(PotPar, &(MT_AtomTypes_CPU[i]));
		Potential_CPU.SetSigma(0.0);
		r[i] = Potential_CPU.AtomicRadius_rms(Dim);
		r[i+NE] = Potential_CPU.AtomicRadius_Cutoff(Dim, Vrl);
		r[i+2*NE] = MT_AtomTypes_CPU[i].ra_e;
	}

	delete [] MT_AtomTypes_CPU; MT_AtomTypes_CPU = 0;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	int PotPar, Dim;
	double *r, Vrl;	
	PotPar = (int)mxGetScalar(prhs[0]);
	Dim = (int)mxGetScalar(prhs[1]);
	Vrl = mxGetScalar(prhs[2]);

	plhs[0] = mxCreateDoubleMatrix(NE, 3, mxREAL);
	r = mxGetPr(plhs[0]);

	GetAtomicRadius(PotPar, Dim, Vrl, r);
}