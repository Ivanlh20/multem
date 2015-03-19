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
#include "hPotential_CPU.h"
#include <mex.h>
	
void mexFunction(int nlhs,mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int Z, PotPar, IntType, Dim;
	double sigma, *r, *Vr, *dVr;
	int m, n, nr;

	PotPar = (int)mxGetScalar(prhs[0]);
	Z = (int)mxGetScalar(prhs[1]);
	sigma = mxGetScalar(prhs[2]);
	IntType = (int)mxGetScalar(prhs[3]);
	Dim = (int)mxGetScalar(prhs[4]);
	r = mxGetPr(prhs[5]);
	m = (int)mxGetM(prhs[5]); 
	n = (int)mxGetN(prhs[5]); 
	nr = m*n;
	plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
	Vr = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
	dVr = mxGetPr(plhs[1]);

	cMT_AtomTypes_CPU MT_AtomTypes_CPU;
	MT_AtomTypes_CPU.SetAtomTypes(Z, PotPar, stVrl, stnR, 0);

	cPotential_CPU Potential_CPU;
	Potential_CPU.SetAtomTypes(PotPar, &MT_AtomTypes_CPU);
	Potential_CPU.SetSigma(sigma);
	Potential_CPU.Vr(IntType, Dim, nr, r, Vr, dVr);
}