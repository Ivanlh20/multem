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
#include "hMT_General_CPU.h"
#include "hMatlab2Cpp.h"
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
int PotPar, nAtomTypes;
	sAtomTypesCPU *AtomTypes=0;

PotPar = (int)mxGetScalar(prhs[0]);

	nAtomTypes = NE;
	AtomTypes = new sAtomTypesCPU[nAtomTypes];
	f_SetAtomTypes(PotPar, 0, stVrl, nAtomTypes, AtomTypes);
	f_AtomTypesCPU2Matlab(nAtomTypes, AtomTypes, plhs[0]);
 
delete [] AtomTypes; AtomTypes = 0;
}