/**
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

#include <cstring>
#include "types.cuh"
#include "input_multislice.cuh"
#include "specimen.hpp"
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int natomsM, iConfFP, nAtoms, nPlanes_u;
	double *atomsM, *Planes_u;
	Input_Multislice input_multislice;

	natomsM = mxGetM(prhs[0]); 
	atomsM = mxGetPr(prhs[0]);
	input_multislice.lx = mxGetScalar(prhs[1]); 
	input_multislice.ly = mxGetScalar(prhs[2]); 
	iConfFP = (int)mxGetScalar(prhs[3]); 
	input_multislice.fp_dim = (int)mxGetScalar(prhs[4]); 
	input_multislice.fp_seed = (int)mxGetScalar(prhs[5]);

	/**************************Input data**************************/
	Specimen specimen;
	specimen.set_input_data(&input_multislice, natomsM, atomsM);
	specimen.move_atoms(iConfFP);

	nPlanes_u = specimen.nPlanes_u;
	/*************************Output data**************************/
	plhs[0] = mxCreateDoubleMatrix(nPlanes_u, 1, mxREAL);
	Planes_u = mxGetPr(plhs[0]);
	memcpy(Planes_u, specimen.Planes_u, nPlanes_u*cSizeofRD);
}