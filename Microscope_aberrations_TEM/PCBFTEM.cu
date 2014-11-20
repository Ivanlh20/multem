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

#include "math.h"
#include "hConstTypes.h"
#include "hMT_General_GPU.h"
#include "hMatlab2Cpp.h"
#include "hTEMIm.h"

#include "cuda.h"
#include "cuda_runtime.h"

#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	double *M2Psi_h;
	sInTEMIm InTEMIm;
	cTEMIm TEMIm;

	Matlab2InTEMIm(prhs[0], InTEMIm);
	TEMIm.SetInputData(InTEMIm);
	plhs[0] = mxCreateDoubleMatrix((mwSize)TEMIm.ny, (mwSize)TEMIm.nx, mxREAL);
	M2Psi_h = mxGetPr(plhs[0]);
	TEMIm.TEMImage(TEMIm.Psirh, TEMIm.Psiih, M2Psi_h);
 	/*************************************************************/
	TEMIm.freeMemory();
}