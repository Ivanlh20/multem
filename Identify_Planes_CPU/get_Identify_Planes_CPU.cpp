/**
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
#include "hConstTypes.h"
#include "hMT_MGP_CPU.h"
#include "hMT_Specimen_CPU.h"
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int nAtomsM, iConfFP, nAtoms, nPlanesu;
	double *AtomsM, *Planesu;
	cMT_MGP_CPU MT_MGP_CPU;

	nAtomsM = (int)mxGetM(prhs[0]); 
	AtomsM = mxGetPr(prhs[0]);
	MT_MGP_CPU.lx = mxGetScalar(prhs[1]); 
	MT_MGP_CPU.ly = mxGetScalar(prhs[2]); 
	iConfFP = (int)mxGetScalar(prhs[3]); 
	MT_MGP_CPU.DimFP = (int)mxGetScalar(prhs[4]); 
	MT_MGP_CPU.SeedFP = (int)mxGetScalar(prhs[5]);

	/**************************Input data**************************/
	cMT_Specimen_CPU MT_Specimen_CPU;
	MT_Specimen_CPU.SetInputData(&MT_MGP_CPU, nAtomsM, AtomsM);
	MT_Specimen_CPU.MoveAtoms(iConfFP);

	nPlanesu = MT_Specimen_CPU.nPlanesu;
	/*************************Output data**************************/
	plhs[0] = mxCreateDoubleMatrix(nPlanesu, 1, mxREAL);
	Planesu = mxGetPr(plhs[0]);
	memcpy(Planesu, MT_Specimen_CPU.Planesu, nPlanesu*cSizeofRD);
}