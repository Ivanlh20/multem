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

#include <cstring>
#include "hConstTypes.h"
#include "hMT_General_CPU.h"
#include "hMT_Specimen_CPU.h"
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int nAtomsM, iConfFP, nAtoms, nPlanesu;
	double *AtomsM, *Planesu;
	sMGP MGP;

	f_sMPG_Init(MGP);

	nAtomsM = (int)mxGetM(prhs[0]); 
	AtomsM = mxGetPr(prhs[0]);
	MGP.lx = mxGetScalar(prhs[1]); 
	MGP.ly = mxGetScalar(prhs[2]); 
	iConfFP = (int)mxGetScalar(prhs[3]); 
	MGP.DimFP = (int)mxGetScalar(prhs[4]); 
	MGP.SeedFP = (int)mxGetScalar(prhs[5]);

	/**************************Input data**************************/
	cMT_Specimen_CPU MT_Specimen_CPU;
	MT_Specimen_CPU.SetInputData(MGP, nAtomsM, AtomsM);
	MT_Specimen_CPU.MoveAtoms(iConfFP);

	nPlanesu = MT_Specimen_CPU.nPlanesu;
	/*************************Output data**************************/
	plhs[0] = mxCreateDoubleMatrix(nPlanesu, 1, mxREAL);
	Planesu = mxGetPr(plhs[0]);
	memcpy(Planesu, MT_Specimen_CPU.Planesu, nPlanesu*cSizeofRD);
}