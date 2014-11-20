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
#include "hMT_General_GPU.h"
#include "hMT_Potential_GPU.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include <device_functions.h>
#include "cufft.h"
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int iConfFP, iSlice;
	int nAtomsM;
	double *AtomsM;
	cMGP MGP;
	sGP GP;

	nAtomsM = (int)mxGetM(prhs[0]); 
	AtomsM = mxGetPr(prhs[0]);
	MGP.gpu = (int)mxGetScalar(prhs[1]); 
	MGP.nx = (int)mxGetScalar(prhs[2]); 
	MGP.ny = (int)mxGetScalar(prhs[3]);
	MGP.lx = mxGetScalar(prhs[4]); 
	MGP.ly = mxGetScalar(prhs[5]);
	MGP.dz = mxGetScalar(prhs[6]);
	iConfFP = (int)mxGetScalar(prhs[7]); 
	MGP.DimFP = (int)mxGetScalar(prhs[8]); 
	MGP.SeedFP = (int)mxGetScalar(prhs[9]);
	iSlice = (int)mxGetScalar(prhs[10])-1; 
	if(iSlice<0)
		iSlice = 0;

	/*************************Output data**************************/
	plhs[0] = mxCreateDoubleMatrix(MGP.ny, MGP.nx, mxREAL);
	double *V0h = mxGetPr(plhs[0]);

	plhs[1] = mxCreateDoubleMatrix(MGP.ny, MGP.nx, mxREAL);
	double *V1h = mxGetPr(plhs[1]);

	cudaSetDevice(MGP.gpu);
	f_sGP_Cal(MGP.nx, MGP.ny, MGP.lx, MGP.ly, MGP.dz, MGP.PBC_xy, MGP.BWL, GP);

	cMT_Potential_GPU MT_Potential_GPU;	
	MT_Potential_GPU.SetInputData(MGP, GP, nAtomsM, AtomsM);
	MT_Potential_GPU.MT_Specimen_CPU.MoveAtoms(iConfFP);
	MT_Potential_GPU.ProjectedPotential(iSlice);

	f_fft2Shift_MD(GP, MT_Potential_GPU.V0, MT_Potential_GPU.V1);
	cudaMemcpy(V0h, MT_Potential_GPU.V0, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
	cudaMemcpy(V1h, MT_Potential_GPU.V1, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);

	//MT_Potential_GPU.freeMemory();
	//cudaDeviceReset();
}