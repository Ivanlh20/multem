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
	cMT_MGP_CPU MT_MGP_CPU;

	nAtomsM = (int)mxGetM(prhs[0]); 
	AtomsM = mxGetPr(prhs[0]);
	MT_MGP_CPU.gpu = (int)mxGetScalar(prhs[1]); 
	MT_MGP_CPU.nx = (int)mxGetScalar(prhs[2]); 
	MT_MGP_CPU.ny = (int)mxGetScalar(prhs[3]);
	MT_MGP_CPU.lx = mxGetScalar(prhs[4]); 
	MT_MGP_CPU.ly = mxGetScalar(prhs[5]);
	MT_MGP_CPU.dz = mxGetScalar(prhs[6]);
	iConfFP = (int)mxGetScalar(prhs[7]); 
	MT_MGP_CPU.DimFP = (int)mxGetScalar(prhs[8]); 
	MT_MGP_CPU.SeedFP = (int)mxGetScalar(prhs[9]);
	iSlice = (int)mxGetScalar(prhs[10])-1; 
	if(iSlice<0) iSlice = 0;

	/************************Output data**************************/
	plhs[0] = mxCreateDoubleMatrix(MT_MGP_CPU.ny, MT_MGP_CPU.nx, mxREAL);
	double *V0h = mxGetPr(plhs[0]);

	plhs[1] = mxCreateDoubleMatrix(MT_MGP_CPU.ny, MT_MGP_CPU.nx, mxREAL);
	double *V1h = mxGetPr(plhs[1]);

	cudaSetDevice(MT_MGP_CPU.gpu);

	cMT_Potential_GPU MT_Potential_GPU;	
	MT_Potential_GPU.SetInputData(&MT_MGP_CPU, nAtomsM, AtomsM);
	MT_Potential_GPU.MoveAtoms(iConfFP);
	MT_Potential_GPU.ProjectedPotential(iSlice);

	f_fft2Shift_MD(MT_Potential_GPU.GP, MT_Potential_GPU.V0, MT_Potential_GPU.V1);
	cudaMemcpy(V0h, MT_Potential_GPU.V0, MT_Potential_GPU.GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
	cudaMemcpy(V1h, MT_Potential_GPU.V1, MT_Potential_GPU.GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);

	//MT_Potential_GPU.freeMemory();
	//cudaDeviceReset();
}