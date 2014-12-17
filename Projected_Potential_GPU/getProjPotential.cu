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
#include "hMatlab2Cpp.h"

#include "cuda.h"
#include "cuda_runtime.h"
#include <device_functions.h>
#include "cufft.h"
#include <mex.h>

/**********************read input Probe*************************/
void f_Matlab2InProjPotential(const mxArray *mxInProjPotential, sInProjPotential &InProjPotential){
	InProjPotential.gpu = ReadValuemxField<int>(mxInProjPotential, 0, "gpu");						// gpu device
	InProjPotential.MulOrder = 2;																	// 1: First order, 2: Second order
	InProjPotential.iConfFP = ReadValuemxField<int>(mxInProjPotential, 0, "iConfFP");				// Frozen phonon configuration
	if(InProjPotential.iConfFP<0) InProjPotential.iConfFP = 0;
	InProjPotential.DimFP = ReadValuemxField<int>(mxInProjPotential, 0, "DimFP");					// Dimensions phonon configurations
	InProjPotential.DistFP = 1;																		// 1: Gaussian (Phonon distribution)
	InProjPotential.SeedFP = ReadValuemxField<int>(mxInProjPotential, 0, "SeedFP");					// Random seed(frozen phonon)
	InProjPotential.PotPar = ReadValuemxField<int>(mxInProjPotential, 0, "PotPar");					// Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
	InProjPotential.PBC_xy = 1;																		// 1: true, 2: false (Peridic boundary contions)

	InProjPotential.nx = ReadValuemxField<int>(mxInProjPotential, 0, "nx");							// Number of pixels in x direction
	InProjPotential.ny = ReadValuemxField<int>(mxInProjPotential, 0, "ny");							// Number of pixels in y direction
	InProjPotential.lx = ReadValuemxField<double>(mxInProjPotential, 0, "lx");						// distance in x direction(Angstroms)
	InProjPotential.ly = ReadValuemxField<double>(mxInProjPotential, 0, "ly");						// distance in y direction(Angstroms)
	InProjPotential.dz = ReadValuemxField<double>(mxInProjPotential, 0, "dz");						// Slice thickness

	mxArray *mxAtomsM = mxGetField(mxInProjPotential, 0, "Atoms");
	InProjPotential.nAtomsM = (int)mxGetM(mxAtomsM);												// Number of Atoms
	InProjPotential.AtomsM = mxGetPr(mxAtomsM);														// Atoms in a matrix form
	InProjPotential.iSlice = ReadValuemxField<int>(mxInProjPotential, 0, "iSlice")-1;				// Slice
	if(InProjPotential.iSlice<0) InProjPotential.iSlice=0;
 }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	sInProjPotential InProjPotential;

	f_Matlab2InProjPotential(prhs[0], InProjPotential);

	cMT_MGP_CPU MT_MGP_CPU;
	MT_MGP_CPU.SetInputData(InProjPotential);

	/************************Output data**************************/
	plhs[0] = mxCreateDoubleMatrix(MT_MGP_CPU.ny, MT_MGP_CPU.nx, mxREAL);
	double *V0h = mxGetPr(plhs[0]);

	plhs[1] = mxCreateDoubleMatrix(MT_MGP_CPU.ny, MT_MGP_CPU.nx, mxREAL);
	double *V1h = mxGetPr(plhs[1]);

	cudaSetDevice(MT_MGP_CPU.gpu);
	cMT_Potential_GPU MT_Potential_GPU;	
	MT_Potential_GPU.SetInputData(&MT_MGP_CPU, InProjPotential.nAtomsM, InProjPotential.AtomsM);
	MT_Potential_GPU.MoveAtoms(InProjPotential.iConfFP);
	MT_Potential_GPU.ProjectedPotential(InProjPotential.iSlice, 2);

	f_fft2Shift_MD(MT_Potential_GPU.GP, MT_Potential_GPU.V0, MT_Potential_GPU.V1);
	cudaMemcpy(V0h, MT_Potential_GPU.V0, MT_Potential_GPU.GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
	cudaMemcpy(V1h, MT_Potential_GPU.V1, MT_Potential_GPU.GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);

	MT_Potential_GPU.freeMemoryReset();
}