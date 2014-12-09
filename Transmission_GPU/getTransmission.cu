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
#include "hMT_General_GPU.h"
#include "hMT_Transmission_GPU.h"
#include "hMatlab2Cpp.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>
#include <mex.h>

/**********************read input Probe*************************/
void f_Matlab2InTransmission(const mxArray *mxInTransmission, sInTransmission &InTransmission){
	InTransmission.gpu = ReadValuemxField<int>(mxInTransmission, 0, "gpu");							// gpu device
	InTransmission.MulOrder = 2;																	// 1: First order, 2: Second order
	InTransmission.iConfFP = ReadValuemxField<int>(mxInTransmission, 0, "iConfFP");					// Frozen phonon configuration
	if(InTransmission.iConfFP<0) InTransmission.iConfFP = 0;
	InTransmission.DimFP = ReadValuemxField<int>(mxInTransmission, 0, "DimFP");						// Dimensions phonon configurations
	InTransmission.DistFP = 1;																		// 1: Gaussian (Phonon distribution)
	InTransmission.SeedFP = ReadValuemxField<int>(mxInTransmission, 0, "SeedFP");					// Random seed(frozen phonon)
	InTransmission.PotPar = ReadValuemxField<int>(mxInTransmission, 0, "PotPar");					// Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
	InTransmission.ApproxModel = ReadValuemxField<int>(mxInTransmission, 0, "ApproxModel");			// 1: Mulstilice, 2: Projection approximation, 3: Phase object approximation, 4: Weak phase object approximation
	if(InTransmission.ApproxModel>2) InTransmission.iConfFP = 0;
	InTransmission.BWL = ReadValuemxField<int>(mxInTransmission, 0, "BWL");							// 1: true, 2: false (bandwidth limited)
	InTransmission.FastCal = ReadValuemxField<int>(mxInTransmission, 0, "FastCal");					// 1: normal mode(low memory consumption), 2: fast calculation(high memory consumption)
	InTransmission.PBC_xy = 1;																		// 1: true, 2: false (Peridic boundary contions)

	/**************************Multislice*************************/
	InTransmission.E0 = ReadValuemxField<double>(mxInTransmission, 0, "E0");						// Acceleration voltage
	InTransmission.theta = ReadValuemxField<double>(mxInTransmission, 0, "theta", deg2rad);			// incident tilt (in spherical coordinates) (degrees-->rad)
	InTransmission.phi = ReadValuemxField<double>(mxInTransmission, 0, "phi", deg2rad);				// incident tilt (in spherical coordinates) (degrees-->rad)
	InTransmission.nx = ReadValuemxField<int>(mxInTransmission, 0, "nx");							// Number of pixels in x direction
	InTransmission.ny = ReadValuemxField<int>(mxInTransmission, 0, "ny");							// Number of pixels in y direction
	InTransmission.lx = ReadValuemxField<double>(mxInTransmission, 0, "lx");						// distance in x direction(Angstroms)
	InTransmission.ly = ReadValuemxField<double>(mxInTransmission, 0, "ly");						// distance in y direction(Angstroms)
	InTransmission.dz = ReadValuemxField<double>(mxInTransmission, 0, "dz");						// Slice thickness

	mxArray *mxAtomsM = mxGetField(mxInTransmission, 0, "Atoms");
	InTransmission.nAtomsM = (int)mxGetM(mxAtomsM);													// Number of Atoms
	InTransmission.AtomsM = mxGetPr(mxAtomsM);														// Atoms in a matrix form
	InTransmission.iSlice = ReadValuemxField<int>(mxInTransmission, 0, "iSlice")-1;					// Slice
	if(InTransmission.iSlice<0) InTransmission.iSlice = 0;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	sInTransmission InTransmission;

	f_Matlab2InTransmission(prhs[0], InTransmission);

	cMT_MGP_CPU MT_MGP_CPU;
	MT_MGP_CPU.SetInputData(InTransmission);

	/************************Output data**************************/
	sComplex Transh;
	plhs[0] = mxCreateDoubleMatrix(MT_MGP_CPU.ny, MT_MGP_CPU.nx, mxCOMPLEX);
	Transh.real = mxGetPr(plhs[0]);
	Transh.imag  = mxGetPi(plhs[0]);

	cudaSetDevice(MT_MGP_CPU.gpu);
	cufftHandle PlanTrans=0;
	cufftPlan2d(&PlanTrans, MT_MGP_CPU.ny, MT_MGP_CPU.nx, CUFFT_Z2Z);

	cMT_Transmission_GPU MT_Transmission_GPU;	
	MT_Transmission_GPU.SetInputData(&MT_MGP_CPU, PlanTrans, InTransmission.nAtomsM, InTransmission.AtomsM);
	MT_Transmission_GPU.MoveAtoms(InTransmission.iConfFP);
	double2 *Trans;
	Trans = MT_Transmission_GPU.getTrans(InTransmission.iSlice, 2);

	f_fft2Shift_MC(MT_Transmission_GPU.GP, Trans);
	/*********************copy data to host************************/
	f_Copy_MCd(MT_Transmission_GPU.GP, Trans, MT_Transmission_GPU.V0, MT_Transmission_GPU.V1, Transh);

	cufftDestroyn(PlanTrans);

	MT_Transmission_GPU.freeMemoryReset();
}