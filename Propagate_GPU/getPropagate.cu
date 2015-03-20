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
#include "cuda.h"
#include "cuda_runtime.h"
#include <device_functions.h>
#include "cufft.h"
#include <mex.h>

class cPropagator{
	private:
		sGP GP;
		sACD Prop_x;
		sACD Prop_y;
		sComplex Psii;
		double2 *Psi;
		cufftHandle PlanPsi;
		double lambda;
	public:
		void freeMemory();
		cPropagator();
		~cPropagator();
		void SetInputData(cMT_MGP_CPU *MT_MGP_CPU_i, sComplex &Psiih);
		void Propagate(double dz, sComplex &Psih);
};

void cPropagator::freeMemory()
{
	cudaDeviceSynchronize(); // wait to finish the work in the GPU
	f_sGP_Init(GP);
	f_sACD_Free_GPU(Prop_x);
	f_sACD_Free_GPU(Prop_y);
	lambda = 0;
	cudaFreen(Psii.real);
	cudaFreen(Psii.imag);
	cudaFreen(Psi);
	cufftDestroyn(PlanPsi);
}

cPropagator::~cPropagator()
{
	freeMemory();
}

cPropagator::cPropagator()
{
	f_sGP_Init(GP);
	f_sACD_Init_GPU(Prop_x);
	f_sACD_Init_GPU(Prop_y);
	lambda = 0;
	Psii.real = 0;
	Psii.imag = 0;
	Psi = 0;
	PlanPsi = 0;
}

void cPropagator::SetInputData(cMT_MGP_CPU *MT_MGP_CPU_i, sComplex &Psiih)
{
	freeMemory();

	cudaSetDevice(MT_MGP_CPU_i->GPU_Device);
	lambda = f_getLambda(MT_MGP_CPU_i->E0);
	f_sGP_SetInputData(MT_MGP_CPU_i, GP);

	f_sACD_Malloc_GPU(GP.nx, Prop_x);
	f_sACD_Malloc_GPU(GP.ny, Prop_y);
	cudaMalloc((void**)&(Psii.real), GP.nxy*cSizeofRD);
	cudaMalloc((void**)&(Psii.imag), GP.nxy*cSizeofRD);
	cudaMalloc((void**)&Psi, GP.nxy*cSizeofCD);
	cufftPlan2d(&PlanPsi, GP.nx, GP.ny, CUFFT_Z2Z);

	cudaMemcpy(Psii.real, Psiih.real, GP.nxy*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(Psii.imag, Psiih.imag, GP.nxy*cSizeofRD, cudaMemcpyHostToDevice);
	f_fft2Shift_MC_GPU(GP, Psii);
}

void cPropagator::Propagate(double dz, sComplex &Psih)
{
	double gxu = 0, gyu = 0;
	f_Set_MC_GPU(GP, Psii, Psi);
	f_Propagate_GPU(PlanPsi, GP, eSReal, gxu, gyu, lambda, dz, Prop_x, Prop_y, Psi);
	// fft2shift 
	f_fft2Shift_MC_GPU(GP, Psi);
	/*********************copy data to host************************/
	f_Copy_MCd(GP, Psi, Psii.real, Psii.imag, Psih);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	cMT_MGP_CPU MT_MGP_CPU;
	sComplex Psiih, Psioh;
	cPropagator Propagator;

	Psiih.real = mxGetPr(prhs[0]);
	Psiih.imag = mxGetPi(prhs[0]);
	MT_MGP_CPU.GPU_Device = (int)mxGetScalar(prhs[1]); 
	MT_MGP_CPU.E0 = mxGetScalar(prhs[2]); 
	MT_MGP_CPU.nx = (int)mxGetScalar(prhs[3]); 
	MT_MGP_CPU.ny = (int)mxGetScalar(prhs[4]);
	MT_MGP_CPU.lx = mxGetScalar(prhs[5]); 
	MT_MGP_CPU.ly = mxGetScalar(prhs[6]);
	MT_MGP_CPU.dz = mxGetScalar(prhs[7]);
	MT_MGP_CPU.BWL = (int)mxGetScalar(prhs[8]);
	/************************Output data**************************/
	plhs[0] = mxCreateDoubleMatrix(MT_MGP_CPU.ny, MT_MGP_CPU.nx, mxCOMPLEX);
	Psioh.real = mxGetPr(plhs[0]);
	Psioh.imag = mxGetPi(plhs[0]);

	Propagator.SetInputData(&MT_MGP_CPU, Psiih);
	Propagator.Propagate(MT_MGP_CPU.dz, Psioh);

	//Propagator.freeMemory();
	//cudaDeviceReset();
}