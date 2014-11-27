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
		sBT BT;
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
		void SetInputData(int gpu, double E0, int nx, int ny, double lx, double ly, bool BWL, sComplex &Psiih);
		void Propagate(double dz, sComplex &Psih);
};

void cPropagator::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU
	f_sGP_Init(GP);
	f_sACD_cudaFree(Prop_x);
	f_sACD_cudaFree(Prop_y);
	lambda = 0;
	cudaFreen(Psii.real);
	cudaFreen(Psii.imag);
	cudaFreen(Psi);
	cufftDestroyn(PlanPsi);
}

cPropagator::~cPropagator(){
	freeMemory();
}

cPropagator::cPropagator(){
	f_sGP_Init(GP);
	f_sACD_cudaInit(Prop_x);
	f_sACD_cudaInit(Prop_y);
	lambda = 0;
	Psii.real = 0;
	Psii.imag = 0;
	Psi = 0;
	PlanPsi = 0;
}

void cPropagator::SetInputData(int gpu, double E0, int nx, int ny, double lx, double ly, bool BWL, sComplex &Psiih){
	freeMemory();

	cudaSetDevice(gpu);
	lambda = f_getLambda(E0);
	f_sGP_Cal(nx, ny, lx, ly, 0, true, BWL, GP);

	f_sACD_cudaMalloc(GP.nx, Prop_x);
	f_sACD_cudaMalloc(GP.ny, Prop_y);
	cudaMalloc((void**)&(Psii.real), GP.nxy*cSizeofRD);
	cudaMalloc((void**)&(Psii.imag), GP.nxy*cSizeofRD);
	cudaMalloc((void**)&Psi, GP.nxy*cSizeofCD);
	cufftPlan2d(&PlanPsi, GP.nx, GP.ny, CUFFT_Z2Z);

	cudaMemcpy(Psii.real, Psiih.real, GP.nxy*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(Psii.imag, Psiih.imag, GP.nxy*cSizeofRD, cudaMemcpyHostToDevice);
	f_fft2Shift_MC(GP, Psii);
}

void cPropagator::Propagate(double dz, sComplex &Psih){
	double gxu = 0, gyu = 0;
	f_Set_MC(GP, Psii, Psi);
	f_Propagate(PlanPsi, GP, eSReal, gxu, gyu, lambda, dz, Prop_x, Prop_y, Psi);
	// fft2shift 
	f_fft2Shift_MC(GP, Psi);
	/*********************copy data to host************************/
	f_Copy_MCd(GP, Psi, Psii.real, Psii.imag, Psih);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int nx, ny, gpu;
	double lx, ly, dz, E0;
	bool BWL;
	sComplex Psiih, Psioh;
	cPropagator Propagator;

	Psiih.real = mxGetPr(prhs[0]);
	Psiih.imag = mxGetPi(prhs[0]);
	gpu = (int)mxGetScalar(prhs[1]); 
	E0 = mxGetScalar(prhs[2]); 
	nx = (int)mxGetScalar(prhs[3]); 
	ny = (int)mxGetScalar(prhs[4]);
	lx = mxGetScalar(prhs[5]); 
	ly = mxGetScalar(prhs[6]);
	dz = mxGetScalar(prhs[7]);
	BWL = ((int)mxGetScalar(prhs[8])==1)?true:false;
	/************************Output data**************************/
	plhs[0] = mxCreateDoubleMatrix(ny, nx, mxCOMPLEX);
	Psioh.real = mxGetPr(plhs[0]);
	Psioh.imag = mxGetPi(plhs[0]);

	Propagator.SetInputData(gpu, E0, nx, ny, lx, ly, BWL, Psiih);
	Propagator.Propagate(dz, Psioh);

	//Propagator.freeMemory();
	//cudaDeviceReset();
}