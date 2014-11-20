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

#include "hConstTypes.h"
#include "hMT_General_CPU.h"
#include "hMT_General_GPU.h"
#include "hMT_IncidentWave_GPU.h"
#include "hMatlab2Cpp.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include <device_functions.h>
#include "cufft.h"
#include <mex.h>

class cProbe{
	private:
		double x0;
		double y0;
		sMGP MGP;
		sGP GP;
		sLens Lens;
		sComplex Psii;
		double2 *Psi;
		cufftHandle PlanPsi;
		cMT_IncidentWave_GPU MT_IncidentWave_GPU;
	public:
		void freeMemory();
		cProbe();
		~cProbe();
		void SetInputData(sInProbe &InProbe);
		void getProbe(sComplex &Psih);
};

void cProbe::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	f_sMPG_Init(MGP);
	f_sGP_Init(GP);
	f_sLens_Init(Lens);

	cudaFreen(Psii.real);
	cudaFreen(Psii.imag);
	cudaFreen(Psi);
	cufftDestroyn(PlanPsi);
}

cProbe::~cProbe(){
	freeMemory();
}

cProbe::cProbe(){
	f_sMPG_Init(MGP);
	f_sGP_Init(GP);
	f_sLens_Init(Lens);
	Psii.real = 0;
	Psii.imag = 0;
	Psi = 0;
	PlanPsi = 0;
}

void cProbe::SetInputData(sInProbe &InProbe){
	freeMemory();

	MGP.gpu = InProbe.gpu;
	MGP.E0 = InProbe.E0;	
	MGP.theta = InProbe.theta;	
	MGP.phi = InProbe.phi;
	MGP.lx = InProbe.lx;
	MGP.ly = InProbe.ly;
	MGP.nx = InProbe.nx;
	MGP.ny = InProbe.ny;
	MGP.BWL = true;

	cudaSetDevice(MGP.gpu);

	f_sGP_Cal(MGP.nx, MGP.ny, MGP.lx, MGP.ly, MGP.dz, MGP.PBC_xy, MGP.BWL, GP);

	Lens.m = InProbe.m;
	Lens.f = InProbe.f;
	Lens.Cs3 = InProbe.Cs3;
	Lens.Cs5 = InProbe.Cs5;
	Lens.mfa2 = InProbe.mfa2;
	Lens.afa2 = InProbe.afa2;
	Lens.mfa3 = InProbe.mfa3;
	Lens.afa3 = InProbe.afa3;
	Lens.aobjl = InProbe.aobjl;
	Lens.aobju = InProbe.aobju;
	Lens.sf = InProbe.sf;
	Lens.nsf = InProbe.nsf;
	Lens.beta = InProbe.beta;
	Lens.nbeta = InProbe.nbeta;
	f_sLens_Cal(MGP.E0, GP, Lens);

	double gxu = sin(MGP.theta)*cos(MGP.phi)/Lens.lambda;
	double gyu = sin(MGP.theta)*sin(MGP.phi)/Lens.lambda;;

	x0 = InProbe.x0;
	y0 = InProbe.y0;

	cudaMalloc((void**)&(Psii.real), GP.nxy*cSizeofRD);
	cudaMalloc((void**)&(Psii.imag), GP.nxy*cSizeofRD);
	cudaMalloc((void**)&Psi, GP.nxy*cSizeofCD);
	cufftPlan2d(&PlanPsi, GP.nx, GP.ny, CUFFT_Z2Z);

	MT_IncidentWave_GPU.SetInputData(GP, Lens, PlanPsi);
}

void cProbe::getProbe(sComplex &Psih){
	MT_IncidentWave_GPU.Psi0(x0, y0, Psi);
	// fft2shift 
	f_fft2Shift_MC(GP, Psi);
	/**********************copy data to host************************/
	f_Copy_MCd(GP, Psi, Psii.real, Psii.imag, Psih);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	sInProbe InProbe;
	sComplex Psih;
	cProbe Probe;

	Matlab2InProbe(prhs[0], InProbe);
	/*************************Output data**************************/
	plhs[0] = mxCreateDoubleMatrix(InProbe.ny, InProbe.nx, mxCOMPLEX);
	Psih.real = mxGetPr(plhs[0]);
	Psih.imag = mxGetPi(plhs[0]);

	Probe.SetInputData(InProbe);
	Probe.getProbe(Psih);
}