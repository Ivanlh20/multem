#include "hConstTypes.h"
#include "hTEMIm.h"
#include "hgeneralCPU.h"
#include "hgeneralGPU.h"
#include "hMicroscopeEffectsGPU.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

cTEMIm::cTEMIm(){
	gpu = 0;
	MEffect = 0;
	Psirh = 0;
	Psiih = 0;
	E0 = 0;
	lx = 0;
	ly = 0;
	nx = 0;
	ny = 0;
	lambda = 0;

	fPsi = 0;
	Psia = 0;
	M2Psis = 0;
	M2Psit = 0;
	PlanPsi = 0;
}

void cTEMIm::freeMemory(){
	gpu = 0;
	MEffect = 0;
	Psirh = 0;
	Psiih = 0;
	E0 = 0;
	lx = 0;
	ly = 0;
	nx = 0;
	ny = 0;
	lambda = 0;

	cudaFreen(fPsi);
	cudaFreen(Psia);
	cudaFreen(M2Psis);
	cudaFreen(M2Psit);
	cufftDestroyn(PlanPsi);

	MicroscopeEffectsGPU.freeMemory();
}

void cTEMIm::GenerateParameters(){
	GP.nx = nx;
	GP.ny = ny;
	GP.nxh = GP.nx/2;
	GP.nyh = GP.ny/2;
	GP.nxy = GP.nx*GP.ny;
	GP.inxy = (GP.nxy==0)?(0.0):(1.0/GP.nxy);
	GP.dRx = (nx==0)?(0.0):(lx/GP.nx);
	GP.dRy = (ny==0)?(0.0):(ly/GP.ny);
	GP.dgx = (lx==0)?(0.0):(1.0/lx);
	GP.dgy = (ly==0)?(0.0):(1.0/ly);

	GP.gmax = MIN(GP.nxh*GP.dgx, GP.nyh*GP.dgy);
	GP.gmax2 = GP.gmax*GP.gmax;
	GP.gmaxl = 2.0*GP.gmax/3.0;
	GP.gmaxl2 = GP.gmaxl*GP.gmaxl;

	lambda = f_getLambda(E0);

	f_sLens_Cal(lambda, GP, Lens);	// Lens coefficients
	f_sBT_Cal(GP, BT);				// Blocks and threads

	cudaFreen(fPsi);
	cudaMalloc((void**)&fPsi, GP.nxy*cSizeofCD);
	cudaFreen(Psia);
	cudaMalloc((void**)&Psia, GP.nxy*cSizeofCD);

	cudaFreen(M2Psis);
	cudaMalloc((void**)&M2Psis, GP.nxy*cSizeofRD);
	cudaFreen(M2Psit);
	cudaMalloc((void**)&M2Psit, GP.nxy*cSizeofRD);

	cufftDestroyn(PlanPsi); 
	cufftPlan2d(&PlanPsi, nx, ny, CUFFT_Z2Z);

	// Microscope parameters
	//MicroscopeEffectsGPU.SetInputData(BT, GP, Lens, PlanPsi, fPsi, Psia, M2Psit);
}

void cTEMIm::SetInputData(sInTEMIm &InTEMIm){
	gpu = InTEMIm.gpu; 
	MEffect = InTEMIm.MEffect;
	Psirh= InTEMIm.Psirh;
	Psiih = InTEMIm.Psiih;
	E0 = InTEMIm.E0;
	lx = InTEMIm.lx;
	ly = InTEMIm.ly;
	nx = InTEMIm.nx;
	ny = InTEMIm.ny;
	Lens.m = InTEMIm.MC_m;
	Lens.f = InTEMIm.MC_f;
	Lens.Cs3 = InTEMIm.MC_Cs3;
	Lens.Cs5 = InTEMIm.MC_Cs5;
	Lens.mfa2 = InTEMIm.MC_mfa2;
	Lens.afa2 = InTEMIm.MC_afa2;
	Lens.mfa3 = InTEMIm.MC_mfa3;
	Lens.afa3 = InTEMIm.MC_afa3;
	Lens.aobjl = InTEMIm.MC_aobjl;
	Lens.aobju = InTEMIm.MC_aobju;
	Lens.sf = InTEMIm.MC_sf;
	Lens.nsf = InTEMIm.MC_nsf;
	Lens.beta = InTEMIm.MC_beta;
	Lens.nbeta = InTEMIm.MC_nbeta;

	cudaSetDevice(gpu);
	cudaDeviceReset();

	GenerateParameters();

 }

// Partially coherent transfer function and Transmission cross coefficient
void cTEMIm::TEMImage(double *Psir_hi, double *Psii_hi, double *M2Psi_ho){
	sComplex Psi;
	Psi.real = M2Psis, Psi.imag = M2Psit;
	// Copy real part of Psi
	cudaMemcpy(Psi.real, Psir_hi, GP.nxy*cSizeofRD, cudaMemcpyHostToDevice);
	// Copy imaginary part of Psi
	cudaMemcpy(Psi.imag, Psii_hi, GP.nxy*cSizeofRD, cudaMemcpyHostToDevice);
	// Set real and imaginary part to Psi
	f_SetRealImagVectorC(BT.Bnxy, BT.Tnxy, GP.nxy, Psi, fPsi);
	// fft2shift
	f_fft2ShiftC(BT.Bhnxny, BT.Thnxny, GP.nxh, GP.nyh, fPsi);
	// Forward fft2
	cufftExecZ2Z(PlanPsi, fPsi, fPsi, CUFFT_FORWARD);

	/**********************Microscope effects**********************/
	//if (MEffect==0)
	//	MicroscopeEffectsGPU.PCLIMWPOTEM(M2Psis);	
	//else
	//	MicroscopeEffectsGPU.PCTCCTEM(M2Psis);
	/*************************************************************/
	// fft2shift
	f_fft2ShiftD(BT.Bhnxny, BT.Thnxny, GP.nxh, GP.nyh, M2Psis);
	// copy M2Psi to the host
	cudaMemcpy(M2Psi_ho, M2Psis, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
}