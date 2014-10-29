#include <memory.h>
#include "hConstTypes.h"
#include "hgeneralCPU.h"
#include "hgeneralGPU.h"
#include "hPotentialCPU.h"
#include "hAtomTypesGPU.h"
#include "math.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>

void CubicPolyCoefHost(int PotPari, sAtomTypesCPU &AtomTypesCPUi, int nR, double *&R, double *&R2, sciVn &ciVR){
	int IntType = 0, Dim = 2;
	double sigma = 0.0, *VR, *dVR;
	double dR2, V, Vn, dV, dVn, m, n;

	VR = new double[nR];
	dVR = new double[nR];
	cPotentialCPU PotentialCPU;
	PotentialCPU.SetAtomTypes(PotPari, AtomTypesCPUi);
	PotentialCPU.SetSigma(sigma);
	PotentialCPU.Vr(IntType, Dim, nR, R, VR, dVR);
	for(int iR=0; iR<nR; iR++){
		ciVR.c0[iR] = VR[iR]/cPotf;
		ciVR.c1[iR] = 0.5*dVR[iR]/(cPotf*R[iR]);
	}

	for(int iR=0; iR<nR-1; iR++){
		dR2 = 1.0/(R2[iR+1]-R2[iR]);
		V = ciVR.c0[iR]; Vn = ciVR.c0[iR+1];
		dV = ciVR.c1[iR]; dVn = ciVR.c1[iR+1];
		m = (Vn-V)*dR2; n = dV+dVn;
		ciVR.c0[iR] = V-ciVR.c0[nR-1];
		ciVR.c2[iR] = (3.0*m-n-dV)*dR2;
		ciVR.c3[iR] = (n-2.0*m)*dR2*dR2;
	}

	delete [] VR; VR = 0;
	delete [] dVR; dVR = 0;
}

// free memory
void cAtomTypesGPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	Z = 0;
	m = 0;
	A = 0;
	rn_e = 0;
	rn_c = 0;
	ra_e = 0;
	ra_c = 0;
	Rmin = 0;
	Rmax = 0;
	Rmin2 = 0;
	Rmax2 = 0;

	f_sCoefPar_Free(cfegh);
	f_sCoefPar_Free(cfxgh);
	f_sCoefPar_Free(cPrh);
	f_sCoefPar_Free(cVrh);
	f_sCoefPar_Free(cVRh);

	f_sCoefPar_cudaFree(cfeg);
	f_sCoefPar_cudaFree(cfxg);
	f_sCoefPar_cudaFree(cPr);
	f_sCoefPar_cudaFree(cVr);
	f_sCoefPar_cudaFree(cVR);

	if(ns>0){		
		for (int i=0; i<ns; i++){
			Vo[i].sigma = 0;\
			f_sCoefPar_cudaFree(Vo[i].cVr);
			f_sCoefPar_cudaFree(Vo[i].cVi);
		}
		delete [] Vo; Vo = 0;
		ns = 0;
	}
	
	nR = 0;

	delete [] Rh; Rh = 0;
	delete [] R2h; R2h = 0;
	f_sciVn_Free(ciVRh);

	cudaFreen(R);
	cudaFreen(R2);
	f_sciVn_cudaFree(ciVR);
};

// Constructor
cAtomTypesGPU::cAtomTypesGPU(){
	Z = 0;
	m = 0;
	A = 0;
	rn_e = 0;
	rn_c = 0;
	ra_e = 0;
	ra_c = 0;
	Rmin = 0;
	Rmax = 0;
	Rmin2 = 0;
	Rmax2 = 0;

	f_sCoefPar_Init(cfegh);
	f_sCoefPar_Init(cfxgh);
	f_sCoefPar_Init(cPrh);
	f_sCoefPar_Init(cVrh);
	f_sCoefPar_Init(cVRh);

	f_sCoefPar_cudaInit(cfeg);
	f_sCoefPar_cudaInit(cfxg);
	f_sCoefPar_cudaInit(cPr);
	f_sCoefPar_cudaInit(cVr);
	f_sCoefPar_cudaInit(cVR);

	ns = 0;
	Vo = 0;

	nR = 0;

	Rh = 0;
	R2h = 0;
	f_sciVn_Init(ciVRh);

	R = 0;
	R2 = 0;
	f_sciVn_cudaInit(ciVR);
}

// Destructor
cAtomTypesGPU::~cAtomTypesGPU(){
	freeMemory(); // clean GPU memory
}

// Set Atom type
void cAtomTypesGPU::SetAtomTypes(int PotPari, sAtomTypesCPU &AtomTypesCPUi, int nRi, double dRmini){
	freeMemory(); // clean GPU memory

	Z = AtomTypesCPUi.Z;
	m = AtomTypesCPUi.m;
	A = AtomTypesCPUi.A;
	rn_e = AtomTypesCPUi.rn_e;
	rn_c = AtomTypesCPUi.rn_c;
	ra_e = AtomTypesCPUi.ra_e;
	ra_c = AtomTypesCPUi.ra_c;
	Rmin = MAX(rn_c, dRmini);
	Rmax = AtomTypesCPUi.Rmax;
	Rmin2 = Rmin*Rmin;
	Rmax2 = Rmax*Rmax;

	/***************************************************************/
	f_sCoefPar_Malloc(6, cfegh);
	memcpy(cfegh.cl, AtomTypesCPUi.cfeg.cl, 6*cSizeofRD);
	memcpy(cfegh.cnl, AtomTypesCPUi.cfeg.cnl, 6*cSizeofRD);

	f_sCoefPar_Malloc(6, cfxgh);
	memcpy(cfxgh.cl, AtomTypesCPUi.cfxg.cl, 6*cSizeofRD);
	memcpy(cfxgh.cnl, AtomTypesCPUi.cfxg.cnl, 6*cSizeofRD);

	f_sCoefPar_Malloc(6, cPrh);
	memcpy(cPrh.cl, AtomTypesCPUi.cPr.cl, 6*cSizeofRD);
	memcpy(cPrh.cnl, AtomTypesCPUi.cPr.cnl, 6*cSizeofRD);

	f_sCoefPar_Malloc(6, cVrh);
	for (int i=0; i<6; i++)
		cVrh.cl[i] = AtomTypesCPUi.cVr.cl[i]/cPotf;
	memcpy(cVrh.cnl, AtomTypesCPUi.cVr.cnl, 6*cSizeofRD);

	f_sCoefPar_Malloc(6, cVRh);
	for (int i=0; i<6; i++)
		cVRh.cl[i] = AtomTypesCPUi.cVR.cl[i]/cPotf;
	memcpy(cVRh.cnl, AtomTypesCPUi.cVR.cnl, 6*cSizeofRD);

	/***************************************************************/
	f_sCoefPar_cudaMalloc(6, cfeg);
	cudaMemcpy(cfeg.cl, cfegh.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(cfeg.cnl, cfegh.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	f_sCoefPar_cudaMalloc(6, cfxg);
	cudaMemcpy(cfxg.cl, cfxgh.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(cfxg.cnl, cfxgh.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	f_sCoefPar_cudaMalloc(6, cPr);
	cudaMemcpy(cPr.cl, cPrh.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(cPr.cnl, cPrh.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	f_sCoefPar_cudaMalloc(6, cVr);
	cudaMemcpy(cVr.cl, cVrh.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(cVr.cnl, cVrh.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	f_sCoefPar_cudaMalloc(6, cVR);
	cudaMemcpy(cVR.cl, cVRh.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(cVR.cnl, cVRh.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	/***********************************************************************************/
	ns = AtomTypesCPUi.ns;
	if(ns>0){
		Vo = new sVoGPU[ns];
		for (int i=0; i<ns; i++){
			Vo[i].sigma = AtomTypesCPUi.Vo[i].sigma;
			f_sCoefPar_cudaMalloc(6, Vo[i].cVr);
			cudaMemcpy(Vo[i].cVr.cl, AtomTypesCPUi.Vo[i].cVr.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
			cudaMemcpy(Vo[i].cVr.cnl, AtomTypesCPUi.Vo[i].cVr.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

			f_sCoefPar_cudaMalloc(6, Vo[i].cVi);
			cudaMemcpy(Vo[i].cVi.cl, AtomTypesCPUi.Vo[i].cVi.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
			cudaMemcpy(Vo[i].cVi.cnl, AtomTypesCPUi.Vo[i].cVi.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);
		}
	}

	nR = nRi;

	double dlnr = log(Rmax/Rmin)/double(nR-1);
	Rh = new double [nR];
	R2h = new double [nR];
	for (int i=0; i<nR; i++){
		Rh[i] = Rmin*exp(double(i)*dlnr);
		R2h[i] = Rh[i]*Rh[i];
	}
	f_sciVn_Malloc(nR, ciVRh);
	CubicPolyCoefHost(PotPari, AtomTypesCPUi, nR, Rh, R2h, ciVRh);

	cudaMalloc((void**)&R, nR*cSizeofRD);
	cudaMemcpy(R, Rh, nR*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&R2, nR*cSizeofRD);
	cudaMemcpy(R2, R2h, nR*cSizeofRD, cudaMemcpyHostToDevice);

	f_sciVn_cudaMalloc(stnR, ciVR);
	cudaMemcpy(ciVR.c0, ciVRh.c0, nR*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(ciVR.c1, ciVRh.c1, nR*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(ciVR.c2, ciVRh.c2, nR*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(ciVR.c3, ciVRh.c3, nR*cSizeofRD, cudaMemcpyHostToDevice);
}