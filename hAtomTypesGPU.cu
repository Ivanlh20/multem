#include <memory.h>
#include "hConstTypes.h"
#include "hAtomTypesGPU.h"
#include "math.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include <device_functions.h>

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
	PotPar = 0;
	occ = 0;

	for (int i=0; i<6; i++){
		fegh.cl[i] = fegh.cnl[i] = 0;
		fxgh.cl[i] = fxgh.cnl[i] = 0;
		Prh.cl[i] = Prh.cnl[i] = 0;
		Vrh.cl[i] = Vrh.cnl[i] = 0;
		VRh.cl[i] = VRh.cnl[i] = 0;
	}

	cudaFreen(feg.cl);
	cudaFreen(feg.cnl);
	cudaFreen(fxg.cl);
	cudaFreen(fxg.cnl);
	cudaFreen(Pr.cl);
	cudaFreen(Pr.cnl);
	cudaFreen(Vr.cl);
	cudaFreen(Vr.cnl);
	cudaFreen(VR.cl);
	cudaFreen(VR.cnl);

	if(ns>0){		
		for (int i=0; i<ns; i++){
			Vo[i].sigma = 0;
			cudaFreen(Vo[i].Vr.cl);
			cudaFreen(Vo[i].Vr.cnl);
			cudaFreen(Vo[i].Vi.cl);
			cudaFreen(Vo[i].Vi.cnl);
		}
		delete [] Vo; Vo = 0;
		ns = 0;
	}
	
	nR = 0;
	cudaFreen(R);
	cudaFreen(R2);
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
	PotPar = 0;
	occ = 0;

	for (int i=0; i<6; i++){
		fegh.cl[i] = fegh.cnl[i] = 0;
		fxgh.cl[i] = fxgh.cnl[i] = 0;
		Prh.cl[i] = Prh.cnl[i] = 0;
		Vrh.cl[i] = Vrh.cnl[i] = 0;
		VRh.cl[i] = VRh.cnl[i] = 0;
	}

	feg.cl = 0;
	feg.cnl = 0;
	fxg.cl = 0;
	fxg.cnl = 0;
	Pr.cl = 0;
	Pr.cnl = 0;
	Vr.cl = 0;
	Vr.cnl = 0;
	VR.cl = 0;
	VR.cnl = 0;

	ns = 0;
	Vo = 0;

	nR = 0;
	R = 0;
	R2 = 0;
}

// Destructor
cAtomTypesGPU::~cAtomTypesGPU(){
	freeMemory(); // clean GPU memory
}

// Set Atom type
void cAtomTypesGPU::SetAtomTypes(sAtomTypesCPU AtomTypesCPUi, int nRi, double dRmini){
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
	PotPar = AtomTypesCPUi.PotPar;
	occ = AtomTypesCPUi.occ;

	/***************************************************************/
	memcpy(fegh.cl, AtomTypesCPUi.feg.cl, 6*cSizeofRD);
	memcpy(fegh.cnl, AtomTypesCPUi.feg.cnl, 6*cSizeofRD);

	memcpy(fxgh.cl, AtomTypesCPUi.fxg.cl, 6*cSizeofRD);
	memcpy(fxgh.cnl, AtomTypesCPUi.fxg.cnl, 6*cSizeofRD);

	memcpy(Prh.cl, AtomTypesCPUi.Pr.cl, 6*cSizeofRD);
	memcpy(Prh.cnl, AtomTypesCPUi.Pr.cnl, 6*cSizeofRD);

	for (int i=0; i<6; i++)
		Vrh.cl[i] = AtomTypesCPUi.Vr.cl[i]/cPotf;
	memcpy(Vrh.cnl, AtomTypesCPUi.Vr.cnl, 6*cSizeofRD);

	for (int i=0; i<6; i++)
		VRh.cl[i] = AtomTypesCPUi.VR.cl[i]/cPotf;
	memcpy(VRh.cnl, AtomTypesCPUi.VR.cnl, 6*cSizeofRD);
	/***************************************************************/

	cudaMalloc((void**)&(feg.cl), 6*cSizeofRD);
	cudaMemcpy(feg.cl, fegh.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&(feg.cnl), 6*cSizeofRD);
	cudaMemcpy(feg.cnl, fegh.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	cudaMalloc((void**)&(fxg.cl), 6*cSizeofRD);
	cudaMemcpy(fxg.cl, fxgh.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&(fxg.cnl), 6*cSizeofRD);
	cudaMemcpy(fxg.cnl, fxgh.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	cudaMalloc((void**)&(Pr.cl), 6*cSizeofRD);
	cudaMemcpy(Pr.cl, Prh.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&(Pr.cnl), 6*cSizeofRD);
	cudaMemcpy(Pr.cnl, Prh.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	cudaMalloc((void**)&(Vr.cl), 6*cSizeofRD);
	cudaMemcpy(Vr.cl, Vrh.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&(Vr.cnl), 6*cSizeofRD);
	cudaMemcpy(Vr.cnl, Vrh.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	cudaMalloc((void**)&(VR.cl), 6*cSizeofRD);
	cudaMemcpy(VR.cl, VRh.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&(VR.cnl), 6*cSizeofRD);
	cudaMemcpy(VR.cnl, VRh.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	/***********************************************************************************/
	ns = AtomTypesCPUi.ns;
	if(ns>0){
		Vo = new sVoGPU[ns];
		for (int i=0; i<ns; i++){
			Vo[i].sigma = AtomTypesCPUi.Vo[i].sigma;
			cudaMalloc((void**)&(Vo[i].Vr.cl), 6*cSizeofRD);
			cudaMemcpy(Vo[i].Vr.cl, AtomTypesCPUi.Vo[i].Vr.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
			cudaMalloc((void**)&(Vo[i].Vr.cnl), 6*cSizeofRD);
			cudaMemcpy(Vo[i].Vr.cnl, AtomTypesCPUi.Vo[i].Vr.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

			cudaMalloc((void**)&(Vo[i].Vi.cl), 6*cSizeofRD);
			cudaMemcpy(Vo[i].Vi.cl, AtomTypesCPUi.Vo[i].Vi.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
			cudaMalloc((void**)&(Vo[i].Vi.cnl), 6*cSizeofRD);
			cudaMemcpy(Vo[i].Vi.cnl, AtomTypesCPUi.Vo[i].Vi.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);
		}
	}

	nR = nRi;
	double dlnr = log(Rmax/Rmin)/double(nR-1);
	double *Rh, *R2h;
	Rh = new double [nR];
	R2h = new double [nR];

	for (int i=0; i<nR; i++){
		Rh[i] = Rmin*exp(double(i)*dlnr);
		R2h[i] = Rh[i]*Rh[i];
	}

	cudaMalloc((void**)&R, nR*cSizeofRD);
	cudaMemcpy(R, Rh, nR*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&R2, nR*cSizeofRD);
	cudaMemcpy(R2, R2h, nR*cSizeofRD, cudaMemcpyHostToDevice);

	delete [] Rh; Rh = 0;
	delete [] R2h; R2h = 0;
}