#include "math.h"
#include <memory>
#include "hConstTypes.h"
#include "hgeneralCPU.h"
#include "hgeneralGPU.h"
#include "hDetectorGPU.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>

void cDetectorGPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	f_sGP_Init(GP);
	BS = dim3();

	nDet = 0;
	f_sDetCir_Free(DetCirh);

	delete [] Tot_h; Tot_h = 0;
	delete [] Coh_h; Coh_h = 0;

	cudaFreen(Tot_d);
	cudaFreen(Coh_d);

	cudaFreen(Totp_d);
	cudaFreen(Cohp_d);
}

cDetectorGPU::~cDetectorGPU(){
	freeMemory();
}

cDetectorGPU::cDetectorGPU(){
	f_sGP_Init(GP);
	BS = dim3();

	nDet = 0;
	f_sDetCir_Init(DetCirh);

	Tot_h = 0;
	Coh_h = 0;

	Tot_d = 0;
	Coh_d = 0;

	Totp_d = 0;
	Cohp_d = 0;
}

/******************************************************************************/
/******************************************************************************/

void cDetectorGPU::SetInputData(sGP &GP_i, sBT &BT_i, int nDeti, sDetCir &DetCirhi){
	freeMemory();

	GP = GP_i;
	BS.x = MIN(BT_i.Bnxny.x, 32); BS.y = MIN(BT_i.Bnxny.y, 32); BS.z = 1;
	TS = BT_i.Tnxny;

	nDet = nDeti;
	f_sDetCir_Malloc(nDet, DetCirh);
	memcpy(DetCirh.g2min, DetCirhi.g2min, nDet*cSizeofRD);
	memcpy(DetCirh.g2max, DetCirhi.g2max, nDet*cSizeofRD);

	Tot_h = new double[nDet];
	Coh_h = new double[nDet];

	cudaMalloc((void**)&Tot_d, nDet*cSizeofRD);
	cudaMalloc((void**)&Coh_d, nDet*cSizeofRD);

	cudaMalloc((void**)&Totp_d, BS.x*BS.y*cSizeofRD);
	cudaMalloc((void**)&Cohp_d, BS.x*BS.y*cSizeofRD);
}

void cDetectorGPU::getDetectorIntensity(double *&aM2Psi, double *&M2aPsi, int ixys, sDetInt *DetInth){
	for(int iDet = 0; iDet<nDet; iDet++)
		f_SumM1M2Det(BS, TS, GP, aM2Psi, M2aPsi, DetCirh.g2min[iDet], DetCirh.g2max[iDet], Totp_d, Cohp_d, iDet, Tot_d, Coh_d);

	cudaMemcpy(Tot_h, Tot_d, nDet*cSizeofRD, cudaMemcpyDeviceToHost);
	cudaMemcpy(Coh_h, Coh_d, nDet*cSizeofRD, cudaMemcpyDeviceToHost);

	for(int iDet = 0; iDet<nDet; iDet++){
		DetInth[iDet].Tot[ixys] = Tot_h[iDet];
		DetInth[iDet].Coh[ixys] = Coh_h[iDet];
	}
}

void cDetectorGPU::getDetectorIntensity(double *&aM2Psi, int ixys, sDetInt *DetInth, bool add){
	for(int iDet = 0; iDet<nDet; iDet++)
		f_SumMDet(BS, TS, GP, aM2Psi, DetCirh.g2min[iDet], DetCirh.g2max[iDet], Totp_d, iDet, Tot_d);

	cudaMemcpy(Tot_h, Tot_d, nDet*cSizeofRD, cudaMemcpyDeviceToHost);

	for(int iDet = 0; iDet<nDet; iDet++){
		if(add)
			DetInth[iDet].Tot[ixys] += Tot_h[iDet];
		else
			DetInth[iDet].Tot[ixys] = Tot_h[iDet];
	}
}
