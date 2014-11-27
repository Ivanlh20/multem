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

#include "math.h"
#include <cstring>
#include "hConstTypes.h"
#include "hMT_General_CPU.h"
#include "hMT_General_GPU.h"
#include "hMT_Detector_GPU.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>

void cMT_Detector_GPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	f_sGP_Init(GP);

	nDet = 0;
	f_sDetCir_Free(DetCirh);

	delete [] Tot_h; Tot_h = 0;
	delete [] Coh_h; Coh_h = 0;

	cudaFreen(Tot_d);
	cudaFreen(Coh_d);

	cudaFreen(M1p_d);
	cudaFreen(M2p_d);
}

cMT_Detector_GPU::cMT_Detector_GPU(){
	f_sGP_Init(GP);

	nDet = 0;
	f_sDetCir_Init(DetCirh);

	Tot_h = 0;
	Coh_h = 0;

	Tot_d = 0;
	Coh_d = 0;

	M1p_d = 0;
	M2p_d = 0;
}

cMT_Detector_GPU::~cMT_Detector_GPU(){
	freeMemory();
}

void cMT_Detector_GPU::SetInputData(sGP &GP_i, int nDeti, sDetCir &DetCirhi){
	freeMemory();

	GP = GP_i;

	nDet = nDeti;
	f_sDetCir_Malloc(nDet, DetCirh);
	memcpy(DetCirh.g2min, DetCirhi.g2min, nDet*cSizeofRD);
	memcpy(DetCirh.g2max, DetCirhi.g2max, nDet*cSizeofRD);

	Tot_h = new double[nDet];
	Coh_h = new double[nDet];

	cudaMalloc((void**)&Tot_d, nDet*cSizeofRD);
	cudaMalloc((void**)&Coh_d, nDet*cSizeofRD);

	cudaMalloc((void**)&M1p_d, 32*32*cSizeofRD);
	cudaMalloc((void**)&M2p_d, 32*32*cSizeofRD);
}

void cMT_Detector_GPU::getDetectorIntensity(double *&aM2Psi, double *&M2aPsi, int ixys, sDetInt *DetInth){
	for(int iDet = 0; iDet<nDet; iDet++)
		f_Sum_MD_Det(GP, aM2Psi, M2aPsi, DetCirh.g2min[iDet], DetCirh.g2max[iDet], M1p_d, M2p_d, iDet, Tot_d, Coh_d);

	cudaMemcpy(Tot_h, Tot_d, nDet*cSizeofRD, cudaMemcpyDeviceToHost);
	cudaMemcpy(Coh_h, Coh_d, nDet*cSizeofRD, cudaMemcpyDeviceToHost);

	for(int iDet = 0; iDet<nDet; iDet++){
		DetInth[iDet].Tot[ixys] = Tot_h[iDet];
		DetInth[iDet].Coh[ixys] = Coh_h[iDet];
	}
}

void cMT_Detector_GPU::getDetectorIntensity(double *&aM2Psi, int ixys, sDetInt *DetInth, bool add){
	for(int iDet = 0; iDet<nDet; iDet++)
		f_Sum_MD_Det(GP, aM2Psi, DetCirh.g2min[iDet], DetCirh.g2max[iDet], M1p_d, iDet, Tot_d);

	cudaMemcpy(Tot_h, Tot_d, nDet*cSizeofRD, cudaMemcpyDeviceToHost);

	for(int iDet = 0; iDet<nDet; iDet++){
		if(add)
			DetInth[iDet].Tot[ixys] += Tot_h[iDet];
		else
			DetInth[iDet].Tot[ixys] = Tot_h[iDet];
	}
}

void cMT_Detector_GPU::getDetectorIntensity(double2 *&aPsi, int ixys, sDetInt *DetInth, bool add){
	for(int iDet = 0; iDet<nDet; iDet++)
		f_Sum_MC_Det(GP, aPsi, DetCirh.g2min[iDet], DetCirh.g2max[iDet], M1p_d, iDet, Tot_d);

	cudaMemcpy(Tot_h, Tot_d, nDet*cSizeofRD, cudaMemcpyDeviceToHost);

	for(int iDet = 0; iDet<nDet; iDet++){
		if(add)
			DetInth[iDet].Tot[ixys] += Tot_h[iDet];
		else
			DetInth[iDet].Tot[ixys] = Tot_h[iDet];
	}
}
