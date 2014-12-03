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
#include "math.h"

#include "hConstTypes.h"
#include "hQuadrature.h"
#include "hMT_General_CPU.h"
#include "hMT_General_GPU.h"
#include "hMT_MGP_CPU.h"
#include "hMT_Specimen_CPU.h"
#include "hMT_AtomTypes_GPU.h"
#include "hMT_Potential_GPU.h"
#include "hMT_Transmission_GPU.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

void cMT_Transmission_GPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	cSynCPU = ccSynCPU;

	if(nSliceM>0)
		if(SliceMTyp==1){
			for(int iSliceM=0; iSliceM<nSliceM; iSliceM++)
				cudaFreen(Trans[iSliceM]);
			delete [] Trans; Trans = 0;
		}else{
			for(int iSliceM=0; iSliceM<nSliceM; iSliceM++)
				cudaFreen(Vpe[iSliceM]);
			delete [] Vpe; Vpe = 0;
		}

	nSliceM = 0;
	SliceMTyp = 0;
}

cMT_Transmission_GPU::cMT_Transmission_GPU()
{
	cSynCPU = ccSynCPU;

	nSliceM = 0;
	SliceMTyp = 0;

	Trans = 0;
	Vpe = 0;
}

cMT_Transmission_GPU::~cMT_Transmission_GPU(){
	freeMemory();
}

void cMT_Transmission_GPU::Potential_Efective(int iSlice, double fPot, float *&Vpe){
	eSlicePos SlicePos = (iSlice==0)?eSPFirst:(iSlice<nSlice)?eSPMedium:eSPLast;

	// Projected potential
	if(iSlice<nSlice) ProjectedPotential(iSlice);

	switch (MT_MGP_CPU->MulOrder)
	{
		case 1:
			f_Potential1(GP, fPot, V0, Vpe);
			break;
		case 2:
			f_Potential2(GP, fPot, V0, V1, V2, SlicePos, Vpe);
			break;
	}
}

void cMT_Transmission_GPU::Transmission(int iSlice, double2 *&Trans){
	eSlicePos SlicePos = (iSlice==0)?eSPFirst:(iSlice<nSlice)?eSPMedium:eSPLast;

	// Projected potential
	if(iSlice<nSlice) ProjectedPotential(iSlice); 

	switch (MT_MGP_CPU->MulOrder)
	{
		case 1:
			if(MT_MGP_CPU->ApproxModel==4) 
				f_TransmissionWPO(PlanPsi, GP, fPot, V0, Trans);
			else 
				f_Transmission1(PlanPsi, GP, fPot, V0, Trans);
			break;
		case 2:
			f_Transmission2(PlanPsi, GP, fPot, MT_Potential_GPU->V0, MT_Potential_GPU->V1, MT_Potential_GPU->V2, SlicePos, Trans);
			break;
	}
}

void cMT_Transmission_GPU::Cal_Trans_Vpe(){
	int iSliceM, nSliceMm = MIN(nSliceM, nSlice);
	if((MT_MGP_CPU->MulOrder==2)&&(nSliceM>nSlice)) nSliceMm++;

	for (iSliceM=0; iSliceM<nSliceMm; iSliceM++){
		if(SliceMTyp==1)
			MT_MulSli_GPU->Transmission(iSliceM, Trans[iSliceM]);
		else
			Potential_Efective(iSliceM, MT_MulSli_GPU->fPot, Vpe[iSliceM]);
	}
	cudaDeviceSynchronize();
}

void cMT_Transmission_GPU::Transmission_Transmit(int iSlice, double2 *&Psi){
	int nSlice = MT_Specimen_CPU->nSlice;
	int nSliceMm = MIN(nSliceM, nSlice);
	if((MT_MGP_CPU->MulOrder==2)&&(nSliceM>nSlice)) nSliceMm++;

	double2 *Transt;
	Transt = (SliceMTyp==1)?Trans[iSlice]:MT_MulSli_GPU->Trans;

	if(iSlice<nSliceMm){
		if(SliceMTyp==2)
			f_Transmission_1_2(MT_MulSli_GPU->PlanPsi, GP, Vpe[iSlice], Transt);
		MT_MulSli_GPU->Transmit(Transt, Psi);
	}else
		MT_MulSli_GPU->Transmission_Transmit(iSlice, Psi);
}

void cMT_Transmission_GPU::SetInputData(cMT_MGP_CPU *MT_MGP_CPU_io, int nAtomsM_i, double *AtomsM_i)
{
	freeMemory();

	cMT_Potential_GPU::SetInputData(MT_MGP_CPU_io, nAtomsM_i, AtomsM_i);

	double Gamma = f_getGamma(MT_MGP_CPU->E0);
	double Lambda = f_getLambda(MT_MGP_CPU->E0);
	fPot = Gamma*Lambda/(cPotf*cos(MT_MGP_CPU->theta));

	int nSliceSigma = ((MT_MGP_CPU->ApproxModel>1)||(MT_MGP_CPU->DimFP%10==0))?0:(int)ceil(6*sigma_max/MT_MGP_CPU->dz);
	int nSliceMax = nSlice + nSliceSigma;
	nSliceMax = (MT_MGP_CPU->MulOrder==1)?nSliceMax:nSliceMax+1;

	size_t SizeFreeM, SizeTotM;
	cudaMemGetInfo(&SizeFreeM, &SizeTotM);
	SizeFreeM = SizeFreeM-10*cMb;
	int nSliceMt;

	if(SizeFreeM/(GP.nxy*cSizeofCD)>=nSliceMax){
		SliceMTyp = 1;
		nSliceMt = SizeFreeM/(GP.nxy*cSizeofCD);
	}else{
		SliceMTyp = 2;
		nSliceMt = SizeFreeM/(GP.nxy*cSizeofRF);
	}

	if((nSliceMt>0)&&(MT_MGP_CPU->SimType==1)&&(MT_MGP_CPU->ApproxModel<=2)){
		nSliceM = MIN(nSliceMt, nSliceMax);
		if(SliceMTyp==1){
			Trans = new double2*[nSliceM];
			for(int iSliceM=0; iSliceM<nSliceM; iSliceM++)
				cudaMalloc((void**)&Trans[iSliceM], GP.nxy*cSizeofCD);
		}else{
			Vpe = new float*[nSliceM];
			for(int iSliceM=0; iSliceM<nSliceM; iSliceM++)
				cudaMalloc((void**)&Vpe[iSliceM], GP.nxy*cSizeofRF);
		}
	}
}
