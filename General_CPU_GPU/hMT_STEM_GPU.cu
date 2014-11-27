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

#include "hMT_STEM_GPU.h"

#include "hConstTypes.h"
#include "hMT_General_CPU.h"
#include "hMT_General_GPU.h"
#include "hMT_Specimen_CPU.h"
#include "hMT_Potential_GPU.h"
#include "hMT_IncidentWave_GPU.h"
#include "hMT_Detector_CPU.h"
#include "hMT_Detector_GPU.h"
#include "hMT_MulSli_GPU.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

void cMT_STEM_GPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	cSynCPU = ccSynCPU;

	MT_MulSli_GPU = 0;
	MT_Specimen_CPU = 0;
	MT_Potential_GPU = 0;
	MT_IncidentWave_GPU = 0;
	delete MT_Detector_GPU; MT_Detector_GPU = 0;

	line = 0;
	FastCal = false;
	ns = 0;	
	x1u = 0;
	y1u = 0;		
	x2u = 0;	
	y2u = 0;	

	f_sDetCir_Free(DetCir);
	for (int iThk = 0; iThk<nThk; iThk++){
		for (int iDet=0; iDet<nDet; iDet++){
			delete [] ImSTEM[iThk].DetInt[iDet].Coh; ImSTEM[iThk].DetInt[iDet].Coh = 0;
			delete [] ImSTEM[iThk].DetInt[iDet].Tot; ImSTEM[iThk].DetInt[iDet].Tot = 0;
		}
		delete [] ImSTEM[iThk].DetInt; ImSTEM[iThk].DetInt = 0;
	}
	delete [] ImSTEM; ImSTEM = 0;
	nDet = 0;

	nThk = 0;
	delete [] Thk; Thk = 0;

	nxs = 0;
	nys = 0;
	delete [] xs; xs = 0;
	delete [] ys; ys = 0;

	nst = 0;
	delete [] xst; xst = 0;
	delete [] yst; yst = 0;

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

cMT_STEM_GPU::cMT_STEM_GPU()
{
	cSynCPU = ccSynCPU;

	MT_MulSli_GPU = 0;
	MT_Specimen_CPU = 0;
	MT_Potential_GPU = 0;
	MT_IncidentWave_GPU = 0;
	MT_Detector_GPU = 0;

	line = 0;
	FastCal = false;
	ns = 0;
	x1u = 0;
	y1u = 0;		
	x2u = 0;	
	y2u = 0;	

	f_sDetCir_Init(DetCir);
	ImSTEM = 0;
	nDet = 0;

	nThk = 0;
	Thk = 0;

	nxs = 0;
	nys = 0;
	xs = 0;
	ys = 0;

	nst = 0;
	xst = 0;
	yst = 0;

	Trans = 0;
	Vpe = 0;

	nSliceM = 0;
	SliceMTyp = 0;
}

cMT_STEM_GPU::~cMT_STEM_GPU(){
	freeMemory();
}

void cMT_STEM_GPU::InitImSTEM()
{
	int iThk, iDet, ist;
	for (iThk = 0; iThk<nThk; iThk++)
		for (iDet=0; iDet<nDet; iDet++)
			for (ist=0; ist<nst; ist++){
				ImSTEM[iThk].DetInt[iDet].Coh[ist] = 0.0;
				ImSTEM[iThk].DetInt[iDet].Tot[ist] = 0.0;
			}
}

void cMT_STEM_GPU::Potential_Efective(int iSlice, double fPot, float *&Vpe){
	int nSlice = MT_Potential_GPU->MT_Specimen_CPU->nSlice;
	eSlicePos SlicePos = (iSlice==0)?eSPFirst:(iSlice<nSlice)?eSPMedium:eSPLast;

	// Projected potential
	if(iSlice<nSlice)
		MT_Potential_GPU->ProjectedPotential(iSlice);

	switch (MT_MGP_CPU.MulOrder)
	{
		case 1:
			f_Potential1(GP, fPot, MT_Potential_GPU->V0, Vpe);
			break;
		case 2:
			f_Potential2(GP, fPot, MT_Potential_GPU->V0, MT_Potential_GPU->V1, MT_Potential_GPU->V2, SlicePos, Vpe);
			break;
	}
}

void cMT_STEM_GPU::Cal_Trans_Vpe(){
	int nSlice = MT_Specimen_CPU->nSlice;
	int iSliceM, nSliceMm = MIN(nSliceM, nSlice);
	if((MT_MGP_CPU.MulOrder==2)&&(nSliceM>nSlice)) nSliceMm++;

	for (iSliceM=0; iSliceM<nSliceMm; iSliceM++){
		if(SliceMTyp==1)
			MT_MulSli_GPU->Transmission(iSliceM, Trans[iSliceM]);
		else
			Potential_Efective(iSliceM, MT_MulSli_GPU->fPot, Vpe[iSliceM]);
	}
	cudaDeviceSynchronize();
}

void cMT_STEM_GPU::Transmission_Transmit(int iSlice, double2 *&Psi){
	int nSlice = MT_Specimen_CPU->nSlice;
	int nSliceMm = MIN(nSliceM, nSlice);
	if((MT_MGP_CPU.MulOrder==2)&&(nSliceM>nSlice)) nSliceMm++;

	double2 *Transt;
	Transt = (SliceMTyp==1)?Trans[iSlice]:MT_MulSli_GPU->Trans;

	if(iSlice<nSliceMm){
		if(SliceMTyp==2)
			f_Transmission_1_2(MT_MulSli_GPU->PlanPsi, GP, Vpe[iSlice], Transt);
		MT_MulSli_GPU->Transmit(Transt, Psi);
	}else
		MT_MulSli_GPU->Transmission_Transmit(iSlice, Psi);
}

void cMT_STEM_GPU::Cal_FAST_STEM_Wavefunction_POA_WPOA(int nConfFP, sDetInt *DetInt){
	int iSlice=0;
	int ist, iThk = 0;
	int iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);
	double nxy2 = pow(double(GP.nxy), 2);

	InitImSTEM();
	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Move atoms
		MT_Specimen_CPU->MoveAtoms(iConf);
		// Transmission
		MT_MulSli_GPU->Transmission(iSlice, MT_MulSli_GPU->Trans);
		for (ist=0; ist<nst; ist++){
			// Plane wave ilumination
			MT_IncidentWave_GPU->Psi0(xst[ist], yst[ist], MT_MulSli_GPU->Psi);
			// Transmit
			MT_MulSli_GPU->Transmit(MT_MulSli_GPU->Trans, MT_MulSli_GPU->Psi);
			// Inclined ilumination
			MT_MulSli_GPU->PhaseMul(MT_MulSli_GPU->Psi);
			// Backward fft2
			cufftExecZ2Z(MT_MulSli_GPU->PlanPsi, MT_MulSli_GPU->Psi, MT_MulSli_GPU->Psi, CUFFT_FORWARD);
			// Add Psi to aM2Psi
			f_Add_wMC2(false, GP, inConfFP/nxy2, MT_MulSli_GPU->Psi, MT_MulSli_GPU->aM2Psi);
			// Detector integration
			MT_Detector_GPU->getDetectorIntensity(MT_MulSli_GPU->aM2Psi, ist, ImSTEM[iThk].DetInt, true);
		}
	}
}

void cMT_STEM_GPU::Cal_FAST_STEM_Wavefunction_MSA(int nConfFP, sDetInt *DetInt){
	int iSlice=0, iSynCPU = 0;
	int ist, iThk = 0;
	int iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);
	double nxy2 = pow(double(GP.nxy), 2);

	InitImSTEM();
	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Move atoms
		MT_Specimen_CPU->MoveAtoms(iConf);

		//Load Trans or Vpe
		Cal_Trans_Vpe();
		for (ist=0; ist<nst; ist++){
			// Plane wave ilumination
			MT_IncidentWave_GPU->Psi0(xst[ist], yst[ist], MT_MulSli_GPU->Psi);
			for (iSlice=0; iSlice<MT_Specimen_CPU->nSlice; iSlice++){
				// Transmission and Transmit
				Transmission_Transmit(iSlice, MT_MulSli_GPU->Psi);
				// Propagate
				MT_MulSli_GPU->Propagate(eSReal, MT_MulSli_GPU->gxu, MT_MulSli_GPU->gyu, MT_Specimen_CPU->get_dz(iSlice), MT_MulSli_GPU->Psi);
				// GPU Synchronize
				f_GPU_Sync_CPU(iSynCPU, cSynCPU); 
			}
			// Last Transmission and Transmit
			if (MT_MGP_CPU.MulOrder==2) Transmission_Transmit(iSlice, MT_MulSli_GPU->Psi);
			// Inclined ilumination
			MT_MulSli_GPU->PhaseMul(MT_MulSli_GPU->Psi);
			// Backward fft2
			cufftExecZ2Z(MT_MulSli_GPU->PlanPsi, MT_MulSli_GPU->Psi, MT_MulSli_GPU->Psi, CUFFT_FORWARD);
			// Add Psi to aM2Psi
			f_Add_wMC2(false, GP, inConfFP/nxy2, MT_MulSli_GPU->Psi, MT_MulSli_GPU->aM2Psi);
			// Detector integration
			MT_Detector_GPU->getDetectorIntensity(MT_MulSli_GPU->aM2Psi, ist, ImSTEM[iThk].DetInt, true);
		}
	}
}

void cMT_STEM_GPU::Cal_STEM(){
	int iThk = 0;

	if(!FastCal){
		for (int ist=0; ist<nst; ist++){
			MT_MulSli_GPU->Image_Convergence_Wave_Illumination(MT_MGP_CPU.nConfFP, eSReciprocal, xst[ist], yst[ist], MT_MulSli_GPU->aPsi, MT_MulSli_GPU->M2aPsi, MT_MulSli_GPU->aM2Psi);
			MT_Detector_GPU->getDetectorIntensity(MT_MulSli_GPU->aM2Psi, MT_MulSli_GPU->M2aPsi, ist, ImSTEM[iThk].DetInt);
		}
	}else{
		if(MT_MGP_CPU.ApproxModel<=2)
			Cal_FAST_STEM_Wavefunction_MSA(MT_MGP_CPU.nConfFP, ImSTEM[iThk].DetInt);
		else
			Cal_FAST_STEM_Wavefunction_POA_WPOA(MT_MGP_CPU.nConfFP, ImSTEM[iThk].DetInt);
	}
}

void cMT_STEM_GPU::SetInputData(cMT_InMULTEM_CPU &MT_InMULTEM_CPU, cMT_MulSli_GPU *MT_MulSli_GPU_i)
{
	freeMemory();

	MT_MulSli_GPU = MT_MulSli_GPU_i;
	MT_Potential_GPU = MT_MulSli_GPU->MT_Potential_GPU;
	MT_Specimen_CPU = MT_Potential_GPU->MT_Specimen_CPU;
	MT_IncidentWave_GPU = MT_MulSli_GPU->MT_IncidentWave_GPU;
	
	MT_MGP_CPU = MT_MulSli_GPU->MT_MGP_CPU;
	GP = MT_MulSli_GPU->GP;

	line = MT_InMULTEM_CPU.STEM_line;
	FastCal = MT_InMULTEM_CPU.STEM_FastCal;
	ns = MT_InMULTEM_CPU.STEM_ns;
	x1u = MT_InMULTEM_CPU.STEM_x1u;	
	y1u = MT_InMULTEM_CPU.STEM_y1u;
	x2u = MT_InMULTEM_CPU.STEM_x2u;
	y2u = MT_InMULTEM_CPU.STEM_y2u;
	f_BuildGrid(line, ns, x1u, y1u, x2u, y2u, nxs, nys, xs, ys);

	nThk = MT_MGP_CPU.nThk;
	if(nThk>0){
		Thk = new double[nThk];
		memcpy(Thk, MT_InMULTEM_CPU.Thickness, nThk*cSizeofRD);
	}

	nDet = MT_InMULTEM_CPU.STEM_nDet;
	double lambda = f_getLambda(MT_InMULTEM_CPU.E0);
	f_sDetCir_Malloc(nDet, DetCir);
	for (int iDet=0; iDet<nDet; iDet++){
		DetCir.g2min[iDet] = pow(MT_InMULTEM_CPU.STEM_DetCir[iDet].InnerAng/lambda, 2);
		DetCir.g2max[iDet] = pow(MT_InMULTEM_CPU.STEM_DetCir[iDet].OuterAng/lambda, 2);
	}

	MT_Detector_GPU = new cMT_Detector_GPU;
	MT_Detector_GPU->SetInputData(GP, nDet, DetCir);

	nst = (line==1)?ns:nxs*nys;
	int ils, ixs, iys, ixys;
	xst = new double[nst];
	yst = new double[nst];
	if(line==1){
		for (ils=0; ils<ns; ils++){
			xst[ils] = xs[ils];
			yst[ils] = ys[ils];
		}
	}else{
		for (ixs=0; ixs<nxs; ixs++)
			for (iys=0; iys<nys; iys++){
				ixys = ixs*nys + iys;
				xst[ixys] = xs[ixs];
				yst[ixys] = ys[iys];
			}
	}

	int iThk, iDet, ist;
	ImSTEM = new sImSTEM[nThk];
	for (iThk = 0; iThk<nThk; iThk++){
		ImSTEM[iThk].DetInt = new sDetInt[nDet];
		for (iDet=0; iDet<nDet; iDet++){
			ImSTEM[iThk].DetInt[iDet].Coh = new double[nst];
			ImSTEM[iThk].DetInt[iDet].Tot = new double[nst];
			for (ist=0; ist<nst; ist++){
				ImSTEM[iThk].DetInt[iDet].Coh[ist] = 0.0;
				ImSTEM[iThk].DetInt[iDet].Tot[ist] = 0.0;
			}
		}
	}

	/****************************************************************/
	int nSliceSigma = ((MT_MGP_CPU.ApproxModel>1)||(MT_MGP_CPU.DimFP%10==0))?0:(int)ceil(6*MT_Specimen_CPU->sigma_max/MT_MGP_CPU.dz);
	int nSliceMax = MT_Specimen_CPU->nSlice + nSliceSigma;
	nSliceMax = (MT_MGP_CPU.MulOrder==1)?nSliceMax:nSliceMax+1;

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

	if((FastCal)&&(nSliceMt>0)&&(MT_MGP_CPU.SimType==1)&&(MT_MGP_CPU.ApproxModel<=2)){
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
