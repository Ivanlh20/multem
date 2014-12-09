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
#include "hMT_InMulSli_CPU.h"
#include "hMT_MGP_CPU.h"
#include "hMT_Transmission_GPU.h"
#include "hMT_IncidentWave_GPU.h"
#include "hMT_MicroscopeEffects_GPU.h"
#include "hMT_STEM_GPU.h"
#include "hMT_MulSli_GPU.h"

#include "math.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

void cMT_MulSli_GPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	cSynCPU = ccSynCPU;

	f_sGP_Init(GP);
	f_sLens_Init(Lens);
	
	gxu = 0;
	gyu = 0;

	/****************************************************************/
	delete STEM;

	CBED.x0 = 0;
	CBED.y0 = 0;

	//HRTEM.xx;

	PED.nrot = 0;
	PED.theta = 0;

	HCI.nrot = 0;
	HCI.theta = 0;
	//HCI.xx;

	//EWRS.xx;

	//EWFS.xx;
	/****************************************************************/

	f_sACD_cudaFree(ExpRg_x);
	f_sACD_cudaFree(ExpRg_x);

	f_sACD_cudaFree(Prop_x);
	f_sACD_cudaFree(Prop_y);

	cudaFreen(Psi);
	cudaFreen(aPsi);

	cudaFreen(M2aPsi);
	cudaFreen(aM2Psi);

	cufftDestroyn(PlanPsi);

	delete MT_Transmission_GPU; MT_Transmission_GPU = 0;
	delete MT_IncidentWave_GPU; MT_IncidentWave_GPU = 0;
	delete MT_MicroscopeEffects_GPU; MT_MicroscopeEffects_GPU = 0;
}

cMT_MulSli_GPU::cMT_MulSli_GPU(){
	cSynCPU = ccSynCPU;

	f_sGP_Init(GP);
	f_sLens_Init(Lens);
	
	gxu = 0;		
	gyu = 0;		

	/****************************************************************/
	STEM = 0;

	CBED.x0 = 0;
	CBED.y0 = 0;

	//HRTEM.xx

	PED.nrot = 0;
	PED.theta = 0;

	HCI.nrot = 0;
	HCI.theta = 0;
	//HCI.xx

	//EWRS.xx

	//EWFS.xx
	/****************************************************************/

	f_sACD_cudaInit(ExpRg_x);
	f_sACD_cudaInit(ExpRg_x);

	f_sACD_cudaInit(Prop_x);
	f_sACD_cudaInit(Prop_y);

	Psi = 0;
	aPsi = 0;

	M2aPsi = 0;
	aM2Psi = 0;

	PlanPsi = 0;

	MT_Transmission_GPU = 0;
	MT_IncidentWave_GPU = 0;
	MT_MicroscopeEffects_GPU = 0;
}

cMT_MulSli_GPU::~cMT_MulSli_GPU(){
	freeMemory();
	cudaDeviceReset();
}

// From Device To Host
void cMT_MulSli_GPU::GPU2CPU(double2 *&Psid, sComplex &Psih){	
	// fft2shift 
	f_fft2Shift_MC(GP, Psid);
	/*********************copy data to host************************/
	f_Copy_MCd(GP, Psid, MT_Transmission_GPU->V0, MT_Transmission_GPU->V1, Psih);
}

// From Device To Host
void cMT_MulSli_GPU::GPU2CPU(double2 *&Psid, double *&M2Psid, sComplex &Psih, double *&M2Psih){
	// fft2shift 
	f_fft2Shift_MC_MD(GP, Psid, M2Psid);
	// Copy wave function squared to the host
	f_Copy_MCd_MDd(GP, Psid, M2Psid, MT_Transmission_GPU->V0, MT_Transmission_GPU->V1, Psih, M2Psih);
}

// From Device To Host
void cMT_MulSli_GPU::GPU2CPU(double2 *&Psid, double *&M2Psi1d, double *&M2Psi2d, sComplex &Psih, double *&M2Psi1h, double *&M2Psi2h){
	// fft2shift 
	f_fft2Shift_MC_MD(GP, Psid, M2Psi1d, M2Psi2d);
	/*********************copy data to host************************/
	f_Copy_MCd_MDd(GP, Psid, M2Psi1d, M2Psi2d, MT_Transmission_GPU->V0, MT_Transmission_GPU->V1, Psih, M2Psi1h, M2Psi2h);
}

void cMT_MulSli_GPU::PhaseMul(double gxu, double gyu, double2 *&Psi){
	if ((gxu != 0.0)||(gyu != 0.0))
		f_PhaseMul(GP, gxu, gyu, ExpRg_x, ExpRg_y, Psi);
}

void cMT_MulSli_GPU::PhaseMul(double2 *&Psi){
	if ((gxu != 0.0)||(gyu != 0.0))
		f_PhaseMul(GP, gxu, gyu, ExpRg_x, ExpRg_y, Psi);
}

void cMT_MulSli_GPU::Propagate(eSpace Space, double gxu, double gyu, double z, double2 *&Psi){
	f_Propagate(PlanPsi, GP, Space, gxu, gyu, Lens.lambda, z, Prop_x, Prop_y, Psi);
}

// Phase object approximation: Space :1 Real and 2 Fourier
void cMT_MulSli_GPU::Cal_Wavefunction_POA_WPOA(eSpace Space, double2 *&Psi){
	int iSlice = 0;
	// get Transmission and Transmit
	MT_Transmission_GPU->Transmit(iSlice, Psi);
	// Inclined ilumination
	PhaseMul(Psi);
	// Backward fft2
	if(Space==eSReciprocal){
		cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
		f_Scale_MC(GP, GP.inxy, Psi);
	}
}

// Exit wave calculation: Space :1 Real and 2 Fourier
void cMT_MulSli_GPU::Cal_Wavefunction_MSA(eSpace Space, double2 *&Psi){
	int iSlice = 0, iSynCPU = 0;

	for (iSlice = 0; iSlice<MT_Transmission_GPU->nSlice; iSlice++){
		// get Transmission and Transmit
		MT_Transmission_GPU->Transmit(iSlice, Psi);
		// Propagate
		Propagate(eSReal, gxu, gyu, MT_Transmission_GPU->get_dz(iSlice), Psi);
		// GPU Synchronize
		f_GPU_Sync_CPU(iSynCPU, cSynCPU); 
	}
	// Last Transmission and Transmit
	if (MT_MGP_CPU.MulOrder==2) MT_Transmission_GPU->Transmit(iSlice, Psi);;
	// Inclined ilumination
	PhaseMul(Psi);
	// Back propagation to our plane reference
	Propagate(Space, gxu, gyu, MT_Transmission_GPU->z_BackProp, Psi);
}

// Exit wave calculation: Space :1 Real and 2 Fourier
void cMT_MulSli_GPU::Cal_Wavefunction(eSpace Space, double2 *&Psi){
	if(MT_MGP_CPU.ApproxModel<=2)
		Cal_Wavefunction_MSA(Space, Psi);
	else
		Cal_Wavefunction_POA_WPOA(Space, Psi);
}

// Get wave function for Plane wave ilumination
void cMT_MulSli_GPU::Image_Plane_Wave_Illumination(int nConfFP, eSpace Space, int MEffect, int STEffect, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi){
	int iConf0 = (nConfFP==0)?0:1;
	bool IsNotSharp = (nConfFP==0)?false:true;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);

	f_Set_MC_MD(GP, 0.0, 0.0, aPsi, 0.0, aM2Psi);

	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Move atoms
		MT_Transmission_GPU->MoveAtoms(iConf);
		// Plane wave ilumination
		MT_IncidentWave_GPU->Psi0(Psi);
		if(MEffect==0){
			// Exit wave
			Cal_Wavefunction(Space, Psi);
			// Add Psi and M2aPsi to aPsi and aM2Psi
			f_Add_wMC_wMD(IsNotSharp, GP, inConfFP, Psi, aPsi, aM2Psi);
		}else{
			// Exit wave(g)
			Cal_Wavefunction(eSReciprocal, Psi);
			// Microscope effects
			MT_MicroscopeEffects_GPU->ApplyMEffects(MEffect, STEffect, Psi, M2aPsi);
			// Backward fft2
			cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_INVERSE);
			// Add Psi and M2aPsi to aPsi and aM2Psi
			f_Add_wMC_wMD(IsNotSharp, GP, inConfFP, Psi, M2aPsi, aPsi, aM2Psi);
		}
	}

	if(MEffect==0)
		f_Add_wMC2(false, GP, 1.0, aPsi, M2aPsi);
	else{
		// Forward fft2
		cufftExecZ2Z(PlanPsi, aPsi, aPsi, CUFFT_FORWARD);
		// Scale vector
		f_Scale_MC(GP, GP.inxy, aPsi);
		// Microscope effects
		MT_MicroscopeEffects_GPU->ApplyMEffects(MEffect, STEffect, aPsi, M2aPsi);
		// Backward fft2
		cufftExecZ2Z(PlanPsi, aPsi, aPsi, CUFFT_INVERSE);
	}

}

// Get wave function for convergent beam ilumination
void cMT_MulSli_GPU::Image_Convergence_Wave_Illumination(int nConfFP, eSpace Space, double xi, double yi, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi){
	int iConf0 = (nConfFP==0)?0:1;
	bool IsNotSharp = (nConfFP==0)?false:true;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);

	if(IsNotSharp)
		f_Set_MC_MD(GP, 0.0, 0.0, aPsi, 0.0, aM2Psi);

	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Move atoms
		MT_Transmission_GPU->MoveAtoms(iConf);
		// Plane wave ilumination
		MT_IncidentWave_GPU->Psi0(xi, yi, Psi);
		// Exit wave
		Cal_Wavefunction(Space, Psi);
		// Add Psi and M2aPsi to aPsi and aM2Psi
		f_Add_wMC_wMD(IsNotSharp, GP, inConfFP, Psi, aPsi, aM2Psi);
	}
	f_Add_wMC2(false, GP, 1.0, aPsi, M2aPsi);
}

// Set Input data
void cMT_MulSli_GPU::SetInputData(cMT_InMulSli_CPU &MT_InMulSli_CPU){
	freeMemory();

	MT_MGP_CPU.SetInputData(MT_InMulSli_CPU);

	cudaSetDevice(MT_MGP_CPU.gpu);

	f_sGP_SetInputData(&MT_MGP_CPU, GP);
	f_sLens_SetInputData(MT_InMulSli_CPU, GP, Lens);

	gxu = sin(MT_MGP_CPU.theta)*cos(MT_MGP_CPU.phi)/Lens.lambda;
	gyu = sin(MT_MGP_CPU.theta)*sin(MT_MGP_CPU.phi)/Lens.lambda;

	f_sACD_cudaMalloc(GP.nx, ExpRg_x);
	f_sACD_cudaMalloc(GP.ny, ExpRg_y);
	f_sACD_cudaMalloc(GP.nx, Prop_x);
	f_sACD_cudaMalloc(GP.ny, Prop_y);

	cudaMalloc((void**)&Psi, GP.nxy*cSizeofCD);
	cudaMalloc((void**)&aPsi, GP.nxy*cSizeofCD);

	cudaMalloc((void**)&M2aPsi, GP.nxy*cSizeofRD);
	cudaMalloc((void**)&aM2Psi, GP.nxy*cSizeofRD);

	cufftPlan2d(&PlanPsi, GP.ny, GP.nx, CUFFT_Z2Z);

	// Transmission parameters
	MT_Transmission_GPU = new cMT_Transmission_GPU ;
	MT_Transmission_GPU->SetInputData(&MT_MGP_CPU, PlanPsi,  MT_InMulSli_CPU.nAtomsM, MT_InMulSli_CPU.AtomsM);

	// Incident wave parameters
	MT_IncidentWave_GPU = new cMT_IncidentWave_GPU;
	MT_IncidentWave_GPU->SetInputData(GP, Lens, PlanPsi);

	// Microscope parameters
	MT_MicroscopeEffects_GPU = new cMT_MicroscopeEffects_GPU;
	MT_MicroscopeEffects_GPU->SetInputData(GP, Lens, PlanPsi, MT_Transmission_GPU->Trans0);

	/***************************************************************************/
	/***************************************************************************/

	switch (MT_InMulSli_CPU.SimType){
		case 11:		// STEM
			STEM = new cMT_STEM_GPU;
			STEM->SetInputData(MT_InMulSli_CPU, this);
			break;
		case 12:		// ISTEM

			break;
		case 21:		// CBED
			CBED.x0 = MT_InMulSli_CPU.CBED_x0;
			CBED.y0 = MT_InMulSli_CPU.CBED_y0;
			break;
		case 22:		// CBEI
			CBEI.x0 = MT_InMulSli_CPU.CBEI_x0;
			CBEI.y0 = MT_InMulSli_CPU.CBEI_y0;
			break;
		case 31:		// ED

			break;
		case 32:		// HRTEM

			break;
		case 41:		// PED
			PED.nrot = MT_InMulSli_CPU.PED_nrot;
			PED.theta = MT_InMulSli_CPU.PED_theta;
			break;
		case 42:		// HCI
			HCI.nrot = MT_InMulSli_CPU.HCI_nrot;
			HCI.theta = MT_InMulSli_CPU.HCI_theta;
			break;
		case 51:	// EW real

			break;
		case 52:	// EW real
			break;
	}
}

void cMT_MulSli_GPU::Cal_ExitWaveRS(sComplex &aPsih, double *&aM2Psih){
	int MEffect = 0, STEffect = 0;
	eSpace Space = eSReal;
	Image_Plane_Wave_Illumination(MT_MGP_CPU.nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	GPU2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_ExitWaveFS(sComplex &aPsih, double *&aM2Psih){
	int MEffect = 0, STEffect = 0;
	eSpace Space = eSReciprocal;
	Image_Plane_Wave_Illumination(MT_MGP_CPU.nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	GPU2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_ED(sComplex &aPsih, double *&aM2Psih){
	int MEffect = 0, STEffect = 0;
	eSpace Space = eSReciprocal;
	Image_Plane_Wave_Illumination(MT_MGP_CPU.nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	GPU2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_HRTEM(sComplex &aPsih, double *&M2aPsih, double *&aM2Psih){
	eSpace Space = eSReal;
	Image_Plane_Wave_Illumination(MT_MGP_CPU.nConfFP, Space, MT_MGP_CPU.MEffect, MT_MGP_CPU.STEffect, aPsi, M2aPsi, aM2Psi);
	GPU2CPU(aPsi, M2aPsi, aM2Psi, aPsih, M2aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_CBED(sComplex &aPsih, double *&aM2Psih){
	Image_Convergence_Wave_Illumination(MT_MGP_CPU.nConfFP, eSReciprocal, CBED.x0, CBED.y0, aPsi, M2aPsi, aM2Psi);
	GPU2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_CBEI(sComplex &aPsih, double *&aM2Psih){
	Image_Convergence_Wave_Illumination(MT_MGP_CPU.nConfFP, eSReal, CBED.x0, CBED.y0, aPsi, M2aPsi, aM2Psi);
	GPU2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

// Hollow cone ilumination
void cMT_MulSli_GPU::CAL_HCI(sComplex &aPsih, double *&aM2Psih){
	int iSlice = 0, iSynCPU = 0;
	int irot, iThk = 0;
	int nConfFP = MT_MGP_CPU.nConfFP, iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);
	double gu, phi, ecos, esin;
	double w = inConfFP/double(HCI.nrot);

	f_Set_MD(GP, 0.0, aM2Psi);
	phi = 2.0*cPi/double(HCI.nrot); gu = sin(HCI.theta)/Lens.lambda;
	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Move atoms
		MT_Transmission_GPU->MoveAtoms(iConf);
		for (irot=0; irot<HCI.nrot; irot++){
			sincos(irot*phi, &esin, &ecos); 
			gxu = -gu*ecos; gyu = -gu*esin;
			// Plane wave ilumination
			MT_IncidentWave_GPU->Psi0(Psi);
			// Exit wave(g)
			Cal_Wavefunction(eSReciprocal, Psi);
			// Microscope effects
			MT_MicroscopeEffects_GPU->ApplyMEffects(MT_MGP_CPU.MEffect, MT_MGP_CPU.STEffect, Psi, M2aPsi);
			// Backward fft2
			cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_INVERSE);
			// Add Psi and M2aPsi to aPsi and aM2Psi
			f_Add_wMC_wMD(true, GP, w, Psi, M2aPsi, aPsi, aM2Psi);
		}
	}
	GPU2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_STEM(){

}

//void cMT_MulSli_GPU::Cal_STEM(){
//	int ist, iThk = 0;
//	double *M2aPsih, *aM2Psih;
//	aM2Psih = new double[GP.nxy];
//	M2aPsih = new double[GP.nxy];
//	cMT_Detector_CPU MT_Detector_CPU;
//	MT_Detector_CPU.SetInputData(GP, STEM.nDet, STEM.DetCir);
//
//	for (ist=0; ist<STEM.nst; ist++){
//		get_Imagefuncion(MT_MGP_CPU.nConfFP, eSReciprocal, STEM.xs[ist], STEM.ys[ist], aPsi, M2aPsi, aM2Psi);
//		cudaMemcpy(aM2Psih, aM2Psi, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
//		cudaMemcpy(M2aPsih, M2aPsi, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
//		MT_Detector_CPU.getDetectorIntensity(aM2Psih, M2aPsih, ist, STEM.ImSTEM[iThk].DetInt);
//	}
//
//	delete [] aM2Psih; aM2Psih = 0;
//	delete [] M2aPsih; M2aPsih = 0;
//}