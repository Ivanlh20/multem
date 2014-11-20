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

#include <cstring>
#include "math.h"
#include "hConstTypes.h"
#include "hMT_General_CPU.h"
#include "hMT_General_GPU.h"
#include "hMT_Potential_GPU.h"
#include "hMT_IncidentWave_GPU.h"
#include "hMT_Detector_CPU.h"
#include "hMT_Detector_GPU.h"
#include "hMT_STEM_GPU.h"
#include "hMT_MulSli_GPU.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

void cMT_MulSli_GPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	cSynCPU = ccSynCPU;

	f_sMPG_Init(MGP);
	f_sGP_Init(GP);
	f_sLens_Init(Lens);
	
	gxu = 0;
	gyu = 0;

	fPot = 0;

	/*****************************************************************/

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
	/*****************************************************************/

	f_sACD_cudaFree(ExpRg_x);
	f_sACD_cudaFree(ExpRg_x);

	f_sACD_cudaFree(Prop_x);
	f_sACD_cudaFree(Prop_y);

	cudaFreen(Trans);

	cudaFreen(Psi);
	cudaFreen(aPsi);

	cudaFreen(M2aPsi);
	cudaFreen(aM2Psi);

	cufftDestroyn(PlanPsi);
}

cMT_MulSli_GPU::cMT_MulSli_GPU(){
	cSynCPU = ccSynCPU;

	f_sMPG_Init(MGP);
	f_sGP_Init(GP);
	f_sLens_Init(Lens);
	
	gxu = 0;		
	gyu = 0;		

	fPot = 0;

	/*****************************************************************/

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
	/*****************************************************************/

	f_sACD_cudaInit(ExpRg_x);
	f_sACD_cudaInit(ExpRg_x);

	f_sACD_cudaInit(Prop_x);
	f_sACD_cudaInit(Prop_y);

	Trans = 0;

	Psi = 0;
	aPsi = 0;

	M2aPsi = 0;
	aM2Psi = 0;

	PlanPsi = 0;
}

cMT_MulSli_GPU::~cMT_MulSli_GPU(){
	freeMemory();
	cudaDeviceReset();
}

// From Device To Host
void cMT_MulSli_GPU::GPU2CPU(double2 *&Psid, sComplex &Psih){	
	// fft2shift 
	f_fft2Shift_MC(GP, Psid);
	/**********************copy data to host************************/
	f_Copy_MCd(GP, Psid, MT_Potential_GPU.V0, MT_Potential_GPU.V1, Psih);
}

// From Device To Host
void cMT_MulSli_GPU::GPU2CPU(double2 *&Psid, double *&M2Psid, sComplex &Psih, double *&M2Psih){
	// fft2shift 
	f_fft2Shift_MC_MD(GP, Psid, M2Psid);
	// Copy wave function squared to the host
	f_Copy_MCd_MDd(GP, Psid, M2Psid, MT_Potential_GPU.V0, MT_Potential_GPU.V1, Psih, M2Psih);
}

// From Device To Host
void cMT_MulSli_GPU::GPU2CPU(double2 *&Psid, double *&M2Psi1d, double *&M2Psi2d, sComplex &Psih, double *&M2Psi1h, double *&M2Psi2h){
	// fft2shift 
	f_fft2Shift_MC_MD(GP, Psid, M2Psi1d, M2Psi2d);
	/**********************copy data to host************************/
	f_Copy_MCd_MDd(GP, Psid, M2Psi1d, M2Psi2d, MT_Potential_GPU.V0, MT_Potential_GPU.V1, Psih, M2Psi1h, M2Psi2h);
}

void cMT_MulSli_GPU::PhaseMul(double gxu, double gyu, double2 *&Psi){
	if ((gxu != 0.0)||(gyu != 0.0))
		f_PhaseMul(GP, gxu, gyu, ExpRg_x, ExpRg_y, Psi);
}

void cMT_MulSli_GPU::Transmission(int iSlice, double fPot, double2 *&Trans){
	int nSlice = MT_Potential_GPU.MT_Specimen_CPU.nSlice;
	eSlicePos SlicePos = (iSlice==0)?eSPFirst:(iSlice<nSlice)?eSPMedium:eSPLast;

	// Projected potential
	if(iSlice<nSlice)
		MT_Potential_GPU.ProjectedPotential(iSlice);

	switch (MGP.MulOrder)
	{
		case 1:
			if(MGP.ApproxModel==4)
				f_TransmissionWPO(PlanPsi, GP, fPot, MT_Potential_GPU.V0, Trans);
			else
				f_Transmission1(PlanPsi, GP, fPot, MT_Potential_GPU.V0, Trans);
			break;
		case 2:
			f_Transmission2(PlanPsi, GP, fPot, MT_Potential_GPU.V0, MT_Potential_GPU.V1, MT_Potential_GPU.V2, SlicePos, Trans);
			break;
	}
}

void cMT_MulSli_GPU::Transmit(double2 *&Trans, double2 *&Psi){
	f_Transmit(GP, Trans, Psi);
}

void cMT_MulSli_GPU::Propagate(eSpace Space, double gxu, double gyu, double z, double2 *&Psi){
	f_Propagate(PlanPsi, GP, Space, gxu, gyu, Lens.lambda, z, Prop_x, Prop_y, Psi);
}
 
/****************************************************************************/
/****************************************************************************/

// Phase object approximation: Space :1 Real and 2 Fourier
void cMT_MulSli_GPU::Cal_Wavefunction_POA_WPOA(eSpace Space, double2 *&Psi){
	int iSlice = 0;
	// Transmission
	Transmission(iSlice, fPot, Trans);
	// Transmit
	Transmit(Trans, Psi);
	// Inclined ilumination
	PhaseMul(gxu, gyu, Psi);
	// Backward fft2
	if(Space==eSReciprocal){
		cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
		f_Scale_MC(GP, GP.inxy, Psi);
	}
}

// Phase object approximation: Space :1 Real and 2 Fourier
void cMT_MulSli_GPU::Cal_Wavefunction_PA(eSpace Space, double2 *&Psi){
	int iSlice = 0, iSynCPU = 0;

	for (iSlice = 0; iSlice<MT_Potential_GPU.MT_Specimen_CPU.nSlice; iSlice++){
		// Transmission
		Transmission(iSlice, fPot, Trans);
		// Transmit
		Transmit(Trans, Psi);
		// Propagate
		if(MT_Potential_GPU.MT_Specimen_CPU.nSlice>1)
			Propagate(eSReal, gxu, gyu, MT_Potential_GPU.MT_Specimen_CPU.get_dz(iSlice), Psi);
		// GPU Synchronize
		f_GPU_Sync_CPU(iSynCPU, cSynCPU); 
	}
	// Inclined ilumination
	PhaseMul(gxu, gyu, Psi);
	// Backward fft2
	if(Space==eSReciprocal){
		cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
		f_Scale_MC(GP, GP.inxy, Psi);
	}
}

// Exit wave calculation: Space :1 Real and 2 Fourier
void cMT_MulSli_GPU::Cal_Wavefunction_MSA(eSpace Space, double2 *&Psi){
	int iSlice = 0, iSynCPU = 0;

	for (iSlice = 0; iSlice<MT_Potential_GPU.MT_Specimen_CPU.nSlice; iSlice++){
		// Transmission
		Transmission(iSlice, fPot, Trans);
		// Transmit
		Transmit(Trans, Psi);
		// Propagate
		Propagate(eSReal, gxu, gyu, MT_Potential_GPU.MT_Specimen_CPU.get_dz(iSlice), Psi);
		// GPU Synchronize
		f_GPU_Sync_CPU(iSynCPU, cSynCPU); 
	}
	// Last transmission function
	if (MGP.MulOrder==2){
		// Transmission
		Transmission(iSlice, fPot, Trans);
		// Transmit
		Transmit(Trans, Psi);
	}
	// Inclined ilumination
	PhaseMul(gxu, gyu, Psi);
	// Back propagation to our plane reference
	Propagate(Space, gxu, gyu, MT_Potential_GPU.MT_Specimen_CPU.z_BackProp, Psi);
}

/****************************************************************************/
/****************************************************************************/

// (Weak) Phase object approximation Space :1 Real and 2 Fourier
void cMT_MulSli_GPU::Cal_FAST_STEM_Wavefunction_POA_WPOA(int nConfFP, sDetInt *DetInt){
	int iSlice = 0;
	int ist, iThk = 0;
	int iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);
	double nxy2 = pow(double(GP.nxy), 2);

	cMT_Detector_GPU MT_Detector_GPU;
	MT_Detector_GPU.SetInputData(GP, STEM.nDet, STEM.DetCir);

	STEM.InitImSTEM();
	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Move atoms
		MT_Potential_GPU.MT_Specimen_CPU.MoveAtoms(iConf);
		// Transmission
		Transmission(iSlice, fPot, Trans);
		for (ist=0; ist<STEM.nst; ist++){
			// Plane wave ilumination
			MT_IncidentWave_GPU.Psi0(STEM.xst[ist], STEM.yst[ist], Psi);
			// Transmit
			Transmit(Trans, Psi);
			// Inclined ilumination
			PhaseMul(gxu, gyu, Psi);
			// Backward fft2
			cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
			// Add Psi to aM2Psi
			f_Add_wMC2(false, GP, inConfFP/nxy2, Psi, aM2Psi);
			// Detector integration
			MT_Detector_GPU.getDetectorIntensity(aM2Psi, ist, STEM.ImSTEM[iThk].DetInt, true);
		}
	}
}

// Phase object approximation: Space :1 Real and 2 Fourier
void cMT_MulSli_GPU::Cal_FAST_STEM_Wavefunction_PA(int nConfFP, sDetInt *DetInt){
	int iSlice = 0, iSynCPU = 0;
	int ist, iThk = 0;
	int iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);
	double nxy2 = pow(double(GP.nxy), 2);
	int iSliceM, nSliceM;

	cMT_Detector_GPU MT_Detector_GPU;
	MT_Detector_GPU.SetInputData(GP, STEM.nDet, STEM.DetCir);

	STEM.InitImSTEM();
	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Move atoms
		MT_Potential_GPU.MT_Specimen_CPU.MoveAtoms(iConf);
		nSliceM = MIN(naTrans, MT_Potential_GPU.MT_Specimen_CPU.nSlice);
		for (iSliceM = 0; iSliceM<nSliceM; iSliceM++)
			Transmission(iSliceM, fPot, aTrans[iSliceM]);
		cudaDeviceSynchronize();

		for (ist=0; ist<STEM.nst; ist++){
			// Plane wave ilumination
			MT_IncidentWave_GPU.Psi0(STEM.xst[ist], STEM.yst[ist], Psi);

			for (iSlice = 0; iSlice<MT_Potential_GPU.MT_Specimen_CPU.nSlice; iSlice++){
				if(iSlice<nSliceM)
					Transmit(aTrans[iSlice], Psi);
				else{
					// Transmission
					Transmission(iSlice, fPot, Trans);
					// Transmit
					Transmit(Trans, Psi);
				}
				if(MT_Potential_GPU.MT_Specimen_CPU.nSlice>1)
					Propagate(eSReal, gxu, gyu, MT_Potential_GPU.MT_Specimen_CPU.get_dz(iSlice), Psi);
				// GPU Synchronize
				f_GPU_Sync_CPU(iSynCPU, cSynCPU); 
			}
			// Inclined ilumination
			PhaseMul(gxu, gyu, Psi);
			// Backward fft2
			cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
			// Add Psi to aM2Psi
			f_Add_wMC2(false, GP, inConfFP/nxy2, Psi, aM2Psi);
			// Detector integration
			MT_Detector_GPU.getDetectorIntensity(aM2Psi, ist, STEM.ImSTEM[iThk].DetInt, true);
		}
	}
}

// Exit wave calculation: Space :1 Real and 2 Fourier
void cMT_MulSli_GPU::Cal_FAST_STEM_Wavefunction_MSA(int nConfFP, sDetInt *DetInt){
	int iSlice = 0, iSynCPU = 0;
	int ist, iThk = 0;
	int iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);
	double nxy2 = pow(double(GP.nxy), 2);
	int iSliceM, nSliceM;

	cMT_Detector_GPU MT_Detector_GPU;
	MT_Detector_GPU.SetInputData(GP, STEM.nDet, STEM.DetCir);

	STEM.InitImSTEM();
	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Move atoms
		MT_Potential_GPU.MT_Specimen_CPU.MoveAtoms(iConf);
		nSliceM = MIN(naTrans, MT_Potential_GPU.MT_Specimen_CPU.nSlice);
		for (iSliceM = 0; iSliceM<nSliceM; iSliceM++)
			Transmission(iSliceM, fPot, aTrans[iSliceM]);
		cudaDeviceSynchronize();

		for (ist=0; ist<STEM.nst; ist++){
			// Plane wave ilumination
			MT_IncidentWave_GPU.Psi0(STEM.xst[ist], STEM.yst[ist], Psi);

			for (iSlice = 0; iSlice<MT_Potential_GPU.MT_Specimen_CPU.nSlice; iSlice++){
				if(iSlice<nSliceM)
					Transmit(aTrans[iSlice], Psi);
				else{
					// Transmission
					Transmission(iSlice, fPot, Trans);
					// Transmit
					Transmit(Trans, Psi);
				}
				// Propagate
				Propagate(eSReal, gxu, gyu, MT_Potential_GPU.MT_Specimen_CPU.get_dz(iSlice), Psi);
				// GPU Synchronize
				f_GPU_Sync_CPU(iSynCPU, cSynCPU); 
			}
			// Last transmission function
			if (MGP.MulOrder==2){
				// Transmission
				Transmission(iSlice, fPot, Trans);
				// Transmit
				Transmit(Trans, Psi);
			}
			// Inclined ilumination
			PhaseMul(gxu, gyu, Psi);
			// Backward fft2
			cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
			// Add Psi to aM2Psi
			f_Add_wMC2(false, GP, inConfFP/nxy2, Psi, aM2Psi);
			// Detector integration
			MT_Detector_GPU.getDetectorIntensity(aM2Psi, ist, STEM.ImSTEM[iThk].DetInt, true);
		}
	}
}

/****************************************************************************/
/****************************************************************************/

// Exit wave calculation: Space :1 Real and 2 Fourier
void cMT_MulSli_GPU::Cal_Wavefunction(eSpace Space, double2 *&Psi){
	switch(MGP.ApproxModel)
	{
		case 1:
			Cal_Wavefunction_MSA(Space, Psi);
			break;
		case 2:
			Cal_Wavefunction_PA(Space, Psi);
			break;
		case 3:
			Cal_Wavefunction_POA_WPOA(Space, Psi);
			break;
		case 4:
			Cal_Wavefunction_POA_WPOA(Space, Psi);
			break;
	}
}

/***************************************************************************/
// Get wave function for inclined ilumination by using frozen phonon
void cMT_MulSli_GPU::get_Imagefuncion(int nConfFP, eSpace Space, int MEffect, int STEffect, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi){
	int iConf0 = (nConfFP==0)?0:1;
	bool IsNotSharp = (nConfFP==0)?false:true;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);

	if(IsNotSharp)
		f_Set_MC_MD(GP, 0.0, 0.0, aPsi, 0.0, aM2Psi);

	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Move atoms
		MT_Potential_GPU.MT_Specimen_CPU.MoveAtoms(iConf);
		// Plane wave ilumination
		MT_IncidentWave_GPU.Psi0(Psi);
		if(MEffect==0){
			// Exit wave
			Cal_Wavefunction(Space, Psi);
			// Add Psi and M2aPsi to aPsi and aM2Psi
			f_Add_wMC_wMD(IsNotSharp, GP, inConfFP, Psi, aPsi, aM2Psi);
		}else{
			// Exit wave(g)
			Cal_Wavefunction(eSReciprocal, Psi);
			// Microscope effects
			MT_MicroscopeEffects_GPU.ApplyMEffects(MEffect, STEffect, Psi, M2aPsi);
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
		MT_MicroscopeEffects_GPU.ApplyMEffects(MEffect, STEffect, aPsi, M2aPsi);
		// Backward fft2
		cufftExecZ2Z(PlanPsi, aPsi, aPsi, CUFFT_INVERSE);
	}

}

// Get wave function for convergent beam and inclined ilumination using frozen phonon
void cMT_MulSli_GPU::get_Imagefuncion(int nConfFP, eSpace Space, double xi, double yi, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi){
	int iConf0 = (nConfFP==0)?0:1;
	bool IsNotSharp = (nConfFP==0)?false:true;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);

	if(IsNotSharp)
		f_Set_MC_MD(GP, 0.0, 0.0, aPsi, 0.0, aM2Psi);

	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Move atoms
		MT_Potential_GPU.MT_Specimen_CPU.MoveAtoms(iConf);
		// Plane wave ilumination
		MT_IncidentWave_GPU.Psi0(xi, yi, Psi);
		// Exit wave
		Cal_Wavefunction(Space, Psi);
		// Add Psi and M2aPsi to aPsi and aM2Psi
		f_Add_wMC_wMD(IsNotSharp, GP, inConfFP, Psi, aPsi, aM2Psi);
	}
	f_Add_wMC2(false, GP, 1.0, aPsi, M2aPsi);
}

/***************************************************************************/
// Set Input data
void cMT_MulSli_GPU::SetInputData(sInMSTEM &InMSTEM){
	freeMemory();

	MGP.gpu = InMSTEM.gpu;
	MGP.SimType = InMSTEM.SimType;	
	MGP.MulOrder = 2;	
	MGP.nConfFP = InMSTEM.nConfFP;		
	MGP.DimFP = InMSTEM.DimFP;	
	MGP.DistFP = 1;
	MGP.SeedFP = InMSTEM.SeedFP;
	MGP.PotPar = InMSTEM.PotPar;
	MGP.MEffect = InMSTEM.MEffect;
	MGP.STEffect = InMSTEM.STEffect;
	MGP.ZeroDefTyp = InMSTEM.ZeroDefTyp;
	MGP.ZeroDefPlane = InMSTEM.ZeroDefPlane;
	MGP.ApproxModel = InMSTEM.ApproxModel;
	MGP.ThkTyp = InMSTEM.ThicknessTyp;
	MGP.nThk = InMSTEM.nThickness;
	memcpy(MGP.Thk, InMSTEM.Thickness, MGP.nThk*cSizeofRD); // change
	MGP.BWL = (InMSTEM.BandwidthLimit==1)?true:false;
	MGP.Vrl = stVrl;
	MGP.E0 = InMSTEM.E0;	
	MGP.theta = InMSTEM.theta;	
	MGP.phi = InMSTEM.phi;
	MGP.lx = InMSTEM.lx;
	MGP.ly = InMSTEM.ly;
	MGP.dz = InMSTEM.dz;
	MGP.nx = InMSTEM.nx;
	MGP.ny = InMSTEM.ny;
	if(MGP.ApproxModel>1){
		MGP.MulOrder = 1;
		MGP.DimFP = MGP.DimFP - MGP.DimFP%10;
	}

	cudaSetDevice(MGP.gpu);

	f_sGP_Cal(MGP.nx, MGP.ny, MGP.lx, MGP.ly, MGP.dz, MGP.PBC_xy, MGP.BWL, GP);

	Lens.m = InMSTEM.MC_m;
	Lens.f = InMSTEM.MC_f;
	Lens.Cs3 = InMSTEM.MC_Cs3;
	Lens.Cs5 = InMSTEM.MC_Cs5;
	Lens.mfa2 = InMSTEM.MC_mfa2;
	Lens.afa2 = InMSTEM.MC_afa2;
	Lens.mfa3 = InMSTEM.MC_mfa3;
	Lens.afa3 = InMSTEM.MC_afa3;
	Lens.aobjl = InMSTEM.MC_aobjl;
	Lens.aobju = InMSTEM.MC_aobju;
	Lens.sf = InMSTEM.MC_sf;
	Lens.nsf = InMSTEM.MC_nsf;
	Lens.beta = InMSTEM.MC_beta;
	Lens.nbeta = InMSTEM.MC_nbeta;
	f_sLens_Cal(MGP.E0, GP, Lens);

	gxu = sin(MGP.theta)*cos(MGP.phi)/Lens.lambda;
	gyu = sin(MGP.theta)*sin(MGP.phi)/Lens.lambda;

	fPot = Lens.gamma*Lens.lambda/cos(MGP.theta);

	/*****************************************************************/
	CBED.x0 = InMSTEM.CBED_x0;
	CBED.y0 = InMSTEM.CBED_y0;
	CBED.Space = (InMSTEM.CBED_Space==1)?eSReal:eSReciprocal;

	// HRTEM.xx = 0;

	PED.nrot = InMSTEM.PED_nrot;
	PED.theta = InMSTEM.PED_theta;

	HCI.nrot = InMSTEM.HCI_nrot;
	HCI.theta = InMSTEM.HCI_theta;
	// HCI.xx = 0;

	// EWRS.xx = 0;

	// EWFS.xx = 0;
	/*****************************************************************/

	f_sACD_cudaMalloc(GP.nx, ExpRg_x);
	f_sACD_cudaMalloc(GP.ny, ExpRg_y);
	f_sACD_cudaMalloc(GP.nx, Prop_x);
	f_sACD_cudaMalloc(GP.ny, Prop_y);

	cudaMalloc((void**)&Trans, GP.nxy*cSizeofCD);

	cudaMalloc((void**)&Psi, GP.nxy*cSizeofCD);
	cudaMalloc((void**)&aPsi, GP.nxy*cSizeofCD);

	cudaMalloc((void**)&M2aPsi, GP.nxy*cSizeofRD);
	cudaMalloc((void**)&aM2Psi, GP.nxy*cSizeofRD);

	cufftPlan2d(&PlanPsi, GP.ny, GP.nx, CUFFT_Z2Z);

	MT_Potential_GPU.SetInputData(MGP, GP, InMSTEM.nAtomsM, InMSTEM.AtomsM);

	// Set parameters for the incident wave
	MT_IncidentWave_GPU.SetInputData(GP, Lens, PlanPsi);

	// Microscope parameters
	MT_MicroscopeEffects_GPU.SetInputData(GP, Lens, PlanPsi, Trans);
}

void cMT_MulSli_GPU::Cal_ExitWaveRS(sComplex &aPsih, double *&aM2Psih){
	int MEffect = 0, STEffect = 0;
	eSpace Space = eSReal;
	get_Imagefuncion(MGP.nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	GPU2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_ExitWaveFS(sComplex &aPsih, double *&aM2Psih){
	int MEffect = 0, STEffect = 0;
	eSpace Space = eSReciprocal;
	get_Imagefuncion(MGP.nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	GPU2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_ED(sComplex &aPsih, double *&aM2Psih){
	int MEffect = 0, STEffect = 0;
	eSpace Space = eSReciprocal;
	get_Imagefuncion(MGP.nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	GPU2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_HRTEM(sComplex &aPsih, double *&M2aPsih, double *&aM2Psih){
	eSpace Space = eSReal;
	get_Imagefuncion(MGP.nConfFP, Space, MGP.MEffect, MGP.STEffect, aPsi, M2aPsi, aM2Psi);
	GPU2CPU(aPsi, M2aPsi, aM2Psi, aPsih, M2aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_CBED(sComplex &aPsih, double *&aM2Psih){
	get_Imagefuncion(MGP.nConfFP, CBED.Space, CBED.x0, CBED.y0, aPsi, M2aPsi, aM2Psi);
	GPU2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_STEM(){
	int iThk = 0;

	if(!STEM.FastCal){
		cMT_Detector_GPU MT_Detector_GPU;
		MT_Detector_GPU.SetInputData(GP, STEM.nDet, STEM.DetCir);
		for (int ist=0; ist<STEM.nst; ist++){
			get_Imagefuncion(MGP.nConfFP, eSReciprocal, STEM.xst[ist], STEM.yst[ist], aPsi, M2aPsi, aM2Psi);
			MT_Detector_GPU.getDetectorIntensity(aM2Psi, M2aPsi, ist, STEM.ImSTEM[iThk].DetInt);
		}
	}else{
		switch(MGP.ApproxModel)
		{
			case 1:
				Cal_FAST_STEM_Wavefunction_MSA(MGP.nConfFP, STEM.ImSTEM[iThk].DetInt);
				break;
			case 2:
				Cal_FAST_STEM_Wavefunction_PA(MGP.nConfFP, STEM.ImSTEM[iThk].DetInt);
				break;
			case 3:
				Cal_FAST_STEM_Wavefunction_POA_WPOA(MGP.nConfFP, STEM.ImSTEM[iThk].DetInt);
				break;
			case 4:
				Cal_FAST_STEM_Wavefunction_POA_WPOA(MGP.nConfFP, STEM.ImSTEM[iThk].DetInt);
				break;
		}
	}
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
//		get_Imagefuncion(MGP.nConfFP, eSReciprocal, STEM.xs[ist], STEM.ys[ist], aPsi, M2aPsi, aM2Psi);
//		cudaMemcpy(aM2Psih, aM2Psi, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
//		cudaMemcpy(M2aPsih, M2aPsi, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
//		MT_Detector_CPU.getDetectorIntensity(aM2Psih, M2aPsih, ist, STEM.ImSTEM[iThk].DetInt);
//	}
//
//	delete [] aM2Psih; aM2Psih = 0;
//	delete [] M2aPsih; M2aPsih = 0;
//}