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

void cMT_MulSli_GPU::freeMemory()
{
	if(IdCall==0) return;

	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	cSynCPU = ccSynCPU;

	f_sGP_Init(GP);
	f_sLens_Init(Lens);
	
	gxu = 0;
	gyu = 0;

	/****************************************************************/
	delete STEM; STEM = 0;

	CBED.x0 = 0;
	CBED.y0 = 0;

	CBEI.x0 = 0;
	CBEI.y0 = 0;

	HRTEM.xx = 0;

	PED.nrot = 0;
	PED.theta = 0;

	HCI.nrot = 0;
	HCI.theta = 0;

	EWRS.xx = 0;

	EWFS.xx = 0;
	/****************************************************************/

	f_sACD_Free_GPU(ExpRg_x);
	f_sACD_Free_GPU(ExpRg_y);

	f_sACD_Free_GPU(Prop_x);
	f_sACD_Free_GPU(Prop_y);

	delete [] MC_h;

	cudaFreen(Psi);
	cudaFreen(aPsi);

	cudaFreen(M2aPsi);
	cudaFreen(aM2Psi);

	cufftDestroyn(PlanPsi);

	delete MT_Transmission_GPU; MT_Transmission_GPU = 0;
	delete MT_IncidentWave_GPU; MT_IncidentWave_GPU = 0;
	delete MT_MicroscopeEffects_GPU; MT_MicroscopeEffects_GPU = 0;
}

cMT_MulSli_GPU::cMT_MulSli_GPU()
{
	IdCall = 0;

	cSynCPU = ccSynCPU;

	f_sGP_Init(GP);
	f_sLens_Init(Lens);
	
	gxu = 0;		
	gyu = 0;		

	/****************************************************************/
	STEM = 0;

	CBED.x0 = 0;
	CBED.y0 = 0;

	CBEI.x0 = 0;
	CBEI.y0 = 0;

	HRTEM.xx = 0;

	PED.nrot = 0;
	PED.theta = 0;

	HCI.nrot = 0;
	HCI.theta = 0;

	EWRS.xx = 0;

	EWFS.xx = 0;
	/****************************************************************/

	f_sACD_Init_GPU(ExpRg_x);
	f_sACD_Init_GPU(ExpRg_y);

	f_sACD_Init_GPU(Prop_x);
	f_sACD_Init_GPU(Prop_y);

	MC_h = 0;

	Psi = 0;
	aPsi = 0;

	M2aPsi = 0;
	aM2Psi = 0;

	PlanPsi = 0;

	MT_Transmission_GPU = 0;
	MT_IncidentWave_GPU = 0;
	MT_MicroscopeEffects_GPU = 0;
}

cMT_MulSli_GPU::~cMT_MulSli_GPU()
{
	freeMemory();
	IdCall = 0;
	cudaDeviceReset();
}

// From Device To Host
void cMT_MulSli_GPU::Gather(double2 *&MC_d_i, sComplex &MC_h_o)
{	
	// fft2shift 
	f_fft2Shift_MC_GPU(GP, MC_d_i);
	/*********************copy data to host************************/
	f_cuDoubleComplex_2_sComplex_GPU(GP, MC_d_i, MC_h, MC_h_o);
}

// From Device To Host
void cMT_MulSli_GPU::Gather(double2 *&MC_d_i, double *&MD_d_i, sComplex &MC_h_o, double *&MD_h_o)
{
	// fft2shift 
	f_fft2Shift_MC_MD_GPU(GP, MC_d_i, MD_d_i);
	/*********************copy data to host************************/
	f_cuDoubleComplex_2_sComplex_GPU(GP, MC_d_i, MC_h, MC_h_o);
	f_cuDouble_2_double_GPU(GP, MD_d_i, MD_h_o);
}

// From Device To Host
void cMT_MulSli_GPU::Gather(double2 *&MC_d_i, double *&MD1_d_i, double *&MD2_d_i, sComplex &MC_h_o, double *&MD1_h_o, double *&MD2_h_o)
{
	// fft2shift 
	f_fft2Shift_MC_MD_GPU(GP, MC_d_i, MD1_d_i, MD2_d_i);
	/*********************copy data to host************************/
	f_cuDoubleComplex_2_sComplex_GPU(GP, MC_d_i, MC_h, MC_h_o);
	f_cuDouble_2_double_GPU(GP, MD1_d_i, MD1_h_o);
	f_cuDouble_2_double_GPU(GP, MD2_d_i, MD2_h_o);
}

void cMT_MulSli_GPU::PhaseMul(double gxu, double gyu, double2 *&Psi)
{
	if(MT_MGP_CPU.ShiftDP) return;

	if((gxu != 0.0)||(gyu != 0.0))
		f_PhaseMul_GPU(GP, gxu, gyu, ExpRg_x, ExpRg_y, Psi);
}

void cMT_MulSli_GPU::Propagate(eSpace Space, double gxu, double gyu, double z, double2 *&Psi)
{
	f_Propagate_GPU(PlanPsi, GP, Space, gxu, gyu, Lens.lambda, z, Prop_x, Prop_y, Psi);
}

// Exit wave calculation: Space :1 Real and 2 Fourier
void cMT_MulSli_GPU::Cal_Wavefunction(eSpace Space, double2 *&Psi)
{
	int iSlice=0, iSynCPU = 0;
	if(MT_MGP_CPU.ApproxModel<=2)
	{
		for(iSlice=0; iSlice<MT_Transmission_GPU->nSlice; iSlice++)
		{
			// get Transmission and Transmit
			MT_Transmission_GPU->Transmit(iSlice, Psi);
			// Propagate
			Propagate(eSReal, gxu, gyu, MT_Transmission_GPU->get_dz(iSlice), Psi);
			// GPU Synchronize
			f_GPU_Sync_CPU(iSynCPU, cSynCPU); 
		}
		// Inclined ilumination
		PhaseMul(gxu, gyu, Psi);
		// Back propagation to our plane reference
		Propagate(Space, gxu, gyu, MT_Transmission_GPU->Thk.zb[0], Psi);
	}
	else
	{
		// get Transmission and Transmit
		MT_Transmission_GPU->Transmit(iSlice, Psi);
		// Inclined ilumination
		PhaseMul(gxu, gyu, Psi);

		if(Space==eSReciprocal)
		{
			cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
			f_Scale_MC_GPU(GP, GP.inxy, Psi);
		}
	}
}

// Get wave function for Plane wave ilumination
void cMT_MulSli_GPU::Image_Plane_Wave_Illumination(int nConfFP, eSpace Space, int MEffect, int STEffect, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi)
{
	int iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);

	f_Set_MC_MD_GPU(GP, 0.0, 0.0, aPsi, 0.0, aM2Psi);
	for(int iConf=iConf0; iConf<=nConfFP; iConf++)
	{
		// Move atoms
		MT_Transmission_GPU->MoveAtoms(iConf);
		// Plane wave ilumination
		MT_IncidentWave_GPU->Psi0(Psi);
		if(MEffect==0)
		{
			// Exit wave
			Cal_Wavefunction(Space, Psi);
			// Add Psi and M2aPsi to aPsi and aM2Psi
			f_Add_wMC_wMD_GPU(GP, inConfFP, Psi, aPsi, aM2Psi);
		}
		else
		{
			// Exit wave(g)
			Cal_Wavefunction(eSReciprocal, Psi);
			// Microscope effects
			MT_MicroscopeEffects_GPU->ApplyMEffects(MEffect, STEffect, Psi, M2aPsi);
			// Backward fft2
			cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_INVERSE);
			// Add Psi and M2aPsi to aPsi and aM2Psi
			f_Add_wMC_wMD_GPU(GP, inConfFP, Psi, M2aPsi, aPsi, aM2Psi);
		}
	}

	if(MEffect==0)
	{
		f_Set_wMC2_GPU(GP, 1.0, aPsi, M2aPsi);
	}
	else
	{
		// Forward fft2
		cufftExecZ2Z(PlanPsi, aPsi, aPsi, CUFFT_FORWARD);
		// Scale vector
		f_Scale_MC_GPU(GP, GP.inxy, aPsi);
		// Microscope effects
		MT_MicroscopeEffects_GPU->ApplyMEffects(MEffect, STEffect, aPsi, M2aPsi);
		// Backward fft2
		cufftExecZ2Z(PlanPsi, aPsi, aPsi, CUFFT_INVERSE);
	}

}

// Get wave function for convergent beam ilumination
void cMT_MulSli_GPU::Image_Convergence_Wave_Illumination(int nConfFP, eSpace Space, double xi, double yi, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi)
{
	int iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);

	f_Set_MC_MD_GPU(GP, 0.0, 0.0, aPsi, 0.0, aM2Psi);
	for(int iConf=iConf0; iConf<=nConfFP; iConf++)
	{
		// Move atoms
		MT_Transmission_GPU->MoveAtoms(iConf);
		// Convergent wave ilumination
		MT_IncidentWave_GPU->Psi0(xi, yi, Psi);
		// Exit wave
		Cal_Wavefunction(Space, Psi);
		// Add Psi and M2aPsi to aPsi and aM2Psi
		f_Add_wMC_wMD_GPU(GP, inConfFP, Psi, aPsi, aM2Psi);
	}
	f_Set_wMC2_GPU(GP, 1.0, aPsi, M2aPsi);
}

// Set Input data
void cMT_MulSli_GPU::SetInputData(cMT_InMulSli_CPU &MT_InMulSli_CPU)
{
	freeMemory();
	IdCall++;

	MT_MGP_CPU.SetInputData(MT_InMulSli_CPU);

	cudaSetDevice(MT_MGP_CPU.GPU_Device);

	f_sGP_SetInputData(&MT_MGP_CPU, GP);
	f_sLens_SetInputData(MT_InMulSli_CPU, GP, Lens);

	gxu = sin(MT_MGP_CPU.theta)*cos(MT_MGP_CPU.phi)/Lens.lambda;
	gyu = sin(MT_MGP_CPU.theta)*sin(MT_MGP_CPU.phi)/Lens.lambda;

	f_sACD_Malloc_GPU(GP.nx, ExpRg_x);
	f_sACD_Malloc_GPU(GP.ny, ExpRg_y);
	f_sACD_Malloc_GPU(GP.nx, Prop_x);
	f_sACD_Malloc_GPU(GP.ny, Prop_y);

	MC_h = new double2[GP.nxy];

	cudaMalloc((void**)&Psi, GP.nxy*cSizeofCD);
	cudaMalloc((void**)&aPsi, GP.nxy*cSizeofCD);

	cudaMalloc((void**)&M2aPsi, GP.nxy*cSizeofRD);
	cudaMalloc((void**)&aM2Psi, GP.nxy*cSizeofRD);

	cufftPlan2d(&PlanPsi, GP.ny, GP.nx, CUFFT_Z2Z);

	// Transmission parameters
	MT_Transmission_GPU = new cMT_Transmission_GPU ;
	MT_Transmission_GPU->SetInputData(&MT_MGP_CPU, PlanPsi, MT_InMulSli_CPU.nAtomsM, MT_InMulSli_CPU.AtomsM);

	// Microscope parameters
	MT_MicroscopeEffects_GPU = new cMT_MicroscopeEffects_GPU;
	MT_MicroscopeEffects_GPU->SetInputData(GP, Lens, PlanPsi, MT_Transmission_GPU->Trans0);

	// Incident wave parameters
	MT_IncidentWave_GPU = new cMT_IncidentWave_GPU;
	MT_IncidentWave_GPU->SetInputData(&MT_MGP_CPU, Lens, PlanPsi, MC_h);

	/***************************************************************************/
	/***************************************************************************/

	switch(MT_InMulSli_CPU.SimType)
	{
		case 11:		// STEM
			STEM = new cMT_STEM_GPU;
			STEM->SetInputData(&MT_InMulSli_CPU, &MT_MGP_CPU, MT_Transmission_GPU->Thk.n);
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

void cMT_MulSli_GPU::Cal_ExitWaveRS(sComplex &aPsih, double *&aM2Psih)
{
	int MEffect = 0, STEffect = 0;
	eSpace Space = eSReal;
	Image_Plane_Wave_Illumination(MT_MGP_CPU.nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	Gather(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_ExitWaveFS(sComplex &aPsih, double *&aM2Psih)
{
	int MEffect = 0, STEffect = 0;
	eSpace Space = eSReciprocal;
	Image_Plane_Wave_Illumination(MT_MGP_CPU.nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	Gather(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_ED(sComplex &aPsih, double *&aM2Psih)
{
	int MEffect = 0, STEffect = 0;
	eSpace Space = eSReciprocal;
	Image_Plane_Wave_Illumination(MT_MGP_CPU.nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	Gather(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_HRTEM(sComplex &aPsih, double *&M2aPsih, double *&aM2Psih)
{
	eSpace Space = eSReal;
	Image_Plane_Wave_Illumination(MT_MGP_CPU.nConfFP, Space, MT_MGP_CPU.MEffect, MT_MGP_CPU.STEffect, aPsi, M2aPsi, aM2Psi);
	Gather(aPsi, M2aPsi, aM2Psi, aPsih, M2aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_CBED(sComplex &aPsih, double *&aM2Psih)
{
	Image_Convergence_Wave_Illumination(MT_MGP_CPU.nConfFP, eSReciprocal, CBED.x0, CBED.y0, aPsi, M2aPsi, aM2Psi);
	Gather(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_CBEI(sComplex &aPsih, double *&aM2Psih)
{
	Image_Convergence_Wave_Illumination(MT_MGP_CPU.nConfFP, eSReal, CBEI.x0, CBEI.y0, aPsi, M2aPsi, aM2Psi);
	Gather(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_PED(sComplex &aPsih, double *&aM2Psih)
{
	int iSlice=0, iSynCPU = 0;
	int irot, iThk = 0;
	int nConfFP = MT_MGP_CPU.nConfFP, iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);
	double gu, phi, ecos, esin;
	double w = inConfFP/double(PED.nrot);

	f_Set_MC_MD_GPU(GP, 0.0, 0.0, aPsi, 0.0, aM2Psi);
	phi = 2.0*cPi/double(PED.nrot); gu = sin(PED.theta)/Lens.lambda;
	for(int iConf=iConf0; iConf<=nConfFP; iConf++)
	{
		// Move atoms
		MT_Transmission_GPU->MoveAtoms(iConf);
		for(irot=0; irot<PED.nrot; irot++)
		{
			sincos(irot*phi, &esin, &ecos); 
			gxu = -gu*ecos; gyu = -gu*esin;
			// Plane wave ilumination
			MT_IncidentWave_GPU->Psi0(Psi);
			// Exit wave(g)
			Cal_Wavefunction(eSReciprocal, Psi);
			// Add Psi and M2aPsi to aPsi and aM2Psi
			f_Add_wMC_wMD_GPU(GP, w, Psi, aPsi, aM2Psi);
		}
	}
	Gather(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_GPU::CAL_HCI(sComplex &aPsih, double *&M2aPsih, double *&aM2Psih)
{
	int iSlice=0, iSynCPU = 0;
	int irot, iThk = 0;
	int nConfFP = MT_MGP_CPU.nConfFP, iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);
	double gu, phi, ecos, esin;
	double w = inConfFP/double(HCI.nrot);

	f_Set_MC_MD_GPU(GP, 0.0, 0.0, aPsi, 0.0, aM2Psi);
	phi = 2.0*cPi/double(HCI.nrot); gu = sin(HCI.theta)/Lens.lambda;
	for(int iConf=iConf0; iConf<=nConfFP; iConf++)
	{
		// Move atoms
		MT_Transmission_GPU->MoveAtoms(iConf);
		for(irot=0; irot<HCI.nrot; irot++)
		{
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
			f_Add_wMC_wMD_GPU(GP, w, Psi, M2aPsi, aPsi, aM2Psi);
		}
	}
	// Forward fft2
	cufftExecZ2Z(PlanPsi, aPsi, aPsi, CUFFT_FORWARD);
	// Scale vector
	f_Scale_MC_GPU(GP, GP.inxy, aPsi);
	// Microscope effects
	MT_MicroscopeEffects_GPU->ApplyMEffects(MT_MGP_CPU.MEffect, MT_MGP_CPU.STEffect, aPsi, M2aPsi);
	// Backward fft2
	cufftExecZ2Z(PlanPsi, aPsi, aPsi, CUFFT_INVERSE);

	Gather(aPsi, M2aPsi, aM2Psi, aPsih, M2aPsih, aM2Psih);
}

void cMT_MulSli_GPU::Cal_STEM()
{
	int ist, iSlice, iThk = 0, iSynCPU = 0;
	int nConfFP = MT_MGP_CPU.nConfFP;
	int iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);

	STEM->InitImSTEM();
	for(int iConf=iConf0; iConf<=nConfFP; iConf++)
	{
		// Move atoms
		MT_Transmission_GPU->MoveAtoms(iConf);
		for(ist=0; ist<STEM->nst; ist++)
		{
			// Convergent wave ilumination
			MT_IncidentWave_GPU->Psi0(STEM->xst[ist], STEM->yst[ist], Psi);
			// Exit wave
			for(iSlice=0; iSlice<MT_Transmission_GPU->nSlice; iSlice++)
			{
				// get Transmission and Transmit
				MT_Transmission_GPU->Transmit(iSlice, Psi);
				// Propagate
				Propagate(eSReciprocal, gxu, gyu, MT_Transmission_GPU->get_dz(iSlice), Psi);
				/***********************************************************************************/
				/***********************************************************************************/
				iThk = MT_Transmission_GPU->IsInThk(iSlice);
				if(iThk>-1) // Detector integration
				{
					STEM->MT_Detector_GPU->getDetectorIntensity(inConfFP, Psi, ist, STEM->ImSTEM[iThk].DetInt, true);
				}
				/***********************************************************************************/
				/***********************************************************************************/
				if(iSlice<MT_Transmission_GPU->nSlice-1)
				{
					cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_INVERSE);
				}
			}
		}
	}
}

//void cMT_MulSli_GPU::Cal_STEM()
//{
//	int ist, iThk = 0;
//	int nConfFP = MT_MGP_CPU.nConfFP;
//	int iConf0 = (nConfFP==0)?0:1;
//	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);
//
//	STEM->InitImSTEM();
//	if(MT_MGP_CPU.FastCal==1)
//	{
//		for(ist=0; ist<STEM->nst; ist++)
//		{
//			f_Set_MC_MD_GPU(GP, 0.0, 0.0, aPsi, 0.0, aM2Psi);
//			for(int iConf=iConf0; iConf<=nConfFP; iConf++)
//			{
//				// Move atoms
//				MT_Transmission_GPU->MoveAtoms(iConf);	
//				// Convergent wave ilumination
//				MT_IncidentWave_GPU->Psi0(STEM->xst[ist], STEM->yst[ist], Psi);
//				// Exit wave
//				Cal_Wavefunction(eSReciprocal, Psi);
//				// Add Psi and M2aPsi to aPsi and aM2Psi
//				f_Add_wMC_wMD_GPU(GP, inConfFP, Psi, aPsi, aM2Psi);	
//			}
//			f_set_wMC2_GPU(GP, 1.0, aPsi, M2aPsi);
//			// Detector integration
//			STEM->MT_Detector_GPU->getDetectorIntensity(aM2Psi, M2aPsi, ist, STEM->ImSTEM[iThk].DetInt);
//		}
//	}
//	else
//	{
//		for(int iConf=iConf0; iConf<=nConfFP; iConf++)
//		{
//			// Move atoms
//			MT_Transmission_GPU->MoveAtoms(iConf);
//			for(ist=0; ist<STEM->nst; ist++)
//			{
//				// Convergent wave ilumination
//				MT_IncidentWave_GPU->Psi0(STEM->xst[ist], STEM->yst[ist], Psi);
//				// Exit wave
//				Cal_Wavefunction(eSReciprocal, Psi);
//				// Add Psi to aM2Psi
//				f_Set_wMC2_GPU(GP, inConfFP, Psi, aM2Psi);
//				// Detector integration
//				STEM->MT_Detector_GPU->getDetectorIntensity(aM2Psi, ist, STEM->ImSTEM[iThk].DetInt, true);
//			}
//		}
//	}
//}
//
//void cMT_MulSli_GPU::Cal_STEM()
//{
//	int ist, iThk = 0;
//	double *M2aPsih, *aM2Psih;
//	aM2Psih = new double[GP.nxy];
//	M2aPsih = new double[GP.nxy];
//	cMT_Detector_CPU MT_Detector_CPU;
//	MT_Detector_CPU.SetInputData(GP, STEM.nDet, STEM.DetCir);
//
//	for(ist=0; ist<STEM.nst; ist++)
//	{
//		get_Imagefuncion(MT_MGP_CPU.nConfFP, eSReciprocal, STEM.xs[ist], STEM.ys[ist], aPsi, M2aPsi, aM2Psi);
//		cudaMemcpy(aM2Psih, aM2Psi, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
//		cudaMemcpy(M2aPsih, M2aPsi, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
//		MT_Detector_CPU.getDetectorIntensity(aM2Psih, M2aPsih, ist, STEM.ImSTEM[iThk].DetInt);
//	}
//
//	delete [] aM2Psih; aM2Psih = 0;
//	delete [] M2aPsih; M2aPsih = 0;
//}
