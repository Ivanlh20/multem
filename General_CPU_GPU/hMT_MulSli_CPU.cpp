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

#include "fftw3.h"
#include "hmathCPU.h"
#include "hConstTypes.h"
#include "hGeneral_CPU.h"
#include "hMT_General_CPU.h"
#include "hMT_InMulSli_CPU.h"
#include "hMT_MGP_CPU.h"
#include "hMT_Transmission_CPU.h"
#include "hMT_IncidentWave_CPU.h"
#include "hMT_MicroscopeEffects_CPU.h"
#include "hMT_STEM_CPU.h"
#include "hMT_MulSli_CPU.h"

void cMT_MulSli_CPU::freeMemory()
{
	if(IdCall==0) return;

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

	f_sACD_Free_CPU(ExpRg_x);
	f_sACD_Free_CPU(ExpRg_y);

	f_sACD_Free_CPU(Prop_x);
	f_sACD_Free_CPU(Prop_y);

	fftw_free(Psi);
	fftw_free(aPsi);

	fftw_free(M2aPsi);
	fftw_free(aM2Psi);

	fftw_destroy_plan(PlanForward);
	fftw_destroy_plan(PlanBackward);

	delete MT_Transmission_CPU; MT_Transmission_CPU = 0;
	delete MT_IncidentWave_CPU; MT_IncidentWave_CPU = 0;
	delete MT_MicroscopeEffects_CPU; MT_MicroscopeEffects_CPU = 0;
}

cMT_MulSli_CPU::cMT_MulSli_CPU()
{
	IdCall = 0;

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

	f_sACD_Init_CPU(ExpRg_x);
	f_sACD_Init_CPU(ExpRg_y);

	f_sACD_Init_CPU(Prop_x);
	f_sACD_Init_CPU(Prop_y);

	Psi = 0;
	aPsi = 0;

	M2aPsi = 0;
	aM2Psi = 0;

	PlanForward = 0;
	PlanBackward = 0;

	MT_Transmission_CPU = 0;
	MT_IncidentWave_CPU = 0;
	MT_MicroscopeEffects_CPU = 0;
}

cMT_MulSli_CPU::~cMT_MulSli_CPU()
{
	freeMemory();
	IdCall = 0;
}


void cMT_MulSli_CPU::Gather(fftw_complex *&Psid, sComplex &Psih)
{	
	// fft2shift 
	f_fft2Shift_MC_CPU(GP, Psid);

	f_fftw_complex_2_sComplex_CPU(GP, Psid, Psih);
}

void cMT_MulSli_CPU::Gather(fftw_complex *&Psid, double *&M2Psid, sComplex &Psih, double *&M2Psih)
{
	// fft2shift 
	f_fft2Shift_MC_MD_CPU(GP, Psid, M2Psid);

	f_fftw_complex_2_sComplex_CPU(GP, Psid, Psih);
	f_fftw_double_2_double_CPU(GP, M2Psid, M2Psih);
}

void cMT_MulSli_CPU::Gather(fftw_complex *&Psid, double *&M2Psi1d, double *&M2Psi2d, sComplex &Psih, double *&M2Psi1h, double *&M2Psi2h)
{
	// fft2shift 
	f_fft2Shift_MC_MD_CPU(GP, Psid, M2Psi1d, M2Psi2d);

	f_fftw_complex_2_sComplex_CPU(GP, Psid, Psih);
	f_fftw_double_2_double_CPU(GP, M2Psi1d, M2Psi1h);
	f_fftw_double_2_double_CPU(GP, M2Psi2d, M2Psi2h);
}

void cMT_MulSli_CPU::PhaseMul(double gxu, double gyu, fftw_complex *&Psi)
{
	if(MT_MGP_CPU.ShiftDP) return;

	if((gxu != 0.0)||(gyu != 0.0))
		f_PhaseMul_CPU(GP, gxu, gyu, ExpRg_x, ExpRg_y, Psi);
}

void cMT_MulSli_CPU::Propagate(eSpace Space, double gxu, double gyu, double z, fftw_complex *&Psi)
{
	f_Propagate_CPU(PlanForward, PlanBackward, GP, Space, gxu, gyu, Lens.lambda, z, Prop_x, Prop_y, Psi);
}

// Exit wave calculation: Space :1 Real and 2 Fourier
void cMT_MulSli_CPU::Cal_Wavefunction(eSpace Space, fftw_complex *&Psi)
{
	int iSlice=0, iSynCPU = 0;
	if(MT_MGP_CPU.ApproxModel<=2)
	{
		for(iSlice=0; iSlice<MT_Transmission_CPU->nSlice; iSlice++)
		{
			// get Transmission and Transmit
			MT_Transmission_CPU->Transmit(iSlice, Psi);
			// Propagate
			Propagate(eSReal, gxu, gyu, MT_Transmission_CPU->get_dz(iSlice), Psi);
		}
		// Last Transmission and Transmit
		if(MT_MGP_CPU.MulOrder==2)
		{
			MT_Transmission_CPU->Transmit(iSlice, Psi);
		}
		// Inclined ilumination
		PhaseMul(gxu, gyu, Psi);
		// Back propagation to our plane reference
		Propagate(Space, gxu, gyu, MT_Transmission_CPU->Thk.zb[0], Psi);
	}
	else
	{
		// get Transmission and Transmit
		MT_Transmission_CPU->Transmit(iSlice, Psi);
		// Inclined ilumination
		PhaseMul(gxu, gyu, Psi);

		if(Space==eSReciprocal)
		{
			fftw_execute_dft(PlanForward, Psi, Psi);
			f_Scale_MC_CPU(GP, GP.inxy, Psi);
		}
	}
}

// Get wave function for Plane wave ilumination
void cMT_MulSli_CPU::Image_Plane_Wave_Illumination(int nConfFP, eSpace Space, int MEffect, int STEffect, fftw_complex *&aPsi, double *&M2aPsi, double *&aM2Psi)
{
	int iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);

	f_Set_MC_MD_CPU(GP, 0.0, 0.0, aPsi, 0.0, aM2Psi);
	for(int iConf=iConf0; iConf<=nConfFP; iConf++)
	{
		// Move atoms
		MT_Transmission_CPU->MoveAtoms(iConf);
		// Plane wave ilumination
		MT_IncidentWave_CPU->Psi0(Psi);
		if(MEffect==0)
		{
			// Exit wave
			Cal_Wavefunction(Space, Psi);
			// Add Psi and M2aPsi to aPsi and aM2Psi
			f_Add_wMC_wMD_CPU(GP, inConfFP, Psi, aPsi, aM2Psi);
		}
		else
		{
			// Exit wave(g)
			Cal_Wavefunction(eSReciprocal, Psi);
			// Microscope effects
			MT_MicroscopeEffects_CPU->ApplyMEffects(MEffect, STEffect, Psi, M2aPsi);
			// Backward fft2
			fftw_execute_dft(PlanBackward, Psi, Psi);
			// Add Psi and M2aPsi to aPsi and aM2Psi
			f_Add_wMC_wMD_CPU(GP, inConfFP, Psi, M2aPsi, aPsi, aM2Psi);
		}
	}

	if(MEffect==0)
	{
		f_Set_wMC2_CPU(GP, 1.0, aPsi, M2aPsi);
	}
	else
	{
		// Forward fft2
		fftw_execute_dft(PlanForward, aPsi, aPsi);
		// Scale vector
		f_Scale_MC_CPU(GP, GP.inxy, aPsi);
		// Microscope effects
		MT_MicroscopeEffects_CPU->ApplyMEffects(MEffect, STEffect, aPsi, M2aPsi);
		// Backward fft2
		fftw_execute_dft(PlanBackward, aPsi, aPsi);
	}

}

// Get wave function for convergent beam ilumination
void cMT_MulSli_CPU::Image_Convergence_Wave_Illumination(int nConfFP, eSpace Space, double xi, double yi, fftw_complex *&aPsi, double *&M2aPsi, double *&aM2Psi)
{
	int iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);

	f_Set_MC_MD_CPU(GP, 0.0, 0.0, aPsi, 0.0, aM2Psi);
	for(int iConf=iConf0; iConf<=nConfFP; iConf++)
	{
		// Move atoms
		MT_Transmission_CPU->MoveAtoms(iConf);
		// Convergent wave ilumination
		MT_IncidentWave_CPU->Psi0(xi, yi, Psi);
		// Exit wave
		Cal_Wavefunction(Space, Psi);
		// Add Psi and M2aPsi to aPsi and aM2Psi
		f_Add_wMC_wMD_CPU(GP, inConfFP, Psi, aPsi, aM2Psi);
	}
	f_Set_wMC2_CPU(GP, 1.0, aPsi, M2aPsi);
}

// Set Input data
void cMT_MulSli_CPU::SetInputData(cMT_InMulSli_CPU &MT_InMulSli_CPU)
{
	freeMemory();
	IdCall++;

	MT_MGP_CPU.SetInputData(MT_InMulSli_CPU);

	f_sGP_SetInputData(&MT_MGP_CPU, GP);
	f_sLens_SetInputData(MT_InMulSli_CPU, GP, Lens);

	gxu = sin(MT_MGP_CPU.theta)*cos(MT_MGP_CPU.phi)/Lens.lambda;
	gyu = sin(MT_MGP_CPU.theta)*sin(MT_MGP_CPU.phi)/Lens.lambda;

	f_sACD_Malloc_CPU(GP.nx, ExpRg_x);
	f_sACD_Malloc_CPU(GP.ny, ExpRg_y);
	f_sACD_Malloc_CPU(GP.nx, Prop_x);
	f_sACD_Malloc_CPU(GP.ny, Prop_y);

	Psi = (fftw_complex*)fftw_malloc(GP.nxy*cSizeofCD);
	aPsi = (fftw_complex*)fftw_malloc(GP.nxy*cSizeofCD);

	M2aPsi = new double[GP.nxy];
	aM2Psi = new double[GP.nxy];

	PlanForward = fftw_plan_dft_2d(GP.ny, GP.nx, Psi, Psi, FFTW_FORWARD, FFTW_ESTIMATE);
	PlanBackward = fftw_plan_dft_2d(GP.ny, GP.nx, Psi, Psi, FFTW_BACKWARD, FFTW_ESTIMATE);

	// Transmission parameters
	MT_Transmission_CPU = new cMT_Transmission_CPU ;
	MT_Transmission_CPU->SetInputData(&MT_MGP_CPU, PlanForward, PlanBackward, MT_InMulSli_CPU.nAtomsM, MT_InMulSli_CPU.AtomsM);

	// Microscope parameters
	MT_MicroscopeEffects_CPU = new cMT_MicroscopeEffects_CPU;
	MT_MicroscopeEffects_CPU->SetInputData(GP, Lens, PlanForward, PlanBackward, MT_Transmission_CPU->Trans0);

	// Incident wave parameters
	MT_IncidentWave_CPU = new cMT_IncidentWave_CPU;
	MT_IncidentWave_CPU->SetInputData(&MT_MGP_CPU, Lens, PlanForward, PlanBackward);

	/***************************************************************************/
	/***************************************************************************/

	switch(MT_InMulSli_CPU.SimType)
	{
		case 11:		// STEM
			STEM = new cMT_STEM_CPU;
			STEM->SetInputData(&MT_InMulSli_CPU, &MT_MGP_CPU, MT_Transmission_CPU->Thk.n);
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

void cMT_MulSli_CPU::Cal_ExitWaveRS(sComplex &aPsih, double *&aM2Psih)
{
	int MEffect = 0, STEffect = 0;
	eSpace Space = eSReal;
	Image_Plane_Wave_Illumination(MT_MGP_CPU.nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	Gather(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_CPU::Cal_ExitWaveFS(sComplex &aPsih, double *&aM2Psih)
{
	int MEffect = 0, STEffect = 0;
	eSpace Space = eSReciprocal;
	Image_Plane_Wave_Illumination(MT_MGP_CPU.nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	Gather(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_CPU::Cal_ED(sComplex &aPsih, double *&aM2Psih)
{
	int MEffect = 0, STEffect = 0;
	eSpace Space = eSReciprocal;
	Image_Plane_Wave_Illumination(MT_MGP_CPU.nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	Gather(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_CPU::Cal_HRTEM(sComplex &aPsih, double *&M2aPsih, double *&aM2Psih)
{
	eSpace Space = eSReal;
	Image_Plane_Wave_Illumination(MT_MGP_CPU.nConfFP, Space, MT_MGP_CPU.MEffect, MT_MGP_CPU.STEffect, aPsi, M2aPsi, aM2Psi);
	Gather(aPsi, M2aPsi, aM2Psi, aPsih, M2aPsih, aM2Psih);
}

void cMT_MulSli_CPU::Cal_CBED(sComplex &aPsih, double *&aM2Psih)
{
	Image_Convergence_Wave_Illumination(MT_MGP_CPU.nConfFP, eSReciprocal, CBED.x0, CBED.y0, aPsi, M2aPsi, aM2Psi);
	Gather(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_CPU::Cal_CBEI(sComplex &aPsih, double *&aM2Psih)
{
	Image_Convergence_Wave_Illumination(MT_MGP_CPU.nConfFP, eSReal, CBEI.x0, CBEI.y0, aPsi, M2aPsi, aM2Psi);
	Gather(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_CPU::Cal_PED(sComplex &aPsih, double *&aM2Psih)
{
	int iSlice=0, iSynCPU = 0;
	int irot, iThk = 0;
	int nConfFP = MT_MGP_CPU.nConfFP, iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);
	double gu, phi;
	double w = inConfFP/double(PED.nrot);

	f_Set_MC_MD_CPU(GP, 0.0, 0.0, aPsi, 0.0, aM2Psi);
	phi = 2.0*cPi/double(PED.nrot); gu = sin(PED.theta)/Lens.lambda;
	for(int iConf=iConf0; iConf<=nConfFP; iConf++)
	{
		// Move atoms
		MT_Transmission_CPU->MoveAtoms(iConf);
		for(irot=0; irot<PED.nrot; irot++)
		{
			gxu = -gu*cos(irot*phi); gyu = -gu*sin(irot*phi);
			// Plane wave ilumination
			MT_IncidentWave_CPU->Psi0(Psi);
			// Exit wave(g)
			Cal_Wavefunction(eSReciprocal, Psi);
			// Add Psi and M2aPsi to aPsi and aM2Psi
			f_Add_wMC_wMD_CPU(GP, w, Psi, aPsi, aM2Psi);
		}
	}
	Gather(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMT_MulSli_CPU::CAL_HCI(sComplex &aPsih, double *&M2aPsih, double *&aM2Psih)
{
	int iSlice=0, iSynCPU = 0;
	int irot, iThk = 0;
	int nConfFP = MT_MGP_CPU.nConfFP, iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);
	double gu, phi;
	double w = inConfFP/double(HCI.nrot);

	f_Set_MC_MD_CPU(GP, 0.0, 0.0, aPsi, 0.0, aM2Psi);
	phi = 2.0*cPi/double(HCI.nrot); gu = sin(HCI.theta)/Lens.lambda;
	for(int iConf=iConf0; iConf<=nConfFP; iConf++)
	{
		// Move atoms
		MT_Transmission_CPU->MoveAtoms(iConf);
		for(irot=0; irot<HCI.nrot; irot++)
		{
			gxu = -gu*cos(irot*phi); gyu = -gu*sin(irot*phi);
			// Plane wave ilumination
			MT_IncidentWave_CPU->Psi0(Psi);
			// Exit wave(g)
			Cal_Wavefunction(eSReciprocal, Psi);
			// Microscope effects
			MT_MicroscopeEffects_CPU->ApplyMEffects(MT_MGP_CPU.MEffect, MT_MGP_CPU.STEffect, Psi, M2aPsi);
			// Backward fft2
			fftw_execute_dft(PlanBackward, Psi, Psi);
			// Add Psi and M2aPsi to aPsi and aM2Psi
			f_Add_wMC_wMD_CPU(GP, w, Psi, M2aPsi, aPsi, aM2Psi);
		}
	}
	// Forward fft2
	fftw_execute_dft(PlanForward, aPsi, aPsi);
	// Scale vector
	f_Scale_MC_CPU(GP, GP.inxy, aPsi);
	// Microscope effects
	MT_MicroscopeEffects_CPU->ApplyMEffects(MT_MGP_CPU.MEffect, MT_MGP_CPU.STEffect, aPsi, M2aPsi);
	// Backward fft2
	fftw_execute_dft(PlanBackward, aPsi, aPsi);

	Gather(aPsi, M2aPsi, aM2Psi, aPsih, M2aPsih, aM2Psih);
}

void cMT_MulSli_CPU::Cal_STEM()
{
	int ist, iSlice, iThk = 0, iSynCPU = 0;
	int nConfFP = MT_MGP_CPU.nConfFP;
	int iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);

	STEM->InitImSTEM();
	for(int iConf=iConf0; iConf<=nConfFP; iConf++)
	{
		// Move atoms
		MT_Transmission_CPU->MoveAtoms(iConf);
		for(ist=0; ist<STEM->nst; ist++)
		{
			// Convergent wave ilumination
			MT_IncidentWave_CPU->Psi0(STEM->xst[ist], STEM->yst[ist], Psi);
			// Exit wave
			for(iSlice=0; iSlice<MT_Transmission_CPU->nSlice; iSlice++)
			{
				// get Transmission and Transmit
				MT_Transmission_CPU->Transmit(iSlice, Psi);
				// Propagate
				Propagate(eSReciprocal, gxu, gyu, MT_Transmission_CPU->get_dz(iSlice), Psi);
				/***********************************************************************************/
				/***********************************************************************************/
				iThk = MT_Transmission_CPU->IsInThk(iSlice);
				if(iThk>-1) // Detector integration
				{
					//STEM->MT_Detector_CPU->getDetectorIntensity(inConfFP, Psi, ist, STEM->ImSTEM[iThk].DetInt, true);
				}
				/***********************************************************************************/
				/***********************************************************************************/
				if(iSlice<MT_Transmission_CPU->nSlice-1)
				{
					fftw_execute_dft(PlanBackward, Psi, Psi);
				}
			}
		}
	}
}

//void cMT_MulSli_CPU::Cal_STEM()
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
//			f_Set_MC_MD_CPU(GP, 0.0, 0.0, aPsi, 0.0, aM2Psi);
//			for(int iConf=iConf0; iConf<=nConfFP; iConf++)
//			{
//				// Move atoms
//				MT_Transmission_CPU->MoveAtoms(iConf);	
//				// Convergent wave ilumination
//				MT_IncidentWave_CPU->Psi0(STEM->xst[ist], STEM->yst[ist], Psi);
//				// Exit wave
//				Cal_Wavefunction(eSReciprocal, Psi);
//				// Add Psi and M2aPsi to aPsi and aM2Psi
//				f_Add_wMC_wMD_CPU(GP, inConfFP, Psi, aPsi, aM2Psi);	
//			}
//			f_set_wMC2_CPU(GP, 1.0, aPsi, M2aPsi);
//			// Detector integration
//			STEM->MT_Detector_CPU->getDetectorIntensity(aM2Psi, M2aPsi, ist, STEM->ImSTEM[iThk].DetInt);
//		}
//	}
//	else
//	{
//		for(int iConf=iConf0; iConf<=nConfFP; iConf++)
//		{
//			// Move atoms
//			MT_Transmission_CPU->MoveAtoms(iConf);
//			for(ist=0; ist<STEM->nst; ist++)
//			{
//				// Convergent wave ilumination
//				MT_IncidentWave_CPU->Psi0(STEM->xst[ist], STEM->yst[ist], Psi);
//				// Exit wave
//				Cal_Wavefunction(eSReciprocal, Psi);
//				// Add Psi to aM2Psi
//				f_Set_wMC2_CPU(GP, inConfFP, Psi, aM2Psi);
//				// Detector integration
//				STEM->MT_Detector_CPU->getDetectorIntensity(aM2Psi, ist, STEM->ImSTEM[iThk].DetInt, true);
//			}
//		}
//	}
//}
//
//void cMT_MulSli_CPU::Cal_STEM()
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
