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

#ifndef hMT_MulSli_CPU_H
#define hMT_MulSli_CPU_H

#include "fftw3.h"
#include "hConstTypes.h"
#include "hMT_General_CPU.h"
#include "hMT_InMulSli_CPU.h"
#include "hMT_MGP_CPU.h"
#include "hMT_Transmission_CPU.h"
#include "hMT_IncidentWave_CPU.h"
#include "hMT_MicroscopeEffects_CPU.h"
#include "hMT_STEM_CPU.h"
#include "hMT_MulSli_CPU.h"

class cMT_STEM_CPU;

class cMT_MulSli_CPU{
	private:
		int IdCall;
	public:
		sLens Lens;													// Aberration parameters
		sGP GP;														// Grid variables

		double gxu;													// incident beam x-tilt in reciprocal units
		double gyu;													// incident beam y-tilt in reciprocal units

		sACD ExpRg_x;												// k_PhaseMul x
		sACD ExpRg_y;												// k_PhaseMul y
		sACD Prop_x;												// Propagator x
		sACD Prop_y;												// Propagator y

		fftw_complex *Psi;											// Wave function
		fftw_complex *aPsi;											// Wave function - temporal

		double *M2aPsi;												// Squared Wave function
		double *aM2Psi;												// Squared Wave function - temporal

		fftw_plan PlanForward;										// Fourier transform's plan Forward
		fftw_plan PlanBackward;										// Fourier transform's plan Backward

		cMT_IncidentWave_CPU *MT_IncidentWave_CPU;					// Incident wave;
		cMT_Transmission_CPU *MT_Transmission_CPU;					// Transmission function
		cMT_MicroscopeEffects_CPU *MT_MicroscopeEffects_CPU;		// Microscope effects

		void PhaseMul(double gxu, double gyu, fftw_complex *&Psi);
		void Propagate(eSpace Space, double gxu, double gyu, double z, fftw_complex *&Psi);
		void Cal_Wavefunction(eSpace Space, fftw_complex *&Psi);

		void Image_Plane_Wave_Illumination(int nConfFP, eSpace Space, int MEffect, int STEffect, fftw_complex *&aPsi, double *&M2aPsi, double *&aM2Psi);
		void Image_Convergence_Wave_Illumination(int nConfFP, eSpace Space, double xi, double yi, fftw_complex *&aPsi, double *&M2aPsi, double *&aM2Psi);

		void GPU2CPU(fftw_complex *&Psid, sComplex &Psih);
		void GPU2CPU(fftw_complex *&Psid, double *&M2Psid, sComplex &Psih, double *&M2Psih);
		void GPU2CPU(fftw_complex *&Psid, double *&M2Psi1d, double *&M2Psi2d, sComplex &Psih, double *&M2Psi1h, double *&M2Psi2h);

		cMT_MGP_CPU MT_MGP_CPU;
		cMT_STEM_CPU *STEM;
		sCBED CBED;
		sCBEI CBEI;
		sHRTEM HRTEM;
		sPED PED;
		sHCI HCI;
		sEWRS EWRS;
		sEWFS EWFS;	
	
		void freeMemory();
		cMT_MulSli_CPU();
		~cMT_MulSli_CPU();

		void SetInputData(cMT_InMulSli_CPU &MT_InMulSli_CPU);
		void Cal_STEM();
		void Cal_ISTEM();

		void Cal_CBED(sComplex &aPsih, double *&aM2Psih);
		void Cal_CBEI(sComplex &aPsih, double *&aM2Psih);

		void Cal_ED(sComplex &aPsih, double *&aM2Psih);
		void Cal_HRTEM(sComplex &aPsih, double *&M2aPsih, double *&aM2Psih);

		void Cal_PED(sComplex &aPsih, double *&aM2Psih);
		void CAL_HCI(sComplex &aPsih, double *&M2aPsih, double *&aM2Psih);

		void Cal_ExitWaveFS(sComplex &aPsih, double *&aM2Psih);
		void Cal_ExitWaveRS(sComplex &aPsih, double *&aM2Psih);
};

#endif