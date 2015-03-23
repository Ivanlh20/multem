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

#ifndef hMT_MulSli_GPU_H
#define hMT_MulSli_GPU_H

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

class cMT_STEM_GPU;

class cMT_MulSli_GPU{
	private:
		int IdCall;
	public:
		int cSynCPU;
		sLens Lens;													// Aberration parameters
		sGP GP;														// Grid variables

		double gxu;													// incident beam x-tilt in reciprocal units
		double gyu;													// incident beam y-tilt in reciprocal units

		sACD ExpRg_x;												// k_PhaseMul x
		sACD ExpRg_y;												// k_PhaseMul y
		sACD Prop_x;												// Propagator x
		sACD Prop_y;												// Propagator y

		double2 *MC_h;												// Host: Complex matrix

		double2 *Psi;												// Wave function
		double2 *aPsi;												// Wave function - temporal

		double *M2aPsi;												// Squared Wave function
		double *aM2Psi;												// Squared Wave function - temporal

		cufftHandle PlanPsi;										// Fourier transform's plan
		
		cMT_IncidentWave_GPU *MT_IncidentWave_GPU;					// Incident wave;
		cMT_Transmission_GPU *MT_Transmission_GPU;					// Transmission function
		cMT_MicroscopeEffects_GPU *MT_MicroscopeEffects_GPU;		// Microscope effects

		void PhaseMul(double gxu, double gyu, double2 *&Psi);
		void Propagate(eSpace Space, double gxu, double gyu, double z, double2 *&Psi);
		void Cal_Wavefunction(eSpace Space, double2 *&Psi);

		void Image_Plane_Wave_Illumination(int nConfFP, eSpace Space, int MEffect, int STEffect, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi);
		void Image_Convergence_Wave_Illumination(int nConfFP, eSpace Space, double xi, double yi, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi);

		void Gather(double2 *&MC_d_i, sComplex &MC_h_o);
		void Gather(double2 *&MC_d_i, double *&MD_d_i, sComplex &MC_h_o, double *&MD_h_o);
		void Gather(double2 *&MC_d_i, double *&MD1_d_i, double *&MD2_d_i, sComplex &MC_h_o, double *&MD1_h_o, double *&MD2_h_o);

		cMT_MGP_CPU MT_MGP_CPU;
		cMT_STEM_GPU *STEM;
		sCBED CBED;
		sCBEI CBEI;
		sHRTEM HRTEM;
		sPED PED;
		sHCI HCI;
		sEWRS EWRS;
		sEWFS EWFS;	
	
		void freeMemory();
		cMT_MulSli_GPU();
		~cMT_MulSli_GPU();

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