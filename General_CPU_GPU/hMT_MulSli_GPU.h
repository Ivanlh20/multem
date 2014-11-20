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

#ifndef hMT_MulSli_GPU_H
#define hMT_MulSli_GPU_H

#include "hConstTypes.h"
#include "hMT_General_CPU.h"
#include "hMT_Specimen_CPU.h"
#include "hMT_AtomTypes_GPU.h"
#include "hMT_Potential_GPU.h"
#include "hMT_IncidentWave_GPU.h"
#include "hMT_MicroscopeEffects_GPU.h"
#include "hMT_STEM_GPU.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

/******************************CBED*****************************/
struct sCBED{
	double x0;
	double y0;
	eSpace Space;
};

/*****************************HRTEM*****************************/
struct sHRTEM{
	int xx;
};

/*************************Precession***************************/
struct sPED {
	int nrot;				// Total number of orientations
	double theta;			// Precession angle in rad
};

/********************Hollow cone ilumination*******************/
struct sHCI{
	int nrot;				// Total number of orientations
	double theta;			// Precession angle in rad
};

/*********************Exit Wave Real Space*********************/
struct sEWRS{
	int xx;
};

/******************Exit Wave Fourier Space*********************/
struct sEWFS{
	int xx;
};

class cMT_MulSli_GPU{
	private:
		int cSynCPU;

		sLens Lens;										// Aberration parameters
		sGP GP;											// Grid variables

		double gxu;										// incident beam x-tilt in reciprocal units
		double gyu;										// incident beam y-tilt in reciprocal units

		double fPot;

		sACD ExpRg_x;									// k_PhaseMul x
		sACD ExpRg_y;									// k_PhaseMul y
		sACD Prop_x;									// Propagator x
		sACD Prop_y;									// Propagator y

		double2 *Trans;									// Transmission function

		double2 *Psi;									// Wave function
		double2 *aPsi;									// Wave function - temporal

		double *M2aPsi;									// Squared Wave function
		double *aM2Psi;									// Squared Wave function - temporal

		cufftHandle PlanPsi;							// Fourier transform's plan

		cMT_Potential_GPU MT_Potential_GPU;						// Potential
		cMT_IncidentWave_GPU MT_IncidentWave_GPU;				// Incident wave;
		cMT_MicroscopeEffects_GPU MT_MicroscopeEffects_GPU;		// Microscope effects

		void PhaseMul(double gxu, double gyu, double2 *&Psi);
		void Transmission(int iSlice, double fPot, double2 *&Trans);
		void Transmit(double2 *&Trans, double2 *&Psi);
		void Propagate(eSpace Space, double gxu, double gyu, double z, double2 *&Psi);
		/***************************************************************************/
		void Cal_Wavefunction_POA_WPOA(eSpace Space, double2 *&Psi);
		void Cal_Wavefunction_PA(eSpace Space, double2 *&Psi);
		void Cal_Wavefunction_MSA(eSpace Space, double2 *&Psi);
		/***************************************************************************/
		void Cal_FAST_STEM_Wavefunction_POA_WPOA(int nConfFP, sDetInt *DetInt);
		void Cal_FAST_STEM_Wavefunction_PA(int nConfFP, sDetInt *DetInt);
		void Cal_FAST_STEM_Wavefunction_MSA(int nConfFP, sDetInt *DetInt);
		/***************************************************************************/
		void Cal_Wavefunction(eSpace Space, double2 *&Psi);

		void get_Imagefuncion(int nConfFP, eSpace Space, int MEffect, int STEffect, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi);
		void get_Imagefuncion(int nConfFP, eSpace Space, double xi, double yi, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi);
		void GPU2CPU(double2 *&Psid, sComplex &Psih);
		void GPU2CPU(double2 *&Psid, double *&M2Psid, sComplex &Psih, double *&M2Psih);
		void GPU2CPU(double2 *&Psid, double *&M2Psi1d, double *&M2Psi2d, sComplex &Psih, double *&M2Psi1h, double *&M2Psi2h);
	public:
		sMGP MGP;

		cSTEM_GPU STEM;
		sCBED CBED;
		sHRTEM HRTEM;
		sPED PED;
		sHCI HCI;
		sEWRS EWRS;
		sEWFS EWFS;	
	
		void freeMemory();
		cMT_MulSli_GPU();
		~cMT_MulSli_GPU();

		void SetInputData(sInMSTEM &InMSTEM);
		void Cal_ExitWaveRS(sComplex &aPsih, double *&aM2Psih);
		void Cal_ExitWaveFS(sComplex &aPsih, double *&aM2Psih);
		void Cal_ED(sComplex &aPsih, double *&aM2Psih);
		void Cal_HRTEM(sComplex &aPsih, double *&M2aPsih, double *&aM2Psih);
		void Cal_CBED(sComplex &aPsih, double *&aM2Psih);
		void Cal_STEM();
};

#endif