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

#ifndef hTEMIm_H
#define hTEMIm_H

#include "hConstTypes.h"
#include "hMT_MicroscopeEffects_GPU.h"
#include <cufft.h>
#include <vector_types.h>

/**************************TEM - imaging*****************************/
class cTEMIm{
	private:
		void GenerateParameters();
	public:
		int gpu;			// gpu card
		int MEffect;		// 0: Exit wave Partial coherente mode, 1: Transmission cross coefficient
		double E0;			// Acceleration volatage in KeV
		double *Psirh;		// Real part of the Wave function
		double *Psiih;		// Imaginary part of the Wave function
		double lx;			// distance in x direction(Angstroms)
		double ly;			// distance in y direction(Angstroms)
		int nx;				// Number of pixels in x direction
		int ny;				// Number of pixels in y direction

		/***********************************************************/
		double lambda;		// wavelenght

		sGP GP;				// Grid variables
		sBT BT;				// Blocks and Threads
		sLens Lens;			// Aberration parameters

		/**********************************************************/
		cufftHandle PlanPsi;	// Fourier transform's plan

		/**********************************************************/
		double2 *fPsi;		// Wave function
		double2 *Psia;		// Wave function - aberrations

		/**********************************************************/
		double *M2Psis;		// Squared Wave function - spatial incoherents
		double *M2Psit;		// Squared Wave function - temporal incoherents

		/**********************************************************/
		cMT_MicroscopeEffects_GPU MT_MicroscopeEffects_GPU;

		/**********************************************************/
		cTEMIm();
		void freeMemory();

		void SetInputData(sInTEMIm &InTEMIm);
		void TEMImage(double *Psir_hi, double *Psii_hi, double *M2Psi_ho);
};

#endif