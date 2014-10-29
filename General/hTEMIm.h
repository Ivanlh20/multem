#ifndef hTEMIm_H
#define hTEMIm_H

#include "hConstTypes.h"
#include "hMicroscopeEffectsGPU.h"
#include <cufft.h>
#include <vector_types.h>

/***************************TEM - imaging*****************************/
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

		/************************************************************/
		double lambda;		// wavelenght

		sGP GP;				// Grid variables
		sBT BT;				// Blocks and Threads
		sLens Lens;			// Aberration parameters

		/***********************************************************/
		cufftHandle PlanPsi;	// Fourier transform's plan

		/***********************************************************/
		double2 *fPsi;		// Wave function
		double2 *Psia;		// Wave function - aberrations

		/***********************************************************/
		double *M2Psis;		// Squared Wave function - spatial incoherents
		double *M2Psit;		// Squared Wave function - temporal incoherents

		/***********************************************************/
		cMicroscopeEffectsGPU MicroscopeEffectsGPU;

		/***********************************************************/
		cTEMIm();
		void freeMemory();

		void SetInputData(sInTEMIm &InTEMIm);
		void TEMImage(double *Psir_hi, double *Psii_hi, double *M2Psi_ho);
};

#endif