#ifndef hDetectorGPU_H
#define hDetectorGPU_H

#include "hConstTypes.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>

/**************************STEM*************************/
class cDetectorGPU{
	private:
		sGP GP;
		dim3 BS;
		dim3 TS;

		double *Totp_d;
		double *Cohp_d;

		int nDet;			// number of Detectors
		sDetCir DetCirh;	// Detectors

		double *Tot_h;		// Total Intensity
		double *Coh_h;		// Coherent Intensity

		double *Tot_d;		// Total Intensity
		double *Coh_d;		// Coherent Intensity
	public:
		void freeMemory();
		cDetectorGPU();
		~cDetectorGPU();

		void SetInputData(sGP &GP_i, sBT &BT_i, int nDeti, sDetCir &DetCirhi);
		void getDetectorIntensity(double *&aM2Psi, double *&M2aPsi, int ixys, sDetInt *DetInth);
		void getDetectorIntensity(double *&aM2Psi, int ixys, sDetInt *DetInth, bool add=false);
};

#endif