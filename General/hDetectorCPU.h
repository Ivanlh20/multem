#ifndef hDetectorCPU_H
#define hDetectorCPU_H

#include "hConstTypes.h"


/**************************STEM*************************/
class cDetectorCPU{
	private:
		sGP GP;

		int nDet;			// number of Detectors
		sDetCir DetCir;	// Detectors

		double *Tot;		// Total Intensity
		double *Coh;		// Coherent Intensity
	public:
		void freeMemory();
		cDetectorCPU();
		~cDetectorCPU();

		void SetInputData(sGP &GP_i, int nDeti, sDetCir &DetCiri);
		void getDetectorIntensity(double *&Toth, double *&Cohh, int ixys, sDetInt *DetInth);
};

#endif