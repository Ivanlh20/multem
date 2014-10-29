#ifndef hIncidentWaveGPU_H
#define hIncidentWaveGPU_H

#include "hConstTypes.h"
#include <cufft.h>
#include <vector_types.h>

/**************************Incident wave*************************/
class cIncidentWaveGPU{
	private:
		sBT BT;
		sGP GP;
		sLens Lens;
		double *Sd;
		cufftHandle PlanPsi;
	public:
		void SetInputData(sBT &BT_i, sGP &GP_i, sLens &Lens_i, cufftHandle &PlanPsi_i);
		void Psi0(double2 *&Psig);
		void Psi0(double x, double y, double2 *&Psig);

		void freeMemory();
		cIncidentWaveGPU();
		~cIncidentWaveGPU();
		
};

#endif