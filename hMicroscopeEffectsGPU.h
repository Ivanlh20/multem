#ifndef hMicroscopeEffectsGPU_H
#define hMicroscopeEffectsGPU_H

#include "hConstTypes.h"
#include "cuda.h"
#include "cufft.h"

/**************************Incident wave*************************/
class cMicroscopeEffectsGPU{
	private:
		int cSynCPU;

		sQ1 Qt;
		int nQs;
		sQ2 Qs;
		sBT BT;
		sGP GP;
		sLens Lens;
		cufftHandle PlanPsi;
		double2 *Psit;

		void ReadSpatialQuadrature(sLens &Lens, int &nQ2, sQ2 &sQ2);
		void PCTCCTEM(int STEffect, double2 *&fPsi, double *&M2PsiM);
		void PCLIMWPOTEM(int STEffect, double2 *&fPsi, double *&M2PsiM);
	public:
		void freeMemory();
		cMicroscopeEffectsGPU();
		~cMicroscopeEffectsGPU();

		void SetInputData(sBT &BT_i, sGP &GP_i, sLens &Lens_i, cufftHandle &PlanPsi_i, double2 *&Psit_i);
		void ApplyMEffects(int MEffect, int STEffect, double2 *&fPsi, double *&M2Psi);
};

#endif