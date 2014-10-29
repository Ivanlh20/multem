#ifndef hPotentialGPU_H
#define hPotentialGPU_H

#include "hConstTypes.h"
#include "hSpecimenCPU.h"
#include "hAtomTypesGPU.h"

/*************************Potential**************************/
class cPotentialGPU{
	private:
		sQ1 Qz;
		scVp *cVph;

		sMGP MGP;						// Multislice general parameters
		sBT BT;							// Blocks and threads
		sGP GP;							// xy-Grid properties

		int nAtomTypesGPU;				// Number of atom types
		cAtomTypesGPU *AtomTypesGPU;	// Atom types

		void SetPotPar(int PotParh);
		int CheckGridLimits(int i, int n);
		void getbn(sGP &GP, double x, double y, double Rmax, sbn &bnx, sbn &bny);
		void setcVp(int ApproxModel, cSpecimenCPU &Specimen, cAtomTypesGPU *&AtomTypesGPU, int iSlice, int iatom, int nsatom, dim3 &BPot, dim3 &TPot, dim3 &BCoef, dim3 &TCoef, dim3 &BEval, dim3 &TEval);
	public:
		cSpecimenCPU Specimen;	// Specimen

		double *V0;				// Zero moment of the potential
		double *V1;				// first moment of the potential
		double *V2;				// Second moment of the potential

		void freeMemory();	
		cPotentialGPU();
		~cPotentialGPU();

		void SetInputData(sMGP &MGP_io, sGP &GP_i, sBT &BT_i, int nAtomsM_i, double *AtomsM_i);
		void ProjectedPotential_Slice(int iSlice);
		void ProjectedPotential_Specimen(int iSlice);
};

#endif