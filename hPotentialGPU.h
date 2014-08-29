#ifndef hPotentialGPU_H
#define hPotentialGPU_H

#include "hConstTypes.h"
#include "hSpecimenCPU.h"
#include "hAtomTypesGPU.h"

/*************************Potential**************************/
class cPotentialGPU{
	private:
		sQ1 Qz;
		scVp cVph[stncVp];

		sBT BT;					// Blocks and threads
		sGP GP;					// xy-Grid properties

		int nAtomTypesGPU;				// Number of atom types
		cAtomTypesGPU *AtomTypesGPU;	// Atom types

		int MulOrder;			// MultiSlice order	

		void cudaInit_sciVn(sciVn &ciVn);
		void cudaFree_sciVn(sciVn &ciVn);
		void cudaMalloc_sciVn(int n, sciVn &ciVn);
		void getbn(sGP &GP, double x, double y, double Rmax, sbn &bnx, sbn &bny);
		void setcVp(cSpecimen &Specimen, cAtomTypesGPU *&AtomTypesGPU, int iSli, int iatom, int nsatom, dim3 &BPot, dim3 &TPot, dim3 &BCoef, dim3 &TCoef, dim3 &BEval, dim3 &TEval);
	public:
		cSpecimen Specimen;		// Specimen

		double *V0;				// Zero moment of the potential
		double *V1;				// first moment of the potential
		double *V2;				// Second moment of the potential

		void freeMemory();	
		cPotentialGPU();
		~cPotentialGPU();

		void SetInputData(sGP &GP_i, sBT &BT_i, int nAtomsMi, double *AtomsMi, double dzi, int PotPari, double Vrli, int DimFPi, int DistFPi, int SeedFPi, int ZeroDefTypi, double ZeroDefPlanei, int MulOrder_i);
		void ProjectedPotential(int iSli);
};

#endif