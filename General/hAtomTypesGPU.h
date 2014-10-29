#ifndef hAtomTypesGPU_H
#define hAtomTypesGPU_H

#include "hConstTypes.h"

class cAtomTypesGPU{
	public:
		int Z;			// Atomic number
		double m;		// Atomic mass
		int A;			// Mass number
		double rn_e;	// Experimental Nuclear radius
		double rn_c;	// Calculated Nuclear radius 
		double ra_e;	// Experimental atomic radius
		double ra_c;	// Calculated atomic radius
		double Rmin;	// Minimum interaction radius-R
		double Rmax;	// Maximum interaction radius-R
		double Rmin2;	// Squared minimum interaction radius-R
		double Rmax2;	// Squared maximum interaction radius-R

		sCoefPar cfegh;	// Electron scattering factor coefficients
		sCoefPar cfxgh;	// X-ray scattering factor coefficients
		sCoefPar cPrh;	// Potential coefficients
		sCoefPar cVrh;	// Potential coefficients
		sCoefPar cVRh;	// Projected potential coefficients

		sCoefPar cfeg;	// Electron scattering factor coefficients
		sCoefPar cfxg;	// X-ray scattering factor coefficients
		sCoefPar cPr;	// Potential coefficients
		sCoefPar cVr;	// Potential coefficients
		sCoefPar cVR;	// Projected potential coefficients

		int ns;			// Number of different sigmas
		sVoGPU *Vo;		// Optical potential coefficients

		int nR;			// Number of grid points

		double *Rh;		// R Grid
		double *R2h;	// R2 Grid
		sciVn ciVRh;	// Look up table - Projected potential coefficients

		double *R;		// R Grid
		double *R2;		// R2 Grid
		sciVn ciVR;		// Look up table - Projected potential coefficients

		void freeMemory();
		cAtomTypesGPU();
		~cAtomTypesGPU();
		void SetAtomTypes(int PotPari, sAtomTypesCPU &AtomTypesCPUi, int nRi, double dRmini);
};

#endif