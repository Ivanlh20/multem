#ifndef hAtomTypesGPU_H
#define hAtomTypesGPU_H

#include "hConstTypes.h"

class cAtomTypesGPU{
	private:
		void freeMemory();
	public:
		int Z;			// Atomic number
		double m;		// Atomic mass
		int A;			// Mass number
		double rn_e;		// Experimental Nuclear radius
		double rn_c;		// Calculated Nuclear radius 
		double ra_e;		// Experimental atomic radius
		double ra_c;		// Calculated atomic radius
		double Rmin;	// Minimum interaction radius-R
		double Rmax;	// Maximum interaction radius-R
		double Rmin2;	// Squared minimum interaction radius-R
		double Rmax2;	// Squared maximum interaction radius-R
		int PotPar;			// Parameterization type
		double occ;		// occupancy

		sCoefCPU fegh;	// Electron scattering factor coefficients
		sCoefCPU fxgh;	// X-ray scattering factor coefficients
		sCoefCPU Prh;	// Potential coefficients
		sCoefCPU Vrh;	// Potential coefficients
		sCoefCPU VRh;	// Projected potential coefficients

		sCoefGPU feg;	// Electron scattering factor coefficients
		sCoefGPU fxg;	// X-ray scattering factor coefficients
		sCoefGPU Pr;	// Potential coefficients
		sCoefGPU Vr;	// Potential coefficients
		sCoefGPU VR;	// Projected potential coefficients

		int ns;			// Number of different sigmas
		sVoGPU *Vo;		// Optical potential coefficients

		int nR;			// Number of grid points
		double *R;		// R Grid
		double *R2;		// R2 Grid
		
		cAtomTypesGPU();
		~cAtomTypesGPU();
		void SetAtomTypes(sAtomTypesCPU AtomTypesCPUi, int nRi, double dRmini);
};

#endif