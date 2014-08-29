#ifndef hMulSliGPU_H
#define hMulSliGPU_H

#include "..\General\hConstTypes.h"
#include "..\General\hSpecimenCPU.h"
#include "..\General\hAtomTypesGPU.h"
#include "..\General\hPotentialGPU.h"
#include "..\General\hIncidentWaveGPU.h"
#include "..\General\hMicroscopeEffectsGPU.h"
#include "..\General\hSTEMGPU.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include <device_functions.h>
#include "cufft.h"

/**************************HRTEM**************************/
struct sHRTEM{
	int MEffect;				// 1: Exit wave Partial coherente mode, 2: Transmission cross coefficient
	int STEffect;				// 1: Spatial and temporal, 2: Temporal, 3: Spatial
	int ZeroDefTyp;				// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
	double ZeroDefPlane;		// Zero defocus plane
};

/******************************Precession****************************/
struct sPED {
	int nrot;				// Total number of orientations
	double theta;			// Precession angle in rad
};

/**********************Hollow cone ilumination***********************/
struct sHCI{
	int nrot;				// Total number of orientations
	double theta;			// Precession angle in rad
	int MEffect;			// 1: Exit wave Partial coherente mode, 2: Transmission cross coefficient
	int STEffect;			// 1: Spatial and temporal, 2: Temporal, 3: Spatial
	int ZeroDefTyp;			// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
	double ZeroDefPlane;	// Zero defocus plane
};

/*************************Exit Wave Real Space************************/
struct sEWRS{
	int ZeroDefTyp;			// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
	double ZeroDefPlane;	// Zero defocus plane
};

/***********************Exit Wave Fourier Space***********************/
struct sEWFS{
	int ZeroDefTyp;			// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
	double ZeroDefPlane;	// Zero defocus plane
};

class cMulSliGPU{
	private:
		int cSynCPU;

		sLens Lens;					// **Aberration parameters
		sGP GP;						// Grid variables
		sBT BT;						// Blocks and Threads

		double gxu;					// incident beam x-tilt in reciprocal units
		double gyu;					// incident beam y-tilt in reciprocal units

		double fPot;
		double fProp;

		sERg ERg;					// Exp2ipiRgu
		sProp Prop;					// Propagator
		double2 *Trans;				// Transmission function

		double2 *Psi;				// Wave function
		double2 *aPsi;				// Wave function - temporal

		double *M2aPsi;				// Squared Wave function
		double *aM2Psi;				// Squared Wave function - temporal

		cufftHandle PlanPsi;		// Fourier transform's plan

		cPotentialGPU PotentialGPU;						// Potential
		cIncidentWaveGPU IncidentWaveGPU;				// Incident wave;
		cMicroscopeEffectsGPU MicroscopeEffectsGPU;		// Microscope effects

		//cSTEMGPU STEMGPU;
		sHRTEM HRTEM;
		sPED PED;
		sHCI HCI;
		sEWRS EWRS;
		sEWFS EWFS;		

		void Transmission(dim3 grid, dim3 threads, int MO, int SlicePos, int nxy, double f, double *Vpg, double *zVpg, double *zVpog, double2 *Transg);
		void Bandlimited2D(dim3 grid, dim3 threads, cufftHandle &PlanPsi, sGP &GP, double2 *Psi);
		void MSExitWaveCal(int iConf, int iSpace, double2 *&Psi);
		void Wavefuncion(int nConfFP, int Space, int MEffect, int STEffect, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi);
		void Gpu2CPU(double2 *&Psid, sComplex &Psih);
		void Gpu2CPU(double2 *&Psid, double *&M2Psid, sComplex &Psih, double *&M2Psih);
		void Gpu2CPU(double2 *&Psid, double *&M2Psi1d, double *M2Psi2d, sComplex &Psih, double *&M2Psi1h, double *&M2Psi2h);
	public:
		int gpu;					// **gpu device
		int SimType;				// **1: STEM, 2: CBED, 3: HRTEM, 4: ED, 5: PED, 6: HCI, ... 10: EW real, 11: EW Fourier
		int MulOrder;				// **1: First order MS, 2: Second Order MS
		int nConfFP;				// **Number of frozen phonon configurations
		int DimFP;					// **Dimensions phonon configurations
		int DistFP;					// **Frozen phonon distribution type 1:Normal, 2:xx
		int SeedFP;					// **Random seed(frozen phonon)
		int PotPar;					// **Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) and 6: Lobato(0-12)

		double E0;					// **Acceleration volatage in KeV
		double theta;				// **Tilt (in spherical coordinates) (rad)
		double phi;					// **Tilt (in spherical coordinates) (rad)

		void freeMemory();
		cMulSliGPU();
		~cMulSliGPU();

		void SetInputData(sInMSTEM &InMSTEM);
		void Cal_ExitWavefuncionRS(sComplex &aPsih, double *&aM2Psih);
		void Cal_ExitWavefuncionFS(sComplex &aPsih, double *&aM2Psih);
		void Cal_ED(sComplex &aPsih, double *&aM2Psih);
		void Cal_HRTEM(sComplex &aPsih, double *&M2aPsih, double *&aM2Psih);
};

#endif