#ifndef hMulSliGPU_H
#define hMulSliGPU_H

#include "hConstTypes.h"
#include "hgeneralCPU.h"
#include "hSpecimenCPU.h"
#include "hAtomTypesGPU.h"
#include "hPotentialGPU.h"
#include "hIncidentWaveGPU.h"
#include "hMicroscopeEffectsGPU.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

/******************************STEM********************************/
struct sSTEM{
	int line;					// 0: Area, 1: Line
	int ns;						// Sampling points
	double x1u;					// Initial scanning position in x
	double y1u;					// Initial scanning in y
	double x2u;					// final scanning position in x
	double y2u;					// final scanning position in y
	int nThk;
	int nDet;					// Number of circular detectors
	sDetCir DetCir;				// Circular detectors
	sImSTEM *ImSTEM;
	bool FastCal;				// 0: normal mode(low memory consumption), 1: fast calculation(high memory consumption)
	int nxs;
	int nys;
	double *xs;
	double *ys;

	int nst;
	double *xst;
	double *yst;

	void Init();
	void free();

	void SetInputData(sInMSTEM &InMSTEM);
	void InitImSTEM();
};

/******************************CBED*****************************/
struct sCBED{
	double x0;
	double y0;
};

/*****************************HRTEM*****************************/
struct sHRTEM{
	int xx;
};

/*************************Precession***************************/
struct sPED {
	int nrot;				// Total number of orientations
	double theta;			// Precession angle in rad
};

/********************Hollow cone ilumination*******************/
struct sHCI{
	int nrot;				// Total number of orientations
	double theta;			// Precession angle in rad
};

/*********************Exit Wave Real Space*********************/
struct sEWRS{
	int xx;
};

/******************Exit Wave Fourier Space*********************/
struct sEWFS{
	int xx;
};

class cMulSliGPU{
	private:
		int cSynCPU;

		sLens Lens;					// Aberration parameters
		sGP GP;						// Grid variables
		sBT BT;						// Blocks and Threads

		double gxu;					// incident beam x-tilt in reciprocal units
		double gyu;					// incident beam y-tilt in reciprocal units

		double fPot;

		sERg ERg;					// k_PhaseMul
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

		sHRTEM HRTEM;
		sPED PED;
		sHCI HCI;
		sEWRS EWRS;
		sEWFS EWFS;	
		sCBED CBED;
		void PhaseMul(double gxu, double gyu, double2 *&Psi);
		void TransmissionWPO(double fPot, double2 *&Trans);
		void Transmission(eSlicePos SlicePos, double fPot, double2 *&Trans);
		void Transmit(double2 *&Trans, double2 *&Psi);
		void Propagate(eSpace Space, double gxu, double gyu, double z, double2 *&Psi);
		/***************************************************************************/
		void Cal_Wavefunction_WPOA(eSpace Space, double2 *&Psi);
		void Cal_Wavefunction_POA(eSpace Space, double2 *&Psi);
		void Cal_Wavefunction_PA(eSpace Space, double2 *&Psi);
		void Cal_Wavefunction_MSA(eSpace Space, double2 *&Psi);
		/***************************************************************************/
		void Cal_FAST_STEM_Wavefunction_WPOA(int nConfFP, sDetInt *DetInt);
		void Cal_FAST_STEM_Wavefunction_POA(int nConfFP, sDetInt *DetInt);
		void Cal_FAST_STEM_Wavefunction_PA(int nConfFP, sDetInt *DetInt);
		void Cal_FAST_STEM_Wavefunction_MSA(int nConfFP, sDetInt *DetInt);
		/***************************************************************************/
		void Cal_Wavefunction(eSpace Space, double2 *&Psi);

		void get_Imagefuncion(int nConfFP, eSpace Space, int MEffect, int STEffect, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi);
		void get_Imagefuncion(int nConfFP, eSpace Space, double xi, double yi, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi);
		void Gpu2CPU(double2 *&Psid, sComplex &Psih);
		void Gpu2CPU(double2 *&Psid, double *&M2Psid, sComplex &Psih, double *&M2Psih);
		void Gpu2CPU(double2 *&Psid, double *&M2Psi1d, double *&M2Psi2d, sComplex &Psih, double *&M2Psi1h, double *&M2Psi2h);
	public:
		sMGP MGP;
		sSTEM STEM;

		void freeMemory();
		cMulSliGPU();
		~cMulSliGPU();

		void SetInputData(sInMSTEM &InMSTEM);
		void Cal_ExitWaveRS(sComplex &aPsih, double *&aM2Psih);
		void Cal_ExitWaveFS(sComplex &aPsih, double *&aM2Psih);
		void Cal_ED(sComplex &aPsih, double *&aM2Psih);
		void Cal_HRTEM(sComplex &aPsih, double *&M2aPsih, double *&aM2Psih);
		void Cal_CBED(sComplex &aPsih, double *&aM2Psih);
		void Cal_STEM();
};

#endif