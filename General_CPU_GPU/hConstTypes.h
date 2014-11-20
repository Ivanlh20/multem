/**
 *  This file is part of MULTEM.
 *  Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
 *
 *  MULTEM is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MULTEM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MULTEM.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef hConstTypes_H
#define hConstTypes_H

#include <cstring>
#include "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v6.5\include\vector_types.h"

#define cHa 27.2113850656389				// Hartree to electron-Volt
#define ca0 0.52917721077817892				// Bohr radius
#define cPotf 47.877645145863056			//
#define c2Pi2a0 10.445539456905012			// 2*pi^2*a0

#define cH 6.62606876e-34					// Planck's constant - J s
#define cbH 1.054571596e-34					// h/(2*pi) - J s
#define cC 2.99792458e+8					// Velocity of light - m s^-1 
#define cQe 1.602176462e-19					// Elementary charge
#define cme 9.10938291e-31					// Electron rest mass [kg]
#define cmp 1.672621637e-27					// Proton rest mass [kg]
#define cKB 1.3806504e-23					// Boltzmann's constant - J K^-1
#define cNa 6.0221415e+23					// Avogadro's Number - mol^-1

#define cE 2.7182818284590452354			// e (base of natural log)

#define cPi 3.141592653589793238463			// pi
#define ciPi 0.3183098861837906715378		// 1.0/pi
#define ci2Pi 1.570796326794896619231		// pi/2
#define ci3Pi 1.047197551196597746154		// pi/3
#define ci4Pi 0.7853981633974483096157		// pi/4
#define c2Pi 6.283185307179586476925		// 2*pi
#define c3Pi 9.424777960769379715388		// 3*pi
#define c4Pi 12.56637061435917295385		// 4*pi
#define cPi2 9.869604401089358618834		// pi^2
#define cPi3 31.00627668029982017548		// pi^3
#define cPi4 97.4090910340024372364			// pi^4
#define cPii2 1.772453850905516027298		// pi^(1/2)
#define cPii3 1.46459188756152326302		// pi^(1/3)
#define cPii4 1.331335363800389712798		// pi^(1/4)

#define c2i2 1.414213562373095048802		// 2^(1/2)
#define c3i2 1.732050807568877293527		// 3^(1/2)
#define c5i2 2.236067977499789696409		// 5^(1/2)
#define c7i2 2.645751311064590590502		// 7^(1/2)	

#define mrad2rad 1.0e-03					// mrad-->rad
#define deg2rad 0.01745329251994329576924	// degrees-->rad
#define mm2Ags 1.0e+07						// mm-->Angstrom

#define eed 1e-14							// error double precision

#define ccSynCPU 5;

#define NE 103
#define stncVp 16
#define stnQz 128
#define stnR 128
#define stngbp 128 // grid blurred potential
#define stnR0c 127
#define stnz0c 256
#define stnr0c 256

#define thrnxny 16
#define thrnxy 256

#define stqh 0.075
#define stqbeta 0.25

#define stVrl 0.015

#define ABS(a)		(((a) < 0)? -(a) : (a))
#define MIN(a,b)	(((a) < (b))?(a) : (b))
#define MAX(a,b)	(((a) < (b))?(b) : (a))
#define SQR(x)		((x)*(x))
#define CUB(x)		((x)*(x)*(x))
#define IsFS(i,nh)	((i<nh)?i:i-2*nh)
#define IsRS(i,nh)	((i<nh)?i+nh:i-nh)

#define cudaFreen(x)	if(x != 0){\
							cudaFree(x);\
							x = 0;\
						}

#define cufftDestroyn(x)	if(x != 0){\
								cufftDestroy(x);\
								x = 0;\
							}

const int cSizeofI = sizeof(int);
const int cSizeofRD = sizeof(double);
const int cSizeofRF = sizeof(float);
const int cSizeofCD = sizeof(double2);
const int cMb = 1024*1024;

enum eMULTEMType{
	eMTT_FLOAT=1, eMTT_DOUBLE=2
};

/**************************Real or Fourier space****************************/
enum eSpace{
	eSReal = 1, eSReciprocal = 2
};

enum eSlicePos{
	eSPFirst = 1, eSPMedium = 2, eSPLast = 3
};

/**************************Complex****************************/
typedef struct sComplex{
	double *real;		// real part
	double *imag;		// imaginary part

} sComplex;

/*************Scattering factors parameterization************/
typedef struct sCoefCPU{
	// Lineal coefficients fep
	double cl[6];
	// Non-Lineal coefficients fep
	double cnl[6];
} sCoefCPU;

typedef struct sCoefPar{
	// Lineal coefficients fep
	double *cl;
	// Non-Lineal coefficients fep
	double *cnl;
} sCoefPar;

typedef struct sfep{	
	sCoefCPU Par[6];
} sfep;

/******************Atom inside the specimen******************/
typedef struct sAtoms {
	int Z;			// Atomic number
	double x;		// x coordinate
	double y;		// y coordinate
	double z;		// z coordinate
	double sigma;	// Isotropic 3D-RMS displacament 
	double occ;		// Occupancy
}sAtoms;

const int cSizeofAtoms = sizeof(sAtoms);

/************************Atom group**************************/
typedef struct sAtomsGroup {
	int nAtoms;
	sAtoms *Atoms;
}sAtomsGroup;

/************************Slice thickness*********************/
typedef struct sSlice{
	double z0;		// Initial z-position
	double ze;		// Final z-position
	int z0_id;		// Index to initial z-position
	int ze_id;		// Index to final z-position
	double z0i;		// Initial interaction z-position
	double zei;		// Final interaction z-position
	int z0i_id;		// Index to initial interaction z-position
	int zei_id;		// Index to final interaction z-position
} sSlice;

const int cSizeofSlice = sizeof(sSlice);

/************************Atomic data*************************/
typedef struct sAtDa{
	// atomic number
	int Z;
	// atomic mass
	double m;
	// mass number
	int A;
	// nuclear radius
	double rn;
	// atomic radius
	double ra;
} sAtDa;

/************Optical potential coefficients CPU**************/
typedef struct sVoCPU{
	double sigma;		// Isotropic standard deviation
	sCoefCPU cVr;		// Real Potential coefficients
	sCoefCPU cVi;		// Imaginary Potential coefficients

	double gr[stngbp];
	double gVr[stngbp];
	double gVi[stngbp];
} sVoCPU;

/************Optical potential coefficients GPU**************/
typedef struct sVoGPU{
	double sigma;		// Isotropic standard deviation
	sCoefPar cVr;		// Real Potential coefficients
	sCoefPar cVi;		// Imaginary Potential coefficients

	//double gr[stngbp];
	//double gVr[stngbp];
	//double gVi[stngbp];
} sVoGPU;

/*************************Atomic type************************/
typedef struct sAtomTypesCPU{
	int Z;			// Atomic number
	double m;		// Atomic mass
	int A;			// Mass number
	double rn_e;	// Experimental Nuclear radius
	double rn_c;	// Calculated Nuclear radius 
	double ra_e;	// Experimental atomic radius
	double ra_c;	// Calculated atomic radius
	double Rmin;	// Minimum interaction radius-R
	double Rmax;	// Maximum interaction radius-R

	sCoefCPU cfeg;	// Electron scattering factor coefficients
	sCoefCPU cfxg;	// X-ray scattering factor coefficients
	sCoefCPU cPr;	// Potential coefficients
	sCoefCPU cVr;	// Potential coefficients
	sCoefCPU cVR;	// Projected potential coefficients

	int ns;			// Number of different sigmas
	sVoCPU *Vo;		// Optical potential coefficients + grid
} sAtomTypesCPU;

/*************************Quadrature*************************/
typedef struct sHt{
	int n;
	int *c;
	double *v;
}sHt;

/*************************Quadrature*************************/
typedef struct sQ1{
	double * __restrict x;
	double * __restrict w;
}sQ1;

typedef struct sQ2{
	double * __restrict x;
	double * __restrict y;
	double * __restrict w;
} sQ2;

/*********************Quadrature Parameters******************/
typedef struct sQp{
	double *r;
	double *rlim;
	double *a;
	double *b;
	double *tmin;
	double *h;
} sQp;

/****************************grid****************************/
typedef struct sGP{
	int nx;				// **Number of pixels in x direction
	int ny;				// **Number of pixels in y direction

	int nxh;			// Half number of pixels in x direction
	int nyh;			// Half number of pixels in y direction

	int nxy;			// nx*ny
	double inxy;		// 1.0/nxy

	double lx;			// **Box size in x direction(Angstroms)
	double ly;			// **Box size in y direction(Angstroms)
	double dz;			// **slice thickness
	bool PBC_xy;		// **Peridic boundary contions
	bool BWL;			// **Band width limit
	double dRx;			// x-sampling in real Space
	double dRy;			// y-sampling in real Space
	double dRmin;		// Minimum sampling in real Space

	double dgx;			// x-sampling in reciprocal Space
	double dgy;			// y-sampling in reciprocal Space
	double dgmin;		// Minimum sampling in reciprocal Space

	double gmax;		// Maximun frequency
	double gmax2;		// Squared of the maximun frequency

	double gmaxl;		// Maximun limited frequency 
	double gmaxl2;		// Squared of the maximun limited frequency 
} sGP;

/*******************Momentum of the potential*********************/
typedef struct smVn{
	double * __restrict Vn;		// int_{z_i}^{z_i+dz}(z-zi)^n V(R, z)
	double * __restrict dVn;	// dVn/dR^2

} smVn; 

/****************Cubic interpolation coefficients*****************/
typedef struct sciVn{
	double * __restrict c0;		// zero coefficient
	double * __restrict c1;		// first coefficient
	double * __restrict c2;		// second coefficient
	double * __restrict c3;		// third coefficient
} sciVn;

/***********************structure borders*************************/
struct sbn{
	int i;
	int n;
};

/******************************scVp*******************************/
struct scVp{
	double x;
	double y;
	double z0;
	double ze;
	bool split;
	double occ;
	double Rmin2;
	double Rmax2;
	double *R2;
	sCoefPar cVr;
	sciVn ciV0;
	sciVn ciV1;
	sbn bnx;
	sbn bny;
};

const int cSizeofcVp = sizeof(scVp);

/*****************************Lens*******************************/
typedef struct sLens{
	double gamma;		// Relativistic factor
	double lambda;		// wavelength(Angstrom)

	int m;				// vortex momentum
	double f;			// defocus
	double Cs3;			// spherical aberration(Angstrom)
	double Cs5;			// spherical aberration(Angstrom)
	double mfa2;		// magnitude 2-fold astigmatism
	double afa2;		// angle 2-fold astigmatism(rad)
	double mfa3;		// magnitude 3-fold astigmatism
	double afa3;		// angle 3-fold astigmatism(rad)
	double aobjl;		// lower objective aperture(rad);
	double aobju;		// upper objective aperture(rad);

	// temporal aberrations
	double sf;			// defocus spread
	int nsf;			// Number of defocus sampling points

	// spatial aberrations
	double beta;		// semi-convergence angle
	int nbeta;			// Number of semi-convergence angle sampling points

	/**************************************************/
	double lambda2;		// wavelength(Angstrom)^2

	double cf;			// pi*f*lambda
	double cCs3;		// -0.5*pi*Cs3*lambda^3
	double cCs5;		// -pi*Cs5*lambda^5/3
	double cmfa2;		// -pi*lambda
	double cmfa3;		// -2*pi*lambda^3/3
	double gmin2;		// aobjl/lambda
	double gmax2;		// aobju/lambda

	double sggs;		// Standard deviation
	int ngxs;			// Number of source sampling points x
	int ngys;			// Number of source sampling points y
	double dgxs;		// source sampling size;
	double dgys;		// source sampling size;
	double gmax2s;		// q maximum square;

} sLens;

/***************************k_PhaseMul***************************/
typedef struct sACD{
	double * __restrict x;
	double * __restrict y;
} sACD;

/*************************STEM Detector**************************/
typedef struct sInDetCir{
	double InnerAng;	// Inner aperture (rad)
	double OuterAng;	// Outer aperture (rad)
} sInDetCir;

typedef struct sDetCir{
	double * __restrict g2min;		// Inner aperture(Angstrom-1)
	double * __restrict g2max;		// Outer aperture(Angstrom-1)
} sDetCir;

typedef struct sDetInt{
	double * __restrict Tot;		// Total Intensity
	double * __restrict Coh;		// Coherent Intensity
} sDetInt;

typedef struct sImSTEM{
	sDetInt *DetInt;
} sImSTEM;

/************************Input - TEM-Image***********************/
typedef struct sInTEMIm{
	int gpu;				// gpu device
	int MEffect;			// 1: Partial coherente mode, 2: Transmission cross coefficient
	int STEffect;			// 1: Spatial and temporal, 2: Temporal, 3: Spatial
	int ZeroDef;			// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
	double ZeroDefPlane;	// Zero defocus plane

	double E0;			// Acceleration volatage in KeV
	double *Psirh;		// Real part of the wavefunction
	double *Psiih;		// Imaginary part of the wavefunction
	int nx;				// Number of pixels in x direction
	int ny;				// Number of pixels in y direction
	double lx;			// distance in x direction(Angstroms)
	double ly;			// distance in y direction(Angstroms)
	
	int MC_m;			// momentum of the vortex
	double MC_f;		// defocus(Angstrom)
	double MC_Cs3;		// spherical aberration(Angstrom)
	double MC_Cs5;		// spherical aberration(Angstrom)
	double MC_mfa2;		// magnitude 2-fold astigmatism(Angstrom)
	double MC_afa2;		// angle 2-fold astigmatism(rad)
	double MC_mfa3;		// magnitude 3-fold astigmatism(Angstrom)
	double MC_afa3;		// angle 3-fold astigmatism(rad)
	double MC_aobjl;	// lower objective aperture(rad)
	double MC_aobju;	// upper objective aperture(rad)
	double MC_sf;		// defocus spread(Angstrom)
	int MC_nsf;			// Number of defocus sampling points
	double MC_beta;		// semi-convergence angle
	int MC_nbeta;			// half number sampling points
} sInTEMIm;

/*************************Input - TEM****************************/
typedef struct sInMSTEM{
	int gpu;					// gpu device
	int SimType;				// 1: STEM, 2: CBED, 3: HRTEM, 4: ED, 5: PED, 6: HCI, ... 10: EW real, 11: EW Fourier
	int nConfFP;				// Number of frozen phonon configurations
	int DimFP;					// Dimensions phonon configurations
	int SeedFP;					// Random seed(frozen phonon)
	int PotPar;					// Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
	int MEffect;				// 1: Partial coherente mode, 2: Transmission cross coefficient
	int STEffect;				// 1: Spatial and temporal, 2: Temporal, 3: Spatial
	int ZeroDefTyp;				// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
	double ZeroDefPlane;		// Zero defocus plane
	int ApproxModel;			// 1: MS, 2: PA, 3: POA, 4:WPOA
	int BandwidthLimit;			// 1: true, 2: false
	int ThicknessTyp;			// 1: Whole specimen, 2: Throught thickness, 3: Through planes
	int nThickness;				// Number of thickness
	double* Thickness;			// Array of thicknesses

	double E0;					// Acceleration volatage in KeV
	double theta;				// incident tilt (in spherical coordinates) (rad)
	double phi;					// incident tilt (in spherical coordinates) (rad)
	int nx;						// Number of pixels in x direction
	int ny;						// Number of pixels in y direction
	double lx;					// distance in x direction(Angstroms)
	double ly;					// distance in y direction(Angstroms)
	double dz;					// slice thickness

	int MC_m;					// momentum of the vortex
	double MC_f;				// defocus(Angstrom)
	double MC_Cs3;				// spherical aberration(Angstrom)
	double MC_Cs5;				// spherical aberration(Angstrom)
	double MC_mfa2;				// magnitude 2-fold astigmatism(Angstrom)
	double MC_afa2;				// angle 2-fold astigmatism(rad)
	double MC_mfa3;				// magnitude 3-fold astigmatism(Angstrom)
	double MC_afa3;				// angle 3-fold astigmatism(rad)
	double MC_aobjl;			// lower objective aperture(rad)
	double MC_aobju;			// upper objective aperture(rad)
	double MC_sf;				// defocus spread(Angstrom)
	int MC_nsf;					// Number of defocus sampling points
	double MC_beta;				// semi-convergence angle
	int MC_nbeta;					// half number sampling points

	int nAtomsM;				// Number of Atoms
	double *AtomsM;				// Atoms in a matrix form

	int STEM_line;				// 0: Area, 1: Line
	bool STEM_FastCal;			// 0: normal mode(low memory consumption), 1: fast calculation(high memory consumption)
	int STEM_ns;				// Sampling points
	double STEM_x1u;			// Initial scanning position in x
	double STEM_y1u;			// Initial scanning in y
	double STEM_x2u;			// final scanning position in x
	double STEM_y2u;			// final scanning position in y
	int STEM_nDet;				// Number of circular detectors
	sInDetCir *STEM_DetCir;		// Circular detectors

	double CBED_x0;				// x position
	double CBED_y0;				// y position
	int CBED_Space;			// y position

	//int HRTEM_xx;

	int PED_nrot;				// Number of orientations
	double PED_theta;			// Precession angle

	int HCI_nrot;				// Number of orientations
	double HCI_theta;			// Precession angle
	//int xx;

	//int EWRS_xx;

	//int EWFS_xx;
} sInMSTEM;

/*************************Input - TEM****************************/
typedef struct sInProbe{
	int gpu;				// gpu device

	double E0;				// Acceleration volatage in KeV
	double theta;			// incident tilt (in spherical coordinates) (rad)
	double phi;				// incident tilt (in spherical coordinates) (rad)
	int nx;					// Number of pixels in x direction
	int ny;					// Number of pixels in y direction
	double lx;				// distance in x direction(Angstroms)
	double ly;				// distance in y direction(Angstroms)

	double x0;				// x position
	double y0;				// y position

	int m;					// momentum of the vortex
	double f;				// defocus(Angstrom)
	double Cs3;				// spherical aberration(Angstrom)
	double Cs5;				// spherical aberration(Angstrom)
	double mfa2;			// magnitude 2-fold astigmatism(Angstrom)
	double afa2;			// angle 2-fold astigmatism(rad)
	double mfa3;			// magnitude 3-fold astigmatism(Angstrom)
	double afa3;			// angle 3-fold astigmatism(rad)
	double aobjl;			// lower objective aperture(rad)
	double aobju;			// upper objective aperture(rad)
	double sf;				// defocus spread(Angstrom)
	int nsf;				// Number of defocus sampling points
	double beta;			// semi-convergence angle
	int nbeta;				// half number sampling points

} sInProbe;

/******************Radial Schrodinger equation*******************/
typedef struct sInRadSchr{
	double E0;					// Acceleration Voltage
	int PotPar;					// Parameterization type
	int n;						// Principal quantum number
	int nr;						// Number of grid points
	int nAtomsM;				// Number of Atoms
	double *AtomsM;				// Atoms
} sInRadSchr;

/****************Multislice general parameters******************/
class cMGP{
	private:
		void freeMemory(){
			SimType = 0;
			MulOrder = 0;
			nConfFP = 0;
			DimFP = 0;
			DistFP = 0;
			SeedFP = 0;
			PotPar = 0;
			MEffect = 0;
			STEffect = 0;
			ZeroDefTyp = 0;
			ZeroDefPlane = 0;
			ApproxModel = 0;
			ThkTyp = 0;
			nThk = 0;
			delete [] Thk; Thk = 0;
			BWL = false;
			PBC_xy = false;
			Vrl = 0;
			E0 = 0;
			theta = 0;	
			phi = 0;
			lx = 0;
			ly = 0;	
			dz = 0;
			nx = 0;
			ny = 0;
		}
	public:
		int gpu;					// gpu device
		int SimType;				// 1: STEM, 2: CBED, 3: HRTEM, 4: ED, 5: PED, 6: HCI, ... 10: EW real, 11: EW Fourier
		int MulOrder;				// 1: First order MS, 2: Second Order MS
		int nConfFP;				// Number of frozen phonon configurations
		int DimFP;					// Dimensions phonon configurations
		int DistFP;					// Frozen phonon distribution type 1:Normal, 2:xx
		int SeedFP;					// Random seed(frozen phonon)
		int PotPar;					// Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) and 6: Lobato(0-12)
		int MEffect;				// 1: Partial coherente mode, 2: Transmission cross coefficient
		int STEffect;				// 1: Spatial and temporal, 2: Temporal, 3: Spatial
		int ZeroDefTyp;				// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
		double ZeroDefPlane;		// Zero defocus plane
		int ApproxModel;			// 1: Mulstilice, 2: Projection approximation, 3: Phase object approximation, 4: Weak phase object approximation
		int ThkTyp;					// 1: Whole specimen, 2: Throught thickness, 3: Through planes
		int nThk;					// Number of thickness
		double* Thk;				// Array of thicknesses	
		bool BWL;					// 1: true, 2: false
		bool PBC_xy;				// Peridic boundary contions	
		double Vrl;					// Atomic potential cut-off
		double E0;					// Acceleration volatage in KeV
		double theta;				// Tilt (in spherical coordinates) (rad)
		double phi;					// Tilt (in spherical coordinates) (rad)
		double lx;					// Box size in x direction(Angstroms)
		double ly;					// Box size in y direction(Angstroms)
		double dz;					// slice thickness
		int nx;						// Number of pixels in x direction
		int ny;						// Number of pixels in y direction

		cMGP(){
			gpu = 0;
			SimType = 10;
			MulOrder = 2;
			nConfFP = 0;
			DimFP = 111;
			DistFP = 1;
			SeedFP = 1983;
			PotPar = 6;
			MEffect = 1;
			STEffect = 1;
			ZeroDefTyp = 3;
			ZeroDefPlane = 0.0;
			ApproxModel = 1;
			ThkTyp = 1;
			nThk = 0;
			Thk = 0;
			BWL = true;
			PBC_xy = true;
			Vrl = stVrl;
			E0 = 300;
			theta = 0;	
			phi = 0;
			lx = 4.078;
			ly = 4.078;	
			dz = 0.25;
			nx = 1024;
			ny = 1024;
		}

		~cMGP(){
			freeMemory();
		}

		cMGP& operator= (const cMGP &MGP){
			freeMemory();

			gpu = MGP.gpu;
			SimType = MGP.SimType;
			MulOrder = MGP.MulOrder;
			nConfFP = MGP.nConfFP;
			DimFP = MGP.DimFP;
			DistFP = MGP.DistFP;
			SeedFP = MGP.SeedFP;
			PotPar = MGP.PotPar;
			MEffect = MGP.MEffect;
			STEffect = MGP.STEffect;
			ZeroDefTyp = MGP.ZeroDefTyp;
			ZeroDefPlane = MGP.ZeroDefPlane;
			ApproxModel = MGP.ApproxModel;
			ThkTyp = MGP.ThkTyp;
			nThk = MGP.nThk;
			Thk = new double[nThk];
			memcpy(Thk, MGP.Thk, nThk*cSizeofRD);
			BWL = MGP.BWL;
			PBC_xy = MGP.PBC_xy;
			Vrl = MGP.Vrl;
			E0 = MGP.E0;
			theta = MGP.theta;	
			phi = MGP.phi;
			lx = MGP.lx;
			ly = MGP.ly;	
			dz = MGP.dz;
			nx = MGP.nx;
			ny = MGP.ny;

			return *this;
		}

		void SetInputData(sInMSTEM &InMSTEM){
			freeMemory();

			gpu = InMSTEM.gpu;
			SimType = InMSTEM.SimType;	
			MulOrder = 2;
			nConfFP = InMSTEM.nConfFP;		
			DimFP = InMSTEM.DimFP;	
			DistFP = 1;
			SeedFP = InMSTEM.SeedFP;
			PotPar = InMSTEM.PotPar;
			MEffect = InMSTEM.MEffect;
			STEffect = InMSTEM.STEffect;
			ZeroDefTyp = InMSTEM.ZeroDefTyp;
			ZeroDefPlane = InMSTEM.ZeroDefPlane;
			ApproxModel = InMSTEM.ApproxModel;
			ThkTyp = InMSTEM.ThicknessTyp;
			nThk = InMSTEM.nThickness;
			memcpy(Thk, InMSTEM.Thickness, nThk*cSizeofRD); // change
			BWL = (InMSTEM.BandwidthLimit==1)?true:false;
			Vrl = stVrl;
			E0 = InMSTEM.E0;	
			theta = InMSTEM.theta;	
			phi = InMSTEM.phi;
			lx = InMSTEM.lx;
			ly = InMSTEM.ly;
			dz = InMSTEM.dz;
			nx = InMSTEM.nx;
			ny = InMSTEM.ny;
			if(ApproxModel>1){
				MulOrder = 1;
				DimFP = DimFP - DimFP%10;
			}
		}
};

#endif