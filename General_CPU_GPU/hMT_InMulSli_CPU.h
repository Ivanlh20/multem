/*
 * This file is part of MULTEM.
 * Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * MULTEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MULTEM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef hMT_InMulSli_CPU_H
#define hMT_InMulSli_CPU_H

#include <cstring>
#include "hConstTypes.h"

class cMT_InMulSli_CPU{
	public:
		int gpu;					// gpu device
		int SimType;				// 11: STEM, 12: ISTEM, 21: CBED, 22: CBEI, 31: ED, 32: HRTEM, 41: PED, 42: HCI, ... 51: EW Fourier, 52: EW real
		int MulOrder;				// 1: First order, 2: Second order
		int nConfFP;				// Number of frozen phonon configurations
		int DimFP;					// Dimensions phonon configurations
		int DistFP;					// 1: Gaussian (Phonon distribution)
		int SeedFP;					// Random seed(frozen phonon)
		int PotPar;					// Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
		int MEffect;				// 1: Partial coherente mode, 2: Transmission cross coefficient
		int STEffect;				// 1: Spatial and temporal, 2: Temporal, 3: Spatial
		int ZeroDefTyp;				// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
		double ZeroDefPlane;		// Zero defocus plane
		int ApproxModel;			// 1: Mulstilice, 2: Projection approximation, 3: Phase object approximation, 4: Weak phase object approximation
		int BWL;					// 1: true, 2: false (bandwidth limited)
		int FastCal;				// 1: normal mode(low memory consumption), 2: fast calculation(high memory consumption)
		int PBC_xy;					// 1: true, 2: false (Peridic boundary contions)
		int ThkTyp;					// 1: Whole specimen, 2: Throught thickness, 3: Through planes
		int nThk;					// Number of thickness
		double *Thk;				// Array of thicknesses
		int Psi0Typ;				// 1: Automatic, 2: User define
		sComplex Psi0; // Input wave

		double E0;					// Acceleration volatage in KeV
		double theta;				// incident tilt (in spherical coordinates) (rad)
		double phi;					// incident tilt (in spherical coordinates) (rad)
		int nx;						// Number of pixels in x direction
		int ny;						// Number of pixels in y direction
		double lx;					// Distance in x direction(Angstroms)
		double ly;					// Distance in y direction(Angstroms)
		double dz;					// Slice thickness

		int MC_m;					// Momentum of the vortex
		double MC_f;				// Defocus(Angstrom)
		double MC_Cs3;				// Spherical aberration(Angstrom)
		double MC_Cs5;				// Spherical aberration(Angstrom)
		double MC_mfa2;				// Magnitude 2-fold astigmatism(Angstrom)
		double MC_afa2;				// Angle 2-fold astigmatism(rad)
		double MC_mfa3;				// Magnitude 3-fold astigmatism(Angstrom)
		double MC_afa3;				// Angle 3-fold astigmatism(rad)
		double MC_aobjl;			// Lower objective aperture(rad)
		double MC_aobju;			// Upper objective aperture(rad)
		double MC_sf;				// Defocus spread(Angstrom)
		int MC_nsf;					// Number of defocus sampling points
		double MC_beta;				// Semi-convergence angle
		int MC_nbeta;				// half number sampling points

		int nAtomsM;				// Number of Atoms
		double *AtomsM;				// Atoms in a matrix form

		int STEM_line;				// 0: Area, 1: Line
		int STEM_ns;				// Sampling points
		double STEM_x1u;			// Initial scanning position in x
		double STEM_y1u;			// Initial scanning in y
		double STEM_x2u;			// final scanning position in x
		double STEM_y2u;			// final scanning position in y
		int STEM_nDet;				// Number of circular detectors
		sInDetCir *STEM_DetCir;		// Circular detectors

		double CBED_x0;				// x position
		double CBED_y0;				// y position

		double CBEI_x0;				// x position
		double CBEI_y0;				// y position

		int PED_nrot;				// Number of orientations
		double PED_theta;			// Precession angle

		int HCI_nrot;				// Number of orientations
		double HCI_theta;			// Precession angle


		void freeMemory();
		void Default_Values();
		cMT_InMulSli_CPU();
		~cMT_InMulSli_CPU();
};

inline void cMT_InMulSli_CPU::freeMemory()
{
	gpu = 0;
	SimType = 0;
	nConfFP = 0;
	DimFP = 0;
	SeedFP = 0;
	PotPar = 0;
	MEffect = 0;		
	STEffect = 0;	
	ZeroDefTyp = 0;
	ZeroDefPlane = 0;
	ApproxModel = 0;
	BWL = 1;
	FastCal = 0;
	ThkTyp = 0;
	nThk = 0;
	delete [] Thk; Thk = 0;
	Psi0Typ = 1;
	Psi0.real = 0;
	Psi0.imag = 0;

	E0 = 0;
	theta = 0;
	phi = 0;
	nx = 0;
	ny = 0;
	lx = 0;
	ly = 0;
	dz = 0;

	MC_m = 0;
	MC_f = 0;
	MC_Cs3 = 0;
	MC_Cs5 = 0;
	MC_mfa2 = 0;
	MC_afa2 = 0;
	MC_mfa3 = 0;
	MC_afa3 = 0;
	MC_aobjl = 0;
	MC_aobju = 0;
	MC_sf = 0;
	MC_nsf = 0;
	MC_beta = 0;
	MC_nbeta = 0;

	nAtomsM = 0;
	AtomsM = 0;

	STEM_line = 0;
	STEM_ns = 0;
	STEM_x1u = 0;
	STEM_y1u = 0;
	STEM_x2u = 0;
	STEM_y2u = 0;
	
	STEM_nDet = 0;
	delete [] STEM_DetCir; STEM_DetCir = 0;

	CBED_x0 = 0;
	CBED_y0 = 0;

	CBEI_x0 = 0;
	CBEI_y0 = 0;

	PED_nrot = 0;
	PED_theta = 0;

	HCI_nrot = 0;
	HCI_theta = 0;
}

inline cMT_InMulSli_CPU::cMT_InMulSli_CPU()
{
	gpu = 0;
	SimType = 0;
	MulOrder = 2;
	nConfFP = 0;
	DimFP = 0;
	SeedFP = 0;
	PotPar = 0;
	MEffect = 0;		
	STEffect = 0;	
	ZeroDefTyp = 0;
	ZeroDefPlane = 0;
	ApproxModel = 0;
	BWL = 0;
	FastCal = 0;
	PBC_xy = 0;
	ThkTyp = 0;
	nThk = 0;
	Thk = 0;
	Psi0Typ = 1;
	Psi0.real = 0;
	Psi0.imag = 0;

	E0 = 0;
	theta = 0;
	phi = 0;
	nx = 0;
	ny = 0;
	lx = 0;
	ly = 0;
	dz = 0;

	MC_m = 0;
	MC_f = 0;
	MC_Cs3 = 0;
	MC_Cs5 = 0;
	MC_mfa2 = 0;
	MC_afa2 = 0;
	MC_mfa3 = 0;
	MC_afa3 = 0;
	MC_aobjl = 0;
	MC_aobju = 0;
	MC_sf = 0;
	MC_nsf = 0;
	MC_beta = 0;
	MC_nbeta = 0;

	nAtomsM = 0;
	AtomsM = 0;

	STEM_line = 0;
	STEM_ns = 0;
	STEM_x1u = 0;
	STEM_y1u = 0;
	STEM_x2u = 0;
	STEM_y2u = 0;
	
	STEM_nDet = 0;
	STEM_DetCir = 0;

	CBED_x0 = 0;
	CBED_y0 = 0;

	CBEI_x0 = 0;
	CBEI_y0 = 0;

	PED_nrot = 0;
	PED_theta = 0;

	HCI_nrot = 0;
	HCI_theta = 0;
}

inline cMT_InMulSli_CPU::~cMT_InMulSli_CPU()
{
	freeMemory();
}

inline void cMT_InMulSli_CPU::Default_Values()
{
	gpu = 0;
 SimType = 32;
	MulOrder = 2;	
	nConfFP = 0;
	DimFP = 111;
	DistFP = 1;
	SeedFP = 1983;
	PotPar = 6;
	MEffect = 1;		
	STEffect = 1;	
	ZeroDefTyp = 3;
	ZeroDefPlane = 0;
	ApproxModel = 2;
	BWL = 1;
	FastCal = 1;
	PBC_xy = 1;
	ThkTyp = 1;
	nThk = 1;
	Thk = new double[nThk];
	Psi0Typ = 1;
	Psi0.real = 0;
	Psi0.imag = 0;

	E0 = 200;
	theta = 0;
	phi = 0;
	nx = 2048;
	ny = 2048;
	lx = 54.307;
	ly = 54.307;
	dz = 1.0195;

	MC_m = 0;
	MC_f = 1110;
	MC_Cs3 = 1.0*mm2Ags;
	MC_Cs5 = 0;
	MC_mfa2 = 0;
	MC_afa2 = 0;
	MC_mfa3 = 0;
	MC_afa3 = 0;
	MC_aobjl = 0;
	MC_aobju = 1000*mrad2rad;
	MC_sf = 88;
	MC_nsf = 10;
	MC_beta = 0.2;
	MC_nbeta = 10;

	nAtomsM = 0;
	AtomsM = 0;

	STEM_line = 1;
	STEM_ns = 10;
	STEM_x1u = 0;
	STEM_y1u = 0;
	STEM_x2u = 1.0;
	STEM_y2u = 1.0;
	
	STEM_nDet = 1;
	STEM_DetCir = 0;

	CBED_x0 = 0;
	CBED_y0 = 0;

	CBEI_x0 = 0;
	CBEI_y0 = 0;

	PED_nrot = 0;
	PED_theta = 0;

	HCI_nrot = 0;
	HCI_theta = 0;
}

#endif
