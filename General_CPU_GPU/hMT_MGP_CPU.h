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

#ifndef hMT_MGP_CPU_H
#define hMT_MGP_CPU_H

#include "hConstTypes.h"
#include "hMT_InMulSli_CPU.h"

class cMT_MGP_CPU{
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

		void freeMemory();
		cMT_MGP_CPU();
		~cMT_MGP_CPU();

		cMT_MGP_CPU& operator= (const cMT_MGP_CPU &MT_MGP_CPU);
		void SetInputData(cMT_InMulSli_CPU &MT_InMulSli_CPU);
};

inline void cMT_MGP_CPU::freeMemory(){
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

inline cMT_MGP_CPU::cMT_MGP_CPU(){
	gpu = 0;
	SimType = 52;
	MulOrder = 2;
	nConfFP = 0;
	DimFP = 110;
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

inline cMT_MGP_CPU::~cMT_MGP_CPU(){
	freeMemory();
}

inline cMT_MGP_CPU& cMT_MGP_CPU::operator= (const cMT_MGP_CPU &MT_MGP_CPU){
	freeMemory();

	gpu = MT_MGP_CPU.gpu;
	SimType = MT_MGP_CPU.SimType;
	MulOrder = MT_MGP_CPU.MulOrder;
	nConfFP = MT_MGP_CPU.nConfFP;
	DimFP = MT_MGP_CPU.DimFP;
	DistFP = MT_MGP_CPU.DistFP;
	SeedFP = MT_MGP_CPU.SeedFP;
	PotPar = MT_MGP_CPU.PotPar;
	MEffect = MT_MGP_CPU.MEffect;
	STEffect = MT_MGP_CPU.STEffect;
	ZeroDefTyp = MT_MGP_CPU.ZeroDefTyp;
	ZeroDefPlane = MT_MGP_CPU.ZeroDefPlane;
	ApproxModel = MT_MGP_CPU.ApproxModel;
	ThkTyp = MT_MGP_CPU.ThkTyp;
	nThk = MT_MGP_CPU.nThk;
	if(nThk>0){
		Thk = new double[nThk];
		memcpy(Thk, MT_MGP_CPU.Thk, nThk*cSizeofRD);
	}
	BWL = MT_MGP_CPU.BWL;
	PBC_xy = MT_MGP_CPU.PBC_xy;
	Vrl = MT_MGP_CPU.Vrl;
	E0 = MT_MGP_CPU.E0;
	theta = MT_MGP_CPU.theta;	
	phi = MT_MGP_CPU.phi;
	lx = MT_MGP_CPU.lx;
	ly = MT_MGP_CPU.ly;	
	dz = MT_MGP_CPU.dz;
	nx = MT_MGP_CPU.nx;
	ny = MT_MGP_CPU.ny;

	return *this;
}

inline void cMT_MGP_CPU::SetInputData(cMT_InMulSli_CPU &MT_InMulSli_CPU){
	freeMemory();

	gpu = MT_InMulSli_CPU.gpu;
	SimType = MT_InMulSli_CPU.SimType;	
	MulOrder = 2;
	nConfFP = MT_InMulSli_CPU.nConfFP;		
	DimFP = MT_InMulSli_CPU.DimFP;	
	DistFP = 1;
	SeedFP = MT_InMulSli_CPU.SeedFP;
	PotPar = MT_InMulSli_CPU.PotPar;
	MEffect = MT_InMulSli_CPU.MEffect;
	STEffect = MT_InMulSli_CPU.STEffect;
	ZeroDefTyp = MT_InMulSli_CPU.ZeroDefTyp;
	ZeroDefPlane = MT_InMulSli_CPU.ZeroDefPlane;
	ApproxModel = MT_InMulSli_CPU.ApproxModel;
	ThkTyp = MT_InMulSli_CPU.ThicknessTyp;
	nThk = (ThkTyp==1)?1:MT_InMulSli_CPU.nThickness;
	Thk = new double[nThk];
	memcpy(Thk, MT_InMulSli_CPU.Thickness, nThk*cSizeofRD); // change
	BWL = (MT_InMulSli_CPU.BandwidthLimit==1)?true:false;
	Vrl = stVrl;
	E0 = MT_InMulSli_CPU.E0;	
	theta = MT_InMulSli_CPU.theta;	
	phi = MT_InMulSli_CPU.phi;
	lx = MT_InMulSli_CPU.lx;
	ly = MT_InMulSli_CPU.ly;
	dz = MT_InMulSli_CPU.dz;
	nx = MT_InMulSli_CPU.nx;
	ny = MT_InMulSli_CPU.ny;
	if(ApproxModel>1){
		MulOrder = 1;
		DimFP = DimFP - DimFP%10;
	}
}

#endif