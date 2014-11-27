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
#include "hMT_InMULTEM_CPU.h"

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
		void SetInputData(cMT_InMULTEM_CPU &MT_InMULTEM_CPU);
};

#endif