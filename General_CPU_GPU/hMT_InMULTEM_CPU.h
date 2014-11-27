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

#ifndef hMT_InMULTEM_CPU_H
#define hMT_InMULTEM_CPU_H

#include <cstring>
#include "hConstTypes.h"

class cMT_InMULTEM_CPU{
	public:
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

		double CBEI_x0;				// x position
		double CBEI_y0;				// y position

		int PED_nrot;				// Number of orientations
		double PED_theta;			// Precession angle

		int HCI_nrot;				// Number of orientations
		double HCI_theta;			// Precession angle


		void freeMemory();
		cMT_InMULTEM_CPU();
		~cMT_InMULTEM_CPU();
};

#endif