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

#ifndef hMT_AtomTypes_GPU_H
#define hMT_AtomTypes_GPU_H

#include "hConstTypes.h"

class cMT_AtomTypes_GPU{
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
		float *R2f;		// R2 Grid
		sciVn ciVR;		// Look up table - Projected potential coefficients

		void freeMemory();
		cMT_AtomTypes_GPU();
		~cMT_AtomTypes_GPU();
		void SetAtomTypes(int PotPari, sAtomTypesCPU &AtomTypesCPUi, int nRi, double dRmini);
};

#endif