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

#ifndef hMT_Detector_GPU_H
#define hMT_Detector_GPU_H

#include "hConstTypes.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>

/*************************STEM*************************/
class cMT_Detector_GPU{
	private:
		sGP GP;

		double *M1p_d;
		double *M2p_d;

		int nDet;			// number of Detectors
		sDetCir DetCirh;	// Detectors

		double *Tot_h;		// Total Intensity
		double *Coh_h;		// Coherent Intensity

		double *Tot_d;		// Total Intensity
		double *Coh_d;		// Coherent Intensity
	public:
		void freeMemory();
		cMT_Detector_GPU();
		~cMT_Detector_GPU();

		void SetInputData(sGP &GP_i, int nDeti, sDetCir &DetCirhi);
		void getDetectorIntensity(double *&aM2Psi, double *&M2aPsi, int ixys, sDetInt *DetInth);
		void getDetectorIntensity(double *&aM2Psi, int ixys, sDetInt *DetInth, bool add=false);
		void getDetectorIntensity(double2 *&aPsi, int ixys, sDetInt *DetInth, bool add=false);
};

#endif