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

#ifndef hMT_Detector_CPU_H
#define hMT_Detector_CPU_H

#include "hConstTypes.h"


/**************************STEM*************************/
class cMT_Detector_CPU{
	private:
		sGP GP;

		int nDet;			// number of Detectors
		sDetCir DetCir;		// Detectors

		double *Tot;		// Total Intensity
		double *Coh;		// Coherent Intensity
	public:
		void freeMemory();
		cMT_Detector_CPU();
		~cMT_Detector_CPU();

		void SetInputData(sGP &GP_i, int nDeti, sDetCir &DetCiri);
		void getDetectorIntensity(double *&Toth, double *&Cohh, int ixys, sDetInt *DetInth);
};

#endif