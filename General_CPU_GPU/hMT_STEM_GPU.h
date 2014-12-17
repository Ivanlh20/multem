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

#ifndef hMT_STEM_GPU_H
#define hMT_STEM_GPU_H

#include "hConstTypes.h"
#include "hMT_General_CPU.h"
#include "hMT_General_GPU.h"
#include "hMT_Specimen_CPU.h"
#include "hMT_Potential_GPU.h"
#include "hMT_IncidentWave_GPU.h"
#include "hMT_Detector_CPU.h"
#include "hMT_Detector_GPU.h"
#include "hMT_MulSli_GPU.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

/*****************************STEM********************************/
class cMT_STEM_GPU{
	public:
		int line;						// 0: Area, 1: Line
		int ns;							// Sampling points
		double x1u;						// Initial scanning position in x
		double y1u;						// Initial scanning in y
		double x2u;						// final scanning position in x
		double y2u;						// final scanning position in y
		int nThk;

		int nDet;						// Number of circular detectors
		sDetCir DetCir;					// Circular detectors
		sImSTEM *ImSTEM;
		int nxs;
		int nys;
		double *xs;
		double *ys;

		int nst;
		double *xst;
		double *yst;

		cMT_Detector_GPU *MT_Detector_GPU;

		void freeMemory();
		cMT_STEM_GPU();
		~cMT_STEM_GPU();

		void SetInputData(cMT_InMulSli_CPU *MT_InMulSli_CPU_i, cMT_MGP_CPU *MT_MGP_CPU_i, int nThk_i);
		void InitImSTEM();

};

#endif