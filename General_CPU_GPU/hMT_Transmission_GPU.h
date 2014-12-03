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

#ifndef hMT_Transmission_GPU_H
#define hMT_Transmission_GPU_H

#include "hConstTypes.h"
#include "hQuadrature.h"
#include "hMT_General_CPU.h"
#include "hMT_General_GPU.h"
#include "hMT_MGP_CPU.h"
#include "hMT_Specimen_CPU.h"
#include "hMT_AtomTypes_GPU.h"
#include "hMT_Potential_GPU.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

class cMT_Transmission_GPU: public cMT_Potential_GPU{
	private:
		int cSynCPU;
		int nSliceM;
		int SliceMTyp;
		double2 **Trans;
		float **Vpe;
		cufftHandle PlanPsi;
		void Potential_Efective(int iSlice, double fPot, float *&Vpe);
		void Transmission(int iSlice, double2 *&Trans);
	public:

		void freeMemory();
		cMT_Transmission_GPU();
		~cMT_Transmission_GPU();

		void SetInputData(cMT_MGP_CPU &MT_MGP_CPU_io, int nAtomsM_i, double *AtomsM_i);
		void getTrans(int iSlice);

};

#endif