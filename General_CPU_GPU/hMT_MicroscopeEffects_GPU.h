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

#ifndef hMT_MicroscopeEffects_GPU_H
#define hMT_MicroscopeEffects_GPU_H

#include "hConstTypes.h"
#include "hQuadrature.h"
#include "hMT_General_CPU.h"
#include "hMT_General_GPU.h"

#include "math.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

/*************************Incident wave*************************/
class cMT_MicroscopeEffects_GPU{
	private:
		int IdCall;
		int cSynCPU;

		sQ1 Qt;
		int nQs;
		sQ2 Qs;
		sGP GP;
		sLens Lens;
		cufftHandle PlanPsi;
		double2 *Psit;

		void ReadTemporalQuadrature(sLens &Lens, sQ1 &Qt);
		void ReadSpatialQuadrature(sLens &Lens, int &nQs, sQ2 &Qs);
		void PCTCCTEM(int STEffect, double2 *&fPsi, double *&M2PsiM);
		void PCLIMWPOTEM(int STEffect, double2 *&fPsi, double *&M2PsiM);
	public:
		void freeMemory();
		cMT_MicroscopeEffects_GPU();
		~cMT_MicroscopeEffects_GPU();

		void SetInputData(sGP &GP_i, sLens &Lens_i, cufftHandle &PlanPsi_i, double2 *&Psit_i);
		void ApplyMEffects(int MEffect, int STEffect, double2 *&fPsi, double *&M2Psi);
};

#endif