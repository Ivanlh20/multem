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

#include "hConstTypes.h"
#include "hMT_General_CPU.h"
#include "hMT_General_GPU.h"
#include "hMT_IncidentWave_GPU.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>


void cMT_IncidentWave_GPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	f_sGP_Init(GP);
	f_sLens_Init(Lens);

	cudaFreen(Mp_d);
	PlanPsi = 0;
}

cMT_IncidentWave_GPU::cMT_IncidentWave_GPU(){
	f_sGP_Init(GP);
	f_sLens_Init(Lens);

	Mp_d = 0;
	PlanPsi = 0;
}

cMT_IncidentWave_GPU::~cMT_IncidentWave_GPU(){
	freeMemory();
}

void cMT_IncidentWave_GPU::SetInputData(sGP &GP_i, sLens &Lens_i, cufftHandle &PlanPsi_i){
	freeMemory();

	GP = GP_i;
	Lens = Lens_i;
	PlanPsi = PlanPsi_i;
	cudaMalloc((void**)&Mp_d, 32*32*cSizeofRD);
}

void cMT_IncidentWave_GPU::Psi0(double2 *&Psi0){
	f_Set_MC(GP, 1.0, 0.0, Psi0);
}

void cMT_IncidentWave_GPU::Psi0(double x, double y, double2 *&Psi0){
	f_Probe_FS(GP, Lens, c2Pi*(0.5*GP.lx-x), c2Pi*(0.5*GP.ly-y), Psi0);
	cufftExecZ2Z(PlanPsi, Psi0, Psi0, CUFFT_INVERSE);
	double Totalsum = f_Sum_MC2(GP, Psi0, Mp_d);
	f_Scale_MC(GP, sqrt(double(GP.nxy)/Totalsum), Psi0);
}