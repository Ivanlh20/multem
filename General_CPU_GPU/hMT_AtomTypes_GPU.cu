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
#include "hMT_General_GPU.h"
#include "hMT_AtomTypes_CPU.h"
#include "hMT_AtomTypes_GPU.h"

#include "math.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>

// free memory
void cMT_AtomTypes_GPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	Z = 0;
	m = 0;
	A = 0;
	rn_e = 0;
	rn_c = 0;
	ra_e = 0;
	ra_c = 0;
	Rmin = 0;
	Rmax = 0;
	Rmin2 = 0;
	Rmax2 = 0;

	f_sCoefPar_cudaFree(cfeg);
	f_sCoefPar_cudaFree(cfxg);
	f_sCoefPar_cudaFree(cPr);
	f_sCoefPar_cudaFree(cVr);
	f_sCoefPar_cudaFree(cVR);

	nR = 0;
	cudaFreen(R);
	cudaFreen(R2);
	f_sciVn_cudaFree(ciVR);
};

// Set Atom type
void cMT_AtomTypes_GPU::SetAtomTypes(cMT_AtomTypes_CPU &MT_AtomTypes_CPU_i){
	freeMemory(); // clean GPU memory

	Z = MT_AtomTypes_CPU_i.Z;
	m = MT_AtomTypes_CPU_i.m;
	A = MT_AtomTypes_CPU_i.A;
	rn_e = MT_AtomTypes_CPU_i.rn_e;
	rn_c = MT_AtomTypes_CPU_i.rn_c;
	ra_e = MT_AtomTypes_CPU_i.ra_e;
	ra_c = MT_AtomTypes_CPU_i.ra_c;
	Rmin = MT_AtomTypes_CPU_i.Rmin;
	Rmin2 = MT_AtomTypes_CPU_i.Rmin2;
	Rmax = MT_AtomTypes_CPU_i.Rmax;
	Rmax2 = MT_AtomTypes_CPU_i.Rmax2;

	f_sCoefPar_cudaMalloc(6, cfeg);
	cudaMemcpy(cfeg.cl, MT_AtomTypes_CPU_i.cfeg.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(cfeg.cnl, MT_AtomTypes_CPU_i.cfeg.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	f_sCoefPar_cudaMalloc(6, cfxg);
	cudaMemcpy(cfxg.cl, MT_AtomTypes_CPU_i.cfxg.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(cfxg.cnl, MT_AtomTypes_CPU_i.cfxg.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	f_sCoefPar_cudaMalloc(6, cPr);
	cudaMemcpy(cPr.cl, MT_AtomTypes_CPU_i.cPr.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(cPr.cnl, MT_AtomTypes_CPU_i.cPr.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	f_sCoefPar_cudaMalloc(6, cVr);
	cudaMemcpy(cVr.cl, MT_AtomTypes_CPU_i.cVr.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(cVr.cnl, MT_AtomTypes_CPU_i.cVr.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	f_sCoefPar_cudaMalloc(6, cVR);
	cudaMemcpy(cVR.cl, MT_AtomTypes_CPU_i.cVR.cl, 6*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(cVR.cnl, MT_AtomTypes_CPU_i.cVR.cnl, 6*cSizeofRD, cudaMemcpyHostToDevice);

	nR = MT_AtomTypes_CPU_i.nR;

	cudaMalloc((void**)&R, nR*cSizeofRD);
	cudaMemcpy(R, MT_AtomTypes_CPU_i.R, nR*cSizeofRD, cudaMemcpyHostToDevice);

	cudaMalloc((void**)&R2, nR*cSizeofRD);
	cudaMemcpy(R2, MT_AtomTypes_CPU_i.R2, nR*cSizeofRD, cudaMemcpyHostToDevice);

	f_sciVn_cudaMalloc(nR, ciVR);
	cudaMemcpy(ciVR.c0, MT_AtomTypes_CPU_i.ciVR.c0, nR*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(ciVR.c1, MT_AtomTypes_CPU_i.ciVR.c1, nR*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(ciVR.c2, MT_AtomTypes_CPU_i.ciVR.c2, nR*cSizeofRD, cudaMemcpyHostToDevice);
	cudaMemcpy(ciVR.c3, MT_AtomTypes_CPU_i.ciVR.c3, nR*cSizeofRD, cudaMemcpyHostToDevice);
}