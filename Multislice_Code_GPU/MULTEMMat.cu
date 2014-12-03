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
#include "hMatlab2Cpp.h"
#include "hMT_InMulSli_CPU.h"
#include "hMT_MulSli_GPU.h"

#include "cuda.h"
#include "cuda_runtime.h"
#include <device_functions.h>
#include "cufft.h"

#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	sComplex aPsih ;
	double *M2aPsih, *aM2Psih;

	cMT_InMulSli_CPU MT_InMulSli_CPU;
	f_Matlab2InMulSli(prhs[0], MT_InMulSli_CPU);

	cMT_MulSli_GPU MulSliGPU;
	MulSliGPU.SetInputData(MT_InMulSli_CPU);

	switch (MT_InMulSli_CPU.SimType){
		case 11:		// STEM
			MulSliGPU.Cal_STEM();
			f_ImSTEM2Matlab(MulSliGPU.STEM->nThk, MulSliGPU.STEM->nDet, MulSliGPU.STEM->line, MulSliGPU.STEM->nxs, MulSliGPU.STEM->nys, MulSliGPU.STEM->ImSTEM, plhs[0]);
			break;
		case 12:		// ISTEM
			MulSliGPU.Cal_STEM();
			f_ImSTEM2Matlab(MulSliGPU.STEM->nThk, MulSliGPU.STEM->nDet, MulSliGPU.STEM->line, MulSliGPU.STEM->nxs, MulSliGPU.STEM->nys, MulSliGPU.STEM->ImSTEM, plhs[0]);
			break;
		case 21:		// CBED
			plhs[0] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[1]);
			MulSliGPU.Cal_CBED(aPsih, aM2Psih);
			break;
		case 22:		// CBEI
			plhs[0] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[1]);
			MulSliGPU.Cal_CBEI(aPsih, aM2Psih);
			break;
		case 31:		// ED					
			plhs[0] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[1]);
			MulSliGPU.Cal_ED(aPsih, aM2Psih);
			break;
		case 32:		// HRTEM						
			plhs[0] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			M2aPsih = mxGetPr(plhs[1]);
			plhs[2] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[2]);
			MulSliGPU.Cal_HRTEM(aPsih, M2aPsih, aM2Psih);
			break;
		case 41:		// PED						

			break;
		case 42:		// HCI						

			break;
		case 51:		// EW Fourier				
			plhs[0] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[1]);
			MulSliGPU.Cal_ExitWaveFS(aPsih, aM2Psih);
			break;
		case 52:		// EW real						
			plhs[0] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[1]);
			MulSliGPU.Cal_ExitWaveRS(aPsih, aM2Psih);
			break;
	}
}