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

#include "math.h"
#include <cstring>
#include "hConstTypes.h"
#include "hMT_General_CPU.h"
#include "hMatlab2Cpp.h"
#include "hMT_MulSli_GPU.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include <device_functions.h>
#include "cufft.h"
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	sInMSTEM InMSTEM;
	sComplex aPsih ;
	double *M2aPsih, *aM2Psih;

	f_InMulSli_Init(InMSTEM);
	Matlab2InMulSli(prhs[0], InMSTEM);

	cMT_MulSli_GPU MulSliGPU;
	MulSliGPU.SetInputData(InMSTEM);

	switch (InMSTEM.SimType){
		case 1:
			MulSliGPU.Cal_STEM();
			ImSTEM2Matlab(MulSliGPU.STEM.nThk, MulSliGPU.STEM.nDet, MulSliGPU.STEM.line, MulSliGPU.STEM.nxs, MulSliGPU.STEM.nys, MulSliGPU.STEM.ImSTEM, plhs[0]);
			break;
		case 2:		
			plhs[0] = mxCreateDoubleMatrix(InMSTEM.ny, InMSTEM.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(InMSTEM.ny, InMSTEM.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[1]);
			MulSliGPU.Cal_CBED(aPsih, aM2Psih);
			break;
		case 3:						
			plhs[0] = mxCreateDoubleMatrix(InMSTEM.ny, InMSTEM.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(InMSTEM.ny, InMSTEM.nx, mxREAL);
			M2aPsih = mxGetPr(plhs[1]);
			plhs[2] = mxCreateDoubleMatrix(InMSTEM.ny, InMSTEM.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[2]);
			MulSliGPU.Cal_HRTEM(aPsih, M2aPsih, aM2Psih);
			break;
		case 4:						
			plhs[0] = mxCreateDoubleMatrix(InMSTEM.ny, InMSTEM.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(InMSTEM.ny, InMSTEM.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[1]);
			MulSliGPU.Cal_ED(aPsih, aM2Psih);
			break;
		case 5:						

			break;
		case 6:						

			break;
		case 10:						
			plhs[0] = mxCreateDoubleMatrix(InMSTEM.ny, InMSTEM.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(InMSTEM.ny, InMSTEM.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[1]);
			MulSliGPU.Cal_ExitWaveRS(aPsih, aM2Psih);
			break;
		case 11:						
			plhs[0] = mxCreateDoubleMatrix(InMSTEM.ny, InMSTEM.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(InMSTEM.ny, InMSTEM.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[1]);
			MulSliGPU.Cal_ExitWaveFS(aPsih, aM2Psih);
			break;
	}

	f_InMulSli_Free(InMSTEM);
}