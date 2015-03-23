/**
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

#include "math.h"
#include "hConstTypes.h"
#include "hMT_General_GPU.h"
#include "hMatlab2Cpp.h"
#include "hTEMIm.h"

#include "cuda.h"
#include "cuda_runtime.h"

#include <mex.h>

/**********************read input TEMim*************************/
void f_Matlab2InTEMIm(const mxArray *mxInTEMIm, sInTEMIm &InTEMIm)
{
	InTEMIm.GPU_Device = ReadValuemxField<int>(mxInTEMIm, 0, "GPU_Device");
	InTEMIm.MEffect = ReadValuemxField<int>(mxInTEMIm, 0, "MEffect");				// 1: Exit wave Partial coherente mode, 2: Transmission cross coefficient
	InTEMIm.STEffect = ReadValuemxField<int>(mxInTEMIm, 0, "STEffect");				// 1: Spatial and temporal, 2: Temporal, 3: Spatial
	InTEMIm.ZeroDef = ReadValuemxField<int>(mxInTEMIm, 0, "ZeroDefTyp");			// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
	InTEMIm.ZeroDefPlane = ReadValuemxField<double>(mxInTEMIm, 0, "ZeroDefPlane");	// Zero defocus plane	

	InTEMIm.E0 = ReadValuemxField<double>(mxInTEMIm, 0, "E0");			// kV
	mxArray *mxPsi = mxGetField(mxInTEMIm, 0, "Psi");					// Exit wave
	InTEMIm.Psirh = mxGetPr(mxPsi);										// real (Exit wave)
	InTEMIm.Psiih = mxGetPi(mxPsi);										// imaginary (Exit wave)
	InTEMIm.lx = ReadValuemxField<double>(mxInTEMIm, 0, "lx");			// distance in x direction(Angstroms)
	InTEMIm.ly = ReadValuemxField<double>(mxInTEMIm, 0, "ly");			// distance in y direction(Angstroms)
	InTEMIm.ny = (int)mxGetM(mxPsi);									// Number of pixels in x direction
	InTEMIm.nx = (int)mxGetN(mxPsi);									// Number of pixels in y direction

	// Microscope parameters
	InTEMIm.MC_m = ReadValuemxField<int>(mxInTEMIm, 0, "m");						// momentum of the vortex
	InTEMIm.MC_f = ReadValuemxField<double>(mxInTEMIm, 0, "f");						// defocus(Angstrom)
	InTEMIm.MC_Cs3 = ReadValuemxField<double>(mxInTEMIm, 0, "Cs3", mm2Ags);			// spherical aberration(mm-->Angstrom)
	InTEMIm.MC_Cs5 = ReadValuemxField<double>(mxInTEMIm, 0, "Cs5", mm2Ags);			// spherical aberration(mm-->Angstrom)
	InTEMIm.MC_mfa2 = ReadValuemxField<double>(mxInTEMIm, 0, "mfa2");				// magnitude 2-fold astigmatism(Angstrom)
	InTEMIm.MC_afa2 = ReadValuemxField<double>(mxInTEMIm, 0, "afa2", deg2rad);		// angle 2-fold astigmatism(degrees-->rad)
	InTEMIm.MC_mfa3 = ReadValuemxField<double>(mxInTEMIm, 0, "mfa3");				// magnitude 3-fold astigmatism(Angstrom)
	InTEMIm.MC_afa3 = ReadValuemxField<double>(mxInTEMIm, 0, "afa3", deg2rad);		// angle 3-fold astigmatism(degrees-->rad)
	InTEMIm.MC_aobjl = ReadValuemxField<double>(mxInTEMIm, 0, "aobjl", mrad2rad);	// lower objective aperture(mrad-->rad)
	InTEMIm.MC_aobju = ReadValuemxField<double>(mxInTEMIm, 0, "aobju", mrad2rad);	// upper objective aperture(mrad-->rad)
	InTEMIm.MC_sf = ReadValuemxField<double>(mxInTEMIm, 0, "sf");					// defocus spread(Angstrom)
	InTEMIm.MC_nsf = ReadValuemxField<int>(mxInTEMIm, 0, "nsf");						// Number of defocus sampling point
	InTEMIm.MC_beta = ReadValuemxField<double>(mxInTEMIm, 0, "beta", mrad2rad);		// semi-convergence angle(mrad-->rad)
	InTEMIm.MC_nbeta = ReadValuemxField<int>(mxInTEMIm, 0, "nbeta");					// half number sampling points
 }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]) 
{
	double *M2Psi_h;
	sInTEMIm InTEMIm;
	cTEMIm TEMIm;

	f_Matlab2InTEMIm(prhs[0], InTEMIm);
	TEMIm.SetInputData(InTEMIm);
	plhs[0] = mxCreateDoubleMatrix((mwSize)TEMIm.ny, (mwSize)TEMIm.nx, mxREAL);
	M2Psi_h = mxGetPr(plhs[0]);
	TEMIm.TEMImage(TEMIm.Psirh, TEMIm.Psiih, M2Psi_h);
 	/*************************************************************/
	TEMIm.freeMemory();
}