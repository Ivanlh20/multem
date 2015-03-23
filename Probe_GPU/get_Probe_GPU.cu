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
#include "hMatlab2Cpp.h"

#include "cuda.h"
#include "cuda_runtime.h"
#include <device_functions.h>
#include "cufft.h"
#include <mex.h>

class cProbe{
	private:
		double x0;
		double y0;
		cMT_MGP_CPU MT_MGP_CPU;
		sGP GP;
		sLens Lens;
		sComplex Psii;
		double2 *Psi;
		cufftHandle PlanPsi;
		cMT_IncidentWave_GPU MT_IncidentWave_GPU;
	public:
		void freeMemory();
		cProbe();
		~cProbe();
		void SetInputData(sInProbe &InProbe);
		void getProbe(sComplex &Psih);
};

void cProbe::freeMemory()
{
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	f_sGP_Init(GP);
	f_sLens_Init(Lens);

	cudaFreen(Psii.real);
	cudaFreen(Psii.imag);
	cudaFreen(Psi);
	cufftDestroyn(PlanPsi);
}

cProbe::~cProbe()
{
	freeMemory();
}

cProbe::cProbe()
{
	f_sGP_Init(GP);
	f_sLens_Init(Lens);
	Psii.real = 0;
	Psii.imag = 0;
	Psi = 0;
	PlanPsi = 0;
}

void cProbe::SetInputData(sInProbe &InProbe)
{
	freeMemory();

	MT_MGP_CPU.SetInputData(InProbe);
	cudaSetDevice(MT_MGP_CPU.GPU_Device);

	f_sGP_SetInputData(&MT_MGP_CPU, GP);

	Lens.m = InProbe.m;
	Lens.f = InProbe.f;
	Lens.Cs3 = InProbe.Cs3;
	Lens.Cs5 = InProbe.Cs5;
	Lens.mfa2 = InProbe.mfa2;
	Lens.afa2 = InProbe.afa2;
	Lens.mfa3 = InProbe.mfa3;
	Lens.afa3 = InProbe.afa3;
	Lens.aobjl = InProbe.aobjl;
	Lens.aobju = InProbe.aobju;
	Lens.sf = InProbe.sf;
	Lens.nsf = InProbe.nsf;
	Lens.beta = InProbe.beta;
	Lens.nbeta = InProbe.nbeta;
	f_sLens_Cal(MT_MGP_CPU.E0, GP, Lens);

	double gxu = sin(MT_MGP_CPU.theta)*cos(MT_MGP_CPU.phi)/Lens.lambda;
	double gyu = sin(MT_MGP_CPU.theta)*sin(MT_MGP_CPU.phi)/Lens.lambda;

	x0 = InProbe.x0;
	y0 = InProbe.y0;

	cudaMalloc((void**)&(Psii.real), GP.nxy*cSizeofRD);
	cudaMalloc((void**)&(Psii.imag), GP.nxy*cSizeofRD);
	cudaMalloc((void**)&Psi, GP.nxy*cSizeofCD);
	cufftPlan2d(&PlanPsi, GP.nx, GP.ny, CUFFT_Z2Z);

	MT_IncidentWave_GPU.SetInputData(&MT_MGP_CPU, Lens, PlanPsi);
}

void cProbe::getProbe(sComplex &Psih)
{
	MT_IncidentWave_GPU.Psi0(x0, y0, Psi);
	// fft2shift 
	f_fft2Shift_MC_GPU(GP, Psi);
	/*********************copy data to host************************/
	f_cuDoubleComplex_2_sComplex_GPU(GP, Psi, Psii.real, Psii.imag, Psih);
}

/**********************read input Probe*************************/
void f_Matlab2InProbe(const mxArray *mxInProbe, sInProbe &InProbe)
{
	InProbe.CPU_GPU = ReadValuemxField<int>(mxInProbe, 0, "CPU_GPU");	
	InProbe.nThread_CPU = ReadValuemxField<int>(mxInProbe, 0, "nThread_CPU");	
	InProbe.GPU_Device = ReadValuemxField<int>(mxInProbe, 0, "GPU_Device");			// GPU_Device device
	InProbe.E0 = ReadValuemxField<double>(mxInProbe, 0, "E0");						// Acceleration voltage
	InProbe.theta = ReadValuemxField<double>(mxInProbe, 0, "theta", deg2rad);		// incident tilt (in spherical coordinates) (degrees-->rad)
	InProbe.phi = ReadValuemxField<double>(mxInProbe, 0, "phi", deg2rad);			// incident tilt (in spherical coordinates) (degrees-->rad)
	InProbe.nx = ReadValuemxField<int>(mxInProbe, 0, "nx");							// Number of pixels in x direction
	InProbe.ny = ReadValuemxField<int>(mxInProbe, 0, "ny");							// Number of pixels in y direction
	InProbe.lx = ReadValuemxField<double>(mxInProbe, 0, "lx");						// distance in x direction(Angstroms)
	InProbe.ly = ReadValuemxField<double>(mxInProbe, 0, "ly");						// distance in y direction(Angstroms)

	InProbe.x0 = ReadValuemxField<double>(mxInProbe, 0, "x0");						// 
	InProbe.y0 = ReadValuemxField<double>(mxInProbe, 0, "y0");						//

	InProbe.m =ReadValuemxField<int>(mxInProbe, 0, "m");							// momentum of the vortex
	InProbe.f = ReadValuemxField<double>(mxInProbe, 0, "f");						// defocus(Angstrom)
	InProbe.Cs3 = ReadValuemxField<double>(mxInProbe, 0, "Cs3", mm2Ags);			// spherical aberration(mm-->Angstrom)
	InProbe.Cs5 = ReadValuemxField<double>(mxInProbe, 0, "Cs5", mm2Ags);			// spherical aberration(mm-->Angstrom)
	InProbe.mfa2 = ReadValuemxField<double>(mxInProbe, 0, "mfa2");					// magnitude 2-fold astigmatism(Angstrom)
	InProbe.afa2 = ReadValuemxField<double>(mxInProbe, 0, "afa2", deg2rad);			// angle 2-fold astigmatism(degrees-->rad)
	InProbe.mfa3 = ReadValuemxField<double>(mxInProbe, 0, "mfa3");					// magnitude 3-fold astigmatism(Angstrom)
	InProbe.afa3 = ReadValuemxField<double>(mxInProbe, 0, "afa3", deg2rad);			// angle 3-fold astigmatism(degrees-->rad)
	InProbe.aobjl = ReadValuemxField<double>(mxInProbe, 0, "aobjl", mrad2rad);		// lower objective aperture(mrad-->rad)
	InProbe.aobju = ReadValuemxField<double>(mxInProbe, 0, "aobju", mrad2rad);		// upper objective aperture(mrad-->rad)
	InProbe.sf = ReadValuemxField<double>(mxInProbe, 0, "sf");						// defocus spread(Angstrom)
	InProbe.nsf = ReadValuemxField<int>(mxInProbe, 0, "nsf");						// Number of defocus sampling point
	InProbe.beta = ReadValuemxField<double>(mxInProbe, 0, "beta", mrad2rad);		// semi-convergence angle(mrad-->rad)
	InProbe.nbeta = ReadValuemxField<int>(mxInProbe, 0, "nbeta");					// half number sampling points
 }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	sInProbe InProbe;
	sComplex Psih;
	cProbe Probe;

	f_Matlab2InProbe(prhs[0], InProbe);
	/************************Output data**************************/
	plhs[0] = mxCreateDoubleMatrix(InProbe.ny, InProbe.nx, mxCOMPLEX);
	Psih.real = mxGetPr(plhs[0]);
	Psih.imag = mxGetPi(plhs[0]);

	Probe.SetInputData(InProbe);
	Probe.getProbe(Psih);
}