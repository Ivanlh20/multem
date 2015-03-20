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
#include "hMT_MulSli_CPU.h"

#include <mex.h>

/**********************read input MulSli************************/
void f_Matlab2InMulSli(const mxArray *mxInMSTEM, cMT_InMulSli_CPU &MT_InMulSli_CPU)
{
	MT_InMulSli_CPU.CPU_GPU = ReadValuemxField<int>(mxInMSTEM, 0, "CPU_GPU");	
	MT_InMulSli_CPU.nThread_CPU = ReadValuemxField<int>(mxInMSTEM, 0, "nThread_CPU");	
	MT_InMulSli_CPU.GPU_Device = ReadValuemxField<int>(mxInMSTEM, 0, "GPU_Device");				// GPU_Device device
	MT_InMulSli_CPU.SimType = ReadValuemxField<int>(mxInMSTEM, 0, "SimType");					// 11: STEM, 12: ISTEM, 21: CBED, 22: CBEI, 31: ED, 32: HRTEM, 41: PED, 42: HCI, ... 51: EW Fourier, 52: EW real
	MT_InMulSli_CPU.MulOrder = 2;																// 1: First order, 2: Second order
	MT_InMulSli_CPU.nConfFP = ReadValuemxField<int>(mxInMSTEM, 0, "nConfFP");					// Number of frozen phonon configurations
	MT_InMulSli_CPU.DimFP = ReadValuemxField<int>(mxInMSTEM, 0, "DimFP");						// Dimensions phonon configurations
	MT_InMulSli_CPU.DistFP = 1;																	// 1: Gaussian (Phonon distribution)
	MT_InMulSli_CPU.SeedFP = ReadValuemxField<int>(mxInMSTEM, 0, "SeedFP");						// Random seed(frozen phonon)
	MT_InMulSli_CPU.PotPar = ReadValuemxField<int>(mxInMSTEM, 0, "PotPar");						// Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
	MT_InMulSli_CPU.MEffect = ReadValuemxField<int>(mxInMSTEM, 0, "MEffect");					// 1: Partial coherente mode, 2: Transmission cross coefficient
	MT_InMulSli_CPU.STEffect = ReadValuemxField<int>(mxInMSTEM, 0, "STEffect");					// 1: Spatial and temporal, 2: Temporal, 3: Spatial
	MT_InMulSli_CPU.ZeroDefTyp = ReadValuemxField<int>(mxInMSTEM, 0, "ZeroDefTyp");				// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
	MT_InMulSli_CPU.ZeroDefPlane = ReadValuemxField<double>(mxInMSTEM, 0, "ZeroDefPlane");		// Zero defocus plane
	MT_InMulSli_CPU.ApproxModel = ReadValuemxField<int>(mxInMSTEM, 0, "ApproxModel");			// 1: Mulstilice, 2: Projection approximation, 3: Phase object approximation, 4: Weak phase object approximation
	MT_InMulSli_CPU.BWL = ReadValuemxField<int>(mxInMSTEM, 0, "BWL");							// 1: true, 2: false (bandwidth limited)
	MT_InMulSli_CPU.FastCal = ReadValuemxField<int>(mxInMSTEM, 0, "FastCal");					// 1: normal mode(low memory consumption), 2: fast calculation(high memory consumption)
	MT_InMulSli_CPU.PBC_xy = 1;																	// 1: true, 2: false (Peridic boundary contions)
	MT_InMulSli_CPU.Psi0Typ = ReadValuemxField<int>(mxInMSTEM, 0, "Psi0Typ");					// 1: Automatic, 2: User define
	if(MT_InMulSli_CPU.Psi0Typ!=1)
	{
		mxArray *mxPsi0 = mxGetField(mxInMSTEM, 0, "Psi0");										// Psi
		MT_InMulSli_CPU.Psi0.real = mxGetPr(mxPsi0);
		MT_InMulSli_CPU.Psi0.imag = mxGetPi(mxPsi0);
	}

	MT_InMulSli_CPU.ThkTyp = ReadValuemxField<int>(mxInMSTEM, 0, "ThkTyp");						// 1: Whole specimen, 2: Throught thickness, 3: Through planes
	mxArray *mxThk; mxThk = mxGetField(mxInMSTEM, 0, "Thk");
	MT_InMulSli_CPU.nThk = mxGetM(mxThk)*mxGetN(mxThk);											// Number of thickness
	MT_InMulSli_CPU.Thk = new double[MT_InMulSli_CPU.nThk];
	ReadValuemxField(mxInMSTEM, 0, "Thk", MT_InMulSli_CPU.nThk, MT_InMulSli_CPU.Thk);			// Array of thicknesses

	/**************************Multislice*************************/
	MT_InMulSli_CPU.E0 = ReadValuemxField<double>(mxInMSTEM, 0, "E0");							// Acceleration voltage
	MT_InMulSli_CPU.theta = ReadValuemxField<double>(mxInMSTEM, 0, "theta", deg2rad);			// incident tilt (in spherical coordinates) (degrees-->rad)
	MT_InMulSli_CPU.phi = ReadValuemxField<double>(mxInMSTEM, 0, "phi", deg2rad);				// incident tilt (in spherical coordinates) (degrees-->rad)
	MT_InMulSli_CPU.nx = ReadValuemxField<int>(mxInMSTEM, 0, "nx");								// Number of pixels in x direction
	MT_InMulSli_CPU.ny = ReadValuemxField<int>(mxInMSTEM, 0, "ny");								// Number of pixels in y direction
	MT_InMulSli_CPU.lx = ReadValuemxField<double>(mxInMSTEM, 0, "lx");							// distance in x direction(Angstroms)
	MT_InMulSli_CPU.ly = ReadValuemxField<double>(mxInMSTEM, 0, "ly");							// distance in y direction(Angstroms)
	MT_InMulSli_CPU.dz = ReadValuemxField<double>(mxInMSTEM, 0, "dz");							// Slice thickness

	/********************Microscope parameters********************/
	mxArray *mxMC= mxGetField(mxInMSTEM, 0, "MC");
	MT_InMulSli_CPU.MC_m = ReadValuemxField<int>(mxMC, 0, "m");									// momentum of the vortex
	MT_InMulSli_CPU.MC_f = ReadValuemxField<double>(mxMC, 0, "f");								// defocus(Angstrom)
	MT_InMulSli_CPU.MC_Cs3 = ReadValuemxField<double>(mxMC, 0, "Cs3", mm2Ags);					// spherical aberration(mm-->Angstrom)
	MT_InMulSli_CPU.MC_Cs5 = ReadValuemxField<double>(mxMC, 0, "Cs5", mm2Ags);					// spherical aberration(mm-->Angstrom)
	MT_InMulSli_CPU.MC_mfa2 = ReadValuemxField<double>(mxMC, 0, "mfa2");						// magnitude 2-fold astigmatism(Angstrom)
	MT_InMulSli_CPU.MC_afa2 = ReadValuemxField<double>(mxMC, 0, "afa2", deg2rad);				// angle 2-fold astigmatism(degrees-->rad)
	MT_InMulSli_CPU.MC_mfa3 = ReadValuemxField<double>(mxMC, 0, "mfa3");						// magnitude 3-fold astigmatism(Angstrom)
	MT_InMulSli_CPU.MC_afa3 = ReadValuemxField<double>(mxMC, 0, "afa3", deg2rad);				// angle 3-fold astigmatism(degrees-->rad)
	MT_InMulSli_CPU.MC_aobjl = ReadValuemxField<double>(mxMC, 0, "aobjl", mrad2rad);			// lower objective aperture(mrad-->rad)
	MT_InMulSli_CPU.MC_aobju = ReadValuemxField<double>(mxMC, 0, "aobju", mrad2rad);			// upper objective aperture(mrad-->rad)
	MT_InMulSli_CPU.MC_sf = ReadValuemxField<double>(mxMC, 0, "sf");							// defocus spread(Angstrom)
	MT_InMulSli_CPU.MC_nsf = ReadValuemxField<int>(mxMC, 0, "nsf");								// Number of defocus sampling point
	MT_InMulSli_CPU.MC_beta = ReadValuemxField<double>(mxMC, 0, "beta", mrad2rad);				// semi-convergence angle(mrad-->rad)
	MT_InMulSli_CPU.MC_nbeta = ReadValuemxField<int>(mxMC, 0, "nbeta");							// half number sampling points

	mxArray *mxAtomsM = mxGetField(mxInMSTEM, 0, "Atoms");
	MT_InMulSli_CPU.nAtomsM = (int)mxGetM(mxAtomsM);											// Number of Atoms
	MT_InMulSli_CPU.AtomsM = mxGetPr(mxAtomsM);													// Atoms in a matrix form

	switch(MT_InMulSli_CPU.SimType)
	{
		case 11:		// STEM
			mxArray *mxSTEM;	
			mxSTEM = mxGetField(mxInMSTEM, 0, "STEM");
			MT_InMulSli_CPU.STEM_line = ReadValuemxField<int>(mxSTEM, 0, "line");
			MT_InMulSli_CPU.STEM_ns = ReadValuemxField<int>(mxSTEM, 0, "ns");
			MT_InMulSli_CPU.STEM_x1u = ReadValuemxField<double>(mxSTEM, 0, "x1u");
			MT_InMulSli_CPU.STEM_y1u = ReadValuemxField<double>(mxSTEM, 0, "y1u");
			MT_InMulSli_CPU.STEM_x2u = ReadValuemxField<double>(mxSTEM, 0, "x2u");
			MT_InMulSli_CPU.STEM_y2u = ReadValuemxField<double>(mxSTEM, 0, "y2u");

			MT_InMulSli_CPU.STEM_nDet = ReadValuemxField<int>(mxSTEM, 0, "nDet");
			if(MT_InMulSli_CPU.STEM_nDet>0)
			{
				MT_InMulSli_CPU.STEM_DetCir = new sInDetCir[MT_InMulSli_CPU.STEM_nDet];
				mxArray *mxDetCir;
				mxDetCir = mxGetField(mxSTEM, 0, "DetCir");
				for(int i=0; i<MT_InMulSli_CPU.STEM_nDet; i++)
				{
					MT_InMulSli_CPU.STEM_DetCir[i].InnerAng = ReadValuemxField<double>(mxDetCir, i, "InnerAng", mrad2rad);	// Inner angle(mrad-->rad)
					MT_InMulSli_CPU.STEM_DetCir[i].OuterAng = ReadValuemxField<double>(mxDetCir, i, "OuterAng", mrad2rad);	// Outer angle(mrad-->rad)
				}
			}
			break;
		case 12:		// ISTEM
			//mxArray *mxISTEM;	
			//mxISTEM = mxGetField(mxInMSTEM, 0, "ISTEM");
			break;
		case 21:		// CBED
			mxArray *mxCBED;
			mxCBED = mxGetField(mxInMSTEM, 0, "CBED");
			MT_InMulSli_CPU.CBED_x0 = ReadValuemxField<double>(mxCBED, 0, "x0");	// 
			MT_InMulSli_CPU.CBED_y0 = ReadValuemxField<double>(mxCBED, 0, "y0");	//
			break;
		case 22:		// CBEI
			mxArray *mxCBEI;
			mxCBEI = mxGetField(mxInMSTEM, 0, "CBEI");
			MT_InMulSli_CPU.CBEI_x0 = ReadValuemxField<double>(mxCBEI, 0, "x0");	// 
			MT_InMulSli_CPU.CBEI_y0 = ReadValuemxField<double>(mxCBEI, 0, "y0");	//
			break;
		case 31:		// ED
			mxArray *mxED;
			mxED = mxGetField(mxInMSTEM, 0, "ED");

			break;
		case 32:		// HRTEM
			mxArray *mxHRTEM;
			mxHRTEM = mxGetField(mxInMSTEM, 0, "HRTEM");

			break;
		case 41:		// PED
			mxArray *mxPED;
			mxPED = mxGetField(mxInMSTEM, 0, "PED");
			MT_InMulSli_CPU.PED_nrot = ReadValuemxField<int>(mxPED, 0, "nrot");					// Number of orientations
			MT_InMulSli_CPU.PED_theta = ReadValuemxField<double>(mxPED, 0, "theta", deg2rad);	// Precession angle(degrees-->rad)
			break;
		case 42:		// HCI
			mxArray *mxHCI;
			mxHCI = mxGetField(mxInMSTEM, 0, "HCI");
			MT_InMulSli_CPU.HCI_nrot = ReadValuemxField<int>(mxHCI, 0, "nrot");					// Number of orientations
			MT_InMulSli_CPU.HCI_theta = ReadValuemxField<double>(mxHCI, 0, "theta", deg2rad);	// Precession angle(degrees-->rad)
			break;
		case 51:		// EW Fourier
			mxArray *mxEWFS;
			mxEWFS = mxGetField(mxInMSTEM, 0, "EWFS");
			break;
		case 52:		// EW real
			mxArray *mxEWRS;
			mxEWRS = mxGetField(mxInMSTEM, 0, "EWRS");
			break;
	}
 }

 // From MT_AtomTypes_CPU to Matlab structure 
void f_ImSTEM2Matlab(int nThk, int nDet, int line, int nxs, int nys, sImSTEM *ImSTEM, mxArray *&mxImSTEM)
{
	const char *field_names_ImSTEM[] = {"DetInt"};
	int number_of_fields_ImSTEM = 1;
	mwSize dims_ImSTEM[2] = {nThk, 1};

	const char *field_names_DetInt[] = {"Tot", "Coh"};
	int number_of_fields_DetInt = 2;
	mwSize dims_DetInt[2] = {nDet, 1};

	mxArray *mxDetInt;
	mxImSTEM = mxCreateStructArray(2, dims_ImSTEM, number_of_fields_ImSTEM, field_names_ImSTEM);
	for(int iThk = 0; iThk<nThk; iThk++)
	{
		mxDetInt = mxCreateStructArray(2, dims_DetInt, number_of_fields_DetInt, field_names_DetInt);
		mxSetField(mxImSTEM, iThk, "DetInt", mxDetInt);
		for(int iDet=0; iDet<nDet; iDet++)
		{
			if(line==1)
			{
				CreateSetValue2mxField(mxDetInt, iDet, "Tot", nxs, ImSTEM[iThk].DetInt[iDet].Tot);
				CreateSetValue2mxField(mxDetInt, iDet, "Coh", nxs, ImSTEM[iThk].DetInt[iDet].Coh);
			}
			else
			{
				CreateSetValue2mxField(mxDetInt, iDet, "Tot", nxs, nys, ImSTEM[iThk].DetInt[iDet].Tot);
				CreateSetValue2mxField(mxDetInt, iDet, "Coh", nxs, nys, ImSTEM[iThk].DetInt[iDet].Coh);
			}
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	sComplex aPsih ;
	double *M2aPsih, *aM2Psih;

	cMT_InMulSli_CPU MT_InMulSli_CPU;
	f_Matlab2InMulSli(prhs[0], MT_InMulSli_CPU);

	cMT_MulSli_CPU MulSliCPU;
	MulSliCPU.SetInputData(MT_InMulSli_CPU);

	switch(MT_InMulSli_CPU.SimType)
	{
		case 11:		// STEM
			MulSliCPU.Cal_STEM();
			f_ImSTEM2Matlab(MulSliCPU.STEM->nThk, MulSliCPU.STEM->nDet, MulSliCPU.STEM->line, MulSliCPU.STEM->nxs, MulSliCPU.STEM->nys, MulSliCPU.STEM->ImSTEM, plhs[0]);
			break;
		case 12:		// ISTEM
			MulSliCPU.Cal_STEM();
			f_ImSTEM2Matlab(MulSliCPU.STEM->nThk, MulSliCPU.STEM->nDet, MulSliCPU.STEM->line, MulSliCPU.STEM->nxs, MulSliCPU.STEM->nys, MulSliCPU.STEM->ImSTEM, plhs[0]);
			break;
		case 21:		// CBED
			plhs[0] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[1]);
			MulSliCPU.Cal_CBED(aPsih, aM2Psih);
			break;
		case 22:		// CBEI
			plhs[0] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[1]);
			MulSliCPU.Cal_CBEI(aPsih, aM2Psih);
			break;
		case 31:		// ED					
			plhs[0] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[1]);
			MulSliCPU.Cal_ED(aPsih, aM2Psih);
			break;
		case 32:		// HRTEM						
			plhs[0] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			M2aPsih = mxGetPr(plhs[1]);
			plhs[2] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[2]);
			MulSliCPU.Cal_HRTEM(aPsih, M2aPsih, aM2Psih);
			break;
		case 41:		// PED						
			plhs[0] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[1]);
			MulSliCPU.Cal_PED(aPsih, aM2Psih);
			break;
		case 42:		// HCI						
			plhs[0] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			M2aPsih = mxGetPr(plhs[1]);
			plhs[2] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[2]);
			MulSliCPU.CAL_HCI(aPsih, M2aPsih, aM2Psih);
			break;
		case 51:		// EW Fourier				
			plhs[0] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[1]);
			MulSliCPU.Cal_ExitWaveFS(aPsih, aM2Psih);
			break;
		case 52:		// EW real						
			plhs[0] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxCOMPLEX);
			aPsih.real = mxGetPr(plhs[0]); aPsih.imag = mxGetPi(plhs[0]);
			plhs[1] = mxCreateDoubleMatrix(MT_InMulSli_CPU.ny, MT_InMulSli_CPU.nx, mxREAL);
			aM2Psih = mxGetPr(plhs[1]);
			MulSliCPU.Cal_ExitWaveRS(aPsih, aM2Psih);
			break;
	}
}