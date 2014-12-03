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

#include <cstring>
#include "hConstTypes.h"
#include "hMT_AtomTypes_CPU.h"
#include "hMT_InMulSli_CPU.h"
#include <mex.h>

/***************************MUlTEM********************************/
// From MT_AtomTypes_CPU to Matlab structure 
void f_ImSTEM2Matlab(int nThk, int nDet, int line, int nxs, int nys, sImSTEM *ImSTEM, mxArray *&mxImSTEM){
	const char *field_names_ImSTEM[] = {"DetInt"};
	int number_of_fields_ImSTEM = 1;
	mwSize dims_ImSTEM[2] = {nThk, 1};

	const char *field_names_DetInt[] = {"Tot", "Coh"};
	int number_of_fields_DetInt = 2;
	mwSize dims_DetInt[2] = {nDet, 1};

	mxArray *mxDetInt;
	mxImSTEM = mxCreateStructArray(2, dims_ImSTEM, number_of_fields_ImSTEM, field_names_ImSTEM);
	for(int iThk = 0; iThk<nThk; iThk++){
		mxDetInt = mxCreateStructArray(2, dims_DetInt, number_of_fields_DetInt, field_names_DetInt);
		mxSetField(mxImSTEM, iThk, "DetInt", mxDetInt);
		for(int iDet=0; iDet<nDet; iDet++){
			if(line==1){
				CreateSetValue2mxField(mxDetInt, iDet, "Tot", nxs, ImSTEM[iThk].DetInt[iDet].Tot);
				CreateSetValue2mxField(mxDetInt, iDet, "Coh", nxs, ImSTEM[iThk].DetInt[iDet].Coh);
			}else{
				CreateSetValue2mxField(mxDetInt, iDet, "Tot", nxs, nys, ImSTEM[iThk].DetInt[iDet].Tot);
				CreateSetValue2mxField(mxDetInt, iDet, "Coh", nxs, nys, ImSTEM[iThk].DetInt[iDet].Coh);
			}
		}
	}
}

/**********************read input TEMim*************************/
void f_Matlab2InTEMIm(const mxArray *mxInTEMIm, sInTEMIm &InTEMIm){
	InTEMIm.gpu = ReadValuemxField<int>(mxInTEMIm, 0, "gpu");
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

/**********************read input MulSli************************/
void f_Matlab2InMulSli(const mxArray *mxInMSTEM, cMT_InMulSli_CPU &MT_InMulSli_CPU){
	MT_InMulSli_CPU.gpu = ReadValuemxField<int>(mxInMSTEM, 0, "gpu");							// gpu device
	MT_InMulSli_CPU.SimType = ReadValuemxField<int>(mxInMSTEM, 0, "SimType");					// 1: STEM, 2: HRTEM, 3: ED, 4: PED, 5: CBED, 6: HCI, ... 10: EW real, 11: EW Fourier
	MT_InMulSli_CPU.nConfFP = ReadValuemxField<int>(mxInMSTEM, 0, "nConfFP");					// Number of frozen phonon configurations
	MT_InMulSli_CPU.DimFP = ReadValuemxField<int>(mxInMSTEM, 0, "DimFP");						// Dimensions phonon configurations
	MT_InMulSli_CPU.SeedFP = ReadValuemxField<int>(mxInMSTEM, 0, "SeedFP");						// Random seed(frozen phonon)
	MT_InMulSli_CPU.PotPar = ReadValuemxField<int>(mxInMSTEM, 0, "PotPar");						// Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
	MT_InMulSli_CPU.MEffect = ReadValuemxField<int>(mxInMSTEM, 0, "MEffect");					// 1: Partial coherente mode, 2: Transmission cross coefficient
	MT_InMulSli_CPU.STEffect = ReadValuemxField<int>(mxInMSTEM, 0, "STEffect");					// 1: Spatial and temporal, 2: Temporal, 3: Spatial
	MT_InMulSli_CPU.ZeroDefTyp = ReadValuemxField<int>(mxInMSTEM, 0, "ZeroDefTyp");				// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
	MT_InMulSli_CPU.ZeroDefPlane = ReadValuemxField<double>(mxInMSTEM, 0, "ZeroDefPlane");		// Zero defocus plane
	MT_InMulSli_CPU.ApproxModel = ReadValuemxField<int>(mxInMSTEM, 0, "ApproxModel");			// 1: MS, 2: PA, 3POA, 4:WPOA
	MT_InMulSli_CPU.BandwidthLimit = ReadValuemxField<int>(mxInMSTEM, 0, "BandwidthLimit");		// 1: true, 2: false

	MT_InMulSli_CPU.ThicknessTyp = ReadValuemxField<int>(mxInMSTEM, 0, "ThicknessTyp");		// 1: Whole specimen, 2: Throught thickness, 3: Through planes
	mxArray *mxThickness;
	mxThickness = mxGetField(mxInMSTEM, 0, "Thickness");
	MT_InMulSli_CPU.nThickness = mxGetM(mxThickness)*mxGetN(mxThickness);						// Number of thickness
	MT_InMulSli_CPU.Thickness = new double[MT_InMulSli_CPU.nThickness];
	ReadValuemxField(mxInMSTEM, 0, "Thickness", MT_InMulSli_CPU.nThickness, MT_InMulSli_CPU.Thickness);	// Array of thicknesses

	/**************************Multislice*************************/
	MT_InMulSli_CPU.E0 = ReadValuemxField<double>(mxInMSTEM, 0, "E0");					// Acceleration voltage
	MT_InMulSli_CPU.theta = ReadValuemxField<double>(mxInMSTEM, 0, "theta", deg2rad);	// incident tilt (in spherical coordinates) (degrees-->rad)
	MT_InMulSli_CPU.phi = ReadValuemxField<double>(mxInMSTEM, 0, "phi", deg2rad);		// incident tilt (in spherical coordinates) (degrees-->rad)
	MT_InMulSli_CPU.nx = ReadValuemxField<int>(mxInMSTEM, 0, "nx");						// Number of pixels in x direction
	MT_InMulSli_CPU.ny = ReadValuemxField<int>(mxInMSTEM, 0, "ny");						// Number of pixels in y direction
	MT_InMulSli_CPU.lx = ReadValuemxField<double>(mxInMSTEM, 0, "lx");					// distance in x direction(Angstroms)
	MT_InMulSli_CPU.ly = ReadValuemxField<double>(mxInMSTEM, 0, "ly");					// distance in y direction(Angstroms)
	MT_InMulSli_CPU.dz = ReadValuemxField<double>(mxInMSTEM, 0, "dz");					// slice thickness

	/********************Microscope parameters********************/
	mxArray *mxMC= mxGetField(mxInMSTEM, 0, "MC");
	MT_InMulSli_CPU.MC_m =ReadValuemxField<int>(mxMC, 0, "m");						// momentum of the vortex
	MT_InMulSli_CPU.MC_f = ReadValuemxField<double>(mxMC, 0, "f");						// defocus(Angstrom)
	MT_InMulSli_CPU.MC_Cs3 = ReadValuemxField<double>(mxMC, 0, "Cs3", mm2Ags);			// spherical aberration(mm-->Angstrom)
	MT_InMulSli_CPU.MC_Cs5 = ReadValuemxField<double>(mxMC, 0, "Cs5", mm2Ags);			// spherical aberration(mm-->Angstrom)
	MT_InMulSli_CPU.MC_mfa2 = ReadValuemxField<double>(mxMC, 0, "mfa2");				// magnitude 2-fold astigmatism(Angstrom)
	MT_InMulSli_CPU.MC_afa2 = ReadValuemxField<double>(mxMC, 0, "afa2", deg2rad);		// angle 2-fold astigmatism(degrees-->rad)
	MT_InMulSli_CPU.MC_mfa3 = ReadValuemxField<double>(mxMC, 0, "mfa3");				// magnitude 3-fold astigmatism(Angstrom)
	MT_InMulSli_CPU.MC_afa3 = ReadValuemxField<double>(mxMC, 0, "afa3", deg2rad);		// angle 3-fold astigmatism(degrees-->rad)
	MT_InMulSli_CPU.MC_aobjl = ReadValuemxField<double>(mxMC, 0, "aobjl", mrad2rad);	// lower objective aperture(mrad-->rad)
	MT_InMulSli_CPU.MC_aobju = ReadValuemxField<double>(mxMC, 0, "aobju", mrad2rad);	// upper objective aperture(mrad-->rad)
	MT_InMulSli_CPU.MC_sf = ReadValuemxField<double>(mxMC, 0, "sf");					// defocus spread(Angstrom)
	MT_InMulSli_CPU.MC_nsf = ReadValuemxField<int>(mxMC, 0, "nsf");						// Number of defocus sampling point
	MT_InMulSli_CPU.MC_beta = ReadValuemxField<double>(mxMC, 0, "beta", mrad2rad);		// semi-convergence angle(mrad-->rad)
	MT_InMulSli_CPU.MC_nbeta = ReadValuemxField<int>(mxMC, 0, "nbeta");					// half number sampling points

	mxArray *mxAtomsM = mxGetField(mxInMSTEM, 0, "Atoms");
	MT_InMulSli_CPU.nAtomsM = (int)mxGetM(mxAtomsM);									// Number of Atoms
	MT_InMulSli_CPU.AtomsM = mxGetPr(mxAtomsM);											// Atoms in a matrix form

	switch (MT_InMulSli_CPU.SimType){
		case 11:		// STEM
			mxArray *mxSTEM;	
			mxSTEM = mxGetField(mxInMSTEM, 0, "STEM");
			MT_InMulSli_CPU.STEM_line = ReadValuemxField<int>(mxSTEM, 0, "line");
			MT_InMulSli_CPU.STEM_FastCal = ReadValuemxField<int>(mxSTEM, 0, "FastCal");
			MT_InMulSli_CPU.STEM_ns = ReadValuemxField<int>(mxSTEM, 0, "ns");
			MT_InMulSli_CPU.STEM_x1u = ReadValuemxField<double>(mxSTEM, 0, "x1u");
			MT_InMulSli_CPU.STEM_y1u = ReadValuemxField<double>(mxSTEM, 0, "y1u");
			MT_InMulSli_CPU.STEM_x2u = ReadValuemxField<double>(mxSTEM, 0, "x2u");
			MT_InMulSli_CPU.STEM_y2u = ReadValuemxField<double>(mxSTEM, 0, "y2u");

			MT_InMulSli_CPU.STEM_nDet = ReadValuemxField<int>(mxSTEM, 0, "nDet");
			if(MT_InMulSli_CPU.STEM_nDet>0){
				MT_InMulSli_CPU.STEM_DetCir = new sInDetCir[MT_InMulSli_CPU.STEM_nDet];
				mxArray *mxDetCir;
				mxDetCir = mxGetField(mxSTEM, 0, "DetCir");
				for (int i=0; i<MT_InMulSli_CPU.STEM_nDet; i++){
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
