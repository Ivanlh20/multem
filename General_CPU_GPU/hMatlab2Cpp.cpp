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
#include "hMT_InMULTEM_CPU.h"
#include <mex.h>

void iScaleVector(int n, double *Vi, double f, double *Vo){
	for(int i=0; i<n; i++)
		Vo[i] = Vi[i]/f;
}

void CreateSetValue2mxField(mxArray *mxB, int p, const char *field_name, double field_value){
	mxArray *mxfield = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(mxfield) = field_value;
	mxSetField(mxB, p, field_name, mxfield);
}

void CreateSetValue2mxField(mxArray *mxB, int p, const char *field_name, int n, double *field_value){
	mxArray *mxfield = mxCreateDoubleMatrix(1, n, mxREAL);
	double *pfield = mxGetPr(mxfield); 
	memcpy(pfield, field_value, n*cSizeofRD);
	mxSetField(mxB, p, field_name, mxfield);
}

void CreateSetValue2mxField(mxArray *mxB, int p, const char *field_name, int nx, int ny, double *field_value){
	mxArray *mxfield = mxCreateDoubleMatrix(ny, nx, mxREAL);
	double *pfield = mxGetPr(mxfield); 
	memcpy(pfield, field_value, ny*nx*cSizeofRD);
	mxSetField(mxB, p, field_name, mxfield);
}

template <class Type>
Type ReadValuemxField(const mxArray *mxB, int p, const char *field_name, double factor=1){
	double val =factor*mxGetScalar(mxGetField(mxB, p, field_name));
	return (Type)val;
}

void ReadValuemxField(const mxArray *mxB, int p, const char *field_name, int n, double *field_value, double factor=1){
	double *pfield = mxGetPr(mxGetField(mxB, p, field_name));
	memcpy(field_value, pfield, n*cSizeofRD);
	for(int i=0; i<n; i++)
		field_value[i] *= factor;
}

/*******************Matlab to layer unit cell*********************/
void f_Matlab2uLayer(const mxArray *mxCrystal, int &na, int &nb, int &nc, double &a, double &b, double &c, int &nuLayer, sAtomsGroup *&uLayer){
int i, j, nAtoms;
	double *AtomsM;

	na = ReadValuemxField<int>(mxCrystal, 0, "na"); 
	nb = ReadValuemxField<int>(mxCrystal, 0, "nb");
	nc = ReadValuemxField<int>(mxCrystal, 0, "nc"); 

	a = ReadValuemxField<double>(mxCrystal, 0, "a");
	b = ReadValuemxField<double>(mxCrystal, 0, "b"); 
	c = ReadValuemxField<double>(mxCrystal, 0, "c");

	nuLayer = ReadValuemxField<int>(mxCrystal, 0, "nuLayer");

	mxArray *mexuLayer, *mexAtoms; 
	mexuLayer = mxGetField(mxCrystal, 0, "uLayer");
	
	uLayer = new sAtomsGroup[nuLayer];
	for (i=0; i<nuLayer; i++){
		mexAtoms = mxGetField(mexuLayer, i, "Atoms");
		AtomsM = mxGetPr(mexAtoms);
		uLayer[i].nAtoms = nAtoms = (int)mxGetM(mexAtoms);
		uLayer[i].Atoms = new sAtoms[nAtoms];
		for(j=0; j<nAtoms; j++){
			uLayer[i].Atoms[j].x = a*AtomsM[0*nAtoms+j];
			uLayer[i].Atoms[j].y = b*AtomsM[1*nAtoms+j];
			uLayer[i].Atoms[j].z = c*AtomsM[2*nAtoms+j];
			uLayer[i].Atoms[j].Z = (int)AtomsM[3*nAtoms+j];
			uLayer[i].Atoms[j].sigma = AtomsM[4*nAtoms+j];
			uLayer[i].Atoms[j].occ = AtomsM[5*nAtoms+j];
		}
	}
}

/*************Matlab to radial Schrodinger equation***************/
void f_Matlab2RadSchr(const mxArray *mxRadSchr, sInRadSchr &InRadSchr){
	InRadSchr.E0 = ReadValuemxField<double>(mxRadSchr, 0, "E0");
	InRadSchr.PotPar = ReadValuemxField<int>(mxRadSchr, 0, "PotPar");
	InRadSchr.n = ReadValuemxField<int>(mxRadSchr, 0, "n");
	InRadSchr.nr = ReadValuemxField<int>(mxRadSchr, 0, "nr");
	InRadSchr.nAtomsM = (int)mxGetM(mxGetField(mxRadSchr, 0, "Atoms"));
	InRadSchr.AtomsM = mxGetPr(mxGetField(mxRadSchr, 0, "Atoms"));
}

/***************************MUlTEM********************************/
// From AtomTypesCPU to Matlab structure 
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

// From AtomTypesCPU to Matlab structure 
void f_AtomTypesCPU2Matlab(int nAtomTypesCPU, sAtomTypesCPU *&AtomTypesCPU, mxArray *&mxAtomTypesCPU){
	const char *field_names[] = {"Z", "m", "A", "rn_e", "rn_c", "ra_e", "ra_c", "Rmin", "Rmax", "cfeg", "cfxg", "cPr", "cVr", "cVR", "ns", "Vo"};
	int number_of_fields = 16;
	mwSize dims[2] = {nAtomTypesCPU, 1};

	const char *field_names_Vpog[] = {"sigma", "Vr", "Vi", "gr", "gVr", "gVi"};
	int number_of_fields_Vpog = 6;
	mwSize dims_Vpog[2] = {1, 1};

	const char *field_names_Coef[] = {"cl", "cnl"};
	int number_of_fields_Coef = 2;
	mwSize dims_Coef[2] = {1, 1};

	mxArray *mxfield, *mxVpog;
	int i, j;

	mxAtomTypesCPU = mxCreateStructArray(2, dims, number_of_fields, field_names);
	for (i=0; i<nAtomTypesCPU; i++){
		CreateSetValue2mxField(mxAtomTypesCPU, i, "Z", AtomTypesCPU[i].Z);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "m", AtomTypesCPU[i].m);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "A", AtomTypesCPU[i].A);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "rn_e", AtomTypesCPU[i].rn_e);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "rn_c", AtomTypesCPU[i].rn_c);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "ra_e", AtomTypesCPU[i].ra_e);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "ra_c", AtomTypesCPU[i].ra_c);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "Rmin", AtomTypesCPU[i].Rmin);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "Rmax", AtomTypesCPU[i].Rmax);

		/*************************fg***************************/
		mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mxAtomTypesCPU, i, "cfeg", mxfield);
		CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cfeg.cl);
		CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cfeg.cnl);

		/*************************fx***************************/
		mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mxAtomTypesCPU, i, "cfxg", mxfield);
		CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cfxg.cl);
		CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cfxg.cnl);

		/*************************Pr***************************/
		mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mxAtomTypesCPU, i, "cPr", mxfield);
		CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cPr.cl);
		CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cPr.cnl);

		/*************************Vr***************************/
		mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mxAtomTypesCPU, i, "cVr", mxfield);
		CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cVr.cl);
		CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cVr.cnl);

		/*************************VR***************************/
		mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mxAtomTypesCPU, i, "cVR", mxfield);
		CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cVR.cl);
		CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cVR.cnl);

		CreateSetValue2mxField(mxAtomTypesCPU, i, "ns", AtomTypesCPU[i].ns);

		if (AtomTypesCPU[i].ns>0){
			dims_Vpog[1] = AtomTypesCPU[i].ns;
			mxVpog = mxCreateStructArray(2, dims_Vpog, number_of_fields_Vpog, field_names_Vpog);
			mxSetField(mxAtomTypesCPU, i, "Vo", mxVpog);
			for (j=0; j<AtomTypesCPU[i].ns; j++){
				CreateSetValue2mxField(mxVpog, j, "sigma", AtomTypesCPU[i].Vo[j].sigma);

				/*************************Vr***************************/
				mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
				mxSetField(mxVpog, j, "cVr", mxfield);
				CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Vo[j].cVr.cl);
				CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Vo[j].cVr.cnl);

				/*************************Vi***************************/
				mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
				mxSetField(mxVpog, j, "cVi", mxfield);
				CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Vo[j].cVi.cl);
				CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Vo[j].cVi.cnl);

				/************************Grid**************************/
				CreateSetValue2mxField(mxVpog, j, "gr", stngbp, AtomTypesCPU[i].Vo[j].gr);
				CreateSetValue2mxField(mxVpog, j, "gVr", stngbp, AtomTypesCPU[i].Vo[j].gVr);
				CreateSetValue2mxField(mxVpog, j, "gVi", stngbp, AtomTypesCPU[i].Vo[j].gVi);
			}
		}
	}
}

// From Matlab structure to AtomTypesCPU
void f_Matlab2AtomTypesCPU(const mxArray *mxAtomTypesCPU, int &nAtomTypesCPU, sAtomTypesCPU *&AtomTypesCPU){
	mxArray *mxfield, *mxVpog;	
	int i, j, ns;

	nAtomTypesCPU = (int)mxGetM(mxAtomTypesCPU)*(int)mxGetN(mxAtomTypesCPU);
	delete [] AtomTypesCPU; AtomTypesCPU = new sAtomTypesCPU[nAtomTypesCPU]; 
	for (i=0; i<nAtomTypesCPU; i++){
		AtomTypesCPU[i].Z = ReadValuemxField<int>(mxAtomTypesCPU, i, "Z");
		AtomTypesCPU[i].m = ReadValuemxField<double>(mxAtomTypesCPU, i, "m");
		AtomTypesCPU[i].A = ReadValuemxField<int>(mxAtomTypesCPU, i, "A");
		AtomTypesCPU[i].rn_e = ReadValuemxField<double>(mxAtomTypesCPU, i, "rn_e");
		AtomTypesCPU[i].rn_c = ReadValuemxField<double>(mxAtomTypesCPU, i, "rn_c");
		AtomTypesCPU[i].ra_e = ReadValuemxField<double>(mxAtomTypesCPU, i, "ra_e");
		AtomTypesCPU[i].ra_c = ReadValuemxField<double>(mxAtomTypesCPU, i, "ra_c");
		AtomTypesCPU[i].Rmin = ReadValuemxField<double>(mxAtomTypesCPU, i, "Rmin");
		AtomTypesCPU[i].Rmax = ReadValuemxField<double>(mxAtomTypesCPU, i, "Rmax");

		/******************************fg***********************************/
		mxfield = mxGetField(mxAtomTypesCPU, i, "cfeg");
		ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cfeg.cl);
		ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cfeg.cnl);

		/******************************fx***********************************/
		mxfield = mxGetField(mxAtomTypesCPU, i, "cfxg");
		ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cfxg.cl);
		ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cfxg.cnl);

		/******************************Pr***********************************/
		mxfield = mxGetField(mxAtomTypesCPU, i, "cPr");
		ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cPr.cl);
		ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cPr.cnl);

		/******************************Vr***********************************/
		mxfield = mxGetField(mxAtomTypesCPU, i, "cVr");
		ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cVr.cl);
		ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cVr.cnl);

		/******************************VR***********************************/
		mxfield = mxGetField(mxAtomTypesCPU, i, "cVR");
		ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cVR.cl);
		ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cVR.cnl);

		mxVpog= mxGetField(mxAtomTypesCPU, i, "Vo");
		ns = (int)mxGetM(mxVpog)*(int)mxGetN(mxVpog);

		if(ns>0){
			AtomTypesCPU[i].Vo = new sVoCPU[ns];
			for (j=0; j<ns; j++){
				AtomTypesCPU[i].Vo[j].sigma = ReadValuemxField<double>(mxVpog, j, "sigma");

				/******************************Vr**********************************/
				mxfield = mxGetField(mxVpog, j, "cVr");
				ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Vo[j].cVr.cl);
				ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Vo[j].cVr.cnl);

				/******************************Vi**********************************/
				mxfield = mxGetField(mxVpog, j, "cVi");
				ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Vo[j].cVi.cl);
				ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Vo[j].cVi.cnl);

				/*****************************Grid*********************************/
				ReadValuemxField(mxVpog, j, "gr", stngbp, AtomTypesCPU[i].Vo[j].gr);
				ReadValuemxField(mxVpog, j, "gVr", stngbp, AtomTypesCPU[i].Vo[j].gVr);
				ReadValuemxField(mxVpog, j, "gVi", stngbp, AtomTypesCPU[i].Vo[j].gVi);
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
void f_Matlab2InMulSli(const mxArray *mxInMSTEM, cMT_InMULTEM_CPU &MT_InMULTEM_CPU){
	MT_InMULTEM_CPU.gpu = ReadValuemxField<int>(mxInMSTEM, 0, "gpu");							// gpu device
	MT_InMULTEM_CPU.SimType = ReadValuemxField<int>(mxInMSTEM, 0, "SimType");					// 1: STEM, 2: HRTEM, 3: ED, 4: PED, 5: CBED, 6: HCI, ... 10: EW real, 11: EW Fourier
	MT_InMULTEM_CPU.nConfFP = ReadValuemxField<int>(mxInMSTEM, 0, "nConfFP");					// Number of frozen phonon configurations
	MT_InMULTEM_CPU.DimFP = ReadValuemxField<int>(mxInMSTEM, 0, "DimFP");						// Dimensions phonon configurations
	MT_InMULTEM_CPU.SeedFP = ReadValuemxField<int>(mxInMSTEM, 0, "SeedFP");						// Random seed(frozen phonon)
	MT_InMULTEM_CPU.PotPar = ReadValuemxField<int>(mxInMSTEM, 0, "PotPar");						// Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
	MT_InMULTEM_CPU.MEffect = ReadValuemxField<int>(mxInMSTEM, 0, "MEffect");					// 1: Partial coherente mode, 2: Transmission cross coefficient
	MT_InMULTEM_CPU.STEffect = ReadValuemxField<int>(mxInMSTEM, 0, "STEffect");					// 1: Spatial and temporal, 2: Temporal, 3: Spatial
	MT_InMULTEM_CPU.ZeroDefTyp = ReadValuemxField<int>(mxInMSTEM, 0, "ZeroDefTyp");				// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
	MT_InMULTEM_CPU.ZeroDefPlane = ReadValuemxField<double>(mxInMSTEM, 0, "ZeroDefPlane");		// Zero defocus plane
	MT_InMULTEM_CPU.ApproxModel = ReadValuemxField<int>(mxInMSTEM, 0, "ApproxModel");			// 1: MS, 2: PA, 3POA, 4:WPOA
	MT_InMULTEM_CPU.BandwidthLimit = ReadValuemxField<int>(mxInMSTEM, 0, "BandwidthLimit");		// 1: true, 2: false

	MT_InMULTEM_CPU.ThicknessTyp = ReadValuemxField<int>(mxInMSTEM, 0, "ThicknessTyp");		// 1: Whole specimen, 2: Throught thickness, 3: Through planes
	mxArray *mxThickness;
	mxThickness = mxGetField(mxInMSTEM, 0, "Thickness");
	MT_InMULTEM_CPU.nThickness = mxGetM(mxThickness)*mxGetN(mxThickness);						// Number of thickness
	MT_InMULTEM_CPU.Thickness = new double[MT_InMULTEM_CPU.nThickness];
	ReadValuemxField(mxInMSTEM, 0, "Thickness", MT_InMULTEM_CPU.nThickness, MT_InMULTEM_CPU.Thickness);	// Array of thicknesses

	/**************************Multislice*************************/
	MT_InMULTEM_CPU.E0 = ReadValuemxField<double>(mxInMSTEM, 0, "E0");					// Acceleration voltage
	MT_InMULTEM_CPU.theta = ReadValuemxField<double>(mxInMSTEM, 0, "theta", deg2rad);	// incident tilt (in spherical coordinates) (degrees-->rad)
	MT_InMULTEM_CPU.phi = ReadValuemxField<double>(mxInMSTEM, 0, "phi", deg2rad);		// incident tilt (in spherical coordinates) (degrees-->rad)
	MT_InMULTEM_CPU.nx = ReadValuemxField<int>(mxInMSTEM, 0, "nx");						// Number of pixels in x direction
	MT_InMULTEM_CPU.ny = ReadValuemxField<int>(mxInMSTEM, 0, "ny");						// Number of pixels in y direction
	MT_InMULTEM_CPU.lx = ReadValuemxField<double>(mxInMSTEM, 0, "lx");					// distance in x direction(Angstroms)
	MT_InMULTEM_CPU.ly = ReadValuemxField<double>(mxInMSTEM, 0, "ly");					// distance in y direction(Angstroms)
	MT_InMULTEM_CPU.dz = ReadValuemxField<double>(mxInMSTEM, 0, "dz");					// slice thickness

	/********************Microscope parameters********************/
	mxArray *mxMC= mxGetField(mxInMSTEM, 0, "MC");
	MT_InMULTEM_CPU.MC_m =ReadValuemxField<int>(mxMC, 0, "m");						// momentum of the vortex
	MT_InMULTEM_CPU.MC_f = ReadValuemxField<double>(mxMC, 0, "f");						// defocus(Angstrom)
	MT_InMULTEM_CPU.MC_Cs3 = ReadValuemxField<double>(mxMC, 0, "Cs3", mm2Ags);			// spherical aberration(mm-->Angstrom)
	MT_InMULTEM_CPU.MC_Cs5 = ReadValuemxField<double>(mxMC, 0, "Cs5", mm2Ags);			// spherical aberration(mm-->Angstrom)
	MT_InMULTEM_CPU.MC_mfa2 = ReadValuemxField<double>(mxMC, 0, "mfa2");				// magnitude 2-fold astigmatism(Angstrom)
	MT_InMULTEM_CPU.MC_afa2 = ReadValuemxField<double>(mxMC, 0, "afa2", deg2rad);		// angle 2-fold astigmatism(degrees-->rad)
	MT_InMULTEM_CPU.MC_mfa3 = ReadValuemxField<double>(mxMC, 0, "mfa3");				// magnitude 3-fold astigmatism(Angstrom)
	MT_InMULTEM_CPU.MC_afa3 = ReadValuemxField<double>(mxMC, 0, "afa3", deg2rad);		// angle 3-fold astigmatism(degrees-->rad)
	MT_InMULTEM_CPU.MC_aobjl = ReadValuemxField<double>(mxMC, 0, "aobjl", mrad2rad);	// lower objective aperture(mrad-->rad)
	MT_InMULTEM_CPU.MC_aobju = ReadValuemxField<double>(mxMC, 0, "aobju", mrad2rad);	// upper objective aperture(mrad-->rad)
	MT_InMULTEM_CPU.MC_sf = ReadValuemxField<double>(mxMC, 0, "sf");					// defocus spread(Angstrom)
	MT_InMULTEM_CPU.MC_nsf = ReadValuemxField<int>(mxMC, 0, "nsf");						// Number of defocus sampling point
	MT_InMULTEM_CPU.MC_beta = ReadValuemxField<double>(mxMC, 0, "beta", mrad2rad);		// semi-convergence angle(mrad-->rad)
	MT_InMULTEM_CPU.MC_nbeta = ReadValuemxField<int>(mxMC, 0, "nbeta");					// half number sampling points

	mxArray *mxAtomsM = mxGetField(mxInMSTEM, 0, "Atoms");
	MT_InMULTEM_CPU.nAtomsM = (int)mxGetM(mxAtomsM);									// Number of Atoms
	MT_InMULTEM_CPU.AtomsM = mxGetPr(mxAtomsM);											// Atoms in a matrix form

	switch (MT_InMULTEM_CPU.SimType){
		case 11:		// STEM
			mxArray *mxSTEM;	
			mxSTEM = mxGetField(mxInMSTEM, 0, "STEM");
			MT_InMULTEM_CPU.STEM_line = ReadValuemxField<int>(mxSTEM, 0, "line");
			MT_InMULTEM_CPU.STEM_FastCal = ReadValuemxField<int>(mxSTEM, 0, "FastCal");
			MT_InMULTEM_CPU.STEM_ns = ReadValuemxField<int>(mxSTEM, 0, "ns");
			MT_InMULTEM_CPU.STEM_x1u = ReadValuemxField<double>(mxSTEM, 0, "x1u");
			MT_InMULTEM_CPU.STEM_y1u = ReadValuemxField<double>(mxSTEM, 0, "y1u");
			MT_InMULTEM_CPU.STEM_x2u = ReadValuemxField<double>(mxSTEM, 0, "x2u");
			MT_InMULTEM_CPU.STEM_y2u = ReadValuemxField<double>(mxSTEM, 0, "y2u");

			MT_InMULTEM_CPU.STEM_nDet = ReadValuemxField<int>(mxSTEM, 0, "nDet");
			if(MT_InMULTEM_CPU.STEM_nDet>0){
				MT_InMULTEM_CPU.STEM_DetCir = new sInDetCir[MT_InMULTEM_CPU.STEM_nDet];
				mxArray *mxDetCir;
				mxDetCir = mxGetField(mxSTEM, 0, "DetCir");
				for (int i=0; i<MT_InMULTEM_CPU.STEM_nDet; i++){
					MT_InMULTEM_CPU.STEM_DetCir[i].InnerAng = ReadValuemxField<double>(mxDetCir, i, "InnerAng", mrad2rad);	// Inner angle(mrad-->rad)
					MT_InMULTEM_CPU.STEM_DetCir[i].OuterAng = ReadValuemxField<double>(mxDetCir, i, "OuterAng", mrad2rad);	// Outer angle(mrad-->rad)
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
			MT_InMULTEM_CPU.CBED_x0 = ReadValuemxField<double>(mxCBED, 0, "x0");	// 
			MT_InMULTEM_CPU.CBED_y0 = ReadValuemxField<double>(mxCBED, 0, "y0");	//
			break;
		case 22:		// CBEI
			mxArray *mxCBEI;
			mxCBEI = mxGetField(mxInMSTEM, 0, "CBEI");
			MT_InMULTEM_CPU.CBEI_x0 = ReadValuemxField<double>(mxCBEI, 0, "x0");	// 
			MT_InMULTEM_CPU.CBEI_y0 = ReadValuemxField<double>(mxCBEI, 0, "y0");	//
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
			MT_InMULTEM_CPU.PED_nrot = ReadValuemxField<int>(mxPED, 0, "nrot");					// Number of orientations
			MT_InMULTEM_CPU.PED_theta = ReadValuemxField<double>(mxPED, 0, "theta", deg2rad);	// Precession angle(degrees-->rad)
			break;
		case 42:		// HCI
			mxArray *mxHCI;
			mxHCI = mxGetField(mxInMSTEM, 0, "HCI");
			MT_InMULTEM_CPU.HCI_nrot = ReadValuemxField<int>(mxHCI, 0, "nrot");					// Number of orientations
			MT_InMULTEM_CPU.HCI_theta = ReadValuemxField<double>(mxHCI, 0, "theta", deg2rad);	// Precession angle(degrees-->rad)
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

/**********************read input Probe*************************/
void f_Matlab2InProbe(const mxArray *mxInProbe, sInProbe &InProbe){
	InProbe.gpu = ReadValuemxField<int>(mxInProbe, 0, "gpu");						// gpu device
	
	/**************************Multislice*************************/
	InProbe.E0 = ReadValuemxField<double>(mxInProbe, 0, "E0");						// Acceleration voltage
	InProbe.theta = ReadValuemxField<double>(mxInProbe, 0, "theta", deg2rad);		// incident tilt (in spherical coordinates) (degrees-->rad)
	InProbe.phi = ReadValuemxField<double>(mxInProbe, 0, "phi", deg2rad);			// incident tilt (in spherical coordinates) (degrees-->rad)
	InProbe.nx = ReadValuemxField<int>(mxInProbe, 0, "nx");							// Number of pixels in x direction
	InProbe.ny = ReadValuemxField<int>(mxInProbe, 0, "ny");							// Number of pixels in y direction
	InProbe.lx = ReadValuemxField<double>(mxInProbe, 0, "lx");						// distance in x direction(Angstroms)
	InProbe.ly = ReadValuemxField<double>(mxInProbe, 0, "ly");						// distance in y direction(Angstroms)

	InProbe.x0 = ReadValuemxField<double>(mxInProbe, 0, "x0");	// 
	InProbe.y0 = ReadValuemxField<double>(mxInProbe, 0, "y0");	//

	InProbe.m =ReadValuemxField<int>(mxInProbe, 0, "m");						// momentum of the vortex
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