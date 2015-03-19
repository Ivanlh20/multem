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
#include "hMT_AtomTypes_CPU.h"
#include "hMatlab2Cpp.h"
#include <mex.h>

// From MT_AtomTypes_CPU to Matlab structure 
void f_AtomTypesCPU2Matlab(int nAtomTypesCPU, cMT_AtomTypes_CPU *&MT_AtomTypes_CPU, mxArray *&mxAtomTypesCPU)
{
	const char *field_names[] = {"Z", "m", "A", "rn_e", "rn_c", "ra_e", "ra_c", "Rmin", "Rmax", "cfeg", "cfxg", "cPr", "cVr", "cVR", "R2", "ciVR"};
	int number_of_fields = 16;
	mwSize dims[2] = {nAtomTypesCPU, 1};

	mxArray *mxfield_CoefPar;
	const char *field_names_CoefPar[] = {"cl", "cnl"};
	int number_of_fields_CoefPar = 2;
	mwSize dims_CoefPar[2] = {1, 1};

	mxArray *mxfield_ciVR;
	const char *field_names_ciVR[] = {"c0", "c1", "c2", "c3"};
	int number_of_fields_ciVR = 4;
	mwSize dims_ciVR[2] = {1, 1};

	int i, j;

	mxAtomTypesCPU = mxCreateStructArray(2, dims, number_of_fields, field_names);
	for(i=0; i<nAtomTypesCPU; i++)
	{
		CreateSetValue2mxField(mxAtomTypesCPU, i, "Z", MT_AtomTypes_CPU[i].Z);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "m", MT_AtomTypes_CPU[i].m);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "A", MT_AtomTypes_CPU[i].A);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "rn_e", MT_AtomTypes_CPU[i].rn_e);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "rn_c", MT_AtomTypes_CPU[i].rn_c);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "ra_e", MT_AtomTypes_CPU[i].ra_e);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "ra_c", MT_AtomTypes_CPU[i].ra_c);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "Rmin", MT_AtomTypes_CPU[i].Rmin);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "Rmax", MT_AtomTypes_CPU[i].Rmax);

		/*************************fg***************************/
		mxfield_CoefPar = mxCreateStructArray(2, dims_CoefPar, number_of_fields_CoefPar, field_names_CoefPar);
		mxSetField(mxAtomTypesCPU, i, "cfeg", mxfield_CoefPar);
		CreateSetValue2mxField(mxfield_CoefPar, 0, "cl", 6, MT_AtomTypes_CPU[i].cfeg.cl);
		CreateSetValue2mxField(mxfield_CoefPar, 0, "cnl", 6, MT_AtomTypes_CPU[i].cfeg.cnl);

		/*************************fx***************************/
		mxfield_CoefPar = mxCreateStructArray(2, dims_CoefPar, number_of_fields_CoefPar, field_names_CoefPar);
		mxSetField(mxAtomTypesCPU, i, "cfxg", mxfield_CoefPar);
		CreateSetValue2mxField(mxfield_CoefPar, 0, "cl", 6, MT_AtomTypes_CPU[i].cfxg.cl);
		CreateSetValue2mxField(mxfield_CoefPar, 0, "cnl", 6, MT_AtomTypes_CPU[i].cfxg.cnl);

		/*************************Pr***************************/
		mxfield_CoefPar = mxCreateStructArray(2, dims_CoefPar, number_of_fields_CoefPar, field_names_CoefPar);
		mxSetField(mxAtomTypesCPU, i, "cPr", mxfield_CoefPar);
		CreateSetValue2mxField(mxfield_CoefPar, 0, "cl", 6, MT_AtomTypes_CPU[i].cPr.cl);
		CreateSetValue2mxField(mxfield_CoefPar, 0, "cnl", 6, MT_AtomTypes_CPU[i].cPr.cnl);

		/*************************Vr***************************/
		mxfield_CoefPar = mxCreateStructArray(2, dims_CoefPar, number_of_fields_CoefPar, field_names_CoefPar);
		mxSetField(mxAtomTypesCPU, i, "cVr", mxfield_CoefPar);
		CreateSetValue2mxField(mxfield_CoefPar, 0, "cl", 6, MT_AtomTypes_CPU[i].cVr.cl);
		CreateSetValue2mxField(mxfield_CoefPar, 0, "cnl", 6, MT_AtomTypes_CPU[i].cVr.cnl);

		/*************************VR***************************/
		mxfield_CoefPar = mxCreateStructArray(2, dims_CoefPar, number_of_fields_CoefPar, field_names_CoefPar);
		mxSetField(mxAtomTypesCPU, i, "cVR", mxfield_CoefPar);
		CreateSetValue2mxField(mxfield_CoefPar, 0, "cl", 6, MT_AtomTypes_CPU[i].cVR.cl);
		CreateSetValue2mxField(mxfield_CoefPar, 0, "cnl", 6, MT_AtomTypes_CPU[i].cVR.cnl);

		/*************************ciVR***************************/
		CreateSetValue2mxField(mxAtomTypesCPU, i, "R2", MT_AtomTypes_CPU[i].nR, MT_AtomTypes_CPU[i].R2);
		mxfield_ciVR = mxCreateStructArray(2, dims_ciVR, number_of_fields_ciVR, field_names_ciVR);
		mxSetField(mxAtomTypesCPU, i, "ciVR", mxfield_ciVR);
		CreateSetValue2mxField(mxfield_ciVR, 0, "c0", MT_AtomTypes_CPU[i].nR, MT_AtomTypes_CPU[i].ciVR.c0);
		CreateSetValue2mxField(mxfield_ciVR, 0, "c1", MT_AtomTypes_CPU[i].nR, MT_AtomTypes_CPU[i].ciVR.c1);
		CreateSetValue2mxField(mxfield_ciVR, 0, "c2", MT_AtomTypes_CPU[i].nR, MT_AtomTypes_CPU[i].ciVR.c3);
		CreateSetValue2mxField(mxfield_ciVR, 0, "c3", MT_AtomTypes_CPU[i].nR, MT_AtomTypes_CPU[i].ciVR.c3);
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int nMT_AtomTypes;
	cMT_AtomTypes_CPU *MT_AtomTypes_CPU;

	int PotPar = (int)mxGetScalar(prhs[0]);

	nMT_AtomTypes = stNAE;
	MT_AtomTypes_CPU = new cMT_AtomTypes_CPU[nMT_AtomTypes];

	for(int i=0; i<nMT_AtomTypes; i++)
		MT_AtomTypes_CPU[i].SetAtomTypes(i+1, PotPar, stVrl, stnR, 0);

	f_AtomTypesCPU2Matlab(nMT_AtomTypes, MT_AtomTypes_CPU, plhs[0]);
 
	delete [] MT_AtomTypes_CPU; MT_AtomTypes_CPU = 0;
}