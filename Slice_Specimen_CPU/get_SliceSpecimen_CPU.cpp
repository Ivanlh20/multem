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
#include "hMT_MGP_CPU.h"
#include "hMT_Specimen_CPU.h"
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int nAtomsM, iConfFP, nAtoms;
	double *AtomsM;
	cMT_MGP_CPU MT_MGP_CPU;

	nAtomsM = (int)mxGetM(prhs[0]); 
	AtomsM = mxGetPr(prhs[0]);
	MT_MGP_CPU.lx = mxGetScalar(prhs[1]); 
	MT_MGP_CPU.ly = mxGetScalar(prhs[2]); 
	MT_MGP_CPU.dz = mxGetScalar(prhs[3]);
	iConfFP = (int)mxGetScalar(prhs[4]); 
	MT_MGP_CPU.DimFP = (int)mxGetScalar(prhs[5]); 
	MT_MGP_CPU.SeedFP = (int)mxGetScalar(prhs[6]);

	/*************************Input data**************************/
	cMT_Specimen_CPU MT_Specimen_CPU;
	MT_Specimen_CPU.SetInputData(&MT_MGP_CPU, nAtomsM, AtomsM);
	MT_Specimen_CPU.MoveAtoms(iConfFP);

	nAtoms = MT_Specimen_CPU.nAtoms;
	/************************Output data**************************/
	plhs[0] = mxCreateDoubleMatrix(nAtoms, 6, mxREAL);
	double *Atoms = mxGetPr(plhs[0]);

	int nSlice = MT_Specimen_CPU.nSlice;
	plhs[1] = mxCreateDoubleMatrix(nSlice, 8, mxREAL);
	double *Slice = mxGetPr(plhs[1]);

	for(int i=0; i<nAtoms; i++)
	{
		Atoms[0*nAtoms+i] = MT_Specimen_CPU.Atoms[i].x;
		Atoms[1*nAtoms+i] = MT_Specimen_CPU.Atoms[i].y;
		Atoms[2*nAtoms+i] = MT_Specimen_CPU.Atoms[i].z;
		Atoms[3*nAtoms+i] = MT_Specimen_CPU.Atoms[i].Z;
		Atoms[4*nAtoms+i] = MT_Specimen_CPU.Atoms[i].sigma;
		Atoms[5*nAtoms+i] = MT_Specimen_CPU.Atoms[i].occ;
	}

	for(int i=0; i<nSlice; i++)
	{
		Slice[0*nSlice+i] = MT_Specimen_CPU.Slice[i].z0;
		Slice[1*nSlice+i] = MT_Specimen_CPU.Slice[i].ze;
		Slice[2*nSlice+i] = MT_Specimen_CPU.Slice[i].z0_id+1;
		Slice[3*nSlice+i] = MT_Specimen_CPU.Slice[i].ze_id+1;
		Slice[4*nSlice+i] = MT_Specimen_CPU.Slice[i].z0i;
		Slice[5*nSlice+i] = MT_Specimen_CPU.Slice[i].zei;
		Slice[6*nSlice+i] = MT_Specimen_CPU.Slice[i].z0i_id+1;
		Slice[7*nSlice+i] = MT_Specimen_CPU.Slice[i].zei_id+1;
	}
}