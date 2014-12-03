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

#include "hMT_Crystal_CPU.h"
#include "hMatlab2Cpp.h"
#include <mex.h>

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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
int na, nb, nc;
	double a, b, c;
	int nuLayer;
	sAtomsGroup *uLayer;
	int nAtomsM;
	double *AtomsM;
	cMT_Crystal_CPU CrystalCPU;
 
	f_Matlab2uLayer(prhs[0], na, nb, nc, a, b, c, nuLayer, uLayer);
	CrystalCPU.SetInputData(na, nb, nc, a, b, c, nuLayer, uLayer, nAtomsM);
	plhs[0] = mxCreateDoubleMatrix(nAtomsM, 6, mxREAL);
	AtomsM = mxGetPr(plhs[0]);
	CrystalCPU.Create3DCrystal(nAtomsM, AtomsM);

	for (int i=0; i<nuLayer; i++){
		uLayer[i].nAtoms = 0;
		delete [] uLayer[i].Atoms; uLayer[i].Atoms = 0;
	}

	nuLayer = 0;
	delete [] uLayer; uLayer = 0;
}