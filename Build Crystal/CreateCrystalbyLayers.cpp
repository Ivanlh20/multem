#include "..\General\hCrystalCPU.h"
#include "..\General\hMatlab2Cpp.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int na, nb, nc;
	double a, b, c;
	int nuLayer;
	sAtomsGroup *uLayer;
	int nAtomsM;
	double *AtomsM;
	cCrystalCPU CrystalCPU;
       
	Matlab2uLayer(prhs[0], na, nb, nc, a, b, c, nuLayer, uLayer);
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