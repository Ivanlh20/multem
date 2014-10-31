#ifndef hMatlab2Cpp_H
#define hMatlab2Cpp_H

#include "hConstTypes.h"

#include <mex.h>

void Matlab2uLayer(const mxArray *mxCrystal, int &na, int &nb, int &nc, double &a, double &b, double &c, int &nuLayer, sAtomsGroup *&uLayer);

/********************************************************************/
void Matlab2RadSchr(const mxArray *mxRadSchr, sInRadSchr &RadSchr);

/********************************************************************/
void ImSTEM2Matlab(int nThk, int nDet, int line, int nxs, int nys, sImSTEM *ImSTEM, mxArray *&mxImSTEM);

void AtomTypesCPU2Matlab(int nAtomTypesCPU, sAtomTypesCPU *&AtomTypesCPU, mxArray *&mxAtomTypesCPU);

void Matlab2AtomTypesCPU(const mxArray *mxAtomTypesCPU, int &nAtomTypesCPU, sAtomTypesCPU *&AtomTypesCPU);

/********************************************************************/

void Matlab2InTEMIm(const mxArray *mxInTEMIm, sInTEMIm &InTEMIm);

void Matlab2InMulSli(const mxArray *mxInMSTEM, sInMSTEM &InMSTEM);

void Matlab2InProbe(const mxArray *mxInProbe, sInProbe &InProbe);

void CheckMulSliPar(sInMSTEM &InMSTEMi, sInMSTEM &InMSTEMo);

#endif