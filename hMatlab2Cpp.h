#ifndef hMatlab2Cpp_H
#define hMatlab2Cpp_H

#include "hConstTypes.h"
#include "mex.h"

void Matlab2RadSchr(const mxArray *mxRadSchr, sRadSchr &RadSchr);

/********************************************************************/

void Matlab2InTEMIm(const mxArray *mxInTEMIm, sInTEMIm &InTEMIm);

void Matlab2InMulSli(const mxArray *mxInMSTEM, sInMSTEM &InMSTEM);

void CheckMulSliPar(sInMSTEM &InMSTEMi, sInMSTEM &InMSTEMo);

/********************************************************************/

void AtomTypesCPU2Matlab(int nAtomTypesCPU, sAtomTypesCPU *&AtomTypesCPU, mxArray *&mxAtomTypesCPU);

void Matlab2AtomTypesCPU(const mxArray *mxAtomTypesCPU, int &nAtomTypesCPU, sAtomTypesCPU *&AtomTypesCPU);

#endif