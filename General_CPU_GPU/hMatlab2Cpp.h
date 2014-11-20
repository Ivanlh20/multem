/**
 *  This file is part of MULTEM.
 *  Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
 *
 *  MULTEM is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MULTEM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MULTEM.  If not, see <http://www.gnu.org/licenses/>.
 */

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