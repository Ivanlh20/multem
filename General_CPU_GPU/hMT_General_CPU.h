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

#ifndef hMT_GeneralCPU_H
#define hMT_GeneralCPU_H

#include "hConstTypes.h"
#include "hMT_AtomTypes_CPU.h"
#include "hMT_MGP_CPU.h"
#include "hMT_inMulSli_CPU.h"

// Input: E0(keV), Output: lambda (electron wave)
double f_getLambda(double E0);

// Input: E0(keV), Output: sigma (Interaction parameter)
double f_getSigma(double E0);

// Input: E0(keV), Output: gamma(relativistic factor)
double f_getGamma(double E0);

// get index (with typ=0: bottom index for equal values and typ=1: upper index for equal values)
int f_getIndex(int ixmin, int ixmax, double *x, int typ, double x0);

// get two dimensional radial distribution for regular grid
void f_get2DRadDist(int nR, double *R, double *fR, int nRl, double *Rl, double *rl, double *frl, double *cfrl, bool reg, int typ=0);

// Set Atoms
void f_AtomsM2Atoms(int nAtomsM_i, double *AtomsM_i, bool PBC_xyi, double lxi, double lyi, int &nAtoms, sAtoms *&Atoms, double &sigma_min, double &sigma_max);

// get 2D maximum interaction distance
double f_getRMax(int nAtoms, sAtoms *&Atoms, cMT_AtomTypes_CPU *&MT_AtomTypes_CPU);

/***************************************************************************/
/***************************************************************************/

// Grid's parameter initialization
void f_sGP_Init(sGP &GP);

// Grid's parameter calculation
void f_sGP_Cal(int nx, int ny, double lx, double ly, double dz, bool PBC_xy, bool BWL, sGP &GP);

// Grid's parameter calculation
void f_sGP_SetInputData(cMT_MGP_CPU *MT_MGP_CPU, sGP &GP);

/***************************************************************************/
/***************************************************************************/

void f_get_BTnxny(sGP &GP, dim3 &B, dim3 &T);

void f_get_BTnxhnyh(sGP &GP, dim3 &B, dim3 &T);

void f_get_BTmnxny(sGP &GP, dim3 &B, dim3 &T);

void f_get_BTnxy(sGP &GP, dim3 &B, dim3 &T);

/***************************************************************************/
/***************************************************************************/

// Block and Thread parameter initialization
void f_sLens_Init(sLens &Lens);

// Lens' parameter calculation
void f_sLens_Cal(double E0, sGP &GP, sLens &Lens);

// Set input data Lens' parameter
void f_sLens_SetInputData(cMT_InMulSli_CPU &MT_InMulSli_CPU, sGP &GP, sLens &Lens);
/***************************************************************************/
/***************************************************************************/

void f_sCoefPar_Free(sCoefPar &CoefPar);

void f_sCoefPar_Init(sCoefPar &CoefPar);

void f_sCoefPar_Malloc(int nCoefPar, sCoefPar &CoefPar);

/***************************************************************************/
/***************************************************************************/

void f_sciVn_Free(sciVn &ciVn);

void f_sciVn_Init(sciVn &ciVn);

void f_sciVn_Malloc(int nciVn, sciVn &ciVn);

/***************************************************************************/
/***************************************************************************/

void f_sDetCir_Free(sDetCir &DetCir);

void f_sDetCir_Init(sDetCir &DetCir);

void f_sDetCir_Malloc(int nDetCir, sDetCir &DetCir);

/***************************************************************************/
/***************************************************************************/

void f_scVp_Init(int ncVp, scVp *cVp);

/***************************************************************************/
/***************************************************************************/

void f_BuildGrid(int line, int ns, double x0, double y0, double xe, double ye, int &nxs, int &nys, double *&xs, double *&ys);

#endif