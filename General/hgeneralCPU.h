#ifndef hgeneralCPU_H
#define hgeneralCPU_H

#include "hConstTypes.h"

// Input: E0(keV), Output: lambda (electron wave)
double f_getLambda(double E0);

// Input: E0(keV), Output: sigma (Interaction parameter)
double f_getSigma(double E0);

// Input: E0(keV), Output: gamma(relativistic factor)
double f_get2DRadDist(double E0);

// get index (with typ=0: bottom index for equal values and typ=1: upper index for equal values)
int f_getIndex(int ixmin, int ixmax, double *x, int typ, double x0);

// get two dimensional radial distribution for regular grid
void f_get2DRadDist(int nR, double *R, double *fR, int nRl, double *Rl, double *rl, double *frl, double *cfrl, bool reg, int typ=0);

// Set atom types
void f_SetAtomTypes(int Z, int PotPar, int ns, double Vrl, sAtomTypesCPU &AtomTypes);

// Set atom types
void f_SetAtomTypes(int PotPar, int ns, double Vrl, int nAtomTypes, sAtomTypesCPU *&AtomTypes);

// Set Atoms
void f_AtomsM2Atoms(int nAtomsM_i, double *AtomsM_i, bool PBC_xyi, double lxi, double lyi, int &nAtoms, sAtoms *&Atoms);

// get 2D maximum interaction distance
double f_getRMax(int nAtoms, sAtoms *&Atoms, sAtomTypesCPU *&AtomTypes);

/*************************************************************************/
// Multislice's general parameters initialization
void f_sMPG_Init(sMGP &MGP);

/*************************************************************************/
// Grid's parameter initialization
void f_sGP_Init(sGP &GP);

// Grid's parameter calculation
void f_sGP_Cal(int nx, int ny, double lx, double ly, double dz, bool PBC_xy, sGP &GP);

/*************************************************************************/
// Grid's parameter initialization
void f_sBT_Init(sBT &BT);

// Block and Thread parameter calculation
void f_sBT_Cal(sGP &GP, sBT &BT);

/*************************************************************************/
// Block and Thread parameter initialization
void f_sLens_Init(sLens &Lens);

// Lens' parameter calculation
void f_sLens_Cal(double E0, sGP &GP, sLens &Lens);

/*************************************************************************/
void f_sCoefPar_Free(sCoefPar &CoefPar);

void f_sCoefPar_Init(sCoefPar &CoefPar);

void f_sCoefPar_Malloc(int nCoefPar, sCoefPar &CoefPar);

/*************************************************************************/
void f_sciVn_Free(sciVn &ciVn);

void f_sciVn_Init(sciVn &ciVn);

void f_sciVn_Malloc(int nciVn, sciVn &ciVn);

/*************************************************************************/
void f_sDetCir_Free(sDetCir &DetCir);

void f_sDetCir_Init(sDetCir &DetCir);

void f_sDetCir_Malloc(int nDetCir, sDetCir &DetCir);

/*************************************************************************/
void f_BuildGrid(int line, int ns, double x0, double y0, double xe, double ye, int &nxs, int &nys, double *&xs, double *&ys);

/*************************************************************************/
void f_InMulSli_Free(sInMSTEM &InMSTEM);

void f_InMulSli_Init(sInMSTEM &InMSTEM);

#endif