#ifndef hgeneralCPU_H
#define hgeneralCPU_H

#include "hConstTypes.h"

// Input: E0(keV), Output: lambda (electron wave)
double fgetLambda(double E0);

// Input: E0(keV), Output: sigma (Interaction parameter)
double fgetSigma(double E0);

// Input: E0(keV), Output: gamma(relativistic factor)
double fgetGamma(double E0);

// get index (with typ=0: bottom index for equal values and typ=1: upper index for equal values)
int fgetIndex(int ixmin, int ixmax, double *x, int typ, double x0);

// get two dimensional radial distribution for regular grid
void fget2DRadDist(int nR, double *R, double *fR, int nRl, double *Rl, double *rl, double *frl, double *cfrl, bool reg, int typ=0);

// Get atom types
void fAtoms2AtomTypes(int nAtoms, double *AtomsM, int PotPar, double Vrl, int &nAtomTypes, sAtomTypesCPU *&AtomTypes);

// Set atom types
void fSetAtomTypes(int Z, double occ, int PotPar, int ns, double Vrl, sAtomTypesCPU &AtomTypes);

// Set atom types
void fSetAtomTypes(int Zi, int Ze, double occ, int PotPar, int ns, double Vrl, int nAtomTypes, sAtomTypesCPU *AtomTypes);

// Grid's parameter initialization
void fsGP_Init(sGP &GP);

// Grid's parameter initialization
void fsBT_Init(sBT &BT);

// Block and Thread parameter initialization
void fsLens_Init(sLens &Lens);

// Grid's parameter calculation
void fsGP_Cal(int nx, int ny, double lx, double ly, double dz, sGP &GP);

// Block and Thread parameter calculation
void fsBT_Cal(sGP &GP, sBT &BT);

// Lens' parameter calculation
void fsLens_Cal(double E0, sGP &GP, sLens &Lens);

#endif