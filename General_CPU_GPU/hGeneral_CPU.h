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

#ifndef hGeneralCPU_H
#define hGeneralCPU_H

#include "hConstTypes.h"

// Input: E0(keV), Output: lambda (electron wave)
double f_getLambda(double E0);

// Input: E0(keV), Output: sigma (Interaction parameter)
double f_getSigma(double E0);

// Input: E0(keV), Output: gamma(relativistic factor)
double f_getGamma(double E0);

// Input: E0(keV), Output: gamma*lambda/cPotf
double f_getfPot(double E0, double theta=0.0);

// get index (with typ=0: bottom index for equal values and typ=1: upper index for equal values)
int f_getIndex(int ixmin, int ixmax, double *x, int typ, double x0);

// get two dimensional Hanning_Filter
void f_getHanning_Filter_2D(int ny, int nx, double dx, double dy, double k, int shift, double *fI);

// get two dimensional Gaussian_Filter
void f_getGaussian_Filter_2D(int ny, int nx, double dx, double dy, double Sigma, int shift, double *fI);

// get two dimensional Butterworth_Filter
void f_getButterworth_Filter_2D(int ny, int nx, double dx, double dy, double Radius, int n, int lpf, int shift, double *fI);

// get two dimensional radial distribution for regular grid
void f_get2DRadDist(int nR, double *R, double *fR, int nRl, double *Rl, double *rl, double *frl, double *cfrl, bool reg, int typ=0);

// get information limit for regular grid
double f_getFFT_InformationLimit_2D(int ny, int nx, int shift, double *fI);

// Set Atoms
void f_AtomsM2Atoms(int nAtomsM_i, double *AtomsM_i, int PBC_xyi, double lxi, double lyi, int &nAtoms, sAtoms *&Atoms, double &sigma_min, double &sigma_max);

// Match vector s_i in x
void f_MatchTwoVectors(int ns_i, double *s_i, int nx, double *x, int &ns_o, double *&s_o);

// Match vector s_i in x
void f_MatchTwoVectors(int ns_i, double *s_i, int nx, double *x, int &ns_o, double *&s_o, int *&is_o);

// Build Grid
void f_BuildGrid(int line, int ns, double x0, double y0, double xe, double ye, int &nxs, int &nys, double *&xs, double *&ys);

#endif