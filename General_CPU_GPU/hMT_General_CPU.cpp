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

#include <cmath>
#include "fftw3.h"
#include "hConstTypes.h"
#include "hQuadrature.h"
#include "hMT_MGP_CPU.h"
#include "hMT_InMulSli_CPU.h"
#include "hMT_AtomTypes_CPU.h"
#include "hMT_General_CPU.h"

// Input: E0(keV), Output: lambda (electron wave)
double f_getLambda(double E0)
{
	double emass = 510.99906;		// electron rest mass in keV
	double hc = 12.3984244;			// Planck's const x speed of light
		
	double lambda = hc/sqrt(E0*(2.0*emass + E0));
	return lambda;
}

// Input: E0(keV), Output: sigma (Interaction parameter)
double f_getSigma(double E0)
{
	double emass = 510.99906;		// electron rest mass in keV
	double hc = 12.3984244;			// Planck's const x speed of light
	double x = (emass + E0)/(2.0*emass + E0);
	
	double lambda = hc/sqrt(E0*(2.0*emass + E0));
	double sigma = 2.0*cPi*x/(lambda*E0);
	return sigma;
}

// Input: E0(keV), Output: gamma(relativistic factor)
double f_getGamma(double E0)
{
	double emass = 510.99906;		// electron rest mass in keV

	double gamma = (1.0 + E0/emass);
	return gamma;
}

// Input: E0(keV), Output: gamma*lambda/cPotf
double f_getfPot(double E0, double theta)
{
	double emass = 510.99906;		// electron rest mass in keV
	double hc = 12.3984244;			// Planck's const x speed of light
		
	double lambda = hc/sqrt(E0*(2.0*emass + E0));
	double gamma = (1.0 + E0/emass);
	double fPot = gamma*lambda/(cPotf*cos(theta));
	return fPot;
}

// get index (with typ=0: bottom index for equal values and typ=1: upper index for equal values)
int f_getIndex(int ixmin, int ixmax, double *x, int typ, double x0)
{
	int ixmid;	
	switch(typ)
	{
		case 0:
			do{
				ixmid = (ixmin + ixmax)>>1;		// divide by 2
				if(x0 <= x[ixmid]) 
					ixmax = ixmid;
				else 
					ixmin = ixmid;
			}while ((ixmax-ixmin)>1);
		case 1:
			do{
				ixmid = (ixmin + ixmax)>>1;		// divide by 2
				if(x0 < x[ixmid]) 
					ixmax = ixmid;
				else 
					ixmin = ixmid;
			}while ((ixmax-ixmin)>1);
	}

	if(x0==x[ixmax])
		return ixmax;
	else
		return ixmin;
}

/***************************************************************************/
/***************************************************************************/
//symmetric coordinates(Fourier Space Coordinates)
int FSC(int i, int nh, int shift=2)
{
	int j;
	if(shift==1)
		j = (i<nh)?i:i-2*nh;
	else
		j = i-nh;
	return j;
}

// get two dimensional Hanning_Filter
void f_getHanning_Filter_2D(int ny, int nx, double dx, double dy, double k, int shift, double *fI)
{	
	int ix, iy;
	int nxh = nx/2, nyh = ny/2;
	double Rx, Ry, *fx, *fy, f;
	double cx = c2Pi/(dx*(nx-1)), cy = c2Pi/(dy*(ny-1));
	fx = new double[nx];
	fy = new double[ny];

	for(ix=0; ix<nx; ix++)
	{
		Rx = FSC(ix, nxh, shift)*dx; 
		fx[ix] = 0.5*(1.0+cos(cx*Rx));
	}
	for(iy=0; iy<ny; iy++)
	{
		Ry = FSC(iy, nyh, shift)*dy; 
		fy[iy] = 0.5*(1.0+cos(cy*Ry));
	}

	for(ix=0; ix<nx; ix++)
		for(iy=0; iy<ny; iy++)
		{
			f = fx[ix]*fy[iy];
			if(f>k)
				fI[ix*ny+iy] = 1;
			else
				fI[ix*ny+iy] = f/k;
		}		

	delete [] fx; fx = 0;
	delete [] fy; fy = 0; 
}

// get two dimensional Gaussian_Filter
void f_getGaussian_Filter_2D(int ny, int nx, double dx, double dy, double Sigma, int shift, double *fI)
{	
	int ix, iy;
	int nxh = nx/2, nyh = ny/2;
	double Rx, Ry, *fx, *fy;
	double c = 0.5/(Sigma*Sigma);
	fx = new double[nx];
	fy = new double[ny];

	for(ix=0; ix<nx; ix++)
	{
		Rx = FSC(ix, nxh, shift)*dx; 
		fx[ix] = exp(-c*Rx*Rx);
	}

	for(iy=0; iy<ny; iy++)
	{
		Ry = FSC(iy, nyh, shift)*dy; 
		fy[iy] = exp(-c*Ry*Ry);
	}

	for(ix=0; ix<nx; ix++)
		for(iy=0; iy<ny; iy++)
		{
			fI[ix*ny+iy] = fx[ix]*fy[iy];
		}

	delete [] fx; fx = 0;
	delete [] fy; fy = 0; 
}

// get two dimensional Butterworth_Filter
void f_getButterworth_Filter_2D(int ny, int nx, double dx, double dy, double Radius, int n, int lpf, int shift, double *fI)
{	
	int ix, iy;
	int nxh = nx/2, nyh = ny/2;
	double R2, Rx, Ry;
	double Radius2 = Radius*Radius;

	for(ix=0; ix<nx; ix++)
	{
		Rx = FSC(ix, nxh, shift)*dx; 
		for(iy=0; iy<ny; iy++)
		{
			Ry = FSC(iy, nyh, shift)*dy;
			R2 = Rx*Rx + Ry*Ry;
			fI[ix*ny+iy] = (lpf==1)?1.0/(1.0+pow(R2/Radius2, n)):(R2<1e-14)?0.0:1.0/(1.0+pow(Radius2/R2, n));
		}
	}
}

/***************************************************************************/
/***************************************************************************/
// get two dimensional radial distribution for regular grid
void f_get2DRadDist(int nR, double *R, double *fR, int nRl, double *Rl, double *rl, double *frl, double *cfrl, bool reg, int typ)
{	
	int i, j;
	double Rlmin = Rl[0], Rlmax = Rl[nRl-1], dRl = Rl[1]-Rl[0];
 
	for(i=0; i<nRl-1; i++)
	{
		rl[i] = 0.5*(Rl[i]+Rl[i+1]);
		frl[i] = 0.0;
		cfrl[i] = 0.0;
	}

	for(i=0; i<nR; i++)
		if((Rlmin<=R[i])&&(R[i]<Rlmax))
		{
			j = (reg)?(int)floor((R[i]-Rlmin)/dRl):f_getIndex(0, nRl-1, Rl, 0, R[i]);
			frl[j] += fR[i];
			cfrl[j] += 1.0;
		}

	if(typ==0)
		for(i=0; i<nRl-1; i++)
			if(cfrl[i]>0)
				frl[i] /= cfrl[i];
}

// get two dimensional radial distribution for regular grid
void f_getCumRadDist_2D(int ny, int nx, int shift, double *fI, int nr, double *r, double *fIr)
{	
	int i, ix, iy;
	int nxh = nx/2, nyh = ny/2;
	double *cfIr;

	cfIr = new double [nr];
	for(i=0; i<nr; i++)
	{
		fIr[i] = cfIr[i] = 0.0;
	}

	double Rx, Ry, R;
	for(ix=0; ix<nx; ix++)
	{
		Rx = FSC(ix, nxh, shift);
		for(iy=0; iy<ny; iy++)
		{
			Ry = FSC(iy, nyh, shift);
			R = sqrt(Rx*Rx + Ry*Ry);
			if((0<=R)&&(R<nr))
			{
				i = (int)floor(R);
				fIr[i] += fI[ix*ny+iy];
				cfIr[i] += 1.0;
			}
		}
	}

	r[0] = 0; 
	fIr[0] /= cfIr[0];
	for(i=1; i<nr; i++)
	{
		r[i] = i;
		if(cfIr[i]>0) 
			fIr[i] /= cfIr[i];
		fIr[i] += fIr[i-1];
	}

	delete [] cfIr; cfIr = 0;
}

// get index to maximun distance
int f_get_dmaxIndex_Point2Line(double x1, double y1, double x2, double y2, int ix1, int ix2, double *x, double *y)
{
	int imax;
	double x0, y0, dmax;
	double d, c1, c2, c3;

	c1 = y2-y1;
	c2 = -(x2-x1); 
	c3 = x2*y1-y2*x1;
	d = sqrt(c1*c1*+c2*c2);
	c1 /= d; c2 /= d; c3 /= d;

	imax = 0; dmax = 0;
	for(int i=ix1; i<=ix2; i++)
	{
		x0 = x[i]; y0 = y[i];
		d = abs(c1*x0+c2*y0+c3);
		if(d>dmax)
		{
			dmax = d;
			imax = i;
		}
	}

	return imax;
}

double f_getLengthCurve(int ix1, int ix2, double *x, double *y)
{
	double x1, y1, x2, y2;
	double d, dx, dy;
	d = 0;
	for(int i=ix1; i<ix2-1; i++)
	{
		x1 = x[i]; y1 = y[i];
		x2 = x[i+1]; y2 = y[i+1];
		dx = x2-x1; dy = y2-y1;
		d = d + sqrt(dx*dx + dy*dy);
	}
	return d;
}

double f_getLengthCurve(int ix1, int ix2, double *x, double *y, double lmax, int &il)
{
	double x1, y1, x2, y2;
	double l, dx, dy;
	if(ix1<ix2)
	{
		l = 0; il = ix2;
		for(int i=ix1; i<ix2-1; i++)
		{
			x1 = x[i]; y1 = y[i];
			x2 = x[i+1]; y2 = y[i+1];
			dx = x2-x1; dy = y2-y1;
			l = l + sqrt(dx*dx + dy*dy);
			if((lmax>0)&&(l>=lmax))
			{
				il = i;
				break;
			}
		}
	}
	else
	{
		l = 0; il = ix2;
		for(int i=ix1; i>ix2; i--)
		{
			x1 = x[i-1]; y1 = y[i-1];
			x2 = x[i]; y2 = y[i];
			dx = x2-x1; dy = y2-y1;
			l = l + sqrt(dx*dx + dy*dy);
			if((lmax>0)&&(l>=lmax))
			{
				il = i;
				break;
			}
		}
	}
	return l;
}

// get information limit for regular grid
double f_getFFT_InformationLimit_2D(int ny, int nx, int shift, double *fI)
{
	int nr = MIN(nx/2,ny/2)-1;
	double *r, *fIr;

	r = new double[nr];
	fIr = new double[nr]; 

	// Cumulative radial integration
	f_getCumRadDist_2D(ny, nx, shift, fI, nr, r, fIr);

	// Shift and Normalize
	double r0 = r[0], fIr0 = fIr[0];
	for(int i=0; i<nr; i++)
	{
		r[i] = (r[i]-r0)/(r[nr-1]-r0);
		fIr[i] = (fIr[i]-fIr0)/(fIr[nr-1]-fIr0);
	}

	int ir1, ir2, irm;
	double x1, y1, x2, y2;

	ir1 = 0; ir2 = nr-1;
	x1 = r[ir1]; y1 = fIr[ir1];
	x2 = r[ir2]; y2 = fIr[ir2];
	
	irm = f_get_dmaxIndex_Point2Line(x1, y1, x2, y2, ir1, ir2, r, fIr);
	double fIr_lim = 0.45*fIr[irm];

	for(int i=0; i<nr; i++)
		if(fIr[i]>fIr_lim)
		{
			irm = i-1;
			break;
		}

	delete [] r;
	delete [] fIr;

	return irm;
}

// Set Atoms
void f_AtomsM2Atoms(int nAtomsM_i, double *AtomsM_i, int PBC_xyi, double lxi, double lyi, int &nAtoms, sAtoms *&Atoms, double &sigma_min, double &sigma_max)
{
	double x, y, dl = 1e-03;
	double lxb = lxi-dl, lyb = lyi-dl;
	double sigma;
	sigma_min = sigma_max = AtomsM_i[4*nAtomsM_i + 0];
	Atoms = new sAtoms[nAtomsM_i];
	nAtoms = 0;
	for(int i=0; i<nAtomsM_i; i++)
	{		
		x = AtomsM_i[0*nAtomsM_i + i];									// x-position
		y = AtomsM_i[1*nAtomsM_i + i];									// y-position
		if((PBC_xyi==2)||((x<lxb)&&(y<lyb)))
		{
			Atoms[nAtoms].x = x;										// x-position
			Atoms[nAtoms].y = y;										// y-position
			Atoms[nAtoms].z = AtomsM_i[2*nAtomsM_i + i];				// z-position
			Atoms[nAtoms].Z = (int)AtomsM_i[3*nAtomsM_i + i];			// Atomic number
			Atoms[nAtoms].sigma = sigma = AtomsM_i[4*nAtomsM_i + i];	// Standard deviation
			Atoms[nAtoms].occ = AtomsM_i[5*nAtomsM_i + i];				// Occupancy
			if(sigma<sigma_min) sigma_min = sigma;
			if(sigma>sigma_max) sigma_max = sigma;
			nAtoms++;
		}
	}
}

// get 2D maximum interaction distance
double f_getRMax(int nAtoms, sAtoms *&Atoms, cMT_AtomTypes_CPU *&MT_AtomTypes_CPU)
{
	double R, Rmax=0;

	for(int iAtoms=0; iAtoms<nAtoms; iAtoms++)
	{
		R = MT_AtomTypes_CPU[Atoms[iAtoms].Z-1].Rmax;
		if(Rmax < R)
			Rmax = R;
	}
	return Rmax;
}

// Match vector s_i in x
void f_MatchTwoVectors(int ns_i, double *s_i, int nx, double *x, int &ns_o, double *&s_o)
{
	int *xc;
	xc = new int[nx];
	for(int ix=0; ix<nx; ix++) xc[ix] = 0;

	int i, ix, imin;
	double dmin, d;

	delete [] s_o; s_o = new double[ns_i];
	ns_o = 0;
	for(i=0; i<ns_i; i++)
	{
		dmin = 1e+10;
		imin = -1;
		for(ix=0; ix<nx; ix++)
		{
			d = abs(s_i[i]-x[ix]);
			if(d<dmin)
			{
				dmin = d;
				imin = ix;
			}
		}
		if((imin>=0)&&(xc[imin]==0))
		{
			xc[imin] = 1;
			s_o[ns_o] = x[imin];
			ns_o++;
		}
	}
	delete [] xc; xc = 0;
}

// Match vector s_i in x
void f_MatchTwoVectors(int ns_i, double *s_i, int nx, double *x, int &ns_o, double *&s_o, int *&is_o)
{
	int *xc;
	xc = new int[nx];
	for(int ix=0; ix<nx; ix++) xc[ix] = 0;

	int i, ix, imin;
	double dmin, d;

	delete [] s_o; s_o = new double[ns_i];
	delete [] is_o; is_o = new int[ns_i];
	ns_o = 0;
	for(i=0; i<ns_i; i++)
	{
		dmin = 1e+10;
		imin = -1;
		for(ix=0; ix<nx; ix++)
		{
			d = abs(s_i[i]-x[ix]);
			if(d<dmin)
			{
				dmin = d;
				imin = ix;
			}
		}
		if((imin>=0)&&(xc[imin]==0))
		{
			xc[imin] = 1;
			s_o[ns_o] = x[imin];
			is_o[ns_o] = imin;
			ns_o++;
		}
	}
	delete [] xc; xc = 0;
}

/***************************************************************************/
/***************************************************************************/

// Grid's parameter initialization
void f_sGP_Init(sGP &GP)
{
	GP.nx = 0;
	GP.ny = 0;
	GP.nxh = 0;
	GP.nyh = 0;
	GP.nxy = 0;
	GP.inxy = 0;
	GP.lx = 0;
	GP.ly = 0;
	GP.dz = 0;
	GP.PBC_xy = 1;
	GP.BWL = 1;
	GP.dRx = 0;
	GP.dRy = 0;
	GP.dRmin = 0;
	GP.dgx = 0;
	GP.dgy = 0;
	GP.dgmin = 0;
	GP.gmax = 0;
	GP.gmax2 = 0;
	GP.gmaxl = 0;
	GP.gmaxl2 = 0;
}

// Set input data Grid's parameter
void f_sGP_Cal(int nx, int ny, double lx, double ly, double dz, int BWL, int PBC_xy, sGP &GP)
{
	GP.nx = nx;
	GP.ny = ny;
	GP.nxh = nx/2;
	GP.nyh = ny/2;
	GP.nxy = GP.nx*GP.ny; 
	GP.inxy = (GP.nxy==0)?0.0:(1.0/double(GP.nxy));
	GP.lx = lx;
	GP.ly = ly;
	GP.dz = dz;
	GP.BWL = BWL;
	GP.PBC_xy = PBC_xy;
	GP.dRx = (GP.nx==0)?0.0:(GP.lx/double(GP.nx));
	GP.dRy = (GP.ny==0)?0.0:(GP.ly/double(GP.ny));
	GP.dRmin = MIN(GP.dRx, GP.dRy);
	GP.dgx = (GP.lx==0)?0.0:(1.0/GP.lx);
	GP.dgy = (GP.ly==0)?0.0:(1.0/GP.ly);
	GP.dgmin = MIN(GP.dgx, GP.dgy);
	GP.gmax = MIN(double(GP.nxh)*GP.dgx, double(GP.nyh)*GP.dgy);
	GP.gmax2 = pow(GP.gmax, 2);
	GP.gmaxl = 2.0*GP.gmax/3.0;
	GP.gmaxl2 = pow(GP.gmaxl, 2);
}

// Grid's parameter calculation
void f_sGP_SetInputData(cMT_MGP_CPU *MT_MGP_CPU, sGP &GP)
{
	f_sGP_Cal(MT_MGP_CPU->nx, MT_MGP_CPU->ny, MT_MGP_CPU->lx, MT_MGP_CPU->ly, MT_MGP_CPU->dz, MT_MGP_CPU->BWL, MT_MGP_CPU->PBC_xy, GP);
}

/***************************************************************************/
/***************************************************************************/

// Block and Thread parameter initialization
void f_sLens_Init(sLens &Lens)
{
	Lens.gamma = 0;
	Lens.lambda = 0;
	Lens.m = 0;
	Lens.f = 0;
	Lens.Cs3 = 0;
	Lens.Cs5 = 0;
	Lens.mfa2 = 0;
	Lens.afa2 = 0;
	Lens.mfa3 = 0;
	Lens.afa3 = 0;
	Lens.aobjl = 0;
	Lens.aobju = 0;
	Lens.sf = 0;
	Lens.nsf = 0;
	Lens.beta = 0;
	Lens.nbeta = 0;
	Lens.lambda2 = 0;
	Lens.cf = 0;
	Lens.cCs3 = 0;
	Lens.cCs5 = 0;
	Lens.cmfa2 = 0;
	Lens.cmfa3 = 0;
	Lens.gmin2 = 0;
	Lens.gmax2 = 0;
	Lens.sggs = 0;
	Lens.ngxs = 0;
	Lens.ngys = 0;
	Lens.dgxs = 0;
	Lens.dgys = 0;
	Lens.gmax2s = 0;
}

// Lens' parameter calculation
void f_sLens_Cal(double E0, sGP &GP, sLens &Lens)
{
	Lens.gamma = f_getGamma(E0);

	Lens.lambda = f_getLambda(E0);
	Lens.lambda2 = pow(Lens.lambda, 2);

	Lens.cf = (Lens.f==0)?0:cPi*Lens.f*Lens.lambda;
	Lens.cCs3 = (Lens.Cs3==0)?0:-cPi*Lens.Cs3*pow(Lens.lambda, 3)/2.0;
	Lens.cCs5 = (Lens.Cs5==0)?0:-cPi*Lens.Cs5*pow(Lens.lambda, 5)/3.0;
	Lens.cmfa2 = (Lens.mfa2==0)?0:-cPi*Lens.mfa2*Lens.lambda;
	Lens.cmfa3 = (Lens.mfa3==0)?0:-2.0*cPi*Lens.mfa3*pow(Lens.lambda, 2)/3.0;
	Lens.gmin2 = (Lens.aobjl==0)?0:pow(Lens.aobjl/Lens.lambda, 2);
	Lens.gmax2 = (Lens.aobju==0)?0:pow(Lens.aobju/Lens.lambda, 2);

	double g0s = Lens.beta/Lens.lambda;
	Lens.sggs = g0s/c2i2;
	double gmaxs = 3.5*Lens.sggs;
	Lens.gmax2s = gmaxs*gmaxs;
	double dgs = gmaxs/Lens.nbeta;
	int nbeta;

	nbeta = (dgs<GP.dgx)?((int)floor(GP.dgx/dgs)+1):(1);
	Lens.ngxs = (int)floor(nbeta*gmaxs/GP.dgx) + 1;
	Lens.dgxs = gmaxs/Lens.ngxs;

	nbeta = (dgs<GP.dgy)?((int)floor(GP.dgy/dgs)+1):(1);
	Lens.ngys = (int)floor(nbeta*gmaxs/GP.dgy) + 1;
	Lens.dgys = gmaxs/Lens.ngys;
}

// Set input data Lens' parameter
void f_sLens_SetInputData(cMT_InMulSli_CPU &MT_InMulSli_CPU, sGP &GP, sLens &Lens)
{
	Lens.m = MT_InMulSli_CPU.MC_m;
	Lens.f = MT_InMulSli_CPU.MC_f;
	Lens.Cs3 = MT_InMulSli_CPU.MC_Cs3;
	Lens.Cs5 = MT_InMulSli_CPU.MC_Cs5;
	Lens.mfa2 = MT_InMulSli_CPU.MC_mfa2;
	Lens.afa2 = MT_InMulSli_CPU.MC_afa2;
	Lens.mfa3 = MT_InMulSli_CPU.MC_mfa3;
	Lens.afa3 = MT_InMulSli_CPU.MC_afa3;
	Lens.aobjl = MT_InMulSli_CPU.MC_aobjl;
	Lens.aobju = MT_InMulSli_CPU.MC_aobju;
	Lens.sf = MT_InMulSli_CPU.MC_sf;
	Lens.nsf = MT_InMulSli_CPU.MC_nsf;
	Lens.beta = MT_InMulSli_CPU.MC_beta;
	Lens.nbeta = MT_InMulSli_CPU.MC_nbeta;
	f_sLens_Cal(MT_InMulSli_CPU.E0, GP, Lens);
}

/***************************************************************************/
/***************************************************************************/

void f_ReadQuadratureCPU(int typ, int nQCPU, sQ1 &QCPU)
{
	if(nQCPU<=0)
	{
		return;
	}

	QCPU.x = new double[nQCPU]; 
	QCPU.w = new double[nQCPU];
	cQuadrature Quad;
	Quad.ReadQuadrature(typ, nQCPU, QCPU);
}

/***************************************************************************/
/***************************************************************************/

void f_sCoefPar_Free(sCoefPar &CoefPar)
{
	delete [] CoefPar.cl; CoefPar.cl = 0;
	delete [] CoefPar.cnl; CoefPar.cnl = 0;
}

void f_sCoefPar_Init(sCoefPar &CoefPar)
{
	CoefPar.cl = 0;
	CoefPar.cnl = 0;
}

void f_sCoefPar_Malloc(int nCoefPar, sCoefPar &CoefPar)
{
	if(nCoefPar<=0)
	{
		return;
	}

	CoefPar.cl = new double[nCoefPar];
	CoefPar.cnl = new double[nCoefPar];
}

/***************************************************************************/
/***************************************************************************/

void f_sciVn_Free(sciVn &ciVn)
{
	delete [] ciVn.c0; ciVn.c0 = 0;
	delete [] ciVn.c1; ciVn.c1 = 0;
	delete [] ciVn.c2; ciVn.c2 = 0;
	delete [] ciVn.c3; ciVn.c3 = 0;
}

void f_sciVn_Init(sciVn &ciVn)
{
	ciVn.c0 = 0;
	ciVn.c1 = 0;
	ciVn.c2 = 0;
	ciVn.c3 = 0;
}

void f_sciVn_Malloc(int nciVn, sciVn &ciVn)
{
	if(nciVn<=0)
	{
		return;
	}

	ciVn.c0 = new double[nciVn];
	ciVn.c1 = new double[nciVn];
	ciVn.c2 = new double[nciVn];
	ciVn.c3 = new double[nciVn];
}

/***************************************************************************/
/***************************************************************************/

void f_sDetCir_Free(sDetCir &DetCir)
{
	delete [] DetCir.g2min; DetCir.g2min = 0;
	delete [] DetCir.g2max; DetCir.g2max = 0;
}

void f_sDetCir_Init(sDetCir &DetCir)
{
	DetCir.g2min = 0;
	DetCir.g2max = 0;
}

void f_sDetCir_Malloc(int nDetCir, sDetCir &DetCir)
{
	if(nDetCir<=0)
	{
		return;
	}

	DetCir.g2min = new double[nDetCir];
	DetCir.g2max = new double[nDetCir];
}

/***************************************************************************/
/***************************************************************************/

void f_scVp_Init(scVp &cVp)
{
	cVp.x = 0;
	cVp.y = 0;
	cVp.z0 = 0;
	cVp.ze = 0;
	cVp.split = false;
	cVp.occ = 0;
	cVp.Rmin2 = 0;
	cVp.Rmax2 = 0;
	cVp.R2 = 0;
	f_sCoefPar_Init(cVp.cVr);
	cVp.bnx.i = 0;
	cVp.bnx.n = 0;
	cVp.bny.i = 0;
	cVp.bny.n = 0;
	f_sciVn_Init(cVp.ciV0);
	f_sciVn_Init(cVp.ciV1);
}

void f_scVp_Init(int ncVp, scVp *cVp)
{
	for(int icVp=0; icVp<ncVp; icVp++)
	{
		f_scVp_Init(cVp[icVp]);
	}
}

/***************************************************************************/
/***************************************************************************/

void f_BuildGrid(int line, int ns, double x0, double y0, double xe, double ye, int &nxs, int &nys, double *&xs, double *&ys)
{
	if(ns<=0)
	{
		nxs = nys = 0;
		xs = 0; ys = 0;
		return;
	}

	int i;
	double lxs = xe-x0, lys = ye-y0, ld = std::sqrt(lxs*lxs+lys*lys);
	double theta = std::atan(lys/lxs), costheta = std::cos(theta), sintheta = std::sin(theta);

	if(line)
	{
			nxs = ns;
			nys = ns;
	}
	else
	{
		nxs = (std::abs(lxs)>std::abs(lys))?ns:(int)std::ceil(ns*std::abs(lxs/lys));
		nys = (std::abs(lxs)>std::abs(lys))?(int)std::ceil(ns*std::abs(lys/lxs)):ns;
	}

	xs = new double [nxs];
	ys = new double [nys];

	double ds;

	if(line)
	{
		ds = ld/ns;
		for(i=0; i<ns; i++)
		{
			xs[i] = x0 + i*ds*costheta;
			ys[i] = y0 + i*ds*sintheta;
		}
	}
	else
	{
		ds = lxs/nxs;
		for(i=0; i<nxs; i++)
		{
			xs[i] = x0 + i*ds;
		}

		ds = lys/nys;
		for(i=0; i<nys; i++)
		{
			ys[i] = y0 + i*ds;
		}
	}

}

/***************************************************************************/
/***************************************************************************/

// Hadamard product
void f_Hadamard_Product_CPU(sGP &GP, double * __restrict A_i, double * __restrict B_io)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			B_io[ixy] = A_i[ixy]*B_io[ixy];
		}
}

// Hadamard product
void f_Hadamard_Product_CPU(sGP &GP, double * __restrict A_i, double * __restrict B_i, double * __restrict C_o)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			C_o[ixy] = A_i[ixy]*B_i[ixy];
		}
}

// Hadamard product
void f_Hadamard_Product_CPU(sGP &GP, fftw_complex * __restrict A_i, fftw_complex * __restrict B_io)
{
	int ix, iy, ixy;
	double z1r, z1i, z2r, z2i, z3r, z3i;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			z1r = A_i[ixy][0];	z1i = A_i[ixy][1];
			z2r = B_io[ixy][0];	z2i = B_io[ixy][1];
			z3r = z1r*z2r-z1i*z2i;	z3i = z1i*z2r+z1r*z2i;
			B_io[ixy][0] = z3r;
			B_io[ixy][1] = z3i;
		}
}

// Hadamard product
void f_Hadamard_Product_CPU(sGP &GP, fftw_complex * __restrict A_i, fftw_complex * __restrict B_i, fftw_complex * __restrict C_o)
{
	int ix, iy, ixy;
	double z1r, z1i, z2r, z2i, z3r, z3i;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			z1r = A_i[ixy][0];	z1i = A_i[ixy][1];
			z2r = B_i[ixy][0];	z2i = B_i[ixy][1];
			z3r = z1r*z2r-z1i*z2i;	z3i = z1i*z2r+z1r*z2i;
			C_o[ixy][0] = z3r;
			C_o[ixy][1] = z3i;
		}
}

/***************************************************************************/
/***************************************************************************/

// From double to float
void f_Double_2_Float_CPU(sGP &GP, double * __restrict V_i, float * __restrict V_o)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			V_o[ixy] = static_cast<float>(V_i[ixy]);
		}
}

// From double to float
void f_Double_2_Float_CPU(sGP &GP, double * __restrict V1_i, double * __restrict V2_i, float * __restrict V1_o, float * __restrict V2_o)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			V1_o[ixy] = static_cast<float>(V1_i[ixy]);
			V2_o[ixy] = static_cast<float>(V2_i[ixy]);
		}
}

// From float to double
void f_Float_2_Double_CPU(sGP &GP, float * __restrict V_i, double * __restrict V_o)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			V_o[ixy] = static_cast<double>(V_i[ixy]);
		}
}

// From float to double
void f_Float_2_Double_CPU(sGP &GP, float * __restrict V1_i, float * __restrict V2_i, double * __restrict V1_o, double * __restrict V2_o)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			V1_o[ixy] = static_cast<double>(V1_i[ixy]);
			V2_o[ixy] = static_cast<double>(V2_i[ixy]);
		}
}

/***************************************************************************/
/***************************************************************************/

void f_Scale_MD_CPU(sGP &GP, double w, double * __restrict MD_io)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			MD_io[ixy] *= w;
		}
}

void f_Scale_MD_CPU(sGP &GP, double w, double * __restrict MD1_io, double * __restrict MD2_io)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			MD1_io[ixy] *= w;
			MD2_io[ixy] *= w;
		}
}

void f_Scale_MC_CPU(sGP &GP, double w, fftw_complex * __restrict MC_io)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			MC_io[ixy][0] *= w;
			MC_io[ixy][1] *= w;
		}
}

void f_Scale_MC_CPU(sGP &GP, double w, fftw_complex * __restrict MC1_io, fftw_complex * __restrict MC2_io)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			MC1_io[ixy][0] *= w;
			MC1_io[ixy][1] *= w;
			MC2_io[ixy][0] *= w;
			MC2_io[ixy][1] *= w;
		}
}

/***************************************************************************/
/***************************************************************************/

// Set value to Double vector:
void f_Set_MD_CPU(sGP &GP, double M, double * __restrict MD_o)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			MD_o[ixy] = M;
		}
}

// Set value to 2 Double vector:
void f_Set_MD_CPU(sGP &GP, double M, double * __restrict MD1_o, double * __restrict MD2_o)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			MD1_o[ixy] = M;
			MD2_o[ixy] = M;
		}
}

// Set value to Complex vector
void f_Set_MC_CPU(sGP &GP, double Mr, double Mi, fftw_complex * __restrict MC_o)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			MC_o[ixy][0] = Mr;
			MC_o[ixy][1] = Mi;
		}
}

// Set value to 2 Double complex vector:
void f_Set_MC_CPU(sGP &GP, double Mr, double Mi, fftw_complex * __restrict MC1_o, fftw_complex * __restrict MC2_o)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			MC1_o[ixy][0] = Mr;
			MC1_o[ixy][1] = Mi;
			MC2_o[ixy][0] = Mr;
			MC2_o[ixy][1] = Mi;
		}
}

// Set value to Complex vector and Double vector
void f_Set_MC_MD_CPU(sGP &GP, double Mr, double Mi, fftw_complex * __restrict MC_o, double M, double * __restrict MD_o)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			MC_o[ixy][0] = Mr;
			MC_o[ixy][1] = Mi;
			MD_o[ixy] = M;
		}
}

// Set Real and Imaginary part of a Complex vector
void f_Set_MC_CPU(sGP &GP, sComplex &MC_i, fftw_complex * __restrict MC_o)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			MC_o[ixy][0] = MC_i.real[ixy];
			MC_o[ixy][1] = MC_i.imag[ixy];
		}
}

// Get Real and Imaginary part of a Complex vector
void f_Get_MC_CPU(sGP &GP, fftw_complex * __restrict MC_i, sComplex &MC_o)
{
	int ix, iy, ixy;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			MC_o.real[ixy] = MC_i[ixy][0];
			MC_o.imag[ixy] = MC_i[ixy][1];
		}
}

/***************************************************************************/
/***************************************************************************/

template <bool add>
void t_Add_Set_wMD(sGP &GP, double w, const double * __restrict MD_i, double * __restrict MD_io)
{
	int ix, iy, ixy;
	double x;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			x = MD_i[ixy];
			if(add)
				MD_io[ixy] += x = w*x;
			else
				MD_io[ixy] = x = w*x;
		}
}

template <bool add>
void t_Add_Set_wMD(sGP &GP, double w, const double * __restrict MD1_i, const double * __restrict MD2_i, double * __restrict MD1_io, double * __restrict MD2_io)
{
	int ix, iy, ixy;
	double x1, x2;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			x1 = MD1_i[ixy]; 
			x2 = MD2_i[ixy];
			if(add)
			{
				MD1_io[ixy] += x1 = w*x1;
				MD2_io[ixy] += x2 = w*x2;
			}
			else
			{
				MD1_io[ixy] = x1 = w*x1;
				MD2_io[ixy] = x2 = w*x2;
			}
		}
}

template <bool add>
void t_Add_Set_wMC(sGP &GP, double w, const fftw_complex * __restrict MC_i, fftw_complex * __restrict MC_io)
{
	int ix, iy, ixy;
	double x, y;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			x = MC_i[ixy][0]; 
			y = MC_i[ixy][1];
			if(add)
			{
				MC_io[ixy][0] += x = w*x;
				MC_io[ixy][1] += y = w*y;
			}
			else
			{
				MC_io[ixy][0] = x = w*x;
				MC_io[ixy][1] = y = w*y;
			}
		}
}

template <bool add>
void t_Add_Set_wMC(sGP &GP, double w, const fftw_complex * __restrict MC1_i, const fftw_complex * __restrict MC2_i, fftw_complex * __restrict MC1_io, fftw_complex * __restrict MC2_io)
{
	int ix, iy, ixy;
	double x1, y1, x2, y2;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			x1 = MC1_i[ixy][0]; 
			y1 = MC1_i[ixy][1];
			x2 = MC2_i[ixy][0]; 
			y2 = MC2_i[ixy][1];
			if(add)
			{
				MC1_io[ixy][0] += x1 = w*x1;
				MC1_io[ixy][1] += y1 = w*y1;
				MC2_io[ixy][0] += x2 = w*x2;
				MC2_io[ixy][1] += y2 = w*y2;
			}
			else
			{
				MC1_io[ixy][0] = x1 = w*x1;
				MC1_io[ixy][1] = y1 = w*y1;
				MC2_io[ixy][0] = x2 = w*x2;
				MC2_io[ixy][1] = y2 = w*y2;
			}
		}
}

template <bool add>
void t_Add_Set_wMC2(sGP &GP, double w, const fftw_complex * __restrict MC_i, double * __restrict MD_io)
{
	int ix, iy, ixy;
	double x, y, z;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			x = MC_i[ixy][0]; 
			y = MC_i[ixy][1]; 
			z = x*x+y*y;
			if(add)
			{
				MD_io[ixy] += z = w*z;
			}
			else
			{
				MD_io[ixy] = z = w*z;
			}
		}
}

template <bool add>
void t_Add_Set_wMC_wMD(sGP &GP, double w, const fftw_complex * __restrict MC_i, fftw_complex * __restrict MC_io, double * __restrict MD_io)
{
	int ix, iy, ixy;
	double x, y, z;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			x = MC_i[ixy][0]; 
			y = MC_i[ixy][1]; 
			z = x*x+y*y;
			if(add)
			{
				MC_io[ixy][0] += x = w*x; 
				MC_io[ixy][1] += y = w*y;
				MD_io[ixy] += z = w*z;
			}
			else
			{
				MC_io[ixy][0] = x = w*x; 
				MC_io[ixy][1] = y = w*y;
				MD_io[ixy] = z = w*z;
			}
		}
}

template <bool add>
void t_Add_Set_wMC_wMD(sGP &GP, double w, const fftw_complex * __restrict MC_i, const double * __restrict MD_i, fftw_complex * __restrict MC_io, double * __restrict MD_io)
{
	int ix, iy, ixy;
	double x, y, z;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			x = MC_i[ixy][0]; 
			y = MC_i[ixy][1]; 
			z = MD_i[ixy];
			if(add)
			{
				MC_io[ixy][0] += x = w*x; 
				MC_io[ixy][1] += y = w*y;
				MD_io[ixy] += z = w*z;
			}
			else
			{
				MC_io[ixy][0] = x = w*x; 
				MC_io[ixy][1] = y = w*y;
				MD_io[ixy] = z = w*z;
			}
		}
}


void f_Set_wMD_CPU(sGP &GP, double w, double * __restrict MD_i, double * __restrict MD_io)
{
	t_Add_Set_wMD<false>(GP, w, MD_i, MD_io);
}

void f_Set_wMD_CPU(sGP &GP, double w, double * __restrict MD1_i, double * __restrict MD2_i, double * __restrict MD1_io, double * __restrict MD2_io)
{
	t_Add_Set_wMD<false>(GP, w, MD1_i, MD2_i, MD1_io, MD2_io);
}

void f_Set_wMC_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, fftw_complex * __restrict MC_io)
{
	t_Add_Set_wMC<false>(GP, w, MC_i, MC_io);
}

void f_Set_wMC_CPU(sGP &GP, double w, fftw_complex * __restrict MC1_i, fftw_complex * __restrict MC2_i, fftw_complex * __restrict MC1_io, fftw_complex * __restrict MC2_io)
{
	t_Add_Set_wMC<false>(GP, w, MC1_i, MC2_i, MC1_io, MC2_io);
}

void f_Set_wMC2_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, double * __restrict MD_io)
{
	t_Add_Set_wMC2<false>(GP, w, MC_i, MD_io);
}

void f_Set_wMC_wMD_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, fftw_complex * __restrict MC_io, double * __restrict MD_io)
{
	t_Add_Set_wMC_wMD<false>(GP, w, MC_i, MC_io, MD_io);
}

void f_Set_wMC_wMD_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, double * __restrict MD_i, fftw_complex * __restrict MC_io, double * __restrict MD_io)
{
	t_Add_Set_wMC_wMD<false>(GP, w, MC_i, MD_i, MC_io, MD_io);
}


void f_Add_wMD_CPU(sGP &GP, double w, double * __restrict MD_i, double * __restrict MD_io)
{
	t_Add_Set_wMD<true>(GP, w, MD_i, MD_io);
}

void f_Add_wMD_CPU(sGP &GP, double w, double * __restrict MD1_i, double * __restrict MD2_i, double * __restrict MD1_io, double * __restrict MD2_io)
{
	t_Add_Set_wMD<true>(GP, w, MD1_i, MD2_i, MD1_io, MD2_io);
}

void f_Add_wMC_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, fftw_complex * __restrict MC_io)
{
	t_Add_Set_wMC<true>(GP, w, MC_i, MC_io);
}

void f_Add_wMC_CPU(sGP &GP, double w, fftw_complex * __restrict MC1_i, fftw_complex * __restrict MC2_i, fftw_complex * __restrict MC1_io, fftw_complex * __restrict MC2_io)
{
	t_Add_Set_wMC<true>(GP, w, MC1_i, MC2_i, MC1_io, MC2_io);
}

void f_Add_wMC2_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, double * __restrict MD_io)
{
	t_Add_Set_wMC2<true>(GP, w, MC_i, MD_io);
}

void f_Add_wMC_wMD_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, fftw_complex * __restrict MC_io, double * __restrict MD_io)
{
	t_Add_Set_wMC_wMD<true>(GP, w, MC_i, MC_io, MD_io);
}

void f_Add_wMC_wMD_CPU(sGP &GP, double w, fftw_complex * __restrict MC_i, double * __restrict MD_i, fftw_complex * __restrict MC_io, double * __restrict MD_io)
{
	t_Add_Set_wMC_wMD<true>(GP, w, MC_i, MD_i, MC_io, MD_io);
}

/***************************************************************************/
/***************************************************************************/
inline void ExcVal(double &v1, double &v2)
{
	 double v = v1;
	 v1 = v2;
	 v2 = v; 
}

inline void ExcVal(fftw_complex &v1, fftw_complex &v2)
{
	 double vx = v1[0], vy = v1[1];
	 v1[0] = v2[0]; v1[1] = v2[1];
	 v2[0] = vx; v2[1] = vy;
}

// Shift Double matrix respect to (nxh, nyh)
void f_fft2Shift_MD_CPU(sGP &GP, double * __restrict MD_io)
{
	int ix, iy, ixy, ixys;
	for(iy=0; iy<GP.nyh; iy++)
		for(ix=0; ix<GP.nxh; ix++)
		{		
			ixy = iy*GP.nx + ix; ixys = (GP.nyh+iy)*GP.nx+(GP.nxh+ix);
			ExcVal(MD_io[ixy], MD_io[ixys]);

			ixy = (GP.nxh+ix) + iy*GP.nx; ixys = (GP.nyh+iy)*GP.nx+ix;
			ExcVal(MD_io[ixy], MD_io[ixys]);
		}
}

// Shift Double matrix respect to (nxh, nyh)
void f_fft2Shift_MD_CPU(sGP &GP, double * __restrict MD1_io, double * __restrict MD2_io)
{
	int ix, iy, ixy, ixys;
	for(iy=0; iy<GP.nyh; iy++)
		for(ix=0; ix<GP.nxh; ix++)
		{		
			ixy = iy*GP.nx + ix; ixys = (GP.nyh+iy)*GP.nx+(GP.nxh+ix);
			ExcVal(MD1_io[ixy], MD1_io[ixys]);
			ExcVal(MD2_io[ixy], MD2_io[ixys]);

			ixy = (GP.nxh+ix) + iy*GP.nx; ixys = (GP.nyh+iy)*GP.nx+ix;
			ExcVal(MD1_io[ixy], MD1_io[ixys]);
			ExcVal(MD2_io[ixy], MD2_io[ixys]);
		}
}

// Shift Complex matrix respect to (nxh, nyh)
void f_fft2Shift_MC_CPU(sGP &GP, fftw_complex * __restrict MC_io)
{
	int ix, iy, ixy, ixys;
	for(iy=0; iy<GP.nyh; iy++)
		for(ix=0; ix<GP.nxh; ix++)
		{		
			ixy = iy*GP.nx + ix; ixys = (GP.nyh+iy)*GP.nx+(GP.nxh+ix);
			ExcVal(MC_io[ixy], MC_io[ixys]);

			ixy = (GP.nxh+ix) + iy*GP.nx; ixys = (GP.nyh+iy)*GP.nx+ix;
			ExcVal(MC_io[ixy], MC_io[ixys]);
		}
}

// Shift Complex matrix respect to (nxh, nyh)
void f_fft2Shift_MC_CPU(sGP &GP, fftw_complex * __restrict MC1_io, fftw_complex * __restrict MC2_io)
{
	int ix, iy, ixy, ixys;
	for(iy=0; iy<GP.nyh; iy++)
		for(ix=0; ix<GP.nxh; ix++)
		{		
			ixy = iy*GP.nx + ix; ixys = (GP.nyh+iy)*GP.nx+(GP.nxh+ix);
			ExcVal(MC1_io[ixy], MC1_io[ixys]);
			ExcVal(MC2_io[ixy], MC2_io[ixys]);

			ixy = (GP.nxh+ix) + iy*GP.nx; ixys = (GP.nyh+iy)*GP.nx+ix;
			ExcVal(MC1_io[ixy], MC1_io[ixys]);
			ExcVal(MC2_io[ixy], MC2_io[ixys]);
		}
}

// Shift Complex matrix respect to (nxh, nyh)
void f_fft2Shift_MC_MD_CPU(sGP& GP, fftw_complex * __restrict MC_io, double * __restrict MD_io)
{
	int ix, iy, ixy, ixys;
	for(iy=0; iy<GP.nyh; iy++)
		for(ix=0; ix<GP.nxh; ix++)
		{		
			ixy = iy*GP.nx + ix; ixys = (GP.nyh+iy)*GP.nx+(GP.nxh+ix);
			ExcVal(MC_io[ixy], MC_io[ixys]);
			ExcVal(MD_io[ixy], MD_io[ixys]);

			ixy = (GP.nxh+ix) + iy*GP.nx; ixys = (GP.nyh+iy)*GP.nx+ix;
			ExcVal(MC_io[ixy], MC_io[ixys]);
			ExcVal(MD_io[ixy], MD_io[ixys]);
		}
}

// Shift Complex matrix respect to (nxh, nyh)
void f_fft2Shift_MC_MD_CPU(sGP &GP, fftw_complex * __restrict MC_io, double * __restrict MD1_io, double * __restrict MD2_io)
{
	int ix, iy, ixy, ixys;
	for(iy=0; iy<GP.nyh; iy++)
		for(ix=0; ix<GP.nxh; ix++)
		{		
			ixy = iy*GP.nx + ix; ixys = (GP.nyh+iy)*GP.nx+(GP.nxh+ix);
			ExcVal(MC_io[ixy], MC_io[ixys]);
			ExcVal(MD1_io[ixy], MD1_io[ixys]);
			ExcVal(MD2_io[ixy], MD2_io[ixys]);

			ixy = (GP.nxh+ix) + iy*GP.nx; ixys = (GP.nyh+iy)*GP.nx+ix;
			ExcVal(MC_io[ixy], MC_io[ixys]);
			ExcVal(MD1_io[ixy], MD1_io[ixys]);
			ExcVal(MD2_io[ixy], MD2_io[ixys]);
		}
}

// Shift Complex matrix respect to (nxh, nyh)
void f_fft2Shift_MC_CPU(sGP &GP, sComplex &MC_io)
{
	int ix, iy, ixy, ixys;
	double zr, zi;
	for(iy=0; iy<GP.nyh; iy++)
		for(ix=0; ix<GP.nxh; ix++)
		{		
			ixy = iy*GP.nx + ix; ixys = (GP.nyh+iy)*GP.nx+(GP.nxh+ix);
			zr = MC_io.real[ixy]; zi = MC_io.imag[ixy];
			MC_io.real[ixy] = MC_io.real[ixys]; MC_io.imag[ixy] = MC_io.imag[ixys];
			MC_io.real[ixys] = zr; MC_io.imag[ixys] = zi;

			ixy = (GP.nxh+ix) + iy*GP.nx; ixys = (GP.nyh+iy)*GP.nx+ix;
			zr = MC_io.real[ixy]; zi = MC_io.imag[ixy];
			MC_io.real[ixy] = MC_io.real[ixys]; MC_io.imag[ixy] = MC_io.imag[ixys];
			MC_io.real[ixys] = zr; MC_io.imag[ixys] = zi;
		}
}

/***************************************************************************/
/***************************************************************************/

double f_Sum_MD_CPU(sGP &GP, double w_i, const double * __restrict MD_i)
{
	int ix, iy, ixy;
	double SumT = 0.0;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			SumT += MD_i[ixy];
		}
	return w_i*SumT;
}

double f_Sum_MC2_CPU(sGP &GP, double w_i, const fftw_complex * __restrict MC_i)
{
	int ix, iy, ixy;
	double x, y, SumT = 0.0;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			x = MC_i[ixy][0]; 
			y = MC_i[ixy][1];
			SumT += x*x+y*y;
		}
	return w_i*SumT;
}

void f_Sum_MD_Det_CPU(sGP &GP, double w_i, const double * __restrict MD_i, double gmin2, double gmax2, int iSum, double * __restrict Sum_MD_o)
{
	int ix, iy, ixy;
	double gx, gy, g2;
	double SumT = 0.0;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			gx = IsFS(ix,GP.nxh)*GP.dgx;
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			ixy = iy*GP.nx+ix;
			if((gmin2<=g2)&&(g2<=gmax2))
				SumT += MD_i[ixy];
		}

	Sum_MD_o[iSum] = w_i*SumT;
}

void f_Sum_MD_Det_CPU(sGP &GP, double w_i, const double * __restrict MD1_i, const double * __restrict MD2_i, double gmin2, double gmax2, int iSum, double * __restrict Sum_MD1_o, double * __restrict Sum_MD2_o)
{
	int ix, iy, ixy;
	double gx, gy, g2;
	double SumT1 = 0.0, SumT2 = 0.0;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			gx = IsFS(ix,GP.nxh)*GP.dgx;
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			ixy = iy*GP.nx+ix;
			if((gmin2<=g2)&&(g2<=gmax2))
			{
				SumT1 += MD1_i[ixy];
				SumT2 += MD2_i[ixy];
			}
		}

	Sum_MD1_o[iSum] = w_i*SumT1;
	Sum_MD2_o[iSum] = w_i*SumT2;
}

void f_Sum_MC_Det_CPU(sGP &GP, double w_i, const fftw_complex * __restrict MC_i, double gmin2, double gmax2, double * __restrict MDp_i, int iSum, double * __restrict Sum_MD_o)
{
	int ix, iy, ixy;
	double x, y;
	double gx, gy, g2;
	double SumT = 0.0;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			gx = IsFS(ix,GP.nxh)*GP.dgx;
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			ixy = iy*GP.nx+ix;
			if((gmin2<=g2)&&(g2<=gmax2))
			{
				x = MC_i[ixy][0]; 
				y = MC_i[ixy][1];
				SumT += x*x+y*y;
			}
		}

	Sum_MD_o[iSum] = w_i*SumT;
}

/***************************************************************************/
/***************************************************************************/

// Phase
void f_Phase_CPU(sGP &GP, double gxu, double gyu, sACD &ExpR_x_o, sACD &ExpR_y_o)
{
	int ix, iy;
	double theta, R;

	for(ix=0; ix<GP.nx; ix++)
	{
		R = IsRS(ix,GP.nxh)*GP.dRx;
		theta = c2Pi*R*gxu;
		ExpR_x_o.x[ix] = cos(theta);
		ExpR_x_o.y[ix] = sin(theta);
	}

	for(iy=0; iy<GP.ny; iy++)
	{
		R = IsRS(iy,GP.nyh)*GP.dRy;
		theta = c2Pi*R*gyu;
		ExpR_y_o.x[ix] = cos(theta);
		ExpR_y_o.y[ix] = sin(theta);
	}
}

// Phase multiplication
void f_PhaseMul_CPU(sGP &GP, double gxu, double gyu, sACD &ExpR_x_o, sACD &ExpR_y_o, fftw_complex *Psi_io)
{
	f_Phase_CPU(GP, gxu, gyu, ExpR_x_o, ExpR_y_o);

	int ix, iy, ixy;
	double z1r, z1i, z2r, z2i, z3r, z3i;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			z1r = ExpR_x_o.x[ix];	z1i = ExpR_x_o.y[ix];
			z2r = ExpR_y_o.x[iy];	z2i = ExpR_y_o.y[iy];
			z3r = z1r*z2r-z1i*z2i;	z3i = z1i*z2r+z1r*z2i;
			z2r = Psi_io[ixy][0];	z2i = Psi_io[ixy][1]; 
			Psi_io[ixy][0] = z1r = z3r*z2r-z3i*z2i; 
			Psi_io[ixy][1] = z1i = z3i*z2r+z3r*z2i; 
		}
}

// BandWidth Limit function
void f_BandwidthLimit2D_CPU(fftw_plan &PlanForward, fftw_plan &PlanBackward, sGP &GP, fftw_complex *MC_io)
{
	if(GP.BWL!=1) return;

	// Forward fft2
	fftw_execute_dft(PlanForward, MC_io, MC_io);
	// AntiAliasing, scale and bandlimited
	int ix, iy, ixy;
	double gx, gy, g2;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			gx = IsFS(ix,GP.nxh)*GP.dgx;
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			g2 = gx*gx + gy*gy;
			ixy = iy*GP.nx+ix;
			if(g2 < GP.gmaxl2)
			{
				MC_io[ixy][0] *= GP.inxy; 
				MC_io[ixy][1] *= GP.inxy; 
			}
			else
			{
 				MC_io[ixy][0] = 0.0; 
				MC_io[ixy][1] = 0.0; 
			}	
		}
	// Backward fft2
	fftw_execute_dft(PlanBackward, MC_io, MC_io);
}

// Build propagator function
void f_Propagator_CPU(sGP &GP, double gxu, double gyu, double scale, sACD &Prop_x_o, sACD &Prop_y_o)
{
	int ix, iy;
	double theta, gx, gy;

	for(ix=0; ix<GP.nx; ix++)
	{
		gx = IsFS(ix,GP.nxh)*GP.dgx;
		theta = scale*(gx+gxu)*(gx+gxu);
		Prop_x_o.x[ix] = cos(theta);
		Prop_x_o.y[ix] = sin(theta);
	}

	for(iy=0; iy<GP.ny; iy++)
	{
		gy = IsFS(iy,GP.nyh)*GP.dgy;
		theta = scale*(gy+gyu)*(gy+gyu);
		Prop_y_o.x[iy] = cos(theta);
		Prop_y_o.y[iy] = sin(theta);
	}
}

// Propagate, scale with cut-off (2/3)gmax
void f_Propagate_CPU(fftw_plan &PlanForward, fftw_plan &PlanBackward, sGP &GP, eSpace Space, double gxu, double gyu, double lambda, double z, sACD &Prop_x_o, sACD &Prop_y_o, fftw_complex *Psi_io)
{
	if(z==0)
	{
		if(Space==eSReal) 
		{
			return;
		}
		else
		{
			fftw_execute_dft(PlanForward, Psi_io, Psi_io);	// Forward fft2
			f_Scale_MC_CPU(GP, GP.inxy, Psi_io);
			return;
		}
	}

	fftw_execute_dft(PlanForward, Psi_io, Psi_io);	// Forward fft2

	f_Propagator_CPU(GP, gxu, gyu, -cPi*lambda*z, Prop_x_o, Prop_y_o);

	int ix, iy, ixy;
	double gx, gy, g2;
	double z1r, z1i, z2r, z2i, z3r, z3i;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			gx = IsFS(ix,GP.nxh)*GP.dgx;
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			g2 = gx*gx + gy*gy;
			ixy = iy*GP.nx+ix;
			if((GP.BWL!=1)||(g2 < GP.gmaxl2))
			{
				z1r = Prop_x_o.x[ix]; z1i = Prop_x_o.y[ix];
				z2r = Prop_y_o.x[iy]; z2i = Prop_y_o.y[iy];
				z3r = z1r*z2r-z1i*z2i; z3i = z1i*z2r+z1r*z2i;
				z2r = Psi_io[ixy][0];		z2i = Psi_io[ixy][1]; 
				Psi_io[ixy][0] = z1r = (z3r*z2r-z3i*z2i)*GP.inxy;
				Psi_io[ixy][1] = z1i = (z3i*z2r+z3r*z2i)*GP.inxy; 
			}
			else
			{
 				Psi_io[ixy][0] = 0.0; 
				Psi_io[ixy][1] = 0.0; 
			}
		}

	// Backward fft2
	if(Space==eSReal)
	{
		fftw_execute_dft(PlanBackward, Psi_io, Psi_io);	// Backward fft2
	}
}

// Probe in Fourier space
void f_Probe_FS_CPU(sGP &GP, sLens Lens, double x, double y, fftw_complex *fPsi_o)
{
	int ix, iy, ixy;
	double chi, phi;
	double gx, gy, g, g2;

	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			gx = IsFS(ix,GP.nxh)*GP.dgx;
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			g2 = gx*gx + gy*gy;
			ixy = iy*GP.nx+ix;
			if((Lens.gmin2 <= g2)&&(g2 < Lens.gmax2))
			{
				chi = x*gx + y*gy + g2*(Lens.cCs5*g2*g2+Lens.cCs3*g2+Lens.cf);
				if((Lens.m!=0)||(Lens.cmfa2!=0)||(Lens.cmfa3!=0))
				{
					g = sqrt(g2);
					phi = atan2(gy, gx);
					chi += Lens.m*phi + Lens.cmfa2*g2*sin(2*(phi-Lens.afa2)) + Lens.cmfa3*g*g2*sin(3*(phi-Lens.afa3));				
				}	
				fPsi_o[ixy][0] = cos(chi); 
				fPsi_o[ixy][1] = sin(chi);	
			}
			else
			{
 				fPsi_o[ixy][0] = 0.0; 
				fPsi_o[ixy][1] = 0.0; 
			}
		}
}

// Apply Coherent transfer function
void f_Apply_CTF_CPU(sGP &GP, sLens &Lens, double gxu, double gyu, fftw_complex *fPsi_io)
{
	int ix, iy, ixy;
	double chix, chiy;
	double chi, phi, Psix, Psiy;
	double gx, gy, g, g2;

	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			gx = IsFS(ix,GP.nxh)*GP.dgx;
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			g2 = gx*gx + gy*gy;
			ixy = iy*GP.nx+ix;
			if((Lens.gmin2 <= g2)&&(g2 < Lens.gmax2))
			{
				g2 = (gx-gxu)*(gx-gxu) + (gy-gyu)*(gy-gyu);
				chi = g2*(Lens.cCs5*g2*g2+Lens.cCs3*g2+Lens.cf);
				if((Lens.cmfa2!=0)||(Lens.cmfa3!=0))
				{
					g = sqrt(g2);
					phi = atan2(gy, gx);
					chi += Lens.cmfa2*g2*sin(2*(phi-Lens.afa2)) + Lens.cmfa3*g*g2*sin(3*(phi-Lens.afa3));				
				}
				Psix = fPsi_io[ixy][0]; Psiy = fPsi_io[ixy][1]; 
				chix = cos(chi); chiy = sin(chi);	
				fPsi_io[ixy][0] = gx = chix*Psix-chiy*Psiy; 
				fPsi_io[ixy][1] = gy = chix*Psiy+chiy*Psix;
			}
			else
			{
 				fPsi_io[ixy][0] = 0.0; 
				fPsi_io[ixy][1] = 0.0; 
			}
		}
}

// Apply Coherent transfer function
void f_Apply_CTF_CPU(sGP &GP, sLens &Lens, double gxu, double gyu, fftw_complex *fPsi_i, fftw_complex *fPsi_o)
{
	int ix, iy, ixy;
	double chix, chiy;
	double chi, phi, Psix, Psiy;
	double gx, gy, g, g2;

	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			gx = IsFS(ix,GP.nxh)*GP.dgx;
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			g2 = gx*gx + gy*gy;
			ixy = iy*GP.nx+ix;
			if((Lens.gmin2 <= g2)&&(g2 < Lens.gmax2))
			{
				g2 = (gx-gxu)*(gx-gxu) + (gy-gyu)*(gy-gyu);
				chi = g2*(Lens.cCs5*g2*g2+Lens.cCs3*g2+Lens.cf);
				if((Lens.cmfa2!=0)||(Lens.cmfa3!=0))
				{
					g = sqrt(g2);
					phi = atan2(gy, gx);
					chi += Lens.cmfa2*g2*sin(2*(phi-Lens.afa2)) + Lens.cmfa3*g*g2*sin(3*(phi-Lens.afa3));				
				}
				Psix = fPsi_i[ixy][0]; Psiy = fPsi_i[ixy][1]; 
				chix = cos(chi); chiy = sin(chi);	
				fPsi_o[ixy][0] = gx = chix*Psix-chiy*Psiy; 
				fPsi_o[ixy][1] = gy = chix*Psiy+chiy*Psix;
			}
			else
			{
 				fPsi_o[ixy][0] = 0.0; 
				fPsi_o[ixy][1] = 0.0; 
			}
		}
}

// Partially coherent transfer function, linear image model and weak phase object
void f_Apply_PCTF_CPU(sGP &GP, sLens &Lens, fftw_complex *fPsi_io)
{
	int ix, iy, ixy;
	double chix, chiy, c, u;
	double sie, tie, sti;
	double chi, phi, Psix, Psiy;
	double gx, gy, g, g2;

	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			gx = IsFS(ix,GP.nxh)*GP.dgx;
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			g2 = gx*gx + gy*gy;
			ixy = iy*GP.nx+ix;
			if((Lens.gmin2 <= g2)&&(g2 < Lens.gmax2))
			{			
				chi = g2*(Lens.cCs3*g2+Lens.cf);
				c = cPi*Lens.beta*Lens.sf;
				u = 1.0 + c*c*g2;

				c = cPi*Lens.sf*Lens.lambda*g2;
				sie = 0.25*c*c;

				c = cPi*Lens.beta*(Lens.Cs3*Lens.lambda2*g2-Lens.f);
				tie = c*c*g2;

				sti = exp(-(sie+tie)/u);

				Psix = fPsi_io[ixy][0]; Psiy = fPsi_io[ixy][1]; 
				chix = cos(chi); chiy = sin(chi);		 	

				fPsi_io[ixy][0] = gx = sti*(chix*Psix-chiy*Psiy); 
				fPsi_io[ixy][1] = gy = sti*(chix*Psiy+chiy*Psix);
			}
			else
			{
 				fPsi_io[ixy][0] = 0.0; 
				fPsi_io[ixy][1] = 0.0; 
			}
		}
}

// Partially coherent transfer function, linear image model and weak phase object
void f_Apply_PCTF_CPU(sGP &GP, sLens &Lens, fftw_complex *fPsi_i, fftw_complex *fPsi_o)
{
	int ix, iy, ixy;
	double chix, chiy, c, u;
	double sie, tie, sti;
	double chi, phi, Psix, Psiy;
	double gx, gy, g, g2;

	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			gx = IsFS(ix,GP.nxh)*GP.dgx;
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			g2 = gx*gx + gy*gy;
			ixy = iy*GP.nx+ix;
			if((Lens.gmin2 <= g2)&&(g2 < Lens.gmax2))
			{			
				chi = g2*(Lens.cCs3*g2+Lens.cf);
				c = cPi*Lens.beta*Lens.sf;
				u = 1.0 + c*c*g2;

				c = cPi*Lens.sf*Lens.lambda*g2;
				sie = 0.25*c*c;

				c = cPi*Lens.beta*(Lens.Cs3*Lens.lambda2*g2-Lens.f);
				tie = c*c*g2;

				sti = exp(-(sie+tie)/u);

				Psix = fPsi_i[ixy][0]; Psiy = fPsi_i[ixy][1]; 
				chix = cos(chi); chiy = sin(chi);		 	

				fPsi_o[ixy][0] = gx = sti*(chix*Psix-chiy*Psiy); 
				fPsi_o[ixy][1] = gy = sti*(chix*Psiy+chiy*Psix);
			}
			else
			{
 				fPsi_o[ixy][0] = 0.0; 
				fPsi_o[ixy][1] = 0.0; 
			}
		}
}

/***************************************************************************/
/***************************************************************************/

void f_Row_2_Column_Format_CPU(sGP &GP, double *&MC_i, double *&MC_o)
{	
	int ix, iy, ixyi, ixyo;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			MC_o[ix*GP.ny+iy] = MC_i[iy*GP.nx+ix];
		}
}

void f_Column_2_Row_Format_CPU(sGP &GP, double *&MC_i, double *&MC_o)
{	
	int ix, iy, ixyi, ixyo;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			MC_o[iy*GP.nx+ix] = MC_i[ix*GP.ny+iy];
		}
}

void f_fftw_complex_2_sComplex_CPU(sGP &GP, fftw_complex *&MC_i, sComplex &MC_o)
{	
	int ix, iy, ixyi, ixyo;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixyi = iy*GP.nx+ix;
			ixyo = ix*GP.ny+iy;
			MC_o.real[ixyo] = MC_i[ixyi][0];
			MC_o.imag[ixyo] = MC_i[ixyi][1];
		}
}

void f_sComplex_2_fftw_complex_CPU(sGP &GP, sComplex &MC_i, fftw_complex *&MC_o)
{	
	int ix, iy, ixyi, ixyo;
	for(iy=0; iy<GP.ny; iy++)
		for(ix=0; ix<GP.nx; ix++)
		{
			ixyi = ix*GP.ny+iy;
			ixyo = iy*GP.nx+ix;
			MC_o[ixyo][0] = MC_i.real[ixyi];
			MC_o[ixyo][1] = MC_i.imag[ixyi];
		}
}
