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
#include "hConstTypes.h"
#include "hGeneral_CPU.h"

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
		d = std::abs(c1*x0+c2*y0+c3);
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
			d = std::abs(s_i[i]-x[ix]);
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
			d = std::abs(s_i[i]-x[ix]);
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

// Build Grid
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
