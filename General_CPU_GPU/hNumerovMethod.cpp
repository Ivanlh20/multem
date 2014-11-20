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

#include "math.h"
#include <cstring>
#include "hConstTypes.h"
#include "hNumerovMethod.h"

void cNumerov::ReadInputdata(int nri, double *ri, double *Vri)
{ 
	nr = nri;
	// allocate Space
	delete [] r;	r = new double [nr];
	delete [] ro;	ro = new double [nr];
	delete [] sqr;	sqr = new double [nr];
	delete [] r2;	r2 = new double [nr];
	delete [] Vr;	Vr = new double [nr];
	delete [] f;	f = new double [nr]; 
	delete [] y;	y = new double [nr]; 
	delete [] Vrt;	Vrt = new double [nr]; 

	// Convert2AtomicUnits:	hb = 1, me = 1, qe = 1.
	 for (int i=0; i<nr; i++){
		ro[i] = ri[i];
		r[i] = ri[i]/ca0;
		sqr[i] = sqrt(r[i]);
		r2[i] = r[i]*r[i];
		Vr[i] = -Vri[i]/cHa;
	 }

	/****************Useful constants****************/
	dx = log(r[nr-1]/r[0])/double(nr-1);
	dxf = dx*dx/12;
	def = dx/(2.0*dxf);
}
 
void cNumerov::SolveRadSchEq(int ni, double gamma, int Dim, double *aEner, double *aPsi) 
{
	int i, j, k, c;
	n = ni;
	c = 0;
	memcpy(aPsi, ro, cSizeofRD*nr);
	for (n=1; n<=ni; n++){
		for (l=0; l<n; l++){
			nodes = n-l-1;
			switch (Dim){
				case 1:
					/******************************************/
					break;
				case 2:
					lh2 = 0.5*l*l;
					// Determination of the wave-function in the first two points
					y0 = (1-r[0]/(l+1.0))*pow(r[0], l);
					y1 = (1-r[1]/(l+1.0))*pow(r[1], l);
					break;
				case 3:
					lh2 = 0.5*(l+0.5)*(l+0.5);
					// Determination of the wave-function in the first two points
					y0 = (1-r[0]/(l+1.0))*pow(r[0], l+0.5);
					y1 = (1-r[1]/(l+1.0))*pow(r[1], l+0.5);
					break;
			}
			// Set initial lower and upper bounds to the eigenvalue
			elw = eup = Vr[nr-1];
			for (i=0; i<nr; i++){
				Vrt[i] = Vr[i] + lh2/r2[i];
				elw = MIN(elw, Vrt[i]);
			}
			e = 0.5*(elw + eup);

			for (k=0; k<Ni; k++){
				icl = -1;
				f[0] = 2.0*dxf*r2[0]*(e-Vrt[0]);
				for (i = 1; i < nr; i++){
					f[i] = 2.0*dxf*r2[i]*(e-Vrt[i]); 

					if (f[i] == 0.0) 
						f[i] = 1e-20; 

					if (f[i] != _copysign(f[i], f[i - 1])) 
						icl = i; 
				}

				if ((icl<0) || (icl>=nr-3)){
					eup = e;
					e = 0.5*(eup+elw);
					continue;
				}

				for (i=0; i<nr; i++)
					f[i] = 1.0 + f[i];

				// Outward integration condition
				y[0] = y0; y[1] = y1;
				// Outward integration, count number of crossings
				ncross = 0;
				for (i=1; i<=icl-1; i++){
					y[i+1] = ((12.0-10.0*f[i])*y[i]-f[i-1]*y[i-1])/f[i+1];
					if (y[i] != _copysign(y[i],y[i+1]))
					++ncross;
				}
				yicl = y[icl];

				// Check number of crossings
				if (ncross != nodes){
					if (ncross > nodes)
						eup = e;
					else
						elw = e;
 
					e = 0.5*(eup + elw);
					continue;
				}
 
				// Determination of the wave-function in the last two points, assuming y(Nr+1) = 0 and y(Nr) = dx
				y[nr-1] = dx;
				y[nr-2] = (12.0-10.0*f[nr-1])*y[nr-1]/f[nr-2];
				// Inward integration
				for (i = nr-2; i >= icl+1; i--){
					y[i-1] = ((12.0-10.0*f[i])*y[i]-f[i+1]*y[i+1])/f[i-1];
					if (y[i-1] > 1e10)
						for (j = nr-1; j >= i-1; j--)
							y[j] /= y[i-1];
				}

				// Rescale function to match at the classical turning point (icl)
				yicl /= y[icl];
				for (i = icl; i < nr; i++)
					y[i] *= yicl;

				// Normalize on the [-xmax,xmax] segment 
				norm = 0.;
				for (i=0; i<nr; i++)
					norm += y[i]*y[i]*r2[i]*dx;

				norm = sqrt(norm);
				for (i=0; i<nr; i++)
					y[i] /= norm;

				// find the value of the cusp at the matching point (icl)
				i = icl;
				ycusp = (y[icl-1]*f[icl-1]+f[icl+1]*y[icl+1]+10.0*f[icl]*y[icl])/12.0;
				dfcusp = f[icl]*(y[icl]/ycusp-1.0);

				// eigenvalue update using perturbation theory
				de = dfcusp*ycusp*ycusp*def;
				if (de > 0.0)
					elw = e;
				else if (de < 0.0)
					eup = e;

				// prevent e to go out of bounds, i.e. e > eup or e < elw (might happen far from convergence)
				e = MAX(MIN(e+de, eup), elw);

				// convergence not achieved
				if (ABS(de) < eps)
					break;
			}
			aEner[c++] = e*cHa/gamma;
			switch (Dim){
				case 1:
					/******************************************/
					break;
				case 2:
					/******************************************/
					break;
				case 3:
					for (i=0; i<nr; i++)
						y[i] /= sqr[i];
					break;
			}

			memcpy(&aPsi[c*nr], y, cSizeofRD*nr);
		}
	}
}

cNumerov::cNumerov()
{ 
	r = 0;
	ro = 0;
	sqr = 0;
	r2 = 0;
	Vr = 0;
	f = 0;
	y = 0;
	Vrt = 0;

	eps = 1e-10;
	Ni = 200;
}

cNumerov::~cNumerov()
{ 
	delete [] r;
	delete [] ro;
	delete [] sqr;
	delete [] r2;
	delete [] Vr;
	delete [] f;
	delete [] y;
	delete [] Vrt;
}