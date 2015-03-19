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

#ifndef hSAGaussianFit_H
#define hSAGaussianFit_H

#include "hmathCPU.h"
#include <cstring>
#include "hConstTypes.h"
#include "hRandGen.h"
#include "lapack.h"

#if !defined(_WIN32)
#define dgelsd dgelsd_
#define dgels dgels_
#endif

#define ngmax 32
#define nxmax 2048

class cSAGF{
	private:
		cRandGen RandGen, RandGena;
		ptrdiff_t rank, lwork, info, *iwork;
		double *work, rcond, *Sv;	
		ptrdiff_t nx, ng, ngh, onei;		
		double x[nxmax], x2[nxmax], y[nxmax], w[nxmax], wb[nxmax], *M, *Mt, *F;
		double clt[ngmax], cnlt[ngmax], cl[ngmax], cnl[ngmax], cln[ngmax], cnln[ngmax], clbc[ngmax], cnlbc[ngmax];
		double cnlmin[ngmax], cnlmax[ngmax], cnll[ngmax], cnld[ngmax], icnl2[ngmax];
		bool bc;
		char *cTn;
		inline void SetInputData(int nxi, double *xi, double *yi, int ngi);
		inline void GetMinMaxcnl();
		inline void RandGenSetLimits();
		inline void RandGenSetLimits(double *cnli, double dn);
		inline void RandGenNext(double *cnlo);

		inline double BoundaryConditions(double *cli, double *cnli, double &eeb, double &ee);
		inline void GetLinealCoeffCostFun(double *cnl, double *cl, bool &bc, double &eeb, double &ee);
		inline void NextCoeff(double *clo, double *cnlo, double &ee);

		inline void GetIniPoint(int ni, double *clo, double *cnlo, double &eeo, double &T);
		inline double dfactor(double p);
	public:
		inline void SimulatedAnnealing(int nxi, double *xi, double *yi, int ngi, int Nt, int Ns, double rT1, double rT2, double *a, double *b);
		cSAGF();
		~cSAGF();		
};

cSAGF::cSAGF()
{
	/***************************************************************************/
	M = 0;
	Mt = 0;
	F = 0;
	/***************************************************************************/
	iwork = 0;
	work = 0;
	Sv = 0;
	rcond = -1;
	cTn = "N";
	/***************************************************************************/
	RandGen.reset();
	/***************************************************************************/
	for(int i=0; i<ngmax; i++)
	{
		clt[i] = cnlt[i] = cl[i] = cnl[i] = cln[i] = cnln[i] = clbc[i] = cnlbc[i] = 0;
		cnlmin[i] = cnlmax[i] = cnll[i] = cnld[i] = 0;
	}
}

cSAGF::~cSAGF()
{
	cTn = "N";
	/***************************************************************************/
	delete [] M; M = 0;
	delete [] Mt; Mt = 0;
	delete [] F; F = 0;
	/***************************************************************************/
	delete [] iwork; iwork = 0;
	delete [] work; work = 0;
	delete [] Sv; Sv = 0;
}

inline void cSAGF::SetInputData(int nxi, double *xi, double *yi, int ngi)
{	
	nx = nxi;
	ng = ngi;
	/***************************************************************/
	for(int i=0; i<nx; i++)
	{
		x[i] = xi[i];
		x2[i] = x[i]*x[i];
		y[i] = yi[i];
	}
	/***************************************************************/
	delete [] F; F = new double[nx*1];
	delete [] M; M = new double[nx*ng];
	delete [] Mt; Mt = new double[ng*nx];
	/***************************************************************/
	lwork = 8*nx*ng;
	delete [] work; work = new double[lwork];
	delete [] iwork; iwork = new ptrdiff_t[8*nx*ng]; 
	delete [] Sv; Sv = new double[ng*ng];
	rcond = -1;
	onei = 1;
	cTn = "N";
	/***************************************************************/
}

inline void cSAGF::GetMinMaxcnl()
{
	int i;
	double t0, t1, t2;
	double dx, smin, smax, ff;

	t0 = t1 = t2 = 0;
	for(i=0; i<(nx-1); i++)
	{
		dx = x[i+1]-x[i];
		t0 += (y[i+1] + y[i])*dx;
		t1 += (x2[i+1]*y[i+1] + x2[i]*y[i])*dx;
		t2 += (SQR(x2[i+1])*y[i+1]+SQR(x2[i])*y[i])*dx;
	}

	smin = sqrt(t1/t0);
	smax = sqrt(t2/t1);
	for(i=0; i<ng; i++)
	{
		cnlmin[i] = smin/30.0;
		cnlmax[i] = 3.0*smax;
	}
	/*******************************************************/
	ff = 1.0/(0.5*(t0+t1));
	for(i=1; i<(nx-1); i++)
		wb[i] = 0.5*(x[i+1]-x[i-1]);

	wb[0] = 10*wb[1]; wb[nx-1] = 5*wb[nx-2];

	for(i=0; i<nx; i++)
	{
		w[i] = ff*wb[i]*sqrt(1.0+x2[i]);
	}
}

inline void cSAGF::RandGenSetLimits()
{
	for(int i=0; i<ng; i++)
	{
		cnll[i] = cnlmin[i];	
		cnld[i] = cnlmax[i]-cnlmin[i];
	}
}

inline void cSAGF::RandGenSetLimits(double *cnli, double dn)
{
	double di;
	for(int i=0; i<ng; i++)
	{
		di = (cnlmax[i]-cnlmin[i])*dn;
		cnll[i] = MAX(cnlmin[i], cnli[i]-di);	
		cnld[i] = MIN(cnlmax[i], cnli[i]+di) - cnll[i];
	}
}

inline void cSAGF::RandGenNext(double *cnlo)
{
	bool Ind;
	do{
		Ind = false;
		cnlo[0] = cnll[0] + cnld[0]*RandGen.randu();
		for(int i=1; i<ng; i++)
		{
			cnlo[i] = cnll[i] + cnld[i]*RandGen.randu();
			if(cnlo[i-1] < cnlo[i])
			{
				Ind = true;
				break;
			}
		}
	} while (Ind);
}

/*inline double cSAGF::BoundaryConditions(double *cli, double *cnli)
{
	int i, j;
	double ya;
	double bc = true;

	for(i= (nx-1); i>=0; i--)
	{
		ya = 0;
		for(j=0; j<ng; j++)
			ya += cli[j]*Mt[i*ng+j];

		if(ya<0)
		{
			bc = false;
			break;
		}
	}

	return bc;
}*/

inline double cSAGF::BoundaryConditions(double *cli, double *cnli, double &eeb, double &ee)
{
	int i, j;
	double yab, ya, t;
	double bc = true;

	ee = eeb = 0;
	yab = 1e+200;

	for(i=0; i<nx; i++)
	{
		ya = 0;
		for(j=0; j<ng; j++)
			ya += cli[j]*Mt[i*ng+j];

		t = wb[i]*(y[i]-ya); 
		eeb += t*t;

		t = w[i]*(y[i]-ya); 
		ee += t*t;

		if(bc)
		{
			if((yab<ya)||(ya<0))
			{
				bc = false;
			}
			else
			{
				yab = ya;
			}
		}
	}

	return bc;
}

inline void cSAGF::GetLinealCoeffCostFun(double *cnl, double *cl, bool &bc, double &eeb, double &ee)
{
	int i, j;
	double t;
	for(j=0; j<ng; j++)
		icnl2[j] = 0.5/SQR(cnl[j]);

	for(i=0; i<nx; i++)
	{
		F[i] = w[i]*y[i];
		for(j=0; j<ng; j++)
		{
			Mt[i*ng + j] = t = exp(-icnl2[j]*x2[i]);
			M[j*nx + i] = w[i]*t;
		}
	}

	// dgelsd(&nx, &ng, &onei, M, &nx, F, &nx, Sv, &rcond, &rank, work, &lwork, iwork, &info); 
	dgels(cTn, &nx, &ng, &onei, M, &nx, F, &nx, work, &lwork, &info); 
	
	for(j=0; j<ng; j++)
		cl[j] = F[j];

	bc = BoundaryConditions(cl, cnl, eeb, ee);
}

inline void cSAGF::NextCoeff(double *clo, double *cnlo, double &ee)
{
	int c=0;
	double eet, eeb, eebmin = 1e+200;
	do
	{
		RandGenNext(cnlo);		
		GetLinealCoeffCostFun(cnlo, clo, bc, eeb, ee);
		if(eeb<eebmin)
		{
			eebmin = eeb;
			eet = ee;
			memcpy(clbc, clo, ng*cSizeofRD);
			memcpy(cnlbc, cnlo, ng*cSizeofRD);
		}
		c++;
	} while ((!bc)&&(c<25));

	if(c==25)
	{
		ee = eet;
		memcpy(clo, clbc, ng*cSizeofRD);
		memcpy(cnlo, cnlbc, ng*cSizeofRD);
	}
}

inline void cSAGF::GetIniPoint(int ni, double *clo, double *cnlo, double &eeo, double &T)
{
	int i, j;
	double een, eea, *eet;

	eet = new double [ni];
	eeo = 1e+200;
	RandGen.reset();
	RandGenSetLimits();
	for(i=0; i<ni; i++)
	{
		NextCoeff(cln, cnln, een);
		if(een < eeo)
		{
			eeo = een;
			memcpy(clo, cln, ng*cSizeofRD);
			memcpy(cnlo, cnln, ng*cSizeofRD);			
		}
		eet[i] = een;
	}

	// Sort by energy
	for(i=0; i<ni; i++)
	{
		for(j=i+1; j<ni; j++)
		{
			if(eet[j] < eet[i])
			{
				eea = eet[i];
				eet[i] = eet[j];
				eet[j] = eea;
			}
		}
	}
	
	T = -(eet[ni/2]-eeo)/log(0.9);

	delete [] eet;
}

inline double cSAGF::dfactor(double p)
{
	double df, c = 2.0, p1 = 0.6, p2 = 0.4;
	if(p>p1)
	{
		df = 1.0 + c*(p-p1)/p2;
	}
	else if(p<p2)
	{
		df = 1.0/(1.0 + c*(p2-p)/p2);
	}
	else
	{
		df = 1;
	}
	return df;
}

inline void cSAGF::SimulatedAnnealing(int nxi, double *xi, double *yi, int ngi, int Nt, int Ns, double rT1, double rT2, double *a, double *b)
{
	int i, j;
	double T0, T, rT;
	double dn, dmax;
	double p, eet, ee, een;	

	RandGena.reset();
	SetInputData(nxi, xi, yi, ngi);
	GetMinMaxcnl();
	GetIniPoint(Nt*Ns, clt, cnlt, eet, T0);

	T = T0;
	dmax = dn = 1.0;

	ee = eet;
	memcpy(cl, clt, ng*cSizeofRD);
	memcpy(cnl, cnlt, ng*cSizeofRD);	
	RandGenSetLimits();
	while (dn>1e-6)
	{	
		RandGen.reset();
		for(i=0; i<Nt; i++)
		{
			p = 0;		
			RandGenSetLimits(cnl, dn);
			for(j=0; j<Ns; j++)
			{
				NextCoeff(cln, cnln, een);
				if(een < ee)
				{
					p = p + 1;
					ee = een;
					memcpy(cl, cln, ng*cSizeofRD);
					memcpy(cnl, cnln, ng*cSizeofRD);
					RandGenSetLimits(cnl, dn);
					if(een < eet)
					{
						eet = een;
						memcpy(clt, cln, ng*cSizeofRD);
						memcpy(cnlt, cnln, ng*cSizeofRD);
					}
				}
				else if(exp(-(een-ee)/T) > RandGena.randu())
				{
					p = p + 1;
					ee = een;
					memcpy(cl, cln, ng*cSizeofRD);
					memcpy(cnl, cnln, ng*cSizeofRD);
					RandGenSetLimits(cnl, dn);
				}
			}	
			dn = MIN(dn*dfactor(p/Ns), dmax);
			if(dn==0)			
			{
				dn = MIN(T, dmax);
			}
		}

		if(dn==dmax)
		{
			rT = rT1;
		}
		else
		{
			rT = rT2;
		}

		T = rT*T;
		ee = eet;
		memcpy(cl, clt, ng*cSizeofRD);
		memcpy(cnl, cnlt, ng*cSizeofRD);
	}
	for(j = 0; j<ng; j++)
	{
		cnlt[j] = 0.5/SQR(cnlt[j]);
	}

	memcpy(a, clt, ng*cSizeofRD);
	memcpy(b, cnlt, ng*cSizeofRD);
}

#endif