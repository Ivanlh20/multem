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

#ifndef hSAGaussianExpFit_H
#define hSAGaussianExpFit_H

#include "hmathCPU.h"
#include <cstring>
#include "hConstTypes.h"
#include "hRandGen.h"
#include "blas.h"
#include "lapack.h"

#if !defined(_WIN32)
#define dgelsd dgelsd_
#endif

class cSAGEF{
	private:
		cRandGen RandGen, RandGena;
		ptrdiff_t rank, lwork, info, *iwork;
		double *work, rcond, *Sv;	
		ptrdiff_t nx, ng, ne, nt, onei;		
		double *x, *x2, *y, *M, *Mt, *F;
		double clt[32], cnlt[32], cl[32], cnl[32], cln[32], cnln[32];
		double cnlmin[32], cnlmax[32], cnll[32], cnld[32], icnl2[32];
		bool bc;
		inline void GetMinMaxcnl();
		inline void RandGenSetLimits();
		inline void RandGenSetLimits(double *cnli, double dn);
		inline void RandGenNext(double *cnlo);
		inline void SetInputData(int nxi, double *xi, double *yi, int ngi, int nei);
		inline void NextCoeff(double *clo, double *cnlo);
		inline void GetLinealCoeff(double *cnl, double *cl, bool &bc);
		inline void CostFun(double *cli, double *cnli, double &ee);
		inline void GetIniPoint(int ni, double *clo, double *cnlo, double &eeo, double &T);
		inline double dfactor(double p);
	public:
		inline void SimulatedAnnealing(int nxi, double *xi, double *yi, int ngi, int nei, int Nt, int Ns, double rT1, double rT2, double *a, double *b);
		cSAGEF();
		~cSAGEF();		
};

cSAGEF::cSAGEF(){
	x = 0;
	x2 = 0;
	y = 0;
	M = 0;
	Mt = 0;
	F = 0;
	/******************************************/
	iwork = 0;
	work = 0;
	Sv = 0;
	rcond = -1;
	/******************************************/
	RandGen.reset();
	/******************************************/
	for (int i=0; i<32; i++){
		clt[i] = cnlt[i] = cl[i] = cnl[i] = cln[i] = cnln[i] = 0;
		cnlmin[i] = cnlmax[i] = cnll[i] = cnld[i] = icnl2[i] = 0;
	}
}

cSAGEF::~cSAGEF(){
	x = 0;
	y = 0;
	/******************************************/
	delete [] x2; x2 = 0;
	delete [] M; M = 0;
	delete [] Mt; Mt = 0;
	delete [] F; F = 0;
	/******************************************/
	delete [] iwork; iwork = 0;
	delete [] work; work = 0;
	delete [] Sv; Sv = 0;
}

inline void cSAGEF::GetMinMaxcnl(){
	double t1, t2;
	double dx, sm;

	t1 = t2 = 0;
	for (int i=0; i<(nx-1); i++){
		dx = x[i+1]-x[i];
		t1 += (y[i+1]+y[i])*dx;
		t2 += (SQR(x[i+1])*y[i+1]+SQR(x[i])*y[i])*dx;
	}
	sm = sqrt(t2/t1);
	for (int i=0; i<nt; i++){
		cnlmin[i] = sm*1e-04;
		cnlmax[i] = 10*sm;
	}
}

inline void cSAGEF::RandGenSetLimits(){
	for (int i=0; i<nt; i++){
		cnll[i] = cnlmin[i];	
		cnld[i] = cnlmax[i]-cnlmin[i];
	}
}

inline void cSAGEF::RandGenSetLimits(double *cnli, double dn){
	double di;
	for (int i=0; i<nt; i++){
		di = (cnlmax[i]-cnlmin[i])*dn;
		cnll[i] = MAX(cnlmin[i], cnli[i]-di);	
		cnld[i] = MIN(cnlmax[i], cnli[i]+di) - cnll[i];
	}
}

inline void cSAGEF::RandGenNext(double *cnlo){
	int i;
	bool Ind;
	do{
		Ind = false;
		cnlo[0] = cnll[0] + cnld[0]*RandGen.randu();
		for (i=1; i<ng; i++){
			cnlo[i] = cnll[i] + cnld[i]*RandGen.randu();
			if (cnlo[i-1] < cnlo[i]){
				Ind = true;
				break;
			}
		}
	} while (Ind);

	do{
		Ind = false;
		cnlo[ng] = cnll[ng] + cnld[ng]*RandGen.randu();
		for (i = ng+1; i<nt; i++){
			cnlo[i] = cnll[i] + cnld[i]*RandGen.randu();
			if (cnlo[i-1] < cnlo[i]){
				Ind = true;
				break;
			}
		}
	} while (Ind);
}

inline void cSAGEF::SetInputData(int nxi, double *xi, double *yi, int ngi, int nei){
	nx = nxi;
	x = xi;
	y = yi;
	ng = ngi;
	ne = nei;
	nt = ng + ne;
	/****************************************************************/
	delete [] x2; x2 = new double[nx*1];
	for (int i=0; i<nx; i++)
		x2[i] = x[i]*x[i];
	/****************************************************************/
	delete [] F; F = new double[nx*1];
	delete [] M; M = new double[nx*nt];
	delete [] Mt; Mt = new double[nt*nx];
	/****************************************************************/
	lwork = 256*256;
	delete [] work; work = new double[lwork];
	delete [] iwork; iwork = new ptrdiff_t[256*nt]; 
	delete [] Sv; Sv = new double[nt*nt];
	rcond = -1;
	onei = 1;
	/****************************************************************/
}

inline void cSAGEF::NextCoeff(double *clo, double *cnlo){
	do {
		RandGenNext(cnlo);		
		GetLinealCoeff(cnlo, clo, bc);
	} while (!bc);
}

inline void cSAGEF::GetLinealCoeff(double *cnl, double *cl, bool &bc){
	int i, j;
	
	for (j=0; j<nt; j++)
		icnl2[j] = (j<ng)?(0.5/SQR(cnl[j])):(c2i2/cnl[j]);

	for (i=0; i<nx; i++){
		F[i] = y[i];
		for (j=0; j<nt; j++)		
			Mt[i*nt + j] = M[j*nx + i] = (j<ng)?(exp(-icnl2[j]*x2[i])):(exp(-icnl2[j]*x[i]));
	}

	dgelsd(&nx, &nt, &onei, M, &nx, F, &nx, Sv, &rcond, &rank, work, &lwork, iwork, &info); 

	bc = true;
	for (i=0; i<nt; i++)
		cl[i] = F[i];
}

inline void cSAGEF::CostFun(double *cli, double *cnli, double &ee){
	int i, j;
	double ya;

	ee = 0;
	for (i=0; i<nx; i++){
		ya = 0;
		for (j=0; j<nt; j++)
			ya += cli[j]*Mt[i*nt+j];
		ee += (y[i]-ya)*(y[i]-ya);		
	}
}

inline void cSAGEF::GetIniPoint(int ni, double *clo, double *cnlo, double &eeo, double &T){
	int i, j;
	double een, eea, *eet;

	eet = new double [ni];
	eeo = 1e+200;
	RandGen.reset();
	RandGenSetLimits();
	for (i=0; i<ni; i++){
		NextCoeff(cln, cnln);	
		CostFun(cln, cnln, een);
		if (een < eeo){
			eeo = een;
			memcpy(clo, cln, nt*cSizeofRD);
			memcpy(cnlo, cnln, nt*cSizeofRD);			
		}
		eet[i] = een;
	}

	// Sort by energy
	for(i=0; i<ni; i++)
		for (j=i+1; j<ni; j++)
			if (eet[j] < eet[i]){
				eea = eet[i];
				eet[i] = eet[j];
				eet[j] = eea;
			}
	
	eea = eet[ni/2];
	
	T = -(eea - eeo)/log(0.9);

	delete [] eet;
}

inline double cSAGEF::dfactor(double p){
	double df, c = 2.0, p1 = 0.6, p2 = 0.4;
	if (p>p1)
		df = 1.0 + c*(p-p1)/p2;
	else if(p<p2)
		df = 1.0/(1.0 + c*(p2-p)/p2);
	else
		df = 1;
	return df;
}

inline void cSAGEF::SimulatedAnnealing(int nxi, double *xi, double *yi, int ngi, int nei, int Nt, int Ns, double rT1, double rT2, double *a, double *b){
	int i, j;
	double T0, Tmin, T, rT;
	double dn, dmax;
	double p, eet, ee, een;	

	RandGena.reset();

	SetInputData(nxi, xi, yi, ngi, nei);

	GetMinMaxcnl();

	GetIniPoint(Nt*Ns, clt, cnlt, eet, T0);

	T = T0; Tmin = 1e-9;
	dmax = dn = 1.0;

	ee = eet;
	memcpy(cl, clt, nt*cSizeofRD);
	memcpy(cnl, cnlt, nt*cSizeofRD);

	RandGenSetLimits();
	while ((T>Tmin)&&(dn>1e-5)){
		RandGen.reset();
		for (i=0; i<Nt; i++){
			p = 0;		
			RandGenSetLimits(cnl, dn);
			for (j=0; j<Ns; j++){
				NextCoeff(cln, cnln);
				CostFun(cln, cnln, een);
				if (een < ee){
					p = p + 1;
					ee = een;
					memcpy(cl, cln, nt*cSizeofRD);
					memcpy(cnl, cnln, nt*cSizeofRD);
					RandGenSetLimits(cnl, dn);
					if (een < eet){
						eet = een;
						memcpy(clt, cln, nt*cSizeofRD);
						memcpy(cnlt, cnln, nt*cSizeofRD);
					}
				}else if (exp(-(een-ee)/T) > RandGena.randu()){
					p = p + 1;
					ee = een;
					memcpy(cl, cln, nt*cSizeofRD);
					memcpy(cnl, cnln, nt*cSizeofRD);
					RandGenSetLimits(cnl, dn);
				}
			}	
			dn = MIN(dn*dfactor(p/Ns), dmax);
			if (dn==0)
				dn = MIN(T, dmax);
		}

		if (dn==dmax)
			rT = rT1;
		else
			rT = rT2;

		T = rT*T;
		ee = eet;
		memcpy(cl, clt, nt*cSizeofRD);
		memcpy(cnl, cnlt, nt*cSizeofRD);
	}
	for (j=0; j<nt; j++)
		cnlt[j] = (j<ng)?(0.5/SQR(cnlt[j])):(c2i2/cnlt[j]);

	memcpy(a, clt, nt*cSizeofRD);
	memcpy(b, cnlt, nt*cSizeofRD);
}

#endif