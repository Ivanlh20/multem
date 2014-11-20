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

#include "hSAHF.h"
#include "hmathCPU.h"
#include <cstring>
#include "hConstTypes.h"
#include "hMT_General_CPU.h"
#include "hRandGen.h"
#include "hfxegTabData.h"
#include "blas.h"
#include "lapack.h"

#if !defined(_WIN32)
#define dgemm dgemm_
#endif

#if !defined(_WIN32)
#define dgesv dgesv_
#endif

#if !defined(_WIN32)
#define dgelsd dgelsd_
#endif

void cSAHF::freeMemory(){
	/******************************************/
	delete [] M; M = 0;
	delete [] Mt; Mt = 0;
	delete [] N; N = 0;
	delete [] Nt; Nt = 0;
	delete [] St; St = 0;
	delete [] Rm; Rm = 0;
	delete [] Rmt; Rmt = 0;
	delete [] Fm; Fm = 0;
	/******************************************/
	rcond = 0;
	rank = 0;
	lwork = 0;
	info = 0;
	Z = 0;
	Zk = 0;
	onei = 1;
	oned = 1;
	/******************************************/
	delete [] iwork; iwork = 0;
	delete [] work; work = 0;
	delete [] Sv; Sv = 0;
	delete [] ipiv; ipiv = 0;
	/******************************************/
	for(int i=0; i<nhmax; i++){
		clmin[i] = clmax[i] = cnlmin[i] = cnlmax[i] = cnll[i] = cnld[i] = 0;
		srcnl[i] = icnl[i] = cnl2[i] = gamma[i] = 0;
		clfx[i] = cnlfx[i] = clPr[i] = cnlPr[i] = clVr[i] = cnlVr[i] = 0;
	}

	for(int i=0; i<i; i++){
		g[i] = g2[i] = fxg[i] = feg[i] = 0;
		wxg[i] = fxgw[i] = weg[i] = fegw[i] = rb[i] = 0;
	}
}

cSAHF::~cSAHF(){
	freeMemory();
}

cSAHF::cSAHF(){
	/******************************************/
	M = 0;
	Mt = 0;
	N = 0;
	Nt = 0;
	St = 0;
	Rm = 0;
	Rmt = 0;
	Fm = 0;
	/******************************************/
	rcond = 0;
	rank = 0;
	lwork = 0;
	info = 0;
	Z = 0;
	Zk = 0;
	onei = 1;
	oned = 1;
	/******************************************/
	iwork = 0;
	work = 0;
	Sv = 0;
	ipiv = 0;
	/******************************************/
	for(int i=0; i<nhmax; i++){
		clmin[i] = clmax[i] = cnlmin[i] = cnlmax[i] = cnll[i] = cnld[i] = 0;
		srcnl[i] = icnl[i] = cnl2[i] = gamma[i] = 0;
		clfx[i] = cnlfx[i] = clPr[i] = cnlPr[i] = clVr[i] = cnlVr[i] = 0;
	}

	for(int i=0; i<i; i++){
		g[i] = g2[i] = fxg[i] = feg[i] = 0;
		wxg[i] = fxgw[i] = weg[i] = fegw[i] = rb[i] = 0;
	}
}

void cSAHF::SetInputData(int Zi, int typi, int nhti, int nhbi, int dgi, double *cnlmini, double *cnlmaxi){	
	freeMemory();

	Z = Zi;
	nht = nhti;	
	nhb = nhbi;
	nhbc = nht-nhb;
	memcpy(cnlmin, cnlmini, nht*cSizeofRD);
	memcpy(cnlmax, cnlmaxi, nht*cSizeofRD);

	// load fxg and feg 
	cfxegTabData fxegTabData;
	int ngt;
	fxegTabData.ReadTabData(Z, typi, dgi, ngt, g, g2, fxg, feg);
	ng = ngt;

	// weight and fegw
	double t;
	for (int i=0; i<ng; i++){	
		t = (double)Z/fxg[i]-1.0;
		weg[i] = sqrt((1.0 + t*t)/(2.0*ng*feg[0]*feg[0]));
		fegw[i] = weg[i]*feg[i];
	}

	// numerical grid for Vr and Pr
	nrb = 32; 
	sAtomTypesCPU AtomTypesCPU;
	f_SetAtomTypes(Z, 6, 0, stVrl, AtomTypesCPU);
	double rmin = AtomTypesCPU.rn_c, rmax = 20;	
	double dlnr = log(rmax/rmin)/(nrb-1);
	for (int i=0; i<nrb; i++)
		rb[i] = rmin*exp(i*dlnr);

	Zk = (double)Z/ca0;

	/******************************************/
	delete [] M; M = new double[nhb*nhb];
	delete [] Mt; Mt = new double[nhb*nhb];
	delete [] N; N = new double[nhb*nhbc];
	delete [] Nt; Nt = new double[nhbc*nhb];
	delete [] St; St = new double[nhb*ng];
	delete [] Rm; Rm = new double[ng*nhbc];
	delete [] Rmt; Rmt = new double[nhbc*ng];
	delete [] Fm; Fm = new double[ng*1];

	gamma[0] = (double)Z/c2Pi2a0;
	gamma[1] = 0.0;
	gamma[2] = feg[0];

	/******************************************/
	lwork = 256*256;
	work = new double[lwork];
	iwork = new ptrdiff_t[256*nht]; 
	Sv = new double[nht*nht];
	ipiv = new ptrdiff_t[nhb];
	rcond = -1;

	/******************************************/
	cTn = "N";
	onei = 1;
	oned = 1.0;
}

void cSAHF::RandGenPReset(){
	RandGenP.reset();
}

void cSAHF::RandGenPSetLimits(){
	for (int i=0; i<nht; i++){
		cnll[i] = cnlmin[i];	
		cnld[i]= cnlmax[i]-cnlmin[i];
	}
}

void cSAHF::RandGenPSetLimits(double *cnli, double dn){
	double di;
	for (int i=0; i<nht; i++){
		di = (cnlmax[i]-cnlmin[i])*dn;
		cnll[i] = MAX(cnlmin[i], cnli[i]-di);	
		cnld[i]= MIN(cnlmax[i], cnli[i]+di) - cnll[i];
	}
}

void cSAHF::RandGenPNext(double *cnlo){
	bool Ind;
	do{
		Ind = false;
		cnlo[0] = cnll[0] + cnld[0]*RandGenP.randu();
		for (int i=1; i<nht; i++){
			cnlo[i] = cnll[i] + cnld[i]*RandGenP.randu();
			if (cnlo[i-1] < cnlo[i]){
				Ind = true;
				break;
			}
		}
	} while (Ind);
}

bool cSAHF::BoundaryConditions(double *cl, double *cnl){
	int i, j;
	bool bc = true;
	double texp;
	double rho, rhob;
	double Vr, Vrb;

	for (j=0; j<nht; j++){
		clfx[j] = cl[j]*icnl[j];
		cnlfx[j] = cnl[j];

		clPr[j] = cl[j]/(cnl2[j]*srcnl[j]);
		cnlPr[j] = c2Pi/srcnl[j];

		clVr[j] = cl[j]/(cnl[j]*srcnl[j]);
		cnlVr[j] = cnlPr[j];
	}

	Vrb = rhob = 1e+20;
	for (i=0; i<nrb; i++){
		rho = Vr = 0;
		for (j=0; j<nht; j++){
			texp = exp(-cnlPr[j]*rb[i]);
			rho += clPr[j]*texp;
			Vr += clVr[j]*(2.0/(cnlVr[j]*rb[i]) + 1.0)*texp;
		}

		if ((Vr < 0)||(Vr > Vrb)||(rho < 0)||(rho > rhob)){
			bc = false;
			break;
		}else{
			Vrb = Vr;
			rhob = rho;
		}
	}

	return bc;
}

void cSAHF::GetLinealCoeff(double *cnl, double *cl, bool &bc){
	int i, j;
	double t;

	for (j=0; j<nht; j++){
		cnl2[j] = cnl[j]*cnl[j];
		srcnl[j] = sqrt(cnl[j]);
		icnl[j] = 1.0/cnl[j];
	}

	for (i=0; i<nhb; i++)
		for (j=0; j<nht; j++){
			switch (i){
				case 0:
					t = icnl[j];
					break;
				case 1:
					t = (cPi/srcnl[j] - Zk)/(cnl2[j]*srcnl[j]);
					break;
				case 2:
					t = 2.0;
					break;
			}
			if(j<nhb)
				M[j*nhb+i] = Mt[i*nhb+j] = t;
			else
				N[(j-nhb)*nhb+i] = Nt[i*nhbc+j-nhb] = -t;
		}

	for (j=0; j<ng; j++)
		for (i=0; i<nht; i++){
			t = 1.0/(1.0 + cnl[i]*g2[j]);
			t = weg[j]*(t + t*t);
			if(i<nhb)
				St[j*nhb+i] = t;
			else
				Rmt[j*nhbc+i-nhb] = t;
		}	
	/*************************************************************/	
	dgesv(&nhb, &ng, Mt, &nhb, ipiv, St, &nhb, &info);
	/*************************************************************/	
	dgemm(cTn, cTn, &nhbc, &ng, &nhb, &oned, Nt, &nhbc, St, &nhb, &oned, Rmt, &nhbc);
	/*************************************************************/
	for (j=0; j<ng; j++){
		t = 0;
		for (i=0; i<nhb; i++)
			t += gamma[i]*St[j*nhb+i];

		Fm[j] = fegw[j] - t;

		for (i=0; i<nhbc; i++)
			Rm[i*ng+j] = Rmt[j*nhbc+i];
	}
	/*************************************************************/		
	dgelsd(&ng, &nhbc, &onei, Rm, &ng, Fm, &ng, Sv, &rcond, &rank, work, &lwork, iwork, &info); 
	/*************************************************************/	

	for (i=0; i<nht; i++)
		cl[i] = (i<nhb)?gamma[i]:Fm[i-nhb];

	/*************************************************************/	
	dgemm(cTn, cTn, &nhb, &onei, &nhbc, &oned, N, &nhb, &cl[nhb], &nhbc, &oned, cl, &nhb);
	/*************************************************************/	
	dgesv(&nhb, &onei, M, &nhb, ipiv, cl, &nhb, &info);
	/*************************************************************/	

	bc = BoundaryConditions(cl, cnl);
}

void cSAHF::NextCoeff(double *clo, double *cnlo){
	bool bc;
	do {
		RandGenPNext(cnlo);		
		GetLinealCoeff(cnlo, clo, bc);
	} while (!bc);
}

void cSAHF::CostFun(double *cli, double *cnli, double &ee, double &sigx, double &sige){
	int i, j;
	double t, tx, te;
	double fxa, fea;

	sigx = sige = ee = 0;
	for (i=0; i<ng; i++){
		fea = fxa = 0;
		for (j=0; j<nht; j++){
			t = 1.0/(1.0 + cnli[j]*g2[i]);
			fxa += c2Pi2a0*cli[j]*t*t/cnli[j];
			fea += cli[j]*(t + t*t);
		}
		tx = (fxg[i]-fxa);	sigx += tx*tx; 
		te = (feg[i]-fea);	sige += te*te;
		ee += te*te*weg[i]*weg[i];		
	}
}

void cSAHF::GetCfegVr(double *cli, double *cnli, double *clfego, double *cnlfego, double *clVro, double *cnlVro){
	for (int i=0; i<nht; i++){			
		clfego[i] = cli[i];
		cnlfego[i] = cnli[i];	
		clVro[i] = cPotf*cPi2*cli[i]/pow(cnli[i], 1.5);
		cnlVro[i] = c2Pi/cnli[i];	
	}
}

void cSAHF::CopyData(int n, double eei, double *cli, double *cnli, double &eeo, double *clo, double *cnlo){
	eeo = eei;
	memcpy(clo, cli, n*cSizeofRD);
	memcpy(cnlo, cnli, n*cSizeofRD);
}

void cSAHF::CopyData(int n, double eei, double sigxi, double sigei, double *cli, double *cnli, double &eeo, double &sigxo, double &sigeo, double *clo, double *cnlo){
	eeo = eei;
	sigxo = sigxi;
	sigeo = sigei;
	memcpy(clo, cli, n*cSizeofRD);
	memcpy(cnlo, cnli, n*cSizeofRD);
}

void cSAHF::GetIniPoint(int ni, double *clfego, double *cnlfego, double *clVro, double *cnlVro, double *eeo, double *sigxo, double *sigeo){
	double cln[nhmax], cnln[nhmax];
	double een, sigxn, sigen;
	RandGenPReset();
	RandGenPSetLimits();
	for (int i=0; i<ni; i++){
		NextCoeff(cln, cnln);	
		CostFun(cln, cnln, een, sigxn, sigen);
		eeo[i] = een;
		sigxo[i] = sqrt(sigxn/ng);
		sigeo[i] = sqrt(sigen/ng);
		GetCfegVr(cln, cnln, &clfego[i*nht], &cnlfego[i*nht], &clVro[i*nht], &cnlVro[i*nht]);	
	}
}

void cSAHF::GetIniPoint(int ni, double *clo, double *cnlo, double &eeo, double &sigxo, double &sigeo, double &To){
	int i, j;
	double cln[nhmax], cnln[nhmax];
	double een, sigxn, sigen;
	double eea, eem, *eet;

	eet = new double [ni];
	eeo = eem = 1e+200;
	RandGenPReset();
	RandGenPSetLimits();
	for (i=0; i<ni; i++){
		NextCoeff(cln, cnln);	
		CostFun(cln, cnln, een, sigxn, sigen);

		if (een < eeo)
			CopyData(nht, een, sigxn, sigen, cln, cnln, eeo, sigxo, sigeo, clo, cnlo);

		eet[i] = een;
		if(een<eem)
			eem = een;
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
	
	To = -(eea-eem)/log(0.9);

	delete [] eet;
}

double cSAHF::dfactor(double p){
	double df, c = 2.0, p1 = 0.6, p2 = 0.4;
	if (p>p1)
		df = 1.0 + c*(p-p1)/p2;
	else if(p<p2)
		df = 1.0/(1.0 + c*(p2-p)/p2);
	else
		df = 1;
	return df;
}

void cSAHF::SimulatedAnnealing(int Z, int nht, int nhb, int dgi, double *cnlmini, double *cnlmaxi, int N, 
double *clfego, double *cnlfego, double *clVro, double *cnlVro, double *eeo, double *sigxo, double *sigeo){
	int typ = 2;
	SetInputData(Z, typ, nht, nhb, dgi, cnlmini, cnlmaxi);
	GetIniPoint(N, clfego, cnlfego, clVro, cnlVro, eeo, sigxo, sigeo);
}

void cSAHF::SimulatedAnnealing(int Z, int nht, int nhb, int dgi, double *cnlmini, double *cnlmaxi, int Nt, int Ns, double rT1, 
	double rT2, double *clfego, double *cnlfego, double *clVro, double *cnlVro, double *eeo, double *sigxo, double *sigeo){
	int i, j, typ = 2;
	double T0, Tmin, T, rT;
	double dn, dmax, p;

	double clopt[nhmax], cnlopt[nhmax];
	double eeopt, sigxopt, sigeopt;

	double cla[nhmax], cnla[nhmax];
	double eea, sigxa, sigea;	

	double cln[nhmax], cnln[nhmax];
	double een, sigxn, sigen;

	RandGenT.reset();
	
	SetInputData(Z, typ, nht, nhb, dgi, cnlmini, cnlmaxi);
	GetIniPoint(Nt*Ns, clopt, cnlopt, eeopt, sigxopt, sigeopt, T0);
	CopyData(nht, eeopt, sigxopt, sigeopt, clopt, cnlopt, eea, sigxa, sigea, cla, cnla);

	T = T0; Tmin = 1e-10;
	dmax = dn = 1.0;

	while ((T>Tmin)&&(dn>1e-8)){
		RandGenPReset();
		for (i=0; i<Nt; i++){
			p = 0;		
			RandGenPSetLimits(cnla, dn);
			for (j=0; j<Ns; j++){
				NextCoeff(cln, cnln);
				CostFun(cln, cnln, een, sigxn, sigen);
				if (een < eea){
					p = p + 1;
					CopyData(nht, een, sigxn, sigen, cln, cnln, eea, sigxa, sigea, cla, cnla);
					RandGenPSetLimits(cnla, dn);					
					if (een < eeopt)
						CopyData(nht, een, sigxn, sigen, cln, cnln, eeopt, sigxopt, sigeopt, clopt, cnlopt);
					
				}else if (exp(-(een-eea)/T) > RandGenT.randu()){
					p = p + 1;
					CopyData(nht, een, sigxn, sigen, cln, cnln, eea, sigxa, sigea, cla, cnla);
					RandGenPSetLimits(cnla, dn);					
				}
			}	
			dn = MIN(dn*dfactor(p/Ns), dmax);
			if (dn==0)
				dn = MIN(T, dmax);
		}

		rT = (dn==dmax)?rT1:rT2;

		T = rT*T;
		CopyData(nht, eeopt, sigxopt, sigeopt, clopt, cnlopt, eea, sigxa, sigea, cla, cnla);
	}
	GetCfegVr(clopt, cnlopt, clfego, cnlfego, clVro, cnlVro);
	sigxopt = sqrt(sigxopt/ng);	*sigxo = sigxopt;
	sigeopt = sqrt(sigeopt/ng);	*sigeo = sigeopt;
	*eeo = eeopt;
}