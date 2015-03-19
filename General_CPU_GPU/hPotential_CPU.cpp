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

#include "hmathCPU.h"
#include <cstring>
#include "hConstTypes.h"
#include "hMT_AtomTypes_CPU.h"
#include "hQuadrature.h"
#include "hPotential_CPU.h"


// Constructor
cPotential_CPU::cPotential_CPU()
{
	PotPar = 0;
	MT_AtomTypes_CPU = 0;
	sigma = 0;
	nQ1 = 127;
	Q1.x = new double [nQ1];
	Q1.w = new double [nQ1];
	cQuadrature Quad;
	Quad.ReadQuadrature(1, nQ1, Q1);	// 1: int_0^infty f(x) dx - ExpSinh quadrature
	cl = cnl = 0;
}

// Destructor
cPotential_CPU::~cPotential_CPU()
{
	PotPar = 0;
	MT_AtomTypes_CPU = 0;
	sigma = 0;
	nQ1 = 0;
	delete [] Q1.x; Q1.x = 0;
	delete [] Q1.w; Q1.w = 0;
	cl = cnl = 0;
}

// Set RMS factor (Angstroms)
void cPotential_CPU::SetSigma(double sigmai)
{
	sigma = sigmai;
}

// Set Atom type
void cPotential_CPU::SetAtomTypes(int PotPar_i, cMT_AtomTypes_CPU *MT_AtomTypes_CPU_i)
{
	PotPar = PotPar_i;
	MT_AtomTypes_CPU = MT_AtomTypes_CPU_i;
}

/***************************************************************************/
// 1D Potential calculation (Vx, dVx) where dVx is the first derivative along x
void cPotential_CPU::Pot1Da(double &r, double *cl, double *cnl, double &f, double &df)
{

}

// 1D Potential calculation (Vx, dVx) where dVx is the first derivative along x
void cPotential_CPU::Pot1Dn(double &r, double *cl, double *cnl, double &f, double &df)
{

}

/***************************************************************************/
// 2D Potential calculation (VR, dVR) where dVR is the first derivative along R
void cPotential_CPU::Pot2Da(double &r, double *cl, double *cnl,double &f, double &df)
{
	int i;
	double icnl;
 double r2 = r*r;
	f = df = 0.0;
	switch(PotPar)
	{
		case 1:
			// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
			for(i=0; i<4; i++)
			{
				f += ft = cl[i]*exp(-cnl[i]*r2);
				df += -2.0*cnl[i]*r*ft;
			}
			break;
		case 2:
			// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
			for(i=0; i<5; i++)
			{
				f += ft = cl[i]*exp(-cnl[i]*r2);
				df += -2.0*cnl[i]*r*ft;
			}
			break;
		case 3:
			// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
			for(i=0; i<5; i++)
			{
				f += ft = cl[i]*exp(-cnl[i]*r2);
				df += -2.0*cnl[i]*r*ft;
			}
			break;
		case 4:
			// 4: Kirkland parameterization - 3 Yukawa + 3 Gaussians - [0, 12]			
			for(i=0; i<3;i++)
			{
				f += cl[i]*besskCPU(0, cnl[i]*r);
				df += -cl[i]*cnl[i]*besskCPU(1, cnl[i]*r);
			}	
			for(i=3; i<6;i++)
			{
				f += ft = cl[i]*exp(-cnl[i]*r2);
				df += -2.0*cnl[i]*r*ft;
			}
			break;
		case 5:
			// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]

			break;
		case 6:
			// 6: lobato parameterization - 5 hydrogen fe - [0, 12]
			for(i=0; i<5; i++)
			{
				k0 = besskCPU(0, cnl[i]*r); k1 = besskCPU(1, cnl[i]*r); icnl = 1.0/cnl[i];
				f += cl[i]*(2.0*icnl*k0 + r*k1);
				df += -cl[i]*(cnl[i]*r*k0 + 2.0*k1);
			}
			break;
	}
}

// 2D Potential calculation (VR, dVR) where dVR is the first derivative along R
void cPotential_CPU::Pot2Dn(double &r, double *cl, double *cnl, double &f, double &df)
{
//	int i, j, nm, np;
//	double M, h, f, alpha, beta;
//	double em, ep, t, dt, ut;
//	double gi, gi2, wgi, Vg, dVg;
//	double zi, wzi, r, r2, Vr, dVr;
//
//	VR = dVR = 0; h = 0.075; beta = 0.25; 
//	for(j=0; j<nQ1; j++){
//		zi = Q1.x[j]; wzi = Q1.w[j]; r2 = R*R + zi*zi; r = sqrt(r2);
//
//		f = c2Pi*r; M = cPi/(f*h); alpha = beta/sqrt(1.0+M*log(1.0+M)/c4Pi);
//	
//		ep = log(1e-10/M); t = log(-ep/beta + 1.0); t = log((2.0*t + log(t) - ep)/beta + 1.0);
//		np = 1 + (int)(log((2.0*t + log(t) - ep)/beta + 1.0)/h);
//
//		em = log(1e-20/M); t = log(-em/alpha + 1.0); t = log((2.0*t + log(t) - em)/alpha + 1.0);
//		nm = -1 - (int)(log((2.0*t + log(t) - em)/alpha + 1.0)/h);
//	
//		Vr = dVr = 0;
//		for(i=nm; i<=np; i++){
//			t = i*h;
//			if(t==0){
//				ut = 1.0/(2.0+alpha+beta);
//				gi = M*ut;
//				wgi = 0.5*M*h*(1.0+(alpha-beta)*ut*ut)*sin(f*gi);
//			}else{
//				em = exp(-t); ep = beta*(1.0/em-1); em = alpha*(em-1); 
//				dt = exp(-2*t+em-ep); ut = 1.0/(1.0-dt);
//				gi = M*t*ut;
//				wgi = M*h*(1.0-(2.0+alpha+beta+em+ep)*t*(ut-1))*ut;
//				if(i>5)
//					wgi *= sin(((i&1)?-1:1)*f*gi*dt);
//				else
//					wgi *= sin(f*gi);
//			}
//			fea(gi, Vg, dVg);
//			if(dw==0){
//				Vr += wgi*gi*Vg;
//				dVr += -wgi*gi*(3*Vg+gi*dVg);
//			}else{
//				gi2 = gi*gi; 
//				Vr += wgi*gi*Vg*exp(-dw*gi2);
//				dVr += -wgi*gi*((3-2*gi2*dw)*Vg+gi*dVg)*exp(-dw*gi2);
//			}
//		}
//
//		VR += Vr*wzi/r;
//		dVR += dVr*wzi/(r*r2);
//	}
//	VR *= 4.0*cPotf;
//	dVR *= 4.0*R*cPotf;
}

/***************************************************************************/
// 3D Potential calculation (Vr, dVr) where dVr is the first derivative along r
void cPotential_CPU::Pot3Da(double &r, double *cl, double *cnl, double &f, double &df)
{
	int i;
	double icnl;
	double ir = 1.0/r, r2 = r*r;

	f = df = 0.0;
	switch(PotPar)
	{
		case 1:
			// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
			for(i=0; i<4; i++)
			{
				f += ft = cl[i]*exp(-cnl[i]*r2);
				df += -2.0*cnl[i]*r*ft;
			}
			break;
		case 2:
			// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
			for(i=0; i<5; i++)
			{
				f += ft = cl[i]*exp(-cnl[i]*r2);
				df += -2.0*cnl[i]*r*ft;
			}
			break;
		case 3:
			// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
			for(i=0; i<5; i++)
			{
				f += ft = cl[i]*exp(-cnl[i]*r2);
				df += -2.0*cnl[i]*r*ft;
			}
			break;
		case 4:
			// 4: Kirkland parameterization - 3 Yukawa + 3 Gaussians - [0, 12]			
			for(i=0; i<3; i++)
			{
				f += ft = cl[i]*exp(-cnl[i]*r)*ir;
				df += -(cnl[i]+ ir)*ft;
			}

			for(i=3; i<6; i++)
			{
				f += ft = cl[i]*exp(-cnl[i]*r2);
				df += -2.0*cnl[i]*r*ft;
			}
			break;
		case 5:
			// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
			for(i=0; i<6; i++)
			{
				f += ft = cl[i]*erfcCPU(cnl[i]*r)*ir;
				df += (-2.0*cl[i]*cnl[i]*exp(-cnl[i]*cnl[i]*r2)/cPii2-ft)*ir;
			}
			break;
		case 6:
			// 6: Lobato parameterization - 5 Hydrogen fe - [0, 12]
			for(i=0; i<5; i++)
			{
				ft = cl[i]*exp(-cnl[i]*r)*ir; icnl = 1.0/cnl[i];
				f += ft*(2.0*icnl + r);
				df += -ft*(2.0*icnl*ir + 2.0 + cnl[i]*r);
			}
			break;
	}
}

// 3D Potential calculation (Vr, dVr) where dVr is the first derivative along r
void cPotential_CPU::Pot3Dn(double &r, double *cl, double *cnl, double &f, double &df)
{
//	int i, nm, np;
//	double M, h, f, alpha, beta;
//	double em, ep, t, dt, ut;
//	double gi, gi2, wgi, Vg, dVg;
//
//	f = c2Pi*r; h = 0.075; M = cPi/(f*h);
//	beta = 0.25; alpha = beta/sqrt(1.0+M*log(1.0+M)/c4Pi);
//	
//	ep = log(1e-10/M); t = log(-ep/beta + 1.0); t = log((2.0*t + log(t) - ep)/beta + 1.0);
//	np = 1 + (int)(log((2.0*t + log(t) - ep)/beta + 1.0)/h);
//
//	em = log(1e-20/M); t = log(-em/alpha + 1.0); t = log((2.0*t + log(t) - em)/alpha + 1.0);
//	nm = -1 - (int)(log((2.0*t + log(t) - em)/alpha + 1.0)/h);
//	
//	Vr = dVr = 0;
//	for(i=nm; i<=np; i++){
//		t = i*h;
//		if(t==0){
//			ut = 1.0/(2.0+alpha+beta);
//			gi = M*ut;
//			wgi = 0.5*(1.0+(alpha-beta)*ut*ut)*sin(f*gi);
//		}else{
//			em = exp(-t); ep = beta*(1.0/em-1); em = alpha*(em-1); 
//			dt = exp(-2*t+em-ep); ut = 1.0/(1.0-dt);
//			gi = M*t*ut;
//			wgi = (1.0-(2.0+alpha+beta+em+ep)*t*(ut-1))*ut;
//			if(i>5)
//				wgi *= sin(((i&1)?-1:1)*f*gi*dt);
//			else
//				wgi *= sin(f*gi);
//		}
//		gi2 = gi*gi;
//		fea(gi, Vg, dVg);
//		if(dw==0){
//			Vr += wgi*gi*Vg;
//			dVr += -wgi*gi*(3*Vg+gi*dVg);
//		}else{
//			Vr += wgi*gi*Vg*exp(-dw*gi2);
//			dVr += -wgi*gi*((3-2*dw*gi2)*Vg+gi*dVg)*exp(-dw*gi2);
//		}
//	}
//	Vr = 2.0*cPotf*M*h*Vr/r;
//	dVr = 2.0*cPotf*M*h*dVr/(r*r);
}

/***************************************************************************/
void cPotential_CPU::Vr(int IntType, int Dim, double r, double &f, double &df)
{
	switch(Dim)
	{
		case 1:
			if(~IntType)
			{
				Pot1Da(r, cl, cnl, f, df);
			}
			else
			{
				Pot1Dn(r, cl, cnl, f, df);
			}
			break;
		case 2:
			cl = MT_AtomTypes_CPU->cVR.cl;
			cnl = MT_AtomTypes_CPU->cVR.cnl;
			if(~IntType)
			{
				Pot2Da(r, cl, cnl, f, df);
			}
			else
			{
				Pot2Dn(r, cl, cnl, f, df);
			}
			break;
		case 3:
			cl = MT_AtomTypes_CPU->cVr.cl;
			cnl = MT_AtomTypes_CPU->cVr.cnl;
			if(~IntType)
			{
				Pot3Da(r, cl, cnl, f, df);
			}
			else
			{
				Pot3Dn(r, cl, cnl, f, df);
			}
			break;
	}

}

/***************************************************************************/
// 2D Potential calculation (VR, dVR) where dVR is the first derivative along R
void cPotential_CPU::Vr(int IntType, int Dim, int nr, double *r, double *f, double *df)
{
	for(int i=0; i<nr; i++)
	{
		Vr(IntType, Dim, r[i], f[i], df[i]);
	}
}

// Calculate the total Potential
void cPotential_CPU::Vr(int PotPar, int nAtoms, sAtoms *Atoms, int nMT_AtomTypes, cMT_AtomTypes_CPU *MT_AtomTypes_CPU, int IntType, int Dim, int nr, double *r, double factor, double *f, double *df)
{
	int i, j;
	double *ft, *dft; 
	double occ;

	ft = new double[nr]; 
	dft = new double[nr];

	// Initiazed potential
	for(i=0; i<nr; i++)
	{
		f[i] = df[i] = 0.0;
	}

	for(i=0; i<nMT_AtomTypes; i++)
	{
		occ = 0.0;
		for(j=0; j<nAtoms; j++)
		{
			if(Atoms[j].Z==MT_AtomTypes_CPU[i].Z)		
			{
				occ += Atoms[j].occ;
			}
		}

		if(occ>0)
		{
			SetAtomTypes(PotPar, &(MT_AtomTypes_CPU[i]));
			SetSigma(0.0);
			Vr(IntType, Dim, nr, r, ft, dft);
			// Add potential
			for(j=0; j<nr; j++)
			{
				f[j] += occ*factor*ft[j];
				df[j] += occ*factor*dft[j];
			}
		}
	}

	delete [] ft;
	delete [] dft;
}

/***************************************************************************/
double cPotential_CPU::AtomicRadius_rms(int Dim)
{
	double iVr, iVrn, Vrt;
	double ri, wi, Vri, dVri;
	double ra, rmin;

	if(MT_AtomTypes_CPU->cVr.cl[0]==0)
	{
		ra = 0.0;
	}
	else
	{
		rmin = MT_AtomTypes_CPU->rn_c;
		Vr(0, Dim, rmin, Vri, dVri);

		iVr = Vri*pow(rmin, Dim)/Dim;
		iVrn = Vri*pow(rmin, Dim+2)/(Dim+2);
		/************************************************************/
		for(int i=0; i<nQ1; i++)
		{
			ri = Q1.x[i] + rmin;
			wi = Q1.w[i];
			Vr(0, Dim, ri, Vri, dVri);		
			iVr += Vrt = wi*Vri*pow(ri, Dim-1);
			iVrn += Vrt*ri*ri;
		}
		ra = sqrt(iVrn/iVr);
	}
	return ra;
}

double cPotential_CPU::AtomicRadius_Cutoff(int Dim, double Vrl)
{
	double eps = 1e-08;
	double rmin, rmax;
	double rm, Vrm, dVrm;

	if((Vrl==0)||(MT_AtomTypes_CPU->cVr.cl[0]==0))
	{
		rmax = 0;
	}
	else
	{
		rmin = MT_AtomTypes_CPU->rn_c; 
		rmax = 50.0;
		Vr(0, Dim, rmin, Vrm, dVrm); 
		Vrm = std::abs(Vrm);
		while (std::abs(Vrm-Vrl)>eps)
		{
			rm = 0.5*(rmin+rmax);
			Vr(0, Dim, rm, Vrm, dVrm);
			Vrm = std::abs(Vrm);
			if(Vrm>Vrl)			
			{
				rmin = rm;
			}
			else
			{
				rmax = rm;
			}
		}
	}
	return rmax;
}
