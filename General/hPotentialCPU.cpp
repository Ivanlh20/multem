#include "hmathCPU.h"
#include <cmath>
#include <cstring>
#include "hConstTypes.h"
#include "hPotentialCPU.h"
#include "hQuadrature.h"

// Constructor
cPotentialCPU::cPotentialCPU(){
	PotPar = 0;
	sigma = 0;
	nQ1 = 127;
	Q1.x = new double [nQ1];
	Q1.w = new double [nQ1];
	cQuadrature Quad;
	Quad.ReadQuadrature(1, nQ1, Q1);	// 1: int_0^infty f(x) dx - ExpSinh quadrature
}

// Destructor
cPotentialCPU::~cPotentialCPU(){
	PotPar = 0;
	sigma = 0;
	nQ1 = 0;
	delete [] Q1.x; Q1.x = 0;
	delete [] Q1.w; Q1.w = 0;
}

// Set RMS factor (Angstroms)
void cPotentialCPU::SetSigma(double sigmai){
	sigma = sigmai;
}

// Set Atom type
void cPotentialCPU::SetAtomTypes(int PotPari, sAtomTypesCPU AtomTypesi){
	PotPar = PotPari;
	std::memcpy(&AtomTypes, &AtomTypesi, sizeof(sAtomTypesCPU));
}

/************************************************************************************/
// 1D Potential calculation (Vx, dVx) where dVx is the first derivative along x
void cPotentialCPU::Pot1Da(double &r, double &f, double &df){

}

// 1D Potential calculation (Vx, dVx) where dVx is the first derivative along x
void cPotentialCPU::Pot1Dn(double &r, double &f, double &df){

}

/************************************************************************************/
// 2D Potential calculation (VR, dVR) where dVR is the first derivative along R
void cPotentialCPU::Pot2Da(double &r, double &f, double &df){
	int i;
	r2 = r*r;
	f = df = 0.0;
	switch (PotPar){
		case 1:
			// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
			for (i=0; i<4; i++){
				cl = AtomTypes.cVR.cl[i]; cnl = AtomTypes.cVR.cnl[i];
				f += ft = cl*exp(-cnl*r2);
				df += -2.0*cnl*r*ft;
			}
			break;
		case 2:
			// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
			for (i=0; i<5; i++){
				cl = AtomTypes.cVR.cl[i]; cnl = AtomTypes.cVR.cnl[i];
				f += ft = cl*exp(-cnl*r2);
				df += -2.0*cnl*r*ft;
			}
			break;
		case 3:
			// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
			for (i=0; i<5; i++){
				cl = AtomTypes.cVR.cl[i]; cnl = AtomTypes.cVR.cnl[i];
				f += ft = cl*exp(-cnl*r2);
				df += -2.0*cnl*r*ft;
			}
			break;
		case 4:
			// 4: Kirkland parameterization - 3 Yukawa + 3 Gaussians - [0, 12]			
			for (i=0; i<3;i++){
				cl = AtomTypes.cVR.cl[i]; cnl = AtomTypes.cVR.cnl[i];
				f += cl*besskCPU(0, cnl*r);
				df += -cl*cnl*besskCPU(1, cnl*r);
			}	
			for (i=3; i<6;i++){
				cl = AtomTypes.cVR.cl[i]; cnl = AtomTypes.cVR.cnl[i];
				f += ft = cl*exp(-cnl*r2);
				df += -2.0*cnl*r*ft;
			}
			break;
		case 5:
			// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]

			break;
		case 6:
			// 6: lobato parameterization - 5 hydrogen fe - [0, 12]
			for (i=0; i<5; i++){
				cl = AtomTypes.cVR.cl[i]; cnl = AtomTypes.cVR.cnl[i]; icnl = 1.0/cnl; 
				k0 = besskCPU(0, cnl*r); k1 = besskCPU(1, cnl*r);
				f += cl*(2.0*icnl*k0 + r*k1);
				df += -cl*(cnl*r*k0 + 2.0*k1);
			}
			break;
	}
}

// 2D Potential calculation (VR, dVR) where dVR is the first derivative along R
void cPotentialCPU::Pot2Dn(double &r, double &f, double &df){
//	int i, j, nm, np;
//	double M, h, f, alpha, beta;
//	double em, ep, t, dt, ut;
//	double gi, gi2, wgi, Vg, dVg;
//	double zi, wzi, r, r2, Vr, dVr;
//
//	VR = dVR = 0; h = 0.075; beta = 0.25; 
//	for (j=0; j<nQ1; j++){
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
//		for (i=nm; i<=np; i++){
//			t = i*h;
//			if (t==0){
//				ut = 1.0/(2.0+alpha+beta);
//				gi = M*ut;
//				wgi = 0.5*M*h*(1.0+(alpha-beta)*ut*ut)*sin(f*gi);
//			}else{
//				em = exp(-t); ep = beta*(1.0/em-1); em = alpha*(em-1); 
//				dt = exp(-2*t+em-ep); ut = 1.0/(1.0-dt);
//				gi = M*t*ut;
//				wgi = M*h*(1.0-(2.0+alpha+beta+em+ep)*t*(ut-1))*ut;
//				if (i>5)
//					wgi *= sin(((i&1)?-1:1)*f*gi*dt);
//				else
//					wgi *= sin(f*gi);
//			}
//			fea(gi, Vg, dVg);
//			if (dw==0){
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

/************************************************************************************/
// 3D Potential calculation (Vr, dVr) where dVr is the first derivative along r
void cPotentialCPU::Pot3Da(double &r, double &f, double &df){
	int i;
	ir = 1.0/r; r2 = r*r;

	f = df = 0.0;
	switch (PotPar){
		case 1:
			// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
			for (i=0; i<4; i++){
				cl = AtomTypes.cVr.cl[i]; cnl = AtomTypes.cVr.cnl[i];
				f += ft = cl*exp(-cnl*r2);
				df += -2.0*cnl*r*ft;
			}
			break;
		case 2:
			// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
			for (i=0; i<5; i++){
				cl = AtomTypes.cVr.cl[i]; cnl = AtomTypes.cVr.cnl[i];
				f += ft = cl*exp(-cnl*r2);
				df += -2.0*cnl*r*ft;
			}
			break;
		case 3:
			// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
			for (i=0; i<5; i++){
				cl = AtomTypes.cVr.cl[i]; cnl = AtomTypes.cVr.cnl[i];
				f += ft = cl*exp(-cnl*r2);
				df += -2.0*cnl*r*ft;
			}
			break;
		case 4:
			// 4: Kirkland parameterization - 3 Yukawa + 3 Gaussians - [0, 12]			
			for (i=0; i<3; i++){
				cl = AtomTypes.cVr.cl[i]; cnl = AtomTypes.cVr.cnl[i];
				f += ft = cl*exp(-cnl*r)*ir;
				df += -(cnl+ ir)*ft;
			}

			for (i=3; i<6; i++){
				cl = AtomTypes.cVr.cl[i]; cnl = AtomTypes.cVr.cnl[i];
				f += ft = cl*exp(-cnl*r2);
				df += -2.0*cnl*r*ft;
			}
			break;
		case 5:
			// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
			for (i=0; i<6; i++){
				cl = AtomTypes.cVr.cl[i]; cnl = AtomTypes.cVr.cnl[i];
				f += ft = cl*erfcCPU(cnl*r)*ir;
				df += (-2.0*cl*cnl*exp(-cnl*cnl*r2)/cPii2-ft)*ir;
			}
			break;
		case 6:
			// 6: Lobato parameterization - 5 Hydrogen fe - [0, 12]
			for (i=0; i<5; i++){
				cl = AtomTypes.cVr.cl[i]; cnl = AtomTypes.cVr.cnl[i];
				ft = cl*exp(-cnl*r)*ir; icnl = 1.0/cnl;
				f += ft*(2.0*icnl + r);
				df += -ft*(2.0*icnl*ir + 2.0 + cnl*r);
			}
			break;
	}
}

// 3D Potential calculation (Vr, dVr) where dVr is the first derivative along r
void cPotentialCPU::Pot3Dn(double &r, double &f, double &df){
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
//	for (i=nm; i<=np; i++){
//		t = i*h;
//		if (t==0){
//			ut = 1.0/(2.0+alpha+beta);
//			gi = M*ut;
//			wgi = 0.5*(1.0+(alpha-beta)*ut*ut)*sin(f*gi);
//		}else{
//			em = exp(-t); ep = beta*(1.0/em-1); em = alpha*(em-1); 
//			dt = exp(-2*t+em-ep); ut = 1.0/(1.0-dt);
//			gi = M*t*ut;
//			wgi = (1.0-(2.0+alpha+beta+em+ep)*t*(ut-1))*ut;
//			if (i>5)
//				wgi *= sin(((i&1)?-1:1)*f*gi*dt);
//			else
//				wgi *= sin(f*gi);
//		}
//		gi2 = gi*gi;
//		fea(gi, Vg, dVg);
//		if (dw==0){
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

/************************************************************************************/
void cPotentialCPU::Vr(int IntType, int Dim, double r, double &f, double &df){
	switch (Dim){
		case 1:
			if(~IntType) 
				Pot1Da(r, f, df);
			else
				Pot1Dn(r, f, df);
			break;
		case 2:
			if(~IntType) 
				Pot2Da(r, f, df);
			else
				Pot2Dn(r, f, df);
			break;
		case 3:
			if(~IntType) 
				Pot3Da(r, f, df);
			else
				Pot3Dn(r, f, df);
			break;
	}

}

/************************************************************************************/
// 2D Potential calculation (VR, dVR) where dVR is the first derivative along R
void cPotentialCPU::Vr(int IntType, int Dim, int nr, double *r, double *f, double *df){
	for (int i=0; i<nr; i++)
		Vr(IntType, Dim, r[i], f[i], df[i]);
}

// Calculate the total Potential
void cPotentialCPU::Vr(int PotPar, int nAtoms, sAtoms *Atoms, int nAtomTypes, sAtomTypesCPU *AtomTypes, int IntType, int Dim, int nr, double *r, double factor, double *f, double *df){
	int i, j;
	double *ft, *dft; 
	double occ;

	ft = new double[nr]; 
	dft = new double[nr];

	// Initiazed potential
	for (i=0; i<nr; i++){
		f[i] = df[i] = 0.0;
	};

	for (i=0; i<nAtomTypes; i++){
		occ = 0.0;
		for (j=0; j<nAtoms; j++)
			if(Atoms[j].Z==AtomTypes[i].Z)
				occ += Atoms[j].occ;

		if(occ>0){
			SetAtomTypes(PotPar, AtomTypes[i]);
			SetSigma(0.0);
			Vr(IntType, Dim, nr, r, ft, dft);
			// Add potential
			for (j=0; j<nr; j++){
				f[j] += occ*factor*ft[j];
				df[j] += occ*factor*dft[j];
			}
		}
	}

	delete [] ft;
	delete [] dft;
}

/************************************************************************************/
double cPotentialCPU::AtomicRadius_rms(int Dim){
	double iVr, iVrn, Vrt;
	double ri, wi, Vri, dVri;
	double ra, rmin;

	if (AtomTypes.cVr.cl[0]==0)
		ra = 0.0;
	else{
		rmin = AtomTypes.rn_c;
		Vr(0, Dim, rmin, Vri, dVri);

		iVr = Vri*pow(rmin, Dim)/Dim;
		iVrn = Vri*pow(rmin, Dim+2)/(Dim+2);
		/*************************************************************/
		for (int i=0; i<nQ1; i++){
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

double cPotentialCPU::AtomicRadius_Cutoff(int Dim, double Vrl){
	double eps = 1e-08;
	double rmin, rmax;
	double rm, Vrm, dVrm;

	if ((Vrl==0)||(AtomTypes.cVr.cl[0]==0))
		rmax = 0;
	else{
		rmin = AtomTypes.rn_c; 
		rmax = 50.0;
		Vr(0, Dim, rmin, Vrm, dVrm); 
		Vrm = std::abs(Vrm);
		while (std::abs(Vrm-Vrl)>eps){
			rm = 0.5*(rmin+rmax);
			Vr(0, Dim, rm, Vrm, dVrm);
			Vrm = std::abs(Vrm);
			if (Vrm>Vrl)
				rmin = rm;
			else
				rmax = rm;
		}
	}
	return rmax;
}