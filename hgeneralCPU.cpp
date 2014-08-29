#include "hmathCPU.h"
#include "hConstTypes.h"
#include "hgeneralCPU.h"
#include "hAtomicData.h"

// Input: E0(keV), Output: lambda (electron wave)
double fgetLambda(double E0){
	double emass = 510.99906;		// electron rest mass in keV
	double hc = 12.3984244;			// Planck's const x speed of light
		
	double lambda = hc/sqrt(E0*(2.0*emass + E0));
	return lambda;
}

// Input: E0(keV), Output: sigma (Interaction parameter)
double fgetSigma(double E0){
	double emass = 510.99906;		// electron rest mass in keV
	double hc = 12.3984244;			// Planck's const x speed of light
	double x = (emass + E0)/(2.0*emass + E0);
	
	double lambda = hc/sqrt(E0*(2.0*emass + E0));
	double sigma = 2.0*cPi*x/(lambda*E0);
	return sigma;
}

// Input: E0(keV), Output: gamma(relativistic factor)
double fgetGamma(double E0){
	double emass = 510.99906;		// electron rest mass in keV

	double gamma = (1.0 + E0/emass);
	return gamma;
}

// get index (with typ=0: bottom index for equal values and typ=1: upper index for equal values)
int fgetIndex(int ixmin, int ixmax, double *x, int typ, double x0){
	int ixmid;	
	if(typ==0){
		do{
			ixmid = (ixmin + ixmax)>>1;		// divide by 2
			if(x0 <= x[ixmid]) 
				ixmax = ixmid;
			else 
				ixmin = ixmid;
		}while ((ixmax-ixmin)>1);
	}else{
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

// get two dimensional radial distribution for regular grid
void fget2DRadDist(int nR, double *R, double *fR, int nRl, double *Rl, double *rl, double *frl, double *cfrl, bool reg, int typ){	
	int i, j;
	double Rlmin = Rl[0], Rlmax = Rl[nRl-1], dRl = Rl[1]-Rl[0];
  
	for (i=0; i<nRl-1; i++){
		rl[i] = 0.5*(Rl[i]+Rl[i+1]);
		frl[i] = 0.0;
		cfrl[i] = 0.0;
	}

	for (i=0; i<nR; i++)
		if ((Rlmin<=R[i])&&(R[i]<Rlmax)){
			j = (reg)?(int)floor((R[i]-Rlmin)/dRl):fgetIndex(0, nRl-1, Rl, 0, R[i]);
			frl[j] += fR[i];
			cfrl[j] += 1.0;
		}

	if (typ==0)
		for (i=0; i<nRl-1; i++)
			if (cfrl[i]>0)
				frl[i] /= cfrl[i];
}

// Get atom types
void fAtoms2AtomTypes(int nAtomsM, double *AtomsM, int PotPar, double Vrl, int &nAtomTypes, sAtomTypesCPU *&AtomTypes){
	int i, j;
	double Z[NE];
	//double sigma[NE];

	// initialize aZ array
	for (i=0; i<NE; i++){
		Z[i] = 0.0;
		//sigma[i] = 0;
	}

	nAtomTypes = 0;
	for (i=0; i<nAtomsM; i++){
		j = (int)AtomsM[3*nAtomsM + i] - 1; // Atom type
		if (Z[j]==0.0)
			nAtomTypes++;
		Z[j] = Z[j] + AtomsM[5*nAtomsM + i];
		//sigma[j] += Atomsi[4*na + i]; // displacement
	}

	// get atom types
	j = 0;
	AtomTypes = new sAtomTypesCPU[nAtomTypes];
	for (i=0; i<NE; i++){
		if (Z[i]>0){
			AtomTypes[j].Z = i+1;
			AtomTypes[j].occ = Z[i];
			AtomTypes[j].PotPar = PotPar;
			AtomTypes[j].ns = 0;
			//AtomTypesCPU[j].sigma = sigma[i]/double(Z[i]); This have to be modified
			j++;
		}
	}

	cAtomicData AtomicData;
	AtomicData.ReadAtomicData(nAtomTypes, AtomTypes, Vrl);
}

// Set atom types
void fSetAtomTypes(int Z, double occ, int PotPar, int ns, double Vrl, sAtomTypesCPU &AtomTypes){
	int nAtomTypes = 1;

	AtomTypes.Z = Z;
	AtomTypes.occ = occ;
	AtomTypes.PotPar = PotPar;
	AtomTypes.ns = ns;

	cAtomicData AtomicData;
	AtomicData.ReadAtomicData(nAtomTypes, &AtomTypes, Vrl);
}

// Set atom types
void fSetAtomTypes(int Zi, int Ze, double occ, int PotPar, int ns, double Vrl, int nAtomTypes, sAtomTypesCPU *AtomTypes){

	for (int i=Zi; i<=Ze; i++){
		AtomTypes[i-Zi].Z = i;
		AtomTypes[i-Zi].occ = occ;
		AtomTypes[i-Zi].PotPar = PotPar;
		AtomTypes[i-Zi].ns = ns;
	}

	cAtomicData AtomicData;
	AtomicData.ReadAtomicData(nAtomTypes, AtomTypes, Vrl);
}

// Grid's parameter initialization
void fsGP_Init(sGP &GP){
	GP.nx = 0;
	GP.ny = 0;
	GP.nxh = 0;
	GP.nyh = 0;
	GP.nxy = 0;
	GP.inxy = 0;
	GP.lx = 0;
	GP.ly = 0;
	GP.dz = 0;
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

// Grid's parameter initialization
void fsBT_Init(sBT &BT){
	BT.Bnxny.x = 0;		BT.Bnxny.y = 0;		BT.Bnxny.z = 0;
	BT.Tnxny.x = 0;		BT.Tnxny.y = 0;		BT.Tnxny.z = 0;
	BT.Bhnxny.x = 0;	BT.Bhnxny.y = 0;	BT.Bhnxny.z = 0;
	BT.Thnxny.x = 0;	BT.Thnxny.y = 0;	BT.Thnxny.z = 0;
	BT.Bonxny.x = 0;	BT.Bonxny.y = 0;	BT.Bonxny.z = 0;
	BT.Tonxny.x = 0;	BT.Tonxny.y = 0;	BT.Tonxny.z = 0;
	BT.Bmnxny.x = 0;	BT.Bmnxny.y = 0;	BT.Bmnxny.z = 0;
	BT.Tmnxny.x = 0;	BT.Tmnxny.y = 0;	BT.Tmnxny.z = 0;
}

// Block and Thread parameter initialization
void fsLens_Init(sLens &Lens){
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

// Grid's parameter calculation
void fsGP_Cal(int nx, int ny, double lx, double ly, double dz, sGP &GP){
	GP.nx = nx;
	GP.ny = ny;
	GP.nxh = nx/2;
	GP.nyh = ny/2;
	GP.nxy = GP.nx*GP.ny;
	GP.inxy = (GP.nxy==0)?(0.0):(1.0/double(GP.nxy));
	GP.lx = lx;
	GP.ly = ly;
	GP.dz = dz;
	GP.dRx = (GP.nx==0)?(0.0):(GP.lx/double(GP.nx));
	GP.dRy = (GP.ny==0)?(0.0):(GP.ly/double(GP.ny));
	GP.dRmin = MIN(GP.dRx, GP.dRy);
	GP.dgx = (lx==0)?(0.0):(1.0/GP.lx);
	GP.dgy = (ly==0)?(0.0):(1.0/GP.ly);
	GP.dgmin = MIN(GP.dgx, GP.dgy);
	GP.gmax = MIN(double(GP.nxh)*GP.dgx, double(GP.nyh)*GP.dgy);
	GP.gmax2 = pow(GP.gmax, 2);
	GP.gmaxl = 2.0*GP.gmax/3.0;
	GP.gmaxl2 = pow(GP.gmaxl, 2);
}

// Block and Thread parameter calculation
void fsBT_Cal(sGP &GP, sBT &BT){
	int nx = GP.nx, ny = GP.ny;
	int nxh = GP.nxh, nyh = GP.nyh;
	double nxy = GP.nxy, mnxny = MAX(nx, ny);

	BT.Bnxny.x = (ny+thrnxny-1)/thrnxny; BT.Bnxny.y = (nx+thrnxny-1)/thrnxny; BT.Bnxny.z = 1;
	BT.Tnxny.x = thrnxny; BT.Tnxny.y = thrnxny; BT.Tnxny.z = 1;

	BT.Bhnxny.x = (nyh+thrnxny-1)/thrnxny; BT.Bhnxny.y = (nxh+thrnxny-1)/thrnxny; BT.Bhnxny.z = 1;
	BT.Thnxny.x = thrnxny; BT.Thnxny.y = thrnxny; BT.Thnxny.z = 1;

	BT.Bonxny.x = 1; BT.Bonxny.y = 1; BT.Bonxny.z = 1;
	BT.Tonxny.x = thrnxny; BT.Tonxny.y = thrnxny; BT.Tonxny.z = 1;

	BT.Bmnxny.x = (mnxny+thrmnxny-1)/thrmnxny; BT.Bmnxny.y = 1; BT.Bmnxny.z = 1;
	BT.Tmnxny.x = thrmnxny; BT.Tmnxny.y = 1; BT.Tmnxny.z = 1;

	BT.Bnxy.x = MIN(64, (nxy+thrnxy-1)/thrnxy); BT.Bnxy.y = 1; BT.Bnxy.z = 1;
	BT.Tnxy.x = thrnxy; BT.Tnxy.y = 1; BT.Tnxy.z = 1;
}

// Lens' parameter calculation
void fsLens_Cal(double E0, sGP &GP, sLens &Lens){
	Lens.gamma = fgetGamma(E0);

	Lens.lambda = fgetLambda(E0);
	Lens.lambda2 = pow(Lens.lambda, 2);

	Lens.cf = cPi*Lens.f*Lens.lambda;
	Lens.cCs3 = -cPi*Lens.Cs3*pow(Lens.lambda, 3)/2.0;
	Lens.cCs5 = -cPi*Lens.Cs5*pow(Lens.lambda, 5)/3.0;
	Lens.cmfa2 = -cPi*Lens.mfa2*Lens.lambda;
	Lens.cmfa3 = -2.0*Lens.mfa3*cPi*pow(Lens.lambda, 2)/3.0;
	Lens.gmin2 = pow(Lens.aobjl/Lens.lambda, 2);
	Lens.gmax2 = pow(Lens.aobju/Lens.lambda, 2);

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
};