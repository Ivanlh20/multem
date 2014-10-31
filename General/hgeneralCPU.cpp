#include <cmath>
#include <cstring>
#include "hmathCPU.h"
#include "hConstTypes.h"
#include "hgeneralCPU.h"
#include "hAtomicData.h"

// Input: E0(keV), Output: lambda (electron wave)
double f_getLambda(double E0){
	double emass = 510.99906;		// electron rest mass in keV
	double hc = 12.3984244;			// Planck's const x speed of light
		
	double lambda = hc/sqrt(E0*(2.0*emass + E0));
	return lambda;
}

// Input: E0(keV), Output: sigma (Interaction parameter)
double f_getSigma(double E0){
	double emass = 510.99906;		// electron rest mass in keV
	double hc = 12.3984244;			// Planck's const x speed of light
	double x = (emass + E0)/(2.0*emass + E0);
	
	double lambda = hc/sqrt(E0*(2.0*emass + E0));
	double sigma = 2.0*cPi*x/(lambda*E0);
	return sigma;
}

// Input: E0(keV), Output: gamma(relativistic factor)
double f_get2DRadDist(double E0){
	double emass = 510.99906;		// electron rest mass in keV

	double gamma = (1.0 + E0/emass);
	return gamma;
}

// get index (with typ=0: bottom index for equal values and typ=1: upper index for equal values)
int f_getIndex(int ixmin, int ixmax, double *x, int typ, double x0){
	int ixmid;	
	switch (typ){
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

// get two dimensional radial distribution for regular grid
void f_get2DRadDist(int nR, double *R, double *fR, int nRl, double *Rl, double *rl, double *frl, double *cfrl, bool reg, int typ){	
	int i, j;
	double Rlmin = Rl[0], Rlmax = Rl[nRl-1], dRl = Rl[1]-Rl[0];
 
	for (i=0; i<nRl-1; i++){
		rl[i] = 0.5*(Rl[i]+Rl[i+1]);
		frl[i] = 0.0;
		cfrl[i] = 0.0;
	}

	for (i=0; i<nR; i++)
		if ((Rlmin<=R[i])&&(R[i]<Rlmax)){
			j = (reg)?(int)floor((R[i]-Rlmin)/dRl):f_getIndex(0, nRl-1, Rl, 0, R[i]);
			frl[j] += fR[i];
			cfrl[j] += 1.0;
		}

	if (typ==0)
		for (i=0; i<nRl-1; i++)
			if (cfrl[i]>0)
				frl[i] /= cfrl[i];
}

// Set atom types
void f_SetAtomTypes(int Z, int PotPar, int ns, double Vrl, sAtomTypesCPU &AtomTypes){
	int nAtomTypes = 1;

	AtomTypes.Z = Z;
	AtomTypes.ns = ns;

	cAtomicData AtomicData;
	AtomicData.ReadAtomicData(PotPar, nAtomTypes, &AtomTypes, Vrl);
}

// Set atom types
void f_SetAtomTypes(int PotPar, int ns, double Vrl, int nAtomTypes, sAtomTypesCPU *&AtomTypes){
	for (int iAtomTypes=0; iAtomTypes<nAtomTypes; iAtomTypes++){
		AtomTypes[iAtomTypes].Z = iAtomTypes+1;
		AtomTypes[iAtomTypes].ns = ns;
		AtomTypes[iAtomTypes].Vo = 0; //Implementation in the future
	}
	cAtomicData AtomicData;
	AtomicData.ReadAtomicData(PotPar, nAtomTypes, AtomTypes, Vrl);
}

// Set Atoms
void f_AtomsM2Atoms(int nAtomsM_i, double *AtomsM_i, bool PBC_xyi, double lxi, double lyi, int &nAtoms, sAtoms *&Atoms){
	double x, y, dl = 1e-03;
	double lxb = lxi-dl, lyb = lyi-dl;

	Atoms = new sAtoms[nAtomsM_i];
	nAtoms = 0;
	for (int i=0; i<nAtomsM_i; i++){		
		x = AtomsM_i[0*nAtomsM_i + i];				// x-position
		y = AtomsM_i[1*nAtomsM_i + i];				// y-position
		if(!PBC_xyi||((x<lxb)&&(y<lyb))){
			Atoms[nAtoms].x = x;					// x-position
			Atoms[nAtoms].y = y;					// y-position
			Atoms[nAtoms].z = AtomsM_i[2*nAtomsM_i + i];		// z-position
			Atoms[nAtoms].Z = (int)AtomsM_i[3*nAtomsM_i + i];		// Atomic number
			Atoms[nAtoms].sigma = AtomsM_i[4*nAtomsM_i + i];		// Standard deviation
			Atoms[nAtoms].occ = AtomsM_i[5*nAtomsM_i + i];			// Occupancy
			nAtoms++;
		}
	}
}

// get 2D maximum interaction distance
double f_getRMax(int nAtoms, sAtoms *&Atoms, sAtomTypesCPU *&AtomTypes){
	double R, Rmax=0;

	for(int iAtoms=0; iAtoms<nAtoms; iAtoms++){
		R = AtomTypes[Atoms[iAtoms].Z-1].Rmax;
		if (Rmax < R)
			Rmax = R;
	}
	return Rmax;
}

/*************************************************************************/
// Multislice's general parameters initialization
void f_sMPG_Init(sMGP &MGP){
	MGP.gpu = 0;
	MGP.SimType = 10;
	MGP.MulOrder = 2;
	MGP.nConfFP = 0;
	MGP.DimFP = 111;
	MGP.DistFP = 1;
	MGP.SeedFP = 1983;
	MGP.PotPar = 6;
	MGP.MEffect = 1;
	MGP.STEffect = 1;
	MGP.ZeroDefTyp = 3;
	MGP.ZeroDefPlane = 0.0;
	MGP.ApproxModel = 1;
	MGP.PBC_xy = true;
	MGP.Vrl = stVrl;
	MGP.E0 = 300;
	MGP.theta = 0;	
	MGP.phi = 0;
	MGP.lx = 4.078;
	MGP.ly = 4.078;	
	MGP.dz = 0.25;
	MGP.nx = 1024;
	MGP.ny = 1024;
	if(MGP.ApproxModel>1) MGP.MulOrder = 1;
}

/*************************************************************************/
// Grid's parameter initialization
void f_sGP_Init(sGP &GP){
	GP.nx = 0;
	GP.ny = 0;
	GP.nxh = 0;
	GP.nyh = 0;
	GP.nxy = 0;
	GP.inxy = 0;
	GP.lx = 0;
	GP.ly = 0;
	GP.dz = 0;
	GP.PBC_xy = true;
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

// Grid's parameter calculation
void f_sGP_Cal(int nx, int ny, double lx, double ly, double dz, bool PBC_xy, sGP &GP){
	GP.nx = nx;
	GP.ny = ny;
	GP.nxh = nx/2;
	GP.nyh = ny/2;
	GP.nxy = GP.nx*GP.ny; 
	GP.inxy = (GP.nxy==0)?(0.0):(1.0/double(GP.nxy));
	GP.lx = lx;
	GP.ly = ly;
	GP.dz = dz;
	GP.PBC_xy = PBC_xy;
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

/*************************************************************************/
// Grid's parameter initialization
void f_sBT_Init(sBT &BT){
	BT.Bnxny.x = 0;		BT.Bnxny.y = 0;		BT.Bnxny.z = 0;
	BT.Tnxny.x = 0;		BT.Tnxny.y = 0;		BT.Tnxny.z = 0;
	BT.Bhnxny.x = 0;	BT.Bhnxny.y = 0;	BT.Bhnxny.z = 0;
	BT.Thnxny.x = 0;	BT.Thnxny.y = 0;	BT.Thnxny.z = 0;
	BT.Bmnxny.x = 0;	BT.Bmnxny.y = 0;	BT.Bmnxny.z = 0;
	BT.Tmnxny.x = 0;	BT.Tmnxny.y = 0;	BT.Tmnxny.z = 0;
}

// Block and Thread parameter calculation
void f_sBT_Cal(sGP &GP, sBT &BT){
	int nx = GP.nx, ny = GP.ny;
	int nxh = GP.nxh, nyh = GP.nyh;
	double nxy = GP.nxy, mnxny = MAX(nx, ny);

	BT.Bnxny.x = (ny+thrnxny-1)/thrnxny; BT.Bnxny.y = (nx+thrnxny-1)/thrnxny; BT.Bnxny.z = 1;
	BT.Tnxny.x = thrnxny; BT.Tnxny.y = thrnxny; BT.Tnxny.z = 1;

	BT.Bhnxny.x = (nyh+thrnxny-1)/thrnxny; BT.Bhnxny.y = (nxh+thrnxny-1)/thrnxny; BT.Bhnxny.z = 1;
	BT.Thnxny.x = thrnxny; BT.Thnxny.y = thrnxny; BT.Thnxny.z = 1;

	BT.Bmnxny.x = (mnxny+thrmnxny-1)/thrmnxny; BT.Bmnxny.y = 1; BT.Bmnxny.z = 1;
	BT.Tmnxny.x = thrmnxny; BT.Tmnxny.y = 1; BT.Tmnxny.z = 1;

	BT.Bnxy.x = MIN(64, (nxy+thrnxy-1)/thrnxy); BT.Bnxy.y = 1; BT.Bnxy.z = 1;
	BT.Tnxy.x = thrnxy; BT.Tnxy.y = 1; BT.Tnxy.z = 1;
}

/*************************************************************************/
// Block and Thread parameter initialization
void f_sLens_Init(sLens &Lens){
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
void f_sLens_Cal(double E0, sGP &GP, sLens &Lens){
	Lens.gamma = f_get2DRadDist(E0);

	Lens.lambda = f_getLambda(E0);
	Lens.lambda2 = pow(Lens.lambda, 2);

	Lens.cf = cPi*Lens.f*Lens.lambda;
	Lens.cCs3 = -cPi*Lens.Cs3*pow(Lens.lambda, 3)/2.0;
	Lens.cCs5 = -cPi*Lens.Cs5*pow(Lens.lambda, 5)/3.0;
	Lens.cmfa2 = -cPi*Lens.mfa2*Lens.lambda;
	Lens.cmfa3 = -2.0*cPi*Lens.mfa3*pow(Lens.lambda, 2)/3.0;
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

/*************************************************************************/
void f_sCoefPar_Free(sCoefPar &CoefPar){
	delete [] CoefPar.cl; CoefPar.cl = 0;
	delete [] CoefPar.cnl; CoefPar.cnl = 0;
}

void f_sCoefPar_Init(sCoefPar &CoefPar){
	CoefPar.cl = 0;
	CoefPar.cnl = 0;
}

void f_sCoefPar_Malloc(int nCoefPar, sCoefPar &CoefPar){
	if(nCoefPar<=0) return;

	CoefPar.cl = new double[nCoefPar];
	CoefPar.cnl = new double[nCoefPar];
}

/*************************************************************************/
void f_sciVn_Free(sciVn &ciVn){
	delete [] ciVn.c0; ciVn.c0 = 0;
	delete [] ciVn.c1; ciVn.c1 = 0;
	delete [] ciVn.c2; ciVn.c2 = 0;
	delete [] ciVn.c3; ciVn.c3 = 0;
}

void f_sciVn_Init(sciVn &ciVn){
	ciVn.c0 = 0;
	ciVn.c1 = 0;
	ciVn.c2 = 0;
	ciVn.c3 = 0;
}

void f_sciVn_Malloc(int nciVn, sciVn &ciVn){
	if(nciVn<=0) return;

	ciVn.c0 = new double[nciVn];
	ciVn.c1 = new double[nciVn];
	ciVn.c2 = new double[nciVn];
	ciVn.c3 = new double[nciVn];
}

/*************************************************************************/
void f_sDetCir_Free(sDetCir &DetCir){
	delete [] DetCir.g2min; DetCir.g2min = 0;
	delete [] DetCir.g2max; DetCir.g2max = 0;
}

void f_sDetCir_Init(sDetCir &DetCir){
	DetCir.g2min = 0;
	DetCir.g2max = 0;
}

void f_sDetCir_Malloc(int nDetCir, sDetCir &DetCir){
	if(nDetCir<=0) return;

	DetCir.g2min = new double[nDetCir];
	DetCir.g2max = new double[nDetCir];
}

/*************************************************************************/
void f_BuildGrid(int line, int ns, double x0, double y0, double xe, double ye, int &nxs, int &nys, double *&xs, double *&ys){
	if(ns<=0){
		nxs = nys = 0;
		xs = 0; ys = 0;
		return;
	}

	int i;
	double lxs = xe-x0, lys = ye-y0, ld = sqrt(lxs*lxs+lys*lys);
	double theta = atan(lys/lxs), costheta = cos(theta), sintheta = sin(theta);

	if (line){
			nxs = ns;
			nys = ns;
	}else{
		nxs = (std::abs(lxs)>std::abs(lys))?ns:(int)ceil(ns*std::abs(lxs/lys));
		nys = (std::abs(lxs)>std::abs(lys))?(int)ceil(ns*std::abs(lys/lxs)):ns;
	}

	xs = new double [nxs];
	ys = new double [nys];

	double ds;

	if (line){
		ds = ld/ns;
		for (i=0; i<ns; i++){
			xs[i] = x0 + i*ds*costheta;
			ys[i] = y0 + i*ds*sintheta;
		}
	}else{
		ds = lxs/nxs;
		for (i=0; i<nxs; i++)
			xs[i] = x0 + i*ds;

		ds = lys/nys;
		for (i=0; i<nys; i++)
			ys[i] = y0 + i*ds;
	}

}

/*************************************************************************/
void f_InMulSli_Init(sInMSTEM &InMSTEM){
	InMSTEM.gpu = 0;
	InMSTEM.SimType = 0;
	InMSTEM.MulOrder = 0;
	InMSTEM.nConfFP = 0;
	InMSTEM.DimFP = 0;
	InMSTEM.SeedFP = 0;
	InMSTEM.PotPar = 0;
	InMSTEM.MEffect = 0;		
	InMSTEM.STEffect = 0;	
	InMSTEM.ZeroDefTyp = 0;
	InMSTEM.ZeroDefPlane = 0;
	InMSTEM.ApproxModel = 0;

	InMSTEM.E0 = 0;
	InMSTEM.theta = 0;
	InMSTEM.phi = 0;
	InMSTEM.nx = 0;
	InMSTEM.ny = 0;
	InMSTEM.lx = 0;
	InMSTEM.ly = 0;
	InMSTEM.dz = 0;

	InMSTEM.MC_m = 0;
	InMSTEM.MC_f = 0;
	InMSTEM.MC_Cs3 = 0;
	InMSTEM.MC_Cs5 = 0;
	InMSTEM.MC_mfa2 = 0;
	InMSTEM.MC_afa2 = 0;
	InMSTEM.MC_mfa3 = 0;
	InMSTEM.MC_afa3 = 0;
	InMSTEM.MC_aobjl = 0;
	InMSTEM.MC_aobju = 0;
	InMSTEM.MC_sf = 0;
	InMSTEM.MC_nsf = 0;
	InMSTEM.MC_beta = 0;
	InMSTEM.MC_nbeta = 0;

	InMSTEM.nAtomsM = 0;
	InMSTEM.AtomsM = 0;

	InMSTEM.STEM_line = 0;
	InMSTEM.STEM_FastCal = false;
	InMSTEM.STEM_ns = 0;
	InMSTEM.STEM_x1u = 0;
	InMSTEM.STEM_y1u = 0;
	InMSTEM.STEM_x2u = 0;
	InMSTEM.STEM_y2u = 0;
	
	InMSTEM.STEM_nDet = 0;
	InMSTEM.STEM_DetCir = 0;

	InMSTEM.CBED_x0 = 0;
	InMSTEM.CBED_y0 = 0;

	//InMSTEM.HRTEM_xx;

	InMSTEM.PED_nrot = 0;
	InMSTEM.PED_theta = 0;

	InMSTEM.HCI_nrot = 0;
	InMSTEM.HCI_theta = 0;
	//InMSTEM.HCI_xx;

	//InMSTEM.EWRS_xx;

	//InMSTEM.EWFS_xx;
}

void f_InMulSli_Free(sInMSTEM &InMSTEM){
	InMSTEM.gpu = 0;
	InMSTEM.SimType = 0;
	InMSTEM.MulOrder = 0;
	InMSTEM.nConfFP = 0;
	InMSTEM.DimFP = 0;
	InMSTEM.SeedFP = 0;
	InMSTEM.PotPar = 0;
	InMSTEM.MEffect = 0;		
	InMSTEM.STEffect = 0;	
	InMSTEM.ZeroDefTyp = 0;
	InMSTEM.ZeroDefPlane = 0;
	InMSTEM.ApproxModel = 0;
	
	InMSTEM.E0 = 0;
	InMSTEM.theta = 0;
	InMSTEM.phi = 0;
	InMSTEM.nx = 0;
	InMSTEM.ny = 0;
	InMSTEM.lx = 0;
	InMSTEM.ly = 0;
	InMSTEM.dz = 0;

	InMSTEM.MC_m = 0;
	InMSTEM.MC_f = 0;
	InMSTEM.MC_Cs3 = 0;
	InMSTEM.MC_Cs5 = 0;
	InMSTEM.MC_mfa2 = 0;
	InMSTEM.MC_afa2 = 0;
	InMSTEM.MC_mfa3 = 0;
	InMSTEM.MC_afa3 = 0;
	InMSTEM.MC_aobjl = 0;
	InMSTEM.MC_aobju = 0;
	InMSTEM.MC_sf = 0;
	InMSTEM.MC_nsf = 0;
	InMSTEM.MC_beta = 0;
	InMSTEM.MC_nbeta = 0;

	InMSTEM.nAtomsM = 0;
	InMSTEM.AtomsM = 0;

	InMSTEM.STEM_line = 0;
	InMSTEM.STEM_FastCal = false;
	InMSTEM.STEM_ns = 0;
	InMSTEM.STEM_x1u = 0;
	InMSTEM.STEM_y1u = 0;
	InMSTEM.STEM_x2u = 0;
	InMSTEM.STEM_y2u = 0;
	
	InMSTEM.STEM_nDet = 0;
	delete [] InMSTEM.STEM_DetCir; InMSTEM.STEM_DetCir = 0;

	InMSTEM.CBED_x0 = 0;
	InMSTEM.CBED_y0 = 0;

	//InMSTEM.HRTEM_xx;

	InMSTEM.PED_nrot = 0;
	InMSTEM.PED_theta = 0;

	InMSTEM.HCI_nrot = 0;
	InMSTEM.HCI_theta = 0;
	//InMSTEM.HCI_xx;

	//InMSTEM.EWRS_xx;

	//InMSTEM.EWFS_xx;
}

void f_InMulSli_Check(sInMSTEM &InMSTEMi, sInMSTEM &InMSTEMo){
	std::memcpy(&InMSTEMo, &InMSTEMi, sizeof(sInMSTEM));

	//if(InMSTEMo.gpu<0)			// gpu device 
	//	InMSTEMo.gpu = 0;

	//if(InMSTEMo.SimType<0)		// 1: STEM, 2: HRTEM, 3: ED, 4: PED, 5: CBED, 6: HCI, ... 10: EW real, 11: EW Fourier
	//	InMSTEMo.SimType = 10;

	//if((InMSTEMo.MulOrder<0)||(InMSTEMo.MulOrder>2)) // 1: First order MS, 2: Second Order MS
	//	InMSTEMo.MulOrder = 1;

	//if(InMSTEMo.nConfFP<0) // Number of frozen phonon configurations
	//	InMSTEMo.nConfFP = 0;

	//if((InMSTEMo.DimFP<1)||(InMSTEMo.DimFP>3)) // Dimensions phonon configurations
	//	InMSTEMo.DimFP = 3;

	//if(InMSTEMo.SeedFP<=0) // Random seed(frozen phonon)
	//	InMSTEMo.SeedFP = 1983;

	//if((InMSTEMo.MEffect<0)||(InMSTEMo.MEffect>2)) // 1: Exit wave Partial coherente mode, 2: Transmission cross coefficient
	//	InMSTEMo.MEffect = 0;

	//if((InMSTEMo.STEffect<0)||(InMSTEMo.STEffect>2)) // 1: Spatial and temporal, 2: Temporal, 3: Spatial
	//	InMSTEMo.STEffect = 0;

	//if((InMSTEMo.ZeroDefTyp<0)||(InMSTEMo.ZeroDefTyp>3)) // 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
	//	InMSTEMo.ZeroDefTyp = 1;

	//if((InMSTEMo.PotPar<1)||(InMSTEMo.PotPar>6)) // Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
	//	InMSTEMo.PotPar = 6;

	///***************************Multislice*************************/
	//if(InMSTEMo.E0<=0) // Acceleration voltage
	//	InMSTEMo.E0 = 300;

	//if((InMSTEMo.theta<0)||(InMSTEMo.theta>cPi)) // incident tilt (in spherical coordinates) (degrees-->rad)
	//	InMSTEMo.theta = 0;

	//if((InMSTEMo.phi<-c2Pi)||(InMSTEMo.phi>c2Pi)) // incident tilt (in spherical coordinates) (degrees-->rad)
	//	InMSTEMo.phi = 0;

	//if(InMSTEMo.nx<=0) // Number of pixels in x direction
	//	InMSTEMo.nx = 128;

	//if(InMSTEMo.ny<=0) // Number of pixels in y direction
	//	InMSTEMo.ny = 128;

	//if(InMSTEMo.lx<0) // distance in x direction(Angstroms)
	//	InMSTEMo.lx = 0;

	//if(InMSTEMo.ly<0) // distance in y direction(Angstroms)
	//	InMSTEMo.ly = 0;

	//if(InMSTEMo.dz<0) // slice thickness
	//	InMSTEMo.dz = 0;

	///*********************Microscope parameters********************/
	//mxArray *mxMC= mxGetField(mxInMSTEM, 0, "MC");
	//InMSTEMo.MC_m =  ReadValuemxField<int>(mxMC, 0, "m");						// momentum of the vortex
	//InMSTEMo.MC_f = ReadValuemxField<double>(mxMC, 0, "f");					// defocus(Angstrom)
	//InMSTEMo.MC_Cs3 = ReadValuemxField<double>(mxMC, 0, "Cs3", mm2Ags);		// spherical aberration(mm-->Angstrom)
	//InMSTEMo.MC_Cs5 = ReadValuemxField<double>(mxMC, 0, "Cs5", mm2Ags);		// spherical aberration(mm-->Angstrom)
	//InMSTEMo.MC_mfa2 = ReadValuemxField<double>(mxMC, 0, "mfa2");				// magnitude 2-fold astigmatism(Angstrom)
	//InMSTEMo.MC_afa2 = ReadValuemxField<double>(mxMC, 0, "afa2", deg2rad);		// angle 2-fold astigmatism(degrees-->rad)
	//InMSTEMo.MC_mfa3 = ReadValuemxField<double>(mxMC, 0, "mfa3");				// magnitude 3-fold astigmatism(Angstrom)
	//InMSTEMo.MC_afa3 = ReadValuemxField<double>(mxMC, 0, "afa3", deg2rad);		// angle 3-fold astigmatism(degrees-->rad)
	//InMSTEMo.MC_aobjl = ReadValuemxField<double>(mxMC, 0, "aobjl", mrad2rad);	// lower objective aperture(mrad-->rad)
	//InMSTEMo.MC_aobju = ReadValuemxField<double>(mxMC, 0, "aobju", mrad2rad);	// upper objective aperture(mrad-->rad)
	//InMSTEMo.MC_sf = ReadValuemxField<double>(mxMC, 0, "sf");					// defocus spread(Angstrom)
	//InMSTEMo.MC_nsf = ReadValuemxField<int>(mxMC, 0, "nsf");						// Number of defocus sampling point
	//InMSTEMo.MC_beta = ReadValuemxField<double>(mxMC, 0, "beta", mrad2rad);	// semi-convergence angle(mrad-->rad)
	//InMSTEMo.MC_nbeta = ReadValuemxField<int>(mxMC, 0, "nbeta");					// half number sampling points

	//mxArray *mxAtomsM = mxGetField(mxInMSTEM, 0, "Atoms");
	//InMSTEMo.nAtomsM = (int)mxGetM(mxAtomsM);	// Number of Atoms
	//InMSTEMo.AtomsM = mxGetPr(mxAtomsM);		// Atoms in a matrix form

	//switch (InMSTEMo.SimType){	
	//	case 0:	// STEM
	//		mxArray *mxSTEM;	
	//		mxSTEM = mxGetField(mxInMSTEM, 0, "STEM");
	//		InMSTEMo.STEM_line = ReadValuemxField<int>(mxSTEM, 0, "line");
	//		InMSTEMo.STEM_GridType = ReadValuemxField<int>(mxSTEM, 0, "GridType");
	//		InMSTEMo.STEM_ns = ReadValuemxField<int>(mxSTEM, 0, "ns");
	//		InMSTEMo.STEM_nucx = ReadValuemxField<int>(mxSTEM, 0, "nucx");
	//		InMSTEMo.STEM_nucy = ReadValuemxField<int>(mxSTEM, 0, "nucy");
	//		InMSTEMo.STEM_x1u = ReadValuemxField<double>(mxSTEM, 0, "x1u");
	//		InMSTEMo.STEM_y1u = ReadValuemxField<double>(mxSTEM, 0, "y1u");
	//		InMSTEMo.STEM_x2u = ReadValuemxField<double>(mxSTEM, 0, "x2u");
	//		InMSTEMo.STEM_y2u = ReadValuemxField<double>(mxSTEM, 0, "y2u");

	//		InMSTEMo.STEM_nDet = ReadValuemxField<int>(mxSTEM, 0, "nDet");
	//		InMSTEMo.STEM_DetCir = new sInDetCir[InMSTEMo.STEM_nDet];
	//		mxArray *mxDetCir;
	//		mxDetCir = mxGetField(mxSTEM, 0, "DetCir");
	//		for (int i=0; i<InMSTEMo.STEM_nDet; i++){
	//			InMSTEMo.STEM_DetCir[i].InnerAng = ReadValuemxField<double>(mxDetCir, i, "InnerAng", mrad2rad);	// Inner angle(mrad-->rad)
	//			InMSTEMo.STEM_DetCir[i].OuterAng = ReadValuemxField<double>(mxDetCir, i, "OuterAng", mrad2rad);	// Outer angle(mrad-->rad)
	//		}
	//		break;
	//	case 1:		// HRTEM
	//		mxArray *mxHRTEM;
	//		mxHRTEM = mxGetField(mxInMSTEM, 0, "HRTEM");

	//		break;
	//	case 2:		// ED
	//		mxArray *mxED;
	//		mxED = mxGetField(mxInMSTEM, 0, "ED");

	//		break;
	//	case 3:		// PED
	//		mxArray *mxPED;
	//		mxPED = mxGetField(mxInMSTEM, 0, "PED");
	//		InMSTEMo.PED_nrot = ReadValuemxField<int>(mxPED, 0, "nrot");				// Number of orientations
	//		InMSTEMo.PED_theta = ReadValuemxField<double>(mxPED, 0, "theta", deg2rad);	// Precession angle(degrees-->rad)
	//		break;
	//	case 4:		// CBED
	//		mxArray *mxCBED;
	//		mxCBED = mxGetField(mxInMSTEM, 0, "CBED");
	//		InMSTEMo.CBED_x0 = ReadValuemxField<double>(mxCBED, 0, "x0");	// Precession angle(degrees-->rad)
	//		InMSTEMo.CBED_y0 = ReadValuemxField<double>(mxCBED, 0, "y0");	// Precession angle(degrees-->rad)
	//		break;
	//	case 5:		// HCI
	//		mxArray *mxHCI;
	//		mxHCI = mxGetField(mxInMSTEM, 0, "HCI");
	//		InMSTEMo.HCI_nrot = ReadValuemxField<int>(mxHCI, 0, "nrot");				// Number of orientations
	//		InMSTEMo.HCI_theta = ReadValuemxField<double>(mxHCI, 0, "theta", deg2rad);	// Precession angle(degrees-->rad)
	//		break;
	//	case 10:	// EW real
	//		mxArray *mxEWr;
	//		mxEWr = mxGetField(mxInMSTEM, 0, "EWr");

	//		break;
	//	case 11:	// EW Fourier
	//		mxArray *mxEWF;
	//		mxEWF = mxGetField(mxInMSTEM, 0, "EWF");

	//		break;
	//}
 }
