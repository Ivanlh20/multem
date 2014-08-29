#include "hConstTypes.h"
#include "hgeneralGPU.h"
#include "hSpecimenCPU.h"
#include "hAtomTypesGPU.h"
#include "hPotentialGPU.h"
#include "hQuadrature.h"
#include "cuda.h"
#include "math.h"
#include "cuda_runtime.h"

__constant__ int PotPar = 0;
__constant__ scVp cVp[stncVp];

// Potential Evaluation (Vr, dVrir)
inline __device__ void Pot3Da(double &r, volatile double *cl, volatile double *cnl, volatile double *icnl, double f, double &Vr, double &dVrir){
	double ir, r2;
	double Vr0, Vr1, Vr2, Vr3, Vr4, Vr5;

	switch (PotPar){
		case 1:
			// 1: doyle and turner parameterization - 4 gaussians - [0, 4]
			r2 = r*r;
			Vr0 = cl[0]*exp(-cnl[0]*r2); Vr1 = cl[1]*exp(-cnl[1]*r2); Vr2 = cl[2]*exp(-cnl[2]*r2); 
			Vr3 = cl[3]*exp(-cnl[3]*r2);
			Vr = Vr0 + Vr1 + Vr2 + Vr3;
			dVrir = 2.0*(cnl[0]*Vr0 + cnl[1]*Vr1 + cnl[2]*Vr2 + cnl[3]*Vr3);
			Vr *= f;
			dVrir *= -f;
			break;
		case 2:
			// 2: peng et al. parameterization - 5 gaussians - [0, 4]
			r2 = r*r;
			Vr0 = cl[0]*exp(-cnl[0]*r2); Vr1 = cl[1]*exp(-cnl[1]*r2); Vr2 = cl[2]*exp(-cnl[2]*r2); 
			Vr3 = cl[3]*exp(-cnl[3]*r2); Vr4 = cl[4]*exp(-cnl[4]*r2);
			Vr = Vr0 + Vr1 + Vr2 + Vr3 + Vr4;
			dVrir = 2.0*(cnl[0]*Vr0 + cnl[1]*Vr1 + cnl[2]*Vr2 + cnl[3]*Vr3 + cnl[4]*Vr4);
			Vr *= f;
			dVrir *= -f;
			break;
		case 3:		
			// 3: peng et al. parameterization - 5 gaussians - [0, 12]
			r2 = r*r;
			Vr0 = cl[0]*exp(-cnl[0]*r2); Vr1 = cl[1]*exp(-cnl[1]*r2); Vr2 = cl[2]*exp(-cnl[2]*r2); 
			Vr3 = cl[3]*exp(-cnl[3]*r2); Vr4 = cl[4]*exp(-cnl[4]*r2);
			Vr = Vr0 + Vr1 + Vr2 + Vr3 + Vr4;
			dVrir = 2.0*(cnl[0]*Vr0 + cnl[1]*Vr1 + cnl[2]*Vr2 + cnl[3]*Vr3 + cnl[4]*Vr4);
			Vr *= f;
			dVrir *= -f;
			break;
		case 4:
			// 4: kirkland parameterization - 3 yukawa + 3 gaussians - [0, 12]
			ir = 1.0/r; r2 = r*r;
			Vr0 = cl[0]*exp(-cnl[0]*r)*ir; Vr1 = cl[1]*exp(-cnl[1]*r)*ir; Vr2 = cl[2]*exp(-cnl[2]*r)*ir; 
			Vr3 = cl[3]*exp(-cnl[3]*r2); Vr4 = cl[4]*exp(-cnl[4]*r2); Vr5 = cl[5]*exp(-cnl[5]*r2);
			Vr = Vr0 + Vr1 + Vr2 + Vr3 + Vr4 + Vr5;
			dVrir = Vr0*(cnl[0]+ir) + Vr1*(cnl[1]+ir) + Vr2*(cnl[2]+ir) + 2.0*r*(cnl[3]*Vr3 + cnl[4]*Vr4 + cnl[5]*Vr5);
			Vr *= f;
			dVrir *= -f/r;
			break;
		case 5:
			// 5: weickenmeier and h.kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
			r2 = r*r;
			Vr0 = cl[0]*erfc(cnl[0]*r); Vr1 = cl[1]*erfc(cnl[1]*r); Vr2 = cl[2]*erfc(cnl[2]*r); 
			Vr3 = cl[3]*erfc(cnl[3]*r); Vr4 = cl[4]*erfc(cnl[4]*r); Vr5 = cl[5]*erfc(cnl[5]*r);
			Vr = (Vr0 + Vr1 + Vr2 + Vr3 + Vr4 + Vr5)/r;
			dVrir = 2.0*(cl[0]*cnl[0]*exp(-cnl[0]*cnl[0]*r2) + cl[1]*cnl[1]*exp(-cnl[1]*cnl[1]*r2) + cl[2]*cnl[2]*exp(-cnl[2]*cnl[2]*r2)+ 
				cl[3]*cnl[3]*exp(-cnl[3]*cnl[3]*r2) + cl[4]*cnl[4]*exp(-cnl[4]*cnl[4]*r2) + cl[5]*cnl[5]*exp(-cnl[5]*cnl[5]*r2))/cPii2 + Vr;
			Vr *= f;
			dVrir *= -f/r2;
			break;
		case 6:
			// 6: lobato parameterization - 5 hydrogen fe - [0, 12]
			r2 = r*r;
			Vr0 = cl[0]*exp(-cnl[0]*r); Vr1 = cl[1]*exp(-cnl[1]*r); Vr2 = cl[2]*exp(-cnl[2]*r); 
			Vr3 = cl[3]*exp(-cnl[3]*r); Vr4 = cl[4]*exp(-cnl[4]*r);
			Vr = Vr0*(2.0*icnl[0]+r) + Vr1*(2.0*icnl[1]+r) + Vr2*(2.0*icnl[2]+r) + Vr3*(2.0*icnl[3]+r) + Vr4*(2.0*icnl[4]+r);
			dVrir = Vr+Vr0*(r+cnl[0]*r2)+Vr1*(r+cnl[1]*r2)+Vr2*(r+cnl[2]*r2)+Vr3*(r+cnl[3]*r2)+Vr4*(r+cnl[4]*r2);
			Vr = f*Vr/r;
			dVrir = -f*dVrir/(r*r2);
			break;
	}
}

// Linear projected potential: V and zV
__global__ void LinearProjAtomicPotentialGPU(sQ1 Qz){	
	__shared__ double V0s[stnQz];
	__shared__ double dV0s[stnQz];
	__shared__ double V1s[stnQz];
	__shared__ double dV1s[stnQz];
	__shared__ double cls[6];
	__shared__ double cnls[6];
	__shared__ double icnls[6];

	unsigned int i = threadIdx.x, iR = blockIdx.x, isatom = blockIdx.y;

	double x = Qz.x[i], w = Qz.w[i], R2 = cVp[isatom].R2[iR];
	double a = (cVp[isatom].split)?-0.5*cVp[isatom].z0:0.5*(cVp[isatom].ze-cVp[isatom].z0);
	double b = (cVp[isatom].split)?0.5*cVp[isatom].z0:0.5*(cVp[isatom].ze+cVp[isatom].z0);
	double z = a*x + b, r = sqrt(z*z + R2);
	double V, dVir;

    if (i < 6){
        cls[i] = cVp[isatom].Vr.cl[i];
        cnls[i] = cVp[isatom].Vr.cnl[i];
        icnls[i] = 1.0/cnls[i];
	}
    __syncthreads();

	Pot3Da(r, cls, cnls, icnls, a*w, V, dVir);
	V0s[i] = V; dV0s[i] = dVir;
	z = z-cVp[isatom].z0;
	V1s[i] = z*V; dV1s[i] = z*dVir;

	if (cVp[isatom].split){
		a = b = 0.5*cVp[isatom].ze;
		z = a*x + b; r = sqrt(z*z + R2);
		Pot3Da(r, cls, cnls, icnls, a*w, V, dVir);
		V0s[i] += V; dV0s[i] += dVir; 
		z = z-cVp[isatom].z0;
		V1s[i] += z*V; dV1s[i] += z*dVir;
	}

	__syncthreads();

    if (i < 64){
        V0s[i] += V0s[i + 64];
		dV0s[i] += dV0s[i + 64];
		V1s[i] += V1s[i + 64];
		dV1s[i] += dV1s[i + 64];
	}
    __syncthreads();

    if (i < 32){
        V0s[i] += V0s[i + 32];
		dV0s[i] += dV0s[i + 32];
		V1s[i] += V1s[i + 32];
		dV1s[i] += dV1s[i + 32];

        V0s[i] += V0s[i + 16];
		dV0s[i] += dV0s[i + 16];
		V1s[i] += V1s[i + 16];
		dV1s[i] += dV1s[i + 16];

        V0s[i] += V0s[i + 8];
		dV0s[i] += dV0s[i + 8];
		V1s[i] += V1s[i + 8];
		dV1s[i] += dV1s[i + 8];

        V0s[i] += V0s[i + 4];
		dV0s[i] += dV0s[i + 4];
		V1s[i] += V1s[i + 4];
		dV1s[i] += dV1s[i + 4];

        V0s[i] += V0s[i + 2];
		dV0s[i] += dV0s[i + 2];
		V1s[i] += V1s[i + 2];
		dV1s[i] += dV1s[i + 2];

        V0s[i] += V0s[i + 1];
		dV0s[i] += dV0s[i + 1];
		V1s[i] += V1s[i + 1];
		dV1s[i] += dV1s[i + 1];

		if(i==0){
			cVp[isatom].ciV0.c0[iR] = cVp[isatom].occ*V0s[0];			// V0
			cVp[isatom].ciV0.c1[iR] = cVp[isatom].occ*0.5*dV0s[0];		// dR2V0

			cVp[isatom].ciV1.c0[iR] = cVp[isatom].occ*V1s[0];			// V1
			cVp[isatom].ciV1.c1[iR] = cVp[isatom].occ*0.5*dV1s[0];		// dR2V1
		}
    }
}

// Get Local interpolation coefficients
__global__ void CubicPolyCoef(void){
	int iR = threadIdx.x, isatom = blockIdx.x;

	if (iR < stnR-1){
		double dx = 1.0/(cVp[isatom].R2[iR+1]-cVp[isatom].R2[iR]), dx2 = dx*dx;
		double V, Vn, dV, dVn, m, n;
		/*********************************************************/
		V = cVp[isatom].ciV0.c0[iR]; Vn = cVp[isatom].ciV0.c0[iR+1];
		dV = cVp[isatom].ciV0.c1[iR]; dVn = cVp[isatom].ciV0.c1[iR+1];
		m = (Vn-V)*dx; n = dV+dVn;
		cVp[isatom].ciV0.c0[iR] = V-cVp[isatom].ciV0.c0[stnR-1];
		cVp[isatom].ciV0.c2[iR] = (3.0*m-n-dV)*dx;
		cVp[isatom].ciV0.c3[iR] = (n-2.0*m)*dx2;
		/*********************************************************/
		V = cVp[isatom].ciV1.c0[iR]; Vn = cVp[isatom].ciV1.c0[iR+1];
		dV = cVp[isatom].ciV1.c1[iR]; dVn = cVp[isatom].ciV1.c1[iR+1];
		m = (Vn-V)*dx; n = dV+dVn;
		cVp[isatom].ciV1.c0[iR] = V-cVp[isatom].ciV1.c0[stnR-1];
		cVp[isatom].ciV1.c2[iR] = (3.0*m-n-dV)*dx;
		cVp[isatom].ciV1.c3[iR] = (n-2.0*m)*dx2;
	}
}

// Cubic polynomial evaluation
__global__ void CubicPolyEval(sGP GP, int nsatom, double *V0g, double *V1g){	
	int j0 = threadIdx.x + blockIdx.x*blockDim.x;
	int i0 = threadIdx.y + blockIdx.y*blockDim.y;

	for(int isatom=0; isatom<nsatom; isatom++){
		if ((i0 < cVp[isatom].bnx.n)&&(j0 < cVp[isatom].bny.n)){	
			int i = cVp[isatom].bnx.i + i0;
			int j = cVp[isatom].bny.i + j0;

			double Rx = (i-GP.nxh)*GP.dRx - cVp[isatom].x;
			double Ry = (j-GP.nyh)*GP.dRy - cVp[isatom].y;
			double R2 = Rx*Rx + Ry*Ry;
			if (R2 < cVp[isatom].Rmax2){
				i = (i<GP.nxh)?(i+GP.nxh):(i-GP.nxh);
				j = (j<GP.nyh)?(j+GP.nyh):(j-GP.nyh);
				int k = i*GP.ny + j, t;
		
				if(R2 <= cVp[isatom].Rmin2){
					i = 0;
					R2 = cVp[isatom].Rmin2;
				}else{
					i = 0; j = stnR-1;
					t = (i + j)>>1;	// divide by 2
					if(R2 < cVp[isatom].R2[t]) j = t; else  i = t; //64
					t = (i + j)>>1;	// divide by 2
					if(R2 < cVp[isatom].R2[t]) j = t; else  i = t; //32
					t = (i + j)>>1;	// divide by 2
					if(R2 < cVp[isatom].R2[t]) j = t; else  i = t; //16
					t = (i + j)>>1;	// divide by 2
					if(R2 < cVp[isatom].R2[t]) j = t; else  i = t; //8
					t = (i + j)>>1;	// divide by 2
					if(R2 < cVp[isatom].R2[t]) j = t; else  i = t; //4
					t = (i + j)>>1;	// divide by 2
					if(R2 < cVp[isatom].R2[t]) j = t; else  i = t; //2
					t = (i + j)>>1;	// divide by 2
					if(R2 < cVp[isatom].R2[t]) j = t; else  i = t; //1
				} 

				double dx = R2 - cVp[isatom].R2[i], dx2 = dx*dx, dx3 = dx2*dx;

				V0g[k] += cVp[isatom].ciV0.c0[i] + cVp[isatom].ciV0.c1[i]*dx + cVp[isatom].ciV0.c2[i]*dx2 + cVp[isatom].ciV0.c3[i]*dx3;
				V1g[k] += cVp[isatom].ciV1.c0[i] + cVp[isatom].ciV1.c1[i]*dx + cVp[isatom].ciV1.c2[i]*dx2 + cVp[isatom].ciV1.c3[i]*dx3;
			}
		}
	}
}

void cPotentialGPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	double PotParh = 0;
	cudaMemcpyToSymbol(PotPar, &PotParh, cSizeofI);

	cudaFreen(Qz.x);
	cudaFreen(Qz.w);

	for(int i=0; i<stncVp; i++){
		cVph[i].x = 0;
		cVph[i].y = 0;
		cVph[i].z = 0;
		cVph[i].occ = 0;
		cVph[i].Rmin2 = 0;
		cVph[i].Rmax2 = 0;

		cVph[i].z0 = 0;
		cVph[i].ze = 0;
		cVph[i].split = false;

		cVph[i].bnx.i = 0;
		cVph[i].bnx.n = 0;
		cVph[i].bnx.i = 0;
		cVph[i].bnx.n = 0;

		cVph[i].R2 = 0;
		cVph[i].Vr.cl = 0;
		cVph[i].Vr.cnl = 0;

		cudaFree_sciVn(cVph[i].ciV0);
		cudaFree_sciVn(cVph[i].ciV1);
	}

	MulOrder = 0;

	nAtomTypesGPU = 0;
	delete [] AtomTypesGPU;  AtomTypesGPU = 0;

	cudaFreen(V0);
	cudaFreen(V1);
	cudaFreen(V2);
}

cPotentialGPU::cPotentialGPU(){
	double PotParh = 0;
	cudaMemcpyToSymbol(PotPar, &PotParh, cSizeofI);

	Qz.x = 0;
	Qz.w = 0;

	for(int i=0; i<stncVp; i++){
		cVph[i].x = 0;
		cVph[i].y = 0;
		cVph[i].z = 0;
		cVph[i].occ = 0;
		cVph[i].Rmin2 = 0;
		cVph[i].Rmax2 = 0;

		cVph[i].z0 = 0;
		cVph[i].ze = 0;
		cVph[i].split = false;

		cVph[i].bnx.i = 0;
		cVph[i].bnx.n = 0;
		cVph[i].bnx.i = 0;
		cVph[i].bnx.n = 0;

		cVph[i].R2 = 0;
		cVph[i].Vr.cl = 0;
		cVph[i].Vr.cnl = 0;

		cudaInit_sciVn(cVph[i].ciV0);
		cudaInit_sciVn(cVph[i].ciV1);
	}

	MulOrder = 0;

	nAtomTypesGPU = 0;
	AtomTypesGPU = 0;

	V0 = 0;
	V1 = 0;
	V2 = 0;
}

cPotentialGPU::~cPotentialGPU(){
	freeMemory();
}

void cPotentialGPU::cudaFree_sciVn(sciVn &ciVn){
	cudaFreen(ciVn.c0);
	cudaFreen(ciVn.c1);
	cudaFreen(ciVn.c2);
	cudaFreen(ciVn.c3);
}

void cPotentialGPU::cudaInit_sciVn(sciVn &ciVn){
	ciVn.c0 = 0;
	ciVn.c1 = 0;
	ciVn.c2 = 0;
	ciVn.c3 = 0;
}

void cPotentialGPU::cudaMalloc_sciVn(int n, sciVn &ciVn){
	cudaMalloc((void**)&(ciVn.c0), n*cSizeofRD);
	cudaMalloc((void**)&(ciVn.c1), n*cSizeofRD);
	cudaMalloc((void**)&(ciVn.c2), n*cSizeofRD);
	cudaMalloc((void**)&(ciVn.c3), n*cSizeofRD);
}

void cPotentialGPU::SetInputData(sGP &GP_i, sBT &BT_i, int nAtomsMi, double *AtomsMi, double dzi, int PotPari, double Vrli, int DimFPi, int DistFPi, int SeedFPi, int ZeroDefTypi, double ZeroDefPlanei, int MulOrder_i){
	freeMemory();

	BT = BT_i;
	GP = GP_i;

	Specimen.SetInputData(nAtomsMi, AtomsMi, dzi, PotPari, Vrli, DimFPi, DistFPi, SeedFPi, ZeroDefTypi, ZeroDefPlanei);

	nAtomTypesGPU = Specimen.nAtomTypes;
	AtomTypesGPU  = new cAtomTypesGPU[nAtomTypesGPU];
	for (int i=0; i<nAtomTypesGPU; i++)
		AtomTypesGPU[i].SetAtomTypes(Specimen.AtomTypes[i], stnR, GP.dRmin);

	MulOrder = MulOrder_i;

	cudaMemcpyToSymbol(PotPar, &PotPari, cSizeofI);

	ReadQuadratureGPU(0, stnQz, Qz);	// TanhSinh

	for(int i=0; i<stncVp; i++){
		cudaMalloc_sciVn(stnR, cVph[i].ciV0);
		cudaMalloc_sciVn(stnR, cVph[i].ciV1);
	}
	cudaMemcpyToSymbol(cVp, cVph, sizeof(scVp)*stncVp);

	cudaMalloc((void**)&V0, GP.nxy*cSizeofRD);
	cudaMalloc((void**)&V1, GP.nxy*cSizeofRD);
	cudaMalloc((void**)&V2, GP.nxy*cSizeofRD);
}

void cPotentialGPU::getbn(sGP &GP, double x, double y, double Rmax, sbn &bnx, sbn &bny){
	double xmin = -GP.nxh*GP.dRx, xmax = (GP.nxh-1)*GP.dRx;
	double ymin = -GP.nyh*GP.dRy, ymax = (GP.nyh-1)*GP.dRy;
	int i0, ie;

	i0 = (x-Rmax<xmin)?0:(int)ceil((x-Rmax-xmin)/GP.dRx);
	if(i0>=GP.nx) i0 = GP.nx-1;
	ie = (x+Rmax>xmax)?GP.nx-1:(int)floor((x+Rmax-xmin)/GP.dRx);
	if(ie<0) ie = 0;
	bnx.i = i0; bnx.n =(i0==ie)?0:ie-i0+1;

	i0 = (y-Rmax<ymin)?0:(int)ceil((y-Rmax-ymin)/GP.dRy);
	if(i0>=GP.ny) i0 = GP.ny-1;
	ie = (y+Rmax>ymax)?GP.ny-1:(int)floor((y+Rmax-ymin)/GP.dRy);
	if(ie<0) ie = 0;
	bny.i = i0; bny.n =(i0==ie)?0:ie-i0+1;
};

void cPotentialGPU::setcVp(cSpecimen &Specimen, cAtomTypesGPU *&AtomTypesGPU, int iSli, int iatom, int nsatom, dim3 &BPot, dim3 &TPot, dim3 &BCoef, dim3 &TCoef, dim3 &BEval, dim3 &TEval){
	int iZ;
	int nxm, nym;

	nxm = nym = 0;
	for(int i=0; i<nsatom; i++){
		iZ = Specimen.Atoms[iatom+i].iZ;
		cVph[i].x = Specimen.Atoms[iatom+i].x;
		cVph[i].y = Specimen.Atoms[iatom+i].y;
		cVph[i].z = Specimen.Atoms[iatom+i].z;
		cVph[i].occ =  Specimen.Atoms[iatom+i].occ;
		cVph[i].Rmin2 = AtomTypesGPU[iZ].Rmin2;
		cVph[i].Rmax2 = AtomTypesGPU[iZ].Rmax2;
		cVph[i].z0 = Specimen.Slice[iSli].z0 - cVph[i].z; 
		cVph[i].ze = Specimen.Slice[iSli].ze - cVph[i].z;
		cVph[i].split = (cVph[i].z0<0)&&(0<cVph[i].ze);
		getbn(GP, cVph[i].x, cVph[i].y, AtomTypesGPU[iZ].Rmax, cVph[i].bnx, cVph[i].bny);
		if(nxm < cVph[i].bnx.n)
			nxm = cVph[i].bnx.n;
		if(nym < cVph[i].bny.n)
			nym = cVph[i].bny.n;
		cVph[i].Vr.cl = AtomTypesGPU[iZ].Vr.cl;
		cVph[i].Vr.cnl = AtomTypesGPU[iZ].Vr.cnl;
		cVph[i].R2 = AtomTypesGPU[iZ].R2;
	}
	cudaMemcpyToSymbolAsync(cVp, cVph, sizeof(scVp)*nsatom);

	TPot.x = stnQz; TPot.y = 1; TPot.z = 1;
	BPot.x = stnR; BPot.y = nsatom; BPot.z = 1;

	TCoef.x = stnR; TCoef.y = 1; TCoef.z = 1;
	BCoef.x = nsatom; BCoef.y = 1; BCoef.z = 1;

	TEval.x = thrnxny; TEval.y = thrnxny; TEval.z = 1;
	BEval.x = (nym+thrnxny-1)/thrnxny; BEval.y = (nxm+thrnxny-1)/thrnxny; BEval.z = 1;	
}

// Projected potential calculation: iSli = slice position
void cPotentialGPU::ProjectedPotential(int iSli){	
	int iatom, nsatom, iSynCPU, cSynCPU = 16;
	double iatom0 = Specimen.Slice[iSli].z0i_id;
	double iatome = Specimen.Slice[iSli].zei_id;
	dim3 BPot, TPot, BCoef, TCoef, BEval, TEval;

	// Set V0 = 0 and V1 = 0 for all its elements
	SetValueVectorD(BT.Bnxy, BT.Tnxy, GP.nxy, 0.0, V0, V1);
	double dz = Specimen.Slice[iSli].ze - Specimen.Slice[iSli].z0;

	iSynCPU = 0;
	iatom = iatom0;
	while (iatom<=iatome){
		nsatom = MIN(stncVp, iatome-iatom+1);
		setcVp(Specimen, AtomTypesGPU, iSli, iatom, nsatom, BPot, TPot, BCoef, TCoef, BEval, TEval);
		LinearProjAtomicPotentialGPU<<<BPot, TPot>>>(Qz);
		CubicPolyCoef<<<BCoef, TCoef>>>();
		CubicPolyEval<<<BEval, TEval>>>(GP, nsatom, V0, V1);	
		iatom += nsatom;
		iSynCPU++;
		if(iSynCPU%cSynCPU==0)
			cudaDeviceSynchronize();
	}

	if (MulOrder==2)
		ScaleVectorD(BT.Bnxy, BT.Tnxy, GP.nxy, 1.0/dz, V1);
}