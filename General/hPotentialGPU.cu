#include "hConstTypes.h"
#include "hgeneralCPU.h"
#include "hgeneralGPU.h"
#include "hSpecimenCPU.h"
#include "hAtomTypesGPU.h"
#include "hPotentialGPU.h"
#include "hQuadrature.h"
#include "math.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>

__constant__ int PotPar = 0;
__constant__ scVp cVp[stncVp];

// Potential Evaluation (VR, dVRiR)
__device__ inline void Pot2Da(double &R, volatile double *cl, volatile double *cnl, volatile double *icnl, double f, double &VR, double &dVRiR){
	double R2;
	double VR0, VR1, VR2, VR3, VR4, VR5;

	switch (PotPar){
		case 1:
			// 1: doyle and turner parameterization - 4 gaussians - [0, 4]
			R2 = R*R;
			VR0 = cl[0]*exp(-cnl[0]*R2); VR1 = cl[1]*exp(-cnl[1]*R2); VR2 = cl[2]*exp(-cnl[2]*R2); 
			VR3 = cl[3]*exp(-cnl[3]*R2);
			VR = VR0 + VR1 + VR2 + VR3;
			dVRiR = 2.0*(cnl[0]*VR0 + cnl[1]*VR1 + cnl[2]*VR2 + cnl[3]*VR3);
			VR *= f;
			dVRiR *= -f;
			break;
		case 2:
			// 2: peng et al. parameterization - 5 gaussians - [0, 4]
			R2 = R*R;
			VR0 = cl[0]*exp(-cnl[0]*R2); VR1 = cl[1]*exp(-cnl[1]*R2); VR2 = cl[2]*exp(-cnl[2]*R2); 
			VR3 = cl[3]*exp(-cnl[3]*R2); VR4 = cl[4]*exp(-cnl[4]*R2);
			VR = VR0 + VR1 + VR2 + VR3 + VR4;
			dVRiR = 2.0*(cnl[0]*VR0 + cnl[1]*VR1 + cnl[2]*VR2 + cnl[3]*VR3 + cnl[4]*VR4);
			VR *= f;
			dVRiR *= -f;
			break;
		case 3:		
			// 3: peng et al. parameterization - 5 gaussians - [0, 12]
			R2 = R*R;
			VR0 = cl[0]*exp(-cnl[0]*R2); VR1 = cl[1]*exp(-cnl[1]*R2); VR2 = cl[2]*exp(-cnl[2]*R2); 
			VR3 = cl[3]*exp(-cnl[3]*R2); VR4 = cl[4]*exp(-cnl[4]*R2);
			VR = VR0 + VR1 + VR2 + VR3 + VR4;
			dVRiR = 2.0*(cnl[0]*VR0 + cnl[1]*VR1 + cnl[2]*VR2 + cnl[3]*VR3 + cnl[4]*VR4);
			VR *= f;
			dVRiR *= -f;
			break;
		case 4:
			// 4: Kirkland parameterization - 3 Yukawa + 3 Gaussians - [0, 12]
			R2 = R*R;
			VR0 = cl[0]*bessk0GPU(cnl[0]*R); VR1 = cl[1]*bessk0GPU(cnl[1]*R); VR2 = cl[2]*bessk0GPU(cnl[2]*R); 
			VR3 = cl[3]*exp(-cnl[3]*R2); VR4 = cl[4]*exp(-cnl[4]*R2); VR5 = cl[5]*exp(-cnl[5]*R2);
			VR = VR0 + VR1 + VR2 + VR3 + VR4 + VR5;
			dVRiR = cl[0]*cnl[0]*bessk1GPU(cnl[0]*R) + cl[1]*cnl[1]*bessk1GPU(cnl[1]*R) + cl[2]*cnl[2]*bessk1GPU(cnl[2]*R) + 2.0*R*(cnl[3]*VR3 + cnl[4]*VR4 + cnl[5]*VR5);
			VR *= f;
			dVRiR *= -f/R;
		case 5:
			// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
			VR = 0;
			dVRiR = 0;
			break;
		case 6:
			// 6: lobato parameterization - 5 hydrogen fe - [0, 12]
			R2 = R*R;
			VR0 = cl[0]*bessk0GPU(cnl[0]*R); VR1 = bessk0GPU(cnl[1]*R); VR2 = bessk0GPU(cnl[2]*R); 
			VR3 = bessk0GPU(cnl[3]*R); VR4 = bessk0GPU(cnl[4]*R);
			VR = 2.0*(cl[0]*icnl[0]*VR0 + cl[1]*icnl[1]*VR1 + cl[2]*icnl[2]*VR2 + cl[3]*icnl[3]*VR3 + cl[4]*icnl[4]*VR4);
			dVRiR = R*(cl[0]*cnl[0]*VR0 + cl[1]*cnl[1]*VR1 + cl[2]*cnl[2]*VR2 + cl[3]*cnl[3]*VR3 + cl[4]*cnl[4]*VR4);

			VR0 = bessk1GPU(cnl[0]*R); VR1 = bessk1GPU(cnl[1]*R); VR2 = bessk1GPU(cnl[2]*R); 
			VR3 = bessk1GPU(cnl[3]*R); VR4 = bessk1GPU(cnl[4]*R);
			VR += R*(cl[0]*VR0 + cl[1]*VR1 + cl[2]*VR2 + cl[3]*VR3 + cl[4]*VR4);
			dVRiR += 2.0*(cl[0]*VR0 + cl[1]*VR1 + cl[2]*VR2 + cl[3]*VR3 + cl[4]*VR4);
			VR = f*VR;
			dVRiR = -f*dVRiR/R;
			break;
	}
}

// Potential Evaluation (Vr, dVrir)
__device__ inline void Pot3Da(double &r, volatile double *cl, volatile double *cnl, volatile double *icnl, double f, double &Vr, double &dVrir){
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

/************************************************************************************/
// Linear projected potential: V and zV
__global__ void LinearProjAtomicPotentialGPU(sQ1 Qz)
{	
	__shared__ double V0s[stnQz];
	__shared__ double dV0s[stnQz];
	__shared__ double V1s[stnQz];
	__shared__ double dV1s[stnQz];
	__shared__ double cls[6];
	__shared__ double cnls[6];
	__shared__ double icnls[6];

	int ix = threadIdx.x, iR = blockIdx.x, isatom = blockIdx.y;

	double x = Qz.x[ix], w = Qz.w[ix], R2 = cVp[isatom].R2[iR];
	double a = (cVp[isatom].split)?(-0.5*cVp[isatom].z0):(0.5*(cVp[isatom].ze-cVp[isatom].z0));
	double b = (cVp[isatom].split)?(0.5*cVp[isatom].z0):(0.5*(cVp[isatom].ze+cVp[isatom].z0));
	double z = a*x + b, r = sqrt(z*z + R2);
	double V, dVir;

	if (ix < 6){
		 cls[ix] = cVp[isatom].cVr.cl[ix];
		 cnls[ix] = cVp[isatom].cVr.cnl[ix];
		 icnls[ix] = (cnls[ix]>0)?1.0/cnls[ix]:0.0;
	}
	__syncthreads();

	Pot3Da(r, cls, cnls, icnls, a*w, V, dVir);
	V0s[ix] = V; dV0s[ix] = dVir;
	z = z-cVp[isatom].z0;
	V1s[ix] = z*V; dV1s[ix] = z*dVir;

	if (cVp[isatom].split){
		a = b = 0.5*cVp[isatom].ze;
		z = a*x + b; r = sqrt(z*z + R2);
		Pot3Da(r, cls, cnls, icnls, a*w, V, dVir);
		V0s[ix] += V; dV0s[ix] += dVir;
		z = z-cVp[isatom].z0;
		V1s[ix] += z*V; dV1s[ix] += z*dVir;
	}

	__syncthreads();

	reduceBlockDouble128(V0s, dV0s, V1s, dV1s, ix);

	if(ix==0){
		cVp[isatom].ciV0.c0[iR] = V0s[0];			// V0
		cVp[isatom].ciV0.c1[iR] = 0.5*dV0s[0];		// dR2V0
		cVp[isatom].ciV1.c0[iR] = V1s[0];			// V1
		cVp[isatom].ciV1.c1[iR] = 0.5*dV1s[0];		// dR2V1
	}
}

void LinearProjAtomicPotentialGPU(dim3 grid, dim3 threads, sQ1 &Qz){
	LinearProjAtomicPotentialGPU<<<grid, threads>>>(Qz);
}

/************************************************************************************/
// Get Local interpolation coefficients
template <int MulOrder>
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
		if(MulOrder>=2){
			V = cVp[isatom].ciV1.c0[iR]; Vn = cVp[isatom].ciV1.c0[iR+1];
			dV = cVp[isatom].ciV1.c1[iR]; dVn = cVp[isatom].ciV1.c1[iR+1];
			m = (Vn-V)*dx; n = dV+dVn;
			cVp[isatom].ciV1.c0[iR] = V-cVp[isatom].ciV1.c0[stnR-1];
			cVp[isatom].ciV1.c2[iR] = (3.0*m-n-dV)*dx;
			cVp[isatom].ciV1.c3[iR] = (n-2.0*m)*dx2;
		}
	}
}

void CubicPolyCoef(dim3 grid, dim3 threads, int MulOrder){
	switch (MulOrder)
	{
		case 1:
			CubicPolyCoef<1><<<grid, threads>>>();
			break;
		case 2:
			CubicPolyCoef<2><<<grid, threads>>>();
			break;
	} 
}

/************************************************************************************/
// Cubic polynomial evaluation
template <int MulOrder>
__global__ void CubicPolyEval(sGP GP, int nsatom, double *V0g, double *V1g){	
	int iy0 = threadIdx.x + blockIdx.x*blockDim.x;
	int ix0 = threadIdx.y + blockIdx.y*blockDim.y;
	int ix, iy, k;
	double Rx, Ry, R2;
	double dx, dx2, dx3;
	double V0, V1;

	for(int isatom=0; isatom<nsatom; isatom++)
		if ((ix0 < cVp[isatom].bnx.n)&&(iy0 < cVp[isatom].bny.n)){	
			ix = cVp[isatom].bnx.i + ix0;
			iy = cVp[isatom].bny.i + iy0;

			Rx = ix*GP.dRx - cVp[isatom].x;
			Ry = iy*GP.dRy - cVp[isatom].y;
			R2 = Rx*Rx + Ry*Ry;
			if (R2 < cVp[isatom].Rmax2){
				ix = ix - GP.nx*(int)floor((ix*GP.dRx)/GP.lx);
				iy = iy - GP.ny*(int)floor((iy*GP.dRy)/GP.ly);
				ix = IsRS(ix, GP.nxh);
				iy = IsRS(iy, GP.nyh);
				k = ix*GP.ny + iy;
		
				if(R2 <= cVp[isatom].Rmin2){
					ix = 0;
					R2 = cVp[isatom].Rmin2;
				}else
					ix = unrolledBinarySearch128(R2, cVp[isatom].R2);

				dx = R2 - cVp[isatom].R2[ix]; dx2 = dx*dx; dx3 = dx2*dx;
				V0 = cVp[isatom].ciV0.c0[ix] + cVp[isatom].ciV0.c1[ix]*dx + cVp[isatom].ciV0.c2[ix]*dx2 + cVp[isatom].ciV0.c3[ix]*dx3;
				atomicAdd(&V0g[k], cVp[isatom].occ*V0);
				if( MulOrder>=2){
					V1 = cVp[isatom].ciV1.c0[ix] + cVp[isatom].ciV1.c1[ix]*dx + cVp[isatom].ciV1.c2[ix]*dx2 + cVp[isatom].ciV1.c3[ix]*dx3;
					atomicAdd(&V1g[k], cVp[isatom].occ*V1);
				}
			}
		}
}

void CubicPolyEval(dim3 grid, dim3 threads, int MulOrder, sGP &GP, int nsatom, double *&V0g, double *&V1g){
	switch (MulOrder)
	{
		case 1:
			CubicPolyEval<1><<<grid, threads>>>(GP, nsatom, V0g, V1g);
			break;
		case 2:
			CubicPolyEval<2><<<grid, threads>>>(GP, nsatom, V0g, V1g);
			break;
	} 
}

/************************************************************************************/
void cPotentialGPU::SetPotPar(int PotParh){
	cudaMemcpyToSymbol(PotPar, &PotParh, cSizeofI);
}

void cPotentialGPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	SetPotPar(0);

	cudaFreen(Qz.x);
	cudaFreen(Qz.w);
 
	for(int i=0; i<stncVp; i++){
		cVph[i].x = 0;
		cVph[i].y = 0;
		cVph[i].z0 = 0;
		cVph[i].ze = 0;
		cVph[i].split = false;
		cVph[i].occ = 0;
		cVph[i].Rmin2 = 0;
		cVph[i].Rmax2 = 0;
		cVph[i].R2 = 0;
		cVph[i].cVr.cl = 0;
		cVph[i].cVr.cnl = 0;

		f_sciVn_cudaFree(cVph[i].ciV0);
		f_sciVn_cudaFree(cVph[i].ciV1);

		cVph[i].bnx.i = 0;
		cVph[i].bnx.n = 0;
		cVph[i].bny.i = 0;
		cVph[i].bny.n = 0;
	}
	delete [] cVph; cVph = 0;

	f_sMPG_Init(MGP);
	f_sGP_Init(GP);
	f_sBT_Init(BT);

	nAtomTypesGPU = 0;
	delete [] AtomTypesGPU; AtomTypesGPU = 0;

	cudaFreen(V0);
	cudaFreen(V1);
	cudaFreen(V2);
}

cPotentialGPU::cPotentialGPU(){
	SetPotPar(0);

	Qz.x = 0;
	Qz.w = 0;

	cVph = new scVp[stncVp];

	for(int i=0; i<stncVp; i++){
		cVph[i].x = 0;
		cVph[i].y = 0;
		cVph[i].z0 = 0;
		cVph[i].ze = 0;
		cVph[i].split = false;
		cVph[i].occ = 0;
		cVph[i].Rmin2 = 0;
		cVph[i].Rmax2 = 0;
		cVph[i].R2 = 0;
		cVph[i].cVr.cl = 0;
		cVph[i].cVr.cnl = 0;

		f_sciVn_cudaInit(cVph[i].ciV0);
		f_sciVn_cudaInit(cVph[i].ciV1);

		cVph[i].bnx.i = 0;
		cVph[i].bnx.n = 0;
		cVph[i].bny.i = 0;
		cVph[i].bny.n = 0;
	}

	f_sMPG_Init(MGP);
	f_sGP_Init(GP);
	f_sBT_Init(BT);

	nAtomTypesGPU = 0;
	AtomTypesGPU = 0;

	V0 = 0;
	V1 = 0;
	V2 = 0;
}

cPotentialGPU::~cPotentialGPU(){
	freeMemory();
}

void cPotentialGPU::SetInputData(sMGP &MGP_io, sGP &GP_i, sBT &BT_i, int nAtomsM_i, double *AtomsM_i){
	freeMemory();

	MGP = MGP_io;
	BT = BT_i;
	GP = GP_i;

	SetPotPar(MGP.PotPar);
	f_ReadQuadratureGPU(0, stnQz, Qz);	// TanhSinh

	Specimen.SetInputData(MGP, nAtomsM_i, AtomsM_i);

	nAtomTypesGPU = Specimen.nAtomTypes;
	AtomTypesGPU = new cAtomTypesGPU[nAtomTypesGPU];
	for (int iAtomTypesGPU=0; iAtomTypesGPU<nAtomTypesGPU; iAtomTypesGPU++)
		AtomTypesGPU[iAtomTypesGPU].SetAtomTypes(MGP.PotPar, Specimen.AtomTypes[iAtomTypesGPU], stnR, GP.dRmin);

	cVph = new scVp[stncVp];
	for(int i=0; i<stncVp; i++){
		cVph[i].x = 0;
		cVph[i].y = 0;
		cVph[i].z0 = 0;
		cVph[i].ze = 0;
		cVph[i].split = false;
		cVph[i].occ = 0;
		cVph[i].Rmin2 = 0;
		cVph[i].Rmax2 = 0;

		cVph[i].R2 = 0;
		f_sCoefPar_Init(cVph[i].cVr);
		f_sciVn_cudaMalloc(stnR, cVph[i].ciV0);
		f_sciVn_cudaMalloc(stnR, cVph[i].ciV1);

		cVph[i].bnx.i = 0;
		cVph[i].bnx.n = 0;
		cVph[i].bny.i = 0;
		cVph[i].bny.n = 0;
	}
	cudaMemcpyToSymbol(cVp, cVph, sizeof(scVp)*stncVp, 0, cudaMemcpyHostToDevice);

	cudaMalloc((void**)&V0, GP.nxy*cSizeofRD);
	cudaMalloc((void**)&V1, GP.nxy*cSizeofRD);
	cudaMalloc((void**)&V2, GP.nxy*cSizeofRD);
	/**************************************************************************/
	MGP_io = MGP;
}

int cPotentialGPU::CheckGridLimits(int i, int n){
	return (i<0)?0:((i>=n)?n-1:i);
}

void cPotentialGPU::getbn(sGP &GP, double x, double y, double Rmax, sbn &bnx, sbn &bny){
	double x0 = x-Rmax, xe = x+Rmax;
	double y0 = y-Rmax, ye = y+Rmax;
	int ix0 = (int)floor(x0/GP.dRx), ixe = (int)ceil(xe/GP.dRx);
	int iy0 = (int)floor(y0/GP.dRy), iye = (int)ceil(ye/GP.dRy);

	if(!GP.PBC_xy){
		ix0 = CheckGridLimits(ix0, GP.nx);
		ixe = CheckGridLimits(ixe, GP.nx);
		iy0 = CheckGridLimits(iy0, GP.ny);
		iye = CheckGridLimits(iye, GP.ny);
	}
	bnx.i = ix0; bnx.n =(ix0==ixe)?0:ixe-ix0+1;
	bny.i = iy0; bny.n =(iy0==iye)?0:iye-iy0+1;
};

void cPotentialGPU::setcVp(int ApproxModel, cSpecimenCPU &Specimen, cAtomTypesGPU *&AtomTypesGPU, int iSlice, int iatom, int nsatom, dim3 &BPot, dim3 &TPot, dim3 &BCoef, dim3 &TCoef, dim3 &BEval, dim3 &TEval){
	int Z;
	int nxmax = 0, nymax = 0;

	for(int i=0; i<nsatom; i++){
		Z = Specimen.Atoms[iatom+i].Z-1;
		cVph[i].x = Specimen.Atoms[iatom+i].x;
		cVph[i].y = Specimen.Atoms[iatom+i].y;
		cVph[i].z0 = Specimen.Slice[iSlice].z0 - Specimen.Atoms[iatom+i].z; 
		cVph[i].ze = Specimen.Slice[iSlice].ze - Specimen.Atoms[iatom+i].z;
		cVph[i].split = (cVph[i].z0<0)&&(0<cVph[i].ze);
		cVph[i].occ = Specimen.Atoms[iatom+i].occ;
		cVph[i].Rmin2 = AtomTypesGPU[Z].Rmin2;
		cVph[i].Rmax2 = AtomTypesGPU[Z].Rmax2;
		cVph[i].cVr.cl = AtomTypesGPU[Z].cVr.cl;
		cVph[i].cVr.cnl = AtomTypesGPU[Z].cVr.cnl;
		cVph[i].R2 = AtomTypesGPU[Z].R2;
		getbn(GP, cVph[i].x, cVph[i].y, AtomTypesGPU[Z].Rmax, cVph[i].bnx, cVph[i].bny);
		if(nxmax<cVph[i].bnx.n) nxmax = cVph[i].bnx.n;
		if(nymax<cVph[i].bny.n) nymax = cVph[i].bny.n;
		if(ApproxModel>1){
			cudaMemcpyAsync(cVph[i].ciV0.c0, AtomTypesGPU[Z].ciVR.c0, stnR*cSizeofRD, cudaMemcpyDeviceToDevice);
			cudaMemcpyAsync(cVph[i].ciV0.c1, AtomTypesGPU[Z].ciVR.c1, stnR*cSizeofRD, cudaMemcpyDeviceToDevice);
			cudaMemcpyAsync(cVph[i].ciV0.c2, AtomTypesGPU[Z].ciVR.c2, stnR*cSizeofRD, cudaMemcpyDeviceToDevice);
			cudaMemcpyAsync(cVph[i].ciV0.c3, AtomTypesGPU[Z].ciVR.c3, stnR*cSizeofRD, cudaMemcpyDeviceToDevice);
		}
	}
	cudaMemcpyToSymbolAsync(cVp, cVph, sizeof(scVp)*stncVp, 0, cudaMemcpyHostToDevice);

	TPot.x = stnQz; TPot.y = 1; TPot.z = 1;
	BPot.x = stnR; BPot.y = nsatom; BPot.z = 1;

	TCoef.x = stnR; TCoef.y = 1; TCoef.z = 1;
	BCoef.x = nsatom; BCoef.y = 1; BCoef.z = 1;

	TEval.x = thrnxny;  TEval.y = thrnxny; TEval.z = 1;
	BEval.x = (nxmax+thrnxny-1)/thrnxny;  BEval.y = (nymax+thrnxny-1)/thrnxny; BEval.z = 1;
}

// Projected potential calculation: iSlice = slice position
void cPotentialGPU::ProjectedPotential_Slice(int iSlice){	
	int iatom, nsatom;
	int iatom0 = Specimen.Slice[iSlice].z0i_id;
	int iatome = Specimen.Slice[iSlice].zei_id;
	dim3 BPot, TPot, BCoef, TCoef, BEval, TEval;

	// Set V0 = 0 and V1 = 0 for all its elements
	f_SetValueVectorD(BT.Bnxy, BT.Tnxy, GP.nxy, 0.0, V0, V1);
	iatom = iatom0;
	while (iatom<=iatome){
		nsatom = MIN(stncVp, iatome-iatom+1);
		setcVp(MGP.ApproxModel, Specimen, AtomTypesGPU, iSlice, iatom, nsatom, BPot, TPot, BCoef, TCoef, BEval, TEval);
		LinearProjAtomicPotentialGPU(BPot, TPot, Qz);
		CubicPolyCoef(BCoef, TCoef, MGP.MulOrder);
		CubicPolyEval(BEval, TEval, MGP.MulOrder, GP, nsatom, V0, V1);
		iatom += nsatom;
	}

	double dz = Specimen.Slice[iSlice].ze - Specimen.Slice[iSlice].z0;
	if (MGP.MulOrder>=2)
		f_ScaleVectorD(BT.Bnxy, BT.Tnxy, GP.nxy, 1.0/dz, V1);
}

// Projected potential calculation: iSlice = slice position
void cPotentialGPU::ProjectedPotential_Specimen(int iSlice){	
	int iatom, nsatom;
	int iatom0 = Specimen.Slice[iSlice].z0i_id;
	int iatome = Specimen.Slice[iSlice].zei_id;
	dim3 BPot, TPot, BCoef, TCoef, BEval, TEval;

	// Set V0 = 0 and V1 = 0 for all its elements
	f_SetValueVectorD(BT.Bnxy, BT.Tnxy, GP.nxy, 0.0, V0, V1);

	iatom = iatom0;
	while (iatom<=iatome){
		nsatom = MIN(stncVp, iatome-iatom+1);
		setcVp(MGP.ApproxModel, Specimen, AtomTypesGPU, iSlice, iatom, nsatom, BPot, TPot, BCoef, TCoef, BEval, TEval);
		CubicPolyEval(BEval, TEval, MGP.MulOrder, GP, nsatom, V0, V1);
		iatom += nsatom;
	}
}