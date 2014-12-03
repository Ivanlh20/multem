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

#include "hConstTypes.h"
#include "hQuadrature.h"
#include "hMT_General_CPU.h"
#include "hMT_General_GPU.h"
#include "hMT_MGP_CPU.h"
#include "hMT_Specimen_CPU.h"
#include "hMT_AtomTypes_GPU.h"
#include "hMT_Potential_GPU.h"

#include "math.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>

__constant__ int PotPar = 0;
__constant__ scVp cVp[stncVp];
scVp cVph[stncVp];

/*------------------------------------------------------------*/
/* PURPOSE: Atomic addition double values					*/
/*------------------------------------------------------------*/
__device__ inline void atomicAdd(double *address, double val)
{
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed,__double_as_longlong(val +
		 __longlong_as_double(assumed)));
	} while (assumed != old);
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=0.*/
/*------------------------------------------------------------*/
__device__ inline double bessi0GPU(double x)
{
	 double ax,ans;
	 double y;

	 if ((ax=fabs(x)) < 3.75){
		 y=x/3.75,y=y*y;
		 ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
		 +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	 }else{
		 y=3.75/ax;
		 ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
		 +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
		 +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
		 +y*0.392377e-2))))))));
	 }
	 return ans;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=1. */
/*------------------------------------------------------------*/
__device__ inline double bessi1GPU(double x)
{
	 double ax,ans;
	 double y;

	 if ((ax=fabs(x)) < 3.75){
		 y=x/3.75,y=y*y;
		 ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
		 +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
	 }else{
		 y=3.75/ax;
		 ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
		 -y*0.420059e-2));
		 ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
		 +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
		 ans *= (exp(ax)/sqrt(ax));
	 }
	 return x < 0.0 ? -ans : ans;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=0.*/
/*------------------------------------------------------------*/
__device__ inline double bessk0GPU(double x)
{
	 double y,ans;

	 if (x <= 2.0){
		 y=x*x/4.0;
		 ans=(-log(x/2.0)*bessi0GPU(x))+(-0.57721566+y*(0.42278420
		 +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
		 +y*(0.10750e-3+y*0.74e-5))))));
	 }else{
		 y=2.0/x;
		 ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
		 +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
		 +y*(-0.251540e-2+y*0.53208e-3))))));
	 }
	 return ans;
}

/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=1.*/
/*------------------------------------------------------------*/
__device__ inline double bessk1GPU(double x)
{
	 double y,ans;

	 if (x <= 2.0){
		 y=x*x/4.0;
		 ans=(log(x/2.0)*bessi1GPU(x))+(1.0/x)*(1.0+y*(0.15443144
		 +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
		 +y*(-0.110404e-2+y*(-0.4686e-4)))))));
	 }else{
		 y=2.0/x;
		 ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
		 +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
		 +y*(0.325614e-2+y*(-0.68245e-3)))))));
	 }
	 return ans;
}

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

__device__ inline int unrolledBinarySearch128(double x0, const double * __restrict x){
	int i0 = 0, ie = stnR-1;
	int im = (i0 + ie)>>1;	// divide by 2
	if(x0 < x[im]) ie = im; else i0 = im; //64
	im = (i0 + ie)>>1;	// divide by 2
	if(x0 < x[im]) ie = im; else i0 = im; //32
	im = (i0 + ie)>>1;	// divide by 2
	if(x0 < x[im]) ie = im; else i0 = im; //16
	im = (i0 + ie)>>1;	// divide by 2
	if(x0 < x[im]) ie = im; else i0 = im; //8
	im = (i0 + ie)>>1;	// divide by 2
	if(x0 < x[im]) ie = im; else i0 = im; //4
	im = (i0 + ie)>>1;	// divide by 2
	if(x0 < x[im]) ie = im; else i0 = im; //2
	im = (i0 + ie)>>1;	// divide by 2
	if(x0 < x[im]) ie = im; else i0 = im; //1
	
	return i0;
}

__device__ inline void reduceBlockDouble128(volatile double *M1, volatile double *M2, int tid)
{ 
	if (tid < 64){
		M1[tid] += M1[tid + 64];
		M2[tid] += M2[tid + 64];
	}
	__syncthreads();

	if (tid < 32){
		M1[tid] += M1[tid + 32];
		M2[tid] += M2[tid + 32];

		M1[tid] += M1[tid + 16];
		M2[tid] += M2[tid + 16];

		M1[tid] += M1[tid + 8];
		M2[tid] += M2[tid + 8];

		M1[tid] += M1[tid + 4];
		M2[tid] += M2[tid + 4];

		M1[tid] += M1[tid + 2];
		M2[tid] += M2[tid + 2];

		M1[tid] += M1[tid + 1];
		M2[tid] += M2[tid + 1];
	}
}

__device__ inline void reduceBlockDouble128(volatile double *M1, volatile double *M2, volatile double *M3, volatile double *M4, int tid)
{ 
	if (tid < 64){
		M1[tid] += M1[tid + 64];
		M2[tid] += M2[tid + 64];
		M3[tid] += M3[tid + 64];
		M4[tid] += M4[tid + 64];
	}
	__syncthreads();

	if (tid < 32){
		M1[tid] += M1[tid + 32];
		M2[tid] += M2[tid + 32];
		M3[tid] += M3[tid + 32];
		M4[tid] += M4[tid + 32];

		M1[tid] += M1[tid + 16];
		M2[tid] += M2[tid + 16];
		M3[tid] += M3[tid + 16];
		M4[tid] += M4[tid + 16];

		M1[tid] += M1[tid + 8];
		M2[tid] += M2[tid + 8];
		M3[tid] += M3[tid + 8];
		M4[tid] += M4[tid + 8];

		M1[tid] += M1[tid + 4];
		M2[tid] += M2[tid + 4];
		M3[tid] += M3[tid + 4];
		M4[tid] += M4[tid + 4];

		M1[tid] += M1[tid + 2];
		M2[tid] += M2[tid + 2];
		M3[tid] += M3[tid + 2];
		M4[tid] += M4[tid + 2];

		M1[tid] += M1[tid + 1];
		M2[tid] += M2[tid + 1];
		M3[tid] += M3[tid + 1];
		M4[tid] += M4[tid + 1];
	}
}

/***************************************************************************/
/***************************************************************************/

// Linear projected potential: V and zV
template <int MulOrder>
__global__ void k_LinearProjAtomicPotentialGPU(sQ1 Qz)
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
	double z = a*x + b, r = sqrt(z*z + R2), zm = 0.5*(cVp[isatom].ze+cVp[isatom].z0);
	double V, dVir;

	if (ix < 6){
		 cls[ix] = cVp[isatom].cVr.cl[ix];
		 cnls[ix] = cVp[isatom].cVr.cnl[ix];
		 icnls[ix] = (cnls[ix]>0)?1.0/cnls[ix]:0.0;
	}
	__syncthreads();

	Pot3Da(r, cls, cnls, icnls, a*w, V, dVir);
	V0s[ix] = V; dV0s[ix] = dVir;
	if(MulOrder>=2){
		z = z-zm; V1s[ix] = z*V; dV1s[ix] = z*dVir;
	}

	if (cVp[isatom].split){
		a = b = 0.5*cVp[isatom].ze;
		z = a*x + b; r = sqrt(z*z + R2);
		Pot3Da(r, cls, cnls, icnls, a*w, V, dVir);
		V0s[ix] += V; dV0s[ix] += dVir;
		if(MulOrder>=2){
			z = z-zm; V1s[ix] += z*V; dV1s[ix] += z*dVir;
		}
	}

	__syncthreads();

	switch(MulOrder){
		case 1:
			reduceBlockDouble128(V0s, dV0s, ix);
			break;
		case 2:
			reduceBlockDouble128(V0s, dV0s, V1s, dV1s, ix);
			break;
	}
		

	if(ix==0){
		cVp[isatom].ciV0.c0[iR] = V0s[0];			// V0
		cVp[isatom].ciV0.c1[iR] = 0.5*dV0s[0];		// dR2V0
		if(MulOrder>=2){
			cVp[isatom].ciV1.c0[iR] = V1s[0];			// V1
			cVp[isatom].ciV1.c1[iR] = 0.5*dV1s[0];		// dR2V1
		}
	}
}

// Get Local interpolation coefficients
template <int MulOrder>
__global__ void k_CubicPolyCoef(void){
	int iR = threadIdx.x, isatom = blockIdx.x;

	if (iR < stnR-1){
		double dx = 1.0/(cVp[isatom].R2[iR+1]-cVp[isatom].R2[iR]), dx2 = dx*dx;
		double V, Vn, dV, dVn, m, n;
		/********************************************************/
		V = cVp[isatom].ciV0.c0[iR]; Vn = cVp[isatom].ciV0.c0[iR+1];
		dV = cVp[isatom].ciV0.c1[iR]; dVn = cVp[isatom].ciV0.c1[iR+1];
		m = (Vn-V)*dx; n = dV+dVn;
		cVp[isatom].ciV0.c0[iR] = V-cVp[isatom].ciV0.c0[stnR-1];
		cVp[isatom].ciV0.c2[iR] = (3.0*m-n-dV)*dx;
		cVp[isatom].ciV0.c3[iR] = (n-2.0*m)*dx2;
		/********************************************************/
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

// Cubic polynomial evaluation
__global__ void k_CubicPolyEval(sGP GP, scVp cVp, double * __restrict V0g, double * __restrict V1g)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if ((ix < cVp.bnx.n)&&(iy < cVp.bny.n)){
		ix += cVp.bnx.i;
		iy += cVp.bny.i;

		double Rx = ix*GP.dRx - cVp.x;
		double Ry = iy*GP.dRy - cVp.y;
		double R2 = Rx*Rx + Ry*Ry;
		if (R2 < cVp.Rmax2){
			double dx, dx2, dx3, V0, V1;
			if(R2 < cVp.Rmin2) R2 = cVp.Rmin2;

			ix = ix-(int)floor(ix*GP.dRx/GP.lx)*GP.nx;
			iy = iy-(int)floor(iy*GP.dRy/GP.ly)*GP.ny;

			ix = IsRS(ix, GP.nxh);
			iy = IsRS(iy, GP.nyh);
			int ixy = ix*GP.ny + iy;

			ix = unrolledBinarySearch128(R2, cVp.R2);

			dx = R2 - cVp.R2[ix]; dx2 = dx*dx; dx3 = dx2*dx;
			V0 = cVp.occ*(cVp.ciV0.c0[ix] + cVp.ciV0.c1[ix]*dx + cVp.ciV0.c2[ix]*dx2 + cVp.ciV0.c3[ix]*dx3);	
			V1 = cVp.occ*(cVp.ciV1.c0[ix] + cVp.ciV1.c1[ix]*dx + cVp.ciV1.c2[ix]*dx2 + cVp.ciV1.c3[ix]*dx3);
			if((ix<GP.nx)&&(iy<GP.ny)){
				V0g[ixy] += V0;
				V1g[ixy] += V1;
			}else{
				atomicAdd(&V0g[ixy], V0);
				atomicAdd(&V1g[ixy], V1);
			}
		}
	}
}

/***************************************************************************/
/***************************************************************************/
void cMT_Potential_GPU::SetPotPar(int PotParh){
	cudaMemcpyToSymbol(PotPar, &PotParh, cSizeofI);
}

void cMT_Potential_GPU::LinearProjAtomicPotentialGPU(dim3 grid, dim3 threads){
	switch(MT_MGP_CPU->MulOrder){
		case 1:
			k_LinearProjAtomicPotentialGPU<1><<<grid, threads>>>(Qz);
			break;
		case 2:
			k_LinearProjAtomicPotentialGPU<2><<<grid, threads>>>(Qz);
			break;
	}
}

void cMT_Potential_GPU::CubicPolyCoef(dim3 grid, dim3 threads){
	switch(MT_MGP_CPU->MulOrder){
		case 1:
			k_CubicPolyCoef<1><<<grid, threads>>>();
			break;
		case 2:
			k_CubicPolyCoef<2><<<grid, threads>>>();
			break;
	}	
}

void cMT_Potential_GPU::CubicPolyEval(int nsatom, double *&V0g, double *&V1g){
	dim3 BEval, TEval(thrnxny, thrnxny); 
	for(int isatom=0; isatom<nsatom; isatom++){
		BEval.x = (cVph[isatom].bny.n+thrnxny-1)/thrnxny; BEval.y = (cVph[isatom].bnx.n+thrnxny-1)/thrnxny;
		k_CubicPolyEval<<<BEval, TEval>>>(GP, cVph[isatom], V0g, V1g);
	}
}

void cMT_Potential_GPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	SetPotPar(0);

	cudaFreen(Qz.x);
	cudaFreen(Qz.w);

	f_sGP_Init(GP);

	f_scVp_Init(stncVp, cVph);
	for(int icVp=0; icVp<stncVp; icVp++){
		f_sciVn_cudaFree(cVph[icVp].ciV0);
		f_sciVn_cudaFree(cVph[icVp].ciV1);
	}

	delete [] MT_AtomTypes_GPU; MT_AtomTypes_GPU = 0;

	cudaFreen(V0);
	cudaFreen(V1);
	cudaFreen(V2);
}

cMT_Potential_GPU::cMT_Potential_GPU(){
	SetPotPar(0);

	Qz.x = 0;
	Qz.w = 0;

	f_sGP_Init(GP);

	f_scVp_Init(stncVp, cVph);
	for(int icVp=0; icVp<stncVp; icVp++){
		f_sciVn_cudaInit(cVph[icVp].ciV0);
		f_sciVn_cudaInit(cVph[icVp].ciV1);
	}

	MT_AtomTypes_GPU = 0;

	V0 = 0;
	V1 = 0;
	V2 = 0;
}

cMT_Potential_GPU::~cMT_Potential_GPU(){
	freeMemory();
}

int cMT_Potential_GPU::CheckGridLimits(int i, int n){
	return (i<0)?0:((i>=n)?n-1:i);
}

void cMT_Potential_GPU::getbn(sGP &GP, double x, double y, double Rmax, sbn &bnx, sbn &bny){
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
	bnx.i = ix0; bnx.n = (ix0==ixe)?0:ixe-ix0+1;
	bny.i = iy0; bny.n = (iy0==iye)?0:iye-iy0+1;
};

void cMT_Potential_GPU::setcVp(int iSlice, int iatom, int nsatom, dim3 &BPot, dim3 &TPot, dim3 &BCoef, dim3 &TCoef){
	int iZ;

	for(int i=0; i<nsatom; i++){
		iZ = Atoms[iatom+i].Z-1;
		cVph[i].x = Atoms[iatom+i].x;
		cVph[i].y = Atoms[iatom+i].y;
		cVph[i].z0 = Slice[iSlice].z0 - Atoms[iatom+i].z; 
		cVph[i].ze = Slice[iSlice].ze - Atoms[iatom+i].z;
		cVph[i].split = (cVph[i].z0<0)&&(0<cVph[i].ze);
		cVph[i].occ = Atoms[iatom+i].occ;
		cVph[i].Rmin2 = MT_AtomTypes_GPU[iZ].Rmin2;
		cVph[i].Rmax2 = MT_AtomTypes_GPU[iZ].Rmax2;
		cVph[i].cVr.cl = MT_AtomTypes_GPU[iZ].cVr.cl;
		cVph[i].cVr.cnl = MT_AtomTypes_GPU[iZ].cVr.cnl;
		cVph[i].R2 = MT_AtomTypes_GPU[iZ].R2;
		getbn(GP, cVph[i].x, cVph[i].y, MT_AtomTypes_GPU[iZ].Rmax, cVph[i].bnx, cVph[i].bny);
		if(MT_MGP_CPU->ApproxModel>1){
			cudaMemcpyAsync(cVph[i].ciV0.c0, MT_AtomTypes_GPU[iZ].ciVR.c0, stnR*cSizeofRD, cudaMemcpyDeviceToDevice);
			cudaMemcpyAsync(cVph[i].ciV0.c1, MT_AtomTypes_GPU[iZ].ciVR.c1, stnR*cSizeofRD, cudaMemcpyDeviceToDevice);
			cudaMemcpyAsync(cVph[i].ciV0.c2, MT_AtomTypes_GPU[iZ].ciVR.c2, stnR*cSizeofRD, cudaMemcpyDeviceToDevice);
			cudaMemcpyAsync(cVph[i].ciV0.c3, MT_AtomTypes_GPU[iZ].ciVR.c3, stnR*cSizeofRD, cudaMemcpyDeviceToDevice);
		}
	}
	cudaMemcpyToSymbolAsync(cVp, cVph, nsatom*cSizeofcVp, 0, cudaMemcpyHostToDevice);

	TPot.x = stnQz; TPot.y = 1; TPot.z = 1;
	BPot.x = stnR; BPot.y = nsatom; BPot.z = 1;

	TCoef.x = stnR; TCoef.y = 1; TCoef.z = 1;
	BCoef.x = nsatom; BCoef.y = 1; BCoef.z = 1;
}

void cMT_Potential_GPU::SetInputData(cMT_MGP_CPU *MT_MGP_CPU_io, int nAtomsM_i, double *AtomsM_i){
	freeMemory();

	f_sGP_SetInputData(MT_MGP_CPU, GP);
	cMT_Specimen_CPU::SetInputData(MT_MGP_CPU_io, nAtomsM_i, AtomsM_i, GP.dRmin);

	SetPotPar(MT_MGP_CPU->PotPar);			// Set Potential parameterization
	f_ReadQuadratureGPU(0, stnQz, Qz);		// TanhSinh

	for(int icVp=0; icVp<stncVp; icVp++){
		f_sciVn_cudaMalloc(stnR, cVph[icVp].ciV0);
		f_sciVn_cudaMalloc(stnR, cVph[icVp].ciV1);
	}

	MT_AtomTypes_GPU = new cMT_AtomTypes_GPU[nMT_AtomTypes];
	for (int i=0; i<nMT_AtomTypes; i++)
		MT_AtomTypes_GPU[i].SetAtomTypes(MT_AtomTypes_CPU[i]);

	cudaMalloc((void**)&V0, GP.nxy*cSizeofRD);
	cudaMalloc((void**)&V1, GP.nxy*cSizeofRD);
	cudaMalloc((void**)&V2, GP.nxy*cSizeofRD);
}

// Projected potential calculation: iSlice = slice position
void cMT_Potential_GPU::ProjectedPotential(int iSlice){	
	int iatom, nsatom;
	int iatom0 = Slice[iSlice].z0i_id;
	int iatome = Slice[iSlice].zei_id;
	dim3 BPot, TPot, BCoef, TCoef;

	f_Set_MD(GP, 0.0, V0, V1);

	if(iatome<iatom0) return;

	iatom = iatom0;
	while (iatom<=iatome){
		nsatom = MIN(stncVp, iatome-iatom+1);
		setcVp(iSlice, iatom, nsatom, BPot, TPot, BCoef, TCoef);
		if(MT_MGP_CPU->ApproxModel==1){
			LinearProjAtomicPotentialGPU(BPot, TPot);
			CubicPolyCoef(BCoef, TCoef);
		}
		CubicPolyEval(nsatom, V0, V1);
		iatom += nsatom;
	}

	if(MT_MGP_CPU->MulOrder==2)
		f_Scale_MD(GP, 1.0/get_dz(iSlice), V1);
}