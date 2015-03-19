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
#include "hMT_MGP_CPU.h"
#include "hMT_Specimen_CPU.h"
#include "hMT_AtomTypes_CPU.h"
#include "hMT_Potential_CPU.h"

#include "hmathCPU.h"

// Potential Evaluation (VR, dVRiR)
inline void cMT_Potential_CPU::Pot2Da(int PotPar, double &R, double *cl, double *cnl, double *icnl, double f, double &VR, double &dVRiR)
{
	double R2;
	double VR0, VR1, VR2, VR3, VR4, VR5;

	switch(PotPar)
	{
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
			VR0 = cl[0]*besskCPU(0,cnl[0]*R); VR1 = cl[1]*besskCPU(0,cnl[1]*R); VR2 = cl[2]*besskCPU(0,cnl[2]*R); 
			VR3 = cl[3]*exp(-cnl[3]*R2); VR4 = cl[4]*exp(-cnl[4]*R2); VR5 = cl[5]*exp(-cnl[5]*R2);
			VR = VR0 + VR1 + VR2 + VR3 + VR4 + VR5;
			dVRiR = cl[0]*cnl[0]*besskCPU(1,cnl[0]*R) + cl[1]*cnl[1]*besskCPU(1,cnl[1]*R) + cl[2]*cnl[2]*besskCPU(1,cnl[2]*R) + 2.0*R*(cnl[3]*VR3 + cnl[4]*VR4 + cnl[5]*VR5);
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
			VR0 = cl[0]*besskCPU(0,cnl[0]*R); VR1 = besskCPU(0,cnl[1]*R); VR2 = besskCPU(0,cnl[2]*R); 
			VR3 = besskCPU(0,cnl[3]*R); VR4 = besskCPU(0,cnl[4]*R);
			VR = 2.0*(cl[0]*icnl[0]*VR0 + cl[1]*icnl[1]*VR1 + cl[2]*icnl[2]*VR2 + cl[3]*icnl[3]*VR3 + cl[4]*icnl[4]*VR4);
			dVRiR = R*(cl[0]*cnl[0]*VR0 + cl[1]*cnl[1]*VR1 + cl[2]*cnl[2]*VR2 + cl[3]*cnl[3]*VR3 + cl[4]*cnl[4]*VR4);

			VR0 = besskCPU(1,cnl[0]*R); VR1 = besskCPU(1,cnl[1]*R); VR2 = besskCPU(1,cnl[2]*R); 
			VR3 = besskCPU(1,cnl[3]*R); VR4 = besskCPU(1,cnl[4]*R);
			VR += R*(cl[0]*VR0 + cl[1]*VR1 + cl[2]*VR2 + cl[3]*VR3 + cl[4]*VR4);
			dVRiR += 2.0*(cl[0]*VR0 + cl[1]*VR1 + cl[2]*VR2 + cl[3]*VR3 + cl[4]*VR4);
			VR = f*VR;
			dVRiR = -f*dVRiR/R;
			break;
	}
}

// Potential Evaluation (Vr, dVrir)
inline void cMT_Potential_CPU::Pot3Da(int PotPar, double &r, double *cl, double *cnl, double *icnl, double f, double &Vr, double &dVrir)
{
	double ir, r2;
	double Vr0, Vr1, Vr2, Vr3, Vr4, Vr5;

	switch(PotPar)
	{
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

//Binary search
inline int cMT_Potential_CPU::unrolledBinarySearch128(double x0, const double * __restrict x)
{
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

// Linear projected potential: V and zV
void cMT_Potential_CPU::getLinearProjAtomicPotential(int PotPar, sQ1 &Qz, scVp &cVp)
{	
	double V0, dV0, V1, dV1, V, dVir;
	double a, b, z, r, zm = 0.5*(cVp.ze+cVp.z0);
	double *cl, *cnl, icnl[6];
	int i, iR, iQz;

	for(i=0; i<6; i++)
		icnl[i] = (cnl[i]>0)?1.0/cnl[i]:0.0;

	a = (cVp.split)?-0.5*cVp.z0:0.5*(cVp.ze-cVp.z0); 
	b = (cVp.split)?0.5*cVp.z0:0.5*(cVp.ze+cVp.z0);
	for(iR=0; iR<stnR; iR++)
	{
		V0 = dV0 = V1 = dV1 = 0.0;
		for(iQz=0; iQz<stnQz; iQz++)
		{
			z = a*Qz.x[iQz] + b, r = sqrt(z*z + cVp.R2[iR]);
			Pot3Da(PotPar, r, cl, cnl, icnl, a*Qz.w[iQz], V, dVir);
			V0 += V; dV0 += dVir;
			V1 += (z-zm)*V; dV1 += (z-zm)*dVir;
		}
		cVp.ciV0.c0[iR] = V0;			// V0
		cVp.ciV0.c1[iR] = 0.5*dV0;		// dR2V0
		cVp.ciV1.c0[iR] = V1;			// V1
		cVp.ciV1.c1[iR] = 0.5*dV1;		// dR2V1
	}

	if(cVp.split)
	{
		a = b = 0.5*cVp.ze;
		for(iR=0; iR<stnR; iR++)
		{
			V0 = dV0 = V1 = dV1 = 0.0;
			for(iQz=0; iQz<stnQz; iQz++)
			{
				z = a*Qz.x[iQz] + b, r = sqrt(z*z + cVp.R2[iR]);
				Pot3Da(PotPar, r, cl, cnl, icnl, a*Qz.w[iQz], V, dVir);
				V0 += V; dV0 += dVir;
				V1 += (z-zm)*V; dV1 += (z-zm)*dVir;
			}
			cVp.ciV0.c0[iR] += V0;			// V0
			cVp.ciV0.c1[iR] += 0.5*dV0;		// dR2V0
			cVp.ciV1.c0[iR] += V1;			// V1
			cVp.ciV1.c1[iR] += 0.5*dV1;		// dR2V1
		}
	}
}

// Get Local interpolation coefficients
void cMT_Potential_CPU::getCubicPolyCoef(scVp &cVp)
{
	double dx, dx2;
	double V, Vn, dV, dVn, m, n;
	for(int iR=0; iR<stnR-1; iR++)
	{
		dx = 1.0/(cVp.R2[iR+1]-cVp.R2[iR]); dx2 = dx*dx;
		/********************************************************/
		V = cVp.ciV0.c0[iR]; Vn = cVp.ciV0.c0[iR+1];
		dV = cVp.ciV0.c1[iR]; dVn = cVp.ciV0.c1[iR+1];
		m = (Vn-V)*dx; n = dV+dVn;
		cVp.ciV0.c0[iR] = V-cVp.ciV0.c0[stnR-1];
		cVp.ciV0.c2[iR] = (3.0*m-n-dV)*dx;
		cVp.ciV0.c3[iR] = (n-2.0*m)*dx2;
		/********************************************************/
		V = cVp.ciV1.c0[iR]; Vn = cVp.ciV1.c0[iR+1];
		dV = cVp.ciV1.c1[iR]; dVn = cVp.ciV1.c1[iR+1];
		m = (Vn-V)*dx; n = dV+dVn;
		cVp.ciV1.c0[iR] = V-cVp.ciV1.c0[stnR-1];
		cVp.ciV1.c2[iR] = (3.0*m-n-dV)*dx;
		cVp.ciV1.c3[iR] = (n-2.0*m)*dx2;
	}
}

// Cubic polynomial evaluation
void cMT_Potential_CPU::evalCubicPoly(sGP &GP, scVp &cVp, double * __restrict V0g, double * __restrict V1g)
{
	int ix0, iy0, ix, iy, ixy;
	double Rx, Ry, R2;
	double dx, dx2, dx3;

	for(iy0=0; iy0<cVp.bny.n; iy0++)
	{
		iy = iy0 + cVp.bny.i;
		Ry = iy*GP.dRy - cVp.y;
		for(ix0=0; ix0<cVp.bnx.n; ix0++)
		{
			ix = ix0 + cVp.bnx.i;
			Rx = ix*GP.dRx - cVp.x;
			R2 = Rx*Rx + Ry*Ry;
			if(R2 < cVp.Rmax2)
			{
				if(R2 < cVp.Rmin2) R2 = cVp.Rmin2;

				ix = ix-(int)floor(ix*GP.dRx/GP.lx)*GP.nx;
				iy = iy-(int)floor(iy*GP.dRy/GP.ly)*GP.ny;

				ix = IsRS(ix, GP.nxh);
				iy = IsRS(iy, GP.nyh);
				ixy = iy*GP.nx + ix;

				ix = unrolledBinarySearch128(R2, cVp.R2);

				dx = R2 - cVp.R2[ix]; dx2 = dx*dx; dx3 = dx2*dx;
				V0g[ixy] += cVp.occ*(cVp.ciV0.c0[ix] + cVp.ciV0.c1[ix]*dx + cVp.ciV0.c2[ix]*dx2 + cVp.ciV0.c3[ix]*dx3);
				V1g[ixy] += cVp.occ*(cVp.ciV1.c0[ix] + cVp.ciV1.c1[ix]*dx + cVp.ciV1.c2[ix]*dx2 + cVp.ciV1.c3[ix]*dx3);	
			}
		}
	}
}

//get Effective Potential
void cMT_Potential_CPU::getV0(sGP &GP, eSlicePos SlicePo, double dz, double * __restrict V0_io, double * __restrict V1_io, double * __restrict V1o_io)
{
	int ix, iy, ixy;
	double V0, V1;

	switch(SlicePo)
	{
		case eSPFirst:		// initial Slice
			for(iy=0; iy<GP.ny; iy++)
				for(ix=0; ix<GP.nx; ix++)
				{
					ixy = iy*GP.nx+ix;
					V0 = V0_io[ixy], V1 = V1_io[ixy]/dz;
					V0_io[ixy] = V0 = V0-V1;
					V1o_io[ixy] = V1;
				}
			break;
		case eSPMedium:		// intermediate Slice
			for(iy=0; iy<GP.ny; iy++)
				for(ix=0; ix<GP.nx; ix++)
				{
					ixy = iy*GP.nx+ix;
					V0 = V0_io[ixy], V1 = V1_io[ixy]/dz;
					V0_io[ixy] = V0 = V0-V1+V1o_io[ixy];
					V1o_io[ixy] = V1;
				}
			break;
		case eSPLast:		// last Slice
			for(iy=0; iy<GP.ny; iy++)
				for(ix=0; ix<GP.nx; ix++)
				{
					V0_io[ixy] = V0 = V1o_io[ixy];
				}
			break;
	}
}

/***************************************************************************/
/***************************************************************************/

void cMT_Potential_CPU::freeMemory()
{
	if(IdCall==0) return;

	PotPar = 0;

	delete [] Qz.x; Qz.x = 0;
	delete [] Qz.w; Qz.w = 0;

	f_sGP_Init(GP);

	f_sciVn_Free_CPU(ciV0);
	f_sciVn_Free_CPU(ciV1);

	f_scVp_Init(cVp);

	delete [] V0; V0 = 0;
	delete [] V1; V1 = 0;
	delete [] V1o; V1o = 0;
}

cMT_Potential_CPU::cMT_Potential_CPU()
{
	IdCall = 0;

	PotPar = 0;

	Qz.x = 0;
	Qz.w = 0;

	f_sGP_Init(GP);

	f_sciVn_Init_CPU(ciV0);
	f_sciVn_Init_CPU(ciV1);

	f_scVp_Init(cVp);

	V0 = 0;
	V1 = 0;
	V1o = 0;
}

cMT_Potential_CPU::~cMT_Potential_CPU()
{
	freeMemory();
	IdCall = 0;
}

inline eSlicePos cMT_Potential_CPU::SlicePos(int iSlice, int nSlice)
{
	return (iSlice==0)?eSPFirst:(iSlice<nSlice)?eSPMedium:eSPLast;
}

inline int cMT_Potential_CPU::CheckGridLimits(int i, int n)
{
	return (i<0)?0:((i>=n)?n-1:i);
}

void cMT_Potential_CPU::getbn(sGP &GP, double x, double y, double Rmax, sbn &bnx, sbn &bny)
{
	double x0 = x-Rmax, xe = x+Rmax;
	double y0 = y-Rmax, ye = y+Rmax;
	int ix0 = (int)floor(x0/GP.dRx), ixe = (int)ceil(xe/GP.dRx);
	int iy0 = (int)floor(y0/GP.dRy), iye = (int)ceil(ye/GP.dRy);

	if(GP.PBC_xy==2)
	{
		ix0 = CheckGridLimits(ix0, GP.nx);
		ixe = CheckGridLimits(ixe, GP.nx);
		iy0 = CheckGridLimits(iy0, GP.ny);
		iye = CheckGridLimits(iye, GP.ny);
	}
	bnx.i = ix0; bnx.n = (ix0==ixe)?0:ixe-ix0+1;
	bny.i = iy0; bny.n = (iy0==iye)?0:iye-iy0+1;
};

void cMT_Potential_CPU::setcVp(int iSlice, int iatom, scVp &cVp)
{
	int iZ = Atoms[iatom].Z-1;
	cVp.x = Atoms[iatom].x;
	cVp.y = Atoms[iatom].y;
	cVp.z0 = Slice[iSlice].z0 - Atoms[iatom].z; 
	cVp.ze = Slice[iSlice].ze - Atoms[iatom].z;
	cVp.split = (cVp.z0<0)&&(0<cVp.ze);
	cVp.occ = Atoms[iatom].occ;
	cVp.Rmin2 = MT_AtomTypes_CPU[iZ].Rmin2;
	cVp.Rmax2 = MT_AtomTypes_CPU[iZ].Rmax2;
	cVp.cVr.cl = MT_AtomTypes_CPU[iZ].cVr.cl;
	cVp.cVr.cnl = MT_AtomTypes_CPU[iZ].cVr.cnl;
	cVp.R2 = MT_AtomTypes_CPU[iZ].R2;
	getbn(GP, cVp.x, cVp.y, MT_AtomTypes_CPU[iZ].Rmax, cVp.bnx, cVp.bny);
	if(MT_MGP_CPU->ApproxModel>1)
	{
		memcpy(cVp.ciV0.c0, MT_AtomTypes_CPU[iZ].ciVR.c0, stnR*cSizeofRD);
		memcpy(cVp.ciV0.c1, MT_AtomTypes_CPU[iZ].ciVR.c1, stnR*cSizeofRD);
		memcpy(cVp.ciV0.c2, MT_AtomTypes_CPU[iZ].ciVR.c2, stnR*cSizeofRD);
		memcpy(cVp.ciV0.c3, MT_AtomTypes_CPU[iZ].ciVR.c3, stnR*cSizeofRD);
	}
}

void cMT_Potential_CPU::addAtomicProjectedPotential(scVp &cVp, double *&V0g, double *&V1g)
{
	if(MT_MGP_CPU->ApproxModel==1)
	{
		getLinearProjAtomicPotential(PotPar, Qz, cVp);
		getCubicPolyCoef(cVp);
	}
	evalCubicPoly(GP, cVp, V0g, V1g);
}

void cMT_Potential_CPU::getV0(int iSlice, double *&V0, int typ)
{
	if(MT_MGP_CPU->MulOrder==1) return;

	 double dz = get_dz(iSlice);
	if(typ==1)
		getV0(GP, SlicePos(iSlice, nSlice), dz, V0, V1, V1o);
	else
		f_Scale_MD_CPU(GP, dz, V1);
}

void cMT_Potential_CPU::SetInputData(cMT_MGP_CPU *MT_MGP_CPU_io, int nAtomsM_i, double *AtomsM_i)
{
	freeMemory();
	IdCall++;

	f_sGP_SetInputData(MT_MGP_CPU_io, GP);
	cMT_Specimen_CPU::SetInputData(MT_MGP_CPU_io, nAtomsM_i, AtomsM_i, GP.dRmin);

	PotPar = MT_MGP_CPU->PotPar;			// Set Potential parameterization
	f_ReadQuadratureCPU(0, stnQz, Qz);		// TanhSinh

	f_sciVn_Malloc_CPU(stnR, ciV0);
	f_sciVn_Malloc_CPU(stnR, ciV1);

	cVp.ciV0.c0 = ciV0.c0;
	cVp.ciV0.c1 = ciV0.c1;
	cVp.ciV0.c2 = ciV0.c2;
	cVp.ciV0.c3 = ciV0.c3;

	cVp.ciV1.c0 = ciV1.c0;
	cVp.ciV1.c1 = ciV1.c1;
	cVp.ciV1.c2 = ciV1.c2;
	cVp.ciV1.c3 = ciV1.c3;

	V0 = new double[GP.nxy];
	V1 = new double[GP.nxy];
	V1o = new double[GP.nxy];
}

// Projected potential calculation: iSlice = Slice position
void cMT_Potential_CPU::ProjectedPotential(int iSlice, int typ)
{
	if(iSlice==nSlice)
	{
		getV0(iSlice, V0);
		return;
	}

	int iatom0 = Slice[iSlice].z0i_id;
	int iatome = Slice[iSlice].zei_id;
	f_Set_MD_CPU(GP, 0.0, V0, V1);

	if(iatome<iatom0) return;

	f_Set_MD_CPU(GP, 0.0, V0, V1);

	for(int iatom=iatom0; iatom<=iatome; iatom++)
	{
		setcVp(iSlice, iatom, cVp);
		addAtomicProjectedPotential(cVp, V0, V1);
	}

	getV0(iSlice, V0, typ);
}