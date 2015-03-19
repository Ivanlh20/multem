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

#ifndef hMT_Potential_CPU_H
#define hMT_Potential_CPU_H

#include "hConstTypes.h"
#include "hQuadrature.h"
#include "hMT_General_CPU.h"
#include "hMT_MGP_CPU.h"
#include "hMT_Specimen_CPU.h"
#include "hMT_AtomTypes_GPU.h"

class cMT_Potential_CPU: public cMT_Specimen_CPU{
	private:
		int IdCall;
		int PotPar;
		sQ1 Qz;
		sciVn ciV0;
		sciVn ciV1;
		scVp cVp;
		void Pot2Da(int PotPar, double &R, double *cl, double *cnl, double *icnl, double f, double &VR, double &dVRiR);
		void Pot3Da(int PotPar, double &r, double *cl, double *cnl, double *icnl, double f, double &Vr, double &dVrir);
		int unrolledBinarySearch128(double x0, const double * __restrict x);
		void getLinearProjAtomicPotential(int PotPar, sQ1 &Qz, scVp &cVp);
		void getCubicPolyCoef(scVp &cVp);
		void evalCubicPoly(sGP &GP, scVp &cVp, double * __restrict V0g, double * __restrict V1g);
		void getV0(sGP &GP, eSlicePos SlicePo, double dz, double * __restrict V0_io, double * __restrict V1_io, double * __restrict V1o_io);

		void addAtomicProjectedPotential(scVp &cVp, double *&V0g, double *&V1g);
		void getV0(int iSlice, double *&V0, int typ=1);
		inline int CheckGridLimits(int i, int n);
		void getbn(sGP &GP, double x, double y, double Rmax, sbn &bnx, sbn &bny);
		void setcVp(int iSlice, int iatom, scVp &cVp);
	public:
		sGP GP;									// xy-Grid properties

		double *V0;								// Zero moment of the potential
		double *V1;								// first moment of the potential
		double *V1o;							// Second moment of the potential

		void freeMemory();

		cMT_Potential_CPU();
		~cMT_Potential_CPU();

		inline eSlicePos SlicePos(int iSlice, int nSlice);
		void SetInputData(cMT_MGP_CPU *MT_MGP_CPU_io, int nAtomsM_i, double *AtomsM_i);
		void ProjectedPotential(int iSlice, int typ=1);
};

#endif