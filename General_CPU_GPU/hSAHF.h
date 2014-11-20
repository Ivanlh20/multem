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

#ifndef hSAHF_H
#define hSAHF_H

#include "hRandGen.h"

#define nhmax 32
#define ngmax 2048

class cSAHF{
	private:
		cRandGen RandGenP, RandGenT;	
		ptrdiff_t rank, lwork, info, *iwork, *ipiv;
		double *work, rcond, *Sv;
		int Z;
		double Zk;

		double clmin[nhmax], clmax[nhmax], cnlmin[nhmax], cnlmax[nhmax], cnll[nhmax], cnld[nhmax];
		double srcnl[nhmax], icnl[nhmax], cnl2[nhmax], gamma[nhmax];
		double clfx[nhmax], cnlfx[nhmax], clPr[nhmax], cnlPr[nhmax], clVr[nhmax], cnlVr[nhmax];

		ptrdiff_t ng, nht, nhb, nhbc, nrb;	
		double g[ngmax], g2[ngmax], fxg[ngmax], feg[ngmax];
		double wxg[ngmax], fxgw[ngmax], weg[ngmax], fegw[ngmax], rb[ngmax];
		double *M, *Mt, *N, *Nt, *St, *Rm, *Rmt, *Fm;

		char *cTn;
		ptrdiff_t onei;
		double oned;

		void RandGenPReset();
		void RandGenPSetLimits();
		void RandGenPSetLimits(double *cnli, double dn);
		void RandGenPNext(double *cnlo);

		bool BoundaryConditions(double *cl, double *cnl);
		void GetLinealCoeff(double *cnli, double *clo, bool &bc);

		void SetInputData(int Zi, int typi, int nhti, int nhbi, int dgi, double *cnlmini, double *cnlmaxi);
		void NextCoeff(double *clo, double *cnlo);
		void CostFun(double *cli, double *cnli, double &ee, double &sigx, double &sige);
		void GetCfegVr(double *cli, double *cnli, double *clfego, double *cnlfego, double *clVro, double *cnlVro);

		void CopyData(int n, double eei, double *cli, double *cnli, double &eeo, double *clo, double *cnlo);
		void CopyData(int n, double eei, double sigxi, double sigei, double *cli, double *cnli, double &eeo, double &sigxo, double &sigeo, double *clo, double *cnlo);

		void GetIniPoint(int ni, double *clfego, double *cnlfego, double *clVro, double *cnlVro, double *eeo, double *sigxo, double *sigeo);
		void GetIniPoint(int ni, double *clo, double *cnlo, double &eeo, double &sigxo, double &sigeo, double &To);
		
		double dfactor(double p);
	public:
		void freeMemory();
		~cSAHF();	
		cSAHF();

		void SimulatedAnnealing(int Z, int nht, int nhb, int dgi, double *cnlmini, double *cnlmaxi, int N, 
		double *clfego, double *cnlfego, double *clVro, double *cnlVro, double *eeo, double *sigxo, double *sigeo);

		void SimulatedAnnealing(int Z, int nht, int nhb, int dgi, double *cnlmini, double *cnlmaxi, int Nt, int Ns, double rT1, 
		double rT2, double *clfego, double *cnlfego, double *clVro, double *cnlVro, double *eeo, double *sigxo, double *sigeo);
};
#endif