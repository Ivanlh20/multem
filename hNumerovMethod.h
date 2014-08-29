#ifndef hNumerovMethod_H
#define hNumerovMethod_H

#include "math.h"

class cNumerov{
	private:
		int n, l, nodes, ncross, icl, Ni;
		double eps, e, de, elw, eup, dx;
		double yicl, norm, ycusp, dfcusp;
		double dxf, def, y0, y1, lh2;
		int nr;
		double *r, *ro, *sqr, *r2, *Vr; 
		double *y, *f, *Vrt; 
	public:
		void ReadInputdata(int nri, double *ri, double *Vri);
		void SolveRadSchEq(int ni, double gamma, int Dim, double *aEner, double *aPsi);
		cNumerov();
		~cNumerov();
};

#endif