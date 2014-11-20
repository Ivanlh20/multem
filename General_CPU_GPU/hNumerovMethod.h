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