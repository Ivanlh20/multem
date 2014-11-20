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

#include "hConstTypes.h"
#include "hQuadrature.h"
#include <mex.h>
	
void mexFunction(int nlhs,mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int QuadType, nQuad;
	double *g, *feg, *dfeg;
	int m, n, ng;

	QuadType = (int)mxGetScalar(prhs[0]);
	nQuad = (int)mxGetScalar(prhs[1]);

	/******************************************************************************/
	sQ1 Q1;
	plhs[0] = mxCreateDoubleMatrix(nQuad, 1, mxREAL);
	Q1.x = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(nQuad, 1, mxREAL);
	Q1.w = mxGetPr(plhs[1]);

	cQuadrature Quad;
	Quad.ReadQuadrature(QuadType, nQuad, Q1);
}