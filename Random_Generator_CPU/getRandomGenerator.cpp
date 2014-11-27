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
#include "hRandGen.h"
#include <mex.h>

void Random(int tp, int n, double *r){
	cRandGen RandGen;	
	RandGen.reset();

		for (int i=0; i<n; i++){
			switch (tp){
				case 0:
					r[i] = RandGen.randu();
					break;
				case 1:
					r[i] = RandGen.randn();
					break;
			}
		}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ]){
	int tp, m, n;
	double *r;

	tp = (int)mxGetScalar(prhs[0]);
	m = (int)mxGetScalar(prhs[1]);
	n = (int)mxGetScalar(prhs[2]);

	plhs[0] = mxCreateDoubleMatrix((mwSize)m, (mwSize)n, mxREAL);
	r = mxGetPr(plhs[0]);
	Random(tp, m*n, r);
}