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

#include <cstring>
#include "hfxegTabData.h"
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int Z, ng, typ;
	double g[1024], g2[1024], fxg[1024], feg[1024];
	cfxegTabData fxegTabData;

	Z = (int)mxGetScalar(prhs[0]); 
	typ = (int)mxGetScalar(prhs[1]); 
	fxegTabData.ReadTabData(Z, typ, 1, ng, g, g2, fxg, feg);
	
	double *go,*fxgo,*fego;
	plhs[0] = mxCreateDoubleMatrix(ng, 1, mxREAL);
	go = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(ng, 1, mxREAL);
	fxgo = mxGetPr(plhs[1]);
 	plhs[2] = mxCreateDoubleMatrix(ng, 1, mxREAL);
	fego = mxGetPr(plhs[2]);

	memcpy(go, g, ng*sizeof(double));
	memcpy(fxgo, fxg, ng*sizeof(double));
	memcpy(fego, feg, ng*sizeof(double));
}