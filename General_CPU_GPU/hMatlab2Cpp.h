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

#ifndef hMatlab2Cpp_H
#define hMatlab2Cpp_H

#include <mex.h>

inline void CreateSetValue2mxField(mxArray *mxB, int p, const char *field_name, double field_value){
	mxArray *mxfield = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(mxfield) = field_value;
	mxSetField(mxB, p, field_name, mxfield);
}

inline void CreateSetValue2mxField(mxArray *mxB, int p, const char *field_name, int n, double *field_value){
	mxArray *mxfield = mxCreateDoubleMatrix(1, n, mxREAL);
	double *pfield = mxGetPr(mxfield); 
	memcpy(pfield, field_value, n*cSizeofRD);
	mxSetField(mxB, p, field_name, mxfield);
}

inline void CreateSetValue2mxField(mxArray *mxB, int p, const char *field_name, int nx, int ny, double *field_value){
	mxArray *mxfield = mxCreateDoubleMatrix(ny, nx, mxREAL);
	double *pfield = mxGetPr(mxfield); 
	memcpy(pfield, field_value, ny*nx*cSizeofRD);
	mxSetField(mxB, p, field_name, mxfield);
}

template <class Type>
inline Type ReadValuemxField(const mxArray *mxB, int p, const char *field_name, double factor=1){
	double val =factor*mxGetScalar(mxGetField(mxB, p, field_name));
	return (Type)val;
}

inline void ReadValuemxField(const mxArray *mxB, int p, const char *field_name, int n, double *field_value, double factor=1){
	double *pfield = mxGetPr(mxGetField(mxB, p, field_name));
	memcpy(field_value, pfield, n*cSizeofRD);
	for(int i=0; i<n; i++)
		field_value[i] *= factor;
}

#endif