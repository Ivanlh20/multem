/*
 * This file is part of MULTEM.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
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
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef MATLAB_MEX_H
#define MATLAB_MEX_H

#include <algorithm>
#include <type_traits>
#include <cmath>

#include "types.cuh"
#include "traits.cuh"
#include "lin_alg_def.cuh"
#include "matlab_types.cuh"
#include "stream.cuh"
#include <mex.h>

using mt::rmatrix_r;
using mt::rmatrix_c;
using mt::r3d;

/**************************************************************************/
template <class TVector>
void mx_print(const TVector &vector)
{
	for(auto i = 0; i<vector.size(); i++)
	{
		mexPrintf("i = %d, value = %8.5f\n", i, double(vector[i]));
	}
}

inline int mx_get_MxN(const mxArray *mxB)
{
	return static_cast<int>(mxGetM(mxB)*mxGetN(mxB));
}

/**************************************************************************/
template <class T>
T mx_get_scalar(const mxArray *mxB)
{
	auto data = static_cast<int>(std::round(mxGetScalar(mxB)));
	return static_cast<T>(data);
}

template <>
int mx_get_scalar<int>(const mxArray *mxB)
{
	return static_cast<int>(std::round(mxGetScalar(mxB)));
}

template <>
double mx_get_scalar<double>(const mxArray *mxB)
{
	return mxGetScalar(mxB);
}

template <>
float mx_get_scalar<float>(const mxArray *mxB)
{
	return static_cast<float>(mxGetScalar(mxB));
}

/**************************************************************************/
template <class T>
T mx_get_matrix(const mxArray *mxB){};

template <>
rmatrix_r mx_get_matrix<rmatrix_r>(const mxArray *mxB)
{
	mt::rmatrix_r matrix;
	matrix.rows = static_cast<int>(mxGetM(mxB));
	matrix.cols = static_cast<int>(mxGetN(mxB));
	matrix.m_size = matrix.rows*matrix.cols;
	matrix.real = mxGetPr(mxB);

	return matrix;
}

template <>
rmatrix_c mx_get_matrix<rmatrix_c>(const mxArray *mxB)
{
	mt::rmatrix_c matrix;
	matrix.rows = static_cast<int>(mxGetM(mxB));
	matrix.cols = static_cast<int>(mxGetN(mxB));
	matrix.m_size = matrix.rows*matrix.cols;
	matrix.real = mxGetPr(mxB);
	matrix.imag = mxGetPi(mxB);

	return matrix;
}

/**************************************************************************/
bool mx_field_exits(const mxArray *mxB, const int &idx, const char *field_name)
{
	const mxArray *fPtr = mxGetField(mxB, idx, field_name);

	return (fPtr!=nullptr);
}

bool mx_field_exits(const mxArray *mxB, const char *field_name)
{
	return mx_field_exits(mxB, 0, field_name);
}

template <class T>
T mx_get_scalar_field(const mxArray *mxB, const int &idx, const char *field_name)
{
	const mxArray *fPtr = mxGetField(mxB, idx, field_name);

	if(fPtr==nullptr)
	{
		return T(0);
	}
	else
	{
		return mx_get_scalar<T>(fPtr);
	}
}

template <class T>
T mx_get_scalar_field(const mxArray *mxB, const char *field_name)
{
	return mx_get_scalar_field<T>(mxB, 0, field_name);
}

/**************************************************************************/
template <class T>
T mx_get_matrix_field(const mxArray *mxB, const int &idx, const char *field_name)
{
	mxArray *mx_matrix = mxGetField(mxB, idx, field_name);

	return mx_get_matrix<T>(mx_matrix);
}

template <class T>
T mx_get_matrix_field(const mxArray *mxB, const char *field_name)
{
	return mx_get_matrix_field<T>(mxB, 0, field_name);
}


template <class T>
r3d<T> mx_get_r3d_field(const mxArray *mxB, const char *field_name)
{
	auto r = mx_get_matrix_field<rmatrix_r>(mxB, field_name);
	return r3d<T>(r[0], r[1], r[2]);
}

/**************************************************************************/
template <class T>
T mx_create_matrix(const int &rows, const int &cols, mxArray *&mx_M);

template <>
rmatrix_r mx_create_matrix<rmatrix_r>(const int &rows, const int &cols, mxArray *&mx_M)
{
	mx_M = mxCreateDoubleMatrix(rows, cols, mxREAL);

	return mx_get_matrix<rmatrix_r>(mx_M);
}

template <>
rmatrix_c mx_create_matrix<rmatrix_c>(const int &rows, const int &cols, mxArray *&mx_M)
{
	mx_M = mxCreateDoubleMatrix(rows, cols, mxCOMPLEX);

	return mx_get_matrix<rmatrix_c>(mx_M);
}

template <class T>
T mx_create_scalar(mxArray *&mx_M)
{
	return mx_create_matrix<T>(1, 1, mx_M);
}

template <class T, class TGrid>
T mx_create_matrix(TGrid &grid_2d, mxArray *&mx_M)
{
	return mx_create_matrix<T>(grid_2d.ny, grid_2d.nx, mx_M);
}

/**************************************************************************/
template <class T>
inline T mx_create_matrix_field(mxArray *mx_struct, const int &idx, const char *field_name, const int &rows, const int &cols)
{
	mxArray *mxfield;
	T rmatrix = mx_create_matrix<T>(rows, cols, mxfield);
	mxSetField(mx_struct, idx, field_name, mxfield);
	return rmatrix;
}

template <class T>
inline T mx_create_matrix_field(mxArray *mx_struct, const char *field_name, const int &rows, const int &cols)
{
	return mx_create_matrix_field<T>(mx_struct, 0, field_name, rows, cols);
}

template <class T>
inline T mx_create_scalar_field(mxArray *mx_struct, const int &idx, const char *field_name)
{
	return mx_create_matrix_field<T>(mx_struct, idx, field_name, 1, 1);
}

template <class T>
inline T mx_create_scalar_field(mxArray *mx_struct, const char *field_name)
{
	return mx_create_scalar_field<T>(mx_struct, 0, field_name);
}

/**************************************************************************/
template <class T, class TVector>
inline void mx_create_set_matrix_field(mxArray *mx_struct, const int &idx, 
const char *field_name, const int &rows, const int &cols, TVector &field_value)
{
	mxArray *mxfield;
	T matrix = mx_create_matrix<T>(rows, cols, mxfield);
	mt::Stream<mt::e_host> stream(1);
	mt::copy_to_host(stream, field_value, matrix);
	mxSetField(mx_struct, idx, field_name, mxfield);
}

template <class T, class TVector>
inline void mx_create_set_matrix_field(mxArray *mx_struct, const char *field_name, 
const int &rows, const int &cols, TVector &field_value)
{
	mx_create_set_matrix_field<T>(mx_struct, 0, field_name, rows, cols, field_value);
}

template <class T>
inline void mx_create_set_scalar_field(mxArray *mx_struct, const int &idx, 
const char *field_name, const double &field_value)
{
	mxArray *mxfield;
	T matrix = mx_create_scalar<T>(mxfield);
	*(matrix.real) = field_value;
	mxSetField(mx_struct, idx, field_name, mxfield);
}

template <class T>
inline void mx_create_set_scalar_field(mxArray *mx_struct, const char *field_name, 
const double &field_value)
{
	 mx_create_set_scalar_field<T>(mx_struct, 0, field_name, field_value);
}

/**************************************************************************/
namespace mt
{
	System_Configuration read_system_conf(const mxArray *mx_input)
	{
		System_Configuration system_conf;

		if(mxIsStruct(mx_input) && mx_field_exits(mx_input, "precision"))
		{
			system_conf.device = mx_get_scalar_field<mt::eDevice>(mx_input, "device"); 
			system_conf.precision = mx_get_scalar_field<mt::ePrecision>(mx_input, "precision");
			system_conf.cpu_ncores = 1; 
			//system_conf.cpu_ncores = mx_get_scalar_field<int>(mx_input, "cpu_ncores"); 
			system_conf.cpu_nthread = mx_get_scalar_field<int>(mx_input, "cpu_nthread"); 
			system_conf.gpu_device = mx_get_scalar_field<int>(mx_input, "gpu_device");
			system_conf.gpu_nstream = 0; 
			//system_conf.gpu_nstream = mx_get_scalar_field<int>(mx_input, "gpu_nstream"); 
			system_conf.active = true;
		}
		else
		{
			system_conf.cpu_nthread = 4; 
			system_conf.active = false;
		}
		system_conf.validate_parameters();
		system_conf.set_device();

		return system_conf;
	}
}

#endif