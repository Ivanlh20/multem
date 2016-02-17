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
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#include "types.cuh"
#include "matlab_types.cuh"
#include "traits.cuh"
#include "stream.cuh"

#include "input_output_superposition.cuh"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "superposition.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using multem::rmatrix_r;

template<class TInput_Superposition>
void read_input_superposition(const mxArray *mx_input_superposition, TInput_Superposition &input_superposition, bool full =true)
{
	using value_type_r = multem::Value_type<TInput_Superposition>;

	input_superposition.precision = mx_get_scalar_field<multem::ePrecision>(mx_input_superposition, "precision");
	input_superposition.device = mx_get_scalar_field<multem::eDevice>(mx_input_superposition, "device");
	input_superposition.cpu_nthread = mx_get_scalar_field<int>(mx_input_superposition, "cpu_nthread"); 
	input_superposition.gpu_device = mx_get_scalar_field<int>(mx_input_superposition, "gpu_device"); 
	input_superposition.gpu_nstream = mx_get_scalar_field<int>(mx_input_superposition, "gpu_nstream"); 
	
	input_superposition.ff_sigma = mx_get_scalar_field<value_type_r>(mx_input_superposition, "ff_sigma");

	bool bwl = false;
	bool pbc_xy = false; 																

	auto nx = mx_get_scalar_field<int>(mx_input_superposition, "nx");
	auto ny = mx_get_scalar_field<int>(mx_input_superposition, "ny");
	auto dx = mx_get_scalar_field<value_type_r>(mx_input_superposition, "dx");
	auto dy = mx_get_scalar_field<value_type_r>(mx_input_superposition, "dy");
	auto lx = nx*dx;
	auto ly = ny*dy;
	auto dz = 0.5; 				

	auto atoms = mx_get_matrix_field<rmatrix_r>(mx_input_superposition, "data");
	if(full)
	{
		input_superposition.atoms.set_Atoms(atoms.rows, atoms.real, lx, ly);
	}
	input_superposition.grid.set_input_data(nx, ny, lx, ly, dz, bwl, pbc_xy);
	input_superposition.validate_parameters();
 }

void set_superposition(const mxArray *mx_input_superposition, mxArray *&mx_output_superposition, multem::Output_Superposition_Matlab &output_superposition)
{
	multem::Input_Superposition<double> input_superposition;
	read_input_superposition(mx_input_superposition, input_superposition, false);
	output_superposition.set_input_data(&input_superposition);

	output_superposition.Im = mx_create_matrix<rmatrix_r>(input_superposition.grid.ny, input_superposition.grid.nx, mx_output_superposition);
}

template<class T, multem::eDevice dev>
void il_superposition(const mxArray *mxB, multem::Output_Superposition_Matlab &output_superposition)
{
	multem::Input_Superposition<T> input_superposition;
	read_input_superposition(mxB, input_superposition);

	multem::Superposition<T, dev> superposition;
	superposition.set_input_data(&input_superposition);
	superposition.run(output_superposition);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	multem::Output_Superposition_Matlab output_superposition;
	set_superposition(prhs[0], plhs[0], output_superposition);

	if(output_superposition.is_float_host())
	{
		il_superposition<float, multem::e_host>(prhs[0], output_superposition);
	}
	else if(output_superposition.is_double_host())
	{
		il_superposition<double, multem::e_host>(prhs[0], output_superposition);
	}
	if(output_superposition.is_float_device())
	{
		il_superposition<float, multem::e_device>(prhs[0], output_superposition);
	}
	else if(output_superposition.is_double_device())
	{
		il_superposition<double, multem::e_device>(prhs[0], output_superposition);
	}
}