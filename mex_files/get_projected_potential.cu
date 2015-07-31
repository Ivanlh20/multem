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

#include "types.hpp"
#include "traits.cuh"
#include "Stream.cuh"
#include "input_multislice.hpp"
#include "output_multislice.hpp"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "atom_data.hpp"
#include "potential.cuh"

#include <mex.h>
#include "mex_matlab.hpp"

using multem::rmatrix_r;

template<class TInput_Multislice>
void read_input_data(const mxArray *mx_input_multislice, TInput_Multislice &input_multislice, bool full=true)
{
	using value_type = multem::Value_type<TInput_Multislice>;

	input_multislice.precision = mx_get_scalar_field<multem::ePrecision>(mx_input_multislice, "precision");
	input_multislice.device = mx_get_scalar_field<multem::eDevice>(mx_input_multislice, "device");
	input_multislice.cpu_ncores = mx_get_scalar_field<int>(mx_input_multislice, "cpu_ncores"); 
	input_multislice.cpu_nthread = mx_get_scalar_field<int>(mx_input_multislice, "cpu_nthread"); 
	input_multislice.gpu_device = mx_get_scalar_field<int>(mx_input_multislice, "gpu_device"); 
	input_multislice.gpu_nstream = mx_get_scalar_field<int>(mx_input_multislice, "gpu_nstream"); 
	
	input_multislice.simulation_type = multem::eST_PPRS;
	input_multislice.phonon_model = mx_get_scalar_field<multem::ePhonon_Model>(mx_input_multislice, "phonon_model"); 
	input_multislice.interaction_model = mx_get_scalar_field<multem::eElec_Spec_Int_Model>(mx_input_multislice, "interaction_model");
	input_multislice.potential_slicing = mx_get_scalar_field<multem::ePotential_Slicing>(mx_input_multislice, "potential_slicing");
	input_multislice.potential_type = mx_get_scalar_field<multem::ePotential_Type>(mx_input_multislice, "potential_type");

	input_multislice.fp_dim.set(mx_get_scalar_field<int>(mx_input_multislice, "fp_dim"));
	input_multislice.fp_seed = mx_get_scalar_field<int>(mx_input_multislice, "fp_seed");
	input_multislice.fp_single_conf = mx_get_scalar_field<bool>(mx_input_multislice, "fp_single_conf");
	input_multislice.fp_nconf = mx_get_scalar_field<int>(mx_input_multislice, "fp_nconf");

	input_multislice.islice = mx_get_scalar_field<int>(mx_input_multislice, "islice")-1;

	bool bwl = false;
	bool pbc_xy = true; 																

	auto nx = mx_get_scalar_field<int>(mx_input_multislice, "nx");
	auto ny = mx_get_scalar_field<int>(mx_input_multislice, "ny");
	auto lx = mx_get_scalar_field<value_type>(mx_input_multislice, "lx");
	auto ly = mx_get_scalar_field<value_type>(mx_input_multislice, "ly");
	auto dz = mx_get_scalar_field<value_type>(mx_input_multislice, "dz"); 				

	auto atoms = mx_get_matrix_field<rmatrix_r>(mx_input_multislice, "atoms");
	if(full)
	{
		input_multislice.atoms.set_Atoms(atoms.rows, atoms.real, lx, ly, dz);
	}
	input_multislice.grid.set_input_data(nx, ny, lx, ly, dz, bwl, pbc_xy);
	input_multislice.validate_parameters();
 }

template<class TOutput_multislice>
void set_output_data(TOutput_multislice &output_multislice, mxArray *&mx_output_multislice)
{
	const char *field_names_output_multislice[] = {"V"};
	int number_of_fields_output_multislice = 1;
	mwSize dims_output_multislice[2] = {1, output_multislice.thickness.size()};

	mx_output_multislice = mxCreateStructArray(2, dims_output_multislice, number_of_fields_output_multislice, field_names_output_multislice);

	for(auto ithk=0; ithk<output_multislice.thickness.size(); ithk++)
	{
		output_multislice.V[ithk] = mx_create_matrix_field<rmatrix_r>(mx_output_multislice, ithk, "V", output_multislice.ny, output_multislice.nx);
	}
}

template<class T, multem::eDevice dev, class TOutput_multislice>
void get_projected_potential(const mxArray *mxB, TOutput_multislice &output_multislice)
{
	multem::Input_Multislice<T, dev> input_multislice;
	read_input_data(mxB, input_multislice);

	multem::Potential<T, dev> potential;
	multem::Stream<T, dev> stream;

	stream.resize(input_multislice.nstream);

	potential.set_input_data(&input_multislice, &stream);
	potential.move_atoms(input_multislice.fp_nconf);
	potential.projected_potential(input_multislice.islice);

	multem::to_host_shift(input_multislice.grid, potential.V0, output_multislice.V[0]);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	multem::Input_Multislice<double, multem::e_Host> input_multislice;
	read_input_data(prhs[0], input_multislice, false);

	multem::Output_Multislice<rmatrix_r, rmatrix_c> output_multislice;
	output_multislice.set_input_data(&input_multislice);

	set_output_data(output_multislice, plhs[0]);

	if(input_multislice.is_float_Host())
	{
		get_projected_potential<float, multem::e_Host>(prhs[0], output_multislice);
	}
	else if(input_multislice.is_double_Host())
	{
		get_projected_potential<double, multem::e_Host>(prhs[0], output_multislice);
	}
	if(input_multislice.is_float_Device())
	{
		get_projected_potential<float, multem::e_Device>(prhs[0], output_multislice);
	}
	else if(input_multislice.is_double_Device())
	{
		get_projected_potential<double, multem::e_Device>(prhs[0], output_multislice);
	}
}