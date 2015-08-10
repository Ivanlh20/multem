/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
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
#include "fft2.cuh"
#include "atom_data.hpp"
#include "input_multislice.hpp"
#include "output_multislice.hpp"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "transmission.cuh"

#include <mex.h>
#include "mex_matlab.hpp"

using multem::rmatrix_r;
using multem::rmatrix_c;

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
	
	input_multislice.simulation_type = multem::eST_TFRS;
	input_multislice.phonon_model = mx_get_scalar_field<multem::ePhonon_Model>(mx_input_multislice, "phonon_model"); 
	input_multislice.interaction_model = mx_get_scalar_field<multem::eElec_Spec_Int_Model>(mx_input_multislice, "interaction_model");
	input_multislice.potential_slicing = mx_get_scalar_field<multem::ePotential_Slicing>(mx_input_multislice, "potential_slicing");
	input_multislice.potential_type = mx_get_scalar_field<multem::ePotential_Type>(mx_input_multislice, "potential_type");

	input_multislice.fp_dim.set(mx_get_scalar_field<int>(mx_input_multislice, "fp_dim"));
	input_multislice.fp_seed = mx_get_scalar_field<int>(mx_input_multislice, "fp_seed");
	input_multislice.fp_single_conf = true;
	input_multislice.fp_nconf = mx_get_scalar_field<int>(mx_input_multislice, "fp_nconf");
	
	input_multislice.islice = mx_get_scalar_field<int>(mx_input_multislice, "islice")-1;

	input_multislice.E_0 = mx_get_scalar_field<value_type>(mx_input_multislice, "E_0");
	input_multislice.theta = mx_get_scalar_field<value_type>(mx_input_multislice, "theta")*multem::c_deg_2_rad;
	input_multislice.phi = mx_get_scalar_field<value_type>(mx_input_multislice, "phi")*multem::c_deg_2_rad;

	bool bwl = mx_get_scalar_field<bool>(mx_input_multislice, "bwl");
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
void set_output_data(const mxArray *mx_input_multislice, mxArray *&mx_output_multislice, TOutput_multislice &output_multislice)
{
	multem::Input_Multislice<double, multem::e_Host> input_multislice;
	read_input_data(mx_input_multislice, input_multislice, false);
	output_multislice.set_input_data(&input_multislice);

	const char *field_names_output_multislice[] = {"dx", "dy", "x", "y", "thickness", "trans"};
	int number_of_fields_output_multislice = 6;
	mwSize dims_output_multislice[2] = {1, 1};

	mx_output_multislice = mxCreateStructArray(2, dims_output_multislice, number_of_fields_output_multislice, field_names_output_multislice);

	mx_create_set_scalar_field<rmatrix_r>(mx_output_multislice, 0, "dx", output_multislice.dx);
	mx_create_set_scalar_field<rmatrix_r>(mx_output_multislice, 0, "dy", output_multislice.dy);
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "x", 1, output_multislice.x.size(), output_multislice.x.data());
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "y", 1, output_multislice.y.size(), output_multislice.y.data());
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "thickness", 1, output_multislice.thickness.size(), output_multislice.thickness.data());
	output_multislice.trans[0] = mx_create_matrix_field<rmatrix_c>(mx_output_multislice, "trans", output_multislice.ny, output_multislice.nx);
}

template<class T, multem::eDevice dev, class TOutput_multislice>
void get_transmission_function(const mxArray *mxB, TOutput_multislice &output_multislice)
{
	multem::Input_Multislice<T, dev> input_multislice;
	read_input_data(mxB, input_multislice);

	multem::Stream<T, dev> stream;
	multem::FFT2<T, dev> fft2;
	multem::Transmission<T, dev> transmission;

	stream.resize(input_multislice.nstream);
	fft2.create_plan(input_multislice.grid.ny, input_multislice.grid.nx, input_multislice.nstream);
	transmission.set_input_data(&input_multislice, &stream, &fft2);

	transmission.move_atoms(input_multislice.fp_nconf);
	transmission.trans(input_multislice.islice);

	multem::copy_to_host(input_multislice.grid, transmission.trans_0, output_multislice.trans[0]);
	multem::fft2_shift(input_multislice.grid, output_multislice.trans[0]);

	fft2.cleanup();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	multem::Output_Multislice<rmatrix_r, rmatrix_c> output_multislice;
	set_output_data(prhs[0], plhs[0], output_multislice);

	if(output_multislice.is_float_Host())
	{
		get_transmission_function<float, multem::e_Host>(prhs[0], output_multislice);
	}
	else if(output_multislice.is_double_Host())
	{
		get_transmission_function<double, multem::e_Host>(prhs[0], output_multislice);
	}
	if(output_multislice.is_float_Device())
	{
		get_transmission_function<float, multem::e_Device>(prhs[0], output_multislice);
	}
	else if(output_multislice.is_double_Device())
	{
		get_transmission_function<double, multem::e_Device>(prhs[0], output_multislice);
	}
}