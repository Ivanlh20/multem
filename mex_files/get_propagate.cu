/**
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
#include "input_multislice.hpp"
#include "atom_data.hpp"
#include "propagator.cuh"

#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"

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

	input_multislice.simulation_type = multem::eST_HRTEM;

	input_multislice.input_wave_type = mx_get_scalar_field<multem::eInput_Wave_Type>(mx_input_multislice, "input_wave_type");
	if(input_multislice.is_user_define_wave() && full)
	{
		auto psi_0 = mx_get_matrix_field<rmatrix_c>(mx_input_multislice, "psi_0");
		multem::assign(psi_0, input_multislice.psi_0);
		multem::fft2_shift(input_multislice.grid, input_multislice.psi_0);
	}

	input_multislice.E_0 = mx_get_scalar_field<value_type>(mx_input_multislice, "E_0");
	input_multislice.theta = mx_get_scalar_field<value_type>(mx_input_multislice, "theta")*multem::c_deg_2_rad;
	input_multislice.phi = mx_get_scalar_field<value_type>(mx_input_multislice, "phi")*multem::c_deg_2_rad;

	bool bwl = true;
	bool pbc_xy = true;

	auto nx = mx_get_scalar_field<int>(mx_input_multislice, "nx");
	auto ny = mx_get_scalar_field<int>(mx_input_multislice, "ny");
	auto lx = mx_get_scalar_field<value_type>(mx_input_multislice, "lx");
	auto ly = mx_get_scalar_field<value_type>(mx_input_multislice, "ly");
	value_type dz = 0.25;

	input_multislice.grid.set_input_data(nx, ny, lx, ly, dz, bwl, pbc_xy);

	input_multislice.lens.f = mx_get_scalar_field<value_type>(mx_input_multislice, "lens_f"); 									// defocus(Angstrom)
	input_multislice.lens.set_input_data(input_multislice.E_0, input_multislice.grid);

	input_multislice.validate_parameters();
 }

template<class T, multem::eDevice dev>
void get_propagate(const mxArray *mxB, rmatrix_c &psi_host)
{
	multem::Input_Multislice<T, dev> input_multislice;
	read_input_data(mxB, input_multislice);

	multem::Stream<T, dev> stream;
	multem::FFT2<T, dev> fft2;
	multem::Propagator<T, dev> propagator;
	multem::Vector<complex<T>, dev> psi(input_multislice.grid.nxy());

	stream.resize(input_multislice.nstream);
	fft2.create_plan(input_multislice.grid.ny, input_multislice.grid.nx);
	propagator.set_input_data(&input_multislice, &fft2);

	T gx_0 = input_multislice.gx_0();
	T gy_0 = input_multislice.gy_0();

	propagator.propagate(multem::eS_Real, gx_0, gy_0, input_multislice.lens.f, input_multislice.psi_0, psi);

	multem::to_host_shift(input_multislice.grid, psi, psi_host);

	fft2.cleanup();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	multem::Input_Multislice<double, multem::e_Host> input_multislice;
	read_input_data(prhs[0], input_multislice, false);

	auto psi = mx_create_matrix<rmatrix_c>(input_multislice.grid, plhs[0]);

	if(input_multislice.is_float_Host())
	{
		get_propagate<float, multem::e_Host>(prhs[0], psi);
	}
	else if(input_multislice.is_double_Host())
	{
		get_propagate<double, multem::e_Host>(prhs[0], psi);
	}
	if(input_multislice.is_float_Device())
	{
		get_propagate<float, multem::e_Device>(prhs[0], psi);
	}
	else if(input_multislice.is_double_Device())
	{
		get_propagate<double, multem::e_Device>(prhs[0], psi);
	}
}