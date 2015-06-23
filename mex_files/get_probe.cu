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

#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"

#include "Probe.cuh"

#include <mex.h>
#include "matlab2cpp.hpp"

using multem::m_matrix_r;
using multem::m_matrix_c;

template<class TInput_Multislice>
void read_input_data(const mxArray *mx_input_multislice, TInput_Multislice &input_multislice)
{
	using value_type = multem::traits::Value_type<TInput_Multislice>;

	input_multislice.precision = mx_get_scalar_field<multem::ePrecision>(mx_input_multislice, "precision");
	input_multislice.device = mx_get_scalar_field<multem::eDevice>(mx_input_multislice, "device"); 
	input_multislice.cpu_ncores = mx_get_scalar_field<int>(mx_input_multislice, "cpu_ncores"); 
	input_multislice.cpu_nthread = mx_get_scalar_field<int>(mx_input_multislice, "cpu_nthread"); 
	input_multislice.gpu_device = mx_get_scalar_field<int>(mx_input_multislice, "gpu_device"); 
	input_multislice.gpu_nstream = mx_get_scalar_field<int>(mx_input_multislice, "gpu_nstream"); 

	input_multislice.simulation_type = multem::eST_ProbeRS; 

	input_multislice.E_0 = mx_get_scalar_field<value_type>(mx_input_multislice, "E_0");

	auto nx = mx_get_scalar_field<int>(mx_input_multislice, "nx");
	auto ny = mx_get_scalar_field<int>(mx_input_multislice, "ny");
	auto lx = mx_get_scalar_field<value_type>(mx_input_multislice, "lx");
	auto ly = mx_get_scalar_field<value_type>(mx_input_multislice, "ly");
	bool bwl = true;
	bool pbc_xy = true;
	auto dz = value_type(0.25);

	input_multislice.grid.set_input_data(nx, ny, lx, ly, dz, bwl, pbc_xy);

	input_multislice.lens.m = mx_get_scalar_field<int>(mx_input_multislice, "lens_m"); 											// momentum of the vortex
	input_multislice.lens.f = mx_get_scalar_field<value_type>(mx_input_multislice, "lens_f"); 									// defocus(Angstrom)
	input_multislice.lens.Cs3 = mx_get_scalar_field<value_type>(mx_input_multislice, "lens_Cs3")*multem::c_mm_2_Ags; 			// spherical aberration(mm-->Angstrom)
	input_multislice.lens.Cs5 = mx_get_scalar_field<value_type>(mx_input_multislice, "lens_Cs5")*multem::c_mm_2_Ags; 			// spherical aberration(mm-->Angstrom)
	input_multislice.lens.mfa2 = mx_get_scalar_field<value_type>(mx_input_multislice, "lens_mfa2"); 							// magnitude 2-fold astigmatism(Angstrom)
	input_multislice.lens.afa2 = mx_get_scalar_field<value_type>(mx_input_multislice, "lens_afa2")*multem::c_deg_2_rad; 		// angle 2-fold astigmatism(degrees-->rad)
	input_multislice.lens.mfa3 = mx_get_scalar_field<value_type>(mx_input_multislice, "lens_mfa3"); 							// magnitude 3-fold astigmatism(Angstrom)
	input_multislice.lens.afa3 = mx_get_scalar_field<value_type>(mx_input_multislice, "lens_afa3")*multem::c_deg_2_rad; 		// angle 3-fold astigmatism(degrees-->rad)
	input_multislice.lens.aobjl = mx_get_scalar_field<value_type>(mx_input_multislice, "lens_aobjl")*multem::c_mrad_2_rad; 		// lower objective aperture(mrad-->rad)
	input_multislice.lens.aobju = mx_get_scalar_field<value_type>(mx_input_multislice, "lens_aobju")*multem::c_mrad_2_rad; 		// upper objective aperture(mrad-->rad)

	input_multislice.conv_beam_wave_x = mx_get_scalar_field<value_type>(mx_input_multislice, "conv_beam_wave_x");
	input_multislice.conv_beam_wave_y = mx_get_scalar_field<value_type>(mx_input_multislice, "conv_beam_wave_y");

	input_multislice.validate_parameters();
}

template<class T, multem::eDevice dev>
void get_probe(const mxArray *mxB, m_matrix_c &host_probe)
{
	multem::Input_Multislice<T, dev> input_multislice;
	read_input_data(mxB, input_multislice);

	multem::Probe<T, dev> probe;

	probe.set_input_data(&input_multislice);
	probe.get(multem::eS_Real, host_probe);

	probe.cleanup();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	multem::Input_Multislice<double, multem::Host> input_multislice;
	read_input_data(prhs[0], input_multislice);

	auto probe = mx_create_matrix<m_matrix_c>(input_multislice.grid, plhs[0]);

	if(input_multislice.is_float_Host())
	{
		get_probe<float, multem::Host>(prhs[0], probe);
	}
	else if(input_multislice.is_double_Host())
	{
		get_probe<double, multem::Host>(prhs[0], probe);
	}
	if(input_multislice.is_float_Device())
	{
		get_probe<float, multem::Device>(prhs[0], probe);
	}
	else if(input_multislice.is_double_Device())
	{
		get_probe<double, multem::Device>(prhs[0], probe);
	}
}