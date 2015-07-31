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

#include "math.cuh"
#include "types.hpp"
#include "traits.cuh"
#include "Stream.cuh"
#include "FFt2.cuh"
#include "input_multislice.hpp"
#include "output_multislice.hpp"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "atom_data.hpp"
#include "wave_function.cuh"

#include <mex.h>
#include "mex_matlab.hpp"

using multem::rmatrix_r;
using multem::rmatrix_c;
using multem::Vector;

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

	input_multislice.simulation_type = mx_get_scalar_field<multem::eSimulation_Type>(mx_input_multislice, "simulation_type"); 
	input_multislice.phonon_model = mx_get_scalar_field<multem::ePhonon_Model>(mx_input_multislice, "phonon_model"); 
	input_multislice.interaction_model = mx_get_scalar_field<multem::eElec_Spec_Int_Model>(mx_input_multislice, "interaction_model");
	input_multislice.potential_slicing = mx_get_scalar_field<multem::ePotential_Slicing>(mx_input_multislice, "potential_slicing");
	input_multislice.potential_type = mx_get_scalar_field<multem::ePotential_Type>(mx_input_multislice, "potential_type");

	input_multislice.fp_dim.set(mx_get_scalar_field<int>(mx_input_multislice, "fp_dim"));
	input_multislice.fp_seed = mx_get_scalar_field<int>(mx_input_multislice, "fp_seed");
	input_multislice.fp_single_conf = mx_get_scalar_field<bool>(mx_input_multislice, "fp_single_conf");
	input_multislice.fp_nconf = mx_get_scalar_field<int>(mx_input_multislice, "fp_nconf");

	input_multislice.zero_defocus_type = mx_get_scalar_field<multem::eZero_Defocus_Type>(mx_input_multislice, "zero_defocus_type");
	input_multislice.zero_defocus_plane = mx_get_scalar_field<value_type>(mx_input_multislice, "zero_defocus_plane");
	
	input_multislice.thickness_type = mx_get_scalar_field<multem::eThickness_Type>(mx_input_multislice, "thickness_type");
	if(input_multislice.thickness_type != multem::eTT_Whole_Specimen && full)
	{
		auto thickness = mx_get_matrix_field<rmatrix_r>(mx_input_multislice, "thickness");
		input_multislice.thickness.resize(thickness.size);
		std::copy(thickness.real, thickness.real + thickness.size, input_multislice.thickness.begin());
	}

	input_multislice.input_wave_type = mx_get_scalar_field<multem::eInput_Wave_Type>(mx_input_multislice, "input_wave_type");
	if(input_multislice.is_user_define_wave() && full)
	{
		auto psi_0 = mx_get_matrix_field<rmatrix_c>(mx_input_multislice, "psi_0");
		input_multislice.psi_0.resize(psi_0.size);
		multem::rmatrix_c_to_complex(psi_0, input_multislice.psi_0);
		multem::fft2_shift(input_multislice.grid, input_multislice.psi_0);
	}

	input_multislice.coherent_contribution = true;

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
	input_multislice.lens.sf = mx_get_scalar_field<value_type>(mx_input_multislice, "lens_sf"); 								// defocus spread(Angstrom)
	input_multislice.lens.nsf = mx_get_scalar_field<int>(mx_input_multislice, "lens_nsf"); 										// Number of defocus sampling point
	input_multislice.lens.beta = mx_get_scalar_field<value_type>(mx_input_multislice, "lens_beta")*multem::c_mrad_2_rad; 		// semi-convergence angle(mrad-->rad)
	input_multislice.lens.nbeta = mx_get_scalar_field<int>(mx_input_multislice, "lens_nbeta"); 									// half number sampling points
	input_multislice.lens.set_input_data(input_multislice.E_0, input_multislice.grid);

	if (input_multislice.is_EWFS())
	{
		input_multislice.ew_fr.space = multem::eS_Reciprocal;
		input_multislice.ew_fr.convergent_beam = mx_get_scalar_field<bool>(mx_input_multislice, "ewfs_convergent_beam");
		input_multislice.ew_fr.x0 = mx_get_scalar_field<value_type>(mx_input_multislice, "ewfs_x0");
		input_multislice.ew_fr.y0 = mx_get_scalar_field<value_type>(mx_input_multislice, "ewfs_y0");
	}
	else if (input_multislice.is_EWRS())
	{
		input_multislice.ew_fr.space = multem::eS_Real;
		input_multislice.ew_fr.convergent_beam = mx_get_scalar_field<bool>(mx_input_multislice, "ewrs_convergent_beam");
		input_multislice.ew_fr.x0 = mx_get_scalar_field<value_type>(mx_input_multislice, "ewrs_x0");
		input_multislice.ew_fr.y0 = mx_get_scalar_field<value_type>(mx_input_multislice, "ewrs_y0");
	}

	input_multislice.validate_parameters();
}

template<class TInput_Multislice, class TOutput_multislice>
void set_output_data(const TInput_Multislice &input_multislice, mxArray *&mx_plhs0, TOutput_multislice &output_multislice)
{
	const char *field_names_psi[] = {"psi"};
	int number_of_fields_psi = 1;
	mwSize dims_psi[2] = {1, input_multislice.thickness.size()};

	output_multislice.resize(input_multislice.thickness.size());
	mx_plhs0 = mxCreateStructArray(2, dims_psi, number_of_fields_psi, field_names_psi);

	for(auto ithk=0; ithk<input_multislice.thickness.size(); ithk++)
	{
		output_multislice[ithk] = mx_create_matrix_field<rmatrix_c>(mx_plhs0, ithk, "psi", input_multislice.grid.ny, input_multislice.grid.nx);
	}
}

template<class T, multem::eDevice dev>
void get_wave_function(const mxArray *mxB, Vector<rmatrix_c, multem::e_Host> &wave)
{
	multem::Input_Multislice<T, dev> input_multislice;
	read_input_data(mxB, input_multislice);

	multem::Stream<T, dev> stream;
	multem::FFT2<T, dev> fft2;
	multem::Wave_Function<T, dev> wave_function;

	stream.resize(input_multislice.nstream);
	fft2.create_plan(input_multislice.grid.ny, input_multislice.grid.nx, input_multislice.nstream);

	wave_function.set_input_data(&input_multislice, &stream, &fft2);
	wave_function.move_atoms(input_multislice.fp_iconf);
	wave_function.psi_0();
	wave_function.psi(input_multislice.ew_fr.space, wave);

	fft2.cleanup();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	multem::Input_Multislice<double, multem::e_Host> input_multislice;
	read_input_data(prhs[0], input_multislice);

	Vector<rmatrix_c, multem::e_Host> wave;
	set_output_data(input_multislice, plhs[0], wave);

	if(input_multislice.is_float_Host())
	{
		get_wave_function<float, multem::e_Host>(prhs[0], wave);
	}
	else if(input_multislice.is_double_Host())
	{
		get_wave_function<double, multem::e_Host>(prhs[0], wave);
	}
	if(input_multislice.is_float_Device())
	{
		get_wave_function<float, multem::e_Device>(prhs[0], wave);
	}
	else if(input_multislice.is_double_Device())
	{
		get_wave_function<double, multem::e_Device>(prhs[0], wave);
	}
}