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
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#include "types.cuh"
#include "matlab_types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "fft2.cuh"
#include "atom_data.hpp"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "wave_function.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using multem::rmatrix_r;
using multem::rmatrix_c;

template<class TInput_Multislice>
void read_input_multislice(const mxArray *mx_input_multislice, TInput_Multislice &input_multislice, bool full =true)
{
	using value_type_r = multem::Value_type<TInput_Multislice>;

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
	input_multislice.fp_single_conf = true;
	input_multislice.fp_nconf = mx_get_scalar_field<int>(mx_input_multislice, "fp_nconf");

	input_multislice.tm_active = mx_get_scalar_field<bool>(mx_input_multislice, "tm_active");
	input_multislice.tm_theta = mx_get_scalar_field<value_type_r>(mx_input_multislice, "tm_theta")*multem::c_deg_2_rad;
	input_multislice.tm_u0 = mx_get_r3d_field<value_type_r>(mx_input_multislice, "tm_u0");
	input_multislice.tm_rot_point_type = mx_get_scalar_field<multem::eRot_Point_Type>(mx_input_multislice, "tm_rot_point_type");
	input_multislice.tm_p0 = mx_get_r3d_field<value_type_r>(mx_input_multislice, "tm_p0");

	input_multislice.thickness_type = mx_get_scalar_field<multem::eThickness_Type>(mx_input_multislice, "thickness_type");
	if(!input_multislice.is_whole_specimen() && full)
	{
		auto thickness = mx_get_matrix_field<rmatrix_r>(mx_input_multislice, "thickness");
		input_multislice.thickness.resize(thickness.m_size);
		std::copy(thickness.real, thickness.real + thickness.m_size, input_multislice.thickness.begin());
	}

	input_multislice.coherent_contribution = true;

	input_multislice.E_0 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "E_0");
	input_multislice.theta = mx_get_scalar_field<value_type_r>(mx_input_multislice, "theta")*multem::c_deg_2_rad;
	input_multislice.phi = mx_get_scalar_field<value_type_r>(mx_input_multislice, "phi")*multem::c_deg_2_rad;

	bool bwl = mx_get_scalar_field<bool>(mx_input_multislice, "bwl");
	bool pbc_xy = true;

	auto nx = mx_get_scalar_field<int>(mx_input_multislice, "nx");
	auto ny = mx_get_scalar_field<int>(mx_input_multislice, "ny");
	auto lx = mx_get_scalar_field<value_type_r>(mx_input_multislice, "lx");
	auto ly = mx_get_scalar_field<value_type_r>(mx_input_multislice, "ly");
	auto dz = mx_get_scalar_field<value_type_r>(mx_input_multislice, "dz"); 				

	auto atoms = mx_get_matrix_field<rmatrix_r>(mx_input_multislice, "atoms");
	if(full)
	{
		input_multislice.atoms.set_Atoms(atoms.rows, atoms.cols, atoms.real, lx, ly);
	}
	input_multislice.grid.set_input_data(nx, ny, lx, ly, dz, bwl, pbc_xy);

	/********************************* Incident wave ********************************/
	auto iw_type = mx_get_scalar_field<multem::eIncident_Wave_Type>(mx_input_multislice, "iw_type");
	input_multislice.set_incident_wave_type(iw_type);

	if(input_multislice.is_user_define_wave() && full)
	{
		auto iw_psi = mx_get_matrix_field<rmatrix_c>(mx_input_multislice, "iw_psi");
		multem::Stream<multem::e_host> stream(input_multislice.cpu_nthread);
		multem::assign(stream, iw_psi, input_multislice.iw_psi);
	}
	input_multislice.iw_x = mx_get_scalar_field<value_type_r>(mx_input_multislice, "iw_x");
	input_multislice.iw_y = mx_get_scalar_field<value_type_r>(mx_input_multislice, "iw_y");

	/****************************** Condenser lens ********************************/
	input_multislice.cond_lens.m = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_m"); 											// momentum of the vortex
	input_multislice.cond_lens.f = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_f"); 									// defocus(Angstrom)
	input_multislice.cond_lens.Cs3 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_Cs3")*multem::c_mm_2_Ags; 			// third order spherical aberration(mm-->Angstrom)
	input_multislice.cond_lens.Cs5 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_Cs5")*multem::c_mm_2_Ags; 			// fifth order aberration(mm-->Angstrom)
	input_multislice.cond_lens.mfa2 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_mfa2"); 							// magnitude 2-fold astigmatism(Angstrom)
	input_multislice.cond_lens.afa2 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_afa2")*multem::c_deg_2_rad; 		// angle 2-fold astigmatism(degrees-->rad)
	input_multislice.cond_lens.mfa3 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_mfa3"); 							// magnitude 3-fold astigmatism(Angstrom)
	input_multislice.cond_lens.afa3 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_afa3")*multem::c_deg_2_rad; 		// angle 3-fold astigmatism(degrees-->rad)
	input_multislice.cond_lens.inner_aper_ang = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_inner_aper_ang")*multem::c_mrad_2_rad; 		// inner aperture(mrad-->rad)
	input_multislice.cond_lens.outer_aper_ang = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_outer_aper_ang")*multem::c_mrad_2_rad; 		// outer aperture(mrad-->rad)
	input_multislice.cond_lens.sf = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_sf"); 								// defocus spread(Angstrom)
	input_multislice.cond_lens.nsf = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_nsf"); 										// Number of integration steps for the defocus Spread
	input_multislice.cond_lens.beta = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_beta")*multem::c_mrad_2_rad; 		// divergence semi-angle(mrad-->rad)
	input_multislice.cond_lens.nbeta = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_nbeta");									// Number of integration steps for the divergence semi-angle
	input_multislice.cond_lens.zero_defocus_type = mx_get_scalar_field<multem::eZero_Defocus_Type>(mx_input_multislice, "cond_lens_zero_defocus_type");
	input_multislice.cond_lens.zero_defocus_plane = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_zero_defocus_plane");	
	input_multislice.cond_lens.set_input_data(input_multislice.E_0, input_multislice.grid);

	/****************************** Objective lens ********************************/
	input_multislice.obj_lens.m = mx_get_scalar_field<int>(mx_input_multislice, "obj_lens_m"); 												// momentum of the vortex
	input_multislice.obj_lens.f = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_f"); 									// defocus(Angstrom)
	input_multislice.obj_lens.Cs3 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_Cs3")*multem::c_mm_2_Ags; 				// third order spherical aberration(mm-->Angstrom)
	input_multislice.obj_lens.Cs5 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_Cs5")*multem::c_mm_2_Ags; 				// fifth order aberration(mm-->Angstrom)
	input_multislice.obj_lens.mfa2 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_mfa2"); 								// magnitude 2-fold astigmatism(Angstrom)
	input_multislice.obj_lens.afa2 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_afa2")*multem::c_deg_2_rad; 			// angle 2-fold astigmatism(degrees-->rad)
	input_multislice.obj_lens.mfa3 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_mfa3"); 								// magnitude 3-fold astigmatism(Angstrom)
	input_multislice.obj_lens.afa3 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_afa3")*multem::c_deg_2_rad; 			// angle 3-fold astigmatism(degrees-->rad)
	input_multislice.obj_lens.inner_aper_ang = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_inner_aper_ang")*multem::c_mrad_2_rad; 		// inner aperture(mrad-->rad)
	input_multislice.obj_lens.outer_aper_ang = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_outer_aper_ang")*multem::c_mrad_2_rad; 		// outer aperture(mrad-->rad)
	input_multislice.obj_lens.sf = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_sf"); 									// defocus spread(Angstrom)
	input_multislice.obj_lens.nsf = mx_get_scalar_field<int>(mx_input_multislice, "obj_lens_nsf"); 											// Number of integration steps for the defocus Spread
	input_multislice.obj_lens.beta = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_beta")*multem::c_mrad_2_rad; 			// divergence semi-angle(mrad-->rad)
	input_multislice.obj_lens.nbeta = mx_get_scalar_field<int>(mx_input_multislice, "obj_lens_nbeta"); 										// Number of integration steps for the divergence semi-angle
	input_multislice.obj_lens.zero_defocus_type = mx_get_scalar_field<multem::eZero_Defocus_Type>(mx_input_multislice, "obj_lens_zero_defocus_type");
	input_multislice.obj_lens.zero_defocus_plane = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_zero_defocus_plane");		
	input_multislice.obj_lens.set_input_data(input_multislice.E_0, input_multislice.grid);

	input_multislice.validate_parameters();
}

void set_output_multislice(const mxArray *mx_input_multislice, mxArray *&mx_output_multislice, multem::Output_Multislice_Matlab &output_multislice)
{
	multem::Input_Multislice<double> input_multislice;
	read_input_multislice(mx_input_multislice, input_multislice);
	output_multislice.set_input_data(&input_multislice);

	const char *field_names_output_multislice[] = {"dx", "dy", "x", "y", "thickness", "data"};
	int number_of_fields_output_multislice = 6;
	mwSize dims_output_multislice[2] = {1, 1};

	mxArray *mx_field_psi;
	const char *field_names_psi[] = {"psi_coh"};
	int number_of_fields_psi = 1;
	mwSize dims_psi[2] = {1, output_multislice.thickness.size()};

	mx_output_multislice = mxCreateStructArray(2, dims_output_multislice, number_of_fields_output_multislice, field_names_output_multislice);

	mx_create_set_scalar_field<rmatrix_r>(mx_output_multislice, "dx", output_multislice.dx);
	mx_create_set_scalar_field<rmatrix_r>(mx_output_multislice, "dy", output_multislice.dy);
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "x", 1, output_multislice.x.size(), output_multislice.x.data());
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "y", 1, output_multislice.y.size(), output_multislice.y.data());
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "thickness", 1, output_multislice.thickness.size(), output_multislice.thickness.data());

	mx_field_psi = mxCreateStructArray(2, dims_psi, number_of_fields_psi, field_names_psi);
	mxSetField(mx_output_multislice, 0, "data", mx_field_psi);

	for(auto ithk = 0; ithk<output_multislice.thickness.size(); ithk++)
	{
		output_multislice.psi_coh[ithk] = mx_create_matrix_field<rmatrix_c>(mx_field_psi, ithk, "psi_coh", output_multislice.ny, output_multislice.nx);
	}
}

template<class T, multem::eDevice dev>
void il_wave_function(const mxArray *mxB, multem::Output_Multislice_Matlab &output_multislice)
{
	multem::Input_Multislice<T> input_multislice;
	read_input_multislice(mxB, input_multislice);

	multem::Stream<dev> stream;
	multem::FFT2<T, dev> fft2;
	multem::Wave_Function<T, dev> wave_function;

	stream.resize(input_multislice.nstream);
	fft2.create_plan(input_multislice.grid.ny, input_multislice.grid.nx, input_multislice.nstream);
	wave_function.set_input_data(&input_multislice, &stream, &fft2);

	wave_function.move_atoms(input_multislice.fp_nconf);
	wave_function.set_incident_wave(wave_function.psi_z);
	wave_function.psi(1.0, wave_function.psi_z, output_multislice);

	output_multislice.shift();
	output_multislice.clear_temporal_data();

	fft2.cleanup();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	multem::Output_Multislice_Matlab output_multislice;
	set_output_multislice(prhs[0], plhs[0], output_multislice);

	if(output_multislice.is_float_host())
	{
		il_wave_function<float, multem::e_host>(prhs[0], output_multislice);
	}
	else if(output_multislice.is_double_host())
	{
		il_wave_function<double, multem::e_host>(prhs[0], output_multislice);
	}
	if(output_multislice.is_float_device())
	{
		il_wave_function<float, multem::e_device>(prhs[0], output_multislice);
	}
	else if(output_multislice.is_double_device())
	{
		il_wave_function<double, multem::e_device>(prhs[0], output_multislice);
	}
}