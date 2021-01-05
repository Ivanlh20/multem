/**
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

#include "types.cuh"
#include "matlab_types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "fft.cuh"
#include "atomic_data_mt.hpp"
#include "input_multislice.cuh"
#include "output_multislice.hpp"

#include "microscope_effects.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::rmatrix_r;
using mt::rmatrix_c;

template <class TInput_Multislice>
void read_input_multislice(const mxArray *mx_input_multislice, TInput_Multislice &input_multislice, bool full = true)
{
	using T_r = mt::Value_type<TInput_Multislice>;

	/************************ simulation type **************************/
	input_multislice.simulation_type = mt::eTEMST_HRTEM;

	/**************************** Specimen *****************************/
	auto lx = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_lx");
	auto ly = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_ly");
	T_r lz = 0;
	T_r dz = 0.25;
	bool pbc_xy = true;

	/************************** xy sampling ****************************/
	auto nx = mx_get_scalar_field<int>(mx_input_multislice, "nx");
	auto ny = mx_get_scalar_field<int>(mx_input_multislice, "ny");
	bool bwl = false;

	input_multislice.grid_2d.set_input_data(nx, ny, lx, ly, dz, bwl, pbc_xy);

	/************************ Incident wave ****************************/
	auto iw_type = mx_get_scalar_field<mt::eIncident_Wave_Type>(mx_input_multislice, "iw_type");
	input_multislice.set_incident_wave_type(iw_type);

	if (input_multislice.is_user_define_wave() && full)
	{
		auto iw_psi = mx_get_matrix_field<rmatrix_c>(mx_input_multislice, "iw_psi");
		mt::assign(iw_psi, input_multislice.iw_psi);
	}

	// read iw_x and iw_y
	auto iw_x = mx_get_matrix_field<rmatrix_r>(mx_input_multislice, "iw_x");
	auto iw_y = mx_get_matrix_field<rmatrix_r>(mx_input_multislice, "iw_y");

	int n_iw_xy = min(iw_x.size(), iw_y.size());
	input_multislice.iw_x.assign(iw_x.begin(), iw_x.begin() + n_iw_xy);
	input_multislice.iw_y.assign(iw_y.begin(), iw_y.begin() + n_iw_xy);

	/********************* Microscope parameter ***********************/
	input_multislice.E_0 = mx_get_scalar_field<T_r>(mx_input_multislice, "E_0");
	input_multislice.theta = mx_get_scalar_field<T_r>(mx_input_multislice, "theta")*mt::c_deg_2_rad;
	input_multislice.phi = mx_get_scalar_field<T_r>(mx_input_multislice, "phi")*mt::c_deg_2_rad;

	/********************* Illumination model *************************/
	input_multislice.illumination_model = mx_get_scalar_field<mt::eIllumination_Model>(mx_input_multislice, "illumination_model");
	input_multislice.temporal_spatial_incoh = mx_get_scalar_field<mt::eTemporal_Spatial_Incoh>(mx_input_multislice, "temporal_spatial_incoh");

	/********************* source spread function *********************/
	input_multislice.cond_lens.si_sigma = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_si_sigma"); 								// standard deviation: For parallel ilumination(Å^-1); otherwise (Å)
	input_multislice.cond_lens.si_rad_npts = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_si_rad_npts");							// # of integration points

	input_multislice.cond_lens.set_input_data(input_multislice.E_0, input_multislice.grid_2d);

	/*********************** Objective lens ***************************/
	input_multislice.obj_lens.m = mx_get_scalar_field<int>(mx_input_multislice, "obj_lens_m");													// momentum of the vortex
	input_multislice.obj_lens.c_10 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_c_10"); 											// defocus (Angstrom)
	input_multislice.obj_lens.c_12 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_c_12"); 											// 2-fold astigmatism (Angstrom)
	input_multislice.obj_lens.phi_12 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_phi_12")*mt::c_deg_2_rad; 						// Azimuthal angle of 2-fold astigmatism (degrees-->rad)

	input_multislice.obj_lens.c_21 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_c_21"); 											// Axial coma (Angstrom)
	input_multislice.obj_lens.phi_21 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_phi_21")*mt::c_deg_2_rad; 						// Azimuthal angle of axial coma (degrees-->rad)
	input_multislice.obj_lens.c_23 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_c_23"); 											// 3-fold astigmatism (Angstrom)
	input_multislice.obj_lens.phi_23 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_phi_23")*mt::c_deg_2_rad; 						// Azimuthal angle of 3-fold astigmatism (degrees-->rad)

	input_multislice.obj_lens.c_30 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_c_30")*mt::c_mm_2_Angs; 							// 3rd order spherical aberration (mm-->Angstrom)
	input_multislice.obj_lens.c_32 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_c_32"); 											// Axial star aberration (Angstrom)
	input_multislice.obj_lens.phi_32 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_phi_32")*mt::c_deg_2_rad; 						// Azimuthal angle of axial star aberration (degrees-->rad)
	input_multislice.obj_lens.c_34 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_c_34"); 											// 4-fold astigmatism (Angstrom)
	input_multislice.obj_lens.phi_34 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_phi_34")*mt::c_deg_2_rad; 						// Azimuthal angle of 4-fold astigmatism (degrees-->rad)

	input_multislice.obj_lens.c_41 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_c_41"); 											// 4th order axial coma (Angstrom)
	input_multislice.obj_lens.phi_41 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_phi_41")*mt::c_deg_2_rad; 						// Azimuthal angle of 4th order axial coma (degrees-->rad)
	input_multislice.obj_lens.c_43 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_c_43"); 											// 3-lobe aberration (Angstrom)
	input_multislice.obj_lens.phi_43 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_phi_43")*mt::c_deg_2_rad; 						// Azimuthal angle of 3-lobe aberration (degrees-->rad)
	input_multislice.obj_lens.c_45 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_c_45"); 											// 5-fold astigmatism (Angstrom)
	input_multislice.obj_lens.phi_45 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_phi_45")*mt::c_deg_2_rad; 						// Azimuthal angle of 5-fold astigmatism (degrees-->rad)

	input_multislice.obj_lens.c_50 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_c_50")*mt::c_mm_2_Angs; 							// 5th order spherical aberration (mm-->Angstrom)
	input_multislice.obj_lens.c_52 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_c_52"); 											// 5th order axial star aberration (Angstrom)
	input_multislice.obj_lens.phi_52 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_phi_52")*mt::c_deg_2_rad; 						// Azimuthal angle of 5th order axial star aberration (degrees-->rad)
	input_multislice.obj_lens.c_54 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_c_54"); 											// 5th order rosette aberration (Angstrom)
	input_multislice.obj_lens.phi_54 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_phi_54")*mt::c_deg_2_rad; 						// Azimuthal angle of 5th order rosette aberration(degrees-->rad)
	input_multislice.obj_lens.c_56 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_c_56"); 											// 6-fold astigmatism (Angstrom)
	input_multislice.obj_lens.phi_56 = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_phi_56")*mt::c_deg_2_rad; 						// Azimuthal angle of 6-fold astigmatism (degrees-->rad)

	input_multislice.obj_lens.inner_aper_ang = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_inner_aper_ang")*mt::c_mrad_2_rad; 		// inner aperture (mrad-->rad)
	input_multislice.obj_lens.outer_aper_ang = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_outer_aper_ang")*mt::c_mrad_2_rad; 		// outer aperture (mrad-->rad)

	/********************* defocus spread function ********************/
	input_multislice.obj_lens.ti_a = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_ti_a"); 											// Height proportion of a normalized Gaussian [0, 1]
	input_multislice.obj_lens.ti_sigma = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_ti_sigma");   									// Standard deviation of the source spread function for the Gaussian component: For parallel ilumination(Å^-1); otherwise (Å)
	input_multislice.obj_lens.ti_beta = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_ti_beta"); 										// Standard deviation of the source spread function for the Exponential component: For parallel ilumination(Å^-1); otherwise (Å)
	input_multislice.obj_lens.ti_npts = mx_get_scalar_field<int>(mx_input_multislice, "obj_lens_ti_npts");										// Number of integration points

	/********************* source spread function *********************/
	input_multislice.obj_lens.si_a = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_si_a"); 											// Height proportion of a normalized Gaussian [0, 1]
	input_multislice.obj_lens.si_sigma = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_si_sigma");   								// Standard deviation of the source spread function for the Gaussian component: For parallel ilumination(Å^-1); otherwise (Å)
	input_multislice.obj_lens.si_beta = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_si_beta"); 									// Standard deviation of the source spread function for the Exponential component: For parallel ilumination(Å^-1); otherwise (Å)
	input_multislice.obj_lens.si_rad_npts = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_si_rad_npts"); 	 						// Number of radial integration points
	input_multislice.obj_lens.si_azm_npts = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_si_azm_npts");    							// Number of azimuth integration points

	/********************* zero defocus reference ********************/
	input_multislice.obj_lens.zero_defocus_type = mx_get_scalar_field<mt::eZero_Defocus_Type>(mx_input_multislice, "obj_lens_zero_defocus_type");	// Zero defocus type
	input_multislice.obj_lens.zero_defocus_plane = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_zero_defocus_plane");					// Zero defocus position

	input_multislice.obj_lens.set_input_data(input_multislice.E_0, input_multislice.grid_2d);

	/********************* select output region *************************/
	input_multislice.output_area.ix_0 = mx_get_scalar_field<int>(mx_input_multislice, "output_area_ix_0")-1;
	input_multislice.output_area.iy_0 = mx_get_scalar_field<int>(mx_input_multislice, "output_area_iy_0")-1;
	input_multislice.output_area.ix_e = mx_get_scalar_field<int>(mx_input_multislice, "output_area_ix_e")-1;
	input_multislice.output_area.iy_e = mx_get_scalar_field<int>(mx_input_multislice, "output_area_iy_e")-1;

	/********************* validate parameters *************************/
	input_multislice.validate_parameters();
}

template<class TOutput_Multislice>
void set_struct_microscope_aberrations(TOutput_Multislice &output_multislice, mxArray *&mx_output_multislice)
{
	const char *field_names_output_multislice[] = { "dx", "dy", "x", "y", "thick", "m2psi" };
	int number_of_fields_output_multislice = 6;
	mwSize dims_output_multislice[2] = { 1, 1 };

	mx_output_multislice = mxCreateStructArray(2, dims_output_multislice, number_of_fields_output_multislice, field_names_output_multislice);

	mx_create_set_scalar_field<rmatrix_r>(mx_output_multislice, 0, "dx", output_multislice.dx);
	mx_create_set_scalar_field<rmatrix_r>(mx_output_multislice, 0, "dy", output_multislice.dy);
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "x", 1, output_multislice.x.size(), output_multislice.x);
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "y", 1, output_multislice.y.size(), output_multislice.y);
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "thick", 1, output_multislice.thick.size(), output_multislice.thick);
	mx_create_set_matrix_field<rmatrix_c>(mx_output_multislice, "m2psi", output_multislice.ny, output_multislice.nx, output_multislice.m2psi_tot[0]);
}

template <class T, mt::eDevice dev>
void run_microscope_aberrations(mt::System_Configuration &system_conf, const mxArray *mx_input_multislice, mxArray *&mx_output_multislice)
{
	mt::Input_Multislice<T> input_multislice;
	read_input_multislice(mx_input_multislice, input_multislice);
	input_multislice.system_conf = system_conf;

	mt::Stream<dev> stream(system_conf.nstream);
	mt::FFT<T, dev> fft_2d;
	fft_2d.create_plan_2d(input_multislice.grid_2d.ny, input_multislice.grid_2d.nx, system_conf.nstream);

	mt::Microscope_Effects<T, dev> microscope_effects;
	microscope_effects.set_input_data(&input_multislice, &stream, &fft_2d);

	mt::Output_Multislice<T> output_multislice;
	output_multislice.set_input_data(&input_multislice);

	microscope_effects(output_multislice);

	stream.synchronize();

	//output_multislice.gather();
	output_multislice.clean_temporal();
	fft_2d.cleanup();

	set_struct_microscope_aberrations(output_multislice, mx_output_multislice);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	auto system_conf = mt::read_system_conf(prhs[0]);
	int idx_0 = (system_conf.active) ? 1 : 0;

	if (system_conf.is_float_host())
	{
		run_microscope_aberrations<float, mt::e_host>(system_conf, prhs[idx_0], plhs[0]);
	}
	else if (system_conf.is_double_host())
	{
		run_microscope_aberrations<double, mt::e_host>(system_conf, prhs[idx_0], plhs[0]);
	}
	else if (system_conf.is_float_device())
	{
		run_microscope_aberrations<float, mt::e_device>(system_conf, prhs[idx_0], plhs[0]);
	}
	else if (system_conf.is_double_device())
	{
		run_microscope_aberrations<double, mt::e_device>(system_conf, prhs[idx_0], plhs[0]);
	}
}