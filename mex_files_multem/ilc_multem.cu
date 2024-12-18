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

#include <algorithm>

#include "types.cuh"
#include "matlab_types.cuh"
#include "traits.cuh"
#include "cgpu_fcns.cuh"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "atomic_data_mt.hpp"
#include "tem_simulation.cuh"
#include "timing.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::rmatrix_r;
using mt::rmatrix_c;

template <class TInput_Multislice>
void read_input_multislice(const mxArray *mx_input_multislice, TInput_Multislice &input_multislice, bool full = true)
{
	using T_r = mt::Value_type<TInput_Multislice>;

	/************************ simulation type **************************/
	input_multislice.simulation_type = mx_get_scalar_field<mt::eTEM_Sim_Type>(mx_input_multislice, "simulation_type");

	/*******************************************************************/
	input_multislice.interaction_model = mx_get_scalar_field<mt::eElec_Spec_Int_Model>(mx_input_multislice, "interaction_model");
	input_multislice.potential_type = mx_get_scalar_field<mt::ePotential_Type>(mx_input_multislice, "potential_type");

	/*******************************************************************/
	input_multislice.operation_mode = mx_get_scalar_field<mt::eOperation_Mode>(mx_input_multislice, "operation_mode");
	input_multislice.reverse_multislice = mx_get_scalar_field<bool>(mx_input_multislice, "reverse_multislice");

	/************** Electron-Phonon interaction model ******************/
	input_multislice.pn_model = mx_get_scalar_field<mt::ePhonon_Model>(mx_input_multislice, "pn_model");
	input_multislice.pn_coh_contrib = mx_get_scalar_field<bool>(mx_input_multislice, "pn_coh_contrib");
	input_multislice.pn_single_conf = mx_get_scalar_field<bool>(mx_input_multislice, "pn_single_conf");
	input_multislice.pn_nconf = mx_get_scalar_field<int>(mx_input_multislice, "pn_nconf");
	input_multislice.pn_dim.set(mx_get_scalar_field<int>(mx_input_multislice, "pn_dim"));
	input_multislice.pn_seed = mx_get_scalar_field<int>(mx_input_multislice, "pn_seed");

	/**************************** Specimen *****************************/
	auto lx = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_lx");
	auto ly = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_ly");
	auto lz = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_lz");
	auto dz = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_dz");
	bool pbc_xy = true;

	if (input_multislice.is_specimen_required())
	{
		auto atoms = mx_get_matrix_field<rmatrix_r>(mx_input_multislice, "spec_atoms");

		auto ct_na = mx_get_scalar_field<int>(mx_input_multislice, "spec_cryst_na");
		auto ct_nb = mx_get_scalar_field<int>(mx_input_multislice, "spec_cryst_nb");
		auto ct_nc = mx_get_scalar_field<int>(mx_input_multislice, "spec_cryst_nc");
		auto ct_a = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_cryst_a");
		auto ct_b = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_cryst_b");
		auto ct_c = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_cryst_c");
		auto ct_x0 = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_cryst_x0");
		auto ct_y0 = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_cryst_y0");

		auto mx_spec_amorp = mxGetField(mx_input_multislice, 0, "spec_amorp");
		mt::Vector<mt::Amorp_Lay_Info<T_r>, mt::e_host> amorp_lay_info(mxGetN(mx_spec_amorp));
		for (auto i = 0; i < amorp_lay_info.size(); i++)
		{
			amorp_lay_info[i].z_0 = mx_get_scalar_field<T_r>(mx_spec_amorp, i, "z_0");
			amorp_lay_info[i].z_e = mx_get_scalar_field<T_r>(mx_spec_amorp, i, "z_e");
			amorp_lay_info[i].dz = mx_get_scalar_field<T_r>(mx_spec_amorp, i, "dz");
		}

		if (full)
		{
			input_multislice.atoms.set_crystal_parameters(ct_na, ct_nb, ct_nc, ct_a, ct_b, ct_c, ct_x0, ct_y0);
			input_multislice.atoms.set_amorphous_parameters(amorp_lay_info);
			input_multislice.atoms.set_atoms(atoms.rows, atoms.cols, atoms.real, lx, ly, lz, dz);
		}

		/************************ Specimen rotation *************************/
		input_multislice.spec_rot_theta = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_rot_theta")*mt::c_deg_2_rad;
		input_multislice.spec_rot_u0 = mx_get_r3d_field<T_r>(mx_input_multislice, "spec_rot_u0");
		input_multislice.spec_rot_u0.normalized();
		input_multislice.spec_rot_center_type = mx_get_scalar_field<mt::eRot_Point_Type>(mx_input_multislice, "spec_rot_center_type");
		input_multislice.spec_rot_center_p = mx_get_r3d_field<T_r>(mx_input_multislice, "spec_rot_center_p");

		/************************ Specimen thickness ***********************/
		input_multislice.thick_type = mx_get_scalar_field<mt::eThick_Type>(mx_input_multislice, "thick_type");
		if (!input_multislice.is_whole_spec() && full)
		{
			auto thick = mx_get_matrix_field<rmatrix_r>(mx_input_multislice, "thick");
			mt::assign(thick, input_multislice.thick);
		}

		/************************ Potential slicing ************************/
		input_multislice.potential_slicing = mx_get_scalar_field<mt::ePotential_Slicing>(mx_input_multislice, "potential_slicing");
	}

	/************************** xy sampling ****************************/
	auto nx = mx_get_scalar_field<int>(mx_input_multislice, "nx");
	auto ny = mx_get_scalar_field<int>(mx_input_multislice, "ny");
	bool bwl = mx_get_scalar_field<bool>(mx_input_multislice, "bwl");

	input_multislice.grid_2d.set_input_data(nx, ny, lx, ly, dz, bwl, pbc_xy);

	/************************ Incident wave ****************************/
	auto iw_type = mx_get_scalar_field<mt::eIncident_Wave_Type>(mx_input_multislice, "iw_type");
	input_multislice.set_incident_wave_type(iw_type);

	if(input_multislice.is_user_define_wave() && full)
	{
		auto iw_psi = mx_get_matrix_field<rmatrix_c>(mx_input_multislice, "iw_psi");
		mt::assign(iw_psi, input_multislice.iw_psi);
	}

	// read iw_x and iw_y
	auto iw_x = mx_get_matrix_field<rmatrix_r>(mx_input_multislice, "iw_x");
	auto iw_y = mx_get_matrix_field<rmatrix_r>(mx_input_multislice, "iw_y");
	
	int n_iw_xy = min(iw_x.size(), iw_y.size()); 
	input_multislice.iw_x.assign(iw_x.begin(), iw_x.begin()+n_iw_xy);
	input_multislice.iw_y.assign(iw_y.begin(), iw_y.begin()+n_iw_xy);

	/********************* Microscope parameter ***********************/
	input_multislice.E_0 = mx_get_scalar_field<T_r>(mx_input_multislice, "E_0");
	input_multislice.theta = mx_get_scalar_field<T_r>(mx_input_multislice, "theta")*mt::c_deg_2_rad;
	input_multislice.phi = mx_get_scalar_field<T_r>(mx_input_multislice, "phi")*mt::c_deg_2_rad;

	/********************* Illumination model *************************/
	input_multislice.illumination_model = mx_get_scalar_field<mt::eIllumination_Model>(mx_input_multislice, "illumination_model");
	input_multislice.temporal_spatial_incoh = mx_get_scalar_field<mt::eTemporal_Spatial_Incoh>(mx_input_multislice, "temporal_spatial_incoh");

	/*********************** Condenser lens ***************************/
	input_multislice.cond_lens.m = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_m");												// momentum of the vortex
	input_multislice.cond_lens.c_10 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_c_10"); 											// defocus (Angstrom)
	input_multislice.cond_lens.c_12 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_c_12"); 											// 2-fold astigmatism (Angstrom)
	input_multislice.cond_lens.phi_12 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_phi_12")*mt::c_deg_2_rad; 						// Azimuthal angle of 2-fold astigmatism (degrees-->rad)

	input_multislice.cond_lens.c_21 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_c_21"); 											// Axial coma (Angstrom)
	input_multislice.cond_lens.phi_21 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_phi_21")*mt::c_deg_2_rad; 						// Azimuthal angle of axial coma (degrees-->rad)
	input_multislice.cond_lens.c_23 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_c_23"); 											// 3-fold astigmatism (Angstrom)
	input_multislice.cond_lens.phi_23 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_phi_23")*mt::c_deg_2_rad; 						// Azimuthal angle of 3-fold astigmatism (degrees-->rad)

	input_multislice.cond_lens.c_30 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_c_30")*mt::c_mm_2_Angs; 							// 3rd order spherical aberration (mm-->Angstrom)
	input_multislice.cond_lens.c_32 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_c_32"); 											// Axial star aberration (Angstrom)
	input_multislice.cond_lens.phi_32 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_phi_32")*mt::c_deg_2_rad; 						// Azimuthal angle of axial star aberration (degrees-->rad)
	input_multislice.cond_lens.c_34 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_c_34"); 											// 4-fold astigmatism (Angstrom)
	input_multislice.cond_lens.phi_34 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_phi_34")*mt::c_deg_2_rad; 						// Azimuthal angle of 4-fold astigmatism (degrees-->rad)

	input_multislice.cond_lens.c_41 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_c_41"); 											// 4th order axial coma (Angstrom)
	input_multislice.cond_lens.phi_41 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_phi_41")*mt::c_deg_2_rad; 						// Azimuthal angle of 4th order axial coma (degrees-->rad)
	input_multislice.cond_lens.c_43 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_c_43"); 											// 3-lobe aberration (Angstrom)
	input_multislice.cond_lens.phi_43 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_phi_43")*mt::c_deg_2_rad; 						// Azimuthal angle of 3-lobe aberration (degrees-->rad)
	input_multislice.cond_lens.c_45 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_c_45"); 											// 5-fold astigmatism (Angstrom)
	input_multislice.cond_lens.phi_45 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_phi_45")*mt::c_deg_2_rad; 						// Azimuthal angle of 5-fold astigmatism (degrees-->rad)

	input_multislice.cond_lens.c_50 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_c_50")*mt::c_mm_2_Angs; 							// 5th order spherical aberration (mm-->Angstrom)
	input_multislice.cond_lens.c_52 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_c_52"); 											// 5th order axial star aberration (Angstrom)
	input_multislice.cond_lens.phi_52 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_phi_52")*mt::c_deg_2_rad; 						// Azimuthal angle of 5th order axial star aberration (degrees-->rad)
	input_multislice.cond_lens.c_54 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_c_54"); 											// 5th order rosette aberration (Angstrom)
	input_multislice.cond_lens.phi_54 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_phi_54")*mt::c_deg_2_rad; 						// Azimuthal angle of 5th order rosette aberration (degrees-->rad)
	input_multislice.cond_lens.c_56 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_c_56"); 											// 6-fold astigmatism (Angstrom)
	input_multislice.cond_lens.phi_56 = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_phi_56")*mt::c_deg_2_rad; 						// Azimuthal angle of 6-fold astigmatism (degrees-->rad)

	input_multislice.cond_lens.inner_aper_ang = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_inner_aper_ang")*mt::c_mrad_2_rad;		// inner aperture (mrad-->rad)
	input_multislice.cond_lens.outer_aper_ang = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_outer_aper_ang")*mt::c_mrad_2_rad;		// outer aperture (mrad-->rad)

	/********************* defocus spread function ********************/
	input_multislice.cond_lens.ti_a = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_ti_a"); 											// Height proportion of a normalized Gaussian [0, 1] 
	input_multislice.cond_lens.ti_sigma = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_ti_sigma");   								// Standard deviation of the source spread function for the Gaussian component: For parallel ilumination(�^-1); otherwise (�)
	input_multislice.cond_lens.ti_beta = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_ti_beta"); 									// Standard deviation of the source spread function for the Exponential component: For parallel ilumination(�^-1); otherwise (�)
	input_multislice.cond_lens.ti_npts = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_ti_npts");									// Number of integration points

	/********************* source spread function *********************/
	input_multislice.cond_lens.si_a = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_si_a"); 											// Height proportion of a normalized Gaussian [0, 1] 
	input_multislice.cond_lens.si_sigma = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_si_sigma");   								// Standard deviation of the source spread function for the Gaussian component: For parallel ilumination(�^-1); otherwise (�)
	input_multislice.cond_lens.si_beta = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_si_beta"); 									// Standard deviation of the source spread function for the Exponential component: For parallel ilumination(�^-1); otherwise (�)
	input_multislice.cond_lens.si_rad_npts = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_si_rad_npts"); 	 						// Number of radial integration points
	input_multislice.cond_lens.si_azm_npts = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_si_azm_npts"); 

	/********************* zero defocus reference ********************/																			/********************* zero defocus reference ********************/
	input_multislice.cond_lens.zero_defocus_type = mx_get_scalar_field<mt::eZero_Defocus_Type>(mx_input_multislice, "cond_lens_zero_defocus_type");	// Zero defocus type
	input_multislice.cond_lens.zero_defocus_plane = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_zero_defocus_plane");					// Zero defocus position

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
	input_multislice.obj_lens.ti_sigma = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_ti_sigma");   									// Standard deviation of the source spread function for the Gaussian component: For parallel ilumination(�^-1); otherwise (�)
	input_multislice.obj_lens.ti_beta = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_ti_beta"); 										// Standard deviation of the source spread function for the Exponential component: For parallel ilumination(�^-1); otherwise (�)
	input_multislice.obj_lens.ti_npts = mx_get_scalar_field<int>(mx_input_multislice, "obj_lens_ti_npts");										// Number of integration points

	/********************* source spread function *********************/
	input_multislice.obj_lens.si_a = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_si_a"); 											// Height proportion of a normalized Gaussian [0, 1]
	input_multislice.obj_lens.si_sigma = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_si_sigma");   								// Standard deviation of the source spread function for the Gaussian component: For parallel ilumination(�^-1); otherwise (�)
	input_multislice.obj_lens.si_beta = mx_get_scalar_field<T_r>(mx_input_multislice, "cond_lens_si_beta"); 									// Standard deviation of the source spread function for the Exponential component: For parallel ilumination(�^-1); otherwise (�)
	input_multislice.obj_lens.si_rad_npts = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_si_rad_npts"); 	 						// Number of radial integration points
	input_multislice.obj_lens.si_azm_npts = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_si_azm_npts");    							// Number of azimuth integration points

	/********************* zero defocus reference ********************/
	input_multislice.obj_lens.zero_defocus_type = mx_get_scalar_field<mt::eZero_Defocus_Type>(mx_input_multislice, "obj_lens_zero_defocus_type");	// Zero defocus type
	input_multislice.obj_lens.zero_defocus_plane = mx_get_scalar_field<T_r>(mx_input_multislice, "obj_lens_zero_defocus_plane");					// Zero defocus position

	input_multislice.obj_lens.set_input_data(input_multislice.E_0, input_multislice.grid_2d);

	/************************** ISTEM/STEM ***************************/
	if (input_multislice.is_scanning())
	{
		input_multislice.scanning.type = mx_get_scalar_field<mt::eScanning_Type>(mx_input_multislice, "scanning_type");
		input_multislice.scanning.pbc = mx_get_scalar_field<bool>(mx_input_multislice, "scanning_periodic");
		input_multislice.scanning.spxs = mx_get_scalar_field<bool>(mx_input_multislice, "scanning_square_pxs");
		input_multislice.scanning.ns = mx_get_scalar_field<int>(mx_input_multislice, "scanning_ns");
		input_multislice.scanning.x0 = mx_get_scalar_field<T_r>(mx_input_multislice, "scanning_x0");
		input_multislice.scanning.y0 = mx_get_scalar_field<T_r>(mx_input_multislice, "scanning_y0");
		input_multislice.scanning.xe = mx_get_scalar_field<T_r>(mx_input_multislice, "scanning_xe");
		input_multislice.scanning.ye = mx_get_scalar_field<T_r>(mx_input_multislice, "scanning_ye");
		input_multislice.scanning.set_grid();
	}

	if (input_multislice.is_STEM())
	{
		T_r lambda = mt::get_lambda(input_multislice.E_0);
		mxArray *mx_detector = mxGetField(mx_input_multislice, 0, "detector");
		input_multislice.detector.type = mx_get_scalar_field<mt::eDetector_Type>(mx_detector, "type");

		switch (input_multislice.detector.type)
		{
		case mt::eDT_Circular:
		{
			mx_detector = mxGetField(mx_detector, 0, "cir");
			int ndetector = mxGetN(mx_detector);
			if (ndetector > 0)
			{
				input_multislice.detector.resize(ndetector);
				for (auto i = 0; i < input_multislice.detector.size(); i++)
				{
					auto inner_ang = mx_get_scalar_field<T_r>(mx_detector, i, "inner_ang")*mt::c_mrad_2_rad;
					input_multislice.detector.g_inner[i] = sin(inner_ang) / lambda;
					auto outer_ang = mx_get_scalar_field<T_r>(mx_detector, i, "outer_ang")*mt::c_mrad_2_rad;
					input_multislice.detector.g_outer[i] = sin(outer_ang) / lambda;
				}
			}
		}
		break;
		case mt::eDT_Radial:
		{
			mx_detector = mxGetField(mx_detector, 0, "radial");
			int ndetector = mxGetN(mx_detector);
			if (ndetector > 0)
			{
				input_multislice.detector.resize(ndetector);
				for (auto i = 0; i < input_multislice.detector.size(); i++)
				{
					// auto x = mx_get_matrix_field<rmatrix_r>(mx_detector, i, "x");
					// mt::assign(x, input_multislice.detector.x[i]);
					// mt::scale(input_multislice.detector.x[i], 1.0/lambda);

					auto fx = mx_get_matrix_field<rmatrix_r>(mx_detector, i, "fx");
					mt::assign(fx, input_multislice.detector.fx[i]);
				}
			}
		}
		break;
		case mt::eDT_Matrix:
		{
			mx_detector = mxGetField(mx_detector, 0, "matrix");
			int ndetector = mxGetN(mx_detector);
			if (ndetector > 0)
			{
				input_multislice.detector.resize(ndetector);
				for (auto i = 0; i < input_multislice.detector.size(); i++)
				{
					// auto R = mx_get_matrix_field<rmatrix_r>(mx_detector, i, "x");
					// mt::assign(R, input_multislice.detector.R[i]);
					// mt::scale(input_multislice.detector.R[i], 1.0/lambda);
					// mt::fft2_shift(input_multislice.grid_2d, input_multislice.detector.R[i]);

					auto fR = mx_get_matrix_field<rmatrix_r>(mx_detector, i, "fR");
					mt::assign(fR, input_multislice.detector.fR[i]);
				}
			}
		}
		break;
		}
	}
	else if (input_multislice.is_PED())
	{
		input_multislice.theta = mx_get_scalar_field<T_r>(mx_input_multislice, "ped_theta")*mt::c_deg_2_rad;
		input_multislice.nrot = mx_get_scalar_field<T_r>(mx_input_multislice, "ped_nrot");
	}
	else if (input_multislice.is_HCTEM())
	{
		input_multislice.theta = mx_get_scalar_field<T_r>(mx_input_multislice, "hci_theta")*mt::c_deg_2_rad;
		input_multislice.nrot = mx_get_scalar_field<T_r>(mx_input_multislice, "hci_nrot");
	}
	else if (input_multislice.is_EELS())
	{
		mt::eSpace space = mt::eS_Reciprocal;
		auto Z = mx_get_scalar_field<int>(mx_input_multislice, "eels_Z");
		auto E_loss = mx_get_scalar_field<T_r>(mx_input_multislice, "eels_E_loss")*mt::c_eV_2_keV;
		auto collection_angle = mx_get_scalar_field<double>(mx_input_multislice, "eels_collection_angle")*mt::c_mrad_2_rad;
		auto m_selection = mx_get_scalar_field<int>(mx_input_multislice, "eels_m_selection");
		auto channelling_type = mx_get_scalar_field<mt::eChannelling_Type>(mx_input_multislice, "eels_channelling_type");

		input_multislice.eels_fr.set_input_data(space, input_multislice.E_0, E_loss, m_selection, collection_angle, channelling_type, Z);
	}
	else if (input_multislice.is_EFTEM())
	{
		mt::eSpace space = mt::eS_Real;
		auto Z = mx_get_scalar_field<int>(mx_input_multislice, "eftem_Z");
		auto E_loss = mx_get_scalar_field<T_r>(mx_input_multislice, "eftem_E_loss")*mt::c_eV_2_keV;
		auto collection_angle = mx_get_scalar_field<double>(mx_input_multislice, "eftem_collection_angle")*mt::c_mrad_2_rad;
		auto m_selection = mx_get_scalar_field<int>(mx_input_multislice, "eftem_m_selection");
		auto channelling_type = mx_get_scalar_field<mt::eChannelling_Type>(mx_input_multislice, "eftem_channelling_type");

		input_multislice.eels_fr.set_input_data(space, input_multislice.E_0, E_loss, m_selection, collection_angle, channelling_type, Z);
	}

	/********************* select output region *************************/
	input_multislice.output_area.ix_0 = mx_get_scalar_field<int>(mx_input_multislice, "output_area_ix_0")-1;
	input_multislice.output_area.iy_0 = mx_get_scalar_field<int>(mx_input_multislice, "output_area_iy_0")-1;
	input_multislice.output_area.ix_e = mx_get_scalar_field<int>(mx_input_multislice, "output_area_ix_e")-1;
	input_multislice.output_area.iy_e = mx_get_scalar_field<int>(mx_input_multislice, "output_area_iy_e")-1;

	/********************* validate parameters *************************/
	input_multislice.validate_parameters();
}

template<class TOutput_Multislice>
void set_struct_multislice(TOutput_Multislice &output_multislice, mxArray *&mx_output_multislice)
{
	const char *field_names_output_multislice[] = {"dx", "dy", "x", "y", "thick", "data"};
	int number_of_fields_output_multislice = 6;
	mwSize dims_output_multislice[2] = {1, 1};

	mx_output_multislice = mxCreateStructArray(2, dims_output_multislice, number_of_fields_output_multislice, field_names_output_multislice);

	mx_create_set_scalar_field<rmatrix_r>(mx_output_multislice, 0, "dx", output_multislice.dx);
	mx_create_set_scalar_field<rmatrix_r>(mx_output_multislice, 0, "dy", output_multislice.dy);
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "x", 1, output_multislice.x.size(), output_multislice.x);
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "y", 1, output_multislice.y.size(), output_multislice.y);
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "thick", 1, output_multislice.thick.size(), output_multislice.thick);

	if (output_multislice.is_STEM() || output_multislice.is_EELS())
	{
		mxArray *mx_field_data;
		const char *field_names_data_full[] = { "image_tot", "image_coh" };
		const char *field_names_data_partial[] = { "image_tot" };
		const char **field_names_data = (output_multislice.pn_coh_contrib) ? field_names_data_full : field_names_data_partial;
		int number_of_fields_data = (output_multislice.pn_coh_contrib) ? 2 : 1;
		mwSize dims_data[2] = { 1, output_multislice.thick.size() };

		mx_field_data = mxCreateStructArray(2, dims_data, number_of_fields_data, field_names_data);
		mxSetField(mx_output_multislice, 0, "data", mx_field_data);

		mxArray *mx_field_detector_tot;
		mxArray *mx_field_detector_coh;
		const char *field_names_detector[] = { "image" };
		int number_of_fields_detector = 1;
		// mwSize dims_detector[2] = {1, output_multislice.ndetector};
		mwSize dims_detector[2];
		dims_detector[0] = 1;
		dims_detector[1] = output_multislice.ndetector;

		int nx = (output_multislice.scanning.is_line()) ? 1 : output_multislice.nx;
		int ny = output_multislice.ny;
		for (auto ithk = 0; ithk < output_multislice.thick.size(); ithk++)
		{
			mx_field_detector_tot = mxCreateStructArray(2, dims_detector, number_of_fields_detector, field_names_detector);
			mxSetField(mx_field_data, ithk, "image_tot", mx_field_detector_tot);
			if (output_multislice.pn_coh_contrib)
			{
				mx_field_detector_coh = mxCreateStructArray(2, dims_detector, number_of_fields_detector, field_names_detector);
				mxSetField(mx_field_data, ithk, "image_coh", mx_field_detector_coh);
			}

			for (auto iDet = 0; iDet < output_multislice.ndetector; iDet++)
			{
				mx_create_set_matrix_field<rmatrix_r>(mx_field_detector_tot, iDet, "image", ny, nx, output_multislice.image_tot[ithk].image[iDet]);
				if (output_multislice.pn_coh_contrib)
				{
					mx_create_set_matrix_field<rmatrix_r>(mx_field_detector_coh, iDet, "image", ny, nx, output_multislice.image_coh[ithk].image[iDet]);
				}
			}
		}
	}
	else if (output_multislice.is_EWFS_EWRS())
	{
		mxArray *mx_field_data;
		const char *field_names_data_full[] = { "m2psi_tot", "psi_coh" };
		const char *field_names_data_partial[] = { "psi_coh" };
		const char **field_names_data = (!output_multislice.is_EWFS_EWRS_SC()) ? field_names_data_full : field_names_data_partial;
		int number_of_fields_data = (!output_multislice.is_EWFS_EWRS_SC()) ? 2 : 1;
		mwSize dims_data[2] = { 1, output_multislice.thick.size() };

		mx_field_data = mxCreateStructArray(2, dims_data, number_of_fields_data, field_names_data);
		mxSetField(mx_output_multislice, 0, "data", mx_field_data);

		for (auto ithk = 0; ithk < output_multislice.thick.size(); ithk++)
		{
			if (!output_multislice.is_EWFS_EWRS_SC())
			{
				mx_create_set_matrix_field<rmatrix_r>(mx_field_data, ithk, "m2psi_tot", output_multislice.ny, output_multislice.nx, output_multislice.m2psi_tot[ithk]);
			}
			mx_create_set_matrix_field<rmatrix_c>(mx_field_data, ithk, "psi_coh", output_multislice.ny, output_multislice.nx, output_multislice.psi_coh[ithk]);
		}
	}
	else
	{
		mxArray *mx_field_data;
		const char *field_names_data_full[] = { "m2psi_tot", "m2psi_coh" };
		const char *field_names_data_partial[] = { "m2psi_tot" };
		const char **field_names_data = (output_multislice.pn_coh_contrib) ? field_names_data_full : field_names_data_partial;
		int number_of_fields_data = (output_multislice.pn_coh_contrib) ? 2 : 1;
		mwSize dims_data[2] = { 1, output_multislice.thick.size() };

		mx_field_data = mxCreateStructArray(2, dims_data, number_of_fields_data, field_names_data);
		mxSetField(mx_output_multislice, 0, "data", mx_field_data);

		for (auto ithk = 0; ithk < output_multislice.thick.size(); ithk++)
		{
			mx_create_set_matrix_field<rmatrix_r>(mx_field_data, ithk, "m2psi_tot", output_multislice.ny, output_multislice.nx, output_multislice.m2psi_tot[ithk]);
			if (output_multislice.pn_coh_contrib)
			{
				mx_create_set_matrix_field<rmatrix_r>(mx_field_data, ithk, "m2psi_coh", output_multislice.ny, output_multislice.nx, output_multislice.m2psi_coh[ithk]);
			}
		}
	}
}

template <class T, mt::eDevice dev>
void run_multislice(mt::System_Configuration &system_conf, const mxArray *mx_input_multislice, 
mxArray *&mx_output_multislice)
{
	mt::Input_Multislice<T> input_multislice;
	read_input_multislice(mx_input_multislice, input_multislice);
	input_multislice.system_conf = system_conf;

	mt::Stream<dev> stream(system_conf.nstream);

	mt::FFT<T, dev> fft_2d;
	fft_2d.create_plan_2d(input_multislice.grid_2d.ny, input_multislice.grid_2d.nx, system_conf.nstream);

	mt::Multislice<T, dev> tem_simulation;
	tem_simulation.set_input_data(&input_multislice, &stream, &fft_2d);

	mt::Output_Multislice<T> output_multislice;
	output_multislice.set_input_data(&input_multislice);

	tem_simulation(output_multislice);
	stream.synchronize();

	output_multislice.gather();
	output_multislice.clean_temporal();
	fft_2d.cleanup();

	set_struct_multislice(output_multislice, mx_output_multislice);

	auto err = cudaGetLastError();
	if (err != cudaSuccess) {
		mexPrintf("CUDA error: %s\n", cudaGetErrorString(err));
	}
	//auto t3 = time.elapsed_ms();
	//mexPrintf("creation, execution, phase grating, deleting time = %7.5f, %7.5f, %7.5f, %7.5f\n", t1, t2, tx, t3);
	//mexPrintf("send device = %d, get device = %d\n", system_conf.gpu_device, system_conf.get_device());
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	auto system_conf = mt::read_system_conf(prhs[0]);
	int idx_0 = (system_conf.active)?1:0;

	// run_multislice<float, mt::e_device>(system_conf, prhs[idx_0], plhs[0]);

	if (system_conf.is_float_host())
	{
		//mexPrintf("cpu - float precision calculation\n");
		run_multislice<float, mt::e_host>(system_conf, prhs[idx_0], plhs[0]);
	}
	else if (system_conf.is_double_host())
	{
		//mexPrintf("cpu - double precision calculation\n");
		run_multislice<double, mt::e_host>(system_conf, prhs[idx_0], plhs[0]);
	}
	if (system_conf.is_float_device())
	{
		//mexPrintf("gpu - float precision calculation\n");
		run_multislice<float, mt::e_device>(system_conf, prhs[idx_0], plhs[0]);
	}
	else if (system_conf.is_double_device())
	{
		//mexPrintf("gpu - double precision calculation\n");
		run_multislice<double, mt::e_device>(system_conf, prhs[idx_0], plhs[0]);
	}
	// //cudaDeviceReset();
}