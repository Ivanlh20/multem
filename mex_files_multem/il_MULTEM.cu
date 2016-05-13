/*
 * This file is part of MULTEM.
 * Copyright 2016 Ivan Lobato <Ivanlh20@gmail.com>
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

#include "math.cuh"
#include "types.cuh"
#include "matlab_types.cuh"
#include "traits.cuh"
#include "host_device_functions.cuh"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "atom_data.hpp"
#include "multislice.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::rmatrix_r;
using mt::rmatrix_c;

template<class TInput_Multislice>
void read_input_multislice(const mxArray *mx_input_multislice, TInput_Multislice &input_multislice, bool full =true)
{
	using value_type_r = mt::Value_type<TInput_Multislice>;

	input_multislice.precision = mx_get_scalar_field<mt::ePrecision>(mx_input_multislice, "precision");
	input_multislice.device = mx_get_scalar_field<mt::eDevice>(mx_input_multislice, "device"); 
	input_multislice.cpu_ncores = mx_get_scalar_field<int>(mx_input_multislice, "cpu_ncores"); 
	input_multislice.cpu_nthread = mx_get_scalar_field<int>(mx_input_multislice, "cpu_nthread"); 
	input_multislice.gpu_device = mx_get_scalar_field<int>(mx_input_multislice, "gpu_device"); 
	input_multislice.gpu_nstream = mx_get_scalar_field<int>(mx_input_multislice, "gpu_nstream"); 

	input_multislice.simulation_type = mx_get_scalar_field<mt::eTEM_Sim_Type>(mx_input_multislice, "simulation_type"); 
	input_multislice.phonon_model = mx_get_scalar_field<mt::ePhonon_Model>(mx_input_multislice, "phonon_model"); 
	input_multislice.interaction_model = mx_get_scalar_field<mt::eElec_Spec_Int_Model>(mx_input_multislice, "interaction_model");
	input_multislice.potential_slicing = mx_get_scalar_field<mt::ePotential_Slicing>(mx_input_multislice, "potential_slicing");
	input_multislice.potential_type = mx_get_scalar_field<mt::ePotential_Type>(mx_input_multislice, "potential_type");

	input_multislice.fp_dim.set(mx_get_scalar_field<int>(mx_input_multislice, "fp_dim"));
	input_multislice.fp_seed = mx_get_scalar_field<int>(mx_input_multislice, "fp_seed");
	input_multislice.fp_single_conf = mx_get_scalar_field<bool>(mx_input_multislice, "fp_single_conf");
	input_multislice.fp_nconf = mx_get_scalar_field<int>(mx_input_multislice, "fp_nconf");

	input_multislice.tm_active = mx_get_scalar_field<bool>(mx_input_multislice, "tm_active");
	input_multislice.tm_theta = mx_get_scalar_field<value_type_r>(mx_input_multislice, "tm_theta")*mt::c_deg_2_rad;
	input_multislice.tm_u0 = mx_get_r3d_field<value_type_r>(mx_input_multislice, "tm_u0");
	input_multislice.tm_rot_point_type = mx_get_scalar_field<mt::eRot_Point_Type>(mx_input_multislice, "tm_rot_point_type");
	input_multislice.tm_p0 = mx_get_r3d_field<value_type_r>(mx_input_multislice, "tm_p0");

	input_multislice.illumination_model = mx_get_scalar_field<mt::eIllumination_Model>(mx_input_multislice, "illumination_model");
	input_multislice.temporal_spatial_incoh = mx_get_scalar_field<mt::eTemporal_Spatial_Incoh>(mx_input_multislice, "temporal_spatial_incoh");

	input_multislice.thickness_type = mx_get_scalar_field<mt::eThickness_Type>(mx_input_multislice, "thickness_type");
	if(!input_multislice.is_whole_specimen() && full)
	{
		auto thickness = mx_get_matrix_field<rmatrix_r>(mx_input_multislice, "thickness");
		input_multislice.thickness.resize(thickness.m_size);
		std::copy(thickness.real, thickness.real + thickness.m_size, input_multislice.thickness.begin());
	}

	input_multislice.operation_mode = mx_get_scalar_field<mt::eOperation_Mode>(mx_input_multislice, "operation_mode");
	input_multislice.coherent_contribution = mx_get_scalar_field<bool>(mx_input_multislice, "coherent_contribution");

	input_multislice.E_0 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "E_0");
	input_multislice.theta = mx_get_scalar_field<value_type_r>(mx_input_multislice, "theta")*mt::c_deg_2_rad;
	input_multislice.phi = mx_get_scalar_field<value_type_r>(mx_input_multislice, "phi")*mt::c_deg_2_rad;

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

	/****************************** Incident wave ********************************/
	auto iw_type = mx_get_scalar_field<mt::eIncident_Wave_Type>(mx_input_multislice, "iw_type");
	input_multislice.set_incident_wave_type(iw_type);

	if(input_multislice.is_user_define_wave() && full)
	{
		auto iw_psi = mx_get_matrix_field<rmatrix_c>(mx_input_multislice, "iw_psi");
		mt::Stream<mt::e_host> stream(input_multislice.cpu_nthread);
		mt::assign(stream, iw_psi, input_multislice.iw_psi);
	}
	input_multislice.iw_x = mx_get_scalar_field<value_type_r>(mx_input_multislice, "iw_x");
	input_multislice.iw_y = mx_get_scalar_field<value_type_r>(mx_input_multislice, "iw_y");

	/****************************** Condenser lens ********************************/
	input_multislice.cond_lens.m = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_m"); 											// momentum of the vortex
	input_multislice.cond_lens.f = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_f"); 									// defocus(Angstrom)
	input_multislice.cond_lens.Cs3 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_Cs3")*mt::c_mm_2_Ags; 			// third order spherical aberration(mm-->Angstrom)
	input_multislice.cond_lens.Cs5 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_Cs5")*mt::c_mm_2_Ags; 			// fifth order aberration(mm-->Angstrom)
	input_multislice.cond_lens.mfa2 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_mfa2"); 							// magnitude 2-fold astigmatism(Angstrom)
	input_multislice.cond_lens.afa2 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_afa2")*mt::c_deg_2_rad; 		// angle 2-fold astigmatism(degrees-->rad)
	input_multislice.cond_lens.mfa3 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_mfa3"); 							// magnitude 3-fold astigmatism(Angstrom)
	input_multislice.cond_lens.afa3 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_afa3")*mt::c_deg_2_rad; 		// angle 3-fold astigmatism(degrees-->rad)
	input_multislice.cond_lens.inner_aper_ang = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_inner_aper_ang")*mt::c_mrad_2_rad; 		// inner aperture(mrad-->rad)
	input_multislice.cond_lens.outer_aper_ang = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_outer_aper_ang")*mt::c_mrad_2_rad; 		// outer aperture(mrad-->rad)
	input_multislice.cond_lens.sf = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_sf"); 								// defocus spread(Angstrom)
	input_multislice.cond_lens.nsf = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_nsf"); 										// Number of integration steps for the defocus Spread
	input_multislice.cond_lens.beta = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_beta")*mt::c_mrad_2_rad; 		// divergence semi-angle(mrad-->rad)
	input_multislice.cond_lens.nbeta = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_nbeta");									// Number of integration steps for the divergence semi-angle
	input_multislice.cond_lens.zero_defocus_type = mx_get_scalar_field<mt::eZero_Defocus_Type>(mx_input_multislice, "cond_lens_zero_defocus_type");
	input_multislice.cond_lens.zero_defocus_plane = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_zero_defocus_plane");	
	input_multislice.cond_lens.set_input_data(input_multislice.E_0, input_multislice.grid);

	//input_multislice.cdl_var_type = mx_get_scalar_field<mt::eLens_Var_Type>(mx_input_multislice, "cdl_var_type");
	//if(!input_multislice.is_whole_specimen() && full)
	//{
	//	auto thickness = mx_get_matrix_field<rmatrix_r>(mx_input_multislice, "thickness");
	//	input_multislice.thickness.resize(thickness.m_size);
	//	std::copy(thickness.real, thickness.real + thickness.m_size, input_multislice.thickness.begin());
	//}

	/****************************** Objective lens ********************************/
	input_multislice.obj_lens.m = mx_get_scalar_field<int>(mx_input_multislice, "obj_lens_m"); 												// momentum of the vortex
	input_multislice.obj_lens.f = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_f"); 									// defocus(Angstrom)
	input_multislice.obj_lens.Cs3 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_Cs3")*mt::c_mm_2_Ags; 				// third order spherical aberration(mm-->Angstrom)
	input_multislice.obj_lens.Cs5 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_Cs5")*mt::c_mm_2_Ags; 				// fifth order aberration(mm-->Angstrom)
	input_multislice.obj_lens.mfa2 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_mfa2"); 								// magnitude 2-fold astigmatism(Angstrom)
	input_multislice.obj_lens.afa2 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_afa2")*mt::c_deg_2_rad; 			// angle 2-fold astigmatism(degrees-->rad)
	input_multislice.obj_lens.mfa3 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_mfa3"); 								// magnitude 3-fold astigmatism(Angstrom)
	input_multislice.obj_lens.afa3 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_afa3")*mt::c_deg_2_rad; 			// angle 3-fold astigmatism(degrees-->rad)
	input_multislice.obj_lens.inner_aper_ang = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_inner_aper_ang")*mt::c_mrad_2_rad; 		// inner aperture(mrad-->rad)
	input_multislice.obj_lens.outer_aper_ang = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_outer_aper_ang")*mt::c_mrad_2_rad; 		// outer aperture(mrad-->rad)
	input_multislice.obj_lens.sf = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_sf"); 									// defocus spread(Angstrom)
	input_multislice.obj_lens.nsf = mx_get_scalar_field<int>(mx_input_multislice, "obj_lens_nsf"); 											// Number of integration steps for the defocus Spread
	input_multislice.obj_lens.beta = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_beta")*mt::c_mrad_2_rad; 			// divergence semi-angle(mrad-->rad)
	input_multislice.obj_lens.nbeta = mx_get_scalar_field<int>(mx_input_multislice, "obj_lens_nbeta"); 										// Number of integration steps for the divergence semi-angle
	input_multislice.obj_lens.zero_defocus_type = mx_get_scalar_field<mt::eZero_Defocus_Type>(mx_input_multislice, "obj_lens_zero_defocus_type");
	input_multislice.obj_lens.zero_defocus_plane = mx_get_scalar_field<value_type_r>(mx_input_multislice, "obj_lens_zero_defocus_plane");		
	input_multislice.obj_lens.set_input_data(input_multislice.E_0, input_multislice.grid);

	if(input_multislice.is_scanning())
	{
		input_multislice.scanning.type = mx_get_scalar_field<mt::eScanning_Type>(mx_input_multislice, "scanning_type");
		input_multislice.scanning.pbc = mx_get_scalar_field<bool>(mx_input_multislice, "scanning_periodic");
		input_multislice.scanning.ns = mx_get_scalar_field<int>(mx_input_multislice, "scanning_ns");
		input_multislice.scanning.x0 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "scanning_x0");
		input_multislice.scanning.y0 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "scanning_y0");
		input_multislice.scanning.xe = mx_get_scalar_field<value_type_r>(mx_input_multislice, "scanning_xe");
		input_multislice.scanning.ye = mx_get_scalar_field<value_type_r>(mx_input_multislice, "scanning_ye");
		input_multislice.scanning.set_grid();
	}

	if(input_multislice.is_STEM())
	{
		value_type_r lambda = mt::get_lambda(input_multislice.E_0);
		mxArray *mx_detector = mxGetField(mx_input_multislice, 0, "detector");
		input_multislice.detector.type = mx_get_scalar_field<mt::eDetector_Type>(mx_detector, "type");

		switch (input_multislice.detector.type)
		{
		case mt::eDT_Circular:
			{
				mx_detector = mxGetField(mx_detector, 0, "cir");
				int ndetector = mxGetN(mx_detector);
				if(ndetector>0)
				{
					input_multislice.detector.resize(ndetector);
					for(auto i = 0; i<input_multislice.detector.size(); i++)
					{
						auto inner_ang = mx_get_scalar_field<value_type_r>(mx_detector, i, "inner_ang")*mt::c_mrad_2_rad;
						input_multislice.detector.g_inner[i] = sin(inner_ang)/lambda;
						auto outer_ang = mx_get_scalar_field<value_type_r>(mx_detector, i, "outer_ang")*mt::c_mrad_2_rad;
						input_multislice.detector.g_outer[i] = sin(outer_ang)/lambda;
					}
				}
			}
			break;
		case mt::eDT_Radial:
			{
				mx_detector = mxGetField(mx_detector, 0, "radial");
				int ndetector = mxGetN(mx_detector);
				if(ndetector>0)
				{
					mt::Stream<mt::e_host> stream(input_multislice.cpu_nthread);
					input_multislice.detector.resize(ndetector);
					for(auto i = 0; i<input_multislice.detector.size(); i++)
					{
						// auto x = mx_get_matrix_field<rmatrix_r>(mx_detector, i, "x");
						// mt::assign(x, input_multislice.detector.x[i]);
						// mt::scale(input_multislice.detector.x[i], 1.0/lambda);

						auto fx = mx_get_matrix_field<rmatrix_r>(mx_detector, i, "fx");
						mt::assign(stream, fx, input_multislice.detector.fx[i]);
					}
				}
			}
			break;
		case mt::eDT_Matrix:
			{
				mx_detector = mxGetField(mx_detector, 0, "matrix");
				int ndetector = mxGetN(mx_detector);
				if(ndetector>0)
				{
					mt::Stream<mt::e_host> stream(input_multislice.cpu_nthread);
					input_multislice.detector.resize(ndetector);
					for(auto i = 0; i<input_multislice.detector.size(); i++)
					{
						// auto R = mx_get_matrix_field<rmatrix_r>(mx_detector, i, "x");
						// mt::assign(R, input_multislice.detector.R[i]);
						// mt::scale(input_multislice.detector.R[i], 1.0/lambda);
						// mt::fft2_shift(input_multislice.grid, input_multislice.detector.R[i]);

						auto fR = mx_get_matrix_field<rmatrix_r>(mx_detector, i, "fR");
						mt::assign(stream, fR, input_multislice.detector.fR[i]);
					}
				}
			}
			break;
		}
	}
	else if (input_multislice.is_PED())
	{
		input_multislice.theta = mx_get_scalar_field<value_type_r>(mx_input_multislice, "ped_theta")*mt::c_deg_2_rad;
		input_multislice.nrot = mx_get_scalar_field<value_type_r>(mx_input_multislice, "ped_nrot");
	}
	else if (input_multislice.is_HCI())
	{
		input_multislice.theta = mx_get_scalar_field<value_type_r>(mx_input_multislice, "hci_theta")*mt::c_deg_2_rad;
		input_multislice.nrot = mx_get_scalar_field<value_type_r>(mx_input_multislice, "hci_nrot");
	}
	else if (input_multislice.is_EELS())
	{
		mt::eSpace space = mt::eS_Reciprocal;
		value_type_r E_loss = mx_get_scalar_field<value_type_r>(mx_input_multislice, "eels_E_loss")*mt::c_meV_2_keV;
		int m_selection = mx_get_scalar_field<int>(mx_input_multislice, "eels_m_selection");
		value_type_r collection_angle = mx_get_scalar_field<double>(mx_input_multislice, "eels_collection_angle")*mt::c_mrad_2_rad;
		mt::eChannelling_Type channelling_type = mx_get_scalar_field<mt::eChannelling_Type>(mx_input_multislice, "eels_channelling_type");
		int Z = mx_get_scalar_field<int>(mx_input_multislice, "eels_Z");

		input_multislice.eels_fr.set_input_data(space, input_multislice.E_0, E_loss, m_selection, collection_angle, channelling_type, Z);
	}
	else if (input_multislice.is_EFTEM())
	{
		mt::eSpace space = mt::eS_Real;
		value_type_r E_loss = mx_get_scalar_field<value_type_r>(mx_input_multislice, "eftem_E_loss")*mt::c_meV_2_keV;
		int m_selection = mx_get_scalar_field<int>(mx_input_multislice, "eftem_m_selection");
		value_type_r collection_angle = mx_get_scalar_field<double>(mx_input_multislice, "eftem_collection_angle")*mt::c_mrad_2_rad;
		mt::eChannelling_Type channelling_type = mx_get_scalar_field<mt::eChannelling_Type>(mx_input_multislice, "eftem_channelling_type");
		int Z = mx_get_scalar_field<int>(mx_input_multislice, "eftem_Z");

		input_multislice.eels_fr.set_input_data(space, input_multislice.E_0, E_loss, m_selection, collection_angle, channelling_type, Z);
	}

	input_multislice.validate_parameters();
 }

void set_output_multislice(const mxArray *mx_input_multislice, mxArray *&mx_output_multislice, mt::Output_Multislice_Matlab &output_multislice)
{
	mt::Input_Multislice<double> input_multislice;
	read_input_multislice(mx_input_multislice, input_multislice);
	output_multislice.set_input_data(&input_multislice);

	const char *field_names_output_multislice[] = {"dx", "dy", "x", "y", "thickness", "data"};
	int number_of_fields_output_multislice = 6;
	mwSize dims_output_multislice[2] = {1, 1};

	mx_output_multislice = mxCreateStructArray(2, dims_output_multislice, number_of_fields_output_multislice, field_names_output_multislice);

	mx_create_set_scalar_field<rmatrix_r>(mx_output_multislice, "dx", output_multislice.dx);
	mx_create_set_scalar_field<rmatrix_r>(mx_output_multislice, "dy", output_multislice.dy);
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "x", 1, output_multislice.x.size(), output_multislice.x.data());
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "y", 1, output_multislice.y.size(), output_multislice.y.data());
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "thickness", 1, output_multislice.thickness.size(), output_multislice.thickness.data());

	if(output_multislice.is_STEM() || output_multislice.is_EELS())
	{
		mxArray *mx_field_data;
		const char *field_names_data_full[] = {"image_tot", "image_coh"};
		const char *field_names_data_partial[] = {"image_tot"};
		const char **field_names_data = (output_multislice.coherent_contribution)?field_names_data_full:field_names_data_partial;
		int number_of_fields_data = (output_multislice.coherent_contribution)?2:1;
		mwSize dims_data[2] = {1, output_multislice.thickness.size()};

		mx_field_data = mxCreateStructArray(2, dims_data, number_of_fields_data, field_names_data);
		mxSetField(mx_output_multislice, 0, "data", mx_field_data);

		mxArray *mx_field_detector_tot;
		mxArray *mx_field_detector_coh;
		const char *field_names_detector[] = {"image"};
		int number_of_fields_detector = 1;
		// mwSize dims_detector[2] = {1, output_multislice.ndetector};
		mwSize dims_detector[2];
		dims_detector[0] = 1;
		dims_detector[1] = output_multislice.ndetector;

		for(auto ithk = 0; ithk<output_multislice.thickness.size(); ithk++)
		{
			mx_field_detector_tot = mxCreateStructArray(2, dims_detector, number_of_fields_detector, field_names_detector);
			mxSetField(mx_field_data, ithk, "image_tot", mx_field_detector_tot);
			if(output_multislice.coherent_contribution)
			{
				mx_field_detector_coh = mxCreateStructArray(2, dims_detector, number_of_fields_detector, field_names_detector);
				mxSetField(mx_field_data, ithk, "image_coh", mx_field_detector_coh);
			}

			for(auto iDet = 0; iDet<output_multislice.ndetector; iDet++)
			{
				output_multislice.image_tot[ithk].image[iDet] = mx_create_matrix_field<rmatrix_r>(mx_field_detector_tot, iDet, "image", output_multislice.ny, output_multislice.nx);
				if(output_multislice.coherent_contribution)
				{
					output_multislice.image_coh[ithk].image[iDet] = mx_create_matrix_field<rmatrix_r>(mx_field_detector_coh, iDet, "image", output_multislice.ny, output_multislice.nx);
				}
			}
		}
	}
	else if(output_multislice.is_EWFS_EWRS())
	{
		mxArray *mx_field_data;
		const char *field_names_data_full[] = {"m2psi_tot", "psi_coh"};
		const char *field_names_data_partial[] = {"psi_coh"};
		const char **field_names_data = (!output_multislice.is_EWFS_EWRS_SC())?field_names_data_full:field_names_data_partial;
		int number_of_fields_data = (!output_multislice.is_EWFS_EWRS_SC())?2:1;
		mwSize dims_data[2] = {1, output_multislice.thickness.size()};

		mx_field_data = mxCreateStructArray(2, dims_data, number_of_fields_data, field_names_data);
		mxSetField(mx_output_multislice, 0, "data", mx_field_data);

		for(auto ithk = 0; ithk<output_multislice.thickness.size(); ithk++)
		{
			if(!output_multislice.is_EWFS_EWRS_SC())
			{
				output_multislice.m2psi_tot[ithk] = mx_create_matrix_field<rmatrix_r>(mx_field_data, ithk, "m2psi_tot", output_multislice.ny, output_multislice.nx);
			}
			output_multislice.psi_coh[ithk] = mx_create_matrix_field<rmatrix_c>(mx_field_data, ithk, "psi_coh", output_multislice.ny, output_multislice.nx);
		}
	}
	else
	{
		mxArray *mx_field_data;
		const char *field_names_data_full[] = {"m2psi_tot", "m2psi_coh"};
		const char *field_names_data_partial[] = {"m2psi_tot"};
		const char **field_names_data = (output_multislice.coherent_contribution)?field_names_data_full:field_names_data_partial;
		int number_of_fields_data = (output_multislice.coherent_contribution)?2:1;
		mwSize dims_data[2] = {1, output_multislice.thickness.size()};

		mx_field_data = mxCreateStructArray(2, dims_data, number_of_fields_data, field_names_data);
		mxSetField(mx_output_multislice, 0, "data", mx_field_data);

		for(auto ithk = 0; ithk<output_multislice.thickness.size(); ithk++)
		{
			output_multislice.m2psi_tot[ithk] = mx_create_matrix_field<rmatrix_r>(mx_field_data, ithk, "m2psi_tot", output_multislice.ny, output_multislice.nx);
			if(output_multislice.coherent_contribution)
			{
				output_multislice.m2psi_coh[ithk] = mx_create_matrix_field<rmatrix_r>(mx_field_data, ithk, "m2psi_coh", output_multislice.ny, output_multislice.nx);
			}
		}
	}
}

template<class T, mt::eDevice dev>
void get_multislice(const mxArray *mxB, mt::Output_Multislice_Matlab &output_multislice)
{
	mt::Input_Multislice<T> input_multislice;
	read_input_multislice(mxB, input_multislice);

	mt::Multislice<T, dev> multislice;
	multislice.set_input_data(&input_multislice);

	multislice.run(output_multislice);

	multislice.cleanup();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mt::Output_Multislice_Matlab output_multislice;
	set_output_multislice(prhs[0], plhs[0], output_multislice);

	if(output_multislice.is_float_host())
	{
		get_multislice<float, mt::e_host>(prhs[0], output_multislice);
	}
	else if(output_multislice.is_double_host())
	{
		get_multislice<double, mt::e_host>(prhs[0], output_multislice);
	}
	if(output_multislice.is_float_device())
	{
		get_multislice<float, mt::e_device>(prhs[0], output_multislice);
	}
	else if(output_multislice.is_double_device())
	{
		get_multislice<double, mt::e_device>(prhs[0], output_multislice);
	}
}