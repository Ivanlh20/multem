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

#include <algorithm>

#include "math.cuh"
#include "types.hpp"
#include "input_multislice.hpp"
#include "atom_data.hpp"
#include "multislice.cuh"

#include <mex.h>
#include "matlab2cpp.hpp"

using multem::m_matrix_r;
using multem::m_matrix_c;

template<class TInput_Multislice>
void read_input_data(const mxArray *mx_input_multislice, TInput_Multislice &input_multislice, bool full=true)
{
	using value_type = multem::traits::Value_type<TInput_Multislice>;

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
	input_multislice.fp_nconf = mx_get_scalar_field<int>(mx_input_multislice, "fp_nconf");
	input_multislice.fp_iconf = mx_get_scalar_field<int>(mx_input_multislice, "fp_iconf");

	input_multislice.microscope_effect = mx_get_scalar_field<multem::eMicroscope_Effect>(mx_input_multislice, "microscope_effect");
	input_multislice.spatial_temporal_effect = mx_get_scalar_field<multem::eSpatial_Temporal_Effect>(mx_input_multislice, "spatial_temporal_effect");

	input_multislice.zero_defocus_type = mx_get_scalar_field<multem::eZero_Defocus_Type>(mx_input_multislice, "zero_defocus_type");
	input_multislice.zero_defocus_plane = mx_get_scalar_field<value_type>(mx_input_multislice, "zero_defocus_plane");

	input_multislice.thickness_type = mx_get_scalar_field<multem::eThickness_Type>(mx_input_multislice, "thickness_type");
	if(input_multislice.thickness_type!=multem::eTT_Whole_Specimen && full)
	{
		auto thickness = mx_get_matrix_field<m_matrix_r>(mx_input_multislice, "thickness");
		input_multislice.thickness.resize(thickness.size);
		std::copy(thickness.real, thickness.real + thickness.size, input_multislice.thickness.begin());
	}

	input_multislice.input_wave_type = mx_get_scalar_field<multem::eInput_Wave_Type>(mx_input_multislice, "input_wave_type");
	if(input_multislice.is_user_define_wave() && full)
	{
		auto psi_0 = mx_get_matrix_field<m_matrix_c>(mx_input_multislice, "psi_0");
		input_multislice.psi_0.resize(psi_0.size);
		multem::scomplex_to_complex(psi_0, input_multislice.psi_0);
		multem::fft2_shift(input_multislice.grid, input_multislice.psi_0);
	}

	input_multislice.operation_mode = mx_get_scalar_field<multem::eOperation_Mode>(mx_input_multislice, "operation_mode");
	input_multislice.coherent_contribution = mx_get_scalar_field<bool>(mx_input_multislice, "coherent_contribution");

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

	auto atoms = mx_get_matrix_field<m_matrix_r>(mx_input_multislice, "atoms");
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

	if(input_multislice.is_scanning())
	{
		input_multislice.scanning.type = mx_get_scalar_field<multem::eScanning_Type>(mx_input_multislice, "scanning_type");
		input_multislice.scanning.ns = mx_get_scalar_field<int>(mx_input_multislice, "scanning_ns");
		input_multislice.scanning.x0 = mx_get_scalar_field<value_type>(mx_input_multislice, "scanning_x0");
		input_multislice.scanning.y0 = mx_get_scalar_field<value_type>(mx_input_multislice, "scanning_y0");
		input_multislice.scanning.xe = mx_get_scalar_field<value_type>(mx_input_multislice, "scanning_xe");
		input_multislice.scanning.ye = mx_get_scalar_field<value_type>(mx_input_multislice, "scanning_ye");
		input_multislice.scanning.set_grid();
	}

	if(input_multislice.is_STEM())
	{
		mxArray *mx_det_cir = mxGetField(mx_input_multislice, 0, "det_cir");
		int ndet_cir = mxGetN(mx_det_cir);
		if(ndet_cir>0)
		{
			input_multislice.det_cir.resize(ndet_cir);
			for(auto i=0; i<input_multislice.det_cir.size(); i++)
			{
				input_multislice.det_cir.ang_inner[i] = mx_get_scalar_field<value_type>(mx_det_cir, i, "ang_inner")*multem::c_mrad_2_rad; // Inner angle(mrad-->rad)
				input_multislice.det_cir.ang_outer[i] = mx_get_scalar_field<value_type>(mx_det_cir, i, "ang_outer")*multem::c_mrad_2_rad; // Outer angle(mrad-->rad)
			}
		}
	}
	else if (input_multislice.is_CBED())
	{
		input_multislice.cbe_fr.space = multem::eS_Reciprocal;
		input_multislice.cbe_fr.x0 = mx_get_scalar_field<value_type>(mx_input_multislice, "cbed_x0");
		input_multislice.cbe_fr.y0 = mx_get_scalar_field<value_type>(mx_input_multislice, "cbed_y0");
	}
	else if (input_multislice.is_CBEI())
	{
		input_multislice.cbe_fr.space = multem::eS_Real;
		input_multislice.cbe_fr.x0 = mx_get_scalar_field<value_type>(mx_input_multislice, "cbei_x0");
		input_multislice.cbe_fr.y0 = mx_get_scalar_field<value_type>(mx_input_multislice, "cbei_y0");
	}
	else if (input_multislice.is_PED())
	{
		input_multislice.pe_fr.space = multem::eS_Reciprocal;
		input_multislice.pe_fr.nrot = mx_get_scalar_field<value_type>(mx_input_multislice, "ped_nrot");
		input_multislice.pe_fr.theta = mx_get_scalar_field<value_type>(mx_input_multislice, "ped_theta")*multem::c_deg_2_rad;
	}
	else if (input_multislice.is_HCI())
	{
		input_multislice.pe_fr.space = multem::eS_Real;
		input_multislice.pe_fr.nrot = mx_get_scalar_field<value_type>(mx_input_multislice, "hci_nrot");
		input_multislice.pe_fr.theta = mx_get_scalar_field<value_type>(mx_input_multislice, "hci_theta")*multem::c_deg_2_rad;
	}
	else if (input_multislice.is_EWFS())
	{
		input_multislice.ew_fr.space = multem::eS_Reciprocal;
		input_multislice.ew_fr.convergent_beam = mx_get_scalar_field<bool>(mx_input_multislice, "ewfs_convergent_beam");
		input_multislice.ew_fr.y0 = mx_get_scalar_field<value_type>(mx_input_multislice, "ewfs_x0");
		input_multislice.ew_fr.y0 = mx_get_scalar_field<value_type>(mx_input_multislice, "ewfs_y0");
	}
	else if (input_multislice.is_EWRS())
	{
		input_multislice.ew_fr.space = multem::eS_Real;
		input_multislice.ew_fr.convergent_beam = mx_get_scalar_field<bool>(mx_input_multislice, "ewrs_convergent_beam");
		input_multislice.ew_fr.y0 = mx_get_scalar_field<value_type>(mx_input_multislice, "ewrs_x0");
		input_multislice.ew_fr.y0 = mx_get_scalar_field<value_type>(mx_input_multislice, "ewrs_y0");
	}
	else if (input_multislice.is_EELS())
	{
		multem::eSpace space = multem::eS_Reciprocal;
		value_type E_loss = mx_get_scalar_field<value_type>(mx_input_multislice, "eels_E_loss")*multem::c_meV_2_keV;
		int m_selection = mx_get_scalar_field<int>(mx_input_multislice, "eels_m_selection");
		value_type collection_angle = mx_get_scalar_field<int>(mx_input_multislice, "eels_collection_angle")*multem::c_mrad_2_rad;
		multem::eChannelling_Type channelling_type = mx_get_scalar_field<multem::eChannelling_Type>(mx_input_multislice, "eels_channelling_type");
		int Z = mx_get_scalar_field<int>(mx_input_multislice, "eels_Z");

		input_multislice.eels_fr.set_input_data(space, input_multislice.E_0, E_loss, m_selection, collection_angle, channelling_type, Z);
	}
	else if (input_multislice.is_EFTEM())
	{
		multem::eSpace space = multem::eS_Real;
		value_type E_loss = mx_get_scalar_field<value_type>(mx_input_multislice, "eftem_E_loss")*multem::c_meV_2_keV;
		int m_selection = mx_get_scalar_field<int>(mx_input_multislice, "eftem_m_selection");
		value_type collection_angle = 0;
		multem::eChannelling_Type channelling_type = mx_get_scalar_field<multem::eChannelling_Type>(mx_input_multislice, "eftem_channelling_type");
		int Z = mx_get_scalar_field<int>(mx_input_multislice, "eftem_Z");

		input_multislice.eels_fr.set_input_data(space, input_multislice.E_0, E_loss, m_selection, collection_angle, channelling_type, Z);
	}

	input_multislice.validate_parameters();
 }

template<class TInput_Multislice>
void set_output_data(const TInput_Multislice &input_multislice, mxArray *&mx_plhs0, mxArray *&mx_plhs1)
{
	if(input_multislice.is_STEM())
	{
		const char *field_names_det_int[] = {"image"};
		int number_of_fields_det_int = 1;
		mwSize dims_det_int[2] = {1, input_multislice.det_cir.size()};

		mx_plhs0 = mxCreateStructArray(2, dims_det_int, number_of_fields_det_int, field_names_det_int);

		if(input_multislice.coherent_contribution)
		{
			mx_plhs1 = mxCreateStructArray(2, dims_det_int, number_of_fields_det_int, field_names_det_int);
		}

		int rows = (input_multislice.scanning.is_line())?1:input_multislice.scanning.ny;

		for(auto iDet=0; iDet<input_multislice.det_cir.size(); iDet++)
		{
			mx_create_matrix_field<m_matrix_r>(mx_plhs0, iDet, "image", rows, input_multislice.scanning.nx);
			if(input_multislice.coherent_contribution)
			{
				mx_create_matrix_field<m_matrix_r>(mx_plhs1, iDet, "image", rows, input_multislice.scanning.nx);
			}
		}
	}
	else if(input_multislice.is_EWFS_EWRS())
	{
		mx_create_matrix<m_matrix_r>(input_multislice.grid, mx_plhs0);
		if(input_multislice.coherent_contribution)
		{
			mx_create_matrix<m_matrix_c>(input_multislice.grid, mx_plhs1);
		}
	}
	else if(input_multislice.is_EELS())
	{
		int rows = (input_multislice.scanning.is_line())?1:input_multislice.scanning.ny;
		mx_create_matrix<m_matrix_r>(rows, input_multislice.scanning.nx, mx_plhs0);
	}
	else
	{
		mx_create_matrix<m_matrix_r>(input_multislice.grid, mx_plhs0);
		if(input_multislice.coherent_contribution)
		{
			mx_create_matrix<m_matrix_r>(input_multislice.grid, mx_plhs1);
		}
	}
}

template<class T, multem::eDevice dev, class TVector_1, class TVector_2>
void get_multislice(const mxArray *mxB, TVector_1 &tot, TVector_2 &coh)
{
	multem::Input_Multislice<T, dev> input_multislice;
	read_input_data(mxB, input_multislice);

	multem::Multislice<T, dev> multislice;
	multislice.set_input_data(&input_multislice);

	if(input_multislice.coherent_contribution)
	{
		multislice.output_matlab(&tot, &coh);
	}
	else
	{
		multislice.output_matlab(&tot);
	}

	multislice.cleanup();
}

template<class TInput_Multislice, class TVector_1, class TVector_2>
void get_simulation(const mxArray *mxB, const TInput_Multislice &input_multislice, TVector_1 &tot, TVector_2 &coh)
{
	if(input_multislice.is_float_Host())
	{
		get_multislice<float, multem::e_Host>(mxB, tot, coh);
	}
	else if(input_multislice.is_double_Host())
	{
		get_multislice<double, multem::e_Host>(mxB, tot, coh);
	}
	if(input_multislice.is_float_Device())
	{
		get_multislice<float, multem::e_Device>(mxB, tot, coh);
		cudaDeviceReset();
	}
	else if(input_multislice.is_double_Device())
	{
		get_multislice<double, multem::e_Device>(mxB, tot, coh);
		cudaDeviceReset();
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	multem::Input_Multislice<double, multem::e_Host> input_multislice;
	read_input_data(prhs[0], input_multislice, false);

	set_output_data(input_multislice, plhs[0], plhs[1]);

	if(input_multislice.is_STEM())
	{
		multem::Det_Int<double, multem::e_Host> det_int_coh;
		multem::Det_Int<double, multem::e_Host> det_int_tot;
		det_int_tot.resize(input_multislice.det_cir.size(), input_multislice.scanning.size());
		if(input_multislice.coherent_contribution)
		{
			det_int_coh.resize(input_multislice.det_cir.size(), input_multislice.scanning.size());
		}

		get_simulation(prhs[0], input_multislice, det_int_tot, det_int_coh);

		for(auto iDet=0; iDet<input_multislice.det_cir.size(); iDet++)
		{
			m_matrix_r m2psi_tot = mx_get_matrix_field<m_matrix_r>(plhs[0], iDet, "image");
			std::copy(det_int_tot.image[iDet].begin(), det_int_tot.image[iDet].end(), m2psi_tot.real);
			if(input_multislice.coherent_contribution)
			{
				m_matrix_r m2psi_coh = mx_get_matrix_field<m_matrix_r>(plhs[1], iDet, "image");
				std::copy(det_int_coh.image[iDet].begin(), det_int_coh.image[iDet].end(), m2psi_coh.real);
			}
		}
	}
	else if(input_multislice.is_EWFS_EWRS())
	{
		m_matrix_r m2psi_tot = mx_get_matrix<m_matrix_r>(plhs[0]);
		m_matrix_c psi_coh;
		if(input_multislice.coherent_contribution)
		{
			psi_coh = mx_get_matrix<m_matrix_c>(plhs[1]);
		}

		get_simulation(prhs[0], input_multislice, m2psi_tot, psi_coh);
	}
	else if(input_multislice.is_EELS())
	{
		m_matrix_r m2psi_tot = mx_get_matrix<m_matrix_r>(plhs[0]);

		multem::Vector<double, multem::e_Host> v_m2psi_tot(m2psi_tot.size);
		multem::Vector<double, multem::e_Host> v_m2psi_coh;
		get_simulation(prhs[0], input_multislice, v_m2psi_tot, v_m2psi_coh);

		std::copy(v_m2psi_tot.begin(), v_m2psi_tot.end(), m2psi_tot.real);
	}
	else
	{
		m_matrix_r m2psi_tot = mx_get_matrix<m_matrix_r>(plhs[0]);
		m_matrix_r m2psi_coh;
		if(input_multislice.coherent_contribution)
		{
			m2psi_coh = mx_get_matrix<m_matrix_r>(plhs[1]);
		}
		get_simulation(prhs[0], input_multislice, m2psi_tot, m2psi_coh);
	}
}