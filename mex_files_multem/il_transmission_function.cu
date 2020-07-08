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

#include "types.cuh"
#include "matlab_types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "fft.cuh"
#include "atomic_data_mt.hpp"
#include "input_multislice.cuh"
#include "output_multislice.hpp"

#include "transmission_function.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::rmatrix_r;
using mt::rmatrix_c;

template <class TInput_Multislice>
void read_input_multislice(const mxArray *mx_input_multislice, TInput_Multislice &input_multislice, bool full = true)
{
	using T_r = mt::Value_type<TInput_Multislice>;

	/************************ simulation type **************************/
	input_multislice.simulation_type = mt::eTEMST_TFRS;

	/*******************************************************************/
	input_multislice.interaction_model = mx_get_scalar_field<mt::eElec_Spec_Int_Model>(mx_input_multislice, "interaction_model");
	input_multislice.potential_type = mx_get_scalar_field<mt::ePotential_Type>(mx_input_multislice, "potential_type");

	/************** Electron-Phonon interaction model ******************/
	input_multislice.pn_model = mx_get_scalar_field<mt::ePhonon_Model>(mx_input_multislice, "pn_model"); 
	input_multislice.pn_coh_contrib = true;
	input_multislice.pn_single_conf = true;
	input_multislice.pn_nconf = mx_get_scalar_field<int>(mx_input_multislice, "pn_nconf");
	input_multislice.pn_dim.set(mx_get_scalar_field<int>(mx_input_multislice, "pn_dim"));
	input_multislice.pn_seed = mx_get_scalar_field<int>(mx_input_multislice, "pn_seed");

	/**************************** Specimen *****************************/
	auto atoms = mx_get_matrix_field<rmatrix_r>(mx_input_multislice, "spec_atoms");

	auto lx = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_lx");
	auto ly = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_ly");
	auto lz = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_lz");
	auto dz = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_dz");
	bool pbc_xy = true; 

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
	for(auto i = 0; i<amorp_lay_info.size(); i++)
	{
		amorp_lay_info[i].z_0 = mx_get_scalar_field<T_r>(mx_spec_amorp, i, "z_0");
		amorp_lay_info[i].z_e = mx_get_scalar_field<T_r>(mx_spec_amorp, i, "z_e");
		amorp_lay_info[i].dz = mx_get_scalar_field<T_r>(mx_spec_amorp, i, "dz");
	}

	if(full)
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

	/************************ Potential slicing ************************/
	input_multislice.potential_slicing = mx_get_scalar_field<mt::ePotential_Slicing>(mx_input_multislice, "potential_slicing");

	/************************** xy sampling ****************************/
	auto nx = mx_get_scalar_field<int>(mx_input_multislice, "nx");
	auto ny = mx_get_scalar_field<int>(mx_input_multislice, "ny");
	bool bwl = mx_get_scalar_field<bool>(mx_input_multislice, "bwl");

	input_multislice.grid_2d.set_input_data(nx, ny, lx, ly, dz, bwl, pbc_xy);

	/********************* Microscope parameter ***********************/
	input_multislice.E_0 = mx_get_scalar_field<T_r>(mx_input_multislice, "E_0");
	input_multislice.theta = mx_get_scalar_field<T_r>(mx_input_multislice, "theta")*mt::c_deg_2_rad;
	input_multislice.phi = mx_get_scalar_field<T_r>(mx_input_multislice, "phi")*mt::c_deg_2_rad;

	input_multislice.islice = mx_get_scalar_field<int>(mx_input_multislice, "islice")-1;

	/********************* select output region *************************/
	input_multislice.output_area.ix_0 = mx_get_scalar_field<int>(mx_input_multislice, "output_area_ix_0")-1;
	input_multislice.output_area.iy_0 = mx_get_scalar_field<int>(mx_input_multislice, "output_area_iy_0")-1;
	input_multislice.output_area.ix_e = mx_get_scalar_field<int>(mx_input_multislice, "output_area_ix_e")-1;
	input_multislice.output_area.iy_e = mx_get_scalar_field<int>(mx_input_multislice, "output_area_iy_e")-1;

	/********************* validate parameters *************************/
	input_multislice.validate_parameters();
}

template<class TOutput_Multislice>
void set_struct_transmission_function(TOutput_Multislice &output_multislice, mxArray *&mx_output_multislice)
{
	const char *field_names_output_multislice[] = {"dx", "dy", "x", "y", "thick", "trans"};
	int number_of_fields_output_multislice = 6;
	mwSize dims_output_multislice[2] = {1, 1};

	mx_output_multislice = mxCreateStructArray(2, dims_output_multislice, number_of_fields_output_multislice, field_names_output_multislice);

	mx_create_set_scalar_field<rmatrix_r>(mx_output_multislice, 0, "dx", output_multislice.dx);
	mx_create_set_scalar_field<rmatrix_r>(mx_output_multislice, 0, "dy", output_multislice.dy);
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "x", 1, output_multislice.x.size(), output_multislice.x);
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "y", 1, output_multislice.y.size(), output_multislice.y);
	mx_create_set_matrix_field<rmatrix_r>(mx_output_multislice, "thick", 1, output_multislice.thick.size(), output_multislice.thick);
	mx_create_set_matrix_field<rmatrix_c>(mx_output_multislice, "trans", output_multislice.ny, output_multislice.nx, output_multislice.trans[0]);
}

template <class T, mt::eDevice dev>
void run_transmission_function(mt::System_Configuration &system_conf, const mxArray *mx_input_multislice, mxArray *&mx_output_multislice)
{
	mt::Input_Multislice<T> input_multislice;
	read_input_multislice(mx_input_multislice, input_multislice);
	input_multislice.system_conf = system_conf;

	mt::Stream<dev> stream(system_conf.nstream);
	mt::FFT<T, dev> fft_2d;
	fft_2d.create_plan_2d(input_multislice.grid_2d.ny, input_multislice.grid_2d.nx, system_conf.nstream);
	
	mt::Transmission_Function<T, dev> transmission_function;
	transmission_function.set_input_data(&input_multislice, &stream, &fft_2d);

	mt::Output_Multislice<T> output_multislice;
	output_multislice.set_input_data(&input_multislice);

	transmission_function.move_atoms(input_multislice.pn_nconf);
	transmission_function.trans(input_multislice.islice, output_multislice);

	stream.synchronize();

	output_multislice.gather();
	output_multislice.clean_temporal();
	fft_2d.cleanup();

	set_struct_transmission_function(output_multislice, mx_output_multislice);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	auto system_conf = mt::read_system_conf(prhs[0]);
	int idx_0 = (system_conf.active)?1:0;

	if(system_conf.is_float_host())
	{
		run_transmission_function<float, mt::e_host>(system_conf, prhs[idx_0], plhs[0]);
	}
	else if(system_conf.is_double_host())
	{
		run_transmission_function<double, mt::e_host>(system_conf, prhs[idx_0], plhs[0]);
	}
	else if(system_conf.is_float_device())
	{
		run_transmission_function<float, mt::e_device>(system_conf, prhs[idx_0], plhs[0]);
	}
	else if(system_conf.is_double_device())
	{
		run_transmission_function<double, mt::e_device>(system_conf, prhs[idx_0], plhs[0]);
	}
}