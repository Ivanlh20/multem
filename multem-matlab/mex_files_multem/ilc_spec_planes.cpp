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
#include "atomic_data_mt.hpp"
#include "input_multislice.cuh"
#include "slicing.hpp"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::rmatrix_r;

template <class TInput_Multislice>
void read_input_multislice(const mxArray *mx_input_multislice, TInput_Multislice &input_multislice)
{
	using T_r = mt::Value_type<TInput_Multislice>;

	input_multislice.interaction_model = mx_get_scalar_field<mt::eElec_Spec_Int_Model>(mx_input_multislice, "interaction_model");
	input_multislice.potential_type = mt::ePT_Lobato_0_12;

	/************** Electron-Phonon interaction model ******************/
	input_multislice.pn_model = mx_get_scalar_field<mt::ePhonon_Model>(mx_input_multislice, "pn_model"); 
	input_multislice.pn_coh_contrib = false;
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
	bool pbc_xy = false; 

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

	input_multislice.atoms.set_crystal_parameters(ct_na, ct_nb, ct_nc, ct_a, ct_b, ct_c, ct_x0, ct_y0);
	input_multislice.atoms.set_amorphous_parameters(amorp_lay_info);
	input_multislice.atoms.set_atoms(atoms.rows, atoms.cols, atoms.real, lx, ly, lz, dz);

	/************************ Specimen rotation *************************/
	input_multislice.spec_rot_theta = mx_get_scalar_field<T_r>(mx_input_multislice, "spec_rot_theta")*mt::c_deg_2_rad;
	input_multislice.spec_rot_u0 = mx_get_r3d_field<T_r>(mx_input_multislice, "spec_rot_u0");
	input_multislice.spec_rot_u0.normalized();
 input_multislice.spec_rot_center_type = mx_get_scalar_field<mt::eRot_Point_Type>(mx_input_multislice, "spec_rot_center_type");
	input_multislice.spec_rot_center_p = mx_get_r3d_field<T_r>(mx_input_multislice, "spec_rot_center_p");

	/************************ Potential slicing ************************/
	input_multislice.potential_slicing = mx_get_scalar_field<mt::ePotential_Slicing>(mx_input_multislice, "potential_slicing");

	/************************** xy sampling ****************************/
	auto nx = 1024;
	auto ny = 1024;
	bool bwl = false;

	input_multislice.grid_2d.set_input_data(nx, ny, lx, ly, dz, bwl, pbc_xy);

	input_multislice.validate_parameters();
 }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
	/*************************Input data**************************/
	mt::Input_Multislice<double> input_multislice;
	read_input_multislice(prhs[0], input_multislice);

	/***************************************************************************/
	mt::Slicing<double> slicing;
	slicing.set_input_data(&input_multislice, &(input_multislice.atoms));

	// /************************Output data**************************/
	auto rplanes = mx_create_matrix<rmatrix_r>(slicing.z_plane.size(), 1, plhs[0]);

	rplanes.assign(slicing.z_plane.begin(), slicing.z_plane.end());
}