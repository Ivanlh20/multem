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

#include "types.cuh"
#include "matlab_types.cuh"
#include "traits.cuh"
#include "atom_data.hpp"
#include "specimen.hpp"
#include "host_device_functions.cuh"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "input_multislice.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using multem::rmatrix_r;

template<class TInput_Multislice>
void read_input_multislice(const mxArray *mx_input_multislice, TInput_Multislice &input_multislice)
{
	using value_type_r = multem::Value_type<TInput_Multislice>;

	input_multislice.phonon_model = mx_get_scalar_field<multem::ePhonon_Model>(mx_input_multislice, "phonon_model"); 
	input_multislice.interaction_model = mx_get_scalar_field<multem::eElec_Spec_Int_Model>(mx_input_multislice, "interaction_model");
	input_multislice.potential_slicing = mx_get_scalar_field<multem::ePotential_Slicing>(mx_input_multislice, "potential_slicing");
	input_multislice.fp_dim.set(mx_get_scalar_field<int>(mx_input_multislice, "fp_dim"));
	input_multislice.fp_seed = mx_get_scalar_field<int>(mx_input_multislice, "fp_seed");
	input_multislice.fp_single_conf = true;
	input_multislice.fp_nconf = mx_get_scalar_field<int>(mx_input_multislice, "fp_nconf");

	input_multislice.tm_active = mx_get_scalar_field<bool>(mx_input_multislice, "tm_active");
	input_multislice.tm_theta = mx_get_scalar_field<value_type_r>(mx_input_multislice, "tm_theta")*multem::c_deg_2_rad;
	input_multislice.tm_u0 = mx_get_r3d_field<value_type_r>(mx_input_multislice, "tm_u0");
	input_multislice.tm_rot_point_type = mx_get_scalar_field<multem::eRot_Point_Type>(mx_input_multislice, "tm_rot_point_type");
	input_multislice.tm_p0 = mx_get_r3d_field<value_type_r>(mx_input_multislice, "tm_p0");

	int nx = 1024;
	int ny = 1024;
	auto lx = mx_get_scalar_field<value_type_r>(mx_input_multislice, "lx");
	auto ly = mx_get_scalar_field<value_type_r>(mx_input_multislice, "ly");
	auto dz = mx_get_scalar_field<value_type_r>(mx_input_multislice, "dz");
	bool bwl = false;
	bool pbc_xy = false; 

	auto atoms = mx_get_matrix_field<rmatrix_r>(mx_input_multislice, "atoms");

	input_multislice.atoms.set_Atoms(atoms.rows, atoms.cols, atoms.real, lx, ly);
	input_multislice.grid.set_input_data(nx, ny, lx, ly, dz, bwl, pbc_xy);
	input_multislice.validate_parameters();
 }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
	/*************************Input data**************************/
	multem::Input_Multislice<double> input_multislice;
	read_input_multislice(prhs[0], input_multislice);

	multem::Specimen<double> specimen;
	specimen.set_input_data(&input_multislice);
	specimen.move_atoms(input_multislice.fp_nconf);

	// /************************Output data**************************/
	auto atomsM = mx_create_matrix<rmatrix_r>(specimen.atoms.size(), 7, plhs[0]);
	auto sliceM = mx_create_matrix<rmatrix_r>(specimen.slice.size(), 7, plhs[1]);

	for(auto i = 0; i<atomsM.rows; i++)
	{
		atomsM.real[0*atomsM.rows+i] = specimen.atoms.Z[i];
		atomsM.real[1*atomsM.rows+i] = specimen.atoms.x[i];
		atomsM.real[2*atomsM.rows+i] = specimen.atoms.y[i];
		atomsM.real[3*atomsM.rows+i] = specimen.atoms.z[i];
		atomsM.real[4*atomsM.rows+i] = specimen.atoms.sigma[i];
		atomsM.real[5*atomsM.rows+i] = specimen.atoms.occ[i];
		atomsM.real[6*atomsM.rows+i] = specimen.atoms.charge[i];
	}

	for(auto i = 0; i<sliceM.rows; i++)
	{
		sliceM.real[0*sliceM.rows+i] = specimen.slice.z_0[i];
		sliceM.real[1*sliceM.rows+i] = specimen.slice.z_e[i];
		sliceM.real[2*sliceM.rows+i] = specimen.slice.z_int_0[i];
		sliceM.real[3*sliceM.rows+i] = specimen.slice.z_int_e[i];
		sliceM.real[4*sliceM.rows+i] = specimen.slice.iatom_0[i]+1;
		sliceM.real[5*sliceM.rows+i] = specimen.slice.iatom_e[i]+1;
	}
}