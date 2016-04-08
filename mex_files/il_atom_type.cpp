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

#include <vector>
#include "types.cuh"
#include "atomic_data.hpp"

#include <mex.h>
#include "matlab_mex.cuh"

using multem::rmatrix_r;
using multem::e_host;
using multem::Atom_Type;

template<class TAtom_Types>
void set_output_data(TAtom_Types &atom_type, mxArray *&mx_atom_type)
{
	const char *field_names[] = {"Z", "m", "A", "rn_e", "rn_c", "ra_e", "ra_c", "coef"};
	int number_of_fields = 8;
	mwSize dims[2] = {1, atom_type.size()};

	mxArray *mx_field_Coef=nullptr;
	const char *field_names_Coef[] = {"charge", "R_min", "R_max", "feg", "fxg", "Pr", "Vr", "VR", "R2", "ciVR"};
	int number_of_fields_Coef = 10;
	mwSize dims_Coef[2] = {1, 1};

	mxArray *mx_field_CoefPar=nullptr;
	const char *field_names_CoefPar[] = {"cl", "cnl"};
	int number_of_fields_CoefPar = 2;
	mwSize dims_CoefPar[2] = {1, 1};

	mxArray *mxfield_ciVR=nullptr;
	const char *field_names_ciVR[] = {"c0", "c1", "c2", "c3"};
	int number_of_fields_ciVR = 4;
	mwSize dims_ciVR[2] = {1, 1};

	mx_atom_type = mxCreateStructArray(2, dims, number_of_fields, field_names);
	for(auto i = 0; i<atom_type.size(); i++)
	{
		mx_create_set_scalar_field<rmatrix_r>(mx_atom_type, i, "Z", atom_type[i].Z);
		mx_create_set_scalar_field<rmatrix_r>(mx_atom_type, i, "m", atom_type[i].m);
		mx_create_set_scalar_field<rmatrix_r>(mx_atom_type, i, "A", atom_type[i].A);
		mx_create_set_scalar_field<rmatrix_r>(mx_atom_type, i, "rn_e", atom_type[i].rn_e);
		mx_create_set_scalar_field<rmatrix_r>(mx_atom_type, i, "rn_c", atom_type[i].rn_c);
		mx_create_set_scalar_field<rmatrix_r>(mx_atom_type, i, "ra_e", atom_type[i].ra_e);
		mx_create_set_scalar_field<rmatrix_r>(mx_atom_type, i, "ra_c", atom_type[i].ra_c);


		mwSize dims_Coef[2] = {1, atom_type[i].coef.size()};
		mx_field_Coef = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mx_atom_type, i, "coef", mx_field_Coef);

		for(auto j= 0; j<atom_type[i].coef.size(); j++)
		{
			auto &coef = atom_type[i].coef[j];
			mx_create_set_scalar_field<rmatrix_r>(mx_field_Coef, j, "charge", coef.charge);
			mx_create_set_scalar_field<rmatrix_r>(mx_field_Coef, j, "R_min", coef.R_min);
			mx_create_set_scalar_field<rmatrix_r>(mx_field_Coef, j, "R_max", coef.R_max);

			/*************************fg***************************/
			mx_field_CoefPar = mxCreateStructArray(2, dims_CoefPar, number_of_fields_CoefPar, field_names_CoefPar);
			mxSetField(mx_field_Coef, j, "feg", mx_field_CoefPar);
			mx_create_set_matrix_field<rmatrix_r>(mx_field_CoefPar, "cl", 1, coef.feg.size(), coef.feg.cl.data());
			mx_create_set_matrix_field<rmatrix_r>(mx_field_CoefPar, "cnl", 1, coef.feg.size(), coef.feg.cnl.data());

			/*************************fx***************************/
			mx_field_CoefPar = mxCreateStructArray(2, dims_CoefPar, number_of_fields_CoefPar, field_names_CoefPar);
			mxSetField(mx_field_Coef, j, "fxg", mx_field_CoefPar);
			mx_create_set_matrix_field<rmatrix_r>(mx_field_CoefPar, "cl", 1, coef.fxg.size(), coef.fxg.cl.data());
			mx_create_set_matrix_field<rmatrix_r>(mx_field_CoefPar, "cnl", 1, coef.fxg.size(), coef.fxg.cnl.data());

			/*************************Pr***************************/
			mx_field_CoefPar = mxCreateStructArray(2, dims_CoefPar, number_of_fields_CoefPar, field_names_CoefPar);
			mxSetField(mx_field_Coef, j, "Pr", mx_field_CoefPar);
			mx_create_set_matrix_field<rmatrix_r>(mx_field_CoefPar, "cl", 1, coef.Pr.size(), coef.Pr.cl.data());
			mx_create_set_matrix_field<rmatrix_r>(mx_field_CoefPar, "cnl", 1, coef.Pr.size(), coef.Pr.cnl.data());

			/*************************Vr***************************/
			mx_field_CoefPar = mxCreateStructArray(2, dims_CoefPar, number_of_fields_CoefPar, field_names_CoefPar);
			mxSetField(mx_field_Coef, j, "Vr", mx_field_CoefPar);
			mx_create_set_matrix_field<rmatrix_r>(mx_field_CoefPar, "cl", 1, coef.Vr.size(), coef.Vr.cl.data());
			mx_create_set_matrix_field<rmatrix_r>(mx_field_CoefPar, "cnl", 1, coef.Vr.size(), coef.Vr.cnl.data());

			/*************************VR***************************/
			mx_field_CoefPar = mxCreateStructArray(2, dims_CoefPar, number_of_fields_CoefPar, field_names_CoefPar);
			mxSetField(mx_field_Coef, j, "VR", mx_field_CoefPar);
			mx_create_set_matrix_field<rmatrix_r>(mx_field_CoefPar, "cl", 1, coef.VR.size(), coef.VR.cl.data());
			mx_create_set_matrix_field<rmatrix_r>(mx_field_CoefPar, "cnl", 1, coef.VR.size(), coef.VR.cnl.data());

			/*************************ciVR***************************/
			mx_create_set_matrix_field<rmatrix_r>(mx_field_Coef, j, "R2", 1, coef.R2.size(), coef.R2.data());
			mxfield_ciVR = mxCreateStructArray(2, dims_ciVR, number_of_fields_ciVR, field_names_ciVR);
			mxSetField(mx_field_Coef, j, "ciVR", mxfield_ciVR);
			mx_create_set_matrix_field<rmatrix_r>(mxfield_ciVR, "c0", 1, coef.ciVR.size(), coef.ciVR.c0.data());
			mx_create_set_matrix_field<rmatrix_r>(mxfield_ciVR, "c1", 1, coef.ciVR.size(), coef.ciVR.c1.data());
			mx_create_set_matrix_field<rmatrix_r>(mxfield_ciVR, "c2", 1, coef.ciVR.size(), coef.ciVR.c2.data());
			mx_create_set_matrix_field<rmatrix_r>(mxfield_ciVR, "c3", 1, coef.ciVR.size(), coef.ciVR.c3.data());
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	auto potential_type = mx_get_scalar<multem::ePotential_Type>(prhs[0]);

	multem::Atomic_Data atomic_data(potential_type);
	multem::Vector<Atom_Type<double, e_host>, e_host> atom_type(multem::c_nAtomsTypes);

	for(auto i = 0; i<atom_type.size(); i++)
	{
		atomic_data.To_atom_type_CPU(i+1, multem::c_Vrl, multem::c_nR, 0.0, atom_type[i]);
	}

	set_output_data(atom_type, plhs[0]);
}