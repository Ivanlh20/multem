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

#include <vector>
#include "types.hpp"
#include "atomic_data.hpp"

#include <mex.h>
#include "matlab2cpp.hpp"

using multem::m_matrix_r;

template<class TAtom_Types>
void set_output_data(TAtom_Types &atom_type, mxArray *&mx_atom_type)
{
	const char *field_names[] = {"Z", "m", "A", "rn_e", "rn_c", "ra_e", "ra_c", "R_min", "R_max", "feg", "fxg", "Pr", "Vr", "VR", "R2", "ciVR"};
	int number_of_fields = 16;
	mwSize dims[2] = {atom_type.size(), 1};

	mxArray *mxfield_CoefPar;
	const char *field_names_CoefPar[] = {"cl", "cnl"};
	int number_of_fields_CoefPar = 2;
	mwSize dims_CoefPar[2] = {1, 1};

	mxArray *mxfield_ciVR;
	const char *field_names_ciVR[] = {"c0", "c1", "c2", "c3"};
	int number_of_fields_ciVR = 4;
	mwSize dims_ciVR[2] = {1, 1};

	mx_atom_type = mxCreateStructArray(2, dims, number_of_fields, field_names);
	for(auto i=0; i<atom_type.size(); i++)
	{
		mx_create_set_scalar_field<m_matrix_r>(mx_atom_type, i, "Z", atom_type[i].Z);
		mx_create_set_scalar_field<m_matrix_r>(mx_atom_type, i, "m", atom_type[i].m);
		mx_create_set_scalar_field<m_matrix_r>(mx_atom_type, i, "A", atom_type[i].A);
		mx_create_set_scalar_field<m_matrix_r>(mx_atom_type, i, "rn_e", atom_type[i].rn_e);
		mx_create_set_scalar_field<m_matrix_r>(mx_atom_type, i, "rn_c", atom_type[i].rn_c);
		mx_create_set_scalar_field<m_matrix_r>(mx_atom_type, i, "ra_e", atom_type[i].ra_e);
		mx_create_set_scalar_field<m_matrix_r>(mx_atom_type, i, "ra_c", atom_type[i].ra_c);
		mx_create_set_scalar_field<m_matrix_r>(mx_atom_type, i, "R_min", atom_type[i].R_min);
		mx_create_set_scalar_field<m_matrix_r>(mx_atom_type, i, "R_max", atom_type[i].R_max);

		/*************************fg***************************/
		mxfield_CoefPar = mxCreateStructArray(2, dims_CoefPar, number_of_fields_CoefPar, field_names_CoefPar);
		mxSetField(mx_atom_type, i, "feg", mxfield_CoefPar);
		mx_create_set_matrix_field<m_matrix_r>(mxfield_CoefPar, "cl", 1, atom_type[i].feg.size(), atom_type[i].feg.cl.data());
		mx_create_set_matrix_field<m_matrix_r>(mxfield_CoefPar, "cnl", 1, atom_type[i].feg.size(), atom_type[i].feg.cnl.data());

		/*************************fx***************************/
		mxfield_CoefPar = mxCreateStructArray(2, dims_CoefPar, number_of_fields_CoefPar, field_names_CoefPar);
		mxSetField(mx_atom_type, i, "fxg", mxfield_CoefPar);
		mx_create_set_matrix_field<m_matrix_r>(mxfield_CoefPar, "cl", 1, atom_type[i].fxg.size(), atom_type[i].fxg.cl.data());
		mx_create_set_matrix_field<m_matrix_r>(mxfield_CoefPar, "cnl", 1, atom_type[i].fxg.size(), atom_type[i].fxg.cnl.data());

		/*************************Pr***************************/
		mxfield_CoefPar = mxCreateStructArray(2, dims_CoefPar, number_of_fields_CoefPar, field_names_CoefPar);
		mxSetField(mx_atom_type, i, "Pr", mxfield_CoefPar);
		mx_create_set_matrix_field<m_matrix_r>(mxfield_CoefPar, "cl", 1, atom_type[i].Pr.size(), atom_type[i].Pr.cl.data());
		mx_create_set_matrix_field<m_matrix_r>(mxfield_CoefPar, "cnl", 1, atom_type[i].Pr.size(), atom_type[i].Pr.cnl.data());

		/*************************Vr***************************/
		mxfield_CoefPar = mxCreateStructArray(2, dims_CoefPar, number_of_fields_CoefPar, field_names_CoefPar);
		mxSetField(mx_atom_type, i, "Vr", mxfield_CoefPar);
		mx_create_set_matrix_field<m_matrix_r>(mxfield_CoefPar, "cl", 1, atom_type[i].Vr.size(), atom_type[i].Vr.cl.data());
		mx_create_set_matrix_field<m_matrix_r>(mxfield_CoefPar, "cnl", 1, atom_type[i].Vr.size(), atom_type[i].Vr.cnl.data());

		/*************************VR***************************/
		mxfield_CoefPar = mxCreateStructArray(2, dims_CoefPar, number_of_fields_CoefPar, field_names_CoefPar);
		mxSetField(mx_atom_type, i, "VR", mxfield_CoefPar);
		mx_create_set_matrix_field<m_matrix_r>(mxfield_CoefPar, "cl", 1, atom_type[i].VR.size(), atom_type[i].VR.cl.data());
		mx_create_set_matrix_field<m_matrix_r>(mxfield_CoefPar, "cnl", 1, atom_type[i].VR.size(), atom_type[i].VR.cnl.data());

		/*************************ciVR***************************/
		mx_create_set_matrix_field<m_matrix_r>(mx_atom_type, i, "R2", 1, atom_type[i].R2.size(), atom_type[i].R2.data());
		mxfield_ciVR = mxCreateStructArray(2, dims_ciVR, number_of_fields_ciVR, field_names_ciVR);
		mxSetField(mx_atom_type, i, "ciVR", mxfield_ciVR);
		mx_create_set_matrix_field<m_matrix_r>(mxfield_ciVR, "c0", 1, atom_type[i].ciVR.size(), atom_type[i].ciVR.c0.data());
		mx_create_set_matrix_field<m_matrix_r>(mxfield_ciVR, "c1", 1, atom_type[i].ciVR.size(), atom_type[i].ciVR.c1.data());
		mx_create_set_matrix_field<m_matrix_r>(mxfield_ciVR, "c2", 1, atom_type[i].ciVR.size(), atom_type[i].ciVR.c2.data());
		mx_create_set_matrix_field<m_matrix_r>(mxfield_ciVR, "c3", 1, atom_type[i].ciVR.size(), atom_type[i].ciVR.c3.data());
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	auto potential_type = mx_get_scalar<multem::ePotential_Type>(prhs[0]);

	multem::Atomic_Data atomic_data;
	atomic_data.Load_Data(potential_type);
	host_vector<multem::Atom_Type<double, multem::e_Host>> atom_type(multem::c_nAtomsTypes);

	for(auto i=0; i<atom_type.size(); i++)
	{
		atomic_data.To_atom_type_CPU(i+1, multem::c_Vrl, multem::c_nR, 0.0, atom_type[i]);
	}

	set_output_data(atom_type, plhs[0]);
}