/*
 * This file is part of Multem.
 * Copyright 2022 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * Multem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version of the License, or
 * (at your option) any later version.
 *
 * Multem is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Multem. If not, see <http:// www.gnu.org/licenses/>.
 */

#define MATLAB_BLAS_LAPACK

#include "atomic_data_mt.cuh"
#include "atomic_fcns_mt.cuh"

#include "mex.h"
#include "matlab_mex.h"

template <class T>
void set_output_data(const mt::Vctr_cpu<mt::Atomic_Info_cpu<T>>& atomic_info, mxArray*& mx_atomic_info)
{
	const char *field_names[] = {"Z", "m", "A", "rn", "ra", "eels_maj_edg", "eels_min_edg", "coef"};
	dt_int32 number_of_fields = 8;
	mwSize dims[2] = {1, atomic_info.size()};

	mxArray* mx_field_coef = nullptr;
	const char *field_names_Coef[] = {"charge", "feg", "fxg", "pr", "vr", "vzp"};
	dt_int32 number_of_fields_Coef = 6;
	mwSize dims_Coef[2] = {1, 1};

	mxArray* mx_field_coef_par = nullptr;
	const char *field_names_coef_par[] = {"cl", "cnl"};
	dt_int32 number_of_fields_coef_par = 2;
	mwSize dims_coef_par[2] = {1, 1};

	mx_atomic_info = mxCreateStructArray(2, dims, number_of_fields, field_names);
	for(auto i = 0; i < atomic_info.size(); i++)
	{
		mex_create_set_num_field<dt_float64, T>(mx_atomic_info, i, "Z", atomic_info[i].Z);
		mex_create_set_num_field<dt_float64, T>(mx_atomic_info, i, "m", atomic_info[i].m);
		mex_create_set_num_field<dt_float64, T>(mx_atomic_info, i, "A", atomic_info[i].A);
		mex_create_set_num_field<dt_float64, T>(mx_atomic_info, i, "rn", atomic_info[i].rn);
		mex_create_set_num_field<dt_float64, T>(mx_atomic_info, i, "ra", atomic_info[i].ra);	

		mex_create_set_pVctr_field<dt_float64, dt_float32>(mx_atomic_info, i, "eels_maj_edg", atomic_info[i].eels_maj_edg);
		mex_create_set_pVctr_field<dt_float64, dt_float32>(mx_atomic_info, i, "eels_min_edg", atomic_info[i].eels_min_edg);


		mwSize dims_Coef[2] = {1, atomic_info[i].coef.size()};
		mx_field_coef = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mx_atomic_info, i, "coef", mx_field_coef);

		for(auto j = 0; j < atomic_info[i].coef.size(); j++)
		{
			auto &coef = atomic_info[i].coef[j];

			/***************************************************************************************/
			mex_create_set_num_field<dt_float64, T>(mx_field_coef, j, "charge", coef.charge);

			/***************************************************************************************/
			mx_field_coef_par = mxCreateStructArray(2, dims_coef_par, number_of_fields_coef_par, field_names_coef_par);
			mxSetField(mx_field_coef, j, "feg", mx_field_coef_par);
			mex_create_set_pVctr_field<dt_float64, T>(mx_field_coef_par, "cl", coef.feg.cl);
			mex_create_set_pVctr_field<dt_float64, T>(mx_field_coef_par, "cnl", coef.feg.cnl);

			/***************************************************************************************/
			mx_field_coef_par = mxCreateStructArray(2, dims_coef_par, number_of_fields_coef_par, field_names_coef_par);
			mxSetField(mx_field_coef, j, "fxg", mx_field_coef_par);
			mex_create_set_pVctr_field<dt_float64, T>(mx_field_coef_par, "cl", coef.fxg.cl);
			mex_create_set_pVctr_field<dt_float64, T>(mx_field_coef_par, "cnl", coef.fxg.cnl);

			/***************************************************************************************/
			mx_field_coef_par = mxCreateStructArray(2, dims_coef_par, number_of_fields_coef_par, field_names_coef_par);
			mxSetField(mx_field_coef, j, "pr", mx_field_coef_par);
			mex_create_set_pVctr_field<dt_float64, T>(mx_field_coef_par, "cl", coef.pr.cl);
			mex_create_set_pVctr_field<dt_float64, T>(mx_field_coef_par, "cnl", coef.pr.cnl);

			/***************************************************************************************/
			mx_field_coef_par = mxCreateStructArray(2, dims_coef_par, number_of_fields_coef_par, field_names_coef_par);
			mxSetField(mx_field_coef, j, "vr", mx_field_coef_par);
			mex_create_set_pVctr_field<dt_float64, T>(mx_field_coef_par, "cl", coef.vr.cl);
			mex_create_set_pVctr_field<dt_float64, T>(mx_field_coef_par, "cnl", coef.vr.cnl);

			/***************************************************************************************/
			mx_field_coef_par = mxCreateStructArray(2, dims_coef_par, number_of_fields_coef_par, field_names_coef_par);
			mxSetField(mx_field_coef, j, "vzp", mx_field_coef_par);
			mex_create_set_pVctr_field<dt_float64, T>(mx_field_coef_par, "cl", coef.vzp.cl);
			mex_create_set_pVctr_field<dt_float64, T>(mx_field_coef_par, "cnl", coef.vzp.cnl);

		}
	}
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	auto atomic_pot_parm_typ = mex_get_enum<mt::eAtomic_Pot_Parm_Typ>(prhs[0]);

	mt::Vctr_cpu<mt::Atomic_Info_cpu<dt_float64>> atomic_info(mt::c_n_atom_typ);
	mt::Atomic_Data atomic_data_mt(atomic_pot_parm_typ);

	for(auto ik = 0; ik < atomic_info.size(); ik++)
	{
		auto Z = ik+1;
		atomic_info[ik] = atomic_data_mt(Z, atomic_pot_parm_typ);
	}

	set_output_data(atomic_info, plhs[0]);
}