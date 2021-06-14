/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
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

#include "types.cuh"
#include "type_traits_gen.cuh"
#include "cgpu_stream.cuh"
#include "cgpu_fft.cuh"
#include "particles.cuh"
#include "slicing.hpp"
#include "in_classes.cuh"
#include "output_multem.hpp"
#include "cpu_fcns.hpp"
#include "gpu_fcns.cuh"
#include "cgpu_fcns.cuh"
#include "wave_function.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::pMLD;
using mt::pMx_c;

template <class TIn_Multislice>
void read_in_multem(const mxArray *mex_in_multem, TIn_Multislice &in_multem, dt_bool full = true)
{
	using T_r = mt::Value_type<TIn_Multislice>;

	/************************ simulation type **************************/
	in_multem.simulation_type = mex_get_num_from_field<mt::eTEM_Sim_Typ>(mex_in_multem, "simulation_type");
	in_multem.simulation_type = (in_multem.is_EWRS())?mt::eTEMST_EWRS:mt::eTEMST_EWFS;

	in_multem.interaction_model = mex_get_num_from_field<mt::eElec_Spec_Int_Model>(mex_in_multem, "interaction_model");
	in_multem.pot_parm_typ = mex_get_num_from_field<mt::ePot_Parm_Typ>(mex_in_multem, "pot_parm_typ");

	/************** Electron-Phonon_Par interaction model **************/
	mex_read_phonon_par(mex_in_multem, in_multem.phonon_par);
	in_multem.phonon_par.coh_contrib = true;
	in_multem.phonon_par.single_conf = true;

	/**************************** Specimen *****************************/
	auto bs_x = mex_get_num_from_field<T_r>(mex_in_multem, "spec_bs_x");
	auto bs_y = mex_get_num_from_field<T_r>(mex_in_multem, "spec_bs_y");
	auto bs_z = mex_get_num_from_field<T_r>(mex_in_multem, "spec_bs_z");
	auto sli_thk = mex_get_num_from_field<T_r>(mex_in_multem, "spec_dz");
	dt_bool pbc_xy = true;

	if (full)
	{
		/************************* atomic positions ************************/
		mex_read_atoms<T_r>(mex_in_multem, bs_x, bs_y, bs_z, sli_thk, in_multem.atoms);
	}

	/************************ Specimen rotation ************************/
	mex_read_rot_par<T_r>(mex_in_multem, in_multem.rot_par);

	/************************ Specimen thickness ***********************/
	in_multem.thick_type = mex_get_num_from_field<mt::eThick_Typ>(mex_in_multem, "thick_type");
	if (!in_multem.is_whole_spec() && full)
	{
		mex_read_spec_thick<T_r>(mex_in_multem, in_multem.thick);
	}

	/************************ Potential slicing ************************/
	in_multem.potential_slicing = mex_get_num_from_field<mt::ePot_Sli_Typ>(mex_in_multem, "potential_slicing");

	/************************** xy sampling ****************************/
	auto nx = mex_get_num_from_field<dt_int32>(mex_in_multem, "nx");
	auto ny = mex_get_num_from_field<dt_int32>(mex_in_multem, "ny");
	dt_bool bwl = mex_get_num_from_field<dt_bool>(mex_in_multem, "bwl");

	in_multem.grid_2d.set_in_data(nx, ny, bs_x, bs_y, sli_thk, bwl, pbc_xy);

	/************************** Incident wave **************************/
	auto iw_type = mex_get_num_from_field<mt::eIncident_Wave_Typ>(mex_in_multem, "iw_type");
	in_multem.set_incident_wave_type(iw_type);

	if (in_multem.is_user_define_wave() && full)
	{
		mex_read_user_define_wave<T_r>(mex_in_multem, in_multem.iw_psi);
	}

	/********************* read beam positions *************************/
	mex_read_beam_pos<T_r>(mex_in_multem, in_multem.beam_pos_i);

	/********************* Microscope parameter ***********************/
	in_multem.E_0 = mex_get_num_from_field<T_r>(mex_in_multem, "E_0");
	in_multem.theta = mex_get_num_from_field<T_r>(mex_in_multem, "theta")*mt::c_deg_2_rad;
	in_multem.phi = mex_get_num_from_field<T_r>(mex_in_multem, "phi")*mt::c_deg_2_rad;

	/********************* Illumination model *************************/
	in_multem.illumination_model = mt::eIM_Coherent;
	in_multem.temporal_spatial_incoh = mt::eTSI_Temporal_Spatial;

	/*********************** Condenser lens ***************************/
	mex_read_cond_lens<T_r>(mex_in_multem, in_multem.cond_lens);
	in_multem.cond_lens.set_in_data(in_multem.E_0, in_multem.grid_2d);

	/*********************** Objective lens ***************************/
	mex_read_obj_lens<T_r>(mex_in_multem, in_multem.obj_lens);
	in_multem.obj_lens.set_in_data(in_multem.E_0, in_multem.grid_2d);

	/******************** select output region ************************/
	mex_read_output_area(mex_in_multem, in_multem.output_area);

	/********************* validate parameters *************************/
	in_multem.validate_parameters();
}

template <class TOutput_Multem>
void set_struct_wave_function(TOutput_Multem &output_multem, mxArray*& mex_output_multem)
{
	const char *field_names_output_multem[] = { "dx", "dy", "x", "y", "thick", "data" };
	dt_int32 number_of_fields_output_multem = 6;
	mwSize dims_output_multem[2] = { 1, 1 };

	mxArray *mex_field_psi;
	const char *field_names_psi[] = { "psi_coh" };
	dt_int32 number_of_fields_psi = 1;
	mwSize dims_psi[2] = { 1, output_multem.thick.size() };

	mex_output_multem = mxCreateStructArray(2, dims_output_multem, number_of_fields_output_multem, field_names_output_multem);

	mex_create_set_num_field<pMLD>(mex_output_multem, 0, "dx", output_multem.dx);
	mex_create_set_num_field<pMLD>(mex_output_multem, 0, "dy", output_multem.dy);
	mex_create_set_pVctr_field<pMLD>(mex_output_multem, "x", 1, output_multem.x.size(), output_multem.x);
	mex_create_set_pVctr_field<pMLD>(mex_output_multem, "y", 1, output_multem.y.size(), output_multem.y);
	mex_create_set_pVctr_field<pMLD>(mex_output_multem, "thick", 1, output_multem.thick.size(), output_multem.thick);

	mex_field_psi = mxCreateStructArray(2, dims_psi, number_of_fields_psi, field_names_psi);
	mxSetField(mex_output_multem, 0, "data", mex_field_psi);

	for(auto ithk = 0; ithk < output_multem.thick.size(); ithk++)
	{
		mex_create_set_pVctr_field<pMx_c>(mex_field_psi, ithk, "psi_coh", output_multem.ny, output_multem.nx, output_multem.psi_coh[ithk]);
	}
}

template <class T, mt::eDev Dev>
void run_wave_function(mt::System_Config &system_config, const mxArray *mex_in_multem, mxArray*& mex_output_multem)
{
	mt::In_Multem<T> in_multem;
	read_in_multem(mex_in_multem, in_multem);
	in_multem.system_config = system_config;

	mt::Stream<Dev> stream(system_config.n_stream);
	mt::FFT<T, Dev> fft_2d;
	fft_2d.create_plan_2d(in_multem.grid_2d.ny, in_multem.grid_2d.nx, system_config.n_stream);

	mt::Wave_Function<T, Dev> wave_function;
	wave_function.set_in_data(&in_multem, &stream, &fft_2d);

	mt::Output_Multem<T> output_multem;
	output_multem.set_in_data(&in_multem);

	wave_function.move_atoms(in_multem.phonon_par.nconf);
	wave_function.set_incident_wave(wave_function.psi_z);
	wave_function.psi(1.0, wave_function.psi_z, output_multem);

	stream.synchronize();

	output_multem.gather();
	output_multem.clean_temporal();

	fft_2d.cleanup();

	set_struct_wave_function(output_multem, mex_output_multem);
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	auto system_config = mex_read_system_config(prhs[0]);
	system_config.set_gpu();

	dt_int32 idx_0 = (system_config.active)?1:0;

	if (system_config.is_float32_cpu())
	{
		// run_wave_function<dt_float32, mt::edev_cpu>(system_config, prhs[idx_0], plhs[0]);
	}
	else if (system_config.is_float64_cpu())
	{
		// run_wave_function<dt_float64, mt::edev_cpu>(system_config, prhs[idx_0], plhs[0]);
	}
	else if (system_config.is_float32_gpu())
	{
		run_wave_function<dt_float32, mt::edev_gpu>(system_config, prhs[idx_0], plhs[0]);
	}
	else if (system_config.is_float64_gpu())
	{
		// run_wave_function<dt_float64, mt::edev_gpu>(system_config, prhs[idx_0], plhs[0]);
	}
}