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

#include "types.cuh"
#include "type_traits_gen.h"
#include "cgpu_stream.cuh"
#include "cgpu_fft.cuh"
#include "particles.cuh"
#include "spec_slic.hpp"
#include "in_classes.cuh"
#include "output_multem.hpp"
#include "fcns_cpu.h"
#include "fcns_gpu.h"
#include "fcns_gpu.h"
#include "wave_function.cuh"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;
using mt::pMx_c;

template <class TIn_Multislice>
void read_multem_in_parm(const mxArray* mex_multem_in_parm, TIn_Multislice &multem_in_parm, dt_bool full = true)
{
	using T_r = mt::Value_type<TIn_Multislice>;

	/************************ simulation type **************************/
	multem_in_parm.em_sim_typ = mex_get_num_from_field<mt::eEM_Sim_Typ>(mex_multem_in_parm, "em_sim_typ");
	multem_in_parm.em_sim_typ = (multem_in_parm.is_EWRS())?mt::eemst_ewrs:mt::eemst_ewfs;

	multem_in_parm.elec_spec_interact_mod = mex_get_num_from_field<mt::eElec_Spec_Interact_Mod>(mex_multem_in_parm, "elec_spec_interact_mod");
	multem_in_parm.atomic_pot_parm_typ = mex_get_num_from_field<mt::eAtomic_Pot_Parm_Typ>(mex_multem_in_parm, "atomic_pot_parm_typ");

	/************** Electron-Atomic_Vib interaction model **************/
	mex_read_atomic_vib(mex_multem_in_parm, multem_in_parm.atomic_vib);
	multem_in_parm.atomic_vib.coh_contrib = true;
	multem_in_parm.atomic_vib.sgl_conf = true;

	/**************************** Specimen *****************************/
	auto bs_x = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_bs_x");
	auto bs_y = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_bs_y");
	auto bs_z = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_bs_z");
	auto sli_thick = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_dz");
	dt_bool pbc_xy = true;

	if (full)
	{
		/************************* atomic positions ************************/
		mex_read_atoms<T_r>(mex_multem_in_parm, bs_x, bs_y, bs_z, sli_thick, multem_in_parm.atoms);
	}

	/************************ Specimen rotation ************************/
	mex_read_rot_in_parm<T_r>(mex_multem_in_parm, multem_in_parm.rot_in_parm);

	/************************ Specimen thickness ***********************/
	multem_in_parm.thick_type = mex_get_num_from_field<mt::eSim_Thick_Typ>(mex_multem_in_parm, "thick_type");
	if (!multem_in_parm.is_sim_whole_spec() && full)
	{
		mex_read_spec_thick<T_r>(mex_multem_in_parm, multem_in_parm.thick);
	}

	/************************ Potential slicing ************************/
	multem_in_parm.spec_slic_typ = mex_get_num_from_field<mt::eSpec_Slic_Typ>(mex_multem_in_parm, "spec_slic_typ");

	/************************** xy sampling ****************************/
	auto nx = mex_get_num_from_field<dt_int32>(mex_multem_in_parm, "nx");
	auto ny = mex_get_num_from_field<dt_int32>(mex_multem_in_parm, "ny");
	dt_bool bwl = mex_get_bool_from_field(mex_multem_in_parm, "bwl");

	multem_in_parm.grid_2d.set_in_data(nx, ny, bs_x, bs_y, sli_thick, bwl, pbc_xy);

	/************************** Incident wave **************************/
	auto iw_type = mex_get_num_from_field<mt::eIncident_Wave_Typ>(mex_multem_in_parm, "iw_type");
	multem_in_parm.set_incident_wave_type(iw_type);

	if (multem_in_parm.is_user_define_wave() && full)
	{
		mex_read_user_define_wave<T_r>(mex_multem_in_parm, multem_in_parm.iw_psi);
	}

	/********************* read beam positions *************************/
	mex_read_beam_pos<T_r>(mex_multem_in_parm, multem_in_parm.beam_pos_2d);

	/********************* Microscope parameter ***********************/
	multem_in_parm.E_0 = mex_get_num_from_field<T_r>(mex_multem_in_parm, "E_0");
	multem_in_parm.theta = mex_get_num_from_field<T_r>(mex_multem_in_parm, "theta")*mt::c_deg_2_rad;
	multem_in_parm.phi = mex_get_num_from_field<T_r>(mex_multem_in_parm, "phi")*mt::c_deg_2_rad;

	/********************* Illumination model *************************/
	multem_in_parm.illum_mod = mt::eim_coherent;
	multem_in_parm.illum_inc = mt::eii_temporal_spatial;

	/*********************** Condenser lens ***************************/
	mex_read_cond_lens<T_r>(mex_multem_in_parm, multem_in_parm.cond_lens);
	multem_in_parm.cond_lens.set_in_data(multem_in_parm.E_0, multem_in_parm.grid_2d);

	/*********************** Objective lens ***************************/
	mex_read_obj_lens<T_r>(mex_multem_in_parm, multem_in_parm.obj_lens);
	multem_in_parm.obj_lens.set_in_data(multem_in_parm.E_0, multem_in_parm.grid_2d);

	/******************** select output region ************************/
	mex_read_output_area(mex_multem_in_parm, multem_in_parm.output_area);

	/********************* validate parameters *************************/
	multem_in_parm.set_dep_var();
}

template <class TOutput_Multem>
void set_struct_wave_function(TOutput_Multem &output_multem, mxArray*& mex_output_multem)
{
	const char *field_names_output_multem[] = { "dx", "dy", "x", "y", "thick", "data" };
	dt_int32 number_of_fields_output_multem = 6;
	mwSize dims_output_multem[2] = { 1, 1 };

	mxArray* mex_field_psi;
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
void run_wave_function(mt::System_Config &system_config, const mxArray* mex_multem_in_parm, mxArray*& mex_output_multem)
{
	mt::Multem_In_Parm<T> multem_in_parm;
	read_multem_in_parm(mex_multem_in_parm, multem_in_parm);
	multem_in_parm.system_config = system_config;

	mt::Stream<Dev> stream(system_config.n_stream);
	mt::FFT<T, Dev> fft_2d;
	fft_2d.create_plan_2d(multem_in_parm.grid_2d.ny, multem_in_parm.grid_2d.nx, system_config.n_stream);

	mt::Wave_Function<T, Dev> wave_function;
	wave_function.set_in_data(&multem_in_parm, &stream, &fft_2d);

	mt::Output_Multem<T> output_multem;
	output_multem.set_in_data(&multem_in_parm);

	wave_function.move_atoms(multem_in_parm.atomic_vib.nconf);
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