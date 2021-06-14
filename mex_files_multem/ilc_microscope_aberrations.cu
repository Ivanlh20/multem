/**
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

#include "types.cuh"
#include "type_traits_gen.cuh"
#include "cgpu_stream.cuh"
#include "cgpu_fft.cuh"
#include "particles.cuh"
#include "in_classes.cuh"
#include "output_multem.hpp"

#include "microscope_effects.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::pMLD;
using mt::pMx_c;

template <class TIn_Multislice>
void read_in_multem(const mxArray *mex_in_multem, TIn_Multislice &in_multem, dt_bool full = true)
{
	using T_r = mt::Value_type<TIn_Multislice>;

	in_multem.simulation_type = mt::eTEMST_HRTEM;

	/**************************** Specimen *****************************/
	auto bs_x = mex_get_num_from_field<T_r>(mex_in_multem, "spec_bs_x");
	auto bs_y = mex_get_num_from_field<T_r>(mex_in_multem, "spec_bs_y");
	T_r bs_z = 0;
	T_r sli_thk = 0.25;
	dt_bool pbc_xy = true;

	/************************** xy sampling ****************************/
	auto nx = mex_get_num_from_field<dt_int32>(mex_in_multem, "nx");
	auto ny = mex_get_num_from_field<dt_int32>(mex_in_multem, "ny");
	dt_bool bwl = false;

	in_multem.grid_2d.set_in_data(nx, ny, bs_x, bs_y, sli_thk, bwl, pbc_xy);

	/************************ Incident wave ****************************/
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
	in_multem.illumination_model = mex_get_num_from_field<mt::eIllumination_Model>(mex_in_multem, "illumination_model");
	in_multem.temporal_spatial_incoh = mex_get_num_from_field<mt::eTemporal_Spatial_Incoh>(mex_in_multem, "temporal_spatial_incoh");

	/*********************** Condenser lens ***************************/
	mex_read_cond_lens<T_r>(mex_in_multem, in_multem.cond_lens);
	in_multem.cond_lens.set_in_data(in_multem.E_0, in_multem.grid_2d);

	/*********************** Objective lens ***************************/
	mex_read_obj_lens<T_r>(mex_in_multem, in_multem.obj_lens);
	in_multem.obj_lens.set_in_data(in_multem.E_0, in_multem.grid_2d);

	/********************* zero defocus reference ********************/
	in_multem.obj_lens.zero_defocus_type = mt::eZDT_Last;
	in_multem.obj_lens.zero_defocus_plane = 0.0;

	/******************** select output region ************************/
	mex_read_output_area(mex_in_multem, in_multem.output_area);

	/********************* validate parameters *************************/
	in_multem.validate_parameters();
}

template <class TOutput_Multem>
void set_struct_microscope_aberrations(TOutput_Multem &output_multem, mxArray*& mex_output_multem)
{
	const char *field_names_output_multem[] = { "dx", "dy", "x", "y", "thick", "m2psi" };
	dt_int32 number_of_fields_output_multem = 6;
	mwSize dims_output_multem[2] = { 1, 1 };

	mex_output_multem = mxCreateStructArray(2, dims_output_multem, number_of_fields_output_multem, field_names_output_multem);

	mex_create_set_num_field<pMLD>(mex_output_multem, 0, "dx", output_multem.dx);
	mex_create_set_num_field<pMLD>(mex_output_multem, 0, "dy", output_multem.dy);
	mex_create_set_pVctr_field<pMLD>(mex_output_multem, "x", 1, output_multem.x.size(), output_multem.x);
	mex_create_set_pVctr_field<pMLD>(mex_output_multem, "y", 1, output_multem.y.size(), output_multem.y);
	mex_create_set_pVctr_field<pMLD>(mex_output_multem, "thick", 1, output_multem.thick.size(), output_multem.thick);
	mex_create_set_pVctr_field<pMx_c>(mex_output_multem, "m2psi", output_multem.ny, output_multem.nx, output_multem.m2psi_tot[0]);
}

template <class T, mt::eDev Dev>
void run_microscope_aberrations(mt::System_Config &system_config, const mxArray *mex_in_multem, mxArray*& mex_output_multem)
{
	mt::In_Multem<T> in_multem;
	read_in_multem(mex_in_multem, in_multem);
	in_multem.system_config = system_config;

	mt::Stream<Dev> stream(system_config.n_stream);
	mt::FFT<T, Dev> fft_2d;
	fft_2d.create_plan_2d(in_multem.grid_2d.ny, in_multem.grid_2d.nx, system_config.n_stream);

	mt::Microscope_Effects<T, Dev> microscope_effects;
	microscope_effects.set_in_data(&in_multem, &stream, &fft_2d);

	mt::Output_Multem<T> output_multem;
	output_multem.set_in_data(&in_multem);

	microscope_effects(output_multem);

	stream.synchronize();

	output_multem.gather();
	output_multem.clean_temporal();
	fft_2d.cleanup();

	set_struct_microscope_aberrations(output_multem, mex_output_multem);
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	auto system_config = mex_read_system_config(prhs[0]);
	system_config.set_gpu();

	dt_int32 idx_0 = (system_config.active)?1:0;

	if (system_config.is_float32_cpu())
	{
		run_microscope_aberrations<dt_float32, mt::edev_cpu>(system_config, prhs[idx_0], plhs[0]);
	}
	else if (system_config.is_float64_cpu())
	{
		run_microscope_aberrations<dt_float64, mt::edev_cpu>(system_config, prhs[idx_0], plhs[0]);
	}
	else if (system_config.is_float32_gpu())
	{
		run_microscope_aberrations<dt_float32, mt::edev_gpu>(system_config, prhs[idx_0], plhs[0]);
	}
	else if (system_config.is_float64_gpu())
	{
		run_microscope_aberrations<dt_float64, mt::edev_gpu>(system_config, prhs[idx_0], plhs[0]);
	}
}