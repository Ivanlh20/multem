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

#include <algorithm>

#include "types.cuh"
#include "type_traits_gen.h"
#include "fcns_gpu.h"
#include "in_classes.cuh"
#include "output_multem.hpp"
#include "particles.cuh"
#include "tem_simulation.cuh"
#include "timing.cuh"

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

	/***************************************************************************************/
	multem_in_parm.elec_spec_interact_mod = mex_get_num_from_field<mt::eElec_Spec_Interact_Mod>(mex_multem_in_parm, "elec_spec_interact_mod");
	multem_in_parm.atomic_pot_parm_typ = mex_get_num_from_field<mt::eAtomic_Pot_Parm_Typ>(mex_multem_in_parm, "atomic_pot_parm_typ");

	/***************************************************************************************/
	multem_in_parm.operation_mode = mex_get_num_from_field<mt::eOperation_Mode>(mex_multem_in_parm, "operation_mode");
	multem_in_parm.reverse_multislice = mex_get_bool_from_field(mex_multem_in_parm, "reverse_multislice");

	/************** Electron-Atomic_Vib interaction model **************/
	mex_read_atomic_vib(mex_multem_in_parm, multem_in_parm.atomic_vib);

	/**************************** Specimen *****************************/
	auto bs_x = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_bs_x");
	auto bs_y = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_bs_y");
	auto bs_z = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_bs_z");
	auto sli_thick = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_dz");
	dt_bool pbc_xy = true;

	if (multem_in_parm.is_specimen_required())
	{
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
	}

	/************************** xy sampling ****************************/
	auto nx = mex_get_num_from_field<dt_int32>(mex_multem_in_parm, "nx");
	auto ny = mex_get_num_from_field<dt_int32>(mex_multem_in_parm, "ny");
	dt_bool bwl = mex_get_bool_from_field(mex_multem_in_parm, "bwl");

	multem_in_parm.grid_2d.set_in_data(nx, ny, bs_x, bs_y, sli_thick, bwl, pbc_xy);

	/************************ Incident wave ****************************/
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
	multem_in_parm.illum_mod = mex_get_num_from_field<mt::eIllum_Mod>(mex_multem_in_parm, "illum_mod");
	multem_in_parm.illum_inc = mex_get_num_from_field<mt::eIllum_Inc>(mex_multem_in_parm, "illum_inc");

	/*********************** Condenser lens ***************************/
	mex_read_cond_lens<T_r>(mex_multem_in_parm, multem_in_parm.cond_lens);
	multem_in_parm.cond_lens.set_in_data(multem_in_parm.E_0, multem_in_parm.grid_2d);

	/*********************** Objective lens ***************************/
	mex_read_obj_lens<T_r>(mex_multem_in_parm, multem_in_parm.obj_lens);
	multem_in_parm.obj_lens.set_in_data(multem_in_parm.E_0, multem_in_parm.grid_2d);

	/************************** ISTEM/STEM ***************************/
	if (multem_in_parm.is_scanning())
	{
		mex_read_scan<T_r>(mex_multem_in_parm, multem_in_parm.scanning);
	}

	if (multem_in_parm.is_STEM())
	{
		mex_read_detector<T_r>(mex_multem_in_parm, multem_in_parm.E_0, multem_in_parm.detector);
	}
	else if (multem_in_parm.is_PED())
	{
		multem_in_parm.theta = mex_get_num_from_field<T_r>(mex_multem_in_parm, "ped_theta")*mt::c_deg_2_rad;
		multem_in_parm.nrot = mex_get_num_from_field<T_r>(mex_multem_in_parm, "ped_nrot");
	}
	else if (multem_in_parm.is_HCTEM())
	{
		multem_in_parm.theta = mex_get_num_from_field<T_r>(mex_multem_in_parm, "hci_theta")*mt::c_deg_2_rad;
		multem_in_parm.nrot = mex_get_num_from_field<T_r>(mex_multem_in_parm, "hci_nrot");
	}
	else if (multem_in_parm.is_STEM_ISTEM_EELS())
	{
		auto Z = mex_get_num_from_field<dt_int32>(mex_multem_in_parm, "eels_Z");
		auto E_loss = mex_get_num_from_field<T_r>(mex_multem_in_parm, "eels_E_loss")*mt::c_eV_2_keV;
		auto coll_angle = mex_get_num_from_field<dt_float64>(mex_multem_in_parm, "eels_coll_angle")*mt::c_mrad_2_rad<dt_float64>;
		auto m_sel = mex_get_num_from_field<dt_int32>(mex_multem_in_parm, "eels_m_sel");
		auto chan_type = mex_get_num_from_field<mt::eChan_Typ>(mex_multem_in_parm, "eels_chan_type");

 		mt::eSpace space = mt::esp_fourier;
		multem_in_parm.eels_fr.set_in_data(space, multem_in_parm.E_0, E_loss, m_sel, coll_angle, chan_type, Z);
	}
	else if (multem_in_parm.is_EFTEM())
	{
		auto Z = mex_get_num_from_field<dt_int32>(mex_multem_in_parm, "eftem_Z");
		auto E_loss = mex_get_num_from_field<T_r>(mex_multem_in_parm, "eftem_E_loss")*mt::c_eV_2_keV;
		auto coll_angle = mex_get_num_from_field<dt_float64>(mex_multem_in_parm, "eftem_coll_angle")*mt::c_mrad_2_rad<dt_float64>;
		auto m_sel = mex_get_num_from_field<dt_int32>(mex_multem_in_parm, "eftem_m_sel");
		auto chan_type = mex_get_num_from_field<mt::eChan_Typ>(mex_multem_in_parm, "eftem_chan_type");

		mt::eSpace space = (multem_in_parm.is_EFTEMRS())?mt::esp_real:mt::esp_fourier;
		multem_in_parm.eels_fr.set_in_data(space, multem_in_parm.E_0, E_loss, m_sel, coll_angle, chan_type, Z);
	}

	/******************** select output region ************************/
	mex_read_output_area(mex_multem_in_parm, multem_in_parm.output_area);

	/********************* validate parameters *************************/
	multem_in_parm.set_dep_var();
}

template <class TOutput_Multem>
void set_struct_multem(TOutput_Multem &output_multem, mxArray*& mex_output_multem)
{
	const char *field_names_output_multem[] = {"dx", "dy", "x", "y", "thick", "data"};
	dt_int32 number_of_fields_output_multem = 6;
	mwSize dims_output_multem[2] = {1, 1};

	mex_output_multem = mxCreateStructArray(2, dims_output_multem, number_of_fields_output_multem, field_names_output_multem);

	mex_create_set_num_field<pMLD>(mex_output_multem, 0, "dx", output_multem.dx);
	mex_create_set_num_field<pMLD>(mex_output_multem, 0, "dy", output_multem.dy);
	mex_create_set_pVctr_field<pMLD>(mex_output_multem, "x", 1, output_multem.x.size(), output_multem.x);
	mex_create_set_pVctr_field<pMLD>(mex_output_multem, "y", 1, output_multem.y.size(), output_multem.y);
	mex_create_set_pVctr_field<pMLD>(mex_output_multem, "thick", 1, output_multem.thick.size(), output_multem.thick);

	if (output_multem.is_STEM() || output_multem.is_STEM_ISTEM_EELS())
	{
		mxArray* mex_field_data;
		const char *field_names_data_full[] = { "image_tot", "image_coh" };
		const char *field_names_data_partial[] = { "image_tot" };
		const char **field_names_data = (output_multem.atomic_vib.coh_contrib)?field_names_data_full:field_names_data_partial;
		dt_int32 number_of_fields_data = (output_multem.atomic_vib.coh_contrib)?2:1;
		mwSize dims_data[2] = { 1, output_multem.thick.size() };

		mex_field_data = mxCreateStructArray(2, dims_data, number_of_fields_data, field_names_data);
		mxSetField(mex_output_multem, 0, "data", mex_field_data);

		mxArray* mex_field_detector_tot;
		mxArray* mex_field_detector_coh;
		const char *field_names_detector[] = { "image" };
		dt_int32 number_of_fields_detector = 1;
		// mwSize dims_detector[2] = {1, output_multem.ndetector};
		mwSize dims_detector[2];
		dims_detector[0] = 1;
		dims_detector[1] = output_multem.ndetector;

		dt_int32 nx = (output_multem.scanning.is_scan_pat_line())?1:output_multem.nx;
		dt_int32 ny = output_multem.ny;
		for(auto ithk = 0; ithk < output_multem.thick.size(); ithk++)
		{
			mex_field_detector_tot = mxCreateStructArray(2, dims_detector, number_of_fields_detector, field_names_detector);
			mxSetField(mex_field_data, ithk, "image_tot", mex_field_detector_tot);
			if (output_multem.atomic_vib.coh_contrib)
			{
				mex_field_detector_coh = mxCreateStructArray(2, dims_detector, number_of_fields_detector, field_names_detector);
				mxSetField(mex_field_data, ithk, "image_coh", mex_field_detector_coh);
			}

			for(auto iDet = 0; iDet < output_multem.ndetector; iDet++)
			{
				mex_create_set_pVctr_field<pMLD>(mex_field_detector_tot, iDet, "image", ny, nx, output_multem.image_tot(ithk, iDet));
				if (output_multem.atomic_vib.coh_contrib)
				{
					mex_create_set_pVctr_field<pMLD>(mex_field_detector_coh, iDet, "image", ny, nx, output_multem.image_coh(ithk, iDet));
				}
			}
		}
	}
	else if (output_multem.is_EWFS_EWRS())
	{
		mxArray* mex_field_data;
		const char *field_names_data_full[] = { "m2psi_tot", "psi_coh" };
		const char *field_names_data_partial[] = { "psi_coh" };
		const char **field_names_data = (!output_multem.is_EWFS_EWRS_SC())?field_names_data_full:field_names_data_partial;
		dt_int32 number_of_fields_data = (!output_multem.is_EWFS_EWRS_SC())?2:1;
		mwSize dims_data[2] = { 1, output_multem.thick.size() };

		mex_field_data = mxCreateStructArray(2, dims_data, number_of_fields_data, field_names_data);
		mxSetField(mex_output_multem, 0, "data", mex_field_data);

		for(auto ithk = 0; ithk < output_multem.thick.size(); ithk++)
		{
			if (!output_multem.is_EWFS_EWRS_SC())
			{
				mex_create_set_pVctr_field<pMLD>(mex_field_data, ithk, "m2psi_tot", output_multem.ny, output_multem.nx, output_multem.m2psi_tot[ithk]);
			}
			mex_create_set_pVctr_field<pMx_c>(mex_field_data, ithk, "psi_coh", output_multem.ny, output_multem.nx, output_multem.psi_coh[ithk]);
		}
	}
	else
	{
		mxArray* mex_field_data;
		const char *field_names_data_full[] = { "m2psi_tot", "m2psi_coh" };
		const char *field_names_data_partial[] = { "m2psi_tot" };
		const char **field_names_data = (output_multem.atomic_vib.coh_contrib)?field_names_data_full:field_names_data_partial;
		dt_int32 number_of_fields_data = (output_multem.atomic_vib.coh_contrib)?2:1;
		mwSize dims_data[2] = { 1, output_multem.nthk() };

		mex_field_data = mxCreateStructArray(2, dims_data, number_of_fields_data, field_names_data);
		mxSetField(mex_output_multem, 0, "data", mex_field_data);

		for(auto ithk = 0; ithk < output_multem.nthk(); ithk++)
		{
			mex_create_set_pVctr_field<pMLD>(mex_field_data, ithk, "m2psi_tot", output_multem.ny, output_multem.nx, output_multem.m2psi_tot[ithk]);
			if (output_multem.atomic_vib.coh_contrib)
			{
				mex_create_set_pVctr_field<pMLD>(mex_field_data, ithk, "m2psi_coh", output_multem.ny, output_multem.nx, output_multem.m2psi_coh[ithk]);
			}
		}
	}
}

template <class T, mt::eDev Dev>
void run_multem(mt::System_Config &system_config, const mxArray* mex_multem_in_parm, mxArray*& mex_output_multem)
{
	mt::Multem_In_Parm<T> multem_in_parm;
	read_multem_in_parm(mex_multem_in_parm, multem_in_parm);
	multem_in_parm.system_config = system_config;

	// multem has to be created first due to memory allocations
	mt::Multem<T, Dev> multem(&multem_in_parm);
	mt::Output_Multem<T> output_multem(&multem_in_parm, true);

	multem(output_multem);

	set_struct_multem(output_multem, mex_output_multem);
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	auto system_config = mex_read_system_config(prhs[0]);
	dt_int32 idx_0 = (system_config.active)?1:0;

	mt::d_grid_blk(1024, 1024, 0, 1);

	if (system_config.is_float32_cpu())
	{
		// mexPrintf("cpu - dt_float32 precision calculation\n");
		run_multem<dt_float32, mt::edev_cpu>(system_config, prhs[idx_0], plhs[0]);
	}
	else if (system_config.is_float64_cpu())
	{
		// mexPrintf("cpu - dt_float64 precision calculation\n");
		run_multem<dt_float64, mt::edev_cpu>(system_config, prhs[idx_0], plhs[0]);
	}
	if (system_config.is_float32_gpu())
	{
		mexPrintf("gpu - dt_float32 precision calculation\n");
		run_multem<dt_float32, mt::edev_gpu>(system_config, prhs[idx_0], plhs[0]);
	}
	else if (system_config.is_float64_gpu())
	{
		// mexPrintf("gpu - dt_float64 precision calculation\n");
		run_multem<dt_float64, mt::edev_gpu>(system_config, prhs[idx_0], plhs[0]);
	}
	// // cudaDeviceReset();
}