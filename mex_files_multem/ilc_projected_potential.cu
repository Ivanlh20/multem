/*
 * This file is part of Multem.
 * Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
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
#include "particles.cuh"
#include "in_classes.cuh"
#include "output_multem.hpp"

#include "projected_potential.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::pMLD;

template <class TIn_Multislice>
void read_in_multem(const mxArray *mex_in_multem, TIn_Multislice &in_multem, dt_bool full = true)
{
	using T_r = mt::Value_type<TIn_Multislice>;

	/************************ simulation type **************************/
	in_multem.simulation_type = mt::eemst_pprs;
	
	/***************************************************************************************/
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

	/************************ Potential slicing ************************/
	in_multem.potential_slicing = mex_get_num_from_field<mt::ePot_Sli_Typ>(mex_in_multem, "potential_slicing");

	/************************** xy sampling ****************************/
	auto nx = mex_get_num_from_field<dt_int32>(mex_in_multem, "nx");
	auto ny = mex_get_num_from_field<dt_int32>(mex_in_multem, "ny");
	dt_bool bwl = mex_get_num_from_field<dt_bool>(mex_in_multem, "bwl");

	in_multem.grid_2d.set_in_data(nx, ny, bs_x, bs_y, sli_thk, bwl, pbc_xy);

	/************************ simulation type **************************/
	in_multem.islice = mex_get_num_from_field<dt_int32>(mex_in_multem, "islice")-1;

	/******************** select output region ************************/
	mex_read_output_area(mex_in_multem, in_multem.output_area);

	/********************* validate parameters *************************/
	in_multem.validate_parameters();
 }

template <class TOutput_Multem>
void set_struct_projected_potential(TOutput_Multem &output_multem, mxArray*& mex_output_multem)
{
	const char *field_names_output_multem[] = {"dx", "dy", "x", "y", "thick", "V"};
	dt_int32 number_of_fields_output_multem = 6;
	mwSize dims_output_multem[2] = {1, 1};

	mex_output_multem = mxCreateStructArray(2, dims_output_multem, number_of_fields_output_multem, field_names_output_multem);

	mex_create_set_num_field<pMLD>(mex_output_multem, 0, "dx", output_multem.dx);
	mex_create_set_num_field<pMLD>(mex_output_multem, 0, "dy", output_multem.dy);
	mex_create_set_pVctr_field<pMLD>(mex_output_multem, "x", 1, output_multem.x.size(), output_multem.x);
	mex_create_set_pVctr_field<pMLD>(mex_output_multem, "y", 1, output_multem.y.size(), output_multem.y);
	mex_create_set_pVctr_field<pMLD>(mex_output_multem, "thick", 1, output_multem.thick.size(), output_multem.thick);
	mex_create_set_pVctr_field<pMLD>(mex_output_multem, "V", output_multem.ny, output_multem.nx, output_multem.V[0]);
}

template <class T, mt::eDev Dev>
void run_projected_potential(mt::System_Config &system_config, const mxArray *mex_in_multem, mxArray*& mex_output_multem)
{
	mt::In_Multem<T> in_multem;
	read_in_multem(mex_in_multem, in_multem);
	in_multem.system_config = system_config;

	mt::Stream<Dev> stream(system_config.n_stream);
	mt::Projected_Potential<T, Dev> projected_potential;
	projected_potential.set_in_data(&in_multem, &stream);

	mt::Output_Multem<T> output_multem;
	output_multem.set_in_data(&in_multem);

	projected_potential.move_atoms(in_multem.phonon_par.nconf);
	projected_potential(in_multem.islice, output_multem);

	stream.synchronize();

	output_multem.gather();
	output_multem.clean_temporal();

	set_struct_projected_potential(output_multem, mex_output_multem);
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	// the output potential is in V - Angstrom
	auto system_config = mex_read_system_config(prhs[0]);
	system_config.set_gpu();

	dt_int32 idx_0 = (system_config.active)?1:0;

	if (system_config.is_float32_cpu())
	{
		run_projected_potential<dt_float32, mt::edev_cpu>(system_config, prhs[idx_0], plhs[0]);
	}
	else if (system_config.is_float64_cpu())
	{
		run_projected_potential<dt_float64, mt::edev_cpu>(system_config, prhs[idx_0], plhs[0]);
	}
	else if (system_config.is_float32_gpu())
	{
		run_projected_potential<dt_float32, mt::edev_gpu>(system_config, prhs[idx_0], plhs[0]);
	}
	else if (system_config.is_float64_gpu())
	{
		run_projected_potential<dt_float64, mt::edev_gpu>(system_config, prhs[idx_0], plhs[0]);
	}
}