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
#include "type_traits_gen.h"
#include "cgpu_stream.cuh"
#include "particles.cuh"
#include "in_classes.cuh"
#include "output_multem.hpp"

#include "projected_potential.cuh"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;

template <class TIn_Multislice>
void read_multem_in_parm(const mxArray* mex_multem_in_parm, TIn_Multislice &multem_in_parm, dt_bool full = true)
{
	using T_r = mt::Value_type<TIn_Multislice>;

	/************************ simulation type **************************/
	multem_in_parm.em_sim_typ = mt::eemst_pprs;
	
	/***************************************************************************************/
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

	/************************ Potential slicing ************************/
	multem_in_parm.spec_slic_typ = mex_get_num_from_field<mt::eSpec_Slic_Typ>(mex_multem_in_parm, "spec_slic_typ");

	/************************** xy sampling ****************************/
	auto nx = mex_get_num_from_field<dt_int32>(mex_multem_in_parm, "nx");
	auto ny = mex_get_num_from_field<dt_int32>(mex_multem_in_parm, "ny");
	dt_bool bwl = mex_get_bool_from_field(mex_multem_in_parm, "bwl");

	multem_in_parm.grid_2d.set_in_data(nx, ny, bs_x, bs_y, sli_thick, bwl, pbc_xy);

	/************************ simulation type **************************/
	multem_in_parm.islice = mex_get_num_from_field<dt_int32>(mex_multem_in_parm, "islice")-1;

	/******************** select output region ************************/
	mex_read_output_area(mex_multem_in_parm, multem_in_parm.output_area);

	/********************* validate parameters *************************/
	multem_in_parm.set_dep_var();
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
void run_projected_potential(mt::System_Config &system_config, const mxArray* mex_multem_in_parm, mxArray*& mex_output_multem)
{
	mt::Multem_In_Parm<T> multem_in_parm;
	read_multem_in_parm(mex_multem_in_parm, multem_in_parm);
	multem_in_parm.system_config = system_config;

	mt::Stream<Dev> stream(system_config.n_stream);
	mt::Projected_Potential<T, Dev> projected_potential;
	projected_potential.set_in_data(&multem_in_parm, &stream);

	mt::Output_Multem<T> output_multem;
	output_multem.set_in_data(&multem_in_parm);

	projected_potential.move_atoms(multem_in_parm.atomic_vib.nconf);
	projected_potential(multem_in_parm.islice, output_multem);

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