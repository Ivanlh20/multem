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
#include "particles.cuh"
#include "multem_in_parm.cuh"
#include "spec_slic.hpp"

#include <mex.h>
#include "matlab_mex.h"
#include "matlab_multem_io.h"

template <class T_r>
void read_multem_in_parm(const mxArray* mex_multem_in_parm, mt::Multem_In_Parm<T_r> &multem_in_parm)
{
	multem_in_parm.elec_spec_interact_mod = mex_get_num_from_field<mt::eElec_Spec_Interact_Mod>(mex_multem_in_parm, "elec_spec_interact_mod");
	multem_in_parm.atomic_pot_parm_typ = mt::eappt_lobato_0_12;

	/************** Electron-Atomic_Vib interaction model **************/
	mex_read_atomic_vib(mex_multem_in_parm, multem_in_parm.atomic_vib);
	multem_in_parm.atomic_vib.coh_contrib = false;
	multem_in_parm.atomic_vib.sgl_conf = true;

	/**************************** Specimen *****************************/
	auto bs_x = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_bs_x");
	auto bs_y = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_bs_y");
	auto bs_z = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_bs_z");
	auto sli_thick = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_dz");
	dt_bool pbc_xy = false;
	dt_bool b_statistic = true;

	/************************* atomic positions ************************/
	mex_read_atoms<T_r>(mex_multem_in_parm, {bs_x, bs_y, bs_z}, pbc_xy, b_statistic, multem_in_parm.atoms);

	/************************ Specimen rotation ************************/
	mex_read_rot_in_parm<T_r>(mex_multem_in_parm, multem_in_parm.rot_in_parm);

	/************************ Potential slicing ************************/
	multem_in_parm.spec_slic_typ = mex_get_enum_from_field<mt::eSpec_Slic_Typ>(mex_multem_in_parm, "spec_slic_typ");

	/************************** xy sampling ****************************/
	auto nx = 1024;
	auto ny = 1024;
	dt_bool bwl = false;

	multem_in_parm.grid_2d.set_in_data(nx, ny, bs_x, bs_y, sli_thick, bwl, pbc_xy);

	multem_in_parm.set_dep_var();
 }

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{	
	/*************************Input data**************************/
	mt::Multem_In_Parm<dt_float64> multem_in_parm;
	read_multem_in_parm(prhs[0], multem_in_parm);

	 /***************************************************************************************/
	mt::Spec_Slic<dt_float64> slicing;
	slicing.set_in_data(&multem_in_parm, &(multem_in_parm.atoms));

	// /************************Output data**************************/
	auto rplanes = mex_create_pVctr<pMLD>(slicing.z_plane.size(), 1, plhs[0]);

	rplanes.assign(slicing.z_plane.begin(), slicing.z_plane.end());
}