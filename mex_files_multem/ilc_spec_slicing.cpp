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
#include "spec.hpp"
#include "in_classes.cuh"

#include <mex.h>
#include "matlab_mex.h"
#include "matlab_multem_io.h"


template <class TIn_Multislice>
void read_multem_in_parm(const mxArray* mex_multem_in_parm, TIn_Multislice& multem_in_parm)
{
	using T_r = mt::Value_type<TIn_Multislice>;

	multem_in_parm.elec_spec_interact_mod = mex_get_num_from_field<mt::eElec_Spec_Interact_Mod>(mex_multem_in_parm, "elec_spec_interact_mod");
	multem_in_parm.atomic_pot_parm_typ = mt::eappt_lobato_0_12;

	/************** Electron-Atomic_Vib interaction model **************/
	mex_read_atomic_vib(mex_multem_in_parm, multem_in_parm.atomic_vib);
	multem_in_parm.atomic_vib.coh_contrib = true;
	multem_in_parm.atomic_vib.sgl_conf = true;

	/**************************** Specimen *****************************/
	auto bs_x = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_bs_x");
	auto bs_y = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_bs_y");
	auto bs_z = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_bs_z");
	auto sli_thick = mex_get_num_from_field<T_r>(mex_multem_in_parm, "spec_dz");
	dt_bool pbc_xy = false;

	/************************* atomic positions ************************/
	mex_read_atoms<T_r>(mex_multem_in_parm, bs_x, bs_y, bs_z, sli_thick, multem_in_parm.atoms);

	/************************ Specimen rotation ************************/
	mex_read_rot_in_parm<T_r>(mex_multem_in_parm, multem_in_parm.rot_in_parm);

	/************************ Potential slicing ************************/
	mex_read_spec_sli_in_parm(mex_multem_in_parm, multem_in_parm.rot_in_parm);

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

	mt::Spec<dt_float64> spec;
	spec.set_in_data(&multem_in_parm);
	spec.move_atoms(multem_in_parm.atomic_vib.nconf);

	// /************************Output data**************************/
	auto atomsM = mex_create_pVctr<pMLD>(spec.atoms.size(), 8, plhs[0]);
	auto sliceM = mex_create_pVctr<pMLD>(spec.slicing.slice.size(), 6, plhs[1]);

	for(auto i = 0; i < atomsM.rows; i++)
	{
		atomsM.real[0*atomsM.rows+i] = spec.atoms.Z[i];
		atomsM.real[1*atomsM.rows+i] = spec.atoms.x[i];
		atomsM.real[2*atomsM.rows+i] = spec.atoms.y[i];
		atomsM.real[3*atomsM.rows+i] = spec.atoms.z[i];
		atomsM.real[4*atomsM.rows+i] = spec.atoms.sigma[i];
		atomsM.real[5*atomsM.rows+i] = spec.atoms.occ[i];
		atomsM.real[6*atomsM.rows+i] = spec.atoms.tag[i];
		atomsM.real[7*atomsM.rows+i] = spec.atoms.charge[i];
	}

	for(auto i = 0; i<sliceM.rows; i++)
	{
		auto slice = spec.slicing.slice[i];

		sliceM.real[0*sliceM.rows+i] = slice.z_0;
		sliceM.real[1*sliceM.rows+i] = slice.z_e;
		sliceM.real[2*sliceM.rows+i] = slice.z_int_0;
		sliceM.real[3*sliceM.rows+i] = slice.z_int_e;
		sliceM.real[4*sliceM.rows+i] = slice.iatom_0+1;
		sliceM.real[5*sliceM.rows+i] = slice.iatom_e+1;
	}
}