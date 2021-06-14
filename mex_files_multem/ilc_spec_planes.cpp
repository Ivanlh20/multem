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
#include "particles.cuh"
#include "in_classes.cuh"
#include "slicing.hpp"

#include <mex.h>
#include "matlab_mex.cuh"
#include "matlab_multem_io.cuh"

using mt::pMLD;

template <class TIn_Multislice>
void read_in_multem(const mxArray *mex_in_multem, TIn_Multislice &in_multem)
{
	using T_r = mt::Value_type<TIn_Multislice>;

	in_multem.interaction_model = mex_get_num_from_field<mt::eElec_Spec_Int_Model>(mex_in_multem, "interaction_model");
	in_multem.pot_parm_typ = mt::ePPT_lobato_0_12;

	/************** Electron-Phonon_Par interaction model **************/
	mex_read_phonon_par(mex_in_multem, in_multem.phonon_par);
	in_multem.phonon_par.coh_contrib = false;
	in_multem.phonon_par.single_conf = true;

	/**************************** Specimen *****************************/
	auto bs_x = mex_get_num_from_field<T_r>(mex_in_multem, "spec_bs_x");
	auto bs_y = mex_get_num_from_field<T_r>(mex_in_multem, "spec_bs_y");
	auto bs_z = mex_get_num_from_field<T_r>(mex_in_multem, "spec_bs_z");
	auto sli_thk = mex_get_num_from_field<T_r>(mex_in_multem, "spec_dz");
	dt_bool pbc_xy = false;

	/************************* atomic positions ************************/
	mex_read_atoms<T_r>(mex_in_multem, bs_x, bs_y, bs_z, sli_thk, in_multem.atoms);

	/************************ Specimen rotation ************************/
	mex_read_rot_par<T_r>(mex_in_multem, in_multem.rot_par);

	/************************ Potential slicing ************************/
	in_multem.potential_slicing = mex_get_num_from_field<mt::ePot_Sli_Typ>(mex_in_multem, "potential_slicing");

	/************************** xy sampling ****************************/
	auto nx = 1024;
	auto ny = 1024;
	dt_bool bwl = false;

	in_multem.grid_2d.set_in_data(nx, ny, bs_x, bs_y, sli_thk, bwl, pbc_xy);

	in_multem.validate_parameters();
 }

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{	
	/*************************Input data**************************/
	mt::In_Multem<dt_float64> in_multem;
	read_in_multem(prhs[0], in_multem);

	 /***************************************************************************************/
	mt::Slicing<dt_float64> slicing;
	slicing.set_in_data(&in_multem, &(in_multem.atoms));

	// /************************Output data**************************/
	auto rplanes = mex_create_pVctr<pMLD>(slicing.z_plane.size(), 1, plhs[0]);

	rplanes.assign(slicing.z_plane.begin(), slicing.z_plane.end());
}