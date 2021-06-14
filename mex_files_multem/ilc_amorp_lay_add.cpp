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
#include "particles.cuh"
#include "amorp_build.hpp"

#include <mex.h>
#include "matlab_mex.cuh"

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	auto ratoms = mex_get_pvctr<pMLD>(prhs[0]);
	auto bs_x = mex_get_num<dt_float64>(prhs[1]);
	auto bs_y = mex_get_num<dt_float64>(prhs[2]);
	auto bs_z = mex_get_num<dt_float64>(prhs[3]);
	auto r_min = mex_get_num<dt_float64>(prhs[4]);
	auto Z = mex_get_num<dt_int32>(prhs[5]);
	auto rms_3d = mex_get_num<dt_float64>(prhs[6]);
	auto rho = mex_get_num<dt_float64>(prhs[7]);
	auto lay_pos = (nrhs>8)?mex_get_num<mt::eSpec_Lay_Typ>(prhs[8]):mt::eSLT_Top;
	dt_int32 seed = (nrhs>9)?mex_get_num<dt_int32>(prhs[9]):300183;

	/***************************************************************************************/
	mt::Ptc_Atom<dt_float32> atoms;
	atoms.set_ptc(ratoms.rows, ratoms.cols, ratoms.real, bs_x, bs_y);

	dt_int32 region = (ratoms.rows==0)?0:atoms.region_max+1;

	atoms.spec_lay_info.resize(1);
	if (lay_pos==mt::eSLT_Top)
	{
		atoms.spec_lay_info[0].z_0 = atoms.z_min-bs_z;
		atoms.spec_lay_info[0].z_e = atoms.z_min;
		atoms.spec_lay_info[0].sli_thk = 2.0;
		atoms.spec_lay_info[0].type = mt::eSLT_Top;
		atoms.spec_lay_info[0].region = region;
	}
	else
	{
		atoms.spec_lay_info[0].z_0 = atoms.z_max;
		atoms.spec_lay_info[0].z_e = atoms.z_max+bs_z;
		atoms.spec_lay_info[0].sli_thk = 2.0;
		atoms.spec_lay_info[0].type = mt::eSLT_Bottom;
		atoms.spec_lay_info[0].region = region;
	}

	mt::Amorp_Build<dt_float32> spec;
	spec.create(atoms, r_min, Z, rms_3d, rho, seed);

	auto atomsM = mex_create_pVctr<pMLD>(atoms.size(), 8, plhs[0]);
	for(auto idx = 0; idx<atoms.size(); idx++)
	{
		atomsM.real[0*atomsM.rows+idx] = atoms.Z[idx];
		atomsM.real[1*atomsM.rows+idx] = atoms.x[idx];
		atomsM.real[2*atomsM.rows+idx] = atoms.y[idx];
		atomsM.real[3*atomsM.rows+idx] = atoms.z[idx];
		atomsM.real[4*atomsM.rows+idx] = atoms.sigma[idx];
		atomsM.real[5*atomsM.rows+idx] = atoms.occ[idx];
		atomsM.real[6*atomsM.rows+idx] = atoms.region[idx];
		atomsM.real[7*atomsM.rows+idx] = atoms.charge[idx];
	}
}