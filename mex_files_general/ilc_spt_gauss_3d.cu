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

#include "cgpu_classes.cuh"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;

template <class TIn_Spt>
void read_in_spt(const mxArray* mex_in_spt, TIn_Spt &in_spt_gauss_3d, dt_bool full = true)
{
	using T_r = mt::Value_type<TIn_Spt>;

	in_spt_gauss_3d.scf_radius = mex_get_num_from_field<T_r>(mex_in_spt, "scf_radius");

	dt_bool bwl = false;
	dt_bool pbc_xyz = false;

	auto nx = mex_get_num_from_field<dt_int32>(mex_in_spt, "nx");
	auto ny = mex_get_num_from_field<dt_int32>(mex_in_spt, "ny");
	auto nz = mex_get_num_from_field<dt_int32>(mex_in_spt, "nz");
	auto dx = mex_get_num_from_field<T_r>(mex_in_spt, "dx");
	auto dy = mex_get_num_from_field<T_r>(mex_in_spt, "dy");
	auto sli_thick = mex_get_num_from_field<T_r>(mex_in_spt, "sli_thick");
	auto bs_x = nx*dx;
	auto bs_y = ny*dy;
	auto bs_z = nz*sli_thick;
	auto dt = 0.5;

	auto atoms = mex_get_pvctr_from_field<pMLD>(mex_in_spt, "data");
	if (full)
	{
		in_spt_gauss_3d.atoms.set_ptc(atoms.rows, atoms.cols, atoms.real, bs_x, bs_y);
	}
	in_spt_gauss_3d.grid_3d.set_in_data(nx, ny, nz, bs_x, bs_y, bs_z, dt, bwl, pbc_xyz);
 }

template <class T, mt::eDev Dev>
void run_spt_gauss_3d(mt::System_Config &system_config, const mxArray* mxA, pMLD &rM_o)
{
	mt::In_Spt_gauss_3d<T> in_spt_gauss_3d;
	read_in_spt(mxA, in_spt_gauss_3d);

	mt::Stream<Dev> stream(system_config.n_stream);

	mt::Vctr<T, Dev> mx_o(rM_o.size());

	mt::Spt_gauss_3d<T, Dev> spt_gauss_3d(&in_spt_gauss_3d, &stream);

	spt_gauss_3d(mx_o);

	thrust::copy(mx_o.begin(), mx_o.end(), rM_o.begin());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{
	auto system_config = mex_read_system_config(prhs[0]);
	system_config.set_gpu();

	dt_int32 idx_0 = (system_config.active)?1:0;

	auto nx = mex_get_num_from_field<dt_int32>(prhs[idx_0+0], "nx");
	auto ny = mex_get_num_from_field<dt_int32>(prhs[idx_0+0], "ny");
	auto nz = mex_get_num_from_field<dt_int32>(prhs[idx_0+0], "nz");
	auto rM_o = mex_create_pVctr<pMLD>(ny, nx, plhs[0]);

	if (system_config.is_float32_cpu())
	{
		run_spt_gauss_3d<dt_float32, mt::edev_cpu>(system_config, prhs[idx_0+0], rM_o);
	}
	else if (system_config.is_float64_cpu())
	{
		run_spt_gauss_3d<dt_float64, mt::edev_cpu>(system_config, prhs[idx_0+0], rM_o);
	}
	if (system_config.is_float32_gpu())
	{
		run_spt_gauss_3d<dt_float32, mt::edev_gpu>(system_config, prhs[idx_0+0], rM_o);
	}
	else if (system_config.is_float64_gpu())
	{
		run_spt_gauss_3d<dt_float64, mt::edev_gpu>(system_config, prhs[idx_0+0], rM_o);
	}
}