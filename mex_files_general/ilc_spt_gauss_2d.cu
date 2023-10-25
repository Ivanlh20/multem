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

#include "in_classes.cuh"
#include "fcns_gpu.h"
#include "cgpu_spt_ptc.cuh"

#include <mex.h>
#include "matlab_mex.h"

template <class T>
void read_in_spt_gauss_2d(const mxArray* mex_in_spt_gauss_2d, mt::In_Spt_Ptc_2d<T>& in_spt_gauss_2d)
{
	in_spt_gauss_2d.scf_radius = mex_get_num_from_field<T>(mex_in_spt_gauss_2d, "scf_radius");

	dt_bool bwl = false;
	dt_bool pbc_xy = false;

	auto nx = mex_get_num_from_field<dt_int32>(mex_in_spt_gauss_2d, "nx");
	auto ny = mex_get_num_from_field<dt_int32>(mex_in_spt_gauss_2d, "ny");
	auto dx = mex_get_num_from_field<T>(mex_in_spt_gauss_2d, "dx");
	auto dy = mex_get_num_from_field<T>(mex_in_spt_gauss_2d, "dy");
	auto bs_x = T(nx*dx);
	auto bs_y = T(ny*dy);
	auto sli_thick = 0.5;
	dt_int64 icol = 0;

	const auto ed_typ = mex_get_data_type_from_field(mex_in_spt_gauss_2d, "data");
	if (mt::is_efloat32(ed_typ))
	{
		auto ptc = mex_get_pvctr_from_field<dt_float32>(mex_in_spt_gauss_2d, "data");
		in_spt_gauss_2d.ptc.set_ptc(ptc, icol, {bs_x, bs_y}, pbc_xy);
	}
	else if (mt::is_efloat64(ed_typ))
	{
		auto ptc = mex_get_pvctr_from_field<dt_float64>(mex_in_spt_gauss_2d, "data");
		in_spt_gauss_2d.ptc.set_ptc(ptc, icol, {bs_x, bs_y}, pbc_xy);
	}

	in_spt_gauss_2d.grid.set_in_data(bs_x, bs_y, nx, ny, T(0), T(0), pbc_xy, bwl, sli_thick);
 }

template <class T, mt::eDev Dev>
void mex_run(mt::System_Config& system_config, dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	mt::In_Spt_Ptc_2d<T> in_spt_gauss_2d;
	read_in_spt_gauss_2d(prhs[system_config.idx_0+0], in_spt_gauss_2d);

	/***************************************************************************************/
	auto pmx_o = mex_create_pVctr<T>(in_spt_gauss_2d.grid.shape(), plhs[0]);

	mt::Stream<Dev> stream(system_config.n_stream);

	mt::Vctr<T, Dev> mx_o(pmx_o.shape());

	mt::Spt_gauss_2d<T, Dev> spt_gauss_2d(&in_spt_gauss_2d, &stream);

	spt_gauss_2d(mx_o);

	mx_o.cpy_to_cpu_ptr(pmx_o.begin(), pmx_o.end());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto system_config = mex_read_set_system_config(prhs[0]);				
			
	if (system_config.is_float32_cpu())
	{
		mex_run<dt_float32, mt::edev_cpu>(system_config, nlhs, plhs, nrhs, prhs);
	}
	else if (system_config.is_float64_cpu())
	{
		mex_run<dt_float64, mt::edev_cpu>(system_config, nlhs, plhs, nrhs, prhs);
	}
	else if (system_config.is_float32_gpu())
	{
		mex_run<dt_float32, mt::edev_gpu>(system_config, nlhs, plhs, nrhs, prhs);
	}				
	else if (system_config.is_float64_gpu())
	{
		mex_run<dt_float64, mt::edev_gpu>(system_config, nlhs, plhs, nrhs, prhs);
	}

	//MEX_RUN_FCN_FLOAT_SYS_CONF(mex_run, 0);
}