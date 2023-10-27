/*
 * This file is part of Multem.
 * Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
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

#include "math_mt.h"
#include "types.cuh"
#include "type_traits_gen.h"
#include "cgpu_stream.cuh"

#include "cgpu_classes.cuh"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;
using mt::pMx_c;

template <class T, mt::eDev Dev>
void run_fd_sft_2d_br(mt::System_Config &system_config, pMLD &rM1_i, pMLD &rM2_i, 
T bs_x, T bs_y, T p, T sigma_g, pMLD &rdx_i, pMLD &rbd, pMLD &rdx_o)
{
	mt::Grid_2d<T> grid_2d(rM1_i.rows, rM2_i.cols, bs_y, bs_x);
	mt::Stream<Dev> stream(system_config.n_stream);
	mt::FFT<T, Dev> fft_2d;
	mt::Region _Rect_2d<T> bd(grid_2d.bs_x, grid_2d.bs_y, rbd.size(), rbd.real, true);

	mt::Vctr<T, Dev> M1(rM1_i.begin(), rM1_i.end());
	mt::Vctr<T, Dev> M2(rM2_i.begin(), rM2_i.end());
	mt::Vctr<T, mt::edev_cpu> dx(grid_2d.nx, T(0));
	if (rdx_i.size()>=grid_2d.nx)
	{
		dx.assign(rdx_i.begin(), rdx_i.end());
	}

	sigma_g = sigma_g*grid_2d.dgy;		// pixels to g

	mt::fcn_trs_2d(stream, grid_2d.nx, grid_2d.ny, M1);
	mt::fcn_trs_2d(stream, grid_2d.nx, grid_2d.ny, M2);

	mt::Fd_Sft_2d_BC<T, Dev> fd_sft_2d_bc(&stream, &fft_2d, grid_2d);
	fd_sft_2d_bc.set_fft_plan();
	auto dx_o = fd_sft_2d_bc(M1, M2, p, sigma_g, dx, bd, 2);
	fd_sft_2d_bc.cleanup();

	thrust::copy(dx_o.begin(), dx_o.end(), rdx_o.begin());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto system_config = mex_read_system_config(prhs[0]);
	system_config.set_gpu();

	dt_int32 idx_0 = (system_config.active)?1:0;

	auto rM_1i = mex_get_pvctr<pMLD>(prhs[idx_0+0]);
	auto rM_2i = mex_get_pvctr<pMLD>(prhs[idx_0+1]);
	auto bs_x = rM_1i.cols*mex_get_num<dt_float64>(prhs[idx_0+2]);
	auto bs_y = rM_1i.rows*mex_get_num<dt_float64>(prhs[idx_0+3]);
	auto p = mex_get_num<dt_float64>(prhs[idx_0+4]);
	auto sigma_g = mex_get_num<dt_float64>(prhs[idx_0+5]);
	auto rdx = (nrhs>idx_0+6)?mex_get_pvctr<pMLD>(prhs[idx_0+6]):pMLD();
	auto rbd = (nrhs>idx_0+7)?mex_get_pvctr<pMLD>(prhs[idx_0+7]):pMLD();

	/***************************************************************************************/
	auto rdx_o = mex_create_pVctr<pMLD>(rM_1i.rows, 1, plhs[0]);

	if (system_config.is_float32_cpu())
	{
		run_fd_sft_2d_br<dt_float32, mt::edev_cpu>(system_config, rM_1i, rM_2i, bs_x, bs_y, p, sigma_g, rdx, rbd, rdx_o);
	}
	else if (system_config.is_float64_cpu())
	{
		run_fd_sft_2d_br<dt_float64, mt::edev_cpu>(system_config, rM_1i, rM_2i, bs_x, bs_y, p, sigma_g, rdx, rbd, rdx_o);
	}
	else if (system_config.is_float32_gpu())
	{
		run_fd_sft_2d_br<dt_float32, mt::edev_gpu>(system_config, rM_1i, rM_2i, bs_x, bs_y, p, sigma_g, rdx, rbd, rdx_o);
	}
	else if (system_config.is_float64_gpu())
	{
		run_fd_sft_2d_br<dt_float64, mt::edev_gpu>(system_config, rM_1i, rM_2i, bs_x, bs_y, p, sigma_g, rdx, rbd, rdx_o);
	}
}