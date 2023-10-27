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
void run_fd_tr_rot_2d(mt::System_Config &system_config, pMLD &rM1_i, 
pMLD &rM2_i, T bs_x, T bs_y, pMLD &raf_i, T p, T sigma_g, T radius, 
pMLD &rbd, mt::eFil_Sel_Typ bg_opt, T bg, dt_int32 nit_nm, pMLD &raf_o)
{
	mt::Grid_2d<T> grid_2d(rM1_i.cols, rM1_i.rows, bs_x, bs_y);
	mt::Stream<Dev> stream(system_config.n_stream);
	mt::FFT<T, Dev> fft_2d;

	T theta = raf_i[0]*mt::c_deg_2_rad;
	mt::R_2d<T> txy(raf_i[1], raf_i[2]);
	mt::Region _Rect_2d<T> bd(grid_2d.bs_x, grid_2d.bs_y, rbd.size(), rbd.real);

	mt::Vctr<T, Dev> M1(rM1_i.begin(), rM1_i.end());
	mt::Vctr<T, Dev> M2(rM2_i.begin(), rM2_i.end());

	sigma_g = sigma_g*grid_2d.dg_min();		// pixels to g
	radius = mt::fcn_is_zero(radius)?mt::fcn_sigma_r_2_sigma_g(sigma_g):radius;

	mt::Fd_Tr_Rot_2d<T, Dev> fd_tr_rot_2d(&stream, &fft_2d, grid_2d, bg_opt, bg);
	fd_tr_rot_2d.set_fft_plan();
	fd_tr_rot_2d(M1, M2, theta, txy, p, sigma_g, radius, bd, nit_nm);

	fft_2d.cleanup();

	raf_o[0] = theta/mt::c_deg_2_rad;
	raf_o[1] = txy.x;
	raf_o[2] = txy.y;
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
	auto raf_i = mex_get_pvctr<pMLD>(prhs[idx_0+4]);
	auto p = mex_get_num<dt_float64>(prhs[idx_0+5]);
	auto sigma_g = mex_get_num<dt_float64>(prhs[idx_0+6]);
	auto rbd = (nrhs>idx_0+7)?mex_get_pvctr<pMLD>(prhs[idx_0+7]):pMLD();
	auto bg_opt = (nrhs>idx_0+8)?mex_get_num<mt::eFil_Sel_Typ>(prhs[idx_0+8]):mt::efst_min;
	auto bg = (nrhs>idx_0+9)? mex_get_num<dt_float64>(prhs[idx_0+9]):0;
	dt_int32 nit_nm = (nrhs>idx_0+10)?mex_get_num<dt_int32>(prhs[idx_0+10]):10;
	auto radius = (nrhs>idx_0+11)? mex_get_num<dt_float64>(prhs[idx_0+11]):0;

	/***************************************************************************************/
	dt_int32 rows_o = (raf_i.rows>raf_i.cols)?3:1;
	dt_int32 cols_o = (raf_i.rows>raf_i.cols)?1:3;

	auto raf_o = mex_create_pVctr<pMLD>(rows_o, cols_o, plhs[0]);

	if (system_config.is_float32_cpu())
	{
		run_fd_tr_rot_2d<dt_float32, mt::edev_cpu>(system_config, rM_1i, rM_2i, bs_x, bs_y, raf_i, p, sigma_g, radius, rbd, bg_opt, bg, nit_nm, raf_o);
	}
	else if (system_config.is_float64_cpu())
	{
		run_fd_tr_rot_2d<dt_float64, mt::edev_cpu>(system_config, rM_1i, rM_2i, bs_x, bs_y, raf_i, p, sigma_g, radius, rbd, bg_opt, bg, nit_nm, raf_o);
	}
	if (system_config.is_float32_gpu())
	{
		run_fd_tr_rot_2d<dt_float32, mt::edev_gpu>(system_config, rM_1i, rM_2i, bs_x, bs_y, raf_i, p, sigma_g, radius, rbd, bg_opt, bg, nit_nm, raf_o);
	}
	else if (system_config.is_float64_gpu())
	{
		run_fd_tr_rot_2d<dt_float64, mt::edev_gpu>(system_config, rM_1i, rM_2i, bs_x, bs_y, raf_i, p, sigma_g, radius, rbd, bg_opt, bg, nit_nm, raf_o);
	}
}