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
#include "vctr_cpu.h"
#include "vctr_gpu.h"
#include "stream_cpu.h"
#include "stream_gpu.h"
#include "fft_cpu.h"
#include "fft_gpu.h"
#include "fcns_cpu.h"
#include "fcns_gpu.h"

#include <mex.h>
#include "matlab_mex.h"

template <class T, mt::eDev Dev>
void mex_run(mt::System_Config& system_config, dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto mx_i = mex_get_vctr<T>(prhs[system_config.idx_0+0]);
	auto bs = mex_get_r_2d_bs<T>(prhs[system_config.idx_0+1], mx_i.shape());
	auto sigma_r = (nrhs>system_config.idx_0+2)?mex_get_num<T>(prhs[system_config.idx_0+2]):T(1);

	/***************************************************************************************/
	auto pmx_o = mex_create_pVctr<T>(mx_i.shape(), plhs[0]);

	mt::Vctr<complex<T>, Dev> mx_io(mx_i);

	auto sigma_g = mt::fcn_sigma_r_2_sigma_g(sigma_r);
	mt::Fcn_Gauss<T> fcn_gauss(sigma_g);
	mt::Stream<Dev> stream(system_config.n_stream);

	if (mx_i.is_1d())
	{
		mt::Grid_1d<T> grid(bs.x, mx_i.size_32());
		mt::FFT<T, Dev> fft(grid.nx, &stream);
		mt::fcn_cv_fs_gauss_1d(grid, fft, fcn_gauss, mx_io, &stream);
		fft.cleanup();
	}
	else
	{
		mt::Grid_2d<T> grid(bs.x, bs.y, mx_i.s1_32(), mx_i.s0_32());
		mt::FFT<T, Dev> fft(grid.ny, grid.nx, &stream);
		mt::fcn_cv_fs_gauss_2d(grid, fft, fcn_gauss, mx_io, &stream);
		fft.cleanup();
	}

	mx_io.cpy_real_to_cpu_ptr(pmx_o.begin(), pmx_o.end());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	MEX_RUN_FCN_FLOAT_SYS_CONF(mex_run, 0);
}