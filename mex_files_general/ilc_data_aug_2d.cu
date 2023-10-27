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

#include <random>
#include <algorithm>

#include "cgpu_classes.cuh"

#include <mex.h>
#include "matlab_mex.h"

using mt::pMLD;
using mt::pMx_c;

/************************Host device configuration************************/
class In_Dat_Aug
{
	public:
		using value_type = dt_float64;

		pMLD M;

		dt_int32 nx;
		dt_int32 ny;

		dt_float64 bs_x;
		dt_float64 bs_y;

		dt_float64 p0_x;
		dt_float64 p0_y;

		dt_int32 nx_c;
		dt_int32 ny_c;

		dt_float64 theta_0;
		dt_float64 theta_e;

		dt_float64 dx_0;
		dt_float64 dx_e;

		dt_float64 dy_0;
		dt_float64 dy_e;

		dt_float64 sx_0;
		dt_float64 sx_e;

		dt_float64 sy_0;
		dt_float64 sy_e;

		dt_float64 sim_0;
		dt_float64 sim_e;

		dt_uint32 seed;

		In_Dat_Aug():nx(0), ny(0), bs_x(0), bs_y(0), nx_c(0), ny_c(0), 
		theta_0(0), theta_e(mt::c_2pi), dx_0(-2), dx_e(2), dy_0(-2), dy_e(2), 
		sx_0(0.75), sx_e(1.25), sy_0(0.75), sy_e(1.25), sim_0(100), sim_e(1000), seed(300183), dist(0, 1) {}

		void set_seed(dt_int32 seed_i)
		{
			seed = (seed_i<1)?rd():seed_i;
			gen.seed(seed);
			dist.reset();
		}

		dt_float64 theta() 
		{
			return theta_0+(theta_e-theta_0)*rand();
		};

		mt::R_2d<dt_float64> ps() 
		{
			auto x = dx_0+(dx_e-dx_0)*rand();
			auto y = dy_0+(dy_e-dy_0)*rand();
			return mt::R_2d<dt_float64>(x, y);
		};

		dt_float64 sx() 
		{
			return sx_0+(sx_e-sx_0)*rand();
		};

		dt_float64 sy() 
		{
			return sy_0+(sy_e-sy_0)*rand();
		};

		dt_float64 sim() 
		{
			return sim_0+(sim_e-sim_0)*rand();
		};

	private:
		std::random_device rd;
		std::mt19937_64 gen;
		std::uniform_real_distribution<dt_float64> dist;

		dt_float64 rand(){ return dist(gen); }
};

void read_in_dat_aug(const mxArray* mex_in_dat_aug, In_Dat_Aug &in_dat_aug)
{
	in_dat_aug.M = mex_get_pvctr_from_field<pMLD>(mex_in_dat_aug, "M");
	in_dat_aug.nx = in_dat_aug.M.cols;
	in_dat_aug.ny = in_dat_aug.M.rows;
	in_dat_aug.bs_x = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "bs_x");
	in_dat_aug.bs_y = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "bs_y");
	in_dat_aug.p0_x = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "p0_x")-1;
	in_dat_aug.p0_y = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "p0_y")-1;
	in_dat_aug.nx_c = mex_get_num_from_field<dt_int32>(mex_in_dat_aug, "nx_c");
	in_dat_aug.ny_c = mex_get_num_from_field<dt_int32>(mex_in_dat_aug, "ny_c");
	in_dat_aug.theta_0 = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "theta_0")*mt::c_deg_2_rad;
	in_dat_aug.theta_e = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "theta_e")*mt::c_deg_2_rad;
	in_dat_aug.dx_0 = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "dx_0");
	in_dat_aug.dx_e = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "dx_e");
	in_dat_aug.dy_0 = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "dy_0");
	in_dat_aug.dy_e = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "dy_e");
	in_dat_aug.sx_0 = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "sx_0");
	in_dat_aug.sx_e = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "sx_e");
	in_dat_aug.sy_0 = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "sy_0");
	in_dat_aug.sy_e = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "sy_e");
	in_dat_aug.sim_0 = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "sim_0");
	in_dat_aug.sim_e = mex_get_num_from_field<dt_float64>(mex_in_dat_aug, "sim_e");
 auto seed = mex_get_num_from_field<dt_int32>(mex_in_dat_aug, "seed");
	in_dat_aug.set_seed(seed);
}

template <class T, mt::eDev Dev>
void run_data_aug_2d(mt::System_Config &system_config, In_Dat_Aug &in_dat_aug, pMLD &rM_o)
{
	mt::Grid_2d<T> grid_2d(in_dat_aug.nx, in_dat_aug.ny, in_dat_aug.bs_x, in_dat_aug.bs_y);
	mt::Stream<Dev> stream(system_config.n_stream);

	mt::Vctr<T, Dev> M(in_dat_aug.M.begin(), in_dat_aug.M.end());
	mt::Vctr<T, Dev> mx_o(rM_o.size());
	mt::R_2d<T> p0(in_dat_aug.p0_x, in_dat_aug.p0_y);

	mt::Data_Aug_2d<T, Dev> data_aug_2d(&stream, grid_2d);
	data_aug_2d(M, in_dat_aug.theta(), p0, in_dat_aug.sx(), in_dat_aug.sy(), in_dat_aug.ps(), in_dat_aug.sim(), in_dat_aug.nx_c, in_dat_aug.ny_c, mx_o);

	thrust::copy(mx_o.begin(), mx_o.end(), rM_o.begin());
}

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[]) 
{
	auto system_config = mex_read_system_config(prhs[0]);
	system_config.set_gpu();

	dt_int32 idx_0 = (system_config.active)?1:0;

	In_Dat_Aug in_dat_aug;
	read_in_dat_aug(prhs[idx_0+0], in_dat_aug);

	/***************************************************************************************/
	auto rM_o = mex_create_pVctr<pMLD>(in_dat_aug.ny_c, in_dat_aug.nx_c, plhs[0]);

	if (system_config.is_float32_cpu())
	{
		run_data_aug_2d<dt_float32, mt::edev_cpu>(system_config, in_dat_aug, rM_o);
	}
	else if (system_config.is_float64_cpu())
	{
		run_data_aug_2d<dt_float64, mt::edev_cpu>(system_config, in_dat_aug, rM_o);
	}
	// else if (system_config.is_float32_gpu())
	// {
	// 	run_data_aug_2d<dt_float32, mt::edev_gpu>(system_config, in_dat_aug, rM_o);
	// }
	// else if (system_config.is_float64_gpu())
	// {
	// 	run_data_aug_2d<dt_float64, mt::edev_gpu>(system_config, in_dat_aug, rM_o);
	// }
}