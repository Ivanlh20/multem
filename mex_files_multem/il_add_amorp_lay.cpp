/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * MULTEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MULTEM is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#include "types.cuh"
#include "matlab_types.cuh"
#include "atomic_data_mt.hpp"
#include "amorp_spec.hpp"

#include <mex.h>
#include "matlab_mex.cuh"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ])
{
	auto ratoms = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto l_x = mx_get_scalar<double>(prhs[1]);
	auto l_y = mx_get_scalar<double>(prhs[2]);
	auto l_z = mx_get_scalar<double>(prhs[3]);
	auto r_min = mx_get_scalar<double>(prhs[4]);
	auto Z = mx_get_scalar<int>(prhs[5]);
	auto rms_3d = mx_get_scalar<double>(prhs[6]);
	auto rho = mx_get_scalar<double>(prhs[7]);
	auto lay_pos = (nrhs>8)?mx_get_scalar<mt::eAmorp_Lay_Type>(prhs[8]):mt::eALT_Top;
	int seed = (nrhs>9)?mx_get_scalar<int>(prhs[9]):300183;

	/*******************************************************************/
	mt::Atom_Data<float> atoms;
	atoms.set_atoms(ratoms.rows, ratoms.cols, ratoms.real, l_x, l_y);

	int region = (ratoms.rows==0)?0:atoms.region_max+1;

	atoms.amorp_lay_info.resize(1);
	if(lay_pos==mt::eALT_Top)
	{
		atoms.amorp_lay_info[0].z_0 = atoms.z_min-l_z;
		atoms.amorp_lay_info[0].z_e = atoms.z_min;
		atoms.amorp_lay_info[0].dz = 2.0;
		atoms.amorp_lay_info[0].type = mt::eALT_Top;
		atoms.amorp_lay_info[0].region = region;
	}
	else
	{
		atoms.amorp_lay_info[0].z_0 = atoms.z_max;
		atoms.amorp_lay_info[0].z_e = atoms.z_max+l_z;
		atoms.amorp_lay_info[0].dz = 2.0;
		atoms.amorp_lay_info[0].type = mt::eALT_Bottom;
		atoms.amorp_lay_info[0].region = region;
	}

	mt::Amorp_Spec<float> spec;
	spec.create(atoms, r_min, Z, rms_3d, rho, seed);

	auto atomsM = mx_create_matrix<rmatrix_r>(atoms.size(), 8, plhs[0]);
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