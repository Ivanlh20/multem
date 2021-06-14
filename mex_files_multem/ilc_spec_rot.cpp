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

#include <mex.h>
#include "matlab_mex.cuh"

using mt::pMLD;

void mexFunction(dt_int32 nlhs, mxArray* plhs[], dt_int32 nrhs, const mxArray* prhs[])
{	
	using T = dt_float64;

	/*************************Input data**************************/
	auto r_atoms = mex_get_pvctr<pMLD>(prhs[0]);
	auto theta = mex_get_num<dt_float64>(prhs[1])*mt::c_deg_2_rad;
	auto u0 = mex_get_r_3d<T>(prhs[2], R_3d<T>(0, 0, 1));
	auto rot_point_type = mex_get_num<mt::eRot_Point_Typ>(prhs[3]);
	auto p0 = mex_get_r_3d<T>(prhs[4], R_3d<T>(0, 0, 0));

	/************************Output data**************************/
	mt::Ptc_Atom<T> atoms;
	atoms.set_ptc(r_atoms.rows, r_atoms.cols, r_atoms.real);
	u0.normalize();

	if (rot_point_type == mt::erpt_Geometric_Center)
	{
		p0 = R_3d<T>(atoms.x_mean, atoms.y_mean, atoms.z_mean);
	}

	mt::rotate_atoms(atoms, theta, u0, p0);

	auto r_atoms_o = mex_create_pVctr<pMLD>(atoms.size(), 8, plhs[0]);

	for(auto i = 0; i<r_atoms_o.rows; i++)
	{
		r_atoms_o.real[0*r_atoms_o.rows+i] = atoms.Z[i];
		r_atoms_o.real[1*r_atoms_o.rows+i] = atoms.x[i];
		r_atoms_o.real[2*r_atoms_o.rows+i] = atoms.y[i];
		r_atoms_o.real[3*r_atoms_o.rows+i] = atoms.z[i];
		r_atoms_o.real[4*r_atoms_o.rows+i] = atoms.sigma[i];
		r_atoms_o.real[5*r_atoms_o.rows+i] = atoms.occ[i];
		r_atoms_o.real[6*r_atoms_o.rows+i] = atoms.region[i];
		r_atoms_o.real[7*r_atoms_o.rows+i] = atoms.charge[i];
	}
}