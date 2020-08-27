/*
 * This file is part of MULTEM.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
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
#include "traits.cuh"
#include "atomic_data_mt.hpp"
#include "input_multislice.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::rmatrix_r;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
	using T = double;

	/*************************Input data**************************/
	auto r_atoms = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto theta = mx_get_scalar<double>(prhs[1])*mt::c_deg_2_rad;
	auto r_u0 = mx_get_matrix<rmatrix_r>(prhs[2]);
	auto rot_point_type = mx_get_scalar<mt::eRot_Point_Type>(prhs[3]);
	auto r_p0 = mx_get_matrix<rmatrix_r>(prhs[4]);

	/************************Output data**************************/
	mt::Atom_Data<T> atoms;
	atoms.set_atoms(r_atoms.rows, r_atoms.cols, r_atoms.real);

	mt::r3d<T> u0 = (r_u0.size()>=3)?mt::r3d<T>(r_u0[0], r_u0[1], r_u0[2]):mt::r3d<T>(0, 0, 1);
	u0.normalized();
	mt::r3d<T> p0 = (r_p0.size()>=3)?mt::r3d<T>(r_p0[0], r_p0[1], r_p0[2]):mt::r3d<T>(0, 0, 0);

	if(rot_point_type == mt::eRPT_geometric_center)
	{
		p0 = r3d<T>(atoms.x_mean, atoms.y_mean, atoms.z_mean);
	}

	mt::rotate_atoms(atoms, theta, u0, p0);

	auto r_atoms_o = mx_create_matrix<rmatrix_r>(atoms.size(), 8, plhs[0]);

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