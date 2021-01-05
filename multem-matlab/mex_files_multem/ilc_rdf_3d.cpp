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
#include "cpu_fcns.hpp"

#include <mex.h>
#include "matlab_mex.cuh"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[ ])
{
	auto ratoms = mx_get_matrix<rmatrix_r>(prhs[0]);
	auto r_max = mx_get_scalar<double>(prhs[1]);
	auto nr = (nrhs>2)?mx_get_scalar<int>(prhs[2]):10;

	/*******************************************************************/
	auto rr_o = mx_create_matrix<rmatrix_r>(nr, 1, plhs[0]);
	auto rrdf_o = mx_create_matrix<rmatrix_r>(nr, 1, plhs[1]);

	std::vector<float> r(nr);
	std::vector<float> rdf(nr);
	mt::Atom_Data<float> atoms; 
	atoms.set_atoms(ratoms.rows, ratoms.cols, ratoms.real);

	mt::rdf_3d(atoms, r_max, nr, r, rdf);

	thrust::copy(r.begin(), r.end(), rr_o.begin());
	thrust::copy(rdf.begin(), rdf.end(), rrdf_o.begin());
}