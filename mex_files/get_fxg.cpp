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
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#include "types.hpp"
#include "atom_cal.hpp"
#include "atomic_data.hpp"

#include <mex.h>
#include "mex_matlab.hpp"

using multem::rmatrix_r;

void mexFunction(int nlhs,mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	auto potential_type = mx_get_scalar<multem::ePotential_Type>(prhs[0]);
	auto Z = mx_get_scalar<int>(prhs[1]);
	auto g = mx_get_matrix<rmatrix_r>(prhs[2]);

	auto fxg = mx_create_matrix<rmatrix_r>(g.rows, g.cols, plhs[0]);
	auto dfxg = mx_create_matrix<rmatrix_r>(g.rows, g.cols, plhs[1]);

	multem::Atom_Type<double, multem::e_Host> atom_type;
	multem::Atomic_Data atomic_data;
	atomic_data.Load_Data(potential_type);
	atomic_data.To_atom_type_CPU(Z, multem::c_Vrl, multem::c_nR, 0.0, atom_type);

	multem::Atom_Cal<double> atom_cal;
	atom_cal.Set_Atom_Type(potential_type, &atom_type);
	atom_cal.fxg_dfxg(g.size, g.real, fxg.real, dfxg.real);
}