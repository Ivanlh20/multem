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
#include "atomic_data.hpp"

#include <mex.h>
#include "matlab2cpp.hpp"

using multem::m_matrix_r;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	auto potential_type = mx_get_scalar<multem::ePotential_Type>(prhs[0]);
	auto Dim = mx_get_scalar<int>(prhs[1]);
	auto Vrl = mx_get_scalar<double>(prhs[2]);

	auto r = mx_create_matrix<m_matrix_r>(multem::c_nAtomsTypes, 3, plhs[0]);

	multem::Atom_Cal<double> atom_cal;
	multem::Atomic_Data atomic_data;
	atomic_data.Load_Data(potential_type);
	multem::Atom_Type<double, multem::Host> atom_type;

	for(auto i=0; i<r.rows; i++)
	{
		atomic_data.To_atom_type_CPU(i+1, multem::c_Vrl, multem::c_nR, 0.0, atom_type);
		atom_cal.Set_Atom_Type(potential_type, &atom_type);
		r.real[i+0*r.rows] = atom_cal.AtomicRadius_rms(Dim);
		r.real[i+1*r.rows] = atom_cal.AtomicRadius_Cutoff(Dim, Vrl);
		r.real[i+2*r.rows] = atom_type.ra_e;
	}
}