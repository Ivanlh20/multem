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

#include <algorithm>
#include "fxeg_tabulated_data.h"

#include <mex.h>
#include "mex_matlab.hpp"

using multem::rmatrix_r;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int Z, ng, type;
	double g[1024], g2[1024], fxg[1024], feg[1024];
	fxeg_Tabulated_Data fxeg_tabulated_data;

	Z = mx_get_scalar<int>(prhs[0]); 
	type = mx_get_scalar<int>(prhs[1]); 
	fxeg_tabulated_data.ReadTabData(Z, type, 1, ng, g, g2, fxg, feg);
	
	rmatrix_r g_o, fxg_o, feg_o;

	g_o = mx_create_matrix<rmatrix_r>(ng, 1, plhs[0]);
	fxg_o = mx_create_matrix<rmatrix_r>(ng, 1, plhs[1]);
	feg_o = mx_create_matrix<rmatrix_r>(ng, 1, plhs[2]);

	std::copy(g, g + ng, g_o.real);
	std::copy(fxg, fxg + ng, fxg_o.real);
	std::copy(feg, feg + ng, feg_o.real);
}