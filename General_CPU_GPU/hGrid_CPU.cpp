/**
 *  This file is part of MULTEM.
 *  Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
 *
 *  MULTEM is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MULTEM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MULTEM.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "hmathCPU.h"
#include "hGrid_CPU.h"

// Build logarithm grid
void cGrid_CPU::RegularGrid(double rmin, double rmax, int nr, double *r){
	double dr = (rmax-rmin)/(nr-1);
	for (int i=0; i<nr; i++)
		r[i] = rmin + i*dr;
}

// Build logarithm grid
void cGrid_CPU::LogarithmicGrid(double rmin, double rmax, int nr, double *r){
	double dlnr = log(rmax/rmin)/double(nr-1);
	for (int i=0; i<nr; i++)
		r[i] = rmin*exp(i*dlnr);
}

// Build logarithm grid
void cGrid_CPU::LogarithmicGridShifted(double rmin, double rmax, int nr, double *r){
	double dlnr = log(rmax/rmin + 1.0)/double(nr-1);
	for (int i=0; i<nr; i++)
		r[i] = rmin*(exp(i*dlnr)-1.0);
}

// Read logarithm grid
void cGrid_CPU::ReadGrid(double rmin, double rmax, int nr, int gridTyp, double *r){
	switch (gridTyp){
		case 0:
			RegularGrid(rmin, rmax, nr, r);
			break;
		case 1:
			LogarithmicGrid(rmin, rmax, nr, r);
			break;
		case 2:
			LogarithmicGridShifted(rmin, rmax, nr, r);
			break;
	}
}