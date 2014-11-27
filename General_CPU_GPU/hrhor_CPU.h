/*
 * This file is part of MULTEM.
 * Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef hrhor_CPU_H
#define hrhor_CPU_H

#include "hConstTypes.h"

class crhor_CPU{
	private:
		int PotPar;
		sAtomTypesCPU AtomTypesCPU;
		double cl, cnl;	
		double ft, t, ir, r2;
		void rhor(double r, double &f, double &df);
	public:
		void SetAtomT(int PotPari, sAtomTypesCPU AtomTypesCPUi);
		void rhor(int nr, double *r, double *f, double *df);
		crhor_CPU();
		~crhor_CPU();
};

#endif