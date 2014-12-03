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

#ifndef hfxg_CPU_H
#define hfxg_CPU_H

#include "hConstTypes.h"
#include "hMT_AtomTypes_CPU.h"

class cfxg_CPU{
	private:
		int PotPar;
		cMT_AtomTypes_CPU *MT_AtomTypes_CPU;
		double cl, cnl;	
		double ft, t, g2;
		void fxg(double g, double &f, double &df);
	public:
		void SetAtomT(int PotPar_i, cMT_AtomTypes_CPU *MT_AtomTypes_CPU_i);
		void fxg(int ng, double *g, double *f, double *df);
		cfxg_CPU();
		~cfxg_CPU();
};

#endif