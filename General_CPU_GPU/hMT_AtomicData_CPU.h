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

#ifndef hMT_AtomicData_CPU_H
#define hMT_AtomicData_CPU_H

#include "hConstTypes.h"
#include "hPotential_CPU.h"
#include "hMT_AtomTypes_CPU.h"

class cMT_AtomicData_CPU{
	private:
		int PotPar;
		cPotential_CPU Potential_CPU;
		sAtomicData *AtomicData;
		sCoefPar *fegPar;

		void Load_Atomic_Data(sAtomicData *AtomicData);
		void Load_feg_Doyle_0_4(sCoefPar *fegPar);
		void Load_feg_Peng_0_4(sCoefPar *fegPar);
		void Load_feg_Peng_0_12(sCoefPar *fegPar);
		void Load_feg_Kirkland_0_12(sCoefPar *fegPar);
		void Load_feg_Weickenmeier_0_12(sCoefPar *fegPar);
		void Load_feg_Lobato_0_12(sCoefPar *fegPar);

		void getCoef(int Z, sCoefPar &cfeg, sCoefPar &cfxg, sCoefPar &cPr, sCoefPar &cVr, sCoefPar &cVR);
public:
		cMT_AtomicData_CPU(int PotPar_i);
		~cMT_AtomicData_CPU();
		void Load_Data(int PotPar_i);	
		void To_MT_AtomTypes_CPU(int Z_i, double Vrl_i, int nR_i, double Rmin_i, cMT_AtomTypes_CPU *MT_AtomTypes_CPU);
};

#endif