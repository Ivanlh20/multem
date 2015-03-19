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

#include "math.h"
#include <cstring>

#include "hConstTypes.h"
#include "hMT_AtomicData_CPU.h"
#include "hMT_General_CPU.h"
#include "hMT_AtomTypes_CPU.h"

// free memory
void cMT_AtomTypes_CPU::freeMemory()
{
	if(IdCall==0)
	{
		return;
	}

	Z = 0;
	m = 0;
	A = 0;
	rn_e = 0;
	rn_c = 0;
	ra_e = 0;
	ra_c = 0;
	Rmin = 0;
	Rmax = 0;
	Rmin2 = 0;
	Rmax2 = 0;

	f_sCoefPar_Free_CPU(cfeg);
	f_sCoefPar_Free_CPU(cfxg);
	f_sCoefPar_Free_CPU(cPr);
	f_sCoefPar_Free_CPU(cVr);
	f_sCoefPar_Free_CPU(cVR);
	
	nR = 0;
	delete [] R; R = 0;
	delete [] R2; R2 = 0;
	f_sciVn_Free_CPU(ciVR);
}

// Set Atom type
void cMT_AtomTypes_CPU::SetAtomTypes(int Z_i, int PotPar_i, double Vrl_i, int nR_i, double Rmin_i)
{
	freeMemory();	// clean CPU memory
	IdCall++;

	f_sCoefPar_Malloc_CPU(6, cfeg);
	f_sCoefPar_Malloc_CPU(6, cfxg);
	f_sCoefPar_Malloc_CPU(6, cPr);
	f_sCoefPar_Malloc_CPU(6, cVr);
	f_sCoefPar_Malloc_CPU(6, cVR);
	R = new double [nR_i];
	R2 = new double [nR_i];
	f_sciVn_Malloc_CPU(nR_i, ciVR);

	cMT_AtomicData_CPU MT_AtomicData_CPU(PotPar_i);
	MT_AtomicData_CPU.To_MT_AtomTypes_CPU(Z_i, Vrl_i, nR_i, Rmin_i, this);
}