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

#include "fftw3.h"

#include <cmath>
#include "hConstTypes.h"
#include "hMT_General_CPU.h"
#include "hMT_IncidentWave_CPU.h"

void cMT_IncidentWave_CPU::freeMemory()
{
	if(IdCall==0) return;

	MT_MGP_CPU = 0;
	f_sGP_Init(GP);
	f_sLens_Init(Lens);
}

cMT_IncidentWave_CPU::cMT_IncidentWave_CPU()
{
	IdCall = 0;

	MT_MGP_CPU = 0;
	f_sGP_Init(GP);
	f_sLens_Init(Lens);
}

cMT_IncidentWave_CPU::~cMT_IncidentWave_CPU()
{
	freeMemory();
	IdCall = 0;
}

void cMT_IncidentWave_CPU::SetInputData(cMT_MGP_CPU *MT_MGP_CPU_i, sLens &Lens_i, fftw_plan &PlanForward_i, fftw_plan &PlanBackward_i)
{
	freeMemory();
	MT_MGP_CPU = MT_MGP_CPU_i;
	f_sGP_SetInputData(MT_MGP_CPU, GP);
	Lens = Lens_i;
	PlanForward = PlanForward_i;
	PlanBackward = PlanBackward_i;
}

void cMT_IncidentWave_CPU::Psi0(fftw_complex *&Psi0)
{
	if(MT_MGP_CPU->Psi0Typ==1)
	{
		f_Set_MC_CPU(GP, 1.0, 0.0, Psi0);
	}
	else
	{
		f_sComplex_2_fftw_complex_CPU(GP, MT_MGP_CPU->Psi0, Psi0);
		// fft2shift 
		f_fft2Shift_MC_CPU(GP, Psi0);	
	}
}

void cMT_IncidentWave_CPU::Psi0(double x, double y, fftw_complex *&Psi0)
{
	f_Probe_FS_CPU(GP, Lens, c2Pi*(0.5*GP.lx-x), c2Pi*(0.5*GP.ly-y), Psi0);
	// Backward fft2
	fftw_execute_dft(PlanBackward, Psi0, Psi0);
	double factor = sqrt(double(GP.nxy)/f_Sum_MC2_CPU(GP, 1.0, Psi0));
	f_Scale_MC_CPU(GP, factor, Psi0);
}