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

#ifndef hMT_IncidentWave_CPU_H
#define hMT_IncidentWave_CPU_H

#include "fftw3.h"
#include "hConstTypes.h"

/*************************Incident wave*************************/
class cMT_IncidentWave_CPU{
	private:
		int IdCall;
		cMT_MGP_CPU *MT_MGP_CPU;
		sGP GP;
		sLens Lens;
		fftw_plan PlanForward;
		fftw_plan PlanBackward;
	public:
		void freeMemory();
		cMT_IncidentWave_CPU();
		~cMT_IncidentWave_CPU();

		void SetInputData(cMT_MGP_CPU *MT_MGP_CPU_i, sLens &Lens_i, fftw_plan &PlanForward_i, fftw_plan &PlanBackward_i);
		void Psi0(fftw_complex *&Psi0);
		void Psi0(double x, double y, fftw_complex *&Psi0);	
};

#endif