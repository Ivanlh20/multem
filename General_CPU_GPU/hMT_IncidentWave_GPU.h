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

#ifndef hMT_IncidentWave_GPU_H
#define hMT_IncidentWave_GPU_H

#include "hConstTypes.h"
#include <cuda.h>
#include <cufft.h>

/*************************Incident wave*************************/
class cMT_IncidentWave_GPU{
	private:
		int IdCall;
		cMT_MGP_CPU *MT_MGP_CPU;
		sGP GP;
		sLens Lens;
		double *Mp_d;
		cufftHandle PlanPsi;
		double2 *MC_h;
	public:
		void freeMemory();
		cMT_IncidentWave_GPU();
		~cMT_IncidentWave_GPU();

		void SetInputData(cMT_MGP_CPU *MT_MGP_CPU_i, sLens &Lens_i, cufftHandle &PlanPsi_i, double2 *MC_h_i=0);
		void Psi0(double2 *&Psi0);
		void Psi0(double x, double y, double2 *&Psi0);	
};

#endif