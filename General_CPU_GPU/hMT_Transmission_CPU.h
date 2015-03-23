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

#ifndef hMT_Transmission_CPU_H
#define hMT_Transmission_CPU_H

#include "fftw3.h"
#include "hConstTypes.h"
#include "hQuadrature.h"
#include "hMT_General_CPU.h"
#include "hMT_MGP_CPU.h"
#include "hMT_Specimen_CPU.h"
#include "hMT_AtomTypes_CPU.h"
#include "hMT_Potential_CPU.h"

class cMT_Transmission_CPU: public cMT_Potential_CPU{
	private:
		int IdCall;

		double fPot;

		int SliceMemTyp;
		int nSliceMem;
		int nSliceMem0;

		fftw_complex **Trans;
		float **Vpe;

		fftw_plan PlanForward;
		fftw_plan PlanBackward;
		void Cal_Trans_or_Vpe();
	public:
		fftw_complex *Trans0;

		void freeMemory();

		cMT_Transmission_CPU();
		~cMT_Transmission_CPU();

		void SetInputData(cMT_MGP_CPU *MT_MGP_CPU_io, fftw_plan &PlanForward_i, fftw_plan &PlanBackward_i, int nAtomsM_i, double *AtomsM_i);
		void MoveAtoms(int iConf);	
		fftw_complex* getTrans(int iSlice);
		void Transmit(int iSlice, fftw_complex *&Psi_io);
};

#endif