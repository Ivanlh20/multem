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

#ifndef hMT_MicroscopeEffects_CPU_H
#define hMT_MicroscopeEffects_CPU_H

#include "fftw3.h"
#include "hConstTypes.h"

/*************************Incident wave*************************/
class cMT_MicroscopeEffects_CPU{
	private:
		int IdCall;

		sQ1 Qt;
		int nQs;
		sQ2 Qs;
		sGP GP;
		sLens Lens;
		fftw_plan PlanForward;
		fftw_plan PlanBackward;
		fftw_complex *Psit;

		void ReadTemporalQuadrature(sLens &Lens, sQ1 &Qt);
		void ReadSpatialQuadrature(sLens &Lens, int &nQs, sQ2 &Qs);
		void PCTCCTEM(int STEffect, fftw_complex *&fPsi, double *&M2PsiM);
		void PCLIMWPOTEM(int STEffect, fftw_complex *&fPsi, double *&M2PsiM);
	public:
		void freeMemory();
		cMT_MicroscopeEffects_CPU();
		~cMT_MicroscopeEffects_CPU();

		void SetInputData(sGP &GP_i, sLens &Lens_i, fftw_plan &PlanForward_i, fftw_plan &PlanBackward_i, fftw_complex *&Psit_i);
		void ApplyMEffects(int MEffect, int STEffect, fftw_complex *&fPsi, double *&M2Psi);
};

#endif