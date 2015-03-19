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
#include "hQuadrature.h"
#include "hMT_General_CPU.h"
#include "hMT_MicroscopeEffects_CPU.h"

/***************************************************************************/
void cMT_MicroscopeEffects_CPU::freeMemory()
{
	if(IdCall==0) return;

	f_sGP_Init(GP);
	f_sLens_Init(Lens);

	delete [] Qt.x; Qt.x = 0;
	delete [] Qt.w; Qt.w = 0;

	nQs = 0;
	delete [] Qs.x; Qs.x = 0;
	delete [] Qs.y; Qs.y = 0;
	delete [] Qs.w; Qs.w = 0;

	Psit = 0;

	PlanForward = 0;
	PlanBackward = 0;
}

cMT_MicroscopeEffects_CPU::cMT_MicroscopeEffects_CPU()
{
	IdCall = 0;

	f_sGP_Init(GP);
	f_sLens_Init(Lens);

	Qt.x = 0; 
	Qt.w = 0;

	nQs = 0;
	Qs.x = 0;
	Qs.y = 0;
	Qs.w = 0;

	Psit = 0;

	PlanForward = 0;
	PlanBackward = 0;
}

cMT_MicroscopeEffects_CPU::~cMT_MicroscopeEffects_CPU()
{	
	freeMemory();
	IdCall = 0;
}

// Partially coherent transfer function and Transmission cross coefficient
void cMT_MicroscopeEffects_CPU::PCTCCTEM(int STEffect, fftw_complex *&fPsi, double *&M2PsiM)
{
	int i, j;
	double f0 = Lens.f;
	double cf0 = Lens.cf;

	f_Set_MD_CPU(GP, 0.0, M2PsiM);
	switch(STEffect)
	{
		case 1:	// Temporal and Spatial
			for(i=0; i<nQs; i++)
			{
				for(j=0; j<Lens.nsf; j++)
				{
					Lens.f = Lens.sf*Qt.x[j]+f0; 
					Lens.cf = cPi*Lens.lambda*Lens.f;
					// Apply Coherent transfer function
					f_Apply_CTF_CPU(GP, Lens, Qs.x[i], Qs.y[i], fPsi, Psit);
					// Backward fft2
					fftw_execute_dft(PlanBackward, Psit, Psit);
					// Apply weighting factor and add to the general sum
					f_Add_wMC2_CPU(GP, Qs.w[i]*Qt.w[j], Psit, M2PsiM);
				}
			}
			break;
		case 2:	// Temporal
			for(j=0; j<Lens.nsf; j++)
			{
				Lens.f = Lens.sf*Qt.x[j]+f0; 
				Lens.cf = cPi*Lens.lambda*Lens.f;
				// Apply Coherent transfer function
				f_Apply_CTF_CPU(GP, Lens, 0.0, 0.0, fPsi, Psit);
				// Backward fft2
				fftw_execute_dft(PlanBackward, Psit, Psit);
				// Apply weighting factor and add to the general sum
				f_Add_wMC2_CPU(GP, Qt.w[j], Psit, M2PsiM);
			}
			break;
		case 3:	// Spatial
			for(i=0; i<nQs; i++)
			{
				// Apply Coherent transfer function
				f_Apply_CTF_CPU(GP, Lens, Qs.x[i], Qs.y[i], fPsi, Psit);
				// Backward fft2
				fftw_execute_dft(PlanBackward, Psit, Psit);
				// Apply weighting factor and add to the general sum
				f_Add_wMC2_CPU(GP, Qs.w[i], Psit, M2PsiM);
			}
			break;
	}

	Lens.f = f0;
	Lens.cf = cf0;
}

// Partially coherent transfer function, linear image model and weak phase object
void cMT_MicroscopeEffects_CPU::PCLIMWPOTEM(int STEffect, fftw_complex *&fPsi, double *&M2PsiM)
{
	double sf = Lens.sf, beta = Lens.beta;

	switch(STEffect)
	{
		case 2:	// Temporal
			Lens.beta = 0;
			break;
		case 3:	// Spatial
			Lens.sf = 0;
			break;
	}

	f_Apply_PCTF_CPU(GP, Lens, fPsi, Psit);
	// Backward fft2
	fftw_execute_dft(PlanBackward, Psit, Psit);
	// Apply weighting factor and add to the general sum
	f_Set_wMC2_CPU(GP, 1.0, Psit, M2PsiM);

	Lens.sf = sf;
	Lens.beta = beta;
}

void cMT_MicroscopeEffects_CPU::ReadTemporalQuadrature(sLens &Lens, sQ1 &Qt)
{
	double df = 6.0*Lens.sf/double(Lens.nsf-1);
	double f, f0 = -3.0*Lens.sf;
	double sumwia = 0.0;

	Qt.x = new double[Lens.nsf]; 
	Qt.w = new double [Lens.nsf];
	for(int i=0; i<=Lens.nsf-1; i++)
	{
		Qt.x[i] = f = f0 + i*df;
		sumwia += Qt.w[i] = exp(-f*f/(Lens.sf*Lens.sf));
	}
	for(int i=0; i<=Lens.nsf-1; i++)
	{
		Qt.w[i] /= sumwia;
	}
}

void cMT_MicroscopeEffects_CPU::ReadSpatialQuadrature(sLens &Lens, int &nQs, sQ2 &Qs)
{
	int i, j;
	double gxs, gys, g2s, sumwia;
	double alpha = 0.5/pow(Lens.sggs, 2);
	sQ2 Qst;

	Qst.x = new double [(2*Lens.ngxs+1)*(2*Lens.ngys+1)];
	Qst.y = new double [(2*Lens.ngxs+1)*(2*Lens.ngys+1)];
	Qst.w = new double [(2*Lens.ngxs+1)*(2*Lens.ngys+1)];
	/**********************************************************************/
	nQs = 0; sumwia = 0.0;
	 for(j=-Lens.ngys; j<=Lens.ngys; j++)
	 {
		 for(i=-Lens.ngxs; i<=Lens.ngxs; i++)
		 {
			 gxs = i*Lens.dgxs; gys = j*Lens.dgys;
			 g2s = gxs*gxs + gys*gys;
			if(g2s < Lens.gmax2s)
			{
				Qst.x[nQs] = gxs;
				Qst.y[nQs] = gys;
				sumwia += Qst.w[nQs] = exp(-alpha*g2s);
				nQs++;
			}
		}
	 }
	/**********************************************************************/
	Qs.x = new double [nQs];
	Qs.y = new double [nQs];
	Qs.w = new double [nQs];

	for(i=0; i<nQs; i++)
	{
		Qs.x[i] = Qst.x[i];
		Qs.y[i] = Qst.y[i];
		Qs.w[i] = Qst.w[i]/sumwia;
	}
	/**********************************************************************/
	delete [] Qst.x; Qst.x = 0;
	delete [] Qst.y; Qst.y = 0;
	delete [] Qst.w; Qst.w = 0;
}

void cMT_MicroscopeEffects_CPU::SetInputData(sGP &GP_i, sLens &Lens_i, fftw_plan &PlanForward_i, fftw_plan &PlanBackward_i, fftw_complex *&Psit_i)
{
	freeMemory();
	IdCall++;

	GP = GP_i;
	Lens = Lens_i;
	PlanForward = PlanForward_i;
	PlanBackward = PlanBackward_i;
	Psit = Psit_i;
	/*********************Temporal quadrature**********************/
	//ReadTemporalQuadrature(Lens, Qt);
	Qt.x = new double[Lens.nsf]; 
	Qt.w = new double [Lens.nsf];
	cQuadrature Quad;
	Quad.ReadQuadrature(8, Lens.nsf, Qt);		// 8: int_-infty^infty f(x) Exp[-x^2] dx
	for(int i=0; i<Lens.nsf; i++)
	{
		Qt.w[i] /= cPii2;
	}
	/*********************Spatial quadrature**********************/
	ReadSpatialQuadrature(Lens, nQs, Qs);
}

// Inclusion of the microscope effect: TypCal 1: PCLIMWPO, 2: PCTCCTEM
void cMT_MicroscopeEffects_CPU::ApplyMEffects(int MEffect, int STEffect, fftw_complex *&fPsi, double *&M2Psi)
{
	switch(MEffect)
	{
		case 1:
			PCLIMWPOTEM(STEffect, fPsi, M2Psi);	
			break;
		case 2:
			PCTCCTEM(STEffect, fPsi, M2Psi);	
			break;
	}
}