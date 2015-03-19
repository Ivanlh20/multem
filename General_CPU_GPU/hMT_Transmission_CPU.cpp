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
#include "hGeneral_CPU.h"
#include "hMT_General_CPU.h"
#include "hMT_MGP_CPU.h"
#include "hMT_AtomTypes_CPU.h"
#include "hMT_Transmission_CPU.h"

// Calculated transmission function
template <class Type>
void t_Transmission(sGP &GP, int ApproxModel, double f, const Type * __restrict V0_i, fftw_complex * __restrict Trans_o)
{
	int ix, iy, ixy;
	double theta;

	for(iy=0; iy<GP.ny; iy++)
	{
		for(ix=0; ix<GP.nx; ix++)
		{
			ixy = iy*GP.nx+ix;
			theta = f*V0_i[ixy];
			Trans_o[ixy][0] = (ApproxModel!=4)?cos(theta):1.0;
			Trans_o[ixy][1] = (ApproxModel!=4)?sin(theta):theta;
		}
	}

}

void cMT_Transmission_CPU::freeMemory()
{
	if(IdCall==0) return;

	fPot = 0.0;

	fftw_free(Trans0); Trans0 = 0;

	if(nSliceMem0>0)
	{
		if(SliceMemTyp==1)
		{
			for(int iSliceMem=0; iSliceMem<nSliceMem0; iSliceMem++)
			{
				fftw_free(Trans[iSliceMem]); 
				Trans[iSliceMem] = 0;
			}
			delete [] Trans; Trans = 0;
		}
		else
		{
			for(int iSliceMem=0; iSliceMem<nSliceMem0; iSliceMem++)
			{
				delete [] Vpe[iSliceMem];
				Vpe[iSliceMem] = 0;
			}
			delete [] Vpe; Vpe = 0;
		}
	}

	SliceMemTyp = 0;
	nSliceMem = 0;
	nSliceMem0 = 0;
}

cMT_Transmission_CPU::cMT_Transmission_CPU()
{
	IdCall = 0;

	fPot = 0.0;

	SliceMemTyp = 0;
	nSliceMem = 0;
	nSliceMem0 = 0;

	Trans0 = 0;
	Trans = 0;
	Vpe = 0;
}

cMT_Transmission_CPU::~cMT_Transmission_CPU()
{
	freeMemory();
	IdCall = 0;
}

void cMT_Transmission_CPU::Cal_Trans_or_Vpe()
{
	for(int iSliceMem=0; iSliceMem<nSliceMem; iSliceMem++)
	{
		ProjectedPotential(iSliceMem);
		if(SliceMemTyp==1)
		{
			t_Transmission<double>(GP, MT_MGP_CPU->ApproxModel, fPot, V0, Trans[iSliceMem]);
			f_BandwidthLimit2D_CPU(PlanForward, PlanBackward, GP, Trans[iSliceMem]);
		}
		else
		{
			f_Double_2_Float_CPU(GP, V0, Vpe[iSliceMem]);
		}
	}
}

void cMT_Transmission_CPU::SetInputData(cMT_MGP_CPU *MT_MGP_CPU_io, fftw_plan &PlanForward_i, fftw_plan &PlanBackward_i, int nAtomsM_i, double *AtomsM_i)
{
	freeMemory();
	IdCall++;

	cMT_Potential_CPU::SetInputData(MT_MGP_CPU_io, nAtomsM_i, AtomsM_i);
	PlanForward = PlanForward_i;
	PlanBackward = PlanBackward_i;

	Trans0 = (fftw_complex*)fftw_malloc(GP.nxy*cSizeofCD);

	fPot = f_getfPot(MT_MGP_CPU->E0, MT_MGP_CPU->theta);

	if((MT_MGP_CPU->FastCal==1)||(MT_MGP_CPU->ApproxModel>2))
	{
		return;
	}

	int nSliceSigma = (MT_MGP_CPU->DimFP%10==0)?0:(int)ceil(6*sigma_max/MT_MGP_CPU->dz);
	int nSliceMax = nSlice + nSliceSigma;
	if(MT_MGP_CPU->MulOrder==2)
	{
		nSliceMax++;
	}

	size_t SizeFreeMem = 6*cGb; // find solutions
	SizeFreeMem = SizeFreeMem-10*cMb;
	int nSliceMemMax = 0;

	if(SizeFreeMem/(GP.nxy*cSizeofCD)>=nSliceMax)
	{
		SliceMemTyp = 1;
		nSliceMemMax = SizeFreeMem/(GP.nxy*cSizeofCD);
	}
	else
	{
		SliceMemTyp = 2;
		nSliceMemMax = SizeFreeMem/(GP.nxy*cSizeofRF);
	}

	if((nSliceMemMax>0)&&(MT_MGP_CPU->ApproxModel<=2))
	{
		nSliceMem0 = MIN(nSliceMemMax, nSliceMax);
		if(SliceMemTyp==1)
		{
			Trans = new fftw_complex*[nSliceMem0];
			for(int iSliceMem=0; iSliceMem<nSliceMem0; iSliceMem++)
			{
				Trans[iSliceMem] = (fftw_complex*)fftw_malloc(GP.nxy*cSizeofCD);
			}
		}
		else
		{
			Vpe = new float*[nSliceMem0];
			for(int iSliceMem=0; iSliceMem<nSliceMem0; iSliceMem++)
			{
				Vpe[iSliceMem] = new float[GP.nxy*cSizeofRF];
			}
		}
	}

	nSliceMem = MIN(nSliceMem0, nSlice);
	if((MT_MGP_CPU->MulOrder==2)&&(nSliceMem0>nSlice))
	{
		nSliceMem++;
	}
}

void cMT_Transmission_CPU::MoveAtoms(int iConf)
{
	cMT_Potential_CPU::MoveAtoms(iConf);
	int nSliceMem = MIN(nSliceMem0, nSlice);
	if((MT_MGP_CPU->MulOrder==2)&&(nSliceMem0>nSlice))
	{
		nSliceMem++;
	}

	if(nSliceMem0>0)
	{
		Cal_Trans_or_Vpe();
	}

	if(MT_MGP_CPU->ApproxModel>2)
	{
		ProjectedPotential(0);
		t_Transmission<double>(GP, MT_MGP_CPU->ApproxModel, fPot, V0, Trans0);
		f_BandwidthLimit2D_CPU(PlanForward, PlanBackward, GP, Trans0);	
	}
}

fftw_complex* cMT_Transmission_CPU::getTrans(int iSlice, int typ)
{
	if(MT_MGP_CPU->ApproxModel>2)
	{
		return Trans0;
	}

	fftw_complex *Trans_o = Trans0;
	if(iSlice<nSliceMem)
	{
		if(SliceMemTyp==1)
		{
			Trans_o = Trans[iSlice];
		}
		else
		{
			t_Transmission<float>(GP, MT_MGP_CPU->ApproxModel, fPot, Vpe[iSlice], Trans_o);
			f_BandwidthLimit2D_CPU(PlanForward, PlanBackward, GP, Trans_o);	
		}
	}
	else
	{
		ProjectedPotential(iSlice, typ);
		t_Transmission<double>(GP, MT_MGP_CPU->ApproxModel, fPot, V0, Trans_o);
		f_BandwidthLimit2D_CPU(PlanForward, PlanBackward, GP, Trans_o);	
	}
	return Trans_o;
}

void cMT_Transmission_CPU::Transmit(int iSlice, fftw_complex *&Psi_io)
{
	fftw_complex *Trans = (MT_MGP_CPU->ApproxModel>2)?Trans0:getTrans(iSlice);
	f_Hadamard_Product_CPU(GP, Trans, Psi_io);
}