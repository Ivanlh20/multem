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

#include "hConstTypes.h"
#include "hQuadrature.h"
#include "hMT_General_CPU.h"
#include "hMT_General_GPU.h"
#include "hMT_MGP_CPU.h"
#include "hMT_AtomTypes_GPU.h"
#include "hMT_Transmission_GPU.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

// Calculated transmission function
template <class Type>
__global__ void k_Transmission(sGP GP, int ApproxModel, double f, const Type * __restrict V0_i, double2 * __restrict Trans_o)
{
	int iy = threadIdx.x + blockIdx.x*blockDim.x;
	int ix = threadIdx.y + blockIdx.y*blockDim.y;

	if((ix < GP.nx)&&(iy < GP.ny))
	{
		int ixy = ix*GP.ny+iy;
		double V0 = V0_i[ixy];
		double theta = f*V0, x = 1.0, y = theta;
		if(ApproxModel!=4)
		{
			sincos(theta, &y , &x);
		}
		Trans_o[ixy].x = x;
		Trans_o[ixy].y = y;
	}
}

void cMT_Transmission_GPU::freeMemory()
{
	if(IdCall==0)
	{
		return;
	}

	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	cSynCPU = ccSynCPU;
	fPot = 0.0;

	cudaFreen(Trans0);

	if(nSliceMem0>0)
		if(SliceMemTyp==1)
		{
			for(int iSliceMem=0; iSliceMem<nSliceMem0; iSliceMem++)
			{
				cudaFreen(Trans[iSliceMem]);
			}
			delete [] Trans; Trans = 0;
		}
		else
		{
			for(int iSliceMem=0; iSliceMem<nSliceMem0; iSliceMem++)
			{
				cudaFreen(Vpe[iSliceMem]);
			}
			delete [] Vpe; Vpe = 0;
		}

	SliceMemTyp = 0;
	nSliceMem = 0;
	nSliceMem0 = 0;

	PlanTrans = 0;
}

void cMT_Transmission_GPU::freeMemoryReset()
{
	freeMemory();
	IdCall = 0;
	cudaDeviceReset();
}

cMT_Transmission_GPU::cMT_Transmission_GPU()
{
	IdCall = 0;

	cSynCPU = ccSynCPU;
	fPot = 0.0;

	SliceMemTyp = 0;
	nSliceMem = 0;
	nSliceMem0 = 0;

	Trans0 = 0;
	Trans = 0;
	Vpe = 0;

	PlanTrans = 0;
}

cMT_Transmission_GPU::~cMT_Transmission_GPU()
{
	freeMemory();
	IdCall = 0;
}

void cMT_Transmission_GPU::Cal_Trans_or_Vpe()
{
	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);

	for(int iSliceMem=0; iSliceMem<nSliceMem; iSliceMem++)
	{
		ProjectedPotential(iSliceMem);
		if(SliceMemTyp==1)
		{
			k_Transmission<double><<<Bnxny, Tnxny>>>(GP, MT_MGP_CPU->ApproxModel, fPot, V0, Trans[iSliceMem]);
			f_BandwidthLimit2D_GPU(PlanTrans, GP, Trans[iSliceMem]);
		}
		else
		{
			f_Double_2_Float_GPU(GP, V0, Vpe[iSliceMem]);
		}
	}
	cudaDeviceSynchronize();
}

void cMT_Transmission_GPU::SetInputData(cMT_MGP_CPU *MT_MGP_CPU_io, cufftHandle &PlanTrans_i, int nAtomsM_i, double *AtomsM_i)
{
	freeMemory();
	IdCall++;

	cMT_Potential_GPU::SetInputData(MT_MGP_CPU_io, nAtomsM_i, AtomsM_i);
	PlanTrans = PlanTrans_i;

	cudaMalloc((void**)&Trans0, GP.nxy*cSizeofCD);

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

	size_t SizeFreeMem, SizeTotMem;
	cudaMemGetInfo(&SizeFreeMem, &SizeTotMem);
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
			Trans = new double2*[nSliceMem0];
			for(int iSliceMem=0; iSliceMem<nSliceMem0; iSliceMem++)
			{
				cudaMalloc((void**)&Trans[iSliceMem], GP.nxy*cSizeofCD);
			}
		}
		else
		{
			Vpe = new float*[nSliceMem0];
			for(int iSliceMem=0; iSliceMem<nSliceMem0; iSliceMem++)
			{
				cudaMalloc((void**)&Vpe[iSliceMem], GP.nxy*cSizeofRF);
			}
		}
	}

	nSliceMem = MIN(nSliceMem0, nSlice);
	if((MT_MGP_CPU->MulOrder==2)&&(nSliceMem0>nSlice))
	{
		nSliceMem++;
	}
}

void cMT_Transmission_GPU::MoveAtoms(int iConf)
{
	cMT_Potential_GPU::MoveAtoms(iConf);
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
		dim3 Bnxny, Tnxny;
		f_get_BTnxny(GP, Bnxny, Tnxny);
		ProjectedPotential(0);
		k_Transmission<double><<<Bnxny, Tnxny>>>(GP, MT_MGP_CPU->ApproxModel, fPot, V0, Trans0);
		f_BandwidthLimit2D_GPU(PlanTrans, GP, Trans0);	
	}
}

double2* cMT_Transmission_GPU::getTrans(int iSlice, int typ)
{
	if(MT_MGP_CPU->ApproxModel>2) return Trans0;

	dim3 Bnxny, Tnxny;
	f_get_BTnxny(GP, Bnxny, Tnxny);

	double2 *Trans_o = Trans0;
	if(iSlice<nSliceMem)
	{
		if(SliceMemTyp==1)
		{
			Trans_o = Trans[iSlice];
		}
		else
		{
			k_Transmission<float><<<Bnxny, Tnxny>>>(GP, MT_MGP_CPU->ApproxModel, fPot, Vpe[iSlice], Trans_o);
			f_BandwidthLimit2D_GPU(PlanTrans, GP, Trans_o);	
		}
	}
	else
	{
		ProjectedPotential(iSlice, typ);
		k_Transmission<double><<<Bnxny, Tnxny>>>(GP, MT_MGP_CPU->ApproxModel, fPot, V0, Trans_o);
		f_BandwidthLimit2D_GPU(PlanTrans, GP, Trans_o);	
	}
	return Trans_o;
}

void cMT_Transmission_GPU::Transmit(int iSlice, double2 *&Psi_io)
{
	double2 *Trans = (MT_MGP_CPU->ApproxModel>2)?Trans0:getTrans(iSlice);
	f_Hadamard_Product_GPU(GP, Trans, Psi_io);
}