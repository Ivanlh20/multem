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

#include <cstring>
#include "math.h"

#include "hConstTypes.h"
#include "hMT_InMulSli_CPU.h"
#include "hMT_InMulSli_CPU.h"
#include "hMT_Detector_GPU.h"
#include "hMT_STEM_GPU.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

void cMT_STEM_GPU::freeMemory()
{
	if(IdCall==0) return;

	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	delete MT_Detector_GPU; MT_Detector_GPU = 0;

	line = 0;
	ns = 0;	
	x1u = 0;
	y1u = 0;		
	x2u = 0;	
	y2u = 0;	

	f_sDetCir_Free(DetCir);
	for(int iThk = 0; iThk<nThk; iThk++)
	{
		for(int iDet=0; iDet<nDet; iDet++)
		{
			delete [] ImSTEM[iThk].DetInt[iDet].Coh; ImSTEM[iThk].DetInt[iDet].Coh = 0;
			delete [] ImSTEM[iThk].DetInt[iDet].Tot; ImSTEM[iThk].DetInt[iDet].Tot = 0;
		}
		delete [] ImSTEM[iThk].DetInt; ImSTEM[iThk].DetInt = 0;
	}
	delete [] ImSTEM; ImSTEM = 0;
	nDet = 0;
	nThk = 0;

	nxs = 0;
	nys = 0;
	delete [] xs; xs = 0;
	delete [] ys; ys = 0;

	nst = 0;
	delete [] xst; xst = 0;
	delete [] yst; yst = 0;
}

cMT_STEM_GPU::cMT_STEM_GPU()
{
	IdCall = 0;
	MT_Detector_GPU = 0;

	line = 0;
	ns = 0;
	x1u = 0;
	y1u = 0;		
	x2u = 0;	
	y2u = 0;	

	f_sDetCir_Init(DetCir);
	ImSTEM = 0;
	nDet = 0;

	nThk = 0;

	nxs = 0;
	nys = 0;
	xs = 0;
	ys = 0;

	nst = 0;
	xst = 0;
	yst = 0;
}

cMT_STEM_GPU::~cMT_STEM_GPU()
{
	freeMemory();
	IdCall = 0;
}

void cMT_STEM_GPU::InitImSTEM()
{
	int iThk, iDet, ist;
	for(iThk = 0; iThk<nThk; iThk++)
		for(iDet=0; iDet<nDet; iDet++)
			for(ist=0; ist<nst; ist++)
			{
				ImSTEM[iThk].DetInt[iDet].Coh[ist] = 0.0;
				ImSTEM[iThk].DetInt[iDet].Tot[ist] = 0.0;
			}
}

void cMT_STEM_GPU::SetInputData(cMT_InMulSli_CPU *MT_InMulSli_CPU_i, cMT_MGP_CPU *MT_MGP_CPU_i, int nThk_i)
{
	freeMemory();
	IdCall++;

	line = MT_InMulSli_CPU_i->STEM_line;
	ns = MT_InMulSli_CPU_i->STEM_ns;
	x1u = MT_InMulSli_CPU_i->STEM_x1u;	
	y1u = MT_InMulSli_CPU_i->STEM_y1u;
	x2u = MT_InMulSli_CPU_i->STEM_x2u;
	y2u = MT_InMulSli_CPU_i->STEM_y2u;
	f_BuildGrid(line, ns, x1u, y1u, x2u, y2u, nxs, nys, xs, ys);

	nDet = MT_InMulSli_CPU_i->STEM_nDet;
	double lambda = f_getLambda(MT_MGP_CPU_i->E0);
	f_sDetCir_Malloc(nDet, DetCir);
	for(int iDet=0; iDet<nDet; iDet++)
	{
		DetCir.g2min[iDet] = pow(MT_InMulSli_CPU_i->STEM_DetCir[iDet].InnerAng/lambda, 2);
		DetCir.g2max[iDet] = pow(MT_InMulSli_CPU_i->STEM_DetCir[iDet].OuterAng/lambda, 2);
	}

	sGP GP;
	f_sGP_SetInputData(MT_MGP_CPU_i, GP);
	MT_Detector_GPU = new cMT_Detector_GPU;
	MT_Detector_GPU->SetInputData(GP, nDet, DetCir);

	nst = (line==1)?ns:nxs*nys;
	int ils, ixs, iys, ixys;
	xst = new double[nst];
	yst = new double[nst];
	if(line==1)
	{
		for(ils=0; ils<ns; ils++)
		{
			xst[ils] = xs[ils];
			yst[ils] = ys[ils];
		}
	}
	else
	{
		for(ixs=0; ixs<nxs; ixs++)
			for(iys=0; iys<nys; iys++)
			{
				ixys = ixs*nys + iys;
				xst[ixys] = xs[ixs];
				yst[ixys] = ys[iys];
			}
	}

	nThk = nThk_i;
	//nThk = 81;
	int iThk, iDet, ist;
	ImSTEM = new sImSTEM[nThk];
	for(iThk = 0; iThk<nThk; iThk++)
	{
		ImSTEM[iThk].DetInt = new sDetInt[nDet];
		for(iDet=0; iDet<nDet; iDet++)
		{
			ImSTEM[iThk].DetInt[iDet].Coh = new double[nst];
			ImSTEM[iThk].DetInt[iDet].Tot = new double[nst];
			for(ist=0; ist<nst; ist++)
			{
				ImSTEM[iThk].DetInt[iDet].Coh[ist] = 0.0;
				ImSTEM[iThk].DetInt[iDet].Tot[ist] = 0.0;
			}
		}
	}
}