/**
 *  This file is part of MULTEM.
 *  Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
 *
 *  MULTEM is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MULTEM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MULTEM.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstring>
#include "math.h"
#include "hConstTypes.h"
#include "hMT_General_CPU.h"
#include "hMT_General_GPU.h"
#include "hMT_Potential_GPU.h"
#include "hMT_IncidentWave_GPU.h"
#include "hMT_Detector_CPU.h"
#include "hMT_Detector_GPU.h"
#include "hMT_STEM_GPU.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

void cSTEM_GPU::freeMemory()
{
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	line = 0;
	FastCal = false;
	ns = 0;	
	x1u = 0;
	y1u = 0;		
	x2u = 0;	
	y2u = 0;	

	f_sDetCir_Free(DetCir);
	for (int iThk = 0; iThk<nThk; iThk++){
		for (int iDet=0; iDet<nDet; iDet++){
			delete [] ImSTEM[iThk].DetInt[iDet].Coh; ImSTEM[iThk].DetInt[iDet].Coh = 0;
			delete [] ImSTEM[iThk].DetInt[iDet].Tot; ImSTEM[iThk].DetInt[iDet].Tot = 0;
		}
		delete [] ImSTEM[iThk].DetInt; ImSTEM[iThk].DetInt = 0;
	}
	delete [] ImSTEM; ImSTEM = 0;
	nDet = 0;

	nThk = 0;
	delete [] Thk; Thk = 0;

	nxs = 0;
	nys = 0;
	delete [] xs; xs = 0;
	delete [] ys; ys = 0;

	nst = 0;
	delete [] xst; xst = 0;
	delete [] yst; yst = 0;

	if(nSliceM>0){
		switch(SliceMTyp)
		{
			case 1:
				for(int iSliceM=0; iSliceM<nSliceM; iSliceM++)
					cudaFreen(Trans[iSliceM]);
				delete [] Trans; Trans = 0;
				break;
			case 2:
				for(int iSliceM=0; iSliceM<nSliceM; iSliceM++)
					cudaFreen(VpD[iSliceM]);
				delete [] VpD; VpD = 0;
				break;
			case 3:
				for(int iSliceM=0; iSliceM<nSliceM; iSliceM++)
					cudaFreen(VpF[iSliceM]);
				delete [] VpF; VpF = 0;
				break;
		}
	}

	nSliceM = 0;
	SliceMTyp = 0;
}

cSTEM_GPU::cSTEM_GPU()
{
	line = 0;
	FastCal = false;
	ns = 0;
	x1u = 0;
	y1u = 0;		
	x2u = 0;	
	y2u = 0;	

	f_sDetCir_Init(DetCir);
	ImSTEM = 0;
	nDet = 0;

	nThk = 0;
	Thk = 0;

	nxs = 0;
	nys = 0;
	xs = 0;
	ys = 0;

	nst = 0;
	xst = 0;
	yst = 0;

	Trans = 0;
	VpD = 0;
	VpF = 0;

	nSliceM = 0;
	SliceMTyp = 0;
}

cSTEM_GPU::~cSTEM_GPU(){
	freeMemory();
}

void cSTEM_GPU::SetInputData(sInMSTEM &InMSTEM, cMT_Specimen_CPU *MT_Specimen_CPU)
{
	freeMemory();

	MGP = MT_Specimen_CPU->MGP;
	f_sGP_Cal(MGP.nx, MGP.ny, MGP.lx, MGP.ly, MGP.dz, MGP.PBC_xy, MGP.BWL, GP);

	line = InMSTEM.STEM_line;
	FastCal = InMSTEM.STEM_FastCal;
	ns = InMSTEM.STEM_ns;
	x1u = InMSTEM.STEM_x1u;	
	y1u = InMSTEM.STEM_y1u;
	x2u = InMSTEM.STEM_x2u;
	y2u = InMSTEM.STEM_y2u;
	f_BuildGrid(line, ns, x1u, y1u, x2u, y2u, nxs, nys, xs, ys);

	nThk = MGP.nThk;
	if(nThk>0){
		Thk = new double[nThk];
		memcpy(Thk, InMSTEM.Thickness, nThk*cSizeofRD);
	}

	nDet = InMSTEM.STEM_nDet;
	double lambda = f_getLambda(InMSTEM.E0);
	f_sDetCir_Malloc(nDet, DetCir);
	for (int iDet=0; iDet<nDet; iDet++){
		DetCir.g2min[iDet] = pow(InMSTEM.STEM_DetCir[iDet].InnerAng/lambda, 2);
		DetCir.g2max[iDet] = pow(InMSTEM.STEM_DetCir[iDet].OuterAng/lambda, 2);
	}

	nst = (line==1)?ns:nxs*nys;
	int ils, ixs, iys, ixys;
	xst = new double[nst];
	yst = new double[nst];
	if(line==1){
		for (ils=0; ils<ns; ils++){
			xst[ils] = xs[ils];
			yst[ils] = ys[ils];
		}
	}else{
		for (ixs=0; ixs<nxs; ixs++)
			for (iys=0; iys<nys; iys++){
				ixys = ixs*nys + iys;
				xst[ixys] = xs[ixs];
				yst[ixys] = ys[iys];
			}
	}

	int iThk, iDet, ist;
	ImSTEM = new sImSTEM[nThk];
	for (iThk = 0; iThk<nThk; iThk++){
		ImSTEM[iThk].DetInt = new sDetInt[nDet];
		for (iDet=0; iDet<nDet; iDet++){
			ImSTEM[iThk].DetInt[iDet].Coh = new double[nst];
			ImSTEM[iThk].DetInt[iDet].Tot = new double[nst];
			for (ist=0; ist<nst; ist++){
				ImSTEM[iThk].DetInt[iDet].Coh[ist] = 0.0;
				ImSTEM[iThk].DetInt[iDet].Tot[ist] = 0.0;
			}
		}
	}

	/*****************************************************************/
	size_t SizeFreeM, SizeTotM;
	cudaMemGetInfo(&SizeFreeM, &SizeTotM);
	size_t SizeTrans = GP.nxy*(SliceMTyp==1)?cSizeofCD:(SliceMTyp==2)?cSizeofRD:cSizeofRF;
	int nSliceMt = (SizeFreeM-10*cMb)/SizeTrans;

	if((FastCal)&&(nSliceMt>0)&&(MGP.SimType==1)&&(MGP.ApproxModel<=2)){
		nSliceM = nSliceMt;
		switch(SliceMTyp)
		{
			case 1:
				Trans = new double2*[nSliceM];
				for(int iSliceM=0; iSliceM<nSliceM; iSliceM++)
					cudaMalloc((void**)&Trans[iSliceM], SizeTrans);
				break;
			case 2:
				VpD = new double*[nSliceM];
				for(int iSliceM=0; iSliceM<nSliceM; iSliceM++)
					cudaMalloc((void**)&VpD[iSliceM], SizeTrans);
				break;
			case 3:
				VpF = new float*[nSliceM];
				for(int iSliceM=0; iSliceM<nSliceM; iSliceM++)
					cudaMalloc((void**)&VpF[iSliceM], SizeTrans);
				break;
		}
	}
}

void cSTEM_GPU::InitImSTEM()
{
	int iThk, iDet, ist;
	for (iThk = 0; iThk<nThk; iThk++)
		for (iDet=0; iDet<nDet; iDet++)
			for (ist=0; ist<nst; ist++){
				ImSTEM[iThk].DetInt[iDet].Coh[ist] = 0.0;
				ImSTEM[iThk].DetInt[iDet].Tot[ist] = 0.0;
			}
}