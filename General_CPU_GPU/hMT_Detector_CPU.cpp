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
#include <cstring>
#include "hConstTypes.h"
#include "hMT_General_CPU.h"
#include "hMT_Detector_CPU.h"

void cMT_Detector_CPU::freeMemory(){
	f_sGP_Init(GP);

	nDet = 0;
	f_sDetCir_Free(DetCir);

	delete [] Tot; Tot = 0;
	delete [] Coh; Coh = 0;
}

cMT_Detector_CPU::~cMT_Detector_CPU(){
	freeMemory();
}

cMT_Detector_CPU::cMT_Detector_CPU(){
	f_sGP_Init(GP);

	nDet = 0;
	f_sDetCir_Init(DetCir);

	Tot = 0;
	Coh = 0;
}

/*****************************************************************************/
/*****************************************************************************/

void cMT_Detector_CPU::SetInputData(sGP &GP_i, int nDeti, sDetCir &DetCiri){
	freeMemory();

	GP = GP_i;
	nDet = nDeti;
	f_sDetCir_Malloc(nDet, DetCir);
	memcpy(DetCir.g2min, DetCiri.g2min, nDet*cSizeofRD);
	memcpy(DetCir.g2max, DetCiri.g2max, nDet*cSizeofRD);

	Tot = new double[nDet];
	Coh = new double[nDet];
}

void cMT_Detector_CPU::getDetectorIntensity(double *&Toth, double *&Cohh, int ixys, sDetInt *DetInth){
	int iDet, ix, iy, ixy;
	double gx, gy, g2;

	for(iDet = 0; iDet<nDet; iDet++)
		Tot[iDet]= Coh[iDet] = 0.0;

	for(ix = 0; ix<GP.nx; ix++)
	{
		gx = IsFS(ix,GP.nxh)*GP.dgx;
		for(iy = 0; iy<GP.ny; iy++)
		{	
			ixy = ix*GP.ny+iy;
			gy = IsFS(iy,GP.nyh)*GP.dgy;
			g2 = gx*gx + gy*gy;
			for(iDet = 0; iDet<nDet; iDet++)
				if((DetCir.g2min[iDet]<=g2)&&(g2<=DetCir.g2max[iDet]))
				{
					Tot[iDet] += Toth[ixy];
					Coh[iDet] += Cohh[ixy];
				}
		}
	}

	for(iDet = 0; iDet<nDet; iDet++){
		DetInth[iDet].Tot[ixys] = Tot[iDet];
		DetInth[iDet].Coh[ixys] = Coh[iDet];
	}
}