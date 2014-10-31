#include <cmath>
#include <memory>
#include "hConstTypes.h"
#include "hgeneralCPU.h"
#include "hDetectorCPU.h"

void cDetectorCPU::freeMemory(){
	f_sGP_Init(GP);

	nDet = 0;
	f_sDetCir_Free(DetCir);

	delete [] Tot; Tot = 0;
	delete [] Coh; Coh = 0;
}

cDetectorCPU::~cDetectorCPU(){
	freeMemory();
}

cDetectorCPU::cDetectorCPU(){
	f_sGP_Init(GP);

	nDet = 0;
	f_sDetCir_Init(DetCir);

	Tot = 0;
	Coh = 0;
}

/******************************************************************************/
/******************************************************************************/

void cDetectorCPU::SetInputData(sGP &GP_i, int nDeti, sDetCir &DetCiri){
	freeMemory();

	GP = GP_i;
	nDet = nDeti;
	f_sDetCir_Malloc(nDet, DetCir);
	memcpy(DetCir.g2min, DetCiri.g2min, nDet*cSizeofRD);
	memcpy(DetCir.g2max, DetCiri.g2max, nDet*cSizeofRD);

	Tot = new double[nDet];
	Coh = new double[nDet];
}

void cDetectorCPU::getDetectorIntensity(double *&Toth, double *&Cohh, int ixys, sDetInt *DetInth){
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