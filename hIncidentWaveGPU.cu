#include "hConstTypes.h"
#include "hgeneralCPU.h"
#include "hgeneralGPU.h"
#include "hIncidentWaveGPU.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "cufft.h"

// Incident wave in Fourier Space
__global__ void IncidentWaveFS(sGP GP, sLens Lens, double x, double y, double2 *Psig){
	int j = threadIdx.x + blockIdx.x*blockDim.x;
	int i = threadIdx.y + blockIdx.y*blockDim.y;

	if ((i < GP.nx)&&(j < GP.ny)){	
		int k = i*GP.ny+j;
		double gx = ((i<GP.nxh)?i:(i-GP.nx))*GP.dgx;
		double gy = ((j<GP.nyh)?j:(j-GP.ny))*GP.dgy;
		double g2 = gx*gx + gy*gy;
		if ((Lens.gmin2 <= g2)&&(g2 < Lens.gmax2)){
			double chi = x*gx + y*gy + g2*(Lens.cCs5*g2*g2+Lens.cCs3*g2+Lens.cf);
			if ((Lens.m!=0)||(Lens.cmfa2!=0)||(Lens.cmfa3!=0)){
				double g = sqrt(g2);
				double phi = atan(gy/gx);
				chi += Lens.m*phi + Lens.cmfa2*g2*sin(2*(phi-Lens.afa2)) + Lens.cmfa3*g*g2*sin(3*(phi-Lens.afa3));				
			}
			sincos(chi, &gy , &gx);		
			Psig[k].x = gx; 
			Psig[k].y = gy;	
		}
		else{
 			Psig[k].x = 0.0; 
			Psig[k].y = 0.0; 
		}
	}
}

void cIncidentWaveGPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	fsGP_Init(GP);
	fsBT_Init(BT);
	fsLens_Init(Lens);

	cudaFreen(Sd);
	PlanPsi = 0;
}

cIncidentWaveGPU::cIncidentWaveGPU(){
	fsGP_Init(GP);
	fsBT_Init(BT);
	fsLens_Init(Lens);

	Sd = 0;
	PlanPsi = 0;
}

cIncidentWaveGPU::~cIncidentWaveGPU(){
	freeMemory();
}

void cIncidentWaveGPU::SetInputData(sBT &BT_i, sGP &GP_i, sLens &Lens_i, cufftHandle &PlanPsi_i){
	freeMemory();

	BT = BT_i;
	GP = GP_i;
	Lens = Lens_i;
	PlanPsi = PlanPsi_i;
	cudaMalloc((void**)&Sd, 2*thrnxy*cSizeofRD);
}

void cIncidentWaveGPU::Psi0(double2 *&Psig){
	SetValueVectorC(BT.Bnxy, BT.Tnxy, GP.nxy, 1.0, 0.0, Psig);
}

void cIncidentWaveGPU::Psi0(double x, double y, double2 *&Psig){
	IncidentWaveFS<<<BT.Bnxny, BT.Tnxny>>>(GP, Lens, c2Pi*x, c2Pi*y, Psig);
	cufftExecZ2Z(PlanPsi, Psig, Psig, CUFFT_INVERSE);
	double Totalsum = SumAc2(BT.Bnxy, BT.Tnxy, GP.nxy, Psig, Sd);
	ScaleVectorC(BT.Bnxy, BT.Tnxy, GP.nxy, sqrt(GP.nxy/Totalsum), Psig);
}