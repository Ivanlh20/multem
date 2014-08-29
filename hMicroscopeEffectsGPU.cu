#include "math.h"
#include <memory>
#include "hConstTypes.h"
#include "hgeneralGPU.h"
#include "hQuadrature.h"
#include "hMicroscopeEffectsGPU.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "cufft.h"
/*********************************************************************/
// Apply Coherent transfer function
__global__ void ApplyCTF(sGP GP, sLens Lens, double gxu, double gyu, double2 *fPsig_i, double2 *fPsig_o){
	int j = threadIdx.x + blockIdx.x*blockDim.x;
	int i = threadIdx.y + blockIdx.y*blockDim.y;

	if ((i < GP.nx)&&(j < GP.ny)){
		int k = i*GP.ny+j;
		double gx = ((i<GP.nxh)?i:(i-GP.nx))*GP.dgx;
		double gy = ((j<GP.nyh)?j:(j-GP.ny))*GP.dgy;
		double g2 = gx*gx + gy*gy;
		if ((Lens.gmin2 <= g2)&&(g2 <= Lens.gmax2)){
			gx -= gxu; gy -= gyu;
			g2 = gx*gx + gy*gy;
			double chi = g2*(Lens.cCs5*g2*g2+Lens.cCs3*g2+Lens.cf);
			if ((Lens.cmfa2!=0)||(Lens.cmfa3!=0)){
				double g = sqrt(g2);
				double phi = atan(gy/gx);
				chi += Lens.cmfa2*g2*sin(2*(phi-Lens.afa2)) + Lens.cmfa3*g*g2*sin(3*(phi-Lens.afa3));				
			}
			double x = fPsig_i[k].x, y = fPsig_i[k].y; 
			sincos(chi, &gy, &gx);		
			fPsig_o[k].x = (gx*x-gy*y); 
			fPsig_o[k].y = (gx*y+gy*x);
		}
		else{
 			fPsig_o[k].x = 0.0; 
			fPsig_o[k].y = 0.0; 
		}
	}
}

// Partially coherent transfer function, linear image model and weak phase object
__global__ void ApplyPCTF(sGP GP, sLens Lens, double2 *fPsig_i, double2 *fPsig_o){
	int j = threadIdx.x + blockIdx.x*blockDim.x;
	int i = threadIdx.y + blockIdx.y*blockDim.y;

	if ((i < GP.nx)&&(j < GP.ny)){
		int k = i*GP.ny+j;
		double gx = ((i<GP.nxh)?i:(i-GP.nx))*GP.dgx;
		double gy = ((j<GP.nyh)?j:(j-GP.ny))*GP.dgy;
		double g2 = gx*gx + gy*gy;

		if ((Lens.gmin2 <= g2)&&(g2 <= Lens.gmax2)){			
			double chi = g2*(Lens.cCs3*g2+Lens.cf);
			double c = cPi*Lens.beta*Lens.sf;
			double u = 1.0 + c*c*g2;

			c = cPi*Lens.sf*Lens.lambda*g2;
			gx = 0.25*c*c;
			c = cPi*Lens.beta*(Lens.Cs3*Lens.lambda2*g2-Lens.f);
			gy = c*c*g2;
			g2 = exp(-(gx+gy)/u);

			double x = fPsig_i[k].x, y = fPsig_i[k].y; 
			sincos(chi, &gy, &gx);	 	

			fPsig_o[k].x = g2*(gx*x-gy*y); 
			fPsig_o[k].y = g2*(gx*y+gy*x);
		}else{
 			fPsig_o[k].x = 0.0; 
			fPsig_o[k].y = 0.0; 
		}
	}
}

/****************************************************************************/
void cMicroscopeEffectsGPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	cSynCPU = ccSynCPU;

	delete [] Qt.x;  Qt.x = 0;
	delete [] Qt.w;  Qt.w = 0;

	nQs = 0;
	delete [] Qs.x;  Qs.x = 0;
	delete [] Qs.y;  Qs.y = 0;
	delete [] Qs.w;  Qs.w = 0;

	Psit = 0;
}

cMicroscopeEffectsGPU::cMicroscopeEffectsGPU(){
	cSynCPU = ccSynCPU;

	Qt.x = 0; 
	Qt.w = 0;

	nQs = 0;
	Qs.x = 0;
	Qs.y = 0;
	Qs.w = 0;

	Psit = 0;
}

cMicroscopeEffectsGPU::~cMicroscopeEffectsGPU(){	
	freeMemory();
}

// Partially coherent transfer function and Transmission cross coefficient
void cMicroscopeEffectsGPU::PCTCCTEM(int STEffect, double2 *&fPsi, double *&M2PsiM){
	int i, iSynCPU;
	double f0 = Lens.f;
	double cf0 = Lens.cf;

	iSynCPU = 0;
	SetValueVectorD(BT.Bnxy, BT.Tnxy, GP.nxy, 0.0, M2PsiM);
	switch(STEffect){
		case 1:	// Temporal and Spatial
			for(i=0; i<nQs; i++){
				for(int j=0; j<Lens.nsf; j++){
					Lens.f = Lens.sf*Qt.x[j]+f0; 
					Lens.cf = cPi*Lens.lambda*Lens.f;
					// Apply Coherent transfer function
					ApplyCTF<<<BT.Bnxny, BT.Tnxny>>>(GP, Lens, Qs.x[i], Qs.y[i], fPsi, Psit);
					// Backward fft2
					cufftExecZ2Z(PlanPsi, Psit, Psit, CUFFT_INVERSE);
					// Apply weighting factor and add to the general sum
					AddwM2CtoD(BT.Bnxy, BT.Tnxy, true, GP.nxy, Qs.w[i]*Qt.w[j], Psit, M2PsiM);

					iSynCPU++;
					if(iSynCPU%cSynCPU==0)
						cudaDeviceSynchronize();
				}
			}
			break;
		case 2:	// Temporal
			for(int j=0; j<Lens.nsf; j++){
				Lens.f = Lens.sf*Qt.x[j]+f0; 
				Lens.cf = cPi*Lens.lambda*Lens.f;
				// Apply Coherent transfer function
				ApplyCTF<<<BT.Bnxny, BT.Tnxny>>>(GP, Lens, 0.0, 0.0, fPsi, Psit);
				// Backward fft2
				cufftExecZ2Z(PlanPsi, Psit, Psit, CUFFT_INVERSE);
				// Apply weighting factor and add to the general sum
				AddwM2CtoD(BT.Bnxy, BT.Tnxy, true, GP.nxy, Qt.w[j], Psit, M2PsiM);

				iSynCPU++;
				if(iSynCPU%cSynCPU==0)
					cudaDeviceSynchronize();
			}
			break;
		case 3:	// Spatial
			for(i=0; i<nQs; i++){
				// Apply Coherent transfer function
				ApplyCTF<<<BT.Bnxny, BT.Tnxny>>>(GP, Lens, Qs.x[i], Qs.y[i], fPsi, Psit);
				// Backward fft2
				cufftExecZ2Z(PlanPsi, Psit, Psit, CUFFT_INVERSE);
				// Apply weighting factor and add to the general sum
				AddwM2CtoD(BT.Bnxy, BT.Tnxy, true, GP.nxy, Qs.w[i], Psit, M2PsiM);

				iSynCPU++;
				if(iSynCPU%cSynCPU==0)
					cudaDeviceSynchronize();
			}
			break;
	}

	Lens.f = f0;
	Lens.cf = cf0;
}

// Partially coherent transfer function, linear image model and weak phase object
void cMicroscopeEffectsGPU::PCLIMWPOTEM(int STEffect, double2 *&fPsi, double *&M2PsiM){
	double sf = Lens.sf, beta = Lens.beta;

	switch(STEffect){
		case 2:	// Temporal
			Lens.beta = 0;
			break;
		case 3:	// Spatial
			Lens.sf = 0;
			break;
	}

	ApplyPCTF<<<BT.Bnxny, BT.Tnxny>>>(GP, Lens, fPsi, Psit);
	// Backward fft2
	cufftExecZ2Z(PlanPsi, Psit, Psit, CUFFT_INVERSE);	
	// Apply weighting factor and add to the general sum
	AddwM2CtoD(BT.Bnxy, BT.Tnxy, false, GP.nxy, 1.0, Psit, M2PsiM);

	Lens.sf = sf;
	Lens.beta = beta;
}

void cMicroscopeEffectsGPU::ReadSpatialQuadrature(sLens &Lens, int &nQ2, sQ2 &Q2){
	int i, j;
	double gxs, gys, g2s, sumwia;
	double alpha = 0.5/pow(Lens.sggs, 2);
	sQ2 Q2t;

	Q2t.x = new double [(2*Lens.ngxs+1)*(2*Lens.ngys+1)];
	Q2t.y = new double [(2*Lens.ngxs+1)*(2*Lens.ngys+1)];
	Q2t.w = new double [(2*Lens.ngxs+1)*(2*Lens.ngys+1)];
	/***********************************************************************/
	nQ2 = 0; sumwia = 0.0;
    for(j=-Lens.ngys; j<=Lens.ngys; j++)
        for(i=-Lens.ngxs; i<=Lens.ngxs; i++){
            gxs = i*Lens.dgxs; gys = j*Lens.dgys;
            g2s = gxs*gxs + gys*gys;
			if (g2s < Lens.gmax2s){
				Q2.x[nQ2] = gxs;
				Q2.y[nQ2] = gys;
				sumwia += Q2.w[nQ2] = exp(-alpha*g2s);
				nQ2++;
			}
		}
	/***********************************************************************/
	Q2.x = new double [nQ2];
	Q2.y = new double [nQ2];
	Q2.w = new double [nQ2];

	for(i=0; i<nQ2; i++){
		Q2.x[i] = Q2t.x[i];
		Q2.y[i] = Q2t.y[i];
		Q2.w[i] = Q2t.w[i]/sumwia;
	}
	/***********************************************************************/
	delete [] Q2t.x; Q2t.x = 0;
	delete [] Q2t.y; Q2t.y = 0;
	delete [] Q2t.w; Q2t.w = 0;
}

void cMicroscopeEffectsGPU::SetInputData(sBT &BT_i, sGP &GP_i, sLens &Lens_i, cufftHandle &PlanPsi_i, double2 *&Psit_i){
	freeMemory();

	BT = BT_i;
	GP = GP_i;
	Lens = Lens_i;
	PlanPsi = PlanPsi_i;
	Psit = Psit_i;
	/**********************Temporal quadrature**********************/
	Qt.x = new double[Lens.nsf]; 
	Qt.w = new double [Lens.nsf];
	cQuadrature Quad;
	Quad.ReadQuadrature(8, Lens.nsf, Qt);	// 8: int_-infty^infty f(x) Exp[-x^2] dx
	for(int i=0; i<Lens.nsf; i++)
		Qt.w[i] /= cPii2;
	/**********************Spatial quadrature**********************/
	ReadSpatialQuadrature(Lens, nQs, Qs);
}

// Inclusion of the microscope effect: TypCal 1: PCLIMWPO, 2: PCTCCTEM
void cMicroscopeEffectsGPU::ApplyMEffects(int MEffect, int STEffect, double2 *&fPsi, double *&M2Psi){
	switch (MEffect){
		case 1:
			PCLIMWPOTEM(STEffect, fPsi, M2Psi);	
			break;
		case 2:
			PCTCCTEM(STEffect, fPsi, M2Psi);	
			break;
	}
}