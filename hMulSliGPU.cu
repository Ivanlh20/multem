#include "math.h"
#include "hConstTypes.h"
#include "hgeneralCPU.h"
#include "hgeneralGPU.h"
#include "hPotentialGPU.h"
#include "hIncidentWaveGPU.h"
#include "hMulSliGPU.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include <device_functions.h>
#include "cufft.h"

/*********************************/
#include <iostream>
#include <cstdio>
#include "mex.h"
/*********************************/

// Calculated transmission function: Exp(i*sigma*Vpg)
__global__ void Transmission1(int nxy, double f, double *Vpg, double2 *Transg){
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	double theta, Tx, Ty;
	while (i < nxy){
		theta = f*Vpg[i];
		sincos(theta, &Ty , &Tx);
		Transg[i].x = Tx;
		Transg[i].y = Ty;
		i += ii;
	}
}

// Calculated transmission function: VR and zVp
template <int zpos>
__global__ void Transmission2(int nxy, double f, double *Vpg, double *zVpg, double *zVpog, double2 *Transg){
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int ii = blockDim.x*gridDim.x;
	double VR, zVp;
	double theta;
	while (i < nxy){
		VR = Vpg[i]; zVp = zVpg[i];
		switch (zpos){
			case 0:
				theta = f*(VR-zVp);
				zVpog[i] = zVp;
				break;
			case 1:
				theta = f*(VR-zVp+zVpog[i]);
				zVpog[i] = zVp;
				break;
			case 2:
				theta = f*zVpog[i];
				break;
		}	
		sincos(theta, &(Transg[i].y) , &(Transg[i].x));
		i += ii;
	}
}

// Anti-Aliasing, scale with cut-off (2/3)gmax
__global__ void AntiAliasing(sGP GP, double2 *Transg){
	int j = threadIdx.x + blockIdx.x*blockDim.x;
	int i = threadIdx.y + blockIdx.y*blockDim.y;

	if ((i < GP.nx)&&(j < GP.ny)){
		int k = i*GP.ny+j;
		double g2, gx, gy;
		gx =  ((i<GP.nxh)?i:i-GP.nx)*GP.dgx;
		gy =  ((j<GP.nyh)?j:j-GP.ny)*GP.dgy;
		g2 = gx*gx + gy*gy;

		if (g2 < GP.gmaxl2){
			Transg[k].x *= GP.inxy;
			Transg[k].y *= GP.inxy;
		}else{
			Transg[k].x = 0.0; 
			Transg[k].y = 0.0; 
		}
	}
}

// Element by element multiplication
__global__ void Transmit(sGP GP, double2 *Transg, double2 *Psig){
	int j = threadIdx.x + blockIdx.x*blockDim.x;
	int i = threadIdx.y + blockIdx.y*blockDim.y;

	if ((i < GP.nx)&&(j < GP.ny)){
		int k = i*GP.ny+j;
		double a = Transg[k].x, b = Transg[k].y;
		double c = Psig[k].x, d = Psig[k].y;
		Psig[k].x = a*c-b*d;
		Psig[k].y = b*c+a*d;
	}
}

// exp(2*1i*pi(rx*gxu+ry*gyu);
__global__ void Exp2ipiRgu(sGP GP, double gxu, double gyu, sERg ERgg){
	int k = threadIdx.x + blockIdx.x*blockDim.x;
	double theta, R;

	if (k < GP.nx){
		R =  ((k<GP.nxh)?k:k-GP.nx)*GP.dRx;
		theta = 2*cPi*R*gxu;
		sincos(theta, &(ERgg.x[k].y), &(ERgg.x[k].x));
	}

	if (k < GP.ny){
		R =  ((k<GP.nyh)?k:k-GP.ny)*GP.dRy;
		theta = 2*cPi*R*gyu;
		sincos(theta, &(ERgg.y[k].y), &(ERgg.y[k].x));
	}
}

// Phase Multiplication
__global__ void PhaseMul(sGP GP, sERg ERgg, double2 *Psig){
	int j = threadIdx.x + blockIdx.x*blockDim.x;
	int i = threadIdx.y + blockIdx.y*blockDim.y;

	if ((i < GP.nx)&&(j < GP.ny)){
		int k = i*GP.ny+j;
		double a, b, c, d, e, f;
		a = ERgg.x[i].x;	b = ERgg.x[i].y;
		c = ERgg.y[j].x;	d = ERgg.y[j].y;
		e = a*c-b*d;		f = b*c+a*d;
		c = Psig[k].x;		d = Psig[k].y; 
		Psig[k].x = a = e*c-f*d; 
		Psig[k].y = b = f*c+e*d; 
	}
}

// Build propagator function
__global__ void Propagator(sGP GP, double gxu, double gyu, double scale, sProp Propg){
	int k = threadIdx.x + blockIdx.x*blockDim.x;
	double theta, g;

	if (k < GP.nx){
		g = ((k<GP.nxh)?k:k-GP.nx)*GP.dgx;
		theta = scale*(g+gxu)*(g+gxu);
		sincos(theta, &(Propg.x[k].y), &(Propg.x[k].x));
	}
	if (k < GP.ny){
		g = ((k<GP.nyh)?k:k-GP.ny)*GP.dgy;
		theta = scale*(g+gyu)*(g+gyu);
		sincos(theta, &(Propg.y[k].y), &(Propg.y[k].x));
	}
}

// Propagate, scale with cut-off (2/3)gmax
__global__ void Propagate(sGP GP, sProp Propg, double2 *Psig)
{
	int j = threadIdx.x + blockIdx.x*blockDim.x;
	int i = threadIdx.y + blockIdx.y*blockDim.y;

	if ((i < GP.nx)&&(j < GP.ny)){
		int k = i*GP.ny+j;
		double gx, gy, g2;

		gx =  ((i<GP.nxh)?i:i-GP.nx)*GP.dgx;
		gy =  ((j<GP.nyh)?j:j-GP.ny)*GP.dgy;
		g2 = gx*gx + gy*gy;

		if (g2 < GP.gmaxl2){
			double a, b, c, d, e, f;
			a = Propg.x[i].x;	b = Propg.x[i].y;
			c = Propg.y[j].x;	d = Propg.y[j].y;
			e = a*c-b*d;		f = b*c+a*d;
			c = Psig[k].x;		d = Psig[k].y; 
			a = e*c-f*d;		b = f*c+e*d;
			Psig[k].x = a*GP.inxy;
			Psig[k].y = b*GP.inxy; 
		}
		else{
 			Psig[k].x = 0.0; 
			Psig[k].y = 0.0; 
		}
	}
}

/*********************************************************************************/
void cMulSliGPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	cSynCPU = ccSynCPU;

	gpu = 0;
	SimType = 0;
	MulOrder = 0;
	nConfFP = 0;
	DimFP = 0;
	DistFP = 0;
	SeedFP = 0;
	PotPar = 0;

	E0 = 0;
	theta = 0;
	phi = 0;

	fsGP_Init(GP);
	fsBT_Init(BT);
	fsLens_Init(Lens);
	
	gxu = 0;		
	gyu = 0;		

	fPot = 0;
	fProp = 0 ;

	/*****************************************************************/
	HRTEM.MEffect = 0;
	HRTEM.STEffect = 0;
	HRTEM.ZeroDefTyp = 0;
	HRTEM.ZeroDefPlane = 0;

	PED.nrot = 0;
	PED.theta = 0;

	HCI.nrot = 0;
	HCI.theta = 0;
	HCI.MEffect = 0;
	HCI.STEffect = 0;
	HCI.ZeroDefTyp = 0;
	HCI.ZeroDefPlane = 0;

	EWRS.ZeroDefTyp = 0;
	EWRS.ZeroDefPlane = 0;

	EWFS.ZeroDefTyp = 0;
	EWFS.ZeroDefPlane = 0;
	/*****************************************************************/

	cudaFreen(ERg.x);
	cudaFreen(ERg.y);

	cudaFreen(Prop.x);
	cudaFreen(Prop.y);

	cudaFreen(Trans);

	cudaFreen(Psi);
	cudaFreen(aPsi);

	cudaFreen(M2aPsi);
	cudaFreen(aM2Psi);

	cufftDestroyn(PlanPsi);
}

cMulSliGPU::cMulSliGPU(){
	cSynCPU = ccSynCPU;

	gpu = 0;
	SimType = 0;
	MulOrder = 0;
	nConfFP = 0;
	DimFP = 0;
	DistFP = 0;
	SeedFP = 0;
	PotPar = 0;

	E0 = 0;
	theta = 0;
	phi = 0;

	fsGP_Init(GP);
	fsBT_Init(BT);
	fsLens_Init(Lens);
	
	gxu = 0;		
	gyu = 0;		

	fPot = 0;
	fProp = 0 ;

	/*****************************************************************/
	HRTEM.MEffect = 0;
	HRTEM.STEffect = 0;
	HRTEM.ZeroDefTyp = 0;
	HRTEM.ZeroDefPlane = 0;

	PED.nrot = 0;
	PED.theta = 0;

	HCI.nrot = 0;
	HCI.theta = 0;
	HCI.MEffect = 0;
	HCI.STEffect = 0;
	HCI.ZeroDefTyp = 0;
	HCI.ZeroDefPlane = 0;

	EWRS.ZeroDefTyp = 0;
	EWRS.ZeroDefPlane = 0;

	EWFS.ZeroDefTyp = 0;
	EWFS.ZeroDefPlane = 0;
	/*****************************************************************/

	ERg.x = 0;
	ERg.y = 0;

	Prop.x = 0;
	Prop.y = 0;

	Trans = 0;

	Psi = 0;
	aPsi = 0;

	M2aPsi = 0;
	aM2Psi = 0;

	PlanPsi = 0;
}

cMulSliGPU::~cMulSliGPU(){
	freeMemory();
}

// Calculate transmission function
void cMulSliGPU::Transmission(dim3 grid, dim3 threads, int MO, int SlicePos, int nxy, double f, double *Vpg, double *zVpg, double *zVpog, double2 *Transg){
	switch (MO){
		case 1:
			Transmission1<<<grid, threads>>>(nxy, f, Vpg, Transg);
			break;
		case 2:	
			switch (SlicePos){
				case 0:
					Transmission2<0><<<grid, threads>>>(nxy, f, Vpg, zVpg, zVpog, Transg);
					break;
				case 1:
					Transmission2<1><<<grid, threads>>>(nxy, f, Vpg, zVpg, zVpog, Transg);
					break;
				case 2:
					Transmission2<2><<<grid, threads>>>(nxy, f, Vpg, zVpg, zVpog, Transg);
					break;
			}	
			break;
	}
}

// AntiAliasing
void cMulSliGPU::Bandlimited2D(dim3 grid, dim3 threads, cufftHandle &PlanPsi, sGP &GP, double2 *Psi){
	// Forward fft2
	cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
	// AntiAliasing, scale and bandlimited the transmission function
	AntiAliasing<<<grid, threads>>>(GP, Psi);
	// Backward fft2
	cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_INVERSE);
}

// From Device To Host
void cMulSliGPU::Gpu2CPU(double2 *&Psid, sComplex &Psih){	
	// fft2shift 
	fft2ShiftC(BT.Bhnxny, BT.Thnxny, GP.nxh, GP.nyh, Psid);
	/**********************copy data to host************************/
	double *Psir, *Psii;
	Psir = PotentialGPU.V0; Psii = PotentialGPU.V1;
	// Get Real and Imaginary part to complex matrix
	GetRealImagVectorC(BT.Bnxy, BT.Tnxy, GP.nxy, Psid, Psir, Psii);
	// Copy real part of the wave function to the host
	cudaMemcpy(Psih.real, Psir, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
	// Copy imaginary part of the wave function to the host
	cudaMemcpy(Psih.imag, Psii, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
}

// From Device To Host
void cMulSliGPU::Gpu2CPU(double2 *&Psid, double *&M2Psid, sComplex &Psih, double *&M2Psih){
	Gpu2CPU(Psid, Psih);
	// fft2shift 
	fft2ShiftD(BT.Bhnxny, BT.Thnxny, GP.nxh, GP.nyh, M2Psid);
	/**********************copy data to host************************/
	// Copy wave function squared to the host
	cudaMemcpy(M2Psih, M2Psid, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
}

// From Device To Host
void cMulSliGPU::Gpu2CPU(double2 *&Psid, double *&M2Psi1d, double *M2Psi2d, sComplex &Psih, double *&M2Psi1h, double *&M2Psi2h){
	Gpu2CPU(Psid, Psih);
	// fft2shift 
	fft2ShiftD(BT.Bhnxny, BT.Thnxny, GP.nxh, GP.nyh, M2Psi1d);
	fft2ShiftD(BT.Bhnxny, BT.Thnxny, GP.nxh, GP.nyh, M2Psi2d);
	/**********************copy data to host************************/
	// Copy wave function squared to the host
	cudaMemcpy(M2Psi1h, M2Psi1d, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
	cudaMemcpy(M2Psi2h, M2Psi2d, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
}

// Exit wave calculation in Fourier space: Space :1 Real and 2 Fourier
void cMulSliGPU::MSExitWaveCal(int iConf, int iSpace, double2 *&Psi){
	int i, iSynCPU;

	/**************************************/
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start);
	/**************************************/

	iSynCPU = 0;
	PotentialGPU.Specimen.MoveAtoms(iConf);
	for (i=0; i<PotentialGPU.Specimen.nz; i++){
		//// Projected potential
		//PotentialGPU.ProjectedPotential(i);
		///***********************************************************************/
		//// Transmission
		//Transmission(BT.Bnxy, BT.Tnxy, MulOrder, (i==0)?0:1, GP.nxy, fPot, PotentialGPU.V0, PotentialGPU.V1, PotentialGPU.V2, Trans);
		//// AntiAliasing
		//Bandlimited2D(BT.Bnxny, BT.Tnxny, PlanPsi, GP, Trans);
		//// Transmit
		Transmit<<<BT.Bnxny, BT.Tnxny>>>(GP, Trans, Psi);
		///***********************************************************************/
		//fProp = cPi*(PotentialGPU.Specimen.Slice[i].ze - PotentialGPU.Specimen.Slice[i].z0)*Lens.lambda/cos(theta);
		//// Forward propagator
		//Propagator<<<BT.Bmnxny, BT.Tmnxny>>>(GP, gxu, gyu, -fProp, Prop);
		//// Forward fft2
		//cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
		//// Propagate, scale and bandlimited the wave function
		//Propagate<<<BT.Bnxny, BT.Tnxny>>>(GP, Prop, Psi);
		//// Backward fft2
		//cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_INVERSE);

		//iSynCPU++;
		//if(iSynCPU%cSynCPU==0)
		//	cudaDeviceSynchronize();
	}

	/**************************************/
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);
	mexPrintf("%4.3f\n", milliseconds);
	/**************************************/

	/**************************last transmission function***********************/
	if (MulOrder==2){
		// Transmission
		Transmission(BT.Bnxy, BT.Tnxy, MulOrder, 2, GP.nxy, fPot, PotentialGPU.V0, PotentialGPU.V1, PotentialGPU.V2, Trans);
		// AntiAliasing
		Bandlimited2D(BT.Bnxny, BT.Tnxny, PlanPsi, GP, Trans);
		// Transmit
		Transmit<<<BT.Bnxny, BT.Tnxny>>>(GP, Trans, Psi);
	}

	/****************************Inclined ilumination***************************/
	if ((gxu != 0.0)||(gyu != 0.0)){
		// exp(2*1i*pi(rx*gxu+ry*gyu)
		Exp2ipiRgu<<<BT.Bmnxny, BT.Tmnxny>>>(GP, gxu, gyu, ERg);
		// Phase multiplication
		PhaseMul<<<BT.Bnxny, BT.Tnxny>>>(GP, ERg, Psi);
	}

	/******************Back propagation to the Reference point****************/
	// Build back propagator(To our reference point(zero defocus))
	Propagator<<<BT.Bmnxny, BT.Tmnxny>>>(GP, gxu, gyu, -cPi*PotentialGPU.Specimen.dzBackProp*Lens.lambda, Prop);
	// Forward fft2
	cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
	// Propagate and scale the wave function
	Propagate<<<BT.Bnxny, BT.Tnxny>>>(GP, Prop, Psi);

	// Backward fft2
	if(iSpace==1)
		cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_INVERSE);
}

// Get wave function for inclined ilumination by using frozen phonon
void cMulSliGPU::Wavefuncion(int nConfFP, int Space, int MEffect, int STEffect, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi){
	int iConf0 = (nConfFP==0)?0:1;
	bool IsNotSharp = (nConfFP==0)?false:true;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);

	if(IsNotSharp)
		SetValueVectorCD(BT.Bnxy, BT.Tnxy, GP.nxy, 0.0, 0.0, aPsi, 0.0, aM2Psi);

	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Plane wave ilumination
		IncidentWaveGPU.Psi0(Psi);

		if(MEffect==0){
			//Exit wave
			MSExitWaveCal(iConf, Space, Psi);
			// Add Psi and M2aPsi to aPsi and aM2Psi
			AddwCwM2CtoCD(BT.Bnxy, BT.Tnxy, IsNotSharp, GP.nxy, inConfFP, Psi, aPsi, aM2Psi);
		}else{
			//Exit wave(g)
			MSExitWaveCal(iConf, 2, Psi);
			// Microscope effects
			MicroscopeEffectsGPU.ApplyMEffects(MEffect, STEffect, Psi, M2aPsi);
			// Add Psi and M2aPsi to aPsi and aM2Psi
			AddwCwDtoCD(BT.Bnxy, BT.Tnxy, IsNotSharp, GP.nxy, inConfFP, Psi, M2aPsi, aPsi, aM2Psi);
		}
	}

	if(MEffect==0)
		 AddwM2CtoD(BT.Bnxy, BT.Tnxy, false, GP.nxy, 1.0, aPsi, M2aPsi);
	else{
		// Forward fft2
		cufftExecZ2Z(PlanPsi, aPsi, aPsi, CUFFT_FORWARD);
		// Scale vector
		ScaleVectorC(BT.Bnxy, BT.Tnxy, GP.nxy, GP.inxy, aPsi);
		// Microscope effects
		MicroscopeEffectsGPU.ApplyMEffects(MEffect, STEffect, aPsi, M2aPsi);
	}

}

// Set Input data
void cMulSliGPU::SetInputData(sInMSTEM &InMSTEM){
	freeMemory();

	gpu = InMSTEM.gpu;
	SimType = InMSTEM.SimType;	
	MulOrder = InMSTEM.MulOrder;	
	nConfFP = InMSTEM.nConfFP;		
	DimFP = InMSTEM.DimFP;	
	DistFP = 1;
	SeedFP = InMSTEM.SeedFP;
	PotPar = InMSTEM.PotPar;
	E0 = InMSTEM.E0;	
	theta = InMSTEM.theta;	
	phi = InMSTEM.phi;

	cudaSetDevice(gpu);

	fsGP_Cal(InMSTEM.nx, InMSTEM.ny, InMSTEM.lx, InMSTEM.ly, InMSTEM.dz, GP);

	fsBT_Cal(GP, BT);

	fsLens_Cal(E0, GP, Lens);

	gxu = sin(theta)*cos(phi)/Lens.lambda;
	gyu = sin(theta)*sin(phi)/Lens.lambda;

	fPot = Lens.gamma*Lens.lambda/cos(theta);
	fProp = cPi*InMSTEM.dz*Lens.lambda/cos(theta);

	/*****************************************************************/
	HRTEM.MEffect = InMSTEM.HRTEM_MEffect;
	HRTEM.STEffect = InMSTEM.HRTEM_STEffect;
	HRTEM.ZeroDefTyp = InMSTEM.HRTEM_ZeroDefTyp;
	HRTEM.ZeroDefPlane = InMSTEM.HRTEM_ZeroDefPlane;

	PED.nrot = InMSTEM.PED_nrot;
	PED.theta = InMSTEM.PED_theta;

	HCI.nrot = InMSTEM.HCI_nrot;
	HCI.theta = InMSTEM.HCI_theta;
	HCI.MEffect = InMSTEM.HCI_MEffect;
	HCI.STEffect = InMSTEM.HCI_STEffect;
	HCI.ZeroDefTyp = InMSTEM.HCI_ZeroDefTyp;
	HCI.ZeroDefPlane = InMSTEM.HCI_ZeroDefPlane;

	EWRS.ZeroDefTyp = InMSTEM.EWRS_ZeroDefTyp;
	EWRS.ZeroDefPlane = InMSTEM.EWRS_ZeroDefPlane;

	EWFS.ZeroDefTyp = InMSTEM.EWFS_ZeroDefTyp;
	EWFS.ZeroDefPlane = InMSTEM.EWFS_ZeroDefPlane;
	/*****************************************************************/

	cudaMalloc((void**)&(ERg.x), GP.nx*cSizeofCD);
	cudaMalloc((void**)&(ERg.y), GP.ny*cSizeofCD);

	cudaMalloc((void**)&(Prop.x), GP.nx*cSizeofCD);
	cudaMalloc((void**)&(Prop.y), GP.ny*cSizeofCD);

	cudaMalloc((void**)&Trans, GP.nxy*cSizeofCD);

	cudaMalloc((void**)&Psi, GP.nxy*cSizeofCD);
	cudaMalloc((void**)&aPsi, GP.nxy*cSizeofCD);

	cudaMalloc((void**)&M2aPsi, GP.nxy*cSizeofRD);
	cudaMalloc((void**)&aM2Psi, GP.nxy*cSizeofRD);

	cufftPlan2d(&PlanPsi, GP.nx, GP.ny, CUFFT_Z2Z);

	// Potential parameters
	int ZeroDefTyp = 3;
	double ZeroDefPlane = 0;
	switch (SimType){
		case 3:
			ZeroDefTyp = HRTEM.ZeroDefTyp;
			ZeroDefPlane = HRTEM.ZeroDefPlane;
			break;
		case 6:
			ZeroDefTyp = HCI.ZeroDefTyp;
			ZeroDefPlane = HCI.ZeroDefPlane;
			break;
		case 10:
			ZeroDefTyp = EWRS.ZeroDefTyp;
			ZeroDefPlane = EWRS.ZeroDefPlane;
			break;
		case 11:
			ZeroDefTyp = EWFS.ZeroDefTyp;
			ZeroDefPlane = EWFS.ZeroDefPlane;
			break;
	}
	PotentialGPU.SetInputData(GP, BT, InMSTEM.nAtomsM, InMSTEM.AtomsM, InMSTEM.dz, PotPar, stVrl, DimFP, DistFP, SeedFP, ZeroDefTyp, ZeroDefPlane, MulOrder);

	// Set parameters for the incident wave
	IncidentWaveGPU.SetInputData(BT, GP, Lens, PlanPsi);

	// Microscope parameters
	MicroscopeEffectsGPU.SetInputData(BT, GP, Lens, PlanPsi, Trans);
}

void cMulSliGPU::Cal_ExitWavefuncionRS(sComplex &aPsih, double *&aM2Psih){
	int MEffect = 0, STEffect = 0, Space = 1;
	Wavefuncion(nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	Gpu2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMulSliGPU::Cal_ExitWavefuncionFS(sComplex &aPsih, double *&aM2Psih){
	int MEffect = 0, STEffect = 0, Space = 2;
	Wavefuncion(nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	Gpu2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMulSliGPU::Cal_ED(sComplex &aPsih, double *&aM2Psih){
	int MEffect = 0, STEffect = 0, Space = 2;
	Wavefuncion(nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	Gpu2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMulSliGPU::Cal_HRTEM(sComplex &aPsih, double *&M2aPsih, double *&aM2Psih){
	int Space = 1;
	Wavefuncion(nConfFP, Space, HRTEM.MEffect, HRTEM.STEffect, aPsi, M2aPsi, aM2Psi);
	Gpu2CPU(aPsi, M2aPsi, aM2Psi, aPsih, M2aPsih, aM2Psih);
}