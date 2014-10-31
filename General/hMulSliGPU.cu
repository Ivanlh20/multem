#include <memory>
#include <cmath>
#include "hConstTypes.h"
#include "hgeneralCPU.h"
#include "hgeneralGPU.h"
#include "hPotentialGPU.h"
#include "hIncidentWaveGPU.h"
#include "hDetectorCPU.h"
#include "hDetectorGPU.h"
#include "hMulSliGPU.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <cufft.h>

void sSTEM::Init()
{
	line = 0;
	FastCal = false;
	ns = 0;
	x1u = 0;
	y1u = 0;		
	x2u = 0;	
	y2u = 0;	
	nThk = 0;

	f_sDetCir_Init(DetCir);
	ImSTEM = 0;
	nDet = 0;
	nxs = 0;
	nys = 0;
	xs = 0;
	ys = 0;

	nst = 0;
	xst = 0;
	yst = 0;
}

void sSTEM::free()
{
	line = 0;
	FastCal = false;
	ns = 0;	
	x1u = 0;
	y1u = 0;		
	x2u = 0;	
	y2u = 0;	
	f_sDetCir_Free(DetCir);

	for (int iThk=0; iThk<nThk; iThk++){
		for (int iDet=0; iDet<nDet; iDet++){
			delete [] ImSTEM[iThk].DetInt[iDet].Coh; ImSTEM[iThk].DetInt[iDet].Coh = 0;
			delete [] ImSTEM[iThk].DetInt[iDet].Tot; ImSTEM[iThk].DetInt[iDet].Tot = 0;
		}
		delete [] ImSTEM[iThk].DetInt; ImSTEM[iThk].DetInt = 0;
	}
	delete [] ImSTEM; ImSTEM = 0;

	nxs = 0;
	nys = 0;
	delete [] xs;  xs = 0;
	delete [] ys;  ys = 0;

	nst = 0;
	delete [] xst;  xst = 0;
	delete [] yst;  yst = 0;
}

void sSTEM::SetInputData(sInMSTEM &InMSTEM)
{
	free();

	line = InMSTEM.STEM_line;
	FastCal = InMSTEM.STEM_FastCal;
	ns = InMSTEM.STEM_ns;
	x1u = InMSTEM.STEM_x1u;	
	y1u = InMSTEM.STEM_y1u;
	x2u = InMSTEM.STEM_x2u;
	y2u = InMSTEM.STEM_y2u;
	f_BuildGrid(line, ns, x1u, y1u, x2u, y2u, nxs, nys, xs, ys);
	nThk = 1;

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

	ImSTEM = new sImSTEM[nThk];
	for (int iThk=0; iThk<nThk; iThk++){
		ImSTEM[iThk].DetInt = new sDetInt[nDet];
		for (int iDet=0; iDet<nDet; iDet++){
			ImSTEM[iThk].DetInt[iDet].Coh = new double[nst];
			ImSTEM[iThk].DetInt[iDet].Tot = new double[nst];
		}
	}
}

void sSTEM::InitImSTEM()
{
	int iThk, iDet, ist;
	for (iThk=0; iThk<nThk; iThk++)
		for (iDet=0; iDet<nDet; iDet++)
			for (ist=0; ist<nst; ist++){
				ImSTEM[iThk].DetInt[iDet].Coh[ist] = 0.0;
				ImSTEM[iThk].DetInt[iDet].Tot[ist] = 0.0;
			}
}

/****************************************************************************************/
void cMulSliGPU::freeMemory(){
	cudaDeviceSynchronize(); // wait to finish the work in the GPU

	cSynCPU = ccSynCPU;

	f_sMPG_Init(MGP);
	f_sGP_Init(GP);
	f_sBT_Init(BT);
	f_sLens_Init(Lens);
	
	gxu = 0;
	gyu = 0;

	fPot = 0;

	/*****************************************************************/
	STEM.free();

	CBED.x0 = 0;
	CBED.y0 = 0;

	//HRTEM.xx;

	PED.nrot = 0;
	PED.theta = 0;

	HCI.nrot = 0;
	HCI.theta = 0;
	//HCI.xx;

	//EWRS.xx;

	//EWFS.xx;
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

	f_sMPG_Init(MGP);
	f_sGP_Init(GP);
	f_sBT_Init(BT);
	f_sLens_Init(Lens);
	
	gxu = 0;		
	gyu = 0;		

	fPot = 0;

	/*****************************************************************/
	STEM.Init();

	CBED.x0 = 0;
	CBED.y0 = 0;

	//HRTEM.xx

	PED.nrot = 0;
	PED.theta = 0;

	HCI.nrot = 0;
	HCI.theta = 0;
	//HCI.xx

	//EWRS.xx

	//EWFS.xx
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
	cudaDeviceReset();
}

// From Device To Host
void cMulSliGPU::Gpu2CPU(double2 *&Psid, sComplex &Psih){	
	// fft2shift 
	f_fft2ShiftC(BT.Bhnxny, BT.Thnxny, GP.nxh, GP.nyh, Psid);
	/**********************copy data to host************************/
	Cd_2_Ch(GP.nxy, PotentialGPU.V0, PotentialGPU.V1, Psid, Psih);
}

// From Device To Host
void cMulSliGPU::Gpu2CPU(double2 *&Psid, double *&M2Psid, sComplex &Psih, double *&M2Psih){
	// fft2shift 
	f_fft2ShiftC(BT.Bhnxny, BT.Thnxny, GP.nxh, GP.nyh, Psid); 
	f_fft2ShiftD(BT.Bhnxny, BT.Thnxny, GP.nxh, GP.nyh, M2Psid);

	// Copy wave function squared to the host
	CdDd_2_ChDh(GP.nxy, PotentialGPU.V0, PotentialGPU.V1, Psid, M2Psid, Psih, M2Psih);
}

// From Device To Host
void cMulSliGPU::Gpu2CPU(double2 *&Psid, double *&M2Psi1d, double *&M2Psi2d, sComplex &Psih, double *&M2Psi1h, double *&M2Psi2h){
	// fft2shift 
	f_fft2ShiftC(BT.Bhnxny, BT.Thnxny, GP.nxh, GP.nyh, Psid); 
	f_fft2ShiftD(BT.Bhnxny, BT.Thnxny, GP.nxh, GP.nyh, M2Psi1d);
	f_fft2ShiftD(BT.Bhnxny, BT.Thnxny, GP.nxh, GP.nyh, M2Psi2d);
	/**********************copy data to host************************/
	CdDdDd_2_ChDhDh(GP.nxy, PotentialGPU.V0, PotentialGPU.V1, Psid, M2Psi1d, M2Psi2d, Psih, M2Psi1h, M2Psi2h);
}

void cMulSliGPU::PhaseMul(double gxu, double gyu, double2 *&Psi){
	if ((gxu != 0.0)||(gyu != 0.0))
		f_PhaseMul(BT.Bmnxny, BT.Tmnxny, GP, gxu, gyu, ERg, Psi);
}

void cMulSliGPU::TransmissionWPO(double fPot, double2 *&Trans){
	f_TransmissionWPO(BT.Bnxny, BT.Tnxny, PlanPsi, GP, fPot, PotentialGPU.V0, Trans);
}

void cMulSliGPU::Transmission(eSlicePos SlicePos, double fPot, double2 *&Trans){
	switch (MGP.MulOrder)
	{
		case 1:
			f_Transmission1(BT.Bnxny, BT.Tnxny, PlanPsi, GP, fPot, PotentialGPU.V0, Trans);
			break;
		case 2:
			f_Transmission2(BT.Bnxny, BT.Tnxny, PlanPsi, GP, fPot, PotentialGPU.V0, PotentialGPU.V1, PotentialGPU.V2, SlicePos, Trans);
			break;
	}
}

void cMulSliGPU::Transmit(double2 *&Trans, double2 *&Psi){
	f_Transmit(BT.Bnxny, BT.Tnxny, GP, Trans, Psi);
}

void cMulSliGPU::Propagate(eSpace Space, double gxu, double gyu, double z, double2 *&Psi){
	f_Propagate(BT.Bnxny, BT.Tnxny, PlanPsi, GP, Space, gxu, gyu, Lens.lambda, z, Prop, Psi);
}

/***************************************************************************/
// Phase object approximation: Space :1 Real and 2 Fourier
void cMulSliGPU::Cal_Wavefunction_WPOA(eSpace Space, double2 *&Psi){
	// Projected potential
	PotentialGPU.ProjectedPotential_Specimen(0);
	// Transmission
	TransmissionWPO(fPot, Trans);
	// Transmit
	Transmit(Trans, Psi);
	// Inclined ilumination
	PhaseMul(gxu, gyu, Psi);
	// Backward fft2
	if(Space==eSReciprocal){
		cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
		f_ScaleVectorC(BT.Bnxy, BT.Tnxy, GP.nxy, GP.inxy, Psi);
	}
}

// Phase object approximation: Space :1 Real and 2 Fourier
void cMulSliGPU::Cal_Wavefunction_POA(eSpace Space, double2 *&Psi){
	// Projected potential
	PotentialGPU.ProjectedPotential_Specimen(0);
	// Transmission
	Transmission(eSPFirst, fPot, Trans);
	// Transmit
	Transmit(Trans, Psi);
	// Inclined ilumination
	PhaseMul(gxu, gyu, Psi);
	// Backward fft2
	if(Space==eSReciprocal){
		cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
		f_ScaleVectorC(BT.Bnxy, BT.Tnxy, GP.nxy, GP.inxy, Psi);
	}
}

// Phase object approximation: Space :1 Real and 2 Fourier
void cMulSliGPU::Cal_Wavefunction_PA(eSpace Space, double2 *&Psi){
	int iSlice, iSynCPU = 0;
	int nSlice = ceil(PotentialGPU.Specimen.Lzu/GP.dz);
	double dz = PotentialGPU.Specimen.Lzu/double(nSlice);
	// Projected potential
	PotentialGPU.ProjectedPotential_Specimen(0);
	// Transmission
	Transmission(eSPFirst, fPot/double(nSlice), Trans);
	// Multislice
	for (iSlice=0; iSlice<nSlice; iSlice++){
		// Transmit
		Transmit(Trans, Psi);
		// Propagate
		Propagate(eSReal, gxu, gyu, dz, Psi);
		iSynCPU++;
		if(iSynCPU%cSynCPU==0)
			cudaDeviceSynchronize();
	}
	// Inclined ilumination
	PhaseMul(gxu, gyu, Psi);
	// Backward fft2
	if(Space==eSReciprocal){
		cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
		f_ScaleVectorC(BT.Bnxy, BT.Tnxy, GP.nxy, GP.inxy, Psi);
	}
}

// Exit wave calculation: Space :1 Real and 2 Fourier
void cMulSliGPU::Cal_Wavefunction_MSA(eSpace Space, double2 *&Psi){
	int iSlice, iSynCPU;
	double dz;
	iSynCPU = 0;
	for (iSlice=0; iSlice<PotentialGPU.Specimen.nSlice; iSlice++){
		// Projected potential
		PotentialGPU.ProjectedPotential_Slice(iSlice);
		// Transmission
		Transmission((iSlice==0)?eSPFirst:eSPMedium, fPot, Trans);
		// Transmit
		Transmit(Trans, Psi);
		// Propagate
		dz = (PotentialGPU.Specimen.Slice[iSlice].ze - PotentialGPU.Specimen.Slice[iSlice].z0)/cos(MGP.theta);
		Propagate(eSReal, gxu, gyu, dz, Psi);
		iSynCPU++;
		if(iSynCPU%cSynCPU==0)
			cudaDeviceSynchronize();
	}
	// Last transmission function
	if (MGP.MulOrder==2){
		// Transmission
		Transmission(eSPLast, fPot, Trans);
		// Transmit
		Transmit(Trans, Psi);
	}
	// Inclined ilumination
	PhaseMul(gxu, gyu, Psi);
	// Back propagation to our plane reference
	Propagate(Space, gxu, gyu, PotentialGPU.Specimen.z_BackProp, Psi);
}

/***************************************************************************/
// Phase object approximation: Space :1 Real and 2 Fourier
void cMulSliGPU::Cal_FAST_STEM_Wavefunction_WPOA(int nConfFP, sDetInt *DetInt){
	int ist, iThk=0;
	int iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);
	double nxy2 = pow(double(GP.nxy), 2);

	cDetectorGPU DetectorGPU;
	DetectorGPU.SetInputData(GP, BT, STEM.nDet, STEM.DetCir);

	STEM.InitImSTEM();
	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Move atoms
		PotentialGPU.Specimen.MoveAtoms(iConf);
		// Projected potential
		PotentialGPU.ProjectedPotential_Specimen(0);
		// Transmission
		TransmissionWPO(fPot, Trans);
		for (ist=0; ist<STEM.nst; ist++){
			// Plane wave ilumination
			IncidentWaveGPU.Psi0(STEM.xst[ist], STEM.yst[ist], Psi);
			// Transmit
			Transmit(Trans, Psi);
			// Inclined ilumination
			PhaseMul(gxu, gyu, Psi);
			// Backward fft2
			cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
			// Add Psi to aM2Psi
			f_AddwM2CtoD(BT.Bnxy, BT.Tnxy, false, GP.nxy, inConfFP/nxy2, Psi, aM2Psi);
			// Detector integration
			DetectorGPU.getDetectorIntensity(aM2Psi, ist, STEM.ImSTEM[iThk].DetInt, true);
		}
	}
}

// Phase object approximation: Space :1 Real and 2 Fourier
void cMulSliGPU::Cal_FAST_STEM_Wavefunction_POA(int nConfFP, sDetInt *DetInt){
	int ist, iThk=0;
	int iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);
	double nxy2 = pow(double(GP.nxy), 2);

	cDetectorGPU DetectorGPU;
	DetectorGPU.SetInputData(GP, BT, STEM.nDet, STEM.DetCir);

	STEM.InitImSTEM();
	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Move atoms
		PotentialGPU.Specimen.MoveAtoms(iConf);
		// Projected potential
		PotentialGPU.ProjectedPotential_Specimen(0);
		// Transmission
		Transmission(eSPFirst, fPot, Trans);
		for (ist=0; ist<STEM.nst; ist++){
			// Plane wave ilumination
			IncidentWaveGPU.Psi0(STEM.xst[ist], STEM.yst[ist], Psi);
			// Transmit
			Transmit(Trans, Psi);
			// Inclined ilumination
			PhaseMul(gxu, gyu, Psi);
			// Backward fft2
			cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
			// Add Psi to aM2Psi
			f_AddwM2CtoD(BT.Bnxy, BT.Tnxy, false, GP.nxy, inConfFP/nxy2, Psi, aM2Psi);
			// Detector integration
			DetectorGPU.getDetectorIntensity(aM2Psi, ist, STEM.ImSTEM[iThk].DetInt, true);
		}
	}
}

// Phase object approximation: Space :1 Real and 2 Fourier
void cMulSliGPU::Cal_FAST_STEM_Wavefunction_PA(int nConfFP, sDetInt *DetInt){
	int iSlice, ist, iSynCPU, iThk=0;
	int iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);
	double nxy2 = pow(double(GP.nxy), 2);
	int nSlice = ceil(PotentialGPU.Specimen.Lzu/GP.dz);
	double dz = PotentialGPU.Specimen.Lzu/double(nSlice);

	cDetectorGPU DetectorGPU;
	DetectorGPU.SetInputData(GP, BT, STEM.nDet, STEM.DetCir);

	STEM.InitImSTEM();
	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Move atoms
		PotentialGPU.Specimen.MoveAtoms(iConf);
		// Projected potential
		PotentialGPU.ProjectedPotential_Specimen(0);
		// Transmission
		Transmission(eSPFirst, fPot/double(nSlice), Trans);
		for (ist=0; ist<STEM.nst; ist++){
			// Plane wave ilumination
			IncidentWaveGPU.Psi0(STEM.xst[ist], STEM.yst[ist], Psi);
			// Multislice
			for (iSlice=0; iSlice<nSlice; iSlice++){
				// Transmit
				Transmit(Trans, Psi);
				// Propagate
				Propagate(eSReal, gxu, gyu, dz, Psi);
				iSynCPU++;
				if(iSynCPU%cSynCPU==0)
					cudaDeviceSynchronize();
			}
			// Inclined ilumination
			PhaseMul(gxu, gyu, Psi);
			// Backward fft2
			cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_FORWARD);
			// Add Psi to aM2Psi
			f_AddwM2CtoD(BT.Bnxy, BT.Tnxy, false, GP.nxy, inConfFP/nxy2, Psi, aM2Psi);
			// Detector integration
			DetectorGPU.getDetectorIntensity(aM2Psi, ist, STEM.ImSTEM[iThk].DetInt, true);
		}
	}
}

// Exit wave calculation: Space :1 Real and 2 Fourier
void cMulSliGPU::Cal_FAST_STEM_Wavefunction_MSA(int nConfFP, sDetInt *DetInt){
	int ist, iThk=0;
	int iConf0 = (nConfFP==0)?0:1;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);

	cDetectorGPU DetectorGPU;
	DetectorGPU.SetInputData(GP, BT, STEM.nDet, STEM.DetCir);

	STEM.InitImSTEM();
	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Move atoms
		PotentialGPU.Specimen.MoveAtoms(iConf);
		for (ist=0; ist<STEM.nst; ist++){
			// Plane wave ilumination
			IncidentWaveGPU.Psi0(STEM.xst[ist], STEM.yst[ist], Psi);
			// Exit wave
			Cal_Wavefunction(eSReciprocal, Psi);
			// Add Psi to aM2Psi
			f_AddwM2CtoD(BT.Bnxy, BT.Tnxy, false, GP.nxy, inConfFP, Psi, aM2Psi);
			// Detector integration
			DetectorGPU.getDetectorIntensity(aM2Psi, ist, STEM.ImSTEM[iThk].DetInt, true);
		}
	}
}

/***************************************************************************/
// Exit wave calculation: Space :1 Real and 2 Fourier
void cMulSliGPU::Cal_Wavefunction(eSpace Space, double2 *&Psi){
	switch(MGP.ApproxModel)
	{
		case 1:
			Cal_Wavefunction_MSA(Space, Psi);
			break;
		case 2:
			Cal_Wavefunction_PA(Space, Psi);
			break;
		case 3:
			Cal_Wavefunction_POA(Space, Psi);
			break;
		case 4:
			Cal_Wavefunction_WPOA(Space, Psi);
			break;
	}
}

/***************************************************************************/
// Get wave function for inclined ilumination by using frozen phonon
void cMulSliGPU::get_Imagefuncion(int nConfFP, eSpace Space, int MEffect, int STEffect, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi){
	int iConf0 = (nConfFP==0)?0:1;
	bool IsNotSharp = (nConfFP==0)?false:true;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);

	if(IsNotSharp)
		f_SetValueVectorCD(BT.Bnxy, BT.Tnxy, GP.nxy, 0.0, 0.0, aPsi, 0.0, aM2Psi);

	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Move atoms
		PotentialGPU.Specimen.MoveAtoms(iConf);
		// Plane wave ilumination
		IncidentWaveGPU.Psi0(Psi);
		if(MEffect==0){
			// Exit wave
			Cal_Wavefunction(Space, Psi);
			// Add Psi and M2aPsi to aPsi and aM2Psi
			f_AddwCwM2CtoCD(BT.Bnxy, BT.Tnxy, IsNotSharp, GP.nxy, inConfFP, Psi, aPsi, aM2Psi);
		}else{
			// Exit wave(g)
			Cal_Wavefunction(eSReciprocal, Psi);
			// Microscope effects
			MicroscopeEffectsGPU.ApplyMEffects(MEffect, STEffect, Psi, M2aPsi);
			// Backward fft2
			cufftExecZ2Z(PlanPsi, Psi, Psi, CUFFT_INVERSE);
			// Add Psi and M2aPsi to aPsi and aM2Psi
			f_AddwCwDtoCD(BT.Bnxy, BT.Tnxy, IsNotSharp, GP.nxy, inConfFP, Psi, M2aPsi, aPsi, aM2Psi);
		}
	}

	if(MEffect==0)
		f_AddwM2CtoD(BT.Bnxy, BT.Tnxy, false, GP.nxy, 1.0, aPsi, M2aPsi);
	else{
		// Forward fft2
		cufftExecZ2Z(PlanPsi, aPsi, aPsi, CUFFT_FORWARD);
		// Scale vector
		f_ScaleVectorC(BT.Bnxy, BT.Tnxy, GP.nxy, GP.inxy, aPsi);
		// Microscope effects
		MicroscopeEffectsGPU.ApplyMEffects(MEffect, STEffect, aPsi, M2aPsi);
		// Backward fft2
		cufftExecZ2Z(PlanPsi, aPsi, aPsi, CUFFT_INVERSE);
	}

}

// Get wave function for convergent beam and inclined ilumination using frozen phonon
void cMulSliGPU::get_Imagefuncion(int nConfFP, eSpace Space, double xi, double yi, double2 *&aPsi, double *&M2aPsi, double *&aM2Psi){
	int iConf0 = (nConfFP==0)?0:1;
	bool IsNotSharp = (nConfFP==0)?false:true;
	double inConfFP = (nConfFP==0)?1.0:1.0/double(nConfFP);

	if(IsNotSharp)
		f_SetValueVectorCD(BT.Bnxy, BT.Tnxy, GP.nxy, 0.0, 0.0, aPsi, 0.0, aM2Psi);

	for (int iConf=iConf0; iConf<=nConfFP; iConf++){
		// Plane wave ilumination
		IncidentWaveGPU.Psi0(xi, yi, Psi);
		// Move atoms
		PotentialGPU.Specimen.MoveAtoms(iConf);
		// Exit wave
		Cal_Wavefunction(Space, Psi);
		// Add Psi and M2aPsi to aPsi and aM2Psi
		f_AddwCwM2CtoCD(BT.Bnxy, BT.Tnxy, IsNotSharp, GP.nxy, inConfFP, Psi, aPsi, aM2Psi);
	}
	f_AddwM2CtoD(BT.Bnxy, BT.Tnxy, false, GP.nxy, 1.0, aPsi, M2aPsi);
}

/***************************************************************************/
// Set Input data
void cMulSliGPU::SetInputData(sInMSTEM &InMSTEM){
	freeMemory();

	MGP.gpu = InMSTEM.gpu;
	MGP.SimType = InMSTEM.SimType;	
	MGP.MulOrder = InMSTEM.MulOrder;	
	MGP.nConfFP = InMSTEM.nConfFP;		
	MGP.DimFP = InMSTEM.DimFP;	
	MGP.DistFP = 1;
	MGP.SeedFP = InMSTEM.SeedFP;
	MGP.PotPar = InMSTEM.PotPar;
	MGP.MEffect = InMSTEM.MEffect;
	MGP.STEffect = InMSTEM.STEffect;
	MGP.ZeroDefTyp = InMSTEM.ZeroDefTyp;
	MGP.ZeroDefPlane = InMSTEM.ZeroDefPlane;
	MGP.ApproxModel = InMSTEM.ApproxModel;
	MGP.Vrl = stVrl;
	MGP.E0 = InMSTEM.E0;	
	MGP.theta = InMSTEM.theta;	
	MGP.phi = InMSTEM.phi;
	MGP.lx = InMSTEM.lx;
	MGP.ly = InMSTEM.ly;
	MGP.dz = InMSTEM.dz;
	MGP.nx = InMSTEM.nx;
	MGP.ny = InMSTEM.ny;
	if(MGP.ApproxModel>1) MGP.MulOrder = 1;

	cudaSetDevice(MGP.gpu);

	f_sGP_Cal(MGP.nx, MGP.ny, MGP.lx, MGP.ly, MGP.dz, MGP.PBC_xy, GP);
	f_sBT_Cal(GP, BT);

	Lens.m = InMSTEM.MC_m;
	Lens.f = InMSTEM.MC_f;
	Lens.Cs3 = InMSTEM.MC_Cs3;
	Lens.Cs5 = InMSTEM.MC_Cs5;
	Lens.mfa2 = InMSTEM.MC_mfa2;
	Lens.afa2 = InMSTEM.MC_afa2;
	Lens.mfa3 = InMSTEM.MC_mfa3;
	Lens.afa3 = InMSTEM.MC_afa3;
	Lens.aobjl = InMSTEM.MC_aobjl;
	Lens.aobju = InMSTEM.MC_aobju;
	Lens.sf = InMSTEM.MC_sf;
	Lens.nsf = InMSTEM.MC_nsf;
	Lens.beta = InMSTEM.MC_beta;
	Lens.nbeta = InMSTEM.MC_nbeta;
	f_sLens_Cal(MGP.E0, GP, Lens);

	gxu = sin(MGP.theta)*cos(MGP.phi)/Lens.lambda;
	gyu = sin(MGP.theta)*sin(MGP.phi)/Lens.lambda;

	fPot = Lens.gamma*Lens.lambda/cos(MGP.theta);

	/*****************************************************************/
	STEM.SetInputData(InMSTEM);

	CBED.x0 = InMSTEM.CBED_x0;
	CBED.y0 = InMSTEM.CBED_y0;

	// HRTEM.xx = 0;

	PED.nrot = InMSTEM.PED_nrot;
	PED.theta = InMSTEM.PED_theta;

	HCI.nrot = InMSTEM.HCI_nrot;
	HCI.theta = InMSTEM.HCI_theta;
	// HCI.xx = 0;

	// EWRS.xx = 0;

	// EWFS.xx = 0;
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

	PotentialGPU.SetInputData(MGP, GP, BT, InMSTEM.nAtomsM, InMSTEM.AtomsM);

	// Set parameters for the incident wave
	IncidentWaveGPU.SetInputData(BT, GP, Lens, PlanPsi);

	// Microscope parameters
	MicroscopeEffectsGPU.SetInputData(BT, GP, Lens, PlanPsi, Trans);
}

void cMulSliGPU::Cal_ExitWaveRS(sComplex &aPsih, double *&aM2Psih){
	int MEffect = 0, STEffect = 0;
	eSpace Space = eSReal;
	get_Imagefuncion(MGP.nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	Gpu2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMulSliGPU::Cal_ExitWaveFS(sComplex &aPsih, double *&aM2Psih){
	int MEffect = 0, STEffect = 0;
	eSpace Space = eSReciprocal;
	get_Imagefuncion(MGP.nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	Gpu2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMulSliGPU::Cal_ED(sComplex &aPsih, double *&aM2Psih){
	int MEffect = 0, STEffect = 0;
	eSpace Space = eSReciprocal;
	get_Imagefuncion(MGP.nConfFP, Space, MEffect, STEffect, aPsi, M2aPsi, aM2Psi);
	Gpu2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMulSliGPU::Cal_HRTEM(sComplex &aPsih, double *&M2aPsih, double *&aM2Psih){
	eSpace Space = eSReal;
	get_Imagefuncion(MGP.nConfFP, Space, MGP.MEffect, MGP.STEffect, aPsi, M2aPsi, aM2Psi);
	Gpu2CPU(aPsi, M2aPsi, aM2Psi, aPsih, M2aPsih, aM2Psih);
}

void cMulSliGPU::Cal_CBED(sComplex &aPsih, double *&aM2Psih){
	eSpace Space = eSReciprocal;
	get_Imagefuncion(MGP.nConfFP, Space, CBED.x0, CBED.y0, aPsi, M2aPsi, aM2Psi);
	Gpu2CPU(aPsi, aM2Psi, aPsih, aM2Psih);
}

void cMulSliGPU::Cal_STEM(){
	int iThk=0;

	if(!STEM.FastCal){
		cDetectorGPU DetectorGPU;
		DetectorGPU.SetInputData(GP, BT, STEM.nDet, STEM.DetCir);
		for (int ist=0; ist<STEM.nst; ist++){
			get_Imagefuncion(MGP.nConfFP, eSReciprocal, STEM.xst[ist], STEM.yst[ist], aPsi, M2aPsi, aM2Psi);
			DetectorGPU.getDetectorIntensity(aM2Psi, M2aPsi, ist, STEM.ImSTEM[iThk].DetInt);
		}
	}else{
		switch(MGP.ApproxModel)
		{
			case 1:
				Cal_FAST_STEM_Wavefunction_MSA(MGP.nConfFP, STEM.ImSTEM[iThk].DetInt);
				break;
			case 2:
				Cal_FAST_STEM_Wavefunction_PA(MGP.nConfFP, STEM.ImSTEM[iThk].DetInt);
				break;
			case 3:
				Cal_FAST_STEM_Wavefunction_POA(MGP.nConfFP, STEM.ImSTEM[iThk].DetInt);
				break;
			case 4:
				Cal_FAST_STEM_Wavefunction_WPOA(MGP.nConfFP, STEM.ImSTEM[iThk].DetInt);
				break;
		}
	}
}

//void cMulSliGPU::Cal_STEM(){
//	int ist, iThk=0;
//	double *M2aPsih, *aM2Psih;
//	aM2Psih = new double[GP.nxy];
//	M2aPsih = new double[GP.nxy];
//	cDetectorCPU DetectorCPU;
//	DetectorCPU.SetInputData(GP, STEM.nDet, STEM.DetCir);
//
//	for (ist=0; ist<STEM.nst; ist++){
//		get_Imagefuncion(MGP.nConfFP, eSReciprocal, STEM.xs[ist], STEM.ys[ist], aPsi, M2aPsi, aM2Psi);
//		cudaMemcpy(aM2Psih, aM2Psi, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
//		cudaMemcpy(M2aPsih, M2aPsi, GP.nxy*cSizeofRD, cudaMemcpyDeviceToHost);
//		DetectorCPU.getDetectorIntensity(aM2Psih, M2aPsih, ils, STEM.ImSTEM[iThk].DetInt);
//	}
//
//	delete [] aM2Psih; aM2Psih = 0;
//	delete [] M2aPsih; M2aPsih = 0;
//}