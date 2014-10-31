#include <cmath>
#include <memory.h>
#include "hConstTypes.h"
#include "hgeneralCPU.h"
#include <mex.h>

void iScaleVector(int n, double *Vi, double f, double *Vo){
	for(int i=0; i<n; i++)
		Vo[i] = Vi[i]/f;
}

void CreateSetValue2mxField(mxArray *mxB, int p, const char *field_name, double field_value){
	mxArray *mxfield = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(mxfield) = field_value;
	mxSetField(mxB, p, field_name, mxfield);
}

void CreateSetValue2mxField(mxArray *mxB, int p, const char *field_name, int n, double *field_value){
	mxArray *mxfield = mxCreateDoubleMatrix(1, n, mxREAL);
	double *pfield = mxGetPr(mxfield); 
	memcpy(pfield, field_value, n*cSizeofRD);
	mxSetField(mxB, p, field_name, mxfield);
}

void CreateSetValue2mxField(mxArray *mxB, int p, const char *field_name, int nx, int ny, double *field_value){
	mxArray *mxfield = mxCreateDoubleMatrix(ny, nx, mxREAL);
	double *pfield = mxGetPr(mxfield); 
	memcpy(pfield, field_value, ny*nx*cSizeofRD);
	mxSetField(mxB, p, field_name, mxfield);
}

template <class Type>
Type ReadValuemxField(const mxArray *mxB, int p, const char *field_name, double factor=1){
	double val =  factor*mxGetScalar(mxGetField(mxB, p, field_name));
	return (Type)val;
}

void ReadValuemxField(const mxArray *mxB, int p, const char *field_name, int n, double *field_value, double factor=1){
	double *pfield = mxGetPr(mxGetField(mxB, p, field_name));
	memcpy(field_value, pfield, n*cSizeofRD);
	for(int i=0; i<n; i++)
		field_value[i] *= factor;
}

/*******************Matlab to layer unit cell*********************/
void Matlab2uLayer(const mxArray *mxCrystal, int &na, int &nb, int &nc, double &a, double &b, double &c, int &nuLayer, sAtomsGroup *&uLayer){
    int i, j, nAtoms;
	double *AtomsM;

    na = ReadValuemxField<int>(mxCrystal, 0, "na"); 
    nb = ReadValuemxField<int>(mxCrystal, 0, "nb");
    nc = ReadValuemxField<int>(mxCrystal, 0, "nc"); 

    a = ReadValuemxField<double>(mxCrystal, 0, "a");
    b = ReadValuemxField<double>(mxCrystal, 0, "b"); 
    c = ReadValuemxField<double>(mxCrystal, 0, "c");

    nuLayer = ReadValuemxField<int>(mxCrystal, 0, "nuLayer");

	mxArray *mexuLayer, *mexAtoms; 
	mexuLayer = mxGetField(mxCrystal, 0, "uLayer");
	
	uLayer = new sAtomsGroup[nuLayer];
	for (i=0; i<nuLayer; i++){
		mexAtoms = mxGetField(mexuLayer, i, "Atoms");
		AtomsM = mxGetPr(mexAtoms);
		uLayer[i].nAtoms = nAtoms = (int)mxGetM(mexAtoms);
		uLayer[i].Atoms = new sAtoms[nAtoms];
		for(j=0; j<nAtoms; j++){
			uLayer[i].Atoms[j].x = a*AtomsM[0*nAtoms+j];
			uLayer[i].Atoms[j].y = b*AtomsM[1*nAtoms+j];
			uLayer[i].Atoms[j].z = c*AtomsM[2*nAtoms+j];
			uLayer[i].Atoms[j].Z = (int)AtomsM[3*nAtoms+j];
			uLayer[i].Atoms[j].sigma = AtomsM[4*nAtoms+j];
			uLayer[i].Atoms[j].occ = AtomsM[5*nAtoms+j];
		}
	}
}

/***************Matlab to radial Schrodinger equation****************/
void Matlab2RadSchr(const mxArray *mxRadSchr, sInRadSchr &InRadSchr){
	InRadSchr.E0 = ReadValuemxField<double>(mxRadSchr, 0, "E0");
	InRadSchr.PotPar = ReadValuemxField<int>(mxRadSchr, 0, "PotPar");
	InRadSchr.n = ReadValuemxField<int>(mxRadSchr, 0, "n");
	InRadSchr.nr = ReadValuemxField<int>(mxRadSchr, 0, "nr");
	InRadSchr.nAtomsM = (int)mxGetM(mxGetField(mxRadSchr, 0, "Atoms"));
	InRadSchr.AtomsM = mxGetPr(mxGetField(mxRadSchr, 0, "Atoms"));
}

/****************************MUlTEM********************************/
// From AtomTypesCPU to Matlab structure 
void AtomTypesCPU2Matlab(int nAtomTypesCPU, sAtomTypesCPU *&AtomTypesCPU, mxArray *&mxAtomTypesCPU){
	const char *field_names[] = {"Z", "m", "A", "rn_e", "rn_c", "ra_e", "ra_c", "Rmin", "Rmax", "cfeg", "cfxg", "cPr", "cVr", "cVR", "ns", "Vo"};
	int number_of_fields = 16;
	mwSize dims[2] = {nAtomTypesCPU, 1};

	const char *field_names_Vpog[] = {"sigma", "Vr", "Vi", "gr", "gVr", "gVi"};
	int number_of_fields_Vpog = 6;
	mwSize dims_Vpog[2] = {1, 1};

	const char *field_names_Coef[] = {"cl", "cnl"};
	int number_of_fields_Coef = 2;
	mwSize dims_Coef[2] = {1, 1};

	mxArray *mxfield, *mxVpog;
	int i, j;

	mxAtomTypesCPU = mxCreateStructArray(2, dims, number_of_fields, field_names);
	for (i=0; i<nAtomTypesCPU; i++){
		CreateSetValue2mxField(mxAtomTypesCPU, i, "Z", AtomTypesCPU[i].Z);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "m", AtomTypesCPU[i].m);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "A", AtomTypesCPU[i].A);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "rn_e", AtomTypesCPU[i].rn_e);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "rn_c", AtomTypesCPU[i].rn_c);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "ra_e", AtomTypesCPU[i].ra_e);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "ra_c", AtomTypesCPU[i].ra_c);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "Rmin", AtomTypesCPU[i].Rmin);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "Rmax", AtomTypesCPU[i].Rmax);

		/**************************fg***************************/
		mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mxAtomTypesCPU, i, "cfeg", mxfield);
		CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cfeg.cl);
		CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cfeg.cnl);

		/**************************fx***************************/
		mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mxAtomTypesCPU, i, "cfxg", mxfield);
		CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cfxg.cl);
		CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cfxg.cnl);

		/**************************Pr***************************/
		mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mxAtomTypesCPU, i, "cPr", mxfield);
		CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cPr.cl);
		CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cPr.cnl);

		/**************************Vr***************************/
		mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mxAtomTypesCPU, i, "cVr", mxfield);
		CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cVr.cl);
		CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cVr.cnl);

		/**************************VR***************************/
		mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mxAtomTypesCPU, i, "cVR", mxfield);
		CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cVR.cl);
		CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cVR.cnl);

		CreateSetValue2mxField(mxAtomTypesCPU, i, "ns", AtomTypesCPU[i].ns);

		if (AtomTypesCPU[i].ns>0){
			dims_Vpog[1] = AtomTypesCPU[i].ns;
			mxVpog = mxCreateStructArray(2, dims_Vpog, number_of_fields_Vpog, field_names_Vpog);
			mxSetField(mxAtomTypesCPU, i, "Vo", mxVpog);
			for (j=0; j<AtomTypesCPU[i].ns; j++){
				CreateSetValue2mxField(mxVpog, j, "sigma", AtomTypesCPU[i].Vo[j].sigma);

				/**************************Vr***************************/
				mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
				mxSetField(mxVpog, j, "cVr", mxfield);
				CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Vo[j].cVr.cl);
				CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Vo[j].cVr.cnl);

				/**************************Vi***************************/
				mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
				mxSetField(mxVpog, j, "cVi", mxfield);
				CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Vo[j].cVi.cl);
				CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Vo[j].cVi.cnl);

				/*************************Grid**************************/
				CreateSetValue2mxField(mxVpog, j, "gr", stngbp, AtomTypesCPU[i].Vo[j].gr);
				CreateSetValue2mxField(mxVpog, j, "gVr", stngbp, AtomTypesCPU[i].Vo[j].gVr);
				CreateSetValue2mxField(mxVpog, j, "gVi", stngbp, AtomTypesCPU[i].Vo[j].gVi);
			}
		}
	}
}

// From AtomTypesCPU to Matlab structure 
void ImSTEM2Matlab(int nThk, int nDet, int line, int nxs, int nys, sImSTEM *ImSTEM, mxArray *&mxImSTEM){
	const char *field_names_ImSTEM[] = {"DetInt"};
	int number_of_fields_ImSTEM = 1;
	mwSize dims_ImSTEM[2] = {nThk, 1};

	const char *field_names_DetInt[] = {"Tot", "Coh"};
	int number_of_fields_DetInt = 2;
	mwSize dims_DetInt[2] = {nDet, 1};

	mxArray *mxDetInt;
	mxImSTEM = mxCreateStructArray(2, dims_ImSTEM, number_of_fields_ImSTEM, field_names_ImSTEM);
	for(int iThk=0; iThk<nThk; iThk++){
		mxDetInt = mxCreateStructArray(2, dims_DetInt, number_of_fields_DetInt, field_names_DetInt);
		mxSetField(mxImSTEM, iThk, "DetInt", mxDetInt);
		for(int iDet=0; iDet<nDet; iDet++){
			if(line==1){
				CreateSetValue2mxField(mxDetInt, iDet, "Tot", nxs, ImSTEM[iThk].DetInt[iDet].Tot);
				CreateSetValue2mxField(mxDetInt, iDet, "Coh", nxs, ImSTEM[iThk].DetInt[iDet].Coh);
			}else{
				CreateSetValue2mxField(mxDetInt, iDet, "Tot", nxs, nys, ImSTEM[iThk].DetInt[iDet].Tot);
				CreateSetValue2mxField(mxDetInt, iDet, "Coh", nxs, nys, ImSTEM[iThk].DetInt[iDet].Coh);
			}
		}
	}
}

// From Matlab structure to AtomTypesCPU
void Matlab2AtomTypesCPU(const mxArray *mxAtomTypesCPU, int &nAtomTypesCPU, sAtomTypesCPU *&AtomTypesCPU){
	mxArray *mxfield, *mxVpog;	
	int i, j, ns;

	nAtomTypesCPU = (int)mxGetM(mxAtomTypesCPU)*(int)mxGetN(mxAtomTypesCPU);
	delete [] AtomTypesCPU; AtomTypesCPU = new sAtomTypesCPU[nAtomTypesCPU]; 
	for (i=0; i<nAtomTypesCPU; i++){
		AtomTypesCPU[i].Z = ReadValuemxField<int>(mxAtomTypesCPU, i, "Z");
		AtomTypesCPU[i].m = ReadValuemxField<double>(mxAtomTypesCPU, i, "m");
		AtomTypesCPU[i].A = ReadValuemxField<int>(mxAtomTypesCPU, i, "A");
		AtomTypesCPU[i].rn_e = ReadValuemxField<double>(mxAtomTypesCPU, i, "rn_e");
		AtomTypesCPU[i].rn_c = ReadValuemxField<double>(mxAtomTypesCPU, i, "rn_c");
		AtomTypesCPU[i].ra_e = ReadValuemxField<double>(mxAtomTypesCPU, i, "ra_e");
		AtomTypesCPU[i].ra_c = ReadValuemxField<double>(mxAtomTypesCPU, i, "ra_c");
		AtomTypesCPU[i].Rmin = ReadValuemxField<double>(mxAtomTypesCPU, i, "Rmin");
		AtomTypesCPU[i].Rmax = ReadValuemxField<double>(mxAtomTypesCPU, i, "Rmax");

		/*******************************fg***********************************/
		mxfield = mxGetField(mxAtomTypesCPU, i, "cfeg");
		ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cfeg.cl);
		ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cfeg.cnl);

		/*******************************fx***********************************/
		mxfield = mxGetField(mxAtomTypesCPU, i, "cfxg");
		ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cfxg.cl);
		ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cfxg.cnl);

		/*******************************Pr***********************************/
		mxfield = mxGetField(mxAtomTypesCPU, i, "cPr");
		ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cPr.cl);
		ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cPr.cnl);

		/*******************************Vr***********************************/
		mxfield = mxGetField(mxAtomTypesCPU, i, "cVr");
		ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cVr.cl);
		ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cVr.cnl);

		/*******************************VR***********************************/
		mxfield = mxGetField(mxAtomTypesCPU, i, "cVR");
		ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].cVR.cl);
		ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].cVR.cnl);

		mxVpog= mxGetField(mxAtomTypesCPU, i, "Vo");
		ns = (int)mxGetM(mxVpog)*(int)mxGetN(mxVpog);

		if(ns>0){
			AtomTypesCPU[i].Vo = new sVoCPU[ns];
			for (j=0; j<ns; j++){
				AtomTypesCPU[i].Vo[j].sigma = ReadValuemxField<double>(mxVpog, j, "sigma");

				/*******************************Vr**********************************/
				mxfield = mxGetField(mxVpog, j, "cVr");
				ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Vo[j].cVr.cl);
				ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Vo[j].cVr.cnl);

				/*******************************Vi**********************************/
				mxfield = mxGetField(mxVpog, j, "cVi");
				ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Vo[j].cVi.cl);
				ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Vo[j].cVi.cnl);

				/******************************Grid*********************************/
				ReadValuemxField(mxVpog, j, "gr", stngbp, AtomTypesCPU[i].Vo[j].gr);
				ReadValuemxField(mxVpog, j, "gVr", stngbp, AtomTypesCPU[i].Vo[j].gVr);
				ReadValuemxField(mxVpog, j, "gVi", stngbp, AtomTypesCPU[i].Vo[j].gVi);
			}
		}
	}
}

/***********************read input TEMim*************************/
void Matlab2InTEMIm(const mxArray *mxInTEMIm, sInTEMIm &InTEMIm){
	InTEMIm.gpu = ReadValuemxField<int>(mxInTEMIm, 0, "gpu");
	InTEMIm.MEffect = ReadValuemxField<int>(mxInTEMIm, 0, "MEffect");				// 1: Exit wave Partial coherente mode, 2: Transmission cross coefficient
	InTEMIm.STEffect = ReadValuemxField<int>(mxInTEMIm, 0, "STEffect");				// 1: Spatial and temporal, 2: Temporal, 3: Spatial
	InTEMIm.ZeroDef = ReadValuemxField<int>(mxInTEMIm, 0, "ZeroDefTyp");			// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
	InTEMIm.ZeroDefPlane = ReadValuemxField<double>(mxInTEMIm, 0, "ZeroDefPlane");	// Zero defocus plane	

	InTEMIm.E0 = ReadValuemxField<double>(mxInTEMIm, 0, "E0");			// kV
	mxArray *mxPsi = mxGetField(mxInTEMIm, 0, "Psi");					// Exit wave
	InTEMIm.Psirh = mxGetPr(mxPsi);										// real (Exit wave)
	InTEMIm.Psiih = mxGetPi(mxPsi);										// imaginary (Exit wave)
	InTEMIm.lx = ReadValuemxField<double>(mxInTEMIm, 0, "lx");			// distance in x direction(Angstroms)
	InTEMIm.ly = ReadValuemxField<double>(mxInTEMIm, 0, "ly");			// distance in y direction(Angstroms)
	InTEMIm.ny = (int)mxGetM(mxPsi);									// Number of pixels in x direction
	InTEMIm.nx = (int)mxGetN(mxPsi);									// Number of pixels in y direction

	// Microscope parameters
	InTEMIm.MC_m = ReadValuemxField<int>(mxInTEMIm, 0, "m");						// momentum of the vortex
	InTEMIm.MC_f = ReadValuemxField<double>(mxInTEMIm, 0, "f");						// defocus(Angstrom)
	InTEMIm.MC_Cs3 = ReadValuemxField<double>(mxInTEMIm, 0, "Cs3", mm2Ags);			// spherical aberration(mm-->Angstrom)
	InTEMIm.MC_Cs5 = ReadValuemxField<double>(mxInTEMIm, 0, "Cs5", mm2Ags);			// spherical aberration(mm-->Angstrom)
	InTEMIm.MC_mfa2 = ReadValuemxField<double>(mxInTEMIm, 0, "mfa2");				// magnitude 2-fold astigmatism(Angstrom)
	InTEMIm.MC_afa2 = ReadValuemxField<double>(mxInTEMIm, 0, "afa2", deg2rad);		// angle 2-fold astigmatism(degrees-->rad)
	InTEMIm.MC_mfa3 = ReadValuemxField<double>(mxInTEMIm, 0, "mfa3");				// magnitude 3-fold astigmatism(Angstrom)
	InTEMIm.MC_afa3 = ReadValuemxField<double>(mxInTEMIm, 0, "afa3", deg2rad);		// angle 3-fold astigmatism(degrees-->rad)
	InTEMIm.MC_aobjl = ReadValuemxField<double>(mxInTEMIm, 0, "aobjl", mrad2rad);	// lower objective aperture(mrad-->rad)
	InTEMIm.MC_aobju = ReadValuemxField<double>(mxInTEMIm, 0, "aobju", mrad2rad);	// upper objective aperture(mrad-->rad)
	InTEMIm.MC_sf = ReadValuemxField<double>(mxInTEMIm, 0, "sf");					// defocus spread(Angstrom)
	InTEMIm.MC_nsf = ReadValuemxField<int>(mxInTEMIm, 0, "nsf");						// Number of defocus sampling point
	InTEMIm.MC_beta = ReadValuemxField<double>(mxInTEMIm, 0, "beta", mrad2rad);		// semi-convergence angle(mrad-->rad)
	InTEMIm.MC_nbeta = ReadValuemxField<int>(mxInTEMIm, 0, "nbeta");					// half number sampling points
 }

/***********************read input MulSli*************************/
void Matlab2InMulSli(const mxArray *mxInMSTEM, sInMSTEM &InMSTEM){
	InMSTEM.gpu = ReadValuemxField<int>(mxInMSTEM, 0, "gpu");						// gpu device
	InMSTEM.SimType = ReadValuemxField<int>(mxInMSTEM, 0, "SimType");				// 1: STEM, 2: HRTEM, 3: ED, 4: PED, 5: CBED, 6: HCI, ... 10: EW real, 11: EW Fourier
	InMSTEM.MulOrder = ReadValuemxField<int>(mxInMSTEM, 0, "MulOrder");				// 1: First order MS, 2: Second Order MS
	InMSTEM.nConfFP = ReadValuemxField<int>(mxInMSTEM, 0, "nConfFP");				// Number of frozen phonon configurations
	InMSTEM.DimFP = ReadValuemxField<int>(mxInMSTEM, 0, "DimFP");					// Dimensions phonon configurations
	InMSTEM.SeedFP = ReadValuemxField<int>(mxInMSTEM, 0, "SeedFP");					// Random seed(frozen phonon)
	InMSTEM.PotPar = ReadValuemxField<int>(mxInMSTEM, 0, "PotPar");					// Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
	InMSTEM.MEffect = ReadValuemxField<int>(mxInMSTEM, 0, "MEffect");				// 1: Partial coherente mode, 2: Transmission cross coefficient
	InMSTEM.STEffect = ReadValuemxField<int>(mxInMSTEM, 0, "STEffect");				// 1: Spatial and temporal, 2: Temporal, 3: Spatial
	InMSTEM.ZeroDefTyp = ReadValuemxField<int>(mxInMSTEM, 0, "ZeroDefTyp");			// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
	InMSTEM.ZeroDefPlane = ReadValuemxField<double>(mxInMSTEM, 0, "ZeroDefPlane");	// Zero defocus plane
	InMSTEM.ApproxModel = ReadValuemxField<int>(mxInMSTEM, 0, "ApproxModel");		// 1: MS, 2: PA, 3POA, 4:WPOA
	
	/***************************Multislice*************************/
	InMSTEM.E0 = ReadValuemxField<double>(mxInMSTEM, 0, "E0");					// Acceleration voltage
	InMSTEM.theta = ReadValuemxField<double>(mxInMSTEM, 0, "theta", deg2rad);	// incident tilt (in spherical coordinates) (degrees-->rad)
	InMSTEM.phi = ReadValuemxField<double>(mxInMSTEM, 0, "phi", deg2rad);		// incident tilt (in spherical coordinates) (degrees-->rad)
	InMSTEM.nx = ReadValuemxField<int>(mxInMSTEM, 0, "nx");						// Number of pixels in x direction
	InMSTEM.ny = ReadValuemxField<int>(mxInMSTEM, 0, "ny");						// Number of pixels in y direction
	InMSTEM.lx = ReadValuemxField<double>(mxInMSTEM, 0, "lx");					// distance in x direction(Angstroms)
	InMSTEM.ly = ReadValuemxField<double>(mxInMSTEM, 0, "ly");					// distance in y direction(Angstroms)
	InMSTEM.dz = ReadValuemxField<double>(mxInMSTEM, 0, "dz");					// slice thickness

	/*********************Microscope parameters********************/
	mxArray *mxMC= mxGetField(mxInMSTEM, 0, "MC");
	InMSTEM.MC_m =  ReadValuemxField<int>(mxMC, 0, "m");						// momentum of the vortex
	InMSTEM.MC_f = ReadValuemxField<double>(mxMC, 0, "f");						// defocus(Angstrom)
	InMSTEM.MC_Cs3 = ReadValuemxField<double>(mxMC, 0, "Cs3", mm2Ags);			// spherical aberration(mm-->Angstrom)
	InMSTEM.MC_Cs5 = ReadValuemxField<double>(mxMC, 0, "Cs5", mm2Ags);			// spherical aberration(mm-->Angstrom)
	InMSTEM.MC_mfa2 = ReadValuemxField<double>(mxMC, 0, "mfa2");				// magnitude 2-fold astigmatism(Angstrom)
	InMSTEM.MC_afa2 = ReadValuemxField<double>(mxMC, 0, "afa2", deg2rad);		// angle 2-fold astigmatism(degrees-->rad)
	InMSTEM.MC_mfa3 = ReadValuemxField<double>(mxMC, 0, "mfa3");				// magnitude 3-fold astigmatism(Angstrom)
	InMSTEM.MC_afa3 = ReadValuemxField<double>(mxMC, 0, "afa3", deg2rad);		// angle 3-fold astigmatism(degrees-->rad)
	InMSTEM.MC_aobjl = ReadValuemxField<double>(mxMC, 0, "aobjl", mrad2rad);	// lower objective aperture(mrad-->rad)
	InMSTEM.MC_aobju = ReadValuemxField<double>(mxMC, 0, "aobju", mrad2rad);	// upper objective aperture(mrad-->rad)
	InMSTEM.MC_sf = ReadValuemxField<double>(mxMC, 0, "sf");					// defocus spread(Angstrom)
	InMSTEM.MC_nsf = ReadValuemxField<int>(mxMC, 0, "nsf");						// Number of defocus sampling point
	InMSTEM.MC_beta = ReadValuemxField<double>(mxMC, 0, "beta", mrad2rad);		// semi-convergence angle(mrad-->rad)
	InMSTEM.MC_nbeta = ReadValuemxField<int>(mxMC, 0, "nbeta");					// half number sampling points

	mxArray *mxAtomsM = mxGetField(mxInMSTEM, 0, "Atoms");
	InMSTEM.nAtomsM = (int)mxGetM(mxAtomsM);									// Number of Atoms
	InMSTEM.AtomsM = mxGetPr(mxAtomsM);											// Atoms in a matrix form

	switch (InMSTEM.SimType){
		case 1:	// STEM
			mxArray *mxSTEM;	
			mxSTEM = mxGetField(mxInMSTEM, 0, "STEM");
			InMSTEM.STEM_line = ReadValuemxField<int>(mxSTEM, 0, "line");
			InMSTEM.STEM_FastCal = ReadValuemxField<int>(mxSTEM, 0, "FastCal");
			InMSTEM.STEM_ns = ReadValuemxField<int>(mxSTEM, 0, "ns");
			InMSTEM.STEM_x1u = ReadValuemxField<double>(mxSTEM, 0, "x1u");
			InMSTEM.STEM_y1u = ReadValuemxField<double>(mxSTEM, 0, "y1u");
			InMSTEM.STEM_x2u = ReadValuemxField<double>(mxSTEM, 0, "x2u");
			InMSTEM.STEM_y2u = ReadValuemxField<double>(mxSTEM, 0, "y2u");

			InMSTEM.STEM_nDet = ReadValuemxField<int>(mxSTEM, 0, "nDet");
			if(InMSTEM.STEM_nDet>0){
				InMSTEM.STEM_DetCir = new sInDetCir[InMSTEM.STEM_nDet];
				mxArray *mxDetCir;
				mxDetCir = mxGetField(mxSTEM, 0, "DetCir");
				for (int i=0; i<InMSTEM.STEM_nDet; i++){
					InMSTEM.STEM_DetCir[i].InnerAng = ReadValuemxField<double>(mxDetCir, i, "InnerAng", mrad2rad);	// Inner angle(mrad-->rad)
					InMSTEM.STEM_DetCir[i].OuterAng = ReadValuemxField<double>(mxDetCir, i, "OuterAng", mrad2rad);	// Outer angle(mrad-->rad)
				}
			}
			break;
		case 2:		// CBED
			mxArray *mxCBED;
			mxCBED = mxGetField(mxInMSTEM, 0, "CBED");
			InMSTEM.CBED_x0 = ReadValuemxField<double>(mxCBED, 0, "x0");	// 
			InMSTEM.CBED_y0 = ReadValuemxField<double>(mxCBED, 0, "y0");	//
			break;
		case 3:		// HRTEM
			mxArray *mxHRTEM;
			mxHRTEM = mxGetField(mxInMSTEM, 0, "HRTEM");
			//InMSTEM.HRTEM_xx

			break;
		case 4:		// ED
			mxArray *mxED;
			mxED = mxGetField(mxInMSTEM, 0, "ED");

			break;
		case 5:		// PED
			mxArray *mxPED;
			mxPED = mxGetField(mxInMSTEM, 0, "PED");
			InMSTEM.PED_nrot = ReadValuemxField<int>(mxPED, 0, "nrot");					// Number of orientations
			InMSTEM.PED_theta = ReadValuemxField<double>(mxPED, 0, "theta", deg2rad);	// Precession angle(degrees-->rad)
			break;
		case 6:		// HCI
			mxArray *mxHCI;
			mxHCI = mxGetField(mxInMSTEM, 0, "HCI");
			InMSTEM.HCI_nrot = ReadValuemxField<int>(mxHCI, 0, "nrot");					// Number of orientations
			InMSTEM.HCI_theta = ReadValuemxField<double>(mxHCI, 0, "theta", deg2rad);	// Precession angle(degrees-->rad)
			//InMSTEM.HCI_xx;

			break;
		case 10:	// EW real
			mxArray *mxEWRS;
			mxEWRS = mxGetField(mxInMSTEM, 0, "EWRS");
			//InMSTEM.EWRS_xx
			break;
		case 11:	// EW Fourier
			mxArray *mxEWFS;
			mxEWFS = mxGetField(mxInMSTEM, 0, "EWFS");
			//InMSTEM.EWFS_xx
			break;
	}
 }

/***********************read input Probe*************************/
void Matlab2InProbe(const mxArray *mxInProbe, sInProbe &InProbe){
	InProbe.gpu = ReadValuemxField<int>(mxInProbe, 0, "gpu");						// gpu device
	
	/***************************Multislice*************************/
	InProbe.E0 = ReadValuemxField<double>(mxInProbe, 0, "E0");						// Acceleration voltage
	InProbe.theta = ReadValuemxField<double>(mxInProbe, 0, "theta", deg2rad);		// incident tilt (in spherical coordinates) (degrees-->rad)
	InProbe.phi = ReadValuemxField<double>(mxInProbe, 0, "phi", deg2rad);			// incident tilt (in spherical coordinates) (degrees-->rad)
	InProbe.nx = ReadValuemxField<int>(mxInProbe, 0, "nx");							// Number of pixels in x direction
	InProbe.ny = ReadValuemxField<int>(mxInProbe, 0, "ny");							// Number of pixels in y direction
	InProbe.lx = ReadValuemxField<double>(mxInProbe, 0, "lx");						// distance in x direction(Angstroms)
	InProbe.ly = ReadValuemxField<double>(mxInProbe, 0, "ly");						// distance in y direction(Angstroms)

	InProbe.x0 = ReadValuemxField<double>(mxInProbe, 0, "x0");	// 
	InProbe.y0 = ReadValuemxField<double>(mxInProbe, 0, "y0");	//

	InProbe.m =  ReadValuemxField<int>(mxInProbe, 0, "m");						// momentum of the vortex
	InProbe.f = ReadValuemxField<double>(mxInProbe, 0, "f");						// defocus(Angstrom)
	InProbe.Cs3 = ReadValuemxField<double>(mxInProbe, 0, "Cs3", mm2Ags);			// spherical aberration(mm-->Angstrom)
	InProbe.Cs5 = ReadValuemxField<double>(mxInProbe, 0, "Cs5", mm2Ags);			// spherical aberration(mm-->Angstrom)
	InProbe.mfa2 = ReadValuemxField<double>(mxInProbe, 0, "mfa2");					// magnitude 2-fold astigmatism(Angstrom)
	InProbe.afa2 = ReadValuemxField<double>(mxInProbe, 0, "afa2", deg2rad);			// angle 2-fold astigmatism(degrees-->rad)
	InProbe.mfa3 = ReadValuemxField<double>(mxInProbe, 0, "mfa3");					// magnitude 3-fold astigmatism(Angstrom)
	InProbe.afa3 = ReadValuemxField<double>(mxInProbe, 0, "afa3", deg2rad);			// angle 3-fold astigmatism(degrees-->rad)
	InProbe.aobjl = ReadValuemxField<double>(mxInProbe, 0, "aobjl", mrad2rad);		// lower objective aperture(mrad-->rad)
	InProbe.aobju = ReadValuemxField<double>(mxInProbe, 0, "aobju", mrad2rad);		// upper objective aperture(mrad-->rad)
	InProbe.sf = ReadValuemxField<double>(mxInProbe, 0, "sf");						// defocus spread(Angstrom)
	InProbe.nsf = ReadValuemxField<int>(mxInProbe, 0, "nsf");						// Number of defocus sampling point
	InProbe.beta = ReadValuemxField<double>(mxInProbe, 0, "beta", mrad2rad);		// semi-convergence angle(mrad-->rad)
	InProbe.nbeta = ReadValuemxField<int>(mxInProbe, 0, "nbeta");					// half number sampling points
 }