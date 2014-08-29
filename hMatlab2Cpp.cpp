#include "math.h"
#include <memory.h>
#include "hConstTypes.h"
#include "mex.h"

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

/***************Matlab to radial Schrodinger equation****************/
void Matlab2RadSchr(const mxArray *mxRadSchr, sRadSchr &RadSchr){
	RadSchr.E0 = ReadValuemxField<double>(mxRadSchr, 0, "E0");
	RadSchr.PotPar = ReadValuemxField<int>(mxRadSchr, 0, "PotPar");
	RadSchr.n = ReadValuemxField<int>(mxRadSchr, 0, "n");
	RadSchr.nr = ReadValuemxField<int>(mxRadSchr, 0, "nr");
	RadSchr.nAtoms = (int)mxGetM(mxGetField(mxRadSchr, 0, "Atoms"));
	RadSchr.AtomsM = mxGetPr(mxGetField(mxRadSchr, 0, "Atoms"));
}

/****************************MUlTEM********************************/
// From AtomTypesCPU to Matlab structure 
void AtomTypesCPU2Matlab(int nAtomTypesCPU, sAtomTypesCPU *&AtomTypesCPU, mxArray *&mxAtomTypesCPU){
	const char *field_names[] = {"Z", "m", "A", "rn_e", "rn_c", "ra_e", "ra_c", "Rmin", "Rmax", "PotPar", "occ", "feg", "fxg", "Pr", "Vr", "VR", "ns", "Vo"};
	int number_of_fields = 18;
	int dims[2] = {nAtomTypesCPU, 1};

	const char *field_names_Vpog[] = {"sigma", "Vr", "Vi", "gr", "gVr", "gVi"};
	int number_of_fields_Vpog = 6;
	int dims_Vpog[2] = {6, 1};

	const char *field_names_Coef[] = {"cl", "cnl"};
	int number_of_fields_Coef = 2;
	int dims_Coef[2] = {1, 1};

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
		CreateSetValue2mxField(mxAtomTypesCPU, i, "PotPar", AtomTypesCPU[i].PotPar);
		CreateSetValue2mxField(mxAtomTypesCPU, i, "occ", AtomTypesCPU[i].occ);

		/**************************fg***************************/
		mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mxAtomTypesCPU, i, "feg", mxfield);
		CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].feg.cl);
		CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].feg.cnl);

		/**************************fx***************************/
		mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mxAtomTypesCPU, i, "fxg", mxfield);
		CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].fxg.cl);
		CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].fxg.cnl);

		/**************************Pr***************************/
		mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mxAtomTypesCPU, i, "Pr", mxfield);
		CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Pr.cl);
		CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Pr.cnl);

		/**************************Vr***************************/
		mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mxAtomTypesCPU, i, "Vr", mxfield);
		CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Vr.cl);
		CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Vr.cnl);

		/**************************VR***************************/
		mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
		mxSetField(mxAtomTypesCPU, i, "VR", mxfield);
		CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].VR.cl);
		CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].VR.cnl);

		CreateSetValue2mxField(mxAtomTypesCPU, i, "ns", AtomTypesCPU[i].ns);
		mxVpog = mxCreateStructArray(2, dims_Vpog, number_of_fields_Vpog, field_names_Vpog);
		mxSetField(mxAtomTypesCPU, i, "Vo", mxVpog);

		if (AtomTypesCPU[i].ns>0){
			for (j=0; j<AtomTypesCPU[i].ns; j++){
				CreateSetValue2mxField(mxVpog, j, "sigma", AtomTypesCPU[i].Vo[j].sigma);

				/**************************Vr***************************/
				mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
				mxSetField(mxVpog, j, "Vr", mxfield);
				CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Vo[j].Vr.cl);
				CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Vo[j].Vr.cnl);

				/**************************Vi***************************/
				mxfield = mxCreateStructArray(2, dims_Coef, number_of_fields_Coef, field_names_Coef);
				mxSetField(mxVpog, j, "Vi", mxfield);
				CreateSetValue2mxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Vo[j].Vi.cl);
				CreateSetValue2mxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Vo[j].Vi.cnl);

				/*************************Grid**************************/
				CreateSetValue2mxField(mxVpog, j, "gr", stngbp, AtomTypesCPU[i].Vo[j].gr);
				CreateSetValue2mxField(mxVpog, j, "gVr", stngbp, AtomTypesCPU[i].Vo[j].gVr);
				CreateSetValue2mxField(mxVpog, j, "gVi", stngbp, AtomTypesCPU[i].Vo[j].gVi);
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
		AtomTypesCPU[i].PotPar = ReadValuemxField<int>(mxAtomTypesCPU, i, "PotPar");
		AtomTypesCPU[i].occ = ReadValuemxField<double>(mxAtomTypesCPU, i, "occ");

		/*******************************fg***********************************/
		mxfield = mxGetField(mxAtomTypesCPU, i, "feg");
		ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].feg.cl);
		ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].feg.cnl);

		/*******************************fx***********************************/
		mxfield = mxGetField(mxAtomTypesCPU, i, "fxg");
		ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].fxg.cl);
		ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].fxg.cnl);

		/*******************************Pr***********************************/
		mxfield = mxGetField(mxAtomTypesCPU, i, "Pr");
		ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Pr.cl);
		ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Pr.cnl);

		/*******************************Vr***********************************/
		mxfield = mxGetField(mxAtomTypesCPU, i, "Vr");
		ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Vr.cl);
		ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Vr.cnl);

		/*******************************VR***********************************/
		mxfield = mxGetField(mxAtomTypesCPU, i, "VR");
		ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].VR.cl);
		ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].VR.cnl);

		mxVpog= mxGetField(mxAtomTypesCPU, i, "Vo");
		ns = (int)mxGetM(mxVpog)*(int)mxGetN(mxVpog);

		if(ns>0){
			AtomTypesCPU[i].Vo = new sVoCPU[ns];
			for (j=0; j<ns; j++){
				AtomTypesCPU[i].Vo[j].sigma = ReadValuemxField<double>(mxVpog, j, "sigma");

				/*******************************Vr**********************************/
				mxfield = mxGetField(mxVpog, j, "Vr");
				ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Vo[j].Vr.cl);
				ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Vo[j].Vr.cnl);

				/*******************************Vi**********************************/
				mxfield = mxGetField(mxVpog, j, "Vi");
				ReadValuemxField(mxfield, 0, "cl", 6, AtomTypesCPU[i].Vo[j].Vi.cl);
				ReadValuemxField(mxfield, 0, "cnl", 6, AtomTypesCPU[i].Vo[j].Vi.cnl);

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

/***************************Init MulSli**************************/
void InMulSli_Init(sInMSTEM &InMSTEM){
	InMSTEM.gpu = 0;
	InMSTEM.SimType = 0;
	InMSTEM.MulOrder = 0;
	InMSTEM.nConfFP = 0;
	InMSTEM.DimFP = 0;
	InMSTEM.SeedFP = 0;
	InMSTEM.PotPar = 0;

	InMSTEM.E0 = 0;
	InMSTEM.theta = 0;
	InMSTEM.phi = 0;
	InMSTEM.nx = 0;
	InMSTEM.ny = 0;
	InMSTEM.lx = 0;
	InMSTEM.ly = 0;
	InMSTEM.dz = 0;

	InMSTEM.MC_m = 0;
	InMSTEM.MC_f = 0;
	InMSTEM.MC_Cs3 = 0;
	InMSTEM.MC_Cs5 = 0;
	InMSTEM.MC_mfa2 = 0;
	InMSTEM.MC_afa2 = 0;
	InMSTEM.MC_mfa3 = 0;
	InMSTEM.MC_afa3 = 0;
	InMSTEM.MC_aobjl = 0;
	InMSTEM.MC_aobju = 0;
	InMSTEM.MC_sf = 0;
	InMSTEM.MC_nsf = 0;
	InMSTEM.MC_beta = 0;
	InMSTEM.MC_nbeta = 0;

	InMSTEM.nAtomsM = 0;
	InMSTEM.AtomsM = 0;

	InMSTEM.STEM_line = 0;
	InMSTEM.STEM_GridType = 0;
	InMSTEM.STEM_ns = 0;
	InMSTEM.STEM_nucx = 0;
	InMSTEM.STEM_nucy = 0;
	InMSTEM.STEM_x1u = 0;
	InMSTEM.STEM_y1u = 0;
	InMSTEM.STEM_x2u = 0;
	InMSTEM.STEM_y2u = 0;
	
	InMSTEM.STEM_nDetCir = 0;
	InMSTEM.STEM_DetCir = 0;

	InMSTEM.CBED_x0 = 0;
	InMSTEM.CBED_y0 = 0;

	InMSTEM.HRTEM_MEffect = 0;
	InMSTEM.HRTEM_STEffect = 0;
	InMSTEM.HRTEM_ZeroDefTyp = 0;
	InMSTEM.HRTEM_ZeroDefPlane = 0;

	InMSTEM.PED_nrot = 0;
	InMSTEM.PED_theta = 0;

	InMSTEM.HCI_nrot = 0;
	InMSTEM.HCI_theta = 0;
	InMSTEM.HCI_MEffect = 0;
	InMSTEM.HCI_STEffect = 0;
	InMSTEM.HCI_ZeroDefTyp = 0;
	InMSTEM.HCI_ZeroDefPlane = 0;

	InMSTEM.EWRS_ZeroDefTyp = 0;
	InMSTEM.EWRS_ZeroDefPlane = 0;

	InMSTEM.EWFS_ZeroDefTyp = 0;
	InMSTEM.EWFS_ZeroDefPlane = 0;
}

/***********************read input MulSli*************************/
void Matlab2InMulSli(const mxArray *mxInMSTEM, sInMSTEM &InMSTEM){
	InMulSli_Init(InMSTEM);

	InMSTEM.gpu = ReadValuemxField<int>(mxInMSTEM, 0, "gpu");				// gpu device
	InMSTEM.SimType = ReadValuemxField<int>(mxInMSTEM, 0, "SimType");		// 1: STEM, 2: HRTEM, 3: ED, 4: PED, 5: CBED, 6: HCI, ... 10: EW real, 11: EW Fourier
	InMSTEM.MulOrder = ReadValuemxField<int>(mxInMSTEM, 0, "MulOrder");		// 1: First order MS, 2: Second Order MS
	InMSTEM.nConfFP = ReadValuemxField<int>(mxInMSTEM, 0, "nConfFP");		// Number of frozen phonon configurations
	InMSTEM.DimFP = ReadValuemxField<int>(mxInMSTEM, 0, "DimFP");			// Dimensions phonon configurations
	InMSTEM.SeedFP = ReadValuemxField<int>(mxInMSTEM, 0, "SeedFP");			// Random seed(frozen phonon)
	InMSTEM.PotPar = ReadValuemxField<int>(mxInMSTEM, 0, "PotPar");			// Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)

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
	InMSTEM.MC_nbeta = ReadValuemxField<int>(mxMC, 0, "nbeta");						// half number sampling points

	mxArray *mxAtomsM = mxGetField(mxInMSTEM, 0, "Atoms");
	InMSTEM.nAtomsM = (int)mxGetM(mxAtomsM);	// Number of Atoms
	InMSTEM.AtomsM = mxGetPr(mxAtomsM);			// Atoms in a matrix form

	switch (InMSTEM.SimType){	
		case 1:	// STEM
			mxArray *mxSTEM;	
			mxSTEM = mxGetField(mxInMSTEM, 0, "STEM");
			InMSTEM.STEM_line = ReadValuemxField<int>(mxSTEM, 0, "line");
			InMSTEM.STEM_GridType = ReadValuemxField<int>(mxSTEM, 0, "GridType");
			InMSTEM.STEM_ns = ReadValuemxField<int>(mxSTEM, 0, "ns");
			InMSTEM.STEM_nucx = ReadValuemxField<int>(mxSTEM, 0, "nucx");
			InMSTEM.STEM_nucy = ReadValuemxField<int>(mxSTEM, 0, "nucy");
			InMSTEM.STEM_x1u = ReadValuemxField<double>(mxSTEM, 0, "x1u");
			InMSTEM.STEM_y1u = ReadValuemxField<double>(mxSTEM, 0, "y1u");
			InMSTEM.STEM_x2u = ReadValuemxField<double>(mxSTEM, 0, "x2u");
			InMSTEM.STEM_y2u = ReadValuemxField<double>(mxSTEM, 0, "y2u");

			InMSTEM.STEM_nDetCir = ReadValuemxField<int>(mxSTEM, 0, "nDetCir");
			InMSTEM.STEM_DetCir = new sInDetCir[InMSTEM.STEM_nDetCir];
			mxArray *mxDetCir;
			mxDetCir = mxGetField(mxSTEM, 0, "DetCir");
			for (int i=0; i<InMSTEM.STEM_nDetCir; i++){
				InMSTEM.STEM_DetCir[i].InnerAng = ReadValuemxField<double>(mxDetCir, i, "InnerAng", mrad2rad);	// Inner angle(mrad-->rad)
				InMSTEM.STEM_DetCir[i].OuterAng = ReadValuemxField<double>(mxDetCir, i, "OuterAng", mrad2rad);	// Outer angle(mrad-->rad)
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
			InMSTEM.HRTEM_MEffect = ReadValuemxField<int>(mxHRTEM, 0, "MEffect");					// 1: Partial coherente mode, 2: Transmission cross coefficient
			InMSTEM.HRTEM_STEffect = ReadValuemxField<int>(mxHRTEM, 0, "STEffect");					// 1: Spatial and temporal, 2: Temporal, 3: Spatial
			InMSTEM.HRTEM_ZeroDefTyp = ReadValuemxField<int>(mxHRTEM, 0, "ZeroDefTyp");				// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
			InMSTEM.HRTEM_ZeroDefPlane = ReadValuemxField<double>(mxHRTEM, 0, "ZeroDefPlane");		// Zero defocus plane
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
			InMSTEM.HCI_nrot = ReadValuemxField<int>(mxHCI, 0, "nrot");							// Number of orientations
			InMSTEM.HCI_theta = ReadValuemxField<double>(mxHCI, 0, "theta", deg2rad);			// Precession angle(degrees-->rad)
			InMSTEM.HCI_MEffect = ReadValuemxField<int>(mxHCI, 0, "MEffect");					// 1: Partial coherente mode, 2: Transmission cross coefficient
			InMSTEM.HCI_STEffect = ReadValuemxField<int>(mxHCI, 0, "STEffect");					// 1: Spatial and temporal, 2: Temporal, 3: Spatial
			InMSTEM.HCI_ZeroDefTyp = ReadValuemxField<int>(mxHCI, 0, "ZeroDefTyp");				// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
			InMSTEM.HCI_ZeroDefPlane = ReadValuemxField<double>(mxHCI, 0, "ZeroDefPlane");		// Zero defocus plane
			break;
		case 10:	// EW real
			mxArray *mxEWRS;
			mxEWRS = mxGetField(mxInMSTEM, 0, "EWRS");
			InMSTEM.EWRS_ZeroDefTyp = ReadValuemxField<int>(mxEWRS, 0, "ZeroDefTyp");				// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
			InMSTEM.EWRS_ZeroDefPlane = ReadValuemxField<double>(mxEWRS, 0, "ZeroDefPlane");	// Zero defocus plane
			break;
		case 11:	// EW Fourier
			mxArray *mxEWFS;
			mxEWFS = mxGetField(mxInMSTEM, 0, "EWFS");
			InMSTEM.EWFS_ZeroDefTyp = ReadValuemxField<int>(mxEWFS, 0, "ZeroDefTyp");				// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
			InMSTEM.EWFS_ZeroDefPlane = ReadValuemxField<double>(mxEWFS, 0, "ZeroDefPlane");	// Zero defocus plane
			break;
	}
 }

/************************Check input TEM*************************/
void CheckMulSliPar(sInMSTEM &InMSTEMi, sInMSTEM &InMSTEMo){
	memcpy(&InMSTEMo, &InMSTEMi, sizeof(sInMSTEM));

	//if(InMSTEMo.gpu<0)			// gpu device 
	//	InMSTEMo.gpu = 0;

	//if(InMSTEMo.SimType<0)		// 1: STEM, 2: HRTEM, 3: ED, 4: PED, 5: CBED, 6: HCI, ... 10: EW real, 11: EW Fourier
	//	InMSTEMo.SimType = 10;

	//if((InMSTEMo.MulOrder<0)||(InMSTEMo.MulOrder>2)) // 1: First order MS, 2: Second Order MS
	//	InMSTEMo.MulOrder = 1;

	//if(InMSTEMo.nConfFP<0) // Number of frozen phonon configurations
	//	InMSTEMo.nConfFP = 0;

	//if((InMSTEMo.DimFP<1)||(InMSTEMo.DimFP>3)) // Dimensions phonon configurations
	//	InMSTEMo.DimFP = 3;

	//if(InMSTEMo.SeedFP<=0) // Random seed(frozen phonon)
	//	InMSTEMo.SeedFP = 1983;

	//if((InMSTEMo.MEffect<0)||(InMSTEMo.MEffect>2)) // 1: Exit wave Partial coherente mode, 2: Transmission cross coefficient
	//	InMSTEMo.MEffect = 0;

	//if((InMSTEMo.STEffect<0)||(InMSTEMo.STEffect>2)) // 1: Spatial and temporal, 2: Temporal, 3: Spatial
	//	InMSTEMo.STEffect = 0;

	//if((InMSTEMo.ZeroDefTyp<0)||(InMSTEMo.ZeroDefTyp>3)) // 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
	//	InMSTEMo.ZeroDefTyp = 1;

	//if((InMSTEMo.PotPar<1)||(InMSTEMo.PotPar>6)) // Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
	//	InMSTEMo.PotPar = 6;

	///***************************Multislice*************************/
	//if(InMSTEMo.E0<=0) // Acceleration voltage
	//	InMSTEMo.E0 = 300;

	//if((InMSTEMo.theta<0)||(InMSTEMo.theta>cPi)) // incident tilt (in spherical coordinates) (degrees-->rad)
	//	InMSTEMo.theta = 0;

	//if((InMSTEMo.phi<-c2Pi)||(InMSTEMo.phi>c2Pi)) // incident tilt (in spherical coordinates) (degrees-->rad)
	//	InMSTEMo.phi = 0;

	//if(InMSTEMo.nx<=0) // Number of pixels in x direction
	//	InMSTEMo.nx = 128;

	//if(InMSTEMo.ny<=0) // Number of pixels in y direction
	//	InMSTEMo.ny = 128;

	//if(InMSTEMo.lx<0) // distance in x direction(Angstroms)
	//	InMSTEMo.lx = 0;

	//if(InMSTEMo.ly<0) // distance in y direction(Angstroms)
	//	InMSTEMo.ly = 0;

	//if(InMSTEMo.dz<0) // slice thickness
	//	InMSTEMo.dz = 0;

	///*********************Microscope parameters********************/
	//mxArray *mxMC= mxGetField(mxInMSTEM, 0, "MC");
	//InMSTEMo.MC_m =  ReadValuemxField<int>(mxMC, 0, "m");						// momentum of the vortex
	//InMSTEMo.MC_f = ReadValuemxField<double>(mxMC, 0, "f");					// defocus(Angstrom)
	//InMSTEMo.MC_Cs3 = ReadValuemxField<double>(mxMC, 0, "Cs3", mm2Ags);		// spherical aberration(mm-->Angstrom)
	//InMSTEMo.MC_Cs5 = ReadValuemxField<double>(mxMC, 0, "Cs5", mm2Ags);		// spherical aberration(mm-->Angstrom)
	//InMSTEMo.MC_mfa2 = ReadValuemxField<double>(mxMC, 0, "mfa2");				// magnitude 2-fold astigmatism(Angstrom)
	//InMSTEMo.MC_afa2 = ReadValuemxField<double>(mxMC, 0, "afa2", deg2rad);		// angle 2-fold astigmatism(degrees-->rad)
	//InMSTEMo.MC_mfa3 = ReadValuemxField<double>(mxMC, 0, "mfa3");				// magnitude 3-fold astigmatism(Angstrom)
	//InMSTEMo.MC_afa3 = ReadValuemxField<double>(mxMC, 0, "afa3", deg2rad);		// angle 3-fold astigmatism(degrees-->rad)
	//InMSTEMo.MC_aobjl = ReadValuemxField<double>(mxMC, 0, "aobjl", mrad2rad);	// lower objective aperture(mrad-->rad)
	//InMSTEMo.MC_aobju = ReadValuemxField<double>(mxMC, 0, "aobju", mrad2rad);	// upper objective aperture(mrad-->rad)
	//InMSTEMo.MC_sf = ReadValuemxField<double>(mxMC, 0, "sf");					// defocus spread(Angstrom)
	//InMSTEMo.MC_nsf = ReadValuemxField<int>(mxMC, 0, "nsf");						// Number of defocus sampling point
	//InMSTEMo.MC_beta = ReadValuemxField<double>(mxMC, 0, "beta", mrad2rad);	// semi-convergence angle(mrad-->rad)
	//InMSTEMo.MC_nbeta = ReadValuemxField<int>(mxMC, 0, "nbeta");					// half number sampling points

	//mxArray *mxAtomsM = mxGetField(mxInMSTEM, 0, "Atoms");
	//InMSTEMo.nAtomsM = (int)mxGetM(mxAtomsM);	// Number of Atoms
	//InMSTEMo.AtomsM = mxGetPr(mxAtomsM);		// Atoms in a matrix form

	//switch (InMSTEMo.SimType){	
	//	case 0:	// STEM
	//		mxArray *mxSTEM;	
	//		mxSTEM = mxGetField(mxInMSTEM, 0, "STEM");
	//		InMSTEMo.STEM_line = ReadValuemxField<int>(mxSTEM, 0, "line");
	//		InMSTEMo.STEM_GridType = ReadValuemxField<int>(mxSTEM, 0, "GridType");
	//		InMSTEMo.STEM_ns = ReadValuemxField<int>(mxSTEM, 0, "ns");
	//		InMSTEMo.STEM_nucx = ReadValuemxField<int>(mxSTEM, 0, "nucx");
	//		InMSTEMo.STEM_nucy = ReadValuemxField<int>(mxSTEM, 0, "nucy");
	//		InMSTEMo.STEM_x1u = ReadValuemxField<double>(mxSTEM, 0, "x1u");
	//		InMSTEMo.STEM_y1u = ReadValuemxField<double>(mxSTEM, 0, "y1u");
	//		InMSTEMo.STEM_x2u = ReadValuemxField<double>(mxSTEM, 0, "x2u");
	//		InMSTEMo.STEM_y2u = ReadValuemxField<double>(mxSTEM, 0, "y2u");

	//		InMSTEMo.STEM_nDetCir = ReadValuemxField<int>(mxSTEM, 0, "nDetCir");
	//		InMSTEMo.STEM_DetCir = new sInDetCir[InMSTEMo.STEM_nDetCir];
	//		mxArray *mxDetCir;
	//		mxDetCir = mxGetField(mxSTEM, 0, "DetCir");
	//		for (int i=0; i<InMSTEMo.STEM_nDetCir; i++){
	//			InMSTEMo.STEM_DetCir[i].InnerAng = ReadValuemxField<double>(mxDetCir, i, "InnerAng", mrad2rad);	// Inner angle(mrad-->rad)
	//			InMSTEMo.STEM_DetCir[i].OuterAng = ReadValuemxField<double>(mxDetCir, i, "OuterAng", mrad2rad);	// Outer angle(mrad-->rad)
	//		}
	//		break;
	//	case 1:		// HRTEM
	//		mxArray *mxHRTEM;
	//		mxHRTEM = mxGetField(mxInMSTEM, 0, "HRTEM");

	//		break;
	//	case 2:		// ED
	//		mxArray *mxED;
	//		mxED = mxGetField(mxInMSTEM, 0, "ED");

	//		break;
	//	case 3:		// PED
	//		mxArray *mxPED;
	//		mxPED = mxGetField(mxInMSTEM, 0, "PED");
	//		InMSTEMo.PED_nrot = ReadValuemxField<int>(mxPED, 0, "nrot");				// Number of orientations
	//		InMSTEMo.PED_theta = ReadValuemxField<double>(mxPED, 0, "theta", deg2rad);	// Precession angle(degrees-->rad)
	//		break;
	//	case 4:		// CBED
	//		mxArray *mxCBED;
	//		mxCBED = mxGetField(mxInMSTEM, 0, "CBED");
	//		InMSTEMo.CBED_x0 = ReadValuemxField<double>(mxCBED, 0, "x0");	// Precession angle(degrees-->rad)
	//		InMSTEMo.CBED_y0 = ReadValuemxField<double>(mxCBED, 0, "y0");	// Precession angle(degrees-->rad)
	//		break;
	//	case 5:		// HCI
	//		mxArray *mxHCI;
	//		mxHCI = mxGetField(mxInMSTEM, 0, "HCI");
	//		InMSTEMo.HCI_nrot = ReadValuemxField<int>(mxHCI, 0, "nrot");				// Number of orientations
	//		InMSTEMo.HCI_theta = ReadValuemxField<double>(mxHCI, 0, "theta", deg2rad);	// Precession angle(degrees-->rad)
	//		break;
	//	case 10:	// EW real
	//		mxArray *mxEWr;
	//		mxEWr = mxGetField(mxInMSTEM, 0, "EWr");

	//		break;
	//	case 11:	// EW Fourier
	//		mxArray *mxEWF;
	//		mxEWF = mxGetField(mxInMSTEM, 0, "EWF");

	//		break;
	//}
 }