#include "hmathCPU.h"
#include <memory.h>
#include "hConstTypes.h"
#include "hgeneralCPU.h"
#include "hSpecimenCPU.h"
#include "hAtomicData.h"
#include "hRandGen.h"

void cSpecimenCPU::freeMemory(){
	f_sMPG_Init(MGP);

	Lzu = 0.0;
	Lztu = 0.0;
	Rmax = 0.0;
	z_BackProp = 0.0;

	nPlanes = 0;
	delete [] Planes; Planes = 0;

	nSlice = 0;
	delete [] Slice; Slice = 0;

	nAtomsu = 0;
	delete [] Atomsu; Atomsu = 0;

	nAtoms = 0;
	delete [] Atoms; Atoms = 0;

	nAtomTypes = 0;
	delete [] AtomTypes; AtomTypes = 0;
}

// Constructor
cSpecimenCPU::cSpecimenCPU(){
	f_sMPG_Init(MGP);

	Lzu = 0.0;
	Lztu = 0.0;
	Rmax = 0.0;
	z_BackProp = 0.0;

	nPlanes = 0;
	Planes = 0;

	nSlice = 0;
	Slice = 0;

	nAtomsu = 0;
	Atomsu = 0;

	nAtoms = 0;
	Atoms = 0;

	nAtomTypes = 0;
	AtomTypes = 0;
}

// Destructor
cSpecimenCPU::~cSpecimenCPU(){
	freeMemory();
}

// random number generator Seed
void cSpecimenCPU::setRandomSeed(unsigned long s, int iConf){	
	int randiu;
	RandGen.seed(s);
	for(int i=0; i<iConf; i++)
		randiu = RandGen.randiu();
	RandGen.seed(randiu);
}

// Sort atoms along z-axis
void cSpecimenCPU::QuickSortAtomsAlongz(sAtoms *&Atoms, int left, int right){
	int i = left, j = right;
	double pivot = Atoms[(left + right)>>1].z;
	sAtoms Atomst;

	/* partition */
	while (i <= j){
		while (Atoms[i].z < pivot)
		i++;
		while (Atoms[j].z > pivot)
		j--;

		if (i <= j){
			memcpy(&Atomst, &Atoms[i], cSizeofAtoms);
			memcpy(&Atoms[i], &Atoms[j], cSizeofAtoms);
			memcpy(&Atoms[j], &Atomst, cSizeofAtoms);
			i++;
			j--;
		}
	};

	/* recursion */
	if (left < j)
		QuickSortAtomsAlongz(Atoms, left, j);
	if (i < right)
		QuickSortAtomsAlongz(Atoms, i, right);
}

// get dimension components
void cSpecimenCPU::getDimCom(int Dim, double &bx, double &by, double &bz){
	bx = by = bz = 0.0;
	switch (Dim){
		case 111:
			bx = 1.0; by = 1.0; bz = 1.0;
			break;
		case 110:
			bx = 1.0; by = 1.0; bz = 0.0;
			break;
		case 101:
			bx = 1.0; by = 0.0; bz = 1.0;
			break;
		case 11:
			bx = 0.0; by = 1.0; bz = 1.0;
			break;
		case 100:
			bx = 1.0; by = 0.0; bz = 0.0;
			break;
		case 10:
			bx = 0.0; by = 1.0; bz = 0.0;
			break;
		case 1:
			bx = 0.0; by = 0.0; bz = 1.0;
			break;
	}
}

// get planes
void cSpecimenCPU::getPlanes(int nAtoms, sAtoms *&Atoms, int &nPlanes, double *&Planes){
	delete [] Planes; Planes = 0;

	double zmin = Atoms[0].z, zmax = Atoms[nAtoms-1].z;
	double Lz = zmax-zmin, dzp = 1e-03;
	int izp, iAtoms;
	sHt Ht;

	Ht.n = (int)ceil((Lz+dzp)/dzp);
	Ht.c = new int [Ht.n];
	Ht.v = new double [Ht.n];

	for(izp=0; izp<Ht.n; izp++)
		Ht.v[izp] = Ht.c[izp] = 0;

	for(iAtoms=0; iAtoms<nAtoms; iAtoms++){
		izp = (int)floor((Atoms[iAtoms].z-zmin)/dzp+0.5);
		if(izp<Ht.n){
			Ht.c[izp]++;
			Ht.v[izp] += Atoms[iAtoms].z;
		}
	}

	nPlanes = 0;
	for(izp=0; izp<Ht.n-1; izp++)
		if(Ht.c[izp]>0){		
			if(Ht.c[izp+1]>0){
				Ht.c[izp] += Ht.c[izp+1];
				Ht.v[izp] += Ht.v[izp+1];
				Ht.v[izp+1] = Ht.c[izp+1] = 0;
			}
			Ht.v[nPlanes++] = Ht.v[izp]/double(Ht.c[izp]);	
		}
	if(Ht.c[Ht.n-1]>0)
		Ht.v[nPlanes++] = Ht.v[Ht.n-1]/double(Ht.c[Ht.n-1]);

	Planes = new double[nPlanes];
	for(izp=0; izp<nPlanes; izp++)
		Planes[izp] = Ht.v[izp];

	Ht.n = 0;
	delete [] Ht.c; Ht.c = 0;
	delete [] Ht.v; Ht.v = 0;
}

//// get atoms in the slice (1:Bot, 2: middle, 3: top) - positive found an atom if not it is negative.
int cSpecimenCPU::getAtomsInSlice(double z0, double ze, int nAtoms, sAtoms *&Atoms, int &z0_id, int &ze_id){
	int i0, ie, im;
	double za_min = Atoms[0].z, za_max = Atoms[nAtoms-1].z;

	z0_id = 1;	ze_id = 0;
	if((z0<za_min)&&(ze<za_min))
		return 1;			// Bottom
	else if((za_max<z0)&&(za_max<ze))
		return 3;			// Top

	// Bottom
	if(z0 <= za_min) 
		z0_id = 0;
	else { 
		i0 = 0; ie = nAtoms-1;
		do{
			im = (i0 + ie)>>1;		// divide by 2
			if(z0 <= Atoms[im].z) 
				ie = im;
			else 
				i0 = im;
		}while ((ie-i0)>1);
		z0_id = (z0==Atoms[ie].z)?ie:(i0+1);
	}

	// Top
	if(ze >= za_max) 
		ze_id = nAtoms-1;
	else { 
		i0 = z0_id; ie = nAtoms-1;
		do{
			im = (i0 + ie)>>1;		// divide by 2
			if(ze < Atoms[im].z) 
				ie = im;
			else 
				i0 = im;
		}while ((ie-i0)>1);
		ze_id = (ze==Atoms[ie].z)?ie:i0;
	}

 if((Atoms[z0_id].z<z0)||(ze<Atoms[ze_id].z)){
 z0_id = 1;	ze_id = 0;
 }

	return 2;
}

// Select atoms respect to z-coordinate
void cSpecimenCPU::Slicing(int ApproxModel, double dz, double Rmax, int nAtoms, sAtoms *&Atoms, int ZeroDefTyp, double ZeroDefPlane, int &nSlice, sSlice *&Slice, double &z_BackProp){
	delete [] Slice; Slice = 0;
	if(ApproxModel>1){
		nSlice = 1;
		Slice = new sSlice[nSlice];
		// Set integration plane
		Slice[0].zm = 0.5*(Atoms[0].z+Atoms[nAtoms-1].z);
		// Get atom's index in the slice
		Slice[0].z0 = Atoms[0].z; Slice[0].ze = Atoms[nAtoms-1].z;
		// Set interaction distance
		Slice[0].Rint = Rmax;
		// Get atom's index in the interaction slice
		Slice[0].z0i = Atoms[0].z - Slice[0].Rint; Slice[0].zei = Atoms[nAtoms-1].z + Slice[0].Rint;
		// Get atom's index in the slice and the interaction slice
		Slice[0].z0_id = 0; Slice[0].ze_id = nAtoms-1;
		Slice[0].z0i_id = 0; Slice[0].zei_id = nAtoms-1;
		z_BackProp = 0.0;
		return;
	}

	double zmin = Atoms[0].z, Lz = Atoms[nAtoms-1].z-zmin;
	nSlice = ceil(Lz/dz);
	if(nSlice*dz<=Lz) nSlice++;
	double RBor = 0.5*(nSlice*dz-Lz);

	if(Rmax > RBor){
		int nRBor = (int)floor((Rmax-RBor)/dz);
		RBor += dz*nRBor;
		nSlice = nSlice + 2*nRBor;
	}
	
	double dzb = dz+Rmax-RBor, dzt;
	double z0 = zmin-Rmax;
	Slice = new sSlice[nSlice];
	for(int iSlice=0; iSlice<nSlice; iSlice++){
		dzt = ((iSlice==0)||(iSlice==nSlice-1))?dzb:dz;
		// Get atom's index in the slice
		Slice[iSlice].z0 = z0; Slice[iSlice].ze = z0 += dzt;
		// Set integration plane
		Slice[iSlice].zm = 0.5*(Slice[iSlice].z0+Slice[iSlice].ze);
		// Get atom's index in the slice
		getAtomsInSlice(Slice[iSlice].z0, Slice[iSlice].ze, nAtoms, Atoms, Slice[iSlice].z0_id, Slice[iSlice].ze_id);
		// Set interaction distance
		Slice[iSlice].Rint = Rmax;
		// Get atom's index in the interaction slice
		Slice[iSlice].z0i = Slice[iSlice].zm - Slice[iSlice].Rint; Slice[iSlice].zei = Slice[iSlice].zm + Slice[iSlice].Rint;
		if(Slice[iSlice].z0i>Slice[iSlice].z0) Slice[iSlice].z0i = Slice[iSlice].z0;
		if(Slice[iSlice].zei<Slice[iSlice].ze) Slice[iSlice].zei = Slice[iSlice].ze;
		getAtomsInSlice(Slice[iSlice].z0i, Slice[iSlice].zei, nAtoms, Atoms, Slice[iSlice].z0i_id, Slice[iSlice].zei_id);
	}
	z_BackProp = ZeroDefPlane - Slice[nSlice-1].ze;
}

/*****************************************************************************/
// Set atoms: Copy input Atoms to the new format Atoms, count number of atoms, Ascending sort by z, Set relative atomic number position, Get maximum interaction distance
void cSpecimenCPU::SetInputData(sMGP &MGP_io, int nAtomsM_i, double *AtomsM_i){
	freeMemory();

	MGP = MGP_io;
	nAtomTypes = NE;
	AtomTypes = new sAtomTypesCPU[nAtomTypes];
	f_SetAtomTypes(MGP.PotPar, 0, MGP.Vrl, nAtomTypes, AtomTypes);
	/**************************************************************************/
	f_AtomsM2Atoms(nAtomsM_i, AtomsM_i, MGP.PBC_xy, MGP.lx, MGP.ly, nAtomsu, Atomsu);
	QuickSortAtomsAlongz(Atomsu, 0, nAtomsu-1);	// Ascending sort by z
	Rmax = f_getRMax(nAtomsu, Atomsu, AtomTypes); 
	/**************************************************************************/
	Lzu = Atomsu[nAtomsu-1].z-Atomsu[0].z;
	Lztu = Lzu + 2.0*Rmax;
	if(Lzu==0) Lzu = Lztu;

	if(Lzu<MGP.dz){
		MGP.MulOrder = 1;
		MGP.dz = Lzu;
		if(MGP.ApproxModel==1)
			MGP.ApproxModel = 3;
	}
	// get Zero defocus plane
	switch(MGP.ZeroDefTyp){
		case 1:
			MGP.ZeroDefPlane = Atomsu[0].z;
			break;
		case 2:
			MGP.ZeroDefPlane = 0.5*(Atomsu[0].z+Atomsu[nAtomsu-1].z);
			break;
		case 3:
			MGP.ZeroDefPlane = Atomsu[nAtomsu-1].z;
			break;
	}
	if(MGP.ApproxModel>1)
		MGP.ZeroDefPlane = 0.0;

	// Slicing procedure
	Slicing(MGP.ApproxModel, MGP.dz, Rmax, nAtomsu, Atomsu, MGP.ZeroDefTyp, MGP.ZeroDefPlane, nSlice, Slice, z_BackProp);
	
	// Copy Atomsu to Atoms
	nAtoms = nAtomsu;
	Atoms = new sAtoms[nAtoms];
	memcpy(Atoms, Atomsu, nAtoms*cSizeofAtoms);

	// get planes
	getPlanes(nAtoms, Atoms, nPlanes, Planes);
	/**************************************************************************/
	MGP_io = MGP;
}

// Move atoms (ramdom distribution will be included in the future)
void cSpecimenCPU::MoveAtoms(int iConf){
	double sigmax, sigmay, sigmaz;
	double bx, by, bz;

	if (iConf<=0)
		return;

	// Get dimension components
	getDimCom(MGP.DimFP, bx, by, bz);
	setRandomSeed(MGP.SeedFP, iConf);
	for (int i=0; i<nAtoms; i++){
		Atoms[i].Z = Atomsu[i].Z;
		sigmax = sigmay = sigmaz = Atomsu[i].sigma;
		Atoms[i].x = Atomsu[i].x + bx*sigmax*RandGen.randn();
		Atoms[i].y = Atomsu[i].y + by*sigmay*RandGen.randn();
		Atoms[i].z = Atomsu[i].z + bz*sigmaz*RandGen.randn();
		Atoms[i].sigma = 0.0;
		Atoms[i].occ = Atomsu[i].occ ;
	}
	// Ascending sort by z
	QuickSortAtomsAlongz(Atoms, 0, nAtoms-1);	
	// Slicing procedure
	Slicing(MGP.ApproxModel, MGP.dz, Rmax, nAtoms, Atoms, MGP.ZeroDefTyp, MGP.ZeroDefPlane, nSlice, Slice, z_BackProp);	
	// get planes
	getPlanes(nAtoms, Atoms, nPlanes, Planes);
}