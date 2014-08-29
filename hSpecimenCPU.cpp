#include "hmathCPU.h"
#include <memory.h>
#include "hConstTypes.h"
#include "hSpecimenCPU.h"
#include "hAtomicData.h"
#include "hRandGen.h"

// Constructor
cSpecimen::cSpecimen(){	
	RandGen.reset();

	DimFP = 0;
	DistFP = 0;
	SeedFP = 0;

	nz = 0;
	Slice = 0;

	nAtoms = 0;
	Atomsu = 0;
	Atoms = 0;

	nAtomTypes = 0;
	AtomTypes = 0;
}

// Destructor
cSpecimen::~cSpecimen(){
	DimFP = 0;
	DistFP = 0;
	SeedFP = 0;

	ZeroDefTyp = 0;
	ZeroDefPlane = 0;
	dzBackProp = 0;

	dz = 0;

	nz = 0;
	delete [] Slice; Slice = 0;

	nAtoms = 0;
	delete [] Atomsu; Atomsu = 0;
	delete [] Atoms; Atoms = 0;

	AtomTypes = 0;
	delete [] AtomTypes; AtomTypes = 0;
}

// Set atoms: Copy input Atoms to the new format Atoms, count number of atoms, Ascending sort by z, Set relative atomic number position, Get maximum interaction distance
void cSpecimen::SetInputData(int nAtomsMi, double *AtomsMi, double dzi, int PotPari, double Vrli, int DimFPi=0, int DistFPi=0, int SeedFPi=0, int ZeroDefTypi=0, double ZeroDefPlanei=0){
	int i, j;
	double Z[NE];

	DimFP = DimFPi;					// Dimension frozen phonon
	DistFP = DistFPi;				// Random distribution generator
	SeedFP = SeedFPi;				// Random SeedFP generator

	ZeroDefTyp = ZeroDefTypi;		// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
	ZeroDefPlane = ZeroDefPlanei;	// Zero defocus plane
	dz = dzi;						// Set slice thickness

	nAtoms = nAtomsMi;
	for (i=0; i<NE; i++)
		Z[i] = 0.0;

	delete [] Atomsu; Atomsu = new sAtoms[nAtoms];
	// Set input Atomsu to the new format sAtoms
	for (i=0; i<nAtoms; i++){		
		Atomsu[i].x = AtomsMi[0*nAtoms + i];			// x-position
		Atomsu[i].y = AtomsMi[1*nAtoms + i];			// y-position
		Atomsu[i].z = AtomsMi[2*nAtoms + i];			// z-position
		Atomsu[i].Z = (int)AtomsMi[3*nAtoms + i];		// Atomic number
		Atomsu[i].sigma = AtomsMi[4*nAtoms + i];		// Standard deviation
		Atomsu[i].occ = AtomsMi[5*nAtoms + i];			// Occupancy
		
		// Get number of atom types
		if (Z[Atomsu[i].Z-1]==0.0)
			nAtomTypes++;

		Z[Atomsu[i].Z-1] += Atomsu[i].occ;
	}

	/**************************************************************************/
	// Get atom types
	j = 0;
	delete [] AtomTypes; AtomTypes = new sAtomTypesCPU[nAtomTypes];
	for (i=0; i<NE; i++){
		if (Z[i]>0){
			AtomTypes[j].Z = i+1;
			AtomTypes[j].occ = Z[i];
			AtomTypes[j].PotPar = PotPari;
			AtomTypes[j].ns = 0;
			//AtomTypesCPU[j].sigma = sigma[i]/double(Z[i]); This have to be modified
			j++;
		}
	}
	cAtomicData AtomicData;
	AtomicData.ReadAtomicData(nAtomTypes, AtomTypes, Vrli);

	/**************************************************************************/
	// z-position
	// Set relative atomic number positionand get the minimun and maximum z-atomic position
	for (i=0; i<nAtoms; i++){
		for (j=0; j<nAtomTypes; j++)
			if (Atomsu[i].Z==AtomTypes[j].Z){
				Atomsu[i].iZ = j;
				break;
			}
	}
	QuickSortAtomsAlongz(Atomsu, 0, nAtoms-1);			// Ascending sort by z

	// get Zero defocus plane
	switch(ZeroDefTyp){
		case 1:
			ZeroDefPlane = Atomsu[0].z;
			break;
		case 2:
			ZeroDefPlane = 0.5*(Atomsu[0].z+Atomsu[nAtoms-1].z);
			break;
		case 3:
			ZeroDefPlane = Atomsu[nAtoms-1].z;
			break;
	}

	Slicing(dz, nAtoms, Atomsu, AtomTypes, ZeroDefTyp, ZeroDefPlane, nz, Slice, dzBackProp);	// Slicing procedure

	// Copy Atomsu to Atomsd
	delete [] Atoms; Atoms = new sAtoms[nAtoms];
	memcpy(Atoms, Atomsu, nAtoms*cSizeofAtoms);
}

/*****************************************************************************/
// Sort atoms along z-axis
void cSpecimen::QuickSortAtomsAlongz(sAtoms *Atoms, int left, int right){
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

// random number generator Seed
void cSpecimen::setRandomSeed(unsigned long s, int iConf){	
	int randiu;
	RandGen.seed(s);
	for(int i=0; i<iConf; i++)
		randiu = RandGen.randiu();
	RandGen.seed(randiu);
}

// get dimension components
void cSpecimen::getDimCom(int Dim, double &bx, double &by, double &bz){
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
		case 011:
			bx = 0.0; by = 1.0; bz = 1.0;
			break;
		case 100:
			bx = 1.0; by = 0.0; bz = 0.0;
			break;
		case 010:
			bx = 0.0; by = 1.0; bz = 0.0;
			break;
		case 001:
			bx = 0.0; by = 0.0; bz = 1.0;
			break;
	}
}

/*****************************************************************************/
//// get atoms in the slice (1: There is an atom in the slice, >1: non atom in the slice 2: Bot, 3: middle, 4: top)
int cSpecimen::getAtomsInSlice(double z0, double ze, int nAtoms, sAtoms *Atoms, int &z0_id, int &ze_id){
	int i0, ie, im;
	double zmin = Atoms[0].z, zmax = Atoms[nAtoms-1].z;

	// Bottom
	if(z0 <= zmin) 
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
	if(ze >= zmax) 
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

	if((z0<=Atoms[z0_id].z)&&(Atoms[ze_id].z<=ze)){
		return 1;			// There is an atom
	}else if((z0<zmin)&&(ze<zmin)){
		z0_id = 1;	ze_id = 0;
		return 2;			// Bottom
	}else if((zmax<z0)&&(zmax<ze)){
		z0_id = 1;	ze_id = 0;
		return 4;			// Top
	}else{
		z0_id = 1;	ze_id = 0;
		return 3;			// Middle
	}
}

// get maximun interaction distance
double cSpecimen::getRint(int i0, int ie, sAtoms *Atoms, sAtomTypesCPU *AtomTypes, int Pos, double RBot, double Rmid, double RTop){
	double R, Rint = 0.0;
	switch(Pos){
		case 1:
			for(int i=i0; i<=ie; i++){
				R = AtomTypes[Atoms[i].iZ].Rmax;
				if (Rint < R)
					Rint = R;
			}
			break;
		case 2:
			Rint = RBot;
			break;
		case 3:
			Rint = Rmid;
			break;
		case 4:
			Rint = RTop;
			break;
	}

	return Rint;
}

// get slicing parameters
void cSpecimen::getSlicingParameters(double dz, int nAtoms, sAtoms *Atoms, sAtomTypesCPU *AtomTypes, int &nz, double &RBot, double &RTop){
	double z0, ze;
	int z0_id, ze_id;
	double zmin = Atoms[0].z, zmax = Atoms[nAtoms-1].z;

	// Get bottom interaction distance
	RBot = AtomTypes[Atoms[0].iZ].Rmax;
	z0 = zmin; ze = zmin + RBot;
	getAtomsInSlice(z0, ze, nAtoms, Atoms, z0_id, ze_id);
	RBot = getRint(z0_id, ze_id, Atoms, AtomTypes);

	// Get top interaction distance
	RTop = AtomTypes[Atoms[nAtoms-1].iZ].Rmax; 
	z0 = zmax - RTop; ze = zmax; 
	getAtomsInSlice(z0, ze, nAtoms, Atoms, z0_id, ze_id);
	RTop = getRint(z0_id, ze_id, Atoms, AtomTypes);
	
	double Lt = RBot + (zmax-zmin) + RTop;
	nz = ceil(Lt/dz);
	double dR = 0.5*(nz*dz-Lt);
	RBot += dR; RTop += dR;
}

// Select atoms respect to z-coordinate
void cSpecimen::Slicing(double dz, int nAtoms, sAtoms *Atoms, sAtomTypesCPU *AtomTypes, int ZeroDefTyp, double ZeroDefPlane, int &nz, sSlice *&Slice, double &dzBackProp){
	int pos;
	double RTop, RBot, zinit;

	// Get slicing parameters
	getSlicingParameters(dz, nAtoms, Atoms, AtomTypes, nz, RTop, RBot);

	delete [] Slice; Slice = new sSlice[nz];
	zinit = Atoms[0].z - RBot;
	for(int i=0; i<nz; i++){
		// Set integration plane
		Slice[i].zm = zinit + (double(i)+0.5)*dz;
		// Get atom's index in the slice
		Slice[i].z0 = Slice[i].zm - 0.5*dz; Slice[i].ze = Slice[i].zm + 0.5*dz;
		// Get atom's index in the slice
		pos = getAtomsInSlice(Slice[i].z0, Slice[i].ze, nAtoms, Atoms, Slice[i].z0_id, Slice[i].ze_id);
		// Get maximum interaction distance from the atoms inside the slice
		Slice[i].Rint = getRint(Slice[i].z0_id, Slice[i].ze_id, Atoms, AtomTypes, pos, RBot, (i==0)?RBot:Slice[i-1].Rint, RTop);
		// Get atom's index in the interaction slice
		Slice[i].z0i = Slice[i].zm - Slice[i].Rint; Slice[i].zei = Slice[i].zm + Slice[i].Rint;
		getAtomsInSlice(Slice[i].z0i, Slice[i].zei, nAtoms, Atoms, Slice[i].z0i_id, Slice[i].zei_id);
	}
	dzBackProp = ZeroDefPlane - Slice[nz-1].ze;
}

/*****************************************************************************/
// Move atoms (ramdom distribution will be included in the future)
void cSpecimen::MoveAtoms(int iConf){
	double sigmax, sigmay, sigmaz;
	double bx, by, bz;

	if (iConf<=0)
		return;

	getDimCom(DimFP, bx, by, bz);		// Get dimension components
	setRandomSeed(SeedFP, iConf);
	for (int i=0; i<nAtoms; i++){
		Atoms[i].Z = Atomsu[i].Z;
		Atoms[i].iZ = Atomsu[i].iZ;
		sigmax = sigmay = sigmaz = Atomsu[i].sigma;
		Atoms[i].x = Atomsu[i].x + bx*sigmax*RandGen.randn();
		Atoms[i].y = Atomsu[i].y + by*sigmay*RandGen.randn();
		Atoms[i].z = Atomsu[i].z + bz*sigmaz*RandGen.randn();
		Atoms[i].sigma = 0.0;
		Atoms[i].occ = Atomsu[i].occ;
	}
	
	QuickSortAtomsAlongz(Atoms, 0, nAtoms-1);				// Ascending sort by z	
	Slicing(dz, nAtoms, Atoms, AtomTypes, ZeroDefTyp, ZeroDefPlane, nz, Slice, dzBackProp);		// Slicing procedure
}