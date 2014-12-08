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

#include <cstring>
#include "hmathCPU.h"
#include "hConstTypes.h"
#include "hMT_MGP_CPU.h"
#include "hRandGen.h"
#include "hMT_General_CPU.h"

#include "hMT_Specimen_CPU.h"

void cMT_Specimen_CPU::freeMemory()
{
	MT_MGP_CPU = 0;

	sigma_min = 0.0;
	sigma_max = 0.0;

	Lzu = 0.0;
	Lztu = 0.0;
	Rmax = 0.0;
	z_BackProp = 0.0;

	nAtomsu = 0;
	delete [] Atomsu; Atomsu = 0;
	
	nSliceu = 0;
	delete [] Sliceu; Sliceu = 0;

	nPlanesu = 0;
	delete [] Planesu; Planesu = 0;

	nAtoms = 0;
	delete [] Atoms; Atoms = 0;

	nSlice = 0;
	delete [] Slice; Slice = 0;

	nMT_AtomTypes = 0;
	delete [] MT_AtomTypes_CPU; MT_AtomTypes_CPU = 0;
}

// Constructor
cMT_Specimen_CPU::cMT_Specimen_CPU()
{
	MT_MGP_CPU = 0;

	sigma_min = 0.0;
	sigma_max = 0.0;

	Lzu = 0.0;
	Lztu = 0.0;
	Rmax = 0.0;
	z_BackProp = 0.0;

	nAtomsu = 0;
	Atomsu = 0;

	nSliceu = 0;
	Sliceu = 0;

	nPlanesu = 0;
	Planesu = 0;

	nAtoms = 0;
	Atoms = 0;

	nSlice = 0;
	Slice = 0;

	nMT_AtomTypes = 0;
	MT_AtomTypes_CPU = 0;
}

// Destructor
cMT_Specimen_CPU::~cMT_Specimen_CPU(){
	freeMemory();
}

// random number generator Seed
void cMT_Specimen_CPU::setRandomSeed(unsigned long s, int iConf){	
	int randiu;
	RandGen.seed(s);
	for(int iSlice = 0; iSlice<iConf; iSlice++)
		randiu = RandGen.randiu();
	RandGen.seed(randiu);
}

// Sort atoms along z-axis
void cMT_Specimen_CPU::QuickSortAtomsAlongz(sAtoms *&Atoms, int left, int right){
	int iSlice = left, j = right;
	double pivot = Atoms[(left + right)>>1].z;
	sAtoms Atomst;

	/* partition */
	while (iSlice <= j){
		while (Atoms[iSlice].z < pivot)
		iSlice++;
		while (Atoms[j].z > pivot)
		j--;

		if (iSlice <= j){
			memcpy(&Atomst, &Atoms[iSlice], cSizeofAtoms);
			memcpy(&Atoms[iSlice], &Atoms[j], cSizeofAtoms);
			memcpy(&Atoms[j], &Atomst, cSizeofAtoms);
			iSlice++;
			j--;
		}
	};

	/* recursion */
	if (left < j)
		QuickSortAtomsAlongz(Atoms, left, j);
	if (iSlice < right)
		QuickSortAtomsAlongz(Atoms, iSlice, right);
}

// get dimension components
void cMT_Specimen_CPU::getDimCom(int Dim, double &bx, double &by, double &bz){
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
void cMT_Specimen_CPU::getPlanes(int nAtoms, sAtoms *&Atoms, int &nPlanesu, double *&Planesu){
	delete [] Planesu; Planesu = 0;

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

	nPlanesu = 0;
	for(izp=0; izp<Ht.n-1; izp++)
		if(Ht.c[izp]>0){		
			if(Ht.c[izp+1]>0){
				Ht.c[izp] += Ht.c[izp+1];
				Ht.v[izp] += Ht.v[izp+1];
				Ht.v[izp+1] = Ht.c[izp+1] = 0;
			}
			Ht.v[nPlanesu++] = Ht.v[izp]/double(Ht.c[izp]);	
		}
	if(Ht.c[Ht.n-1]>0)
		Ht.v[nPlanesu++] = Ht.v[Ht.n-1]/double(Ht.c[Ht.n-1]);

	Planesu = new double[nPlanesu];
	for(izp=0; izp<nPlanesu; izp++)
		Planesu[izp] = Ht.v[izp];

	Ht.n = 0;
	delete [] Ht.c; Ht.c = 0;
	delete [] Ht.v; Ht.v = 0;
}

// get border
void cMT_Specimen_CPU::getnSlice(double z10, double z1i, double z20, double z2i, double dzi, int &nSlice, double &dz0, double &dze){
	double z0 = z1i, ze = z2i;
	double Lzt = z2i-z1i;

	z10=-z10; z1i=-z1i;

	dz0 = (z1i-z10)-(int)floor((z1i-z10)/dzi)*dzi;
	if(abs(dz0)<((dzi>2.0)?0.25:0.5)*dzi)
		dz0 += dzi;
	/*********************************************************/
	dze = (z2i-z20)-(int)floor((z2i-z20)/dzi)*dzi;
	if(abs(dze)<((dzi>2.0)?0.25:0.5)*dzi)
		dze += dzi;
	/*********************************************************/
	nSlice=1;
	double z=z0+dz0;
	while((z<ze+eed)&&(z-dzi+dze<ze+eed)){
		nSlice++;
		z+=dzi;
	}
	if(dz0+dze+(nSlice-2)*dzi>Lzt+eed)
		nSlice--;
}

// get atoms in the Slice (1:Bot, 2: middle, 3: top) - positive found an atom if not it is negative.
int cMT_Specimen_CPU::getAtomsInSlice(double z0, double ze, int nAtoms, sAtoms *&Atoms, int &z0_id, int &ze_id){
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

// get border
int cMT_Specimen_CPU::getBorderSlicing(double z0, double ze, double zi, double dzi, double &dzb){

	if(zi<z0){
		z0=-z0; ze=-ze; zi=-zi;
	};

	int nzb=0;
	double z = z0;
	while(z<zi){
		z+=dzi;
		nzb++;
	}
	nzb=(nzb>2)?nzb-2:0;
	dzb = zi-(z-dzi);
	if(dzb>=((dzi>2.0)?0.25:0.5)*dzi)
		nzb++;
	else
		dzb = dzb+dzi;

	return nzb;
}

// Select atoms respect to z-coordinate
void cMT_Specimen_CPU::Slicing(double Rmax, int nAtoms, sAtoms *&Atoms, int &nSlice, sSlice *&Slice, double &z_BackProp){
	delete [] Slice; Slice = 0;
	double zmin = Atoms[0].z, zmax = Atoms[nAtoms-1].z;
	double Lz = zmax - zmin;

	if(MT_MGP_CPU->ApproxModel>1){
		nSlice = 1;
		Slice = new sSlice[nSlice];
		// Get atom's index in the Slice
		Slice[0].z0 = zmin - Rmax; Slice[0].ze = zmax + Rmax;
		// Get atom's index in the interaction Slice
		Slice[0].z0i = zmin - Rmax; Slice[0].zei = zmax + Rmax;
		// Get atom's index in the Slice and the interaction Slice
		Slice[0].z0_id = 0; Slice[0].ze_id = nAtoms-1;
		Slice[0].z0i_id = 0; Slice[0].zei_id = nAtoms-1;
		z_BackProp = 0.0;
		return;
	}

	nSlice = (int)ceil(Lz/MT_MGP_CPU->dz);
	if(nSlice*MT_MGP_CPU->dz<=Lz) nSlice++;

	double dz0, dze;
	double z10 = zmin-0.5*(nSlice*MT_MGP_CPU->dz-Lz)+MT_MGP_CPU->dz, z20 = zmax+0.5*(nSlice*MT_MGP_CPU->dz-Lz)-MT_MGP_CPU->dz;
	getnSlice(z10, zmin-Rmax, z20, zmax+Rmax, MT_MGP_CPU->dz, nSlice, dz0, dze);
	if(nSlice==1) dz0 = dze = 2.0*Rmax;

	double dz, zm, z0 = zmin-Rmax;
	Slice = new sSlice[nSlice];
	for(int iSlice = 0; iSlice<nSlice; iSlice++){
		dz = (iSlice==0)?dz0:(iSlice==nSlice-1)?dze:MT_MGP_CPU->dz;
		// Get atom's index in the Slice
		Slice[iSlice].z0 = z0; Slice[iSlice].ze = z0 += dz;
		// Get atom's index in the Slice
		getAtomsInSlice(Slice[iSlice].z0, Slice[iSlice].ze, nAtoms, Atoms, Slice[iSlice].z0_id, Slice[iSlice].ze_id);
		// Set integration plane
		zm = 0.5*(Slice[iSlice].z0+Slice[iSlice].ze);
		// Get atom's index in the interaction Slice
		Slice[iSlice].z0i = zm - Rmax; Slice[iSlice].zei = zm + Rmax;
		if(Slice[iSlice].z0i>Slice[iSlice].z0) Slice[iSlice].z0i = Slice[iSlice].z0;
		if(Slice[iSlice].zei<Slice[iSlice].ze) Slice[iSlice].zei = Slice[iSlice].ze;
		getAtomsInSlice(Slice[iSlice].z0i, Slice[iSlice].zei, nAtoms, Atoms, Slice[iSlice].z0i_id, Slice[iSlice].zei_id);
	}
	z_BackProp = 0.0;
}

// Select atoms respect to z-coordinate
void cMT_Specimen_CPU::Slicing(double Rmax, int nSliceu, sSlice *Sliceu, int nAtoms, sAtoms *&Atoms, int &nSlice, sSlice *&Slice, double &z_BackProp){
	delete [] Slice; Slice = 0;
	double zmin = Atoms[0].z, zmax = Atoms[nAtoms-1].z;
	double Lzt = zmax-zmin+2.0*Rmax;

	if(nSliceu==1){
		nSlice = 1;
		Slice = new sSlice[nSlice];
		// Get atom's index in the Slice
		Slice[0].z0 = zmin - Rmax; Slice[0].ze = zmax + Rmax;
		// Get atom's index in the Slice and the interaction Slice
		Slice[0].z0_id = 0; Slice[0].ze_id = nAtoms-1;
		// Get atom's index in the interaction Slice
		Slice[0].z0i = zmin - Rmax; Slice[0].zei = zmax + Rmax;
		Slice[0].z0i_id = 0; Slice[0].zei_id = nAtoms-1;
		z_BackProp = 0.0;
		return;
	}

	double dz0, dze;
	double z10 = Sliceu[0].ze, z20 = Sliceu[nSliceu-1].z0;
	getnSlice(z10, zmin-Rmax, z20, zmax+Rmax, MT_MGP_CPU->dz, nSlice, dz0, dze);

	double dz, zm, z0 = zmin-Rmax;
	Slice = new sSlice[nSlice];
	for(int iSlice = 0; iSlice<nSlice; iSlice++){
		dz = (iSlice==0)?dz0:(iSlice==nSlice-1)?dze:MT_MGP_CPU->dz;
		// Get atom's index in the Slice
		Slice[iSlice].z0 = z0; Slice[iSlice].ze = z0 += dz;
		// Get atom's index in the Slice
		getAtomsInSlice(Slice[iSlice].z0, Slice[iSlice].ze, nAtoms, Atoms, Slice[iSlice].z0_id, Slice[iSlice].ze_id);
		// Set integration plane
		zm = 0.5*(Slice[iSlice].z0+Slice[iSlice].ze);
		// Get atom's index in the interaction Slice
		Slice[iSlice].z0i = zm - Rmax; Slice[iSlice].zei = zm + Rmax;
		if(Slice[iSlice].z0i>Slice[iSlice].z0) Slice[iSlice].z0i = Slice[iSlice].z0;
		if(Slice[iSlice].zei<Slice[iSlice].ze) Slice[iSlice].zei = Slice[iSlice].ze;
		getAtomsInSlice(Slice[iSlice].z0i, Slice[iSlice].zei, nAtoms, Atoms, Slice[iSlice].z0i_id, Slice[iSlice].zei_id);
	}
	z_BackProp = MT_MGP_CPU->ZeroDefPlane - Slice[nSlice-1].ze;
}

// Select atoms respect to z-coordinate
void cMT_Specimen_CPU::Slicing(int nPlanesu, double *Planesu, int nAtoms, sAtoms *&Atoms, int &nSlice, sSlice *&Slice, double &z_BackProp){
	nSlice = nPlanesu;
	delete [] Slice; Slice = 0;
	double dzh;
	Slice = new sSlice[nPlanesu];
	for(int iSlice = 0; iSlice<nSlice; iSlice++){
		// Get atom's index in the Slice
		if(iSlice<nSlice-1){
			dzh = 0.5*(Planesu[iSlice+1]-Planesu[iSlice]);
			Slice[iSlice].z0 = Planesu[iSlice]-dzh; Slice[iSlice].ze = Planesu[iSlice]+dzh;
		}
		// Get atom's index in the interaction Slice
		Slice[iSlice].z0i = Planesu[iSlice]-0.001; Slice[iSlice].zei = Planesu[iSlice]+0.001;
		getAtomsInSlice(Slice[iSlice].z0i, Slice[iSlice].zei, nAtoms, Atoms, Slice[iSlice].z0i_id, Slice[iSlice].zei_id);
		Slice[iSlice].z0_id = Slice[iSlice].z0i_id; Slice[iSlice].ze_id = Slice[iSlice].zei_id;
	}
	if(nSlice>1)
		Slice[nSlice-1].z0 = Planesu[nSlice-1]-dzh; Slice[nSlice-1].ze = Planesu[nSlice-1]+dzh;

	z_BackProp = MT_MGP_CPU->ZeroDefPlane - Planesu[nSlice-1];
}

// Set atoms: Copy input Atoms to the new format Atoms, count number of atoms, Ascending sort by z, Set relative atomic number position, Get maximum interaction distance
void cMT_Specimen_CPU::SetInputData(cMT_MGP_CPU *MT_MGP_CPU_io, int nAtomsM_i, double *AtomsM_i, double dRmin){
	freeMemory();

	MT_MGP_CPU = MT_MGP_CPU_io;
	nMT_AtomTypes = NE;
	MT_AtomTypes_CPU = new cMT_AtomTypes_CPU[nMT_AtomTypes];
	for(int i=0; i<nMT_AtomTypes; i++)
		MT_AtomTypes_CPU[i].SetAtomTypes(i+1, MT_MGP_CPU->PotPar, MT_MGP_CPU->Vrl, stnR, dRmin);
	/*************************************************************************/
	f_AtomsM2Atoms(nAtomsM_i, AtomsM_i, MT_MGP_CPU->PBC_xy, MT_MGP_CPU->lx, MT_MGP_CPU->ly, nAtomsu, Atomsu, sigma_min, sigma_max);
	QuickSortAtomsAlongz(Atomsu, 0, nAtomsu-1);	// Ascending sort by z
	Rmax = f_getRMax(nAtomsu, Atomsu, MT_AtomTypes_CPU); 
	/*************************************************************************/
	Lzu = Atomsu[nAtomsu-1].z-Atomsu[0].z;
	Lztu = Lzu + 2.0*Rmax;
	if(Lzu==0) Lzu = Lztu;

	if(Lzu<MT_MGP_CPU->dz){
		MT_MGP_CPU->MulOrder = 1;
		MT_MGP_CPU->dz = Lzu;
		if(MT_MGP_CPU->ApproxModel==1)
			MT_MGP_CPU->ApproxModel = 3;
	}
	// get Zero defocus plane
	switch(MT_MGP_CPU->ZeroDefTyp){
		case 1:
			MT_MGP_CPU->ZeroDefPlane = Atomsu[0].z;
			break;
		case 2:
			MT_MGP_CPU->ZeroDefPlane = 0.5*(Atomsu[0].z+Atomsu[nAtomsu-1].z);
			break;
		case 3:
			MT_MGP_CPU->ZeroDefPlane = Atomsu[nAtomsu-1].z;
			break;
	}

	// get planes
	getPlanes(nAtomsu, Atomsu, nPlanesu, Planesu);

	// Slicing procedure
	if(MT_MGP_CPU->ApproxModel==1)
		Slicing(Rmax, nAtomsu, Atomsu, nSliceu, Sliceu, z_BackProp);
	else
		Slicing(nPlanesu, Planesu, nAtomsu, Atomsu, nSliceu, Sliceu, z_BackProp);
	
	/*************************************************************************/
	// Copy Atomsu to Atoms
	nAtoms = nAtomsu;
	Atoms = new sAtoms[nAtoms];
	memcpy(Atoms, Atomsu, nAtoms*cSizeofAtoms);

	// Copy Atomsu to Atoms
	nSlice = nSliceu;
	Slice = new sSlice[nSlice];
	memcpy(Slice, Sliceu, nSlice*cSizeofSlice);
}

// Move atoms (ramdom distribution will be included in the future)
void cMT_Specimen_CPU::MoveAtoms(int iConf){
	if (iConf<=0) return;

	double sigmax, sigmay, sigmaz;
	double bx, by, bz;

	// Get dimension components
	getDimCom(MT_MGP_CPU->DimFP, bx, by, bz);
	setRandomSeed(MT_MGP_CPU->SeedFP, iConf);
	for (int iSlice = 0; iSlice<nAtoms; iSlice++){
		Atoms[iSlice].Z = Atomsu[iSlice].Z;
		sigmax = sigmay = sigmaz = Atomsu[iSlice].sigma;
		Atoms[iSlice].x = Atomsu[iSlice].x + bx*sigmax*RandGen.randn();
		Atoms[iSlice].y = Atomsu[iSlice].y + by*sigmay*RandGen.randn();
		Atoms[iSlice].z = Atomsu[iSlice].z + bz*sigmaz*RandGen.randn();
		Atoms[iSlice].sigma = 0.0;
		Atoms[iSlice].occ = Atomsu[iSlice].occ ;
	}
	if(MT_MGP_CPU->ApproxModel==1){
		// Ascending sort by z
		QuickSortAtomsAlongz(Atoms, 0, nAtoms-1);	
		// Slicing procedure
		Slicing(Rmax, nSliceu, Sliceu, nAtoms, Atoms, nSlice, Slice, z_BackProp);
	}
}

// Get dz
double cMT_Specimen_CPU::get_dz(int iSlice){
	double dz = (iSlice<nSlice)?(Slice[iSlice].ze - Slice[iSlice].z0)/cos(MT_MGP_CPU->theta):0.0;
	return dz;
}