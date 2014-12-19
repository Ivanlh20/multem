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

#include "hMT_Slicing_CPU.h"

void cMT_Slicing_CPU::freeMemory()
{
	MT_MGP_CPU = 0;

	Rmax = 0;

	nAtomsu = 0;
	Atomsu = 0;
	
	nSliceu = 0;
	delete [] Sliceu; Sliceu = 0;

	nPlanesu = 0;
	delete [] Planesu; Planesu = 0;

	nAtoms = 0;
	Atoms = 0;

	nSlice = 0;
	delete [] Slice; Slice = 0;

	Thk.n = 0;
	delete [] Thk.h; Thk.h = 0;
	delete [] Thk.zb; Thk.zb = 0;
	delete [] Thk.iSlice; Thk.iSlice=0;
	delete [] Thk.iAtom; Thk.iAtom = 0;
}

// Constructor
cMT_Slicing_CPU::cMT_Slicing_CPU(){
	MT_MGP_CPU = 0;

	Rmax = 0;

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

	Thk.n = 0;
	Thk.h = 0;
	Thk.zb = 0;
	Thk.iSlice=0;
	Thk.iAtom = 0;
}

// Destructor
cMT_Slicing_CPU::~cMT_Slicing_CPU(){
	freeMemory();
}

// typ 1: bottom, 2: top
int cMT_Slicing_CPU::binary_search(double z0, int iAtom0, int iAtome, sAtoms *&Atoms, int typ){
	int i0, ie, im;
	double eps = 1e-12;
	double z0s = (typ==1)?z0-eps:z0+eps;

	if(z0 < Atoms[iAtom0].z) 
		i0 = iAtom0;
	else if(Atoms[iAtome].z<z0)
		i0 = iAtome;
	else{ 
		i0 = iAtom0; ie = iAtome;
		do{
			im = (i0 + ie)>>1;		// divide by 2
			if(z0s < Atoms[im].z)
				ie = im;
			else
				i0 = im;
		}while ((ie-i0)>1);

        if(i0>0){
			if(typ==1)
				i0++;
			else if(z0==Atoms[ie].z) 
				i0++;
		}
	}

	return i0;
}

// get thickness
void cMT_Slicing_CPU::getThk(int nAtoms, sAtoms *Atoms, int nSlice, sSlice *Slice, sThk &Thk){
	if((nSlice==1)||(MT_MGP_CPU->ThkTyp==1)){
		Thk.n = 1;
		delete [] Thk.h; Thk.h = new double[Thk.n];
		delete [] Thk.zb; Thk.zb = new double[Thk.n];
		delete [] Thk.iSlice; Thk.iSlice = new int[Thk.n];
		delete [] Thk.iAtom; Thk.iAtom = new int[Thk.n];

		Thk.h[0] = Atoms[nAtoms-1].z;
		Thk.zb[0] = MT_MGP_CPU->ZeroDefPlane - Slice[nSlice-1].ze;
		Thk.iSlice[0] = nSlice-1;
		Thk.iAtom[0] = nAtoms-1;
		return;
	}
	
	int nz_i, nz;
	double *z_i = 0, *z = 0;

	if(MT_MGP_CPU->ThkTyp==3){ 
		f_MatchTwoVectors(MT_MGP_CPU->nThk_i, MT_MGP_CPU->Thk_i, nPlanesu, Planesu, nz_i, z_i);
	}else{
		nz_i = MT_MGP_CPU->nThk_i;
		z_i = new double [nz_i];
		memcpy(z_i, MT_MGP_CPU->Thk_i, nz_i*cSizeofRD); 
	}

	nz = nSlice;
	z = new double [nz];
	for(int iz=0; iz<nz; iz++) 
		z[iz] = 0.5*(Slice[iz].z0+Slice[iz].ze);

	f_MatchTwoVectors(nz_i, z_i, nz, z, Thk.n, Thk.h, Thk.iSlice);
	nz = 0; delete [] z; z = 0;
	nz_i = 0; delete [] z_i; z_i = 0;

	delete [] Thk.zb; Thk.zb = new double[Thk.n];
	delete [] Thk.iAtom; Thk.iAtom = new int[Thk.n];
	for(int iThk=0; iThk<Thk.n; iThk++){
		Thk.iAtom[iThk] = binary_search(Thk.h[iThk], 0, nAtoms-1, Atoms, 2);
		Thk.zb[iThk] = MT_MGP_CPU->ZeroDefPlane - (Atoms[Thk.iAtom[iThk]].z+Rmax);
	}
}

// get planes
void cMT_Slicing_CPU::getPlanes(int nAtoms, sAtoms *&Atoms, int &nPlanesu, double *&Planesu){
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
void cMT_Slicing_CPU::get_nSlice(double z10, double z1i, double z20, double z2i, double dzi, int &nSlice, double &dz0, double &dze){
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
void cMT_Slicing_CPU::getAtomsInSlice(double z0, double ze, int nAtoms, sAtoms *&Atoms, int &z0_id, int &ze_id){
	z0_id = binary_search(z0, 0, nAtoms-1, Atoms, 1);
	ze_id = binary_search(ze, z0_id, nAtoms-1, Atoms, 2);

	 if((Atoms[z0_id].z<z0)||(ze<Atoms[ze_id].z)){
		z0_id = 1;	ze_id = 0;
	 }
}

// get dzh between planes
void cMT_Slicing_CPU::get_dzh_Planes(int nPlane, double *Plane, int iPlane, double &dzhd, double &dzhu){
	if(nPlane==1){
		dzhd = dzhu = 0.0;
		return;
	}

	if(iPlane==0){
		dzhd = dzhu = 0.5*(Plane[iPlane+1]-Plane[iPlane]);
	}else if(iPlane==nPlane-1){
		dzhd = dzhu = 0.5*(Plane[nPlane-1]-Plane[nPlane-2]);
	}else{
		dzhd = 0.5*(Plane[iPlane]-Plane[iPlane-1]);
		dzhu = 0.5*(Plane[iPlane+1]-Plane[iPlane]);
	}
}

// Set input data
void cMT_Slicing_CPU::SetInputData(cMT_MGP_CPU *MT_MGP_CPU_i, int nAtomsu_i, sAtoms *Atomsu_i, int nAtoms_i, sAtoms *Atoms_i, double Rmax_i){
	freeMemory();

	MT_MGP_CPU = MT_MGP_CPU_i;

	nAtomsu = nAtomsu_i;
	Atomsu = Atomsu_i;

	nAtoms = nAtoms_i; 
	Atoms = Atoms_i;

	Rmax = Rmax_i;

	// get planes
	getPlanes(nAtomsu, Atomsu, nPlanesu, Planesu);

	// Slicing procedure
	if(MT_MGP_CPU->ApproxModel!=2)
		Slicing_u(nAtomsu, Atomsu, nSliceu, Sliceu, Thk);
	else
		Slicing_by_Planes(nPlanesu, Planesu, nAtomsu, Atomsu, nSliceu, Sliceu, Thk);

	// Copy Atomsu to Atoms
	nSlice = nSliceu;
	Slice = new sSlice[nSlice];
	memcpy(Slice, Sliceu, nSlice*cSizeofSlice);
}

// Select atoms respect to z-coordinate
void cMT_Slicing_CPU::Slicing_u(int nAtoms, sAtoms *&Atoms, int &nSlice, sSlice *&Slice, sThk &Thk){
	delete [] Slice; Slice = 0;
	double zmin = Atoms[0].z, zmax = Atoms[nAtoms-1].z;
	double Lz = zmax - zmin;

	if(MT_MGP_CPU->ApproxModel>2){
		nSlice = 1;
		Slice = new sSlice[nSlice];
		// Get atom's index in the Slice
		Slice[0].z0 = zmin - Rmax; Slice[0].ze = zmax + Rmax;
		// Get atom's index in the Slice and the interaction Slice
		Slice[0].z0_id = 0; Slice[0].ze_id = nAtoms-1;
		// Get atom's index in the interaction Slice
		Slice[0].z0i = zmin - Rmax; Slice[0].zei = zmax + Rmax;
		// Get atom's index in the Slice and the interaction Slice
		Slice[0].z0i_id = 0; Slice[0].zei_id = nAtoms-1;
		// Get thickness
		getThk(nAtoms, Atoms, nSlice, Slice, Thk);
		return;
	}

	nSlice = (int)ceil(Lz/MT_MGP_CPU->dz);
	if(nSlice*MT_MGP_CPU->dz<=Lz) nSlice++;

	double dz0, dze;
	double z10 = zmin-0.5*(nSlice*MT_MGP_CPU->dz-Lz)+MT_MGP_CPU->dz, z20 = zmax+0.5*(nSlice*MT_MGP_CPU->dz-Lz)-MT_MGP_CPU->dz;
	get_nSlice(z10, zmin-Rmax, z20, zmax+Rmax, MT_MGP_CPU->dz, nSlice, dz0, dze);
	if(nSlice==1) dz0 = dze = 2.0*Rmax;

	double dz, zm, z0 = zmin-Rmax;
	Slice = new sSlice[nSlice];
	for(int iSlice=0; iSlice<nSlice; iSlice++){
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
	// Get thickness
	getThk(nAtoms, Atoms, nSlice, Slice, Thk);
}

// Select atoms respect to z-coordinate
void cMT_Slicing_CPU::Slicing_d(int nSliceu, sSlice *Sliceu, int nAtoms, sAtoms *&Atoms, int &nSlice, sSlice *&Slice, sThk &Thk){
	delete [] Slice; Slice = 0;
	double zmin = Atoms[0].z, zmax = Atoms[nAtoms-1].z;

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
		// Get thickness
		getThk(nAtoms, Atoms, nSlice, Slice, Thk);
		return;
	}

	double dz0, dze;
	double z10 = Sliceu[0].ze, z20 = Sliceu[nSliceu-1].z0;
	get_nSlice(z10, zmin-Rmax, z20, zmax+Rmax, MT_MGP_CPU->dz, nSlice, dz0, dze);

	double dz, zm, z0 = zmin-Rmax;
	Slice = new sSlice[nSlice];
	for(int iSlice=0; iSlice<nSlice; iSlice++){
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
	// Get thickness
	getThk(nAtoms, Atoms, nSlice, Slice, Thk);
}

// Select atoms respect to z-coordinate
void cMT_Slicing_CPU::Slicing_by_Planes(int nPlanesu, double *Planesu, int nAtoms, sAtoms *&Atoms, int &nSlice, sSlice *&Slice, sThk &Thk){	
	nSlice = nPlanesu;
	delete [] Slice; Slice = 0;
	Slice = new sSlice[nSlice];

	if(nSlice==1){
		Slice = new sSlice[nSlice];
		// Get atom's index in the Slice
		Slice[0].z0 = Planesu[0] - Rmax; Slice[0].ze = Planesu[0] + Rmax;
		// Get atom's index in the Slice and the interaction Slice
		Slice[0].z0_id = 0; Slice[0].ze_id = nAtoms-1;
		// Get atom's index in the interaction Slice
		Slice[0].z0i = Planesu[0] - Rmax; Slice[0].zei = Planesu[0] + Rmax;
		Slice[0].z0i_id = 0; Slice[0].zei_id = nAtoms-1;
		// Get thickness
		getThk(nAtoms, Atoms, nSlice, Slice, Thk);
		return;
	}

	double dzhu, dzhd, dz;
	for(int iSlice=0; iSlice<nSlice; iSlice++){
		dz = (iSlice<iSlice-1)?Planesu[iSlice+1]-Planesu[iSlice]:Planesu[nSlice-1]-Planesu[nSlice-2];
		Slice[iSlice].z0 = Planesu[iSlice]-0.5*dz; Slice[iSlice].ze = Planesu[iSlice]+0.5*dz;
		get_dzh_Planes(nPlanesu, Planesu, iSlice, dzhd, dzhu);
		Slice[iSlice].z0i = Planesu[iSlice]-dzhd; Slice[iSlice].zei = Planesu[iSlice]+dzhu;
		getAtomsInSlice(Slice[iSlice].z0i, Slice[iSlice].zei, nAtoms, Atoms, Slice[iSlice].z0i_id, Slice[iSlice].zei_id);

		Slice[iSlice].z0_id = Slice[iSlice].z0i_id; Slice[iSlice].ze_id = Slice[iSlice].zei_id;
	}
	// Get thickness
	getThk(nAtoms, Atoms, nSlice, Slice, Thk);
}

int cMT_Slicing_CPU::IsInThk(int iSlice){
	int iThk = -1;
	for(int i=0; i<Thk.n; i++){
		if(Thk.iSlice[i]==iSlice){
			iThk = i;
			break;
		}
	}

	return iThk;
}
