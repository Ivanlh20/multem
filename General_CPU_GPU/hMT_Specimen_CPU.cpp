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

	nAtomsu = 0;
	delete [] Atomsu; Atomsu = 0;

	nAtoms = 0;
	delete [] Atoms; Atoms = 0;

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

	nAtomsu = 0;
	Atomsu = 0;

	nAtoms = 0;
	Atoms = 0;

	nMT_AtomTypes = 0;
	MT_AtomTypes_CPU = 0;
}

// Destructor
cMT_Specimen_CPU::~cMT_Specimen_CPU()
{
	freeMemory();
}

// random number generator Seed
void cMT_Specimen_CPU::setRandomSeed(unsigned long s, int iConf)
{	
	int randiu;
	RandGen.seed(s);
	for(int iSlice=0; iSlice<iConf; iSlice++)
	{
		randiu = RandGen.randiu();
	}
	RandGen.seed(randiu);
}

// Sort atoms along z-axis
void cMT_Specimen_CPU::QuickSortAtomsAlongz(sAtoms *&Atoms, int left, int right)
{
	int iSlice = left, j = right;
	double pivot = Atoms[(left + right)>>1].z;
	sAtoms Atomst;

	/* partition */
	while (iSlice <= j)
	{
		while (Atoms[iSlice].z < pivot)
		{
			iSlice++;
		}

		while (Atoms[j].z > pivot)
		{
			j--;
		}

		if(iSlice <= j)
		{
			memcpy(&Atomst, &Atoms[iSlice], cSizeofAtoms);
			memcpy(&Atoms[iSlice], &Atoms[j], cSizeofAtoms);
			memcpy(&Atoms[j], &Atomst, cSizeofAtoms);
			iSlice++;
			j--;
		}
	};

	/* recursion */
	if(left < j)
	{
		QuickSortAtomsAlongz(Atoms, left, j);
	}

	if(iSlice < right)
	{
		QuickSortAtomsAlongz(Atoms, iSlice, right);
	}
}

// get dimension components
void cMT_Specimen_CPU::getDimCom(int Dim, double &bx, double &by, double &bz)
{
	bx = by = bz = 0.0;
	switch(Dim)
	{
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

// Set atoms: Copy input Atoms to the new format Atoms, count number of atoms, Ascending sort by z, Set relative atomic number position, Get maximum interaction distance
void cMT_Specimen_CPU::SetInputData(cMT_MGP_CPU *MT_MGP_CPU_io, int nAtomsM_i, double *AtomsM_i, double dRmin)
{
	freeMemory();

	MT_MGP_CPU = MT_MGP_CPU_io;
	nMT_AtomTypes = stNAE;
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

	if(Lzu<MT_MGP_CPU->dz)
	{
		MT_MGP_CPU->MulOrder = 1;
		MT_MGP_CPU->dz = Lzu;
		if(MT_MGP_CPU->ApproxModel==1)
			MT_MGP_CPU->ApproxModel = 3;
	}
	// get Zero defocus plane
	switch(MT_MGP_CPU->ZeroDefTyp)
	{
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

	// Copy Atomsu to Atoms
	nAtoms = nAtomsu;
	Atoms = new sAtoms[nAtoms];
	memcpy(Atoms, Atomsu, nAtoms*cSizeofAtoms);

	cMT_Slicing_CPU::SetInputData(MT_MGP_CPU, nAtomsu, Atomsu, nAtoms, Atoms, Rmax);
}

// Move atoms (ramdom distribution will be included in the future)
void cMT_Specimen_CPU::MoveAtoms(int iConf)
{
	if(iConf<=0) return;

	double sigmax, sigmay, sigmaz;
	double bx, by, bz;

	// Get dimension components
	getDimCom(MT_MGP_CPU->DimFP, bx, by, bz);
	setRandomSeed(MT_MGP_CPU->SeedFP, iConf);
	for(int iAtoms = 0; iAtoms<nAtoms; iAtoms++)
	{
		Atoms[iAtoms].Z = Atomsu[iAtoms].Z;
		sigmax = sigmay = sigmaz = Atomsu[iAtoms].sigma;
		Atoms[iAtoms].x = Atomsu[iAtoms].x + bx*sigmax*RandGen.randn();
		Atoms[iAtoms].y = Atomsu[iAtoms].y + by*sigmay*RandGen.randn();
		Atoms[iAtoms].z = Atomsu[iAtoms].z + bz*sigmaz*RandGen.randn();
		Atoms[iAtoms].sigma = 0.0;
		Atoms[iAtoms].occ = Atomsu[iAtoms].occ ;
	}
	if(MT_MGP_CPU->ApproxModel==1)
	{
		// Ascending sort by z
		QuickSortAtomsAlongz(Atoms, 0, nAtoms-1);	
		// Slicing procedure
		Slicing_d(nSliceu, Sliceu, nAtoms, Atoms, nSlice, Slice, Thk);
	}
}

// Get dz
double cMT_Specimen_CPU::get_dz(int iSlice)
{
	double dz;
	if(nSlice==1)
		dz = Rmax;
	else
		dz = (iSlice<nSlice)?(Slice[iSlice].ze - Slice[iSlice].z0)/cos(MT_MGP_CPU->theta):0.0;
	return dz;
}