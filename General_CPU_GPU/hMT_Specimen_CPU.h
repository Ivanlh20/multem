/**
 *  This file is part of MULTEM.
 *  Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
 *
 *  MULTEM is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MULTEM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MULTEM.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef hMT_Specimen_CPU_H
#define hMT_Specimen_CPU_H

#include "hConstTypes.h"
#include "hRandGen.h"

class cMT_Specimen_CPU{
	private:
		cRandGen RandGen;			// Ramdom generator

		int nAtomsu;
		sAtoms *Atomsu;				// Undisplaced Atoms

		int nSliceu;				// Number of slices
		sSlice *Sliceu;				// Slicing procedure

		void QuickSortAtomsAlongz(sAtoms *&Atoms, int left, int right);
		void setRandomSeed(unsigned long s, int iConf);
		void getDimCom(int Dim, double &bx, double &by, double &bz);
		void getPlanes(int nAtoms, sAtoms *&Atoms, int &nPlanes, double *&Planes);
		void getnSlice(double z10, double z1i, double z20, double z2i, double dzi, int &nSlice, double &dz0, double &dze);
		int getBorderSlicing(double z0, double ze, double zi, double dzi, double &dzb);
		int getAtomsInSlice(double z0, double ze, int nAtoms, sAtoms *&Atoms, int &z0_id, int &ze_id);
		void Slicing(sMGP &MGP, double Rmax,int nAtoms, sAtoms *&Atoms, int &nSlice, sSlice *&Slice, double &z_BackProp);
		void Slicing(sMGP &MGP, double Rmax, int nSliceu, sSlice *Sliceu, int nAtoms, sAtoms *&Atoms, int &nSlice, sSlice *&Slice, double &z_BackProp);
		void Slicing(sMGP &MGP, int nPlanesu, double *Planesu, int nAtoms, sAtoms *&Atoms, int &nSlice, sSlice *&Slice, double &z_BackProp);
public:
		sMGP MGP;					// Multislice general parameters

		double Lzu;
		double Lztu;
		double Rmax;
		double z_BackProp;			// Back propagator

		int nAtomTypes;				// Number of atom types
		sAtomTypesCPU *AtomTypes;	// Atom types

		int nAtoms;					// Number of atoms in the specimen
		sAtoms *Atoms;				// Displaced Atoms

		int nSlice;					// Number of slices
		sSlice *Slice;				// Slicing procedure

		int nPlanesu;				// Number of Planes
		double *Planesu;			// Planes

		void freeMemory();
		cMT_Specimen_CPU();
		~cMT_Specimen_CPU();
		void SetInputData(sMGP &MGP_io, int nAtomsM_i, double *AtomsM_i);
		void MoveAtoms(int iConf);
		double get_dz(int iSlice);
};

#endif