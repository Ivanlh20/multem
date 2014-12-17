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

#ifndef hMT_Slicing_CPU_H
#define hMT_Slicing_CPU_H

#include <cstring>
#include "hConstTypes.h"
#include "hMT_MGP_CPU.h"
#include "hRandGen.h"
#include "hMT_General_CPU.h"

class cMT_Slicing_CPU{
	private:
		cMT_MGP_CPU *MT_MGP_CPU;	// Multislice general parameters

		double Rmax;				// Maximum interaction distance

		int nAtomsu;				// Number of atoms in the specimen
		sAtoms *Atomsu;				// Undisplaced Atoms

		int nAtoms;					// Number of atoms in the specimen
		sAtoms *Atoms;				// Displaced Atoms

		int binary_search(double z0, int iAtom0, int iAtome, sAtoms *&Atoms, int typ);
		void getThk(int nAtoms, sAtoms *Atoms, int nSlice, sSlice *Slice, sThk &Thk);
		void getPlanes(int nAtoms, sAtoms *&Atoms, int &nPlanesu, double *&Planesu);
		void get_nSlice(double z10, double z1i, double z20, double z2i, double dzi, int &nSlice, double &dz0, double &dze);
		void getAtomsInSlice(double z0, double ze, int nAtoms, sAtoms *&Atoms, int &z0_id, int &ze_id);
		void get_dzh_Planes(int nPlane, double *Plane, int iPlane, double &dzd, double &dzu);
	public:
		int nSliceu;				// Number of slices
		sSlice *Sliceu;				// Slicing procedure

		int nSlice;					// Number of slices
		sSlice *Slice;				// Slicing procedure

		int nPlanesu;				// Number of Planes
		double *Planesu;			// Planes

		sThk Thk;					//Thicknesses

		void freeMemory();
		cMT_Slicing_CPU();
		~cMT_Slicing_CPU();

		void SetInputData(cMT_MGP_CPU *MT_MGP_CPU_i, int nAtomsu_i, sAtoms *Atomsu_i, int nAtoms_i, sAtoms *Atoms_i, double Rmax_i);

		void Slicing_u(int nAtoms, sAtoms *&Atoms, int &nSlice, sSlice *&Slice, sThk &Thk);
		void Slicing_d(int nSliceu, sSlice *Sliceu, int nAtoms, sAtoms *&Atoms, int &nSlice, sSlice *&Slice, sThk &Thk);
		void Slicing_by_Planes(int nPlanesu, double *Planesu, int nAtoms, sAtoms *&Atoms, int &nSlice, sSlice *&Slice, sThk &Thk);
		int IsInThk(int iSlice);
};

#endif