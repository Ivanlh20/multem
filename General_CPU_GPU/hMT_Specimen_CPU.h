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

#ifndef hMT_Specimen_CPU_H
#define hMT_Specimen_CPU_H

#include <cstring>
#include "hConstTypes.h"
#include "hMT_MGP_CPU.h"
#include "hRandGen.h"
#include "hMT_General_CPU.h"
#include "hMT_Slicing_CPU.h"

class cMT_Specimen_CPU: public cMT_Slicing_CPU{
	private:
		cRandGen RandGen;						// Ramdom generator

		void QuickSortAtomsAlongz(sAtoms *&Atoms, int left, int right);
		void setRandomSeed(unsigned long s, int iConf);
		void getDimCom(int Dim, double &bx, double &by, double &bz);
	protected:
		int nAtomsu;
		sAtoms *Atomsu;							// Undisplaced Atoms

	public:
		cMT_MGP_CPU *MT_MGP_CPU;				// Multislice general parameters

		double sigma_min;
		double sigma_max;

		double Lzu;
		double Lztu;
		double Rmax;

		int nMT_AtomTypes;						// Number of atom types
		cMT_AtomTypes_CPU *MT_AtomTypes_CPU;	// Atom types

		int nAtoms;								// Number of atoms in the specimen
		sAtoms *Atoms;							// Displaced Atoms

		void freeMemory();
		cMT_Specimen_CPU();
		~cMT_Specimen_CPU();

		void SetInputData(cMT_MGP_CPU *MT_MGP_CPU_io, int nAtomsM_i, double *AtomsM_i, double dRmin=0);
		void MoveAtoms(int iConf);
		double get_dz(int iSlice);
};

#endif