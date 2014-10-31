#ifndef hSpecimenCPU_H
#define hSpecimenCPU_H

#include "hConstTypes.h"
#include "hRandGen.h"

class cSpecimenCPU{
	private:
		int nAtomsu;
		sAtoms *Atomsu;				// Undisplaced Atoms
		cRandGen RandGen;			// Ramdom generator
		void QuickSortAtomsAlongz(sAtoms *&Atoms, int left, int right);
		void setRandomSeed(unsigned long s, int iConf);
		void getDimCom(int Dim, double &bx, double &by, double &bz);
		void getPlanes(int nAtoms, sAtoms *&Atoms, int &nPlanes, double *&Planes);
		int getAtomsInSlice(double z0, double ze, int nAtoms, sAtoms *&Atoms, int &z0_id, int &ze_id);
		void Slicing(int ApproxModel, double dz, double Rmax, int nAtoms, sAtoms *&Atoms, int ZeroDefTyp, double ZeroDefPlane, int &nSlice, sSlice *&Slice, double &z_BackProp);
	public:
		sMGP MGP;					// Multislice general parameters

		double Lzu;
		double Lztu;
		double Rmax;
		double z_BackProp;			// Back propagator

		int nPlanes;
		double *Planes;

		int nAtomTypes;				// Number of atom types
		sAtomTypesCPU *AtomTypes;	// Atom types

		int nAtoms;					// Number of atoms in the specimen
		sAtoms *Atoms;				// Displaced Atoms

		int nSlice;					// Number of slices
		sSlice *Slice;				// Slicing procedure

		void freeMemory();
		cSpecimenCPU();
		~cSpecimenCPU();
		void SetInputData(sMGP &MGP_io, int nAtomsM_i, double *AtomsM_i);
		void MoveAtoms(int iConf);
};

#endif