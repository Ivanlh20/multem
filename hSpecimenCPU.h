#ifndef hSpecimenCPU_H
#define hSpecimenCPU_H

#include "hConstTypes.h"
#include "hRandGen.h"

class cSpecimen{
	private:	
		sAtoms *Atomsu;			// Undisplaced Atoms
		cRandGen RandGen;		// Ramdom generator
		int DimFP;				// Dimension frozen phonon
		int DistFP;				// Random distribution generator
		int SeedFP;				// Random seed generator
		int ZeroDefTyp;			// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
		double ZeroDefPlane;	// Zero defocus plane
		/*****************************************************************************/
		void QuickSortAtomsAlongz(sAtoms *Atoms, int left, int right);
		void setRandomSeed(unsigned long s, int iConf);
		void getDimCom(int Dim, double &bx, double &by, double &bz);
		/*****************************************************************************/
		int getAtomsInSlice(double z0, double ze, int nAtoms, sAtoms *Atoms, int &z0_id, int &ze_id);
		double getRint(int i0, int ie, sAtoms *Atoms, sAtomTypesCPU *AtomTypes, int Pos=1, double RBot=0, double Rmid=0, double RTop=0);
		void getSlicingParameters(double dz, int nAtoms, sAtoms *Atoms, sAtomTypesCPU *AtomTypes, int &nz, double &RTop, double &RBot);	
		void Slicing(double dz, int nAtoms, sAtoms *Atoms, sAtomTypesCPU *AtomTypes, int ZeroDefTyp, double ZeroDefPlane, int &nz, sSlice *&Slice, double &dzBackProp);
	public:
		int nAtomTypes;				// Number of atom types
		sAtomTypesCPU *AtomTypes;	// Atom types

		int nAtoms;					// Number of atoms in the specimen
		sAtoms *Atoms;				// Displaced Atoms

		int nz;						// Number of slices
		sSlice *Slice;				// Slicing procedure

		double dz;					// Slice thickness
		double dzBackProp;			//

		cSpecimen();
		~cSpecimen();
		void SetInputData(int nAtomsMi, double *AtomsMi, double dzi, int PotPari, double Vrli, int DimFPi, int DistFPi, int SeedFPi, int ZeroDefTypi, double ZeroDefPlanei);
		void MoveAtoms(int iConf);
};

#endif