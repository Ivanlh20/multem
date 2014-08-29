#ifndef hAtomicData_H
#define hAtomicData_H

#include "hConstTypes.h"

class cAtomicData{
	private:
		sAtDa AtDa[NE];
		sfep fep[NE];
		void AtomicData();
		void AtomicParameterization();	
	public:
		void ReadAtomicData(int nAtomTypes, sAtomTypesCPU *AtomTypes, double Vrl=0);
};

#endif