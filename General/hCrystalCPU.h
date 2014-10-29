#ifndef hCrystalCPU_H
#define hCrystalCPU_H

#include "hConstTypes.h"

class cCrystalCPU{
	private:
		int na;
		int nb;
		int nc;

		double a;
		double b;
		double c;

		int nuLayer;
		sAtomsGroup *uLayer;
		sAtomsGroup *Layers;

		void uLayer2Layer(int na, int nb, double a, double b, sAtomsGroup &uLayer, sAtomsGroup &Layer);
		int StackLayers(int na, int nb, double a, double b, int nuLayer, sAtomsGroup *&uLayer, sAtomsGroup *&Layers);
	public:
		void freeMemory();
		~cCrystalCPU();
		cCrystalCPU();
		void SetInputData(int nai, int nbi, int nci, double ai, double bi, double ci, int nuLayeri, sAtomsGroup *uLayeri, int &nAtoms);
		void Create3DCrystal(int nAtomsM, double *&AtomsM);
};

#endif
