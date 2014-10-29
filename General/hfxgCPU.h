#ifndef hfxgCPU_H
#define hfxgCPU_H

#include "hConstTypes.h"

class cfxgCPU{
	private:
		int PotPar;
		sAtomTypesCPU AtomTypesCPU;
		double cl, cnl;	
		double ft, t, g2;
		void fxg(double g, double &f, double &df);
	public:
		void SetAtomT(int PotPari, sAtomTypesCPU AtomTypesCPUi);
		void fxg(int ng, double *g, double *f, double *df);
		cfxgCPU();
		~cfxgCPU();
};

#endif