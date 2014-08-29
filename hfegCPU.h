#ifndef hfegCPU_H
#define hfegCPU_H

#include "hConstTypes.h"

class cfegCPU{
	private:
		sAtomTypesCPU AtomTypesCPU;
		double cl, cnl;	
		double ft, t, g2;
		void feg(double g, double &f, double &df);
	public:
		void SetAtomT(sAtomTypesCPU AtomTypesCPUi);
		void feg(int ng, double *g, double *f, double *df);
		cfegCPU();
		~cfegCPU();
};

#endif