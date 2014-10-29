#ifndef hrhorCPU_H
#define hrhorCPU_H

#include "hConstTypes.h"

class crhorCPU{
	private:
		int PotPar;
		sAtomTypesCPU AtomTypesCPU;
		double cl, cnl;	
		double ft, t, ir, r2;
		void rhor(double r, double &f, double &df);
	public:
		void SetAtomT(int PotPari, sAtomTypesCPU AtomTypesCPUi);
		void rhor(int nr, double *r, double *f, double *df);
		crhorCPU();
		~crhorCPU();
};

#endif