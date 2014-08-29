#ifndef hPotentialCPU_H
#define hPotentialCPU_H

#include "hConstTypes.h"

class cPotentialCPU{
	private:
		sAtomTypesCPU AtomTypes;
		double cl, cnl, icnl;	
		double ft, t, ir, r2, g2;
		double k0, k1;
		double sigma;

		int nQ1;
		sQ1 Q1;

		void Pot1Da(double &r, double &f, double &df);
		void Pot1Dn(double &r, double &f, double &df);

		void Pot2Da(double &r, double &f, double &df);
		void Pot2Dn(double &r, double &f, double &df);

		void Pot3Da(double &r, double &f, double &df);
		void Pot3Dn(double &r, double &f, double &df);
	public:
		cPotentialCPU();
		~cPotentialCPU();

		void SetSigma(double sigmai);
		void SetAtomTypes(sAtomTypesCPU AtomTypesi);

		void Vr(int IntType, int Dim, double r, double &f, double &df);
		void Vr(int IntType, int Dim, int nr, double *r, double *f, double *df);
		void Vr(int nAtomTypes, sAtomTypesCPU *AtomTypes, int IntType, int Dim, int nr, double *r, double factor, double *f, double *df);

		double AtomicRadius_rms(int Dim);
		double AtomicRadius_Cutoff(int Dim, double Vrl);
};

#endif