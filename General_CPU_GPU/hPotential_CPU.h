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

#ifndef hPotential_CPU_H
#define hPotential_CPU_H

#include "hConstTypes.h"

class cPotential_CPU{
	private:
		int PotPar;
		sAtomTypesCPU AtomTypes;
		double *cl, *cnl;	
		double ft, t, g2;
		double k0, k1;
		double sigma;

		int nQ1;
		sQ1 Q1;

		void Pot1Da(double &r, double *cl, double *cnl, double &f, double &df);
		void Pot1Dn(double &r, double *cl, double *cnl, double &f, double &df);

		void Pot2Da(double &r, double *cl, double *cnl, double &f, double &df);
		void Pot2Dn(double &r, double *cl, double *cnl, double &f, double &df);

		void Pot3Da(double &r, double *cl, double *cnl, double &f, double &df);
		void Pot3Dn(double &r, double *cl, double *cnl, double &f, double &df);
	public:
		cPotential_CPU();
		~cPotential_CPU();

		void SetSigma(double sigmai);
		void SetAtomTypes(int PotPari, sAtomTypesCPU AtomTypesi);

		void Vr(int IntType, int Dim, double r, double &f, double &df);
		void Vr(int IntType, int Dim, int nr, double *r, double *f, double *df);
		void Vr(int PotPar, int nAtoms, sAtoms *Atoms, int nAtomTypes, sAtomTypesCPU *AtomTypes, int IntType, int Dim, int nr, double *r, double factor, double *f, double *df);

		double AtomicRadius_rms(int Dim);
		double AtomicRadius_Cutoff(int Dim, double Vrl);
};

#endif