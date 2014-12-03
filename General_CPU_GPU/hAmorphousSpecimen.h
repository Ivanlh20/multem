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

#ifndef hAmorphousSpecimen_H
#define hAmorphousSpecimen_H

#include "hConstTypes.h"
#include "hRandGen.h"

class cAmorphousSpecimen{
	private:
		const int ng = 200;	// Number of trials 
		double lx;			// size in x
		double ly;			// size in y
		double lz;			// size in z
		double Z;			// Atomic number

		int na;				// Number of cube
		double *x;			// x-position
		double *y;			// y-position
		double *z;			// z-position
		double dmin;		// minimum separatin distance
		double dmin2;		// 
		double rmin;		// minimum radius
		double Vsmin;		// minimum volumen

		int ncxyz;			// Number of cube
		int ncxy;			// 
		int ncx;			// 
		int ncy;			// 
		int ncz;			// 
		double lc;			// size of the cube
		double ilc;			// 1.0/lc
		double Vcmin;		// minimum volumen of the cube
		int *c3d;			// array of indexes

		cRandGen RandGen;	// random generator

		inline void getLimits(int k, int kd, int nk, int &k0, int &ke);
		inline int CheckDistance(double &xn, double &yn, double &zn);
		inline int GeneratePoint(double &xn, double &yn, double &zn);
	public:
		cAmorphousSpecimen();
		~cAmorphousSpecimen();
		void freeMemory();
		void SetInputData(double lxi, double lyi, double lzi, double dmini, int Zi);	
		int CreateAmorphous();
		void Amorphous2Atoms(double *Atoms);
};

#endif