/**
 *  This file is part of MULTEM.
 *  Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
 *
 *  MULTEM is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MULTEM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MULTEM.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef hMT_Crystal_CPU_H
#define hMT_Crystal_CPU_H

#include "hConstTypes.h"

class cMT_Crystal_CPU{
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
		~cMT_Crystal_CPU();
		cMT_Crystal_CPU();
		void SetInputData(int nai, int nbi, int nci, double ai, double bi, double ci, int nuLayeri, sAtomsGroup *uLayeri, int &nAtoms);
		void Create3DCrystal(int nAtomsM, double *&AtomsM);
};

#endif
