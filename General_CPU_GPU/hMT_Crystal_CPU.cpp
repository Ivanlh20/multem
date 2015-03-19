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

#include <cstring>
#include "hConstTypes.h"
#include "hMT_Crystal_CPU.h"

void cMT_Crystal_CPU::freeMemory()
{
	na = 0;
	nb = 0;
	nc = 0;

	a = 0;
	b = 0;
	c = 0;

	for(int i=0; i<nuLayer; i++)
	{
		uLayer[i].nAtoms = 0;
		delete [] uLayer[i].Atoms; uLayer[i].Atoms = 0;
		Layers[i].nAtoms = 0;
		delete [] Layers[i].Atoms; Layers[i].Atoms = 0;
	}

	nuLayer = 0;
	delete [] uLayer; uLayer = 0;
	delete [] Layers; Layers = 0;	
}

cMT_Crystal_CPU::~cMT_Crystal_CPU()
{
	freeMemory();
}

cMT_Crystal_CPU::cMT_Crystal_CPU()
{
	na = 0;
	nb = 0;
	nc = 0;

	a = 0;
	b = 0;
	c = 0;

	nuLayer = 0;
	uLayer = 0;
	Layers = 0;	
}

void cMT_Crystal_CPU::uLayer2Layer(int na, int nb, double a, double b, sAtomsGroup &uLayer, sAtomsGroup &Layer)
{
 int i, j, k, l;
	double xmin, ymin, xmax, ymax;
	double x, y;

	xmin = 0.0 - 1e-05; xmax = na*a + 1e-05;
	ymin = 0.0 - 1e-05; ymax = nb*b + 1e-05;

	l = 0;
	for(j=0; j<=nb; j++)
		for(i=0; i<=na; i++)	
			for(k=0; k<uLayer.nAtoms; k++)
			{
				x = i*a + uLayer.Atoms[k].x;
				y = j*b + uLayer.Atoms[k].y;				
				if((xmin<=x)&&(x<=xmax)&&(ymin<=y)&&(y<=ymax))
				{
					Layer.Atoms[l].x = x;
					Layer.Atoms[l].y = y;
					Layer.Atoms[l].z = uLayer.Atoms[k].z;
					Layer.Atoms[l].Z = uLayer.Atoms[k].Z;
					Layer.Atoms[l].sigma = uLayer.Atoms[k].sigma;
					Layer.Atoms[l].occ = uLayer.Atoms[k].occ;
					l++;
				}
			}
	Layer.nAtoms = l;
}

int cMT_Crystal_CPU::StackLayers(int na, int nb, double a, double b, int nuLayer, sAtomsGroup *&uLayer, sAtomsGroup *&Layers)
{
	int nAtomsLayers = 0;
	for(int i=0; i<nuLayer; i++)
	{
		uLayer2Layer(na, nb, a, b, uLayer[i], Layers[i]);
		nAtomsLayers += Layers[i].nAtoms;
	}
	return nAtomsLayers;
}

void cMT_Crystal_CPU::SetInputData(int nai, int nbi, int nci, double ai, double bi, double ci, int nuLayeri, sAtomsGroup *uLayeri, int &nAtoms)
{
	freeMemory();

	na = nai;
	nb = nbi;
	nc = nci;

	a = ai;
	b = bi;
	c = ci;

	nuLayer = nuLayeri;

	uLayer = new sAtomsGroup[nuLayer];
	for(int i=0; i<nuLayer; i++)
	{
		uLayer[i].nAtoms = uLayeri[i].nAtoms;
		uLayer[i].Atoms = new sAtoms[uLayer[i].nAtoms];
		memcpy(uLayer[i].Atoms, uLayeri[i].Atoms, uLayer[i].nAtoms*sizeof(sAtoms));
	}

	Layers = new sAtomsGroup[nuLayer]; 	
	for(int i=0; i<nuLayer; i++)
	{
		Layers[i].nAtoms = uLayer[i].nAtoms*(na+1)*(nb+1);
		Layers[i].Atoms = new sAtoms[Layers[i].nAtoms];
	}

	int nAtomsLayers;
	nAtomsLayers = StackLayers(na, nb, a, b, nuLayer, uLayer, Layers);
	nAtoms = nc*nAtomsLayers + Layers[0].nAtoms;
}

void cMT_Crystal_CPU::Create3DCrystal(int nAtomsM, double *&AtomsM)
{
	int i, j, k, l;

	l = 0;
	for(k=0; k<nc; k++)
	{
		for(i=0; i<nuLayer; i++)
			for(j=0; j<Layers[i].nAtoms; j++)
			{
				AtomsM[0*nAtomsM+l] = Layers[i].Atoms[j].x;
				AtomsM[1*nAtomsM+l] = Layers[i].Atoms[j].y;
				AtomsM[2*nAtomsM+l] = Layers[i].Atoms[j].z + c*k;
				AtomsM[3*nAtomsM+l] = Layers[i].Atoms[j].Z;
				AtomsM[4*nAtomsM+l] = Layers[i].Atoms[j].sigma;
				AtomsM[5*nAtomsM+l] = Layers[i].Atoms[j].occ;
				l++;
			}
	}

	// Last layer
	for(j=0; j<Layers[0].nAtoms; j++)
	{
		AtomsM[0*nAtomsM+l] = Layers[0].Atoms[j].x;
		AtomsM[1*nAtomsM+l] = Layers[0].Atoms[j].y;
		AtomsM[2*nAtomsM+l] = Layers[0].Atoms[j].z + c*nc;
		AtomsM[3*nAtomsM+l] = Layers[0].Atoms[j].Z;
		AtomsM[4*nAtomsM+l] = Layers[0].Atoms[j].sigma;
		AtomsM[5*nAtomsM+l] = Layers[0].Atoms[j].occ;
		l++;
	}	
}