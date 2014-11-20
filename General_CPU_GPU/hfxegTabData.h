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

#ifndef hfxegTabData_H
#define hfxegTabData_H

class cfxegTabData{
	private:
		int nmax;
		int ng;
		double *g;
		double *g2;
		double *feg;
		double *fxg;
		void fxegActaCrys(int Z);
		void fxegRez(int Z);
		void fxegKirkland(int Z);
		void fxegLobato(int Z);
	public:
		cfxegTabData();
		~cfxegTabData();
		void ReadTabData(int Zi, int Typi, int dni, int &no, double *go, double *g2o, double *fxo, double *feo);
};

#endif