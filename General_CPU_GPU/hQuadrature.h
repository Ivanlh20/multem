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

#ifndef hQuadrature_H
#define hQuadrature_H

#include "hConstTypes.h"

class cQuadrature{
	private:
		int i;
		double h;
		sQ1 Q;
		void CoefTanhSinh(int n, double xmin, double xmax);		// 0: int_-1^1 f(x) dx
		void CoefExpSinh(int n, double xmin, double xmax);		// 1: int_0^infty f(x) dx
		void CoefExpExp(int n, double xmin, double xmax);		// 2: int_0^infty f(x)exp(-x) dx
		void CoefSinhSinh(int n, double xmin, double xmax);		// 3: int_-infty^infty f(x) dx
		void CoefFourierTypeSin(int n, double ta);				// 4: int_0^infty f(x)sin(wx) dx
		void CoefFourierTypeCos(int n, double ta);				// 5: int_0^infty f(x)Cos(wx) dx
		void CoefGaussLegrendre(int n);							// 6: int_-1^1 f(x) dx
		void CoefGaussHermitezero2pinfty(int n);				// 7: int_0^infty f(x) Exp[-x^2] dx
		void CoefGaussHermiteninfty2pinfty(int n);				// 8: int_-infty^infty f(x) Exp[-x^2] dx
		void CoefGaussLaguerrex0zero2pinfty(int n);				// 9: int_0^infty f(x) Exp[-x] dx
		void CoefGaussLaguerrexnhzero2pinfty(int n);			// 10: int_0^infty f(x) Exp[-x]/Sqrt[x] dx
	public:
		cQuadrature();
		~cQuadrature();
		void ReadQuadrature(int qii, int nqi, sQ1 &Qo, double tai=4.0);

};

#endif