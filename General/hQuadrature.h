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