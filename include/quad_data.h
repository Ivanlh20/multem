/*
* This file is part of Multem.
* Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
*
* Multem is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version of the License, or
* (at your option) any later version.
*
* Multem is distributed in the hope that it will be useful, 
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Multem. If not, see <http:// www.gnu.org/licenses/>.
*/

#pragma once

#include "quad_coef_1d.h"
#include "quad_coef_2d.h"

/* template definition */
namespace mt
{
	template <class T> class Quad_Data_T;
}

/* derived class */
namespace mt
{
	using Quad_Data = Quad_Data_T<dt_float64>;
}

namespace mt
{
	template <class T>
	class Quad_Data_T
	{
	public:
		using value_type = T;
		using size_type = dt_int32;

		Quad_Data_T() = default;

		template <class U, eDev Dev>
		Quad_Data_T(const eQuad_Typ& quad_typ, const dt_int32& n_quad, Quad_Coef_1d<U, Dev>& quad, T alpha = 0, T beta = 0, T a = 0, T b = 1);

		template <class U, eDev Dev>
		void operator()(const eQuad_Typ& quad_typ, const dt_int32& n_quad, Quad_Coef_1d<U, Dev>& quad, T alpha=0, T beta=0, T a=0, T b=1);

	private:
		// 1: int_-1^1 f(x) dx
		void nw_tanh_sinh_int_n1_p1(const dt_int32& nx, const T& x_min, const T& x_max, T* px, T* pw);

		// 2: int_0^infty f(x) dx
		void nw_exp_sinh_int_0_pinfty(const dt_int32& nx, const T& x_min, const T& x_max, T* px, T* pw);

		// 3: int_0^infty f(x)exp(-x) dx
		void nw_exp_exp_int_0_pinfty(const dt_int32& nx, const T& x_min, const T& x_max, T* px, T* pw);

		// 4: int_-infty^infty f(x) dx
		void nw_sinh_sinh_int_ninfty_pinfty(const dt_int32& nx, const T& x_min, const T& x_max, T* px, T* pw);

		//1. Ooura, T. & Mori, M. A robust double exponential formula for Fourier-type integrals. J. Comput. Appl. Math. 112, 229ï¿½241 (1999).
		// 5: int_0^infty f(x)sin(wx) dx
		void nw_fourier_sin_int_0_pinfty(const dt_int32& nx, const T& ta, T* px, T* pw);

		// 6: int_0^infty f(x)Cos(wx) dx
		void nw_fourier_cos_int_0_pinfty(const dt_int32& nx, const T& ta, T* px, T* pw);

		// 7: int_-1^1 f(x) dx
		void nw_gauss_legendre_int_n1_p1(const dt_int32& nx, T* px, T* pw);

		// 8: int_-infty^infty f(x) x^0 Exp[-x^2] dx
		void nw_gauss_hermite_x0_int_ninfty_pinfty(const dt_int32& nx, T* px, T* pw);

		// 9: int_-infty^infty f(x) |x|^1 Exp[-x^2] dx
		void nw_gauss_hermite_x1_int_ninfty_pinfty(const dt_int32& nx, T* px, T* pw);

		// 9: int_-infty^infty f(x) |x|^2 Exp[-x^2] dx
		void nw_gauss_hermite_x2_int_ninfty_pinfty(const dt_int32& nx, T* px, T* pw);

		// 11: int_0^infty f(x) x^0 Exp[-x] dx
		void nw_gauss_laguerre_x0_int_0_pinfty(const dt_int32& nx, T* px, T* pw);

		// 12: int_0^infty f(x) x^1 Exp[-x] dx
		void nw_gauss_laguerre_x1_int_0_pinfty(const dt_int32& nx, T* px, T* pw);

		// 13: int_0^infty f(x) x^2 Exp[-x] dx
		void nw_gauss_laguerre_x2_int_0_pinfty(const dt_int32& nx, T* px, T* pw);

		// 14: int_0^infty f(x) Exp[-x]/Sqrt[x] dx
		void nw_gauss_laguerre_xi2_int_0_pinfty(const dt_int32& nx, T* px, T* pw);
 
		/***************************************************************************************/
		// in, dt_int32 KIND, the rule.
		// 1, legendre, (a, b) 1.0
		// 2, chebyshev, (a, b) ((b-x)*(x-a))^(-0.5)
		// 3, gegenbauer, (a, b) ((b-x)*(x-a))^alpha
		// 4, jacobi, (a, b) (b-x)^alpha*(x-a)^beta
		// 5, generalized laguerre, (a, inf) (x-a)^alpha*exp(-b*(x-a))
		// 6, generalized hermite, (-inf, inf) |x-a|^alpha*exp(-b*(x-a)^2)
		// 7, exponential, (a, b) |x-(a+b)/2.0|^alpha
		// 8, rational, (a, inf) (x-a)^alpha*(x+b)^beta

		void cgqf(dt_int32 kind, dt_int32 nt, T alpha, T beta, T a, T b, T t[], T wts[]);

		void cdgqf(dt_int32 kind, dt_int32 nt, T alpha, T beta, T t[], T wts[]);

		T class_matrix(dt_int32 kind, dt_int32 m, T alpha, T beta, T aj[], T bj[]);

		void imtqlx(dt_int32 n, T d[], T e[], T z[]);

		void parchk(dt_int32 kind, dt_int32 m, T alpha, T beta);

		void scqf(dt_int32 nt, T t[], dt_int32 mlt[], T wts[], dt_int32 nwts, dt_int32 ndx[], 
			T swts[], T st[], dt_int32 kind, T alpha, T beta, T a, T b);

		void sgqf(dt_int32 nt, T aj[], T bj[], T zemu, T t[], T wts[]);

		T r8_epsilon();

		T r8_huge();

		T r8_sign(T x);
	};
}

#include "../src/quad_data.inl"