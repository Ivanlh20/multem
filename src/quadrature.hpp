/*
 * This file is part of Multem.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * Multem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
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

#ifndef QUADRATURE_H
#define QUADRATURE_H

#ifdef _MSC_VER
#pragma once
#endif // _MSC_VER

#include "types.cuh"
#include "traits.cuh"

namespace mt
{
	class Quadrature{
		public:
			template <class TQ1>
			void operator()(int Quad_Type, int nQuad, TQ1 &Q, double alpha=0, double beta=0, double a=0, double b=1)
			{
				using value_type = Value_type<TQ1>;

				if (nQuad <= 1)
				{
					Q.resize(1);
					Q.x[0] = 0;
					Q.w[0] = 1;
					return;
				}

				value_type x_min, x_max, r_max = 225;
				double tai = 4.0;

				cQ.resize(nQuad);

				switch(Quad_Type)
				{
					case 0: // 0: int_-1^1 f(x) dx
						x_min = -1.0+Epsilon<float>::eps; 
						x_max = 1.0-Epsilon<float>::eps;
						nw_TanhSinh(cQ, x_min, x_max);
						break;
					case 1: // 1: int_0^infty f(x) dx
						x_min = Epsilon<float>::eps;
						x_max = r_max;
						nw_ExpSinh(cQ, x_min, x_max);
						break;
					case 2: // 2: int_0^infty f(x)exp(-x) dx
						x_min = Epsilon<float>::eps; 
						x_max = r_max;
						nw_ExpExp(cQ, x_min, x_max);
						break;
					case 3: // 3: int_-infty^infty f(x) dx
						x_min = -r_max; 
						x_max = r_max;
						nw_SinhSinh(cQ, x_min, x_max);
						break;
					case 4: // 4: int_0^infty f(x)sin(wx) dx
						nw_Fourier_Sin(cQ, tai); 	// tai = 2.6
						break;
					case 5: // 5: int_0^infty f(x)Cos(wx) dx
						nw_Fourier_Cos(cQ, tai); 	// tai = 2.6
						break;
					case 6: // 6: int_-1^1 f(x) dx
						nw_Gauss_Legendre_int_n1_p1(cQ);
						break;
					case 7: // 7: int_-infty^infty f(x) x^0 Exp[-x^2] dx
						nw_Gauss_Hermite_x0_int_ninfty_pinfty(cQ);
						break;
					case 8: // 8: int_-infty^infty f(x) |x|^1 Exp[-x^2] dx
						nw_Gauss_Hermite_x1_int_ninfty_pinfty(cQ);
						break;
					case 9: // 9: int_-infty^infty f(x) |x|^2 Exp[-x^2] dx
						nw_Gauss_Hermite_x2_int_ninfty_pinfty(cQ);
						break;
					case 10: // 10: int_0^infty f(x) x^0 Exp[-x] dx
						nw_Gauss_Laguerre_x0_int_0_pinfty(cQ);
						break;
					case 11: // 11: int_0^infty f(x) x^1 Exp[-x] dx
						nw_Gauss_Laguerre_x1_int_0_pinfty(cQ);
						break;
					case 12: // 12: int_0^infty f(x) x^2 Exp[-x] dx
						nw_Gauss_Laguerre_x2_int_0_pinfty(cQ);
						break;
					case 13: // 12: int_0^infty f(x) Exp[-x]/Sqrt[x] dx
						nw_Gauss_Laguerre_xi2_int_0_pinfty(cQ);
						break;
					case 21: //    1, Legendre,             (a,b)       1.0
						cgqf(cQ.size(), 1, alpha, beta, a, b, cQ.x.data(), cQ.w.data());
						break;
					case 22: //    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
						cgqf(cQ.size(), 2, alpha, beta, a, b, cQ.x.data(), cQ.w.data());
						break;
					case 23: //    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
						cgqf(cQ.size(), 3, alpha, beta, a, b, cQ.x.data(), cQ.w.data());
						break;
					case 24: //    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
						cgqf(cQ.size(), 4, alpha, beta, a, b, cQ.x.data(), cQ.w.data());
						break;
					case 25: //    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
						cgqf(cQ.size(), 5, alpha, beta, a, b, cQ.x.data(), cQ.w.data());
						break;
					case 26: //    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
						cgqf(cQ.size(), 6, alpha, beta, a, b, cQ.x.data(), cQ.w.data());
						break;
					case 27: //    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
						cgqf(cQ.size(), 7, alpha, beta, a, b, cQ.x.data(), cQ.w.data());
						break;
					case 28: //    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
						cgqf(cQ.size(), 8, alpha, beta, a, b, cQ.x.data(), cQ.w.data());
						break;
				}

				Q.assign(cQ);
			}

		private:
			// 0: int_-1^1 f(x) dx
			void nw_TanhSinh(Q1<double, e_host> &cQ, const double &x_min, const double &x_max)
			{
				auto Aprox = [](double x)->double{ return asinh(atanh(x)/c_i2Pi); };
				double tmin = Aprox(x_min);
				double tmax = Aprox(x_max);

				double ti, tti, xi, wi;
				double h = (tmax-tmin)/(cQ.size()-1);

				for(auto i = 0; i < cQ.size(); i++)
				{
					ti = tmin + i*h;
					tti = c_i2Pi*sinh(ti);
					xi = tanh(tti);
					wi = h*c_i2Pi*cosh(ti)/pow(cosh(tti), 2);
					cQ.x[i] = xi;
					cQ.w[i] = wi;
				}
			}

			// 1: int_0^infty f(x) dx
			void nw_ExpSinh(Q1<double, e_host> &cQ, const double &x_min, const double &x_max)
			{
				auto getLimit = [](double x)->double{ return asinh(log(x)/c_i2Pi); };
				double t_min = getLimit(x_min);
				double t_max = getLimit(x_max);
				int N = cQ.size();
				double h = (t_max-t_min)/(N-1);
				double f = 0.5;

				for(auto ik = 0; ik < N; ik++)
				{
					auto t_ik = t_min + ik*h;
					auto x_i = exp(c_i2Pi*sinh(t_ik));
					auto w_i = f*h*c_i2Pi*cosh(t_ik)*x_i;
					cQ.x[ik] = x_i;
					cQ.w[ik] = w_i;
					f = (ik==N-1)?0.5:1.0;
				}
			}

			// 2: int_0^infty f(x)exp(-x) dx
			void nw_ExpExp(Q1<double, e_host> &cQ, const double &x_min, const double &x_max)
			{
				auto getLimit = [](double x, double t, int n)->double
				{
					auto t_s = log(x);
					for (auto ik = 0; ik < n; ik++)
					{
						auto exp_t = exp(-t);
						t = t - (t - exp_t - t_s)/(1.0 + exp_t);
					}
					return t;
				};

				double t_min = -log(-log(x_min));
				t_min = getLimit(x_min, t_min, 5);
				double t_max = log(x_max);
				t_max = getLimit(x_max, t_max, 5);

				int N = cQ.size();
				double h = (t_max-t_min)/(N-1);
				double f = 0.5;

				for(auto ik = 0; ik < N; ik++)
				{
					auto t_ik = t_min + ik*h;
					auto exp_t_ik = exp(-t_ik);
					auto x_i = exp(t_ik-exp_t_ik);
					auto w_i = f*h*x_i*(1.0+exp_t_ik);
					cQ.x[ik] = x_i;
					cQ.w[ik] = w_i;
					f = (ik==N-1)?0.5:1.0;
				}
			}

			// 3: int_-infty^infty f(x) dx
			void nw_SinhSinh(Q1<double, e_host> &cQ, const double &x_min, const double &x_max)
			{
				auto getLimit = [](double x)->double{ return asinh(asin(x)/c_i2Pi); };
				double tmin = getLimit(x_min);
				double tmax = getLimit(x_max);

				double ti, tti, xi, wi;
				double h = (tmax-tmin)/(cQ.size()-1);

				for(auto i = 0; i < cQ.size(); i++)
				{
					ti = tmin + i*h;
					tti = c_i2Pi*sinh(ti);
					xi = sinh(tti);
					wi = h*c_i2Pi*cosh(ti)*cosh(tti);
					cQ.x[i] = xi;
					cQ.w[i] = wi;
				}
			}

			// 4: int_0^infty f(x)sin(wx) dx
			void nw_Fourier_Sin(Q1<double, e_host> &cQ, const double &ta)
			{
				double M, h, k, xi, wi;
				double ti, ut, phi, dphi;

				h = 2.0*ta/cQ.size(); 
				M = c_Pi/h; k = 6;

				for(auto i = 0; i < cQ.size(); i++)
				{
					ti = -ta + (i+1)*h;
					if(ti == 0 )
					{
						phi = 1.0/k; dphi = 0.5;
						xi = M*phi; 
						wi = c_Pi*dphi*sin(M*phi);
					}
					else
					{
						ut = 1.0-exp(-k*sinh(ti));
						phi = ti/ut; 
						dphi = (1.0+(1.0+k*ti*cosh(ti))*(ut-1))/(ut*ut);
						xi = M*phi; 
						wi = c_Pi*dphi*sin(M*phi);
					}
					cQ.x[i] = xi;
					cQ.w[i] = wi;
				}
			}

			// 5: int_0^infty f(x)Cos(wx) dx
			void nw_Fourier_Cos(Q1<double, e_host> &cQ, const double &ta)
			{
				double M, h, k, xi, wi;
				double ti, ut, phi, dphi;

				h = 2.0*ta/cQ.size(); 
				M = c_Pi/h; k = 6.0;

				for(auto i = 0; i < cQ.size(); i++)
				{
					ti = -ta + (i+0.5)*h;
					if(ti == 0 )
					{
						phi = 1.0/k; 
						dphi = 0.5;
						xi = M*phi; 
						wi = c_Pi*dphi*cos(M*phi);
					}
					else
					{
						ut = 1.0-exp(-k*sinh(ti));
						phi = ti/ut; 
						dphi = (1.0+(1.0+k*ti*cosh(ti))*(ut-1))/(ut*ut);
						xi = M*phi; 
						wi = c_Pi*dphi*cos(M*phi);
					}
					cQ.x[i] = xi;
					cQ.w[i] = wi;
				}
			}

			// 6: int_-1^1 f(x) dx
			void nw_Gauss_Legendre_int_n1_p1(Q1<double, e_host> &cQ)
			{
				//    Input, int KIND, the rule.
				//    1, Legendre,             (a,b)       1.0
				//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
				//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
				//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
				//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
				//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
				//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
				//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

				int n_quad = cQ.size();
				int kind = 1;
				double alpha = 0;
				double beta = 0;
				double a = -1.0;
				double b = 1.0;

				cgqf(n_quad, kind, alpha, beta, a, b, cQ.x.data(), cQ.w.data());
			}

			// 7: int_-infty^infty f(x) x^0 Exp[-x^2] dx
			void nw_Gauss_Hermite_x0_int_ninfty_pinfty(Q1<double, e_host> &cQ)
			{
				//    Input, int KIND, the rule.
				//    1, Legendre,             (a,b)       1.0
				//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
				//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
				//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
				//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
				//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
				//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
				//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

				int n_quad = cQ.size();
				int kind = 6;
				double alpha = 0;
				double beta = 0;
				double a = 0.0;
				double b = 1.0;

				cgqf(n_quad, kind, alpha, beta, a, b, cQ.x.data(), cQ.w.data());
			}

			// 8: int_-infty^infty f(x) |x|^1 Exp[-x^2] dx
			void nw_Gauss_Hermite_x1_int_ninfty_pinfty(Q1<double, e_host> &cQ)
			{
				//    Input, int KIND, the rule.
				//    1, Legendre,             (a,b)       1.0
				//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
				//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
				//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
				//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
				//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
				//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
				//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

				int n_quad = cQ.size();
				int kind = 6;
				double alpha = 1;
				double beta = 0;
				double a = 0.0;
				double b = 1.0;

				cgqf(n_quad, kind, alpha, beta, a, b, cQ.x.data(), cQ.w.data());
			}

			// 9: int_-infty^infty f(x) |x|^2 Exp[-x^2] dx
			void nw_Gauss_Hermite_x2_int_ninfty_pinfty(Q1<double, e_host> &cQ)
			{
				//    Input, int KIND, the rule.
				//    1, Legendre,             (a,b)       1.0
				//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
				//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
				//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
				//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
				//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
				//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
				//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

				int n_quad = cQ.size();
				int kind = 6;
				double alpha = 2;
				double beta = 0;
				double a = 0.0;
				double b = 1.0;

				cgqf(n_quad, kind, alpha, beta, a, b, cQ.x.data(), cQ.w.data());
			}

			// 10: int_0^infty f(x) x^0 Exp[-x] dx
			void nw_Gauss_Laguerre_x0_int_0_pinfty(Q1<double, e_host> &cQ)
			{
				//    Input, int KIND, the rule.
				//    1, Legendre,             (a,b)       1.0
				//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
				//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
				//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
				//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
				//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
				//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
				//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

				int n_quad = cQ.size();
				int kind = 5;
				double alpha = 0;
				double beta = 0;
				double a = 0.0;
				double b = 1.0;

				cgqf(n_quad, kind, alpha, beta, a, b, cQ.x.data(), cQ.w.data());
			}

			// 11: int_0^infty f(x) x^1 Exp[-x] dx
			void nw_Gauss_Laguerre_x1_int_0_pinfty(Q1<double, e_host> &cQ)
			{
				//    Input, int KIND, the rule.
				//    1, Legendre,             (a,b)       1.0
				//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
				//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
				//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
				//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
				//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
				//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
				//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

				int n_quad = cQ.size();
				int kind = 5;
				double alpha = 1;
				double beta = 0;
				double a = 0.0;
				double b = 1.0;

				cgqf(n_quad, kind, alpha, beta, a, b, cQ.x.data(), cQ.w.data());
			}

			// 12: int_0^infty f(x) x^2 Exp[-x] dx
			void nw_Gauss_Laguerre_x2_int_0_pinfty(Q1<double, e_host> &cQ)
			{
				//    Input, int KIND, the rule.
				//    1, Legendre,             (a,b)       1.0
				//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
				//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
				//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
				//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
				//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
				//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
				//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

				int n_quad = cQ.size();
				int kind = 5;
				double alpha = 2;
				double beta = 0;
				double a = 0.0;
				double b = 1.0;

				cgqf(n_quad, kind, alpha, beta, a, b, cQ.x.data(), cQ.w.data());
			}

			// 13: int_0^infty f(x) Exp[-x]/Sqrt[x] dx
			void nw_Gauss_Laguerre_xi2_int_0_pinfty(Q1<double, e_host> &cQ)
			{
				switch(cQ.size())
				{
					case 2:
						cQ.x[0] = 2.7525512860841095e-1; 	cQ.w[0] = 1.6098281800110257e0;
						cQ.x[1] = 2.7247448713915890e0; 	cQ.w[1] = 1.6262567089449035e-1;
						break;
					case 4:
						cQ.x[0] = 1.4530352150331709e-1; 	cQ.w[0] = 1.3222940251164826e0;
						cQ.x[1] = 1.3390972881263614e0; 	cQ.w[1] = 4.1560465162978376e-1;
						cQ.x[2] = 3.9269635013582872e0; 	cQ.w[2] = 3.4155966014826951e-2;
						cQ.x[3] = 8.5886356890120343e0; 	cQ.w[3] = 3.9920814442273524e-4;
						break;
					case 8:
						cQ.x[0] = 7.4791882596818270e-2; 	cQ.w[0] = 1.0158589580332275e0;
						cQ.x[1] = 6.7724908764928915e-1; 	cQ.w[1] = 5.6129491705706735e-1;
						cQ.x[2] = 1.9051136350314284e0; 	cQ.w[2] = 1.6762008279797166e-1;
						cQ.x[3] = 3.8094763614849071e0; 	cQ.w[3] = 2.5760623071019947e-2;
						cQ.x[4] = 6.4831454286271704e0; 	cQ.w[4] = 1.8645680172483611e-3;
						cQ.x[5] = 1.0093323675221343e1; 	cQ.w[5] = 5.4237201850757630e-5;
						cQ.x[6] = 1.4972627088426393e1; 	cQ.w[6] = 4.6419616897304213e-7;
						cQ.x[7] = 2.1984272840962651e1; 	cQ.w[7] = 5.3096149480223645e-10;
						break;
					case 16:
						cQ.x[0] = 3.7962914575313455e-2; 	cQ.w[0] = 7.5047670518560479e-1;
						cQ.x[1] = 3.4220015601094768e-1; 	cQ.w[1] = 5.5491628460505980e-1;
						cQ.x[2] = 9.5355315539086550e-1; 	cQ.w[2] = 3.0253946815328497e-1;
						cQ.x[3] = 1.8779315076960743e0; 	cQ.w[3] = 1.2091626191182523e-1;
						cQ.x[4] = 3.1246010507021443e0; 	cQ.w[4] = 3.5106857663146861e-2;
						cQ.x[5] = 4.7067267076675872e0; 	cQ.w[5] = 7.3097806533088562e-3;
						cQ.x[6] = 6.6422151797414440e0; 	cQ.w[6] = 1.0725367310559441e-3;
						cQ.x[7] = 8.9550013377233902e0; 	cQ.w[7] = 1.0833168123639965e-4;
						cQ.x[8] = 1.1677033673975957e1; 	cQ.w[8] = 7.3011702591247521e-6;
						cQ.x[9] = 1.4851431341801250e1; 	cQ.w[9] = 3.1483355850911881e-7;
						cQ.x[10] = 1.8537743178606694e1; 	cQ.w[10] = 8.1976643295417932e-9;
						cQ.x[11] = 2.2821300693525208e1; 	cQ.w[11] = 1.1866582926793277e-10;
						cQ.x[12] = 2.7831438211328676e1; 	cQ.w[12] = 8.4300204226528951e-13;
						cQ.x[13] = 3.3781970488226166e1; 	cQ.w[13] = 2.3946880341856973e-15;
						cQ.x[14] = 4.1081666525491202e1; 	cQ.w[14] = 1.8463473073036584e-18;
						cQ.x[15] = 5.0777223877537080e1; 	cQ.w[15] = 1.4621352854768325e-22;
						break;
					case 24:
						cQ.x[0] = 2.5437996585689359e-2; 	cQ.w[0] = 6.2200206075592616e-1;
						cQ.x[1] = 2.2910231649262433e-1; 	cQ.w[1] = 5.0792308532951820e-1;
						cQ.x[2] = 6.3729027873266879e-1; 	cQ.w[2] = 3.3840894389128221e-1;
						cQ.x[3] = 1.2517406323627464e0; 	cQ.w[3] = 1.8364459415857036e-1;
						cQ.x[4] = 2.0751129098523806e0; 	cQ.w[4] = 8.0959353969207698e-2;
						cQ.x[5] = 3.1110524551477130e0; 	cQ.w[5] = 2.8889923149962199e-2;
						cQ.x[6] = 4.3642830769353062e0; 	cQ.w[6] = 8.3060098239551049e-3;
						cQ.x[7] = 5.8407332713236080e0; 	cQ.w[7] = 1.9127846396388306e-3;
						cQ.x[8] = 7.5477046800234544e0; 	cQ.w[8] = 3.5030086360234566e-4;
						cQ.x[9] = 9.4940953300264876e0; 	cQ.w[9] = 5.0571980554969778e-5;
						cQ.x[10] = 1.1690695926056073e1; 	cQ.w[10] = 5.6945173834696962e-6;
						cQ.x[11] = 1.4150586187285759e1; 	cQ.w[11] = 4.9373179873395010e-7;
						cQ.x[12] = 1.6889671928527108e1; 	cQ.w[12] = 3.2450282717915397e-8;
						cQ.x[13] = 1.9927425875242462e1; 	cQ.w[13] = 1.5860934990330765e-9;
						cQ.x[14] = 2.3287932824879917e1; 	cQ.w[14] = 5.6305930756763382e-11;
						cQ.x[15] = 2.7001406056472356e1; 	cQ.w[15] = 1.4093865163091778e-12;
						cQ.x[16] = 3.1106464709046565e1; 	cQ.w[16] = 2.3951797309583587e-14;
						cQ.x[17] = 3.5653703516328212e1; 	cQ.w[17] = 2.6303192453168170e-16;
						cQ.x[18] = 4.0711598185543107e1; 	cQ.w[18] = 1.7460319202373353e-18;
						cQ.x[19] = 4.6376979557540133e1; 	cQ.w[19] = 6.3767746470102769e-21;
						cQ.x[20] = 5.2795432527283630e1; 	cQ.w[20] = 1.1129154937804570e-23;
						cQ.x[21] = 6.0206666963057223e1; 	cQ.w[21] = 7.3700721603013398e-27;
						cQ.x[22] = 6.9068601975304369e1; 	cQ.w[22] = 1.1969225386627757e-30;
						cQ.x[23] = 8.0556280819950406e1; 	cQ.w[23] = 1.5871103023654473e-35;
						break;
					case 32:
						cQ.x[0] = 1.9127510968446856e-2; 	cQ.w[0] = 5.4275484988260796e-1;
						cQ.x[1] = 1.7221572414539558e-1; 	cQ.w[1] = 4.6598957212535609e-1;
						cQ.x[2] = 4.7875647727748885e-1; 	cQ.w[2] = 3.4337168469816740e-1;
						cQ.x[3] = 9.3948321450073428e-1; 	cQ.w[3] = 2.1699669861237368e-1;
						cQ.x[4] = 1.5555082314789380e0; 	cQ.w[4] = 1.1747996392819887e-1;
						cQ.x[5] = 2.3283376682103970e0; 	cQ.w[5] = 5.4406257907377837e-2;
						cQ.x[6] = 3.2598922564569419e0; 	cQ.w[6] = 2.1512081019758274e-2;
						cQ.x[7] = 4.3525345293301410e0; 	cQ.w[7] = 7.2451739570689175e-3;
						cQ.x[8] = 5.6091034574961513e0; 	cQ.w[8] = 2.0726581990151553e-3;
						cQ.x[9] = 7.0329577982838936e0; 	cQ.w[9] = 5.0196739702612497e-4;
						cQ.x[10] = 8.6280298574059291e0; 	cQ.w[10] = 1.0251858271572549e-4;
						cQ.x[11] = 1.0398891905552624e1; 	cQ.w[11] = 1.7576998461700718e-5;
						cQ.x[12] = 1.2350838217714770e1; 	cQ.w[12] = 2.5166805020623692e-6;
						cQ.x[13] = 1.4489986690780274e1; 	cQ.w[13] = 2.9910658734544941e-7;
						cQ.x[14] = 1.6823405362953694e1; 	cQ.w[14] = 2.9302506329522187e-8;
						cQ.x[15] = 1.9359271087268714e1; 	cQ.w[15] = 2.3472334846430987e-9;
						cQ.x[16] = 2.2107070382206007e1; 	cQ.w[16] = 1.5230434500290903e-10;
						cQ.x[17] = 2.5077856544198053e1; 	cQ.w[17] = 7.9183555338954479e-12;
						cQ.x[18] = 2.8284583194970531e1; 	cQ.w[18] = 3.2566814614194407e-13;
						cQ.x[19] = 3.1742543790616606e1; 	cQ.w[19] = 1.0437247453181695e-14;
						cQ.x[20] = 3.5469961396173283e1; 	cQ.w[20] = 2.5601867826448761e-16;
						cQ.x[21] = 3.9488797123368127e1; 	cQ.w[21] = 4.7037694213516382e-18;
						cQ.x[22] = 4.3825886369903902e1; 	cQ.w[22] = 6.3045091330075628e-20;
						cQ.x[23] = 4.8514583867416048e1; 	cQ.w[23] = 5.9657255685597023e-22;
						cQ.x[24] = 5.3597231826148512e1; 	cQ.w[24] = 3.8234137666012857e-24;
						cQ.x[25] = 5.9129027934391951e1; 	cQ.w[25] = 1.5723595577851821e-26;
						cQ.x[26] = 6.5184426376135782e1; 	cQ.w[26] = 3.8582071909299337e-29;
						cQ.x[27] = 7.1868499359551422e1; 	cQ.w[27] = 5.0993217982259985e-32;
						cQ.x[28] = 7.9339086528823201e1; 	cQ.w[28] = 3.1147812492595276e-35;
						cQ.x[29] = 8.7856119943133525e1; 	cQ.w[29] = 6.8422760225114810e-39;
						cQ.x[30] = 9.7916716426062762e1; 	cQ.w[30] = 3.3594959802163184e-43;
						cQ.x[31] = 1.1079926894707576e2; 	cQ.w[31] = 1.1088644990767160e-48;
						break;
					case 40:
						cQ.x[0] = 1.5325663331507189e-2; 	cQ.w[0] = 4.8767170761442366e-1;
						cQ.x[1] = 1.3796600174102882e-1; 	cQ.w[1] = 4.3154989254866095e-1;
						cQ.x[2] = 3.8343384139288653e-1; 	cQ.w[2] = 3.3787593855180386e-1;
						cQ.x[3] = 7.5210508353159711e-1; 	cQ.w[3] = 2.3396146760859843e-1;
						cQ.x[4] = 1.2445475511131953e0; 	cQ.w[4] = 1.4320149538368678e-1;
						cQ.x[5] = 1.8615258453171583e0; 	cQ.w[5] = 7.7417199829675957e-2;
						cQ.x[6] = 2.6040079765974233e0; 	cQ.w[6] = 3.6931502308855776e-2;
						cQ.x[7] = 3.4731739113527927e0; 	cQ.w[7] = 1.5528078810799283e-2;
						cQ.x[8] = 4.4704262206423379e0; 	cQ.w[8] = 5.7463950661794487e-3;
						cQ.x[9] = 5.5974030703416056e0; 	cQ.w[9] = 1.8686353739223537e-3;
						cQ.x[10] = 6.8559938552704945e0; 	cQ.w[10] = 5.3295577338670161e-4;
						cQ.x[11] = 8.2483578563697614e0; 	cQ.w[11] = 1.3303482483630581e-4;
						cQ.x[12] = 9.7769463941849785e0; 	cQ.w[12] = 2.8992888382259350e-5;
						cQ.x[13] = 1.1444529069317641e1; 	cQ.w[13] = 5.5014492548024867e-6;
						cQ.x[14] = 1.3254224828602629e1; 	cQ.w[14] = 9.0610711317226041e-7;
						cQ.x[15] = 1.5209538784718065e1; 	cQ.w[15] = 1.2908925900239279e-7;
						cQ.x[16] = 1.7314405960671086e1; 	cQ.w[16] = 1.5845846578300512e-8;
						cQ.x[17] = 1.9573243448531832e1; 	cQ.w[17] = 1.6686111852173423e-9;
						cQ.x[18] = 2.1991012891259404e1; 	cQ.w[18] = 1.4999439829316549e-10;
						cQ.x[19] = 2.4573295756577260e1; 	cQ.w[19] = 1.1446559437004013e-11;
						cQ.x[20] = 2.7326384629341011e1; 	cQ.w[20] = 7.3696996403182414e-13;
						cQ.x[21] = 3.0257394787336121e1; 	cQ.w[21] = 3.9750416693969613e-14;
						cQ.x[22] = 3.3374401770413884e1; 	cQ.w[22] = 1.7818922081332433e-15;
						cQ.x[23] = 3.6686612696144794e1; 	cQ.w[23] = 6.5783240287291484e-17;
						cQ.x[24] = 4.0204582016183319e1; 	cQ.w[24] = 1.9793057979240198e-18;
						cQ.x[25] = 4.3940486724615585e1; 	cQ.w[25] = 4.7956805512357826e-20;
						cQ.x[26] = 4.7908482506508098e1; 	cQ.w[26] = 9.2269421491078646e-22;
						cQ.x[27] = 5.2125172273284642e1; 	cQ.w[27] = 1.3868239639057759e-23;
						cQ.x[28] = 5.6610234272577454e1; 	cQ.w[28] = 1.5970374396346345e-25;
						cQ.x[29] = 6.1387282639433491e1; 	cQ.w[29] = 1.3766914087267507e-27;
						cQ.x[30] = 6.6485076698107068e1; 	cQ.w[30] = 8.6356346860439419e-30;
						cQ.x[31] = 7.1939271994415582e1; 	cQ.w[31] = 3.8059812200638324e-32;
						cQ.x[32] = 7.7795048292119587e1; 	cQ.w[32] = 1.1274385558954499e-34;
						cQ.x[33] = 8.4111230042098140e1; 	cQ.w[33] = 2.1190293856555058e-37;
						cQ.x[34] = 9.0967109281657703e1; 	cQ.w[34] = 2.3384153996394935e-40;
						cQ.x[35] = 9.8474564473065897e1; 	cQ.w[35] = 1.3584866032309069e-43;
						cQ.x[36] = 1.0680170504629888e2; 	cQ.w[36] = 3.5284114086777935e-47;
						cQ.x[37] = 1.1622557258080335e2; 	cQ.w[37] = 3.1343129875498926e-51;
						cQ.x[38] = 1.2727662323235295e2; 	cQ.w[38] = 5.7648218162351852e-52;
						cQ.x[39] = 1.4132130003237776e2; 	cQ.w[39] = 5.6256581525587782e-50;
						break;
					case 48:
						cQ.x[0] = 1.2784572337116694e-2; 	cQ.w[0] = 4.4654004699506952e-1;
						cQ.x[1] = 1.1508148358837804e-1; 	cQ.w[1] = 4.0322602695846852e-1;
						cQ.x[2] = 3.1978387827943194e-1; 	cQ.w[2] = 3.2875926119699403e-1;
						cQ.x[3] = 6.2710953585675593e-1; 	cQ.w[3] = 2.4196623372503065e-1;
						cQ.x[4] = 1.0373867221540579e0; 	cQ.w[4] = 1.6070790016128415e-1;
						cQ.x[5] = 1.5510561307962492e0; 	cQ.w[5] = 9.6279813462180107e-2;
						cQ.x[6] = 2.1686735147948251e0; 	cQ.w[6] = 5.2000680542488898e-2;
						cQ.x[7] = 2.8909130466944233e0; 	cQ.w[7] = 2.5302750225600317e-2;
						cQ.x[8] = 3.7185714571096889e0; 	cQ.w[8] = 1.1083262578877953e-2;
						cQ.x[9] = 4.6525730143474026e0; 	cQ.w[9] = 4.3662639521934365e-3;
						cQ.x[10] = 5.6939754224427003e0; 	cQ.w[10] = 1.5453916433633038e-3;
						cQ.x[11] = 6.8439767318335151e0; 	cQ.w[11] = 4.9083618546105323e-4;
						cQ.x[12] = 8.1039233766453468e0; 	cQ.w[12] = 1.3970864889083082e-4;
						cQ.x[13] = 9.4753194758919079e0; 	cQ.w[13] = 3.5583611836440579e-5;
						cQ.x[14] = 1.0959837563734646e1; 	cQ.w[14] = 8.0964663575545199e-6;
						cQ.x[15] = 1.2559330947448595e1; 	cQ.w[15] = 1.6427089095665616e-6;
						cQ.x[16] = 1.4275847932399555e1; 	cQ.w[16] = 2.9659388100100142e-7;
						cQ.x[17] = 1.6111648203063945e1; 	cQ.w[17] = 4.7547324295148571e-8;
						cQ.x[18] = 1.8069221710408724e1; 	cQ.w[18] = 6.7511694622554086e-9;
						cQ.x[19] = 2.0151310492061653e1; 	cQ.w[19] = 8.4671875670220449e-10;
						cQ.x[20] = 2.2360933946965237e1; 	cQ.w[20] = 9.3520207719701412e-11;
						cQ.x[21] = 2.4701418206395586e1; 	cQ.w[21] = 9.0666067491677764e-12;
						cQ.x[22] = 2.7176430396128092e1; 	cQ.w[22] = 7.6873634497532554e-13;
						cQ.x[23] = 2.9790018780757460e1; 	cQ.w[23] = 5.6775347903386773e-14;
						cQ.x[24] = 3.2546660035349492e1; 	cQ.w[24] = 3.6363382609864670e-15;
						cQ.x[25] = 3.5451315222092618e1; 	cQ.w[25] = 2.0098102983149713e-16;
						cQ.x[26] = 3.8509496489187991e1; 	cQ.w[26] = 9.5336623423732135e-18;
						cQ.x[27] = 4.1727347097011937e1; 	cQ.w[27] = 3.8577586536742773e-19;
						cQ.x[28] = 4.5111738172331784e1; 	cQ.w[28] = 1.3225902352962134e-20;
						cQ.x[29] = 4.8670386683149351e1; 	cQ.w[29] = 3.8125203773751006e-22;
						cQ.x[30] = 5.2412000646782474e1; 	cQ.w[30] = 9.1612021970998299e-24;
						cQ.x[31] = 5.6346459734243123e1; 	cQ.w[31] = 1.8171931294097668e-25;
						cQ.x[32] = 6.0485042530508732e1; 	cQ.w[32] = 2.9424831099483641e-27;
						cQ.x[33] = 6.4840716257286064e1; 	cQ.w[33] = 3.8399414135897572e-29;
						cQ.x[34] = 6.9428511589042573e1; 	cQ.w[34] = 3.9790937885390546e-31;
						cQ.x[35] = 7.4266015688703556e1; 	cQ.w[35] = 3.2177464158769137e-33;
						cQ.x[36] = 7.9374033184791669e1; 	cQ.w[36] = 1.9893600118442518e-35;
						cQ.x[37] = 8.4777491893280191e1; 	cQ.w[37] = 9.1748033317163899e-38;
						cQ.x[38] = 9.0506715916613337e1; 	cQ.w[38] = 3.0635979294223916e-40;
						cQ.x[39] = 9.6599269661943114e1; 	cQ.w[39] = 7.1378771235620131e-43;
						cQ.x[40] = 1.0310272648925973e2; 	cQ.w[40] = 1.1074129639565051e-45;
						cQ.x[41] = 1.1007901167333772e2; 	cQ.w[41] = 1.0766445013940194e-48;
						cQ.x[42] = 1.1761159733441177e2; 	cQ.w[42] = 6.0445878630042245e-52;
						cQ.x[43] = 1.2581828912230779e2; 	cQ.w[43] = 1.7467465671838054e-55;
						cQ.x[44] = 1.3487618875613396e2; 	cQ.w[44] = 2.1867649475035571e-59;
						cQ.x[45] = 1.4507736943128927e2; 	cQ.w[45] = 8.9374048433626366e-64;
						cQ.x[46] = 1.5698162310364024e2; 	cQ.w[46] = 6.9616827754363648e-69;
						cQ.x[47] = 1.7203286674516623e2; 	cQ.w[47] = 7.3843030686139940e-50;
						break;
					case 56:
						cQ.x[0] = 1.0966296964922414e-2; 	cQ.w[0] = 4.1431861461592488e-1;
						cQ.x[1] = 9.8709503981500921e-2; 	cQ.w[1] = 3.7958824335151592e-1;
						cQ.x[2] = 2.7426441281900728e-1; 	cQ.w[2] = 3.1859560878006560e-1;
						cQ.x[3] = 5.3776830760580984e-1; 	cQ.w[3] = 2.4493804875366549e-1;
						cQ.x[4] = 8.8942785389626524e-1; 	cQ.w[4] = 1.7245385388352993e-1;
						cQ.x[5] = 1.3295199945387192e0; 	cQ.w[5] = 1.1116552769641440e-1;
						cQ.x[6] = 1.8583931587330957e0; 	cQ.w[6] = 6.5583934984821083e-2;
						cQ.x[7] = 2.4764687971240665e0; 	cQ.w[7] = 3.5397526052325080e-2;
						cQ.x[8] = 3.1842432594503771e0; 	cQ.w[8] = 1.7469565741779008e-2;
						cQ.x[9] = 3.9822900352751187e0; 	cQ.w[9] = 7.8791146628173459e-3;
						cQ.x[10] = 4.8712623827436266e0; 	cQ.w[10] = 3.2454637275858797e-3;
						cQ.x[11] = 5.8518963752570206e0; 	cQ.w[11] = 1.2200069198095973e-3;
						cQ.x[12] = 6.9250144015294926e0; 	cQ.w[12] = 4.1819509658049854e-4;
						cQ.x[13] = 8.0915291608577938e0; 	cQ.w[13] = 1.3059682676732009e-4;
						cQ.x[14] = 9.3524482027415579e0; 	cQ.w[14] = 3.7118324819779832e-5;
						cQ.x[15] = 1.0708879068458080e1; 	cQ.w[15] = 9.5910914884881867e-6;
						cQ.x[16] = 1.2162035102064631e1; 	cQ.w[16] = 2.2503384827035532e-6;
						cQ.x[17] = 1.3713242009881661e1; 	cQ.w[17] = 4.7880157280422767e-7;
						cQ.x[18] = 1.5363945261179596e1; 	cQ.w[18] = 9.2250198493335919e-8;
						cQ.x[19] = 1.7115718439020667e1; 	cQ.w[19] = 1.6069558346408703e-8;
						cQ.x[20] = 1.8970272669583436e1; 	cQ.w[20] = 2.5265449989883177e-9;
						cQ.x[21] = 2.0929467281561854e1; 	cQ.w[21] = 3.5787741726578897e-10;
						cQ.x[22] = 2.2995321875320932e1; 	cQ.w[22] = 4.5577902726623018e-11;
						cQ.x[23] = 2.5170030015603701e1; 	cQ.w[23] = 5.2076385947795911e-12;
						cQ.x[24] = 2.7455974803255672e1; 	cQ.w[24] = 5.3255665746212511e-13;
						cQ.x[25] = 2.9855746632650704e1; 	cQ.w[25] = 4.8619672412277693e-14;
						cQ.x[26] = 3.2372163504856039e1; 	cQ.w[26] = 3.9515202575705826e-15;
						cQ.x[27] = 3.5008294345466858e1; 	cQ.w[27] = 2.8503559634838643e-16;
						cQ.x[28] = 3.7767485874979349e1; 	cQ.w[28] = 1.8187535954103432e-17;
						cQ.x[29] = 4.0653393704581200e1; 	cQ.w[29] = 1.0228519632208071e-18;
						cQ.x[30] = 4.3670018489449905e1; 	cQ.w[30] = 5.0500039302829468e-20;
						cQ.x[31] = 4.6821748176145590e1; 	cQ.w[31] = 2.1793149108595921e-21;
						cQ.x[32] = 5.0113407645739421e1; 	cQ.w[32] = 8.1812395526908028e-23;
						cQ.x[33] = 5.3550317401224539e1; 	cQ.w[33] = 2.6576488518963704e-24;
						cQ.x[34] = 5.7138363406574332e1; 	cQ.w[34] = 7.4271309204470600e-26;
						cQ.x[35] = 6.0884080798556528e1; 	cQ.w[35] = 1.7740931536888198e-27;
						cQ.x[36] = 6.4794755023535346e1; 	cQ.w[36] = 3.5960724740047246e-29;
						cQ.x[37] = 6.8878545092128855e1; 	cQ.w[37] = 6.1357429171822799e-31;
						cQ.x[38] = 7.3144635233014682e1; 	cQ.w[38] = 8.7325878464563983e-33;
						cQ.x[39] = 7.7603423474903333e1; 	cQ.w[39] = 1.0260913087634106e-34;
						cQ.x[40] = 8.2266758923044609e1; 	cQ.w[40] = 9.8379556343138002e-37;
						cQ.x[41] = 8.7148244251462459e1; 	cQ.w[41] = 7.5937692877807736e-39;
						cQ.x[42] = 9.2263627069792685e1; 	cQ.w[42] = 4.6460504504603365e-41;
						cQ.x[43] = 9.7631314803869374e1; 	cQ.w[43] = 2.2125228051355289e-43;
						cQ.x[44] = 1.0327306509488168e2; 	cQ.w[44] = 8.0267804054671947e-46;
						cQ.x[45] = 1.0921493206709584e2; 	cQ.w[45] = 2.1621413990444480e-48;
						cQ.x[46] = 1.1548859679336965e2; 	cQ.w[46] = 4.1913798221383726e-51;
						cQ.x[47] = 1.2213329501416645e2; 	cQ.w[47] = 5.6258187814957665e-54;
						cQ.x[48] = 1.2919871245898944e2; 	cQ.w[48] = 4.9791444365660730e-57;
						cQ.x[49] = 1.3674952821727593e2; 	cQ.w[49] = 2.7270187769624205e-60;
						cQ.x[50] = 1.4487294472877377e2; 	cQ.w[50] = 8.4854893562244169e-64;
						cQ.x[51] = 1.5369207574377554e2; 	cQ.w[51] = 1.3300105642340421e-67;
						cQ.x[52] = 1.6339209491703519e2; 	cQ.w[52] = 8.7671924329217860e-72;
						cQ.x[53] = 1.7427858611626866e2; 	cQ.w[53] = 1.8070041919711762e-76;
						cQ.x[54] = 1.8693771869007120e2; 	cQ.w[54] = 1.0229121536329777e-48;
						cQ.x[55] = 2.0288303763687224e2; 	cQ.w[55] = 7.4032472491198754e-53;
						break;
					case 64:
						cQ.x[0] = 9.6008293650696286e-3; 	cQ.w[0] = 3.8819522372817551e-1;
						cQ.x[1] = 8.6416074005085619e-2; 	cQ.w[1] = 3.5954616781559853e-1;
						cQ.x[2] = 2.4009251326551299e-1; 	cQ.w[2] = 3.0842087059670877e-1;
						cQ.x[3] = 4.7072219845762405e-1; 	cQ.w[3] = 2.4500654632827139e-1;
						cQ.x[4] = 7.7844358612953627e-1; 	cQ.w[4] = 1.8021735675289784e-1;
						cQ.x[5] = 1.1634419969552807e0; 	cQ.w[5] = 1.2272144200898013e-1;
						cQ.x[6] = 1.6259502332928078e0; 	cQ.w[6] = 7.7347909621273805e-2;
						cQ.x[7] = 2.1662493604094500e0; 	cQ.w[7] = 4.5108620335648845e-2;
						cQ.x[8] = 2.7846696577606368e0; 	cQ.w[8] = 2.4333837728938679e-2;
						cQ.x[9] = 3.4815917481914229e0; 	cQ.w[9] = 1.2137724813851775e-2;
						cQ.x[10] = 4.2574479145349120e0; 	cQ.w[10] = 5.5956788032115785e-3;
						cQ.x[11] = 5.1127236148341415e0; 	cQ.w[11] = 2.3831276289143345e-3;
						cQ.x[12] = 6.0479592093454125e0; 	cQ.w[12] = 9.3710307561682273e-4;
						cQ.x[13] = 7.0637519146269628e0; 	cQ.w[13] = 3.4002817652561882e-4;
						cQ.x[14] = 8.1607580024185598e0; 	cQ.w[14] = 1.1377487520080482e-4;
						cQ.x[15] = 9.3396952637232506e0; 	cQ.w[15] = 3.5080971696187810e-5;
						cQ.x[16] = 1.0601345761568931e1; 	cQ.w[16] = 9.9598490651819740e-6;
						cQ.x[17] = 1.1946558899421827e1; 	cQ.w[17] = 2.6014940064763985e-6;
						cQ.x[18] = 1.3376254836226536e1; 	cQ.w[18] = 6.2457459723578062e-7;
						cQ.x[19] = 1.4891428283653888e1; 	cQ.w[19] = 1.3769162244308700e-7;
						cQ.x[20] = 1.6493152726464029e1; 	cQ.w[20] = 2.7843814305870358e-8;
						cQ.x[21] = 1.8182585113077513e1; 	cQ.w[21] = 5.1587944588527896e-9;
						cQ.x[22] = 1.9960971070661513e1; 	cQ.w[22] = 8.7463733196968069e-10;
						cQ.x[23] = 2.1829650707488929e1; 	cQ.w[23] = 1.3551560975549108e-10;
						cQ.x[24] = 2.3790065075269629e1; 	cQ.w[24] = 1.9160633017471715e-11;
						cQ.x[25] = 2.5843763375899103e1; 	cQ.w[25] = 2.4684289732011134e-12;
						cQ.x[26] = 2.7992411011009711e1; 	cQ.w[26] = 2.8926946423808331e-13;
						cQ.x[27] = 3.0237798589328330e1; 	cQ.w[27] = 3.0780994607070916e-14;
						cQ.x[28] = 3.2581852026749626e1; 	cQ.w[28] = 2.9684476750277130e-15;
						cQ.x[29] = 3.5026643897992604e1; 	cQ.w[29] = 2.5890963186787431e-16;
						cQ.x[30] = 3.7574406227691570e1; 	cQ.w[30] = 2.0378664608465850e-17;
						cQ.x[31] = 4.0227544944021652e1; 	cQ.w[31] = 1.4440214633856584e-18;
						cQ.x[32] = 4.2988656261067658e1; 	cQ.w[32] = 9.1880153546594432e-20;
						cQ.x[33] = 4.5860545309175844e1; 	cQ.w[33] = 5.2349151678696223e-21;
						cQ.x[34] = 4.8846247398170421e1; 	cQ.w[34] = 2.6627357180671792e-22;
						cQ.x[35] = 5.1949052380103696e1; 	cQ.w[35] = 1.2051968064012909e-23;
						cQ.x[36] = 5.5172532680822943e1; 	cQ.w[36] = 4.8368069192953299e-25;
						cQ.x[37] = 5.8520575699338077e1; 	cQ.w[37] = 1.7145660967538707e-26;
						cQ.x[38] = 6.1997421439210377e1; 	cQ.w[38] = 5.3458472401161465e-28;
						cQ.x[39] = 6.5607706448472386e1; 	cQ.w[39] = 1.4593090013536809e-29;
						cQ.x[40] = 6.9356515419807817e1; 	cQ.w[40] = 3.4702060405641224e-31;
						cQ.x[41] = 7.3249442163003752e1; 	cQ.w[41] = 7.1487577917588421e-33;
						cQ.x[42] = 7.7292662138289202e1; 	cQ.w[42] = 1.2679827052732978e-34;
						cQ.x[43] = 8.1493019376821402e1; 	cQ.w[43] = 1.9233416135935014e-36;
						cQ.x[44] = 8.5858131478224402e1; 	cQ.w[44] = 2.4761711159527361e-38;
						cQ.x[45] = 9.0396517560549849e1; 	cQ.w[45] = 2.6829949635287387e-40;
						cQ.x[46] = 9.5117755689162181e1; 	cQ.w[46] = 2.4235590682611838e-42;
						cQ.x[47] = 1.0003267864790140e2; 	cQ.w[47] = 1.8056080277573288e-44;
						cQ.x[48] = 1.0515362028215433e2; 	cQ.w[48] = 1.0960434057957593e-46;
						cQ.x[49] = 1.1049472958850639e2; 	cQ.w[49] = 5.3454875034721357e-49;
						cQ.x[50] = 1.1607237715014811e2; 	cQ.w[50] = 2.0609725041113895e-51;
						cQ.x[51] = 1.2190568994074218e2; 	cQ.w[51] = 6.1641547666785974e-54;
						cQ.x[52] = 1.2801726858944196e2; 	cQ.w[52] = 1.3986145848103912e-56;
						cQ.x[53] = 1.3443417070003274e2; 	cQ.w[53] = 2.3439570024259610e-59;
						cQ.x[54] = 1.4118929376119792e2; 	cQ.w[54] = 2.8089345154809745e-62;
						cQ.x[55] = 1.4832337939847587e2; 	cQ.w[55] = 2.3123103281927504e-65;
						cQ.x[56] = 1.5588802451891946e2; 	cQ.w[56] = 1.2428488323660627e-68;
						cQ.x[57] = 1.6395040789497031e2; 	cQ.w[57] = 4.0831715944879700e-72;
						cQ.x[58] = 1.7260112638093499e2; 	cQ.w[58] = 7.5024317376094500e-76;
						cQ.x[59] = 1.8196813220934213e2; 	cQ.w[59] = 6.8024601738732743e-80;
						cQ.x[60] = 1.9224396472236875e2; 	cQ.w[60] = 2.5224989666770768e-84;
						cQ.x[61] = 2.0374654228943264e2; 	cQ.w[61] = 9.5158011667265254e-49;
						cQ.x[62] = 2.1708611403654403e2; 	cQ.w[62] = 5.3192427849433605e-52;
						cQ.x[63] = 2.3383975178282573e2; 	cQ.w[63] = 2.9535559171937292e-66;
						break;
					case 72:
						cQ.x[0] = 8.5377530337449248e-3; 	cQ.w[0] = 3.6646133805754739e-1;
						cQ.x[1] = 7.6845831749819249e-2; 	cQ.w[1] = 3.4230533863506727e-1;
						cQ.x[2] = 2.1349429706226695e-1; 	cQ.w[2] = 2.9865564028363514e-1;
						cQ.x[3] = 4.1854784885205123e-1; 	cQ.w[3] = 2.4337240372556322e-1;
						cQ.x[4] = 6.9210374777076850e-1; 	cQ.w[4] = 1.8521384480327197e-1;
						cQ.x[5] = 1.0342920697248060e0; 	cQ.w[5] = 1.3161967742383411e-1;
						cQ.x[6] = 1.4452760476859431e0; 	cQ.w[6] = 8.7325860377519644e-2;
						cQ.x[7] = 1.9252525030070945e0; 	cQ.w[7] = 5.4082218044787982e-2;
						cQ.x[8] = 2.4744523690138993e0; 	cQ.w[8] = 3.1257614600552892e-2;
						cQ.x[9] = 3.0931413102677093e0; 	cQ.w[9] = 1.6855126430119594e-2;
						cQ.x[10] = 3.7816204415611087e0; 	cQ.w[10] = 8.4772081114631611e-3;
						cQ.x[11] = 4.5402271514220959e0; 	cQ.w[11] = 3.9753232176891744e-3;
						cQ.x[12] = 5.3693360356771178e0; 	cQ.w[12] = 1.7375161592345435e-3;
						cQ.x[13] = 6.2693599474670913e0; 	cQ.w[13] = 7.0752789404236957e-4;
						cQ.x[14] = 7.2407511710366363e0; 	cQ.w[14] = 2.6830022946769146e-4;
						cQ.x[15] = 8.2840027276389276e0; 	cQ.w[15] = 9.4699597895333223e-5;
						cQ.x[16] = 9.3996498230328796e0; 	cQ.w[16] = 3.1095207369238697e-5;
						cQ.x[17] = 1.0588271447314284e1; 	cQ.w[17] = 9.4930705560516880e-6;
						cQ.x[18] = 1.1850492139239455e1; 	cQ.w[18] = 2.6928856004713719e-6;
						cQ.x[19] = 1.3186983928793866e1; 	cQ.w[19] = 7.0931074479099719e-7;
						cQ.x[20] = 1.4598468473558433e1; 	cQ.w[20] = 1.7336070598064455e-7;
						cQ.x[21] = 1.6085719406466767e1; 	cQ.w[21] = 3.9284921495770124e-8;
						cQ.x[22] = 1.7649564914868501e1; 	cQ.w[22] = 8.2471566241554238e-9;
						cQ.x[23] = 1.9290890573464524e1; 	cQ.w[23] = 1.6025182143590107e-9;
						cQ.x[24] = 2.1010642456716749e1; 	cQ.w[24] = 2.8794742979351125e-10;
						cQ.x[25] = 2.2809830559825846e1; 	cQ.w[25] = 4.7796759804694583e-11;
						cQ.x[26] = 2.4689532561396773e1; 	cQ.w[26] = 7.3213811680310381e-12;
						cQ.x[27] = 2.6650897965571812e1; 	cQ.w[27] = 1.0337134692368372e-12;
						cQ.x[28] = 2.8695152666822676e1; 	cQ.w[28] = 1.3436626222686719e-13;
						cQ.x[29] = 3.0823603986900308e1; 	cQ.w[29] = 1.6058255460718379e-14;
						cQ.x[30] = 3.3037646240818195e1; 	cQ.w[30] = 1.7620660009112015e-15;
						cQ.x[31] = 3.5338766897405687e1; 	cQ.w[31] = 1.7726340324641576e-16;
						cQ.x[32] = 3.7728553410174551e1; 	cQ.w[32] = 1.6323127659605682e-17;
						cQ.x[33] = 4.0208700806318587e1; 	cQ.w[33] = 1.3735451099194697e-18;
						cQ.x[34] = 4.2781020136014409e1; 	cQ.w[34] = 1.0542790113729066e-19;
						cQ.x[35] = 4.5447447901312514e1; 	cQ.w[35] = 7.3672431203789817e-21;
						cQ.x[36] = 4.8210056604429630e1; 	cQ.w[36] = 4.6773118965550795e-22;
						cQ.x[37] = 5.1071066579968236e1; 	cQ.w[37] = 2.6919741928296506e-23;
						cQ.x[38] = 5.4032859305501852e1; 	cQ.w[38] = 1.4012014167408308e-24;
						cQ.x[39] = 5.7097992421357986e1; 	cQ.w[39] = 6.5793246486772861e-26;
						cQ.x[40] = 6.0269216734953006e1; 	cQ.w[40] = 2.7792387816429709e-27;
						cQ.x[41] = 6.3549495539818125e1; 	cQ.w[41] = 1.0530654909535471e-28;
						cQ.x[42] = 6.6942026647284918e1; 	cQ.w[42] = 3.5677140725415228e-30;
						cQ.x[43] = 7.0450267613328067e1; 	cQ.w[43] = 1.0770556703499718e-31;
						cQ.x[44] = 7.4077964749136873e1; 	cQ.w[44] = 2.8865806752335222e-33;
						cQ.x[45] = 7.7829186638083087e1; 	cQ.w[45] = 6.8402497184153889e-35;
						cQ.x[46] = 8.1708363052613987e1; 	cQ.w[46] = 1.4268974543918404e-36;
						cQ.x[47] = 8.5720330384148150e1; 	cQ.w[47] = 2.6077204118485067e-38;
						cQ.x[48] = 8.9870384983719634e1; 	cQ.w[48] = 4.1533050523669125e-40;
						cQ.x[49] = 9.4164346183819770e1; 	cQ.w[49] = 5.7316992987700021e-42;
						cQ.x[50] = 9.8608631264982208e1; 	cQ.w[50] = 6.8102796020749233e-44;
						cQ.x[51] = 1.0321034529045834e2; 	cQ.w[51] = 6.9179617811502756e-46;
						cQ.x[52] = 1.0797738962609581e2; 	cQ.w[52] = 5.9610219491215519e-48;
						cQ.x[53] = 1.1291859418950553e2; 	cQ.w[53] = 4.3190839650906910e-50;
						cQ.x[54] = 1.1804388018177253e2; 	cQ.w[54] = 2.6056833265450411e-52;
						cQ.x[55] = 1.2336446247428126e2; 	cQ.w[55] = 1.2944535536043246e-54;
						cQ.x[56] = 1.2889310430876951e2; 	cQ.w[56] = 5.2287531735186084e-57;
						cQ.x[57] = 1.3464444208967572e2; 	cQ.w[57] = 1.6926455724514031e-59;
						cQ.x[58] = 1.4063540573754781e2; 	cQ.w[58] = 4.3183461553903832e-62;
						cQ.x[59] = 1.4688577190557120e2; 	cQ.w[59] = 8.5145486452100040e-65;
						cQ.x[60] = 1.5341890608220506e2; 	cQ.w[60] = 1.2678660729891205e-67;
						cQ.x[61] = 1.6026278017056644e2; 	cQ.w[61] = 1.3869468532834784e-70;
						cQ.x[62] = 1.6745140389442310e2; 	cQ.w[62] = 1.0778332837412283e-73;
						cQ.x[63] = 1.7502689981511292e2; 	cQ.w[63] = 5.7084939874269743e-77;
						cQ.x[64] = 1.8304262155336948e2; 	cQ.w[64] = 1.9550707098857624e-80;
						cQ.x[65] = 1.9156804971323859e2; 	cQ.w[65] = 4.0440344494888973e-84;
						cQ.x[66] = 2.0069691105825289e2; 	cQ.w[66] = 4.6082741548145242e-88;
						cQ.x[67] = 2.1056162324265798e2; 	cQ.w[67] = 2.5411102164420924e-92;
						cQ.x[68] = 2.2136152670104394e2; 	cQ.w[68] = 5.5818464378910382e-97;
						cQ.x[69] = 2.3342593022383359e2; 	cQ.w[69] = 3.8018252344681360e-50;
						cQ.x[70] = 2.4738731475402643e2; 	cQ.w[70] = 2.0795745925390798e-53;
						cQ.x[71] = 2.6488137073545847e2; 	cQ.w[71] = 5.9029492520141047e-68;
						break;
					case 80:
						cQ.x[0] = 7.6866318444038650e-3; 	cQ.w[0] = 3.4801121234964061e-1;
						cQ.x[1] = 6.9184104726922287e-2; 	cQ.w[1] = 3.2728554809172511e-1;
						cQ.x[2] = 1.9220262418784241e-1; 	cQ.w[2] = 2.8945683894699143e-1;
						cQ.x[3] = 3.7678938735594810e-1; 	cQ.w[3] = 2.4073751764332210e-1;
						cQ.x[4] = 6.2301531452046157e-1; 	cQ.w[4] = 1.8826779554366350e-1;
						cQ.x[5] = 9.3097519935847685e-1; 	cQ.w[5] = 1.3843298565975990e-1;
						cQ.x[6] = 1.3007879104323092e0; 	cQ.w[6] = 9.5693623713172557e-2;
						cQ.x[7] = 1.7325966449946741e0; 	cQ.w[7] = 6.2179086370712325e-2;
						cQ.x[8] = 2.2265692364173354e0; 	cQ.w[8] = 3.7970908655009312e-2;
						cQ.x[9] = 2.7828985168491514e0; 	cQ.w[9] = 2.1788109306695115e-2;
						cQ.x[10] = 3.4018027370152287e0; 	cQ.w[10] = 1.1745069517438239e-2;
						cQ.x[11] = 4.0835260453933245e0; 	cQ.w[11] = 5.9463917421429804e-3;
						cQ.x[12] = 4.8283390293501774e0; 	cQ.w[12] = 2.8268077494244823e-3;
						cQ.x[13] = 5.6365393211929032e0; 	cQ.w[13] = 1.2614066531994309e-3;
						cQ.x[14] = 6.5084522724931878e0; 	cQ.w[14] = 5.2818895724163281e-4;
						cQ.x[15] = 7.4444317004794741e0; 	cQ.w[15] = 2.0746547336341577e-4;
						cQ.x[16] = 8.4448607107699771e0; 	cQ.w[16] = 7.6411540013651783e-5;
						cQ.x[17] = 9.5101526012431861e0; 	cQ.w[17] = 2.6378489263127714e-5;
						cQ.x[18] = 1.0640751852419332e1; 	cQ.w[18] = 8.5315214346510367e-6;
						cQ.x[19] = 1.1837135210363891e1; 	cQ.w[19] = 2.5839402017229379e-6;
						cQ.x[20] = 1.3099812868831478e1; 	cQ.w[20] = 7.3248237980755977e-7;
						cQ.x[21] = 1.4429329758155656e1; 	cQ.w[21] = 1.9423835172306337e-7;
						cQ.x[22] = 1.5826266949269123e1; 	cQ.w[22] = 4.8155353038663436e-8;
						cQ.x[23] = 1.7291243182222957e1; 	cQ.w[23] = 1.1154712143820550e-8;
						cQ.x[24] = 1.8824916529679155e1; 	cQ.w[24] = 2.4126359656654055e-9;
						cQ.x[25] = 2.0427986207095807e1; 	cQ.w[25] = 4.8690306206928333e-10;
						cQ.x[26] = 2.2101194542730669e1; 	cQ.w[26] = 9.1619835022640038e-11;
						cQ.x[27] = 2.3845329122181760e1; 	cQ.w[27] = 1.6061727543126440e-11;
						cQ.x[28] = 2.5661225123992617e1; 	cQ.w[28] = 2.6211365629406585e-12;
						cQ.x[29] = 2.7549767864910021e1; 	cQ.w[29] = 3.9783127491115241e-13;
						cQ.x[30] = 2.9511895575734651e1; 	cQ.w[30] = 5.6106664386726106e-14;
						cQ.x[31] = 3.1548602431399443e1; 	cQ.w[31] = 7.3452477074096906e-15;
						cQ.x[32] = 3.3660941862004724e1; 	cQ.w[32] = 8.9170156740924202e-16;
						cQ.x[33] = 3.5850030175103469e1; 	cQ.w[33] = 1.0027003519531224e-16;
						cQ.x[34] = 3.8117050523647834e1; 	cQ.w[34] = 1.0431576805686864e-17;
						cQ.x[35] = 4.0463257258780339e1; 	cQ.w[35] = 1.0027987156316190e-18;
						cQ.x[36] = 4.2889980712201399e1; 	cQ.w[36] = 8.8958541452840810e-20;
						cQ.x[37] = 4.5398632459316878e1; 	cQ.w[37] = 7.2721294867218350e-21;
						cQ.x[38] = 4.7990711121944935e1; 	cQ.w[38] = 5.4700089499572571e-22;
						cQ.x[39] = 5.0667808778260416e1; 	cQ.w[39] = 3.7798889589128719e-23;
						cQ.x[40] = 5.3431618058147838e1; 	cQ.w[40] = 2.3955393929296669e-24;
						cQ.x[41] = 5.6283940014554144e1; 	cQ.w[41] = 1.3898975276677134e-25;
						cQ.x[42] = 5.9226692876193983e1; 	cQ.w[42] = 7.3686649638685550e-27;
						cQ.x[43] = 6.2261921804579421e1; 	cQ.w[43] = 3.5623617964926723e-28;
						cQ.x[44] = 6.5391809799470289e1; 	cQ.w[44] = 1.5670671394405450e-29;
						cQ.x[45] = 6.8618689922287050e1; 	cQ.w[45] = 6.2579202126857896e-31;
						cQ.x[46] = 7.1945059037830997e1; 	cQ.w[46] = 2.2630134382743232e-32;
						cQ.x[47] = 7.5373593312139102e1; 	cQ.w[47] = 7.3910051432137950e-34;
						cQ.x[48] = 7.8907165750162510e1; 	cQ.w[48] = 2.1738965407946768e-35;
						cQ.x[49] = 8.2548866113398157e1; 	cQ.w[49] = 5.7406301865352811e-37;
						cQ.x[50] = 8.6302023627490826e1; 	cQ.w[50] = 1.3565280389341554e-38;
						cQ.x[51] = 9.0170232976927819e1; 	cQ.w[51] = 2.8582178093556767e-40;
						cQ.x[52] = 9.4157384193266820e1; 	cQ.w[52] = 5.3491017715110987e-42;
						cQ.x[53] = 9.8267697181550113e1; 	cQ.w[53] = 8.8545200956130098e-44;
						cQ.x[54] = 1.0250576180568381e2; 	cQ.w[54] = 1.2905303360994193e-45;
						cQ.x[55] = 1.0687658467989708e2; 	cQ.w[55] = 1.6479024314408629e-47;
						cQ.x[56] = 1.1138564410689640e2; 	cQ.w[56] = 1.8335510486003441e-49;
						cQ.x[57] = 1.1603895498763616e2; 	cQ.w[57] = 1.7670959327066303e-51;
						cQ.x[58] = 1.2084314603612786e2; 	cQ.w[58] = 1.4654679930191558e-53;
						cQ.x[59] = 1.2580555231319460e2; 	cQ.w[59] = 1.0382012259739465e-55;
						cQ.x[60] = 1.3093432701495951e2; 	cQ.w[60] = 6.2325264712640765e-58;
						cQ.x[61] = 1.3623857771756320e2; 	cQ.w[61] = 3.1419744937690600e-60;
						cQ.x[62] = 1.4172853404292507e2; 	cQ.w[62] = 1.3167216299744360e-62;
						cQ.x[63] = 1.4741575620660480e2; 	cQ.w[63] = 4.5348612513885167e-65;
						cQ.x[64] = 1.5331339750560076e2; 	cQ.w[64] = 1.2669383082047862e-67;
						cQ.x[65] = 1.5943653908891510e2; 	cQ.w[65] = 2.8286860158479256e-70;
						cQ.x[66] = 1.6580262329068657e2; 	cQ.w[66] = 4.9608576653838684e-73;
						cQ.x[67] = 1.7243202402097285e2; 	cQ.w[67] = 6.6976362308197039e-76;
						cQ.x[68] = 1.7934881203693207e2; 	cQ.w[68] = 6.7974784724244015e-79;
						cQ.x[69] = 1.8658180447952091e2; 	cQ.w[69] = 5.0405265141181345e-82;
						cQ.x[70] = 1.9416604151140286e2; 	cQ.w[70] = 2.6380873260356514e-85;
						cQ.x[71] = 2.0214492732686698e2; 	cQ.w[71] = 9.3369329897351393e-89;
						cQ.x[72] = 2.1057344821155215e2; 	cQ.w[72] = 2.1169241731154879e-92;
						cQ.x[73] = 2.1952322632251781e2; 	cQ.w[73] = 2.8655194931946936e-96;
						cQ.x[74] = 2.2909090256827742e2; 	cQ.w[74] = 2.1061692586895640e-100;
						cQ.x[75] = 2.3941305411072180e2; 	cQ.w[75] = 7.6547423878179816e-105;
						cQ.x[76] = 2.5069535780481856e2; 	cQ.w[76] = 4.5949033294697062e-51;
						cQ.x[77] = 2.6327773413988054e2; 	cQ.w[77] = 6.4123032211485360e-55;
						cQ.x[78] = 2.7781337017790295e2; 	cQ.w[78] = 2.6355344341971242e-68;
						cQ.x[79] = 2.9599252372507157e2; 	cQ.w[79] = 4.2626228488835649e-81;
						break;
					case 88:
						cQ.x[0] = 6.9898229051547411e-3; 	cQ.w[0] = 3.3209344438409325e-1;
						cQ.x[1] = 6.2911728289785692e-2; 	cQ.w[1] = 3.1405674898440063e-1;
						cQ.x[2] = 1.7477326359221082e-1; 	cQ.w[2] = 2.8086402575909276e-1;
						cQ.x[3] = 3.4260990879552213e-1; 	cQ.w[3] = 2.3752476979765092e-1;
						cQ.x[4] = 5.6647496129652395e-1; 	cQ.w[4] = 1.8994306166644856e-1;
						cQ.x[5] = 8.4643962917804553e-1; 	cQ.w[5] = 1.4361800366234363e-1;
						cQ.x[6] = 1.1825931561845759e0; 	cQ.w[6] = 1.0266599963917014e-1;
						cQ.x[7] = 1.5750429789323641e0; 	cQ.w[7] = 6.9379346479544336e-2;
						cQ.x[8] = 2.0239149170256255e0; 	cQ.w[8] = 4.4316467105237175e-2;
						cQ.x[9] = 2.5293533968962740e0; 	cQ.w[9] = 2.6752780145097017e-2;
						cQ.x[10] = 3.0915217103368590e0; 	cQ.w[10] = 1.5260575426000781e-2;
						cQ.x[11] = 3.7106023088564040e0; 	cQ.w[11] = 8.2241712678339878e-3;
						cQ.x[12] = 4.3867971351580078e0; 	cQ.w[12] = 4.1864435774652050e-3;
						cQ.x[13] = 5.1203279932168699e0; 	cQ.w[13] = 2.0124922623294288e-3;
						cQ.x[14] = 5.9114369586294842e0; 	cQ.w[14] = 9.1338594776452588e-4;
						cQ.x[15] = 6.7603868311109156e0; 	cQ.w[15] = 3.9128377994473102e-4;
						cQ.x[16] = 7.6674616312393195e0; 	cQ.w[16] = 1.5816991465628900e-4;
						cQ.x[17] = 8.6329671437874178e0; 	cQ.w[17] = 6.0313998121755884e-5;
						cQ.x[18] = 9.6572315102419724e0; 	cQ.w[18] = 2.1688652721856188e-5;
						cQ.x[19] = 1.0740605873397186e1; 	cQ.w[19] = 7.3521662189359001e-6;
						cQ.x[20] = 1.1883465077219540e1; 	cQ.w[20] = 2.3485735170662032e-6;
						cQ.x[21] = 1.3086208425523377e1; 	cQ.w[21] = 7.0668575490339494e-7;
						cQ.x[22] = 1.4349260503372554e1; 	cQ.w[22] = 2.0021570985018051e-7;
						cQ.x[23] = 1.5673072065538298e1; 	cQ.w[23] = 5.3385656195939403e-8;
						cQ.x[24] = 1.7058120996802116e1; 	cQ.w[24] = 1.3390564636227622e-8;
						cQ.x[25] = 1.8504913349401250e1; 	cQ.w[25] = 3.1579273648368550e-9;
						cQ.x[26] = 2.0013984463479417e1; 	cQ.w[26] = 6.9984639734477468e-10;
						cQ.x[27] = 2.1585900177035253e1; 	cQ.w[27] = 1.4566522851312418e-10;
						cQ.x[28] = 2.3221258132563997e1; 	cQ.w[28] = 2.8457916194634976e-11;
						cQ.x[29] = 2.4920689188374744e1; 	cQ.w[29] = 5.2152057446839083e-12;
						cQ.x[30] = 2.6684858943448156e1; 	cQ.w[30] = 8.9592785232829909e-13;
						cQ.x[31] = 2.8514469385691634e1; 	cQ.w[31] = 1.4417973875868243e-13;
						cQ.x[32] = 3.0410260674566840e1; 	cQ.w[32] = 2.1719253795631441e-14;
						cQ.x[33] = 3.2373013070326923e1; 	cQ.w[33] = 3.0602546728763147e-15;
						cQ.x[34] = 3.4403549023529919e1; 	cQ.w[34] = 4.0298283180374064e-16;
						cQ.x[35] = 3.6502735440116345e1; 	cQ.w[35] = 4.9551512528428224e-17;
						cQ.x[36] = 3.8671486139183441e1; 	cQ.w[36] = 5.6842543746643022e-18;
						cQ.x[37] = 4.0910764522691699e1; 	cQ.w[37] = 6.0774264422552206e-19;
						cQ.x[38] = 4.3221586478743604e1; 	cQ.w[38] = 6.0500123323576681e-20;
						cQ.x[39] = 4.5605023542830377e1; 	cQ.w[39] = 5.6016998994982451e-21;
						cQ.x[40] = 4.8062206344609867e1; 	cQ.w[40] = 4.8186081446164731e-22;
						cQ.x[41] = 5.0594328371429292e1; 	cQ.w[41] = 3.8463215773862006e-23;
						cQ.x[42] = 5.3202650084026280e1; 	cQ.w[42] = 2.8454232488178569e-24;
						cQ.x[43] = 5.5888503424734049e1; 	cQ.w[43] = 1.9482724641523455e-25;
						cQ.x[44] = 5.8653296764206622e1; 	cQ.w[44] = 1.2329494549301778e-26;
						cQ.x[45] = 6.1498520339319444e1; 	cQ.w[45] = 7.2009452095564265e-28;
						cQ.x[46] = 6.4425752242674225e1; 	cQ.w[46] = 3.8752556618667624e-29;
						cQ.x[47] = 6.7436665033270514e1; 	cQ.w[47] = 1.9184700627092652e-30;
						cQ.x[48] = 7.0533033048677930e1; 	cQ.w[48] = 8.7214115374171135e-32;
						cQ.x[49] = 7.3716740511795214e1; 	cQ.w[49] = 3.6339711712303267e-33;
						cQ.x[50] = 7.6989790540440737e1; 	cQ.w[50] = 1.3850762473413167e-34;
						cQ.x[51] = 8.0354315186114501e1; 	cQ.w[51] = 4.8188380976506615e-36;
						cQ.x[52] = 8.3812586649969540e1; 	cQ.w[52] = 1.5268836098287974e-37;
						cQ.x[53] = 8.7367029850170202e1; 	cQ.w[53] = 4.3955893390762546e-39;
						cQ.x[54] = 9.1020236546461036e1; 	cQ.w[54] = 1.1467178326064145e-40;
						cQ.x[55] = 9.4774981266282465e1; 	cQ.w[55] = 2.7034989251018274e-42;
						cQ.x[56] = 9.8634239323895773e1; 	cQ.w[56] = 5.7430338090994212e-44;
						cQ.x[57] = 1.0260120728198211e2; 	cQ.w[57] = 1.0957767689227991e-45;
						cQ.x[58] = 1.0667932627700808e2; 	cQ.w[58] = 1.8714704533178499e-47;
						cQ.x[59] = 1.1087230871918132e2; 	cQ.w[59] = 2.8505049132116738e-49;
						cQ.x[60] = 1.1518416899019180e2; 	cQ.w[60] = 3.8566146717396706e-51;
						cQ.x[61] = 1.1961925890402067e2; 	cQ.w[61] = 4.6148622579363381e-53;
						cQ.x[62] = 1.2418230887717554e2; 	cQ.w[62] = 4.8611462511088885e-55;
						cQ.x[63] = 1.2887847598742982e2; 	cQ.w[63] = 4.4845808683072384e-57;
						cQ.x[64] = 1.3371340040194700e2; 	cQ.w[64] = 3.6030988117615386e-59;
						cQ.x[65] = 1.3869327205088259e2; 	cQ.w[65] = 2.5057229445134590e-61;
						cQ.x[66] = 1.4382490994553306e2; 	cQ.w[66] = 1.4981434219579486e-63;
						cQ.x[67] = 1.4911585724002039e2; 	cQ.w[67] = 7.6434036433722766e-66;
						cQ.x[68] = 1.5457449608379523e2; 	cQ.w[68] = 3.3000635664912984e-68;
						cQ.x[69] = 1.6021018761433266e2; 	cQ.w[69] = 1.1946162455653543e-70;
						cQ.x[70] = 1.6603344425357073e2; 	cQ.w[70] = 3.5882216086453966e-73;
						cQ.x[71] = 1.7205614404012268e2; 	cQ.w[71] = 8.8381795575894759e-76;
						cQ.x[72] = 1.7829180043052003e2; 	cQ.w[72] = 1.7614304684723708e-78;
						cQ.x[73] = 1.8475590644174345e2; 	cQ.w[73] = 2.7972091009100451e-81;
						cQ.x[74] = 1.9146638017638173e2; 	cQ.w[74] = 3.4772790784565961e-84;
						cQ.x[75] = 1.9844415134556974e2; 	cQ.w[75] = 3.3145008476708899e-87;
						cQ.x[76] = 2.0571394830191976e2; 	cQ.w[76] = 2.3639811252252232e-90;
						cQ.x[77] = 2.1330537759040182e2; 	cQ.w[77] = 1.2252220679209011e-93;
						cQ.x[78] = 2.2125444306253341e2; 	cQ.w[78] = 4.4534409258766820e-97;
						cQ.x[79] = 2.2960574884915661e2; 	cQ.w[79] = 1.0863817141918454e-100;
						cQ.x[80] = 2.3841581114398289e2; 	cQ.w[80] = 1.6822680293726604e-104;
						cQ.x[81] = 2.4775826014248507e2; 	cQ.w[81] = 1.5380603201854840e-108;
						cQ.x[82] = 2.5773247037077271e2; 	cQ.w[82] = 7.5305256625793262e-113;
						cQ.x[83] = 2.6847892170878440e2; 	cQ.w[83] = 1.7204742760900275e-117;
						cQ.x[84] = 2.8020923651927812e2; 	cQ.w[84] = 2.4712355922645459e-52;
						cQ.x[85] = 2.9327329106785954e2; 	cQ.w[85] = 2.0946397303280510e-50;
						cQ.x[86] = 3.0834369296097073e2; 	cQ.w[86] = 9.0428276327931568e-65;
						cQ.x[87] = 3.2716185523285862e2; 	cQ.w[87] = 1.9175767652062620e-79;
						break;
					case 96:
						cQ.x[0] = 6.4088479693086173e-3; 	cQ.w[0] = 3.1817720207314518e-1;
						cQ.x[1] = 5.7682192384603139e-2; 	cQ.w[1] = 3.0229450239014958e-1;
						cQ.x[2] = 1.6024254224241131e-1; 	cQ.w[2] = 2.7286429856942241e-1;
						cQ.x[3] = 3.1411723962617867e-1; 	cQ.w[3] = 2.3399502554972005e-1;
						cQ.x[4] = 5.1934734780450853e-1; 	cQ.w[4] = 1.9063041737616373e-1;
						cQ.x[5] = 7.7598771160231075e-1; 	cQ.w[5] = 1.4752978061991784e-1;
						cQ.x[6] = 1.0841070382287637e0; 	cQ.w[6] = 1.0845245228344696e-1;
						cQ.x[7] = 1.4437879988503862e0; 	cQ.w[7] = 7.5724332195484464e-2;
						cQ.x[8] = 1.8551273512731066e0; 	cQ.w[8] = 5.0214141133917932e-2;
						cQ.x[9] = 2.3182360841752035e0; 	cQ.w[9] = 3.1620089370995634e-2;
						cQ.x[10] = 2.8332395834139107e0; 	cQ.w[10] = 1.8905681570520517e-2;
						cQ.x[11] = 3.4002778210128949e0; 	cQ.w[11] = 1.0731321018794248e-2;
						cQ.x[12] = 4.0195055675263257e0; 	cQ.w[12] = 5.7820073860565274e-3;
						cQ.x[13] = 4.6910926285685173e0; 	cQ.w[13] = 2.9566179089772504e-3;
						cQ.x[14] = 5.4152241063968240e0; 	cQ.w[14] = 1.4345733126805138e-3;
						cQ.x[15] = 6.1921006875403740e0; 	cQ.w[15] = 6.6035034249772304e-4;
						cQ.x[16] = 7.0219389575791665e0; 	cQ.w[16] = 2.8830782856393540e-4;
						cQ.x[17] = 7.9049717442979328e0; 	cQ.w[17] = 1.1936254041747660e-4;
						cQ.x[18] = 8.8414484905679871e0; 	cQ.w[18] = 4.6849104586421198e-5;
						cQ.x[19] = 9.8316356584491658e0; 	cQ.w[19] = 1.7427730160295739e-5;
						cQ.x[20] = 1.0875817166154111e1; 	cQ.w[20] = 6.1427630640308088e-6;
						cQ.x[21] = 1.1974294859679961e1; 	cQ.w[21] = 2.0508762851384060e-6;
						cQ.x[22] = 1.3127389021089503e1; 	cQ.w[22] = 6.4837954271472021e-7;
						cQ.x[23] = 1.4335438915616705e1; 	cQ.w[23] = 1.9403784928607733e-7;
						cQ.x[24] = 1.5598803379982271e1; 	cQ.w[24] = 5.4948470407217337e-8;
						cQ.x[25] = 1.6917861454535469e1; 	cQ.w[25] = 1.4718810972300389e-8;
						cQ.x[26] = 1.8293013062091624e1; 	cQ.w[26] = 3.7279041518397285e-9;
						cQ.x[27] = 1.9724679736612835e1; 	cQ.w[27] = 8.9237873517832415e-10;
						cQ.x[28] = 2.1213305405186056e1; 	cQ.w[28] = 2.0180598524189019e-10;
						cQ.x[29] = 2.2759357227090996e1; 	cQ.w[29] = 4.3094025836606105e-11;
						cQ.x[30] = 2.4363326494124491e1; 	cQ.w[30] = 8.6853173406217016e-12;
						cQ.x[31] = 2.6025729596762688e1; 	cQ.w[31] = 1.6512638032485222e-12;
						cQ.x[32] = 2.7747109061202697e1; 	cQ.w[32] = 2.9598835329453310e-13;
						cQ.x[33] = 2.9528034662837436e1; 	cQ.w[33] = 4.9993510953075103e-14;
						cQ.x[34] = 3.1369104622288073e1; 	cQ.w[34] = 7.9519688037928626e-15;
						cQ.x[35] = 3.3270946890755629e1; 	cQ.w[35] = 1.1903810906669537e-15;
						cQ.x[36] = 3.5234220532166221e1; 	cQ.w[36] = 1.6759557500075327e-16;
						cQ.x[37] = 3.7259617210383481e1; 	cQ.w[37] = 2.2177113226367451e-17;
						cQ.x[38] = 3.9347862790659313e1; 	cQ.w[38] = 2.7561280628581385e-18;
						cQ.x[39] = 4.1499719065504354e1; 	cQ.w[39] = 3.2145199936756972e-19;
						cQ.x[40] = 4.3715985616298953e1; 	cQ.w[40] = 3.5156800159203487e-20;
						cQ.x[41] = 4.5997501823253394e1; 	cQ.w[41] = 3.6025990545465830e-21;
						cQ.x[42] = 4.8345149037785020e1; 	cQ.w[42] = 3.4558483324761214e-22;
						cQ.x[43] = 5.0759852933036353e1; 	cQ.w[43] = 3.1004532672815943e-23;
						cQ.x[44] = 5.3242586050143369e1; 	cQ.w[44] = 2.5990033771037196e-24;
						cQ.x[45] = 5.5794370560013523e1; 	cQ.w[45] = 2.0335626660738988e-25;
						cQ.x[46] = 5.8416281262832383e1; 	cQ.w[46] = 1.4835853103395276e-26;
						cQ.x[47] = 6.1109448850337470e1; 	cQ.w[47] = 1.0080516345867198e-27;
						cQ.x[48] = 6.3875063459139563e1; 	cQ.w[48] = 6.3716760141675274e-29;
						cQ.x[49] = 6.6714378547108836e1; 	cQ.w[49] = 3.7418231588440592e-30;
						cQ.x[50] = 6.9628715129163805e1; 	cQ.w[50] = 2.0389207136140223e-31;
						cQ.x[51] = 7.2619466413811364e1; 	cQ.w[51] = 1.0294436657754091e-32;
						cQ.x[52] = 7.5688102887614222e1; 	cQ.w[52] = 4.8089946872814419e-34;
						cQ.x[53] = 7.8836177901563518e1; 	cQ.w[53] = 2.0753148425343429e-35;
						cQ.x[54] = 8.2065333821298500e1; 	cQ.w[54] = 8.2600326269673452e-37;
						cQ.x[55] = 8.5377308812473848e1; 	cQ.w[55] = 3.0268879909029589e-38;
						cQ.x[56] = 8.8773944343613222e1; 	cQ.w[56] = 1.0193701477126852e-39;
						cQ.x[57] = 9.2257193501856653e1; 	cQ.w[57] = 3.1487953377266536e-41;
						cQ.x[58] = 9.5829130232545822e1; 	cQ.w[58] = 8.9030283898555212e-43;
						cQ.x[59] = 9.9491959632139456e1; 	cQ.w[59] = 2.2991067830207375e-44;
						cQ.x[60] = 1.0324802944619364e2; 	cQ.w[60] = 5.4099632013102152e-46;
						cQ.x[61] = 1.0709984295093964e2; 	cQ.w[61] = 1.1570760541533550e-47;
						cQ.x[62] = 1.1105007342943756e2; 	cQ.w[62] = 2.2434027937343383e-49;
						cQ.x[63] = 1.1510158049277213e2; 	cQ.w[63] = 3.9318356145233163e-51;
						cQ.x[64] = 1.1925742854508090e2; 	cQ.w[64] = 6.2101615333023767e-53;
						cQ.x[65] = 1.2352090775068523e2; 	cQ.w[65] = 8.8106771041998731e-55;
						cQ.x[66] = 1.2789555793525788e2; 	cQ.w[66] = 1.1188872437683772e-56;
						cQ.x[67] = 1.3238519594478634e2; 	cQ.w[67] = 1.2670287078247711e-58;
						cQ.x[68] = 1.3699394710135137e2; 	cQ.w[68] = 1.2741754871478997e-60;
						cQ.x[69] = 1.4172628154048866e2; 	cQ.w[69] = 1.1328848065466888e-62;
						cQ.x[70] = 1.4658705640065854e2; 	cQ.w[70] = 8.8625846243039034e-65;
						cQ.x[71] = 1.5158156507410060e2; 	cQ.w[71] = 6.0683384585205647e-67;
						cQ.x[72] = 1.5671559503799880e2; 	cQ.w[72] = 3.6159126920762266e-69;
						cQ.x[73] = 1.6199549619039953e2; 	cQ.w[73] = 1.8632314937659433e-71;
						cQ.x[74] = 1.6742826215197814e2; 	cQ.w[74] = 8.2452068773611101e-74;
						cQ.x[75] = 1.7302162771302658e2; 	cQ.w[75] = 3.1094858299677929e-76;
						cQ.x[76] = 1.7878418657827880e2; 	cQ.w[76] = 9.9088925545678630e-79;
						cQ.x[77] = 1.8472553489864371e2; 	cQ.w[77] = 2.6428813129457523e-81;
						cQ.x[78] = 1.9085644794111835e2; 	cQ.w[78] = 5.8372225441875102e-84;
						cQ.x[79] = 1.9718909988484638e2; 	cQ.w[79] = 1.0548005520682382e-86;
						cQ.x[80] = 2.0373734053112526e2; 	cQ.w[80] = 1.5381913402837248e-89;
						cQ.x[81] = 2.1051704829955498e2; 	cQ.w[81] = 1.7819545894948543e-92;
						cQ.x[82] = 2.1754658727070438e2; 	cQ.w[82] = 1.6104009725265265e-95;
						cQ.x[83] = 2.2484740894787622e2; 	cQ.w[83] = 1.1114881452557864e-98;
						cQ.x[84] = 2.3244485984497173e2; 	cQ.w[84] = 5.7137871484250277e-102;
						cQ.x[85] = 2.4036928938394857e2; 	cQ.w[85] = 2.1230825957651528e-105;
						cQ.x[86] = 2.4865760911996251e2; 	cQ.w[86] = 5.4979966939612256e-109;
						cQ.x[87] = 2.5735555421768465e2; 	cQ.w[87] = 9.4850943127338216e-113;
						cQ.x[88] = 2.6652108371132283e2; 	cQ.w[88] = 1.0296088521328023e-116;
						cQ.x[89] = 2.7622972228180378e2; 	cQ.w[89] = 6.4585972823499812e-121;
						cQ.x[90] = 2.8658342409205652e2; 	cQ.w[90] = 1.7529595043807457e-117;
						cQ.x[91] = 2.9772635225359436e2; 	cQ.w[91] = 3.1105665306320562e-130;
						cQ.x[92] = 3.0987574005931242e2; 	cQ.w[92] = 2.5102657022454895e-50;
						cQ.x[93] = 3.2339085777326039e2; 	cQ.w[93] = 1.1598746411740026e-54;
						cQ.x[94] = 3.3896263242081578e2; 	cQ.w[94] = 1.8735722301688603e-68;
						cQ.x[95] = 3.5838071121369978e2; 	cQ.w[95] = 3.6413912508706140e-80;
						break;
					case 104:
						cQ.x[0] = 5.9170399902629313e-3; 	cQ.w[0] = 3.0587540675027521e-1;
						cQ.x[1] = 5.3255375119357608e-2; 	cQ.w[1] = 2.9174971091826186e-1;
						cQ.x[2] = 1.4794279594504654e-1; 	cQ.w[2] = 2.6542237352386527e-1;
						cQ.x[3] = 2.9000081703620763e-1; 	cQ.w[3] = 2.3031197872444814e-1;
						cQ.x[4] = 4.7946174387139034e-1; 	cQ.w[4] = 1.9060462767137745e-1;
						cQ.x[5] = 7.1636871330281893e-1; 	cQ.w[5] = 1.5044217126994679e-1;
						cQ.x[6] = 1.0007757476986113e0; 	cQ.w[6] = 1.1324045142276644e-1;
						cQ.x[7] = 1.3327478229276096e0; 	cQ.w[7] = 8.1283295752663497e-2;
						cQ.x[8] = 1.7123609503940247e0; 	cQ.w[8] = 5.5633427651614899e-2;
						cQ.x[9] = 2.1397022733730653e0; 	cQ.w[9] = 3.6305080311286526e-2;
						cQ.x[10] = 2.6148701779441033e0; 	cQ.w[10] = 2.2586734959646651e-2;
						cQ.x[11] = 3.1379744188649923e0; 	cQ.w[11] = 1.3395115370766690e-2;
						cQ.x[12] = 3.7091362607801872e0; 	cQ.w[12] = 7.5717088395729501e-3;
						cQ.x[13] = 4.3284886352066046e0; 	cQ.w[13] = 4.0788547539583762e-3;
						cQ.x[14] = 4.9961763137950506e0; 	cQ.w[14] = 2.0937112294538550e-3;
						cQ.x[15] = 5.7123560984218524e0; 	cQ.w[15] = 1.0239096611160486e-3;
						cQ.x[16] = 6.4771970287254551e0; 	cQ.w[16] = 4.7697940575163359e-4;
						cQ.x[17] = 7.2908806077665737e0; 	cQ.w[17] = 2.1161778351897556e-4;
						cQ.x[18] = 8.1536010465584800e0; 	cQ.w[18] = 8.9399513582116538e-5;
						cQ.x[19] = 9.0655655282866201e0; 	cQ.w[19] = 3.5954968465819779e-5;
						cQ.x[20] = 1.0026994493114549e1; 	cQ.w[20] = 1.3763471639247608e-5;
						cQ.x[21] = 1.1038121944556697e1; 	cQ.w[21] = 5.0135007842022910e-6;
						cQ.x[22] = 1.2099195778488418e1; 	cQ.w[22] = 1.7373642266442212e-6;
						cQ.x[23] = 1.3210478135960773e1; 	cQ.w[23] = 5.7261651269671511e-7;
						cQ.x[24] = 1.4372245781092441e1; 	cQ.w[24] = 1.7944867221251103e-7;
						cQ.x[25] = 1.5584790505424801e1; 	cQ.w[25] = 5.3455660388561925e-8;
						cQ.x[26] = 1.6848419560249673e1; 	cQ.w[26] = 1.5131821405551822e-8;
						cQ.x[27] = 1.8163456118553402e1; 	cQ.w[27] = 4.0690485475645235e-9;
						cQ.x[28] = 1.9530239768367304e1; 	cQ.w[28] = 1.0390864670042617e-9;
						cQ.x[29] = 2.0949127039474128e1; 	cQ.w[29] = 2.5189139018781987e-10;
						cQ.x[30] = 2.2420491965594827e1; 	cQ.w[30] = 5.7944986530631928e-11;
						cQ.x[31] = 2.3944726684371185e1; 	cQ.w[31] = 1.2644148134587503e-11;
						cQ.x[32] = 2.5522242077669684e1; 	cQ.w[32] = 2.6161135978115952e-12;
						cQ.x[33] = 2.7153468454962594e1; 	cQ.w[33] = 5.1301572545741211e-13;
						cQ.x[34] = 2.8838856282796113e1; 	cQ.w[34] = 9.5305252540215033e-14;
						cQ.x[35] = 3.0578876963635225e1; 	cQ.w[35] = 1.6765288843561763e-14;
						cQ.x[36] = 3.2374023667684009e1; 	cQ.w[36] = 2.7912570841914264e-15;
						cQ.x[37] = 3.4224812221621939e1; 	cQ.w[37] = 4.3960232489011666e-16;
						cQ.x[38] = 3.6131782058575501e1; 	cQ.w[38] = 6.5457287737442617e-17;
						cQ.x[39] = 3.8095497234064709e1; 	cQ.w[39] = 9.2097916367430484e-18;
						cQ.x[40] = 4.0116547513131397e1; 	cQ.w[40] = 1.2237146175390552e-18;
						cQ.x[41] = 4.2195549534376423e1; 	cQ.w[41] = 1.5345540324461345e-19;
						cQ.x[42] = 4.4333148057213359e1; 	cQ.w[42] = 1.8150011522883905e-20;
						cQ.x[43] = 4.6530017299294814e1; 	cQ.w[43] = 2.0233572060917363e-21;
						cQ.x[44] = 4.8786862371793691e1; 	cQ.w[44] = 2.1245357575136396e-22;
						cQ.x[45] = 5.1104420821036147e1; 	cQ.w[45] = 2.0995790452740547e-23;
						cQ.x[46] = 5.3483464285898315e1; 	cQ.w[46] = 1.9513862776804163e-24;
						cQ.x[47] = 5.5924800281409502e1; 	cQ.w[47] = 1.7043065745685083e-25;
						cQ.x[48] = 5.8429274120167462e1; 	cQ.w[48] = 1.3975904770003455e-26;
						cQ.x[49] = 6.0997770984486184e1; 	cQ.w[49] = 1.0751199957435543e-27;
						cQ.x[50] = 6.3631218163686521e1; 	cQ.w[50] = 7.7513543786225607e-29;
						cQ.x[51] = 6.6330587472631968e1; 	cQ.w[51] = 5.2326478018346559e-30;
						cQ.x[52] = 6.9096897869537813e1; 	cQ.w[52] = 3.3040599143163473e-31;
						cQ.x[53] = 7.1931218293279340e1; 	cQ.w[53] = 1.9493708377000895e-32;
						cQ.x[54] = 7.4834670742938227e1; 	cQ.w[54] = 1.0734353038941787e-33;
						cQ.x[55] = 7.7808433625208567e1; 	cQ.w[55] = 5.5103919609823204e-35;
						cQ.x[56] = 8.0853745398597933e1; 	cQ.w[56] = 2.6337753192506890e-36;
						cQ.x[57] = 8.3971908547179854e1; 	cQ.w[57] = 1.1705805971768911e-37;
						cQ.x[58] = 8.7164293921072095e1; 	cQ.w[58] = 4.8312277173311803e-39;
						cQ.x[59] = 9.0432345485938809e1; 	cQ.w[59] = 1.8489404848533130e-40;
						cQ.x[60] = 9.3777585529775213e1; 	cQ.w[60] = 6.5514799851540860e-42;
						cQ.x[61] = 9.7201620382189983e1; 	cQ.w[61] = 2.1459149564917285e-43;
						cQ.x[62] = 1.0070614670954704e2; 	cQ.w[62] = 6.4864747764771309e-45;
						cQ.x[63] = 1.0429295845890204e2; 	cQ.w[63] = 1.8061451197301461e-46;
						cQ.x[64] = 1.0796395453496108e2; 	cQ.w[64] = 4.6240601087164368e-48;
						cQ.x[65] = 1.1172114730766026e2; 	cQ.w[65] = 1.0863042774855916e-49;
						cQ.x[66] = 1.1556667206386112e2; 	cQ.w[66] = 2.3367590047772692e-51;
						cQ.x[67] = 1.1950279753563604e2; 	cQ.w[67] = 4.5923018372463456e-53;
						cQ.x[68] = 1.2353193766037873e2; 	cQ.w[68] = 8.2254183241711378e-55;
						cQ.x[69] = 1.2765666475539859e2; 	cQ.w[69] = 1.3393229733255628e-56;
						cQ.x[70] = 1.3187972432286344e2; 	cQ.w[70] = 1.9770897715441433e-58;
						cQ.x[71] = 1.3620405174137066e2; 	cQ.w[71] = 2.6382340385810416e-60;
						cQ.x[72] = 1.4063279114988850e2; 	cQ.w[72] = 3.1724019756972273e-62;
						cQ.x[73] = 1.4516931689069355e2; 	cQ.w[73] = 3.4260518574631090e-64;
						cQ.x[74] = 1.4981725795333762e2; 	cQ.w[74] = 3.3110336284062607e-66;
						cQ.x[75] = 1.5458052595568164e2; 	cQ.w[75] = 2.8523873910804969e-68;
						cQ.x[76] = 1.5946334731603723e2; 	cQ.w[76] = 2.1812446026410382e-70;
						cQ.x[77] = 1.6447030041968208e2; 	cQ.w[77] = 1.4739245558990479e-72;
						cQ.x[78] = 1.6960635877321616e2; 	cQ.w[78] = 8.7573978000268771e-75;
						cQ.x[79] = 1.7487694138470448e2; 	cQ.w[79] = 4.5505714966004102e-77;
						cQ.x[80] = 1.8028797192464818e2; 	cQ.w[80] = 2.0558561295409975e-79;
						cQ.x[81] = 1.8584594863812376e2; 	cQ.w[81] = 8.0232998204915563e-82;
						cQ.x[82] = 1.9155802752806165e2; 	cQ.w[82] = 2.6857246723649169e-84;
						cQ.x[83] = 1.9743212206533502e2; 	cQ.w[83] = 7.6508024905351989e-87;
						cQ.x[84] = 2.0347702367824598e2; 	cQ.w[84] = 1.8386578992442195e-89;
						cQ.x[85] = 2.0970254864305275e2; 	cQ.w[85] = 3.6915967746713164e-92;
						cQ.x[86] = 2.1611971890494698e2; 	cQ.w[86] = 6.1249132399226474e-95;
						cQ.x[87] = 2.2274098706028915e2; 	cQ.w[87] = 8.2946217677700358e-98;
						cQ.x[88] = 2.2958051962428511e2; 	cQ.w[88] = 9.0408055705882247e-101;
						cQ.x[89] = 2.3665455843056904e2; 	cQ.w[89] = 7.8044833919279015e-104;
						cQ.x[90] = 2.4398188860512600e2; 	cQ.w[90] = 5.2375692667980381e-107;
						cQ.x[91] = 2.5158445479007743e2; 	cQ.w[91] = 2.6738203410490589e-110;
						cQ.x[92] = 2.5948818823677792e2; 	cQ.w[92] = 1.0120749135417976e-113;
						cQ.x[93] = 2.6772414159918294e2; 	cQ.w[93] = 2.7544888677231600e-117;
						cQ.x[94] = 2.7633008621234766e2; 	cQ.w[94] = 5.1929686996596160e-121;
						cQ.x[95] = 2.8535282906339674e2; 	cQ.w[95] = 6.4755792717472659e-125;
						cQ.x[96] = 2.9485169696486133e2; 	cQ.w[96] = 5.0218621637969880e-129;
						cQ.x[97] = 3.0490401093647767e2; 	cQ.w[97] = 1.9758139579198723e-133;
						cQ.x[98] = 3.1561417142961915e2; 	cQ.w[98] = 2.6170200395803803e-134;
						cQ.x[99] = 3.2712983445040887e2; 	cQ.w[99] = 2.0394890255436422e-51;
						cQ.x[100] = 3.3967355383586156e2; 	cQ.w[100] = 1.9410589620744451e-56;
						cQ.x[101] = 3.5361350610111344e2; 	cQ.w[101] = 1.9061593758307076e-70;
						cQ.x[102] = 3.6965798176808293e2; 	cQ.w[102] = 9.4662691701560053e-80;
						cQ.x[103] = 3.8964232774217858e2; 	cQ.w[103] = 1.1104725813752798e-82;
						break;
					case 112:
						cQ.x[0] = 5.4953341803301609e-3; 	cQ.w[0] = 2.9489826922398073e-1;
						cQ.x[1] = 4.9459621920962637e-2; 	cQ.w[1] = 2.8222796178145113e-1;
						cQ.x[2] = 1.3739680892391444e-1; 	cQ.w[2] = 2.5849487746587741e-1;
						cQ.x[3] = 2.6932412751335894e-1; 	cQ.w[3] = 2.2657942282144552e-1;
						cQ.x[4] = 4.4526744940170549e-1; 	cQ.w[4] = 1.9006161252137982e-1;
						cQ.x[5] = 6.6526131362871422e-1; 	cQ.w[5] = 1.5256636998998856e-1;
						cQ.x[6] = 9.2934896392688816e-1; 	cQ.w[6] = 1.1719118584424322e-1;
						cQ.x[7] = 1.2375823956108998e0; 	cQ.w[7] = 8.6135126997093899e-2;
						cQ.x[8] = 1.5900224121141247e0; 	cQ.w[8] = 6.0574317417714514e-2;
						cQ.x[9] = 1.9867386913212612e0; 	cQ.w[9] = 4.0755803971739961e-2;
						cQ.x[10] = 2.4278098618726480e0; 	cQ.w[10] = 2.6233066082736348e-2;
						cQ.x[11] = 2.9133235896433742e0; 	cQ.w[11] = 1.6152112426163087e-2;
						cQ.x[12] = 3.4433766746287615e0; 	cQ.w[12] = 9.5123620147754134e-3;
						cQ.x[13] = 4.0180751584974267e0; 	cQ.w[13] = 5.3577229622045945e-3;
						cQ.x[14] = 4.6375344431040623e0; 	cQ.w[14] = 2.8857253991993310e-3;
						cQ.x[15] = 5.3018794202864780e0; 	cQ.w[15] = 1.4861359153491582e-3;
						cQ.x[16] = 6.0112446133054996e0; 	cQ.w[16] = 7.3169919995809469e-4;
						cQ.x[17] = 6.7657743303222200e0; 	cQ.w[17] = 3.4436155411986205e-4;
						cQ.x[18] = 7.5656228303450516e0; 	cQ.w[18] = 1.5489520162150497e-4;
						cQ.x[19] = 8.4109545021192616e0; 	cQ.w[19] = 6.6578086192592747e-5;
						cQ.x[20] = 9.3019440564744304e0; 	cQ.w[20] = 2.7341291782963178e-5;
						cQ.x[21] = 1.0238776732690826e1; 	cQ.w[21] = 1.0725583190203260e-5;
						cQ.x[22] = 1.1221648519494320e1; 	cQ.w[22] = 4.0183900973494323e-6;
						cQ.x[23] = 1.2250766391341530e1; 	cQ.w[23] = 1.4375499392088491e-6;
						cQ.x[24] = 1.3326348560712631e1; 	cQ.w[24] = 4.9095191796603443e-7;
						cQ.x[25] = 1.4448624747189284e1; 	cQ.w[25] = 1.6002973367034795e-7;
						cQ.x[26] = 1.5617836464159565e1; 	cQ.w[26] = 4.9774112595986853e-8;
						cQ.x[27] = 1.6834237324061373e1; 	cQ.w[27] = 1.4768543773235768e-8;
						cQ.x[28] = 1.8098093363150867e1; 	cQ.w[28] = 4.1791497540829771e-9;
						cQ.x[29] = 1.9409683386863716e1; 	cQ.w[29] = 1.1275442160204611e-9;
						cQ.x[30] = 2.0769299336924940e1; 	cQ.w[30] = 2.8996675589492684e-10;
						cQ.x[31] = 2.2177246681458568e1; 	cQ.w[31] = 7.1055803679543191e-11;
						cQ.x[32] = 2.3633844829452103e1; 	cQ.w[32] = 1.6586292637190786e-11;
						cQ.x[33] = 2.5139427571043608e1; 	cQ.w[33] = 3.6868289748704489e-12;
						cQ.x[34] = 2.6694343545222255e1; 	cQ.w[34] = 7.8011677084046385e-13;
						cQ.x[35] = 2.8298956736667379e1; 	cQ.w[35] = 1.5707664005052362e-13;
						cQ.x[36] = 2.9953647003597728e1; 	cQ.w[36] = 3.0084617930185426e-14;
						cQ.x[37] = 3.1658810638663093e1; 	cQ.w[37] = 5.4788162867425923e-15;
						cQ.x[38] = 3.3414860965086388e1; 	cQ.w[38] = 9.4832975673093346e-16;
						cQ.x[39] = 3.5222228970457193e1; 	cQ.w[39] = 1.5594657436863047e-16;
						cQ.x[40] = 3.7081363980789887e1; 	cQ.w[40] = 2.4352390869127641e-17;
						cQ.x[41] = 3.8992734377692841e1; 	cQ.w[41] = 3.6095571137034836e-18;
						cQ.x[42] = 4.0956828361752369e1; 	cQ.w[42] = 5.0757642386344912e-19;
						cQ.x[43] = 4.2974154765518962e1; 	cQ.w[43] = 6.7680555851412924e-20;
						cQ.x[44] = 4.5045243919797104e1; 	cQ.w[44] = 8.5528682924618337e-21;
						cQ.x[45] = 4.7170648577287264e1; 	cQ.w[45] = 1.0237779829567897e-21;
						cQ.x[46] = 4.9350944898013731e1; 	cQ.w[46] = 1.1601029671831009e-22;
						cQ.x[47] = 5.1586733501399526e1; 	cQ.w[47] = 1.2437241818417510e-23;
						cQ.x[48] = 5.3878640590325199e1; 	cQ.w[48] = 1.2607163390079573e-24;
						cQ.x[49] = 5.6227319153038116e1; 	cQ.w[49] = 1.2075186187234281e-25;
						cQ.x[50] = 5.8633450249370054e1; 	cQ.w[50] = 1.0920906578859209e-26;
						cQ.x[51] = 6.1097744388381858e1; 	cQ.w[51] = 9.3197543545708300e-28;
						cQ.x[52] = 6.3620943005294023e1; 	cQ.w[52] = 7.4991199467470311e-29;
						cQ.x[53] = 6.6203820046392479e1; 	cQ.w[53] = 5.6851415515118275e-30;
						cQ.x[54] = 6.8847183671532207e1; 	cQ.w[54] = 4.0573989053601805e-31;
						cQ.x[55] = 7.1551878084912556e1; 	cQ.w[55] = 2.7237314698332196e-32;
						cQ.x[56] = 7.4318785505984439e1; 	cQ.w[56] = 1.7183408111964119e-33;
						cQ.x[57] = 7.7148828293691083e1; 	cQ.w[57] = 1.0178497763141794e-34;
						cQ.x[58] = 8.0042971238764368e1; 	cQ.w[58] = 5.6554871509106191e-36;
						cQ.x[59] = 8.3002224040525493e1; 	cQ.w[59] = 2.9446351279717004e-37;
						cQ.x[60] = 8.6027643986604472e1; 	cQ.w[60] = 1.4351943180018244e-38;
						cQ.x[61] = 8.9120338856236006e1; 	cQ.w[61] = 6.5407638433990210e-40;
						cQ.x[62] = 9.2281470070355184e1; 	cQ.w[62] = 2.7840843888220453e-41;
						cQ.x[63] = 9.5512256114659083e1; 	cQ.w[63] = 1.1054661469676639e-42;
						cQ.x[64] = 9.8813976265183997e1; 	cQ.w[64] = 4.0894445809546643e-44;
						cQ.x[65] = 1.0218797464984955e2; 	cQ.w[65] = 1.4075301928520080e-45;
						cQ.x[66] = 1.0563566468393243e2; 	cQ.w[66] = 4.5010556553905317e-47;
						cQ.x[67] = 1.0915853392266498e2; 	cQ.w[67] = 1.3353375262485678e-48;
						cQ.x[68] = 1.1275814938024115e2; 	cQ.w[68] = 3.6695275610607160e-50;
						cQ.x[69] = 1.1643616337161760e2; 	cQ.w[69] = 9.3251810053654335e-52;
						cQ.x[70] = 1.2019431994181835e2; 	cQ.w[70] = 2.1876527359441781e-53;
						cQ.x[71] = 1.2403446195723057e2; 	cQ.w[71] = 4.7290816436201351e-55;
						cQ.x[72] = 1.2795853894491414e2; 	cQ.w[72] = 9.4017744393225470e-57;
						cQ.x[73] = 1.3196861577960725e2; 	cQ.w[73] = 1.7154770732070957e-58;
						cQ.x[74] = 1.3606688233435014e2; 	cQ.w[74] = 2.8665242731015198e-60;
						cQ.x[75] = 1.4025566423003982e2; 	cQ.w[75] = 4.3763989573541065e-62;
						cQ.x[76] = 1.4453743484248376e2; 	cQ.w[76] = 6.0897571742385102e-64;
						cQ.x[77] = 1.4891482875354178e2; 	cQ.w[77] = 7.7031170730276684e-66;
						cQ.x[78] = 1.5339065686687570e2; 	cQ.w[78] = 8.8328749793303574e-68;
						cQ.x[79] = 1.5796792345012664e2; 	cQ.w[79] = 9.1539586505379100e-70;
						cQ.x[80] = 1.6264984541588618e2; 	cQ.w[80] = 8.5466832132929860e-72;
						cQ.x[81] = 1.6743987421605101e2; 	cQ.w[81] = 7.1643610696800221e-74;
						cQ.x[82] = 1.7234172080122024e2; 	cQ.w[82] = 5.3721243887414992e-76;
						cQ.x[83] = 1.7735938419287581e2; 	cQ.w[83] = 3.5890274005464772e-78;
						cQ.x[84] = 1.8249718433670430e2; 	cQ.w[84] = 2.1271822587291959e-80;
						cQ.x[85] = 1.8775980005795738e2; 	cQ.w[85] = 1.1132951771397071e-82;
						cQ.x[86] = 1.9315231313418438e2; 	cQ.w[86] = 5.1191820200621500e-85;
						cQ.x[87] = 1.9868025975060599e2; 	cQ.w[87] = 2.0567744853280919e-87;
						cQ.x[88] = 2.0434969092759117e2; 	cQ.w[88] = 7.1772526649596168e-90;
						cQ.x[89] = 2.1016724393431546e2; 	cQ.w[89] = 2.1609888032718971e-92;
						cQ.x[90] = 2.1614022726467649e2; 	cQ.w[90] = 5.5733711965802891e-95;
						cQ.x[91] = 2.2227672250383794e2; 	cQ.w[91] = 1.2214353570041591e-97;
						cQ.x[92] = 2.2858570743324025e2; 	cQ.w[92] = 2.2544498870469226e-100;
						cQ.x[93] = 2.3507720612202630e2; 	cQ.w[93] = 3.4698364983238705e-103;
						cQ.x[94] = 2.4176247370399097e2; 	cQ.w[94] = 4.4037992320193089e-106;
						cQ.x[95] = 2.4865422630218547e2; 	cQ.w[95] = 4.5511410854658978e-109;
						cQ.x[96] = 2.5576693054575381e2; 	cQ.w[96] = 3.7753744367507348e-112;
						cQ.x[97] = 2.6311717297716247e2; 	cQ.w[97] = 2.4729177063662401e-115;
						cQ.x[98] = 2.7072413844178640e2; 	cQ.w[98] = 1.2549250228724467e-118;
						cQ.x[99] = 2.7861024009039942e2; 	cQ.w[99] = 4.8255791873729851e-122;
						cQ.x[100] = 2.8680196505406753e2; 	cQ.w[100] = 1.3697149157065236e-125;
						cQ.x[101] = 2.9533103485728541e2; 	cQ.w[101] = 2.7803242111394275e-129;
						cQ.x[102] = 3.0423603893997241e2; 	cQ.w[102] = 5.7192111543195200e-133;
						cQ.x[103] = 3.1356480447692787e2; 	cQ.w[103] = 2.7224139503928083e-135;
						cQ.x[104] = 3.2337796045230511e2; 	cQ.w[104] = 1.4985049959339292e-124;
						cQ.x[105] = 3.3375453828748506e2; 	cQ.w[105] = 2.4300435849929846e-132;
						cQ.x[106] = 3.4480126780658105e2; 	cQ.w[106] = 1.2259326780033627e-133;
						cQ.x[107] = 3.5666913087724963e2; 	cQ.w[107] = 2.4281860924963857e-51;
						cQ.x[108] = 3.6958574691724203e2; 	cQ.w[108] = 1.3766733748904064e-55;
						cQ.x[109] = 3.8392777037966436e2; 	cQ.w[109] = 5.9376811746370755e-69;
						cQ.x[110] = 4.0042001600670709e2; 	cQ.w[110] = 2.6857350954106024e-82;
						cQ.x[111] = 4.2094130637188611e2; 	cQ.w[111] = 2.9834955778444689e-81;
						break;
					case 120:
						cQ.x[0] = 5.1297391669498222e-3; 	cQ.w[0] = 2.8502395317547204e-1;
						cQ.x[1] = 4.6168965558965583e-2; 	cQ.w[1] = 2.7357523293518048e-1;
						cQ.x[2] = 1.2825442268012189e-1; 	cQ.w[2] = 2.5203716182373969e-1;
						cQ.x[3] = 2.5140012578184619e-1; 	cQ.w[3] = 2.2286337164152000e-1;
						cQ.x[4] = 4.1562711419913549e-1; 	cQ.w[4] = 1.8914288655924106e-1;
						cQ.x[5] = 6.2096347114179658e-1; 	cQ.w[5] = 1.5406584277742807e-1;
						cQ.x[6] = 8.6744435015264118e-1; 	cQ.w[6] = 1.2044049267022821e-1;
						cQ.x[7] = 1.1551120082929025e0; 	cQ.w[7] = 9.0358744219355181e-2;
						cQ.x[8] = 1.4840158461306763e0; 	cQ.w[8] = 6.5054591666537773e-2;
						cQ.x[9] = 1.8542124546240439e0; 	cQ.w[9] = 4.4943958036852311e-2;
						cQ.x[10] = 2.2657656690067761e0; 	cQ.w[10] = 2.9793575056951512e-2;
						cQ.x[11] = 2.7187466298012132e0; 	cQ.w[11] = 1.8949643702576734e-2;
						cQ.x[12] = 3.2132338511001456e0; 	cQ.w[12] = 1.1563055121878621e-2;
						cQ.x[13] = 3.7493132962773567e0; 	cQ.w[13] = 6.7686183318686482e-3;
						cQ.x[14] = 4.3270784613050146e0; 	cQ.w[14] = 3.8005184078890596e-3;
						cQ.x[15] = 4.9466304658754044e0; 	cQ.w[15] = 2.0467143746931182e-3;
						cQ.x[16] = 5.6080781525446617e0; 	cQ.w[16] = 1.0570513354181031e-3;
						cQ.x[17] = 6.3115381941373083e0; 	cQ.w[17] = 5.2349052163887461e-4;
						cQ.x[18] = 7.0571352096725975e0; 	cQ.w[18] = 2.4856586787767428e-4;
						cQ.x[19] = 7.8450018890970696e0; 	cQ.w[19] = 1.1314532553553562e-4;
						cQ.x[20] = 8.6752791271324176e0; 	cQ.w[20] = 4.9366594298931427e-5;
						cQ.x[21] = 9.5481161665738945e0; 	cQ.w[21] = 2.0642700128290312e-5;
						cQ.x[22] = 1.0463670751402208e1; 	cQ.w[22] = 8.2711958203927816e-6;
						cQ.x[23] = 1.1422109290101294e1; 	cQ.w[23] = 3.1751663936728944e-6;
						cQ.x[24] = 1.2423607029605695e1; 	cQ.w[24] = 1.1675760027782005e-6;
						cQ.x[25] = 1.3468348240334714e1; 	cQ.w[25] = 4.1119185727034140e-7;
						cQ.x[26] = 1.4556526412806197e1; 	cQ.w[26] = 1.3866304189070192e-7;
						cQ.x[27] = 1.5688344466361026e1; 	cQ.w[27] = 4.4765756397348091e-8;
						cQ.x[28] = 1.6864014970570331e1; 	cQ.w[28] = 1.3832734889608059e-8;
						cQ.x[29] = 1.8083760379941373e1; 	cQ.w[29] = 4.0902585285705370e-9;
						cQ.x[30] = 1.9347813282585301e1; 	cQ.w[30] = 1.1571064421471916e-9;
						cQ.x[31] = 2.0656416663560776e1; 	cQ.w[31] = 3.1309092427946476e-10;
						cQ.x[32] = 2.2009824183662284e1; 	cQ.w[32] = 8.1008838545649928e-11;
						cQ.x[33] = 2.3408300474481054e1; 	cQ.w[33] = 2.0037543940941961e-11;
						cQ.x[34] = 2.4852121450630391e1; 	cQ.w[34] = 4.7368238655988873e-12;
						cQ.x[35] = 2.6341574640096350e1; 	cQ.w[35] = 1.0698826247166949e-12;
						cQ.x[36] = 2.7876959533749559e1; 	cQ.w[36] = 2.3081364934421758e-13;
						cQ.x[37] = 2.9458587955135198e1; 	cQ.w[37] = 4.7547562946209405e-14;
						cQ.x[38] = 3.1086784451746349e1; 	cQ.w[38] = 9.3496621251942795e-15;
						cQ.x[39] = 3.2761886709081786e1; 	cQ.w[39] = 1.7543580097101897e-15;
						cQ.x[40] = 3.4484245988893703e1; 	cQ.w[40] = 3.1400971154418774e-16;
						cQ.x[41] = 3.6254227593144597e1; 	cQ.w[41] = 5.3593522678761522e-17;
						cQ.x[42] = 3.8072211355316697e1; 	cQ.w[42] = 8.7188945766817631e-18;
						cQ.x[43] = 3.9938592160852983e1; 	cQ.w[43] = 1.3515106234226978e-18;
						cQ.x[44] = 4.1853780498657202e1; 	cQ.w[44] = 1.9953021921306626e-19;
						cQ.x[45] = 4.3818203045742856e1; 	cQ.w[45] = 2.8044282474471764e-20;
						cQ.x[46] = 4.5832303287299388e1; 	cQ.w[46] = 3.7508906971614997e-21;
						cQ.x[47] = 4.7896542174639558e1; 	cQ.w[47] = 4.7717681635147569e-22;
						cQ.x[48] = 5.0011398823707265e1; 	cQ.w[48] = 5.7712738310619112e-23;
						cQ.x[49] = 5.2177371257062108e1; 	cQ.w[49] = 6.6327836836945849e-24;
						cQ.x[50] = 5.4394977192518312e1; 	cQ.w[50] = 7.2398182117815891e-25;
						cQ.x[51] = 5.6664754881904218e1; 	cQ.w[51] = 7.5012587727931265e-26;
						cQ.x[52] = 5.8987264003727612e1; 	cQ.w[52] = 7.3734885029706274e-27;
						cQ.x[53] = 6.1363086613885492e1; 	cQ.w[53] = 6.8721443844602975e-28;
						cQ.x[54] = 6.3792828158948708e1; 	cQ.w[54] = 6.0691779183687126e-29;
						cQ.x[55] = 6.6277118556987132e1; 	cQ.w[55] = 5.0759007367076133e-30;
						cQ.x[56] = 6.8816613351385211e1; 	cQ.w[56] = 4.0175221438220706e-31;
						cQ.x[57] = 7.1411994943637213e1; 	cQ.w[57] = 3.0072530813297540e-32;
						cQ.x[58] = 7.4063973911713684e1; 	cQ.w[58] = 2.1273584216828318e-33;
						cQ.x[59] = 7.6773290421263832e1; 	cQ.w[59] = 1.4211829797498313e-34;
						cQ.x[60] = 7.9540715737672631e1; 	cQ.w[60] = 8.9590957177462511e-36;
						cQ.x[61] = 8.2367053847837527e1; 	cQ.w[61] = 5.3251918761543396e-37;
						cQ.x[62] = 8.5253143201480774e1; 	cQ.w[62] = 2.9819503451853612e-38;
						cQ.x[63] = 8.8199858582884791e1; 	cQ.w[63] = 1.5717412988204925e-39;
						cQ.x[64] = 9.1208113125147039e1; 	cQ.w[64] = 7.7908008505381426e-41;
						cQ.x[65] = 9.4278860480418384e1; 	cQ.w[65] = 3.6281971412784443e-42;
						cQ.x[66] = 9.7413097161138711e1; 	cQ.w[66] = 1.5859015177178711e-43;
						cQ.x[67] = 1.0061186506904395e2; 	cQ.w[67] = 6.4996104375459879e-45;
						cQ.x[68] = 1.0387625423072283e2; 	cQ.w[68] = 2.4948970341280962e-46;
						cQ.x[69] = 1.0720740576078862e2; 	cQ.w[69] = 8.9593901038397884e-48;
						cQ.x[70] = 1.1060651507634780e2; 	cQ.w[70] = 3.0064061869829462e-49;
						cQ.x[71] = 1.1407483538944766e2; 	cQ.w[71] = 9.4149529210783543e-51;
						cQ.x[72] = 1.1761368150763627e2; 	cQ.w[72] = 2.7480284196099378e-52;
						cQ.x[73] = 1.2122443397674619e2; 	cQ.w[73] = 7.4655144872656008e-54;
						cQ.x[74] = 1.2490854360461544e2; 	cQ.w[74] = 1.8849778633354656e-55;
						cQ.x[75] = 1.2866753640979513e2; 	cQ.w[75] = 4.4167295028366118e-57;
						cQ.x[76] = 1.3250301904550293e2; 	cQ.w[76] = 9.5884544111130115e-59;
						cQ.x[77] = 1.3641668475632833e2; 	cQ.w[77] = 1.9253846157160579e-60;
						cQ.x[78] = 1.4041031993368383e2; 	cQ.w[78] = 3.5697275921210464e-62;
						cQ.x[79] = 1.4448581134597201e2; 	cQ.w[79] = 6.0993788510648131e-64;
						cQ.x[80] = 1.4864515413120591e2; 	cQ.w[80] = 9.5853191608176058e-66;
						cQ.x[81] = 1.5289046065375615e2; 	cQ.w[81] = 1.3825627736640795e-67;
						cQ.x[82] = 1.5722397034346682e2; 	cQ.w[82] = 1.8262174400320119e-69;
						cQ.x[83] = 1.6164806065516658e2; 	cQ.w[83] = 2.2038494574755515e-71;
						cQ.x[84] = 1.6616525931033000e2; 	cQ.w[84] = 2.4237108166704223e-73;
						cQ.x[85] = 1.7077825801123656e2; 	cQ.w[85] = 2.4226238499885409e-75;
						cQ.x[86] = 1.7548992785259882e2; 	cQ.w[86] = 2.1946080740100573e-77;
						cQ.x[87] = 1.8030333669777762e2; 	cQ.w[87] = 1.7962565586872283e-79;
						cQ.x[88] = 1.8522176883828634e2; 	cQ.w[88] = 1.3240410388018822e-81;
						cQ.x[89] = 1.9024874731879046e2; 	cQ.w[89] = 8.7585968525902169e-84;
						cQ.x[90] = 1.9538805938846759e2; 	cQ.w[90] = 5.1800274294960421e-86;
						cQ.x[91] = 2.0064378563766165e2; 	cQ.w[91] = 2.7279222749641292e-88;
						cQ.x[92] = 2.0602033350188195e2; 	cQ.w[92] = 1.2735953172719079e-90;
						cQ.x[93] = 2.1152247597090631e2; 	cQ.w[93] = 5.2465045442237418e-93;
						cQ.x[94] = 2.1715539653923255e2; 	cQ.w[94] = 1.8971890062657397e-95;
						cQ.x[95] = 2.2292474168927554e2; 	cQ.w[95] = 5.9884763458065071e-98;
						cQ.x[96] = 2.2883668252968495e2; 	cQ.w[96] = 1.6399282759538974e-100;
						cQ.x[97] = 2.3489798764468268e2; 	cQ.w[97] = 3.8700481929789787e-103;
						cQ.x[98] = 2.4111610978413709e2; 	cQ.w[98] = 7.8122606315689856e-106;
						cQ.x[99] = 2.4749928979224946e2; 	cQ.w[99] = 1.3379893855632987e-108;
						cQ.x[100] = 2.5405668221374821e2; 	cQ.w[100] = 1.9266285323004829e-111;
						cQ.x[101] = 2.6079850844627440e2; 	cQ.w[101] = 2.3089140588011642e-114;
						cQ.x[102] = 2.6773624530027240e2; 	cQ.w[102] = 2.2768662844496908e-117;
						cQ.x[103] = 2.7488285964960641e2; 	cQ.w[103] = 1.8239060669490364e-120;
						cQ.x[104] = 2.8225310392367891e2; 	cQ.w[104] = 1.1696325620584706e-123;
						cQ.x[105] = 2.8986389317085953e2; 	cQ.w[105] = 5.9089203696565755e-127;
						cQ.x[106] = 2.9773479340583953e2; 	cQ.w[106] = 2.3982393520466646e-130;
						cQ.x[107] = 3.0588866478394910e2; 	cQ.w[107] = 1.3336608988109749e-133;
						cQ.x[108] = 3.1435252503773286e2; 	cQ.w[108] = 1.8554662333955921e-131;
						cQ.x[109] = 3.2315873437660234e2; 	cQ.w[109] = 5.6745552433723977e-133;
						cQ.x[110] = 3.3234666364534706e2; 	cQ.w[110] = 6.4148023503115847e-132;
						cQ.x[111] = 3.4196511464296833e2; 	cQ.w[111] = 1.7690546779468858e-132;
						cQ.x[112] = 3.5207596053730912e2; 	cQ.w[112] = 2.3559424075535023e-133;
						cQ.x[113] = 3.6275986709694166e2; 	cQ.w[113] = 7.0346860953952109e-133;
						cQ.x[114] = 3.7412578995211861e2; 	cQ.w[114] = 8.2581860285237120e-132;
						cQ.x[115] = 3.8632788878085670e2; 	cQ.w[115] = 3.8298599230421540e-49;
						cQ.x[116] = 3.9959862252409562e2; 	cQ.w[116] = 8.6830434227075531e-52;
						cQ.x[117] = 4.1432274380983757e2; 	cQ.w[117] = 7.8569306274869450e-66;
						cQ.x[118] = 4.3124084763817151e2; 	cQ.w[118] = 1.0876766815421296e-79;
						cQ.x[119] = 4.5227326181239169e2; 	cQ.w[119] = 1.7576416430685541e-80;
						break;
					case 128:
						cQ.x[0] = 4.8097546269833893e-3; 	cQ.w[0] = 2.7607937255104182e-1;
						cQ.x[1] = 4.3288873981462912e-2; 	cQ.w[1] = 2.6566783614880718e-1;
						cQ.x[2] = 1.2025288615546319e-1; 	cQ.w[2] = 2.4600646227023164e-1;
						cQ.x[3] = 2.3571334284486281e-1; 	cQ.w[3] = 2.1920556026860463e-1;
						cQ.x[4] = 3.8968758351746751e-1; 	cQ.w[4] = 1.8795185561410056e-1;
						cQ.x[5] = 5.8219874974842322e-1; 	cQ.w[5] = 1.5506780698139096e-1;
						cQ.x[6] = 8.1327580437858389e-1; 	cQ.w[6] = 1.2310171404959776e-1;
						cQ.x[7] = 1.0829535555341720e0; 	cQ.w[7] = 9.4028368962612467e-2;
						cQ.x[8] = 1.3912726855559166e0; 	cQ.w[8] = 6.9101651733143485e-2;
						cQ.x[9] = 1.7382797848958800e0; 	cQ.w[9] = 4.8857640815463552e-2;
						cQ.x[10] = 2.1240273910504205e0; 	cQ.w[10] = 3.3232907130961684e-2;
						cQ.x[11] = 2.5485740326082378e0; 	cQ.w[11] = 2.1745560021797622e-2;
						cQ.x[12] = 3.0119842785032310e0; 	cQ.w[12] = 1.3687095820753643e-2;
						cQ.x[13] = 3.5143287925730293e0; 	cQ.w[13] = 8.2862847528073521e-3;
						cQ.x[14] = 4.0556843935355590e0; 	cQ.w[14] = 4.8248404263837572e-3;
						cQ.x[15] = 4.6361341205079486e0; 	cQ.w[15] = 2.7017467389205526e-3;
						cQ.x[16] = 5.2557673042044822e0; 	cQ.w[16] = 1.4548097345670083e-3;
						cQ.x[17] = 5.9146796439632562e0; 	cQ.w[17] = 7.5322752691836494e-4;
						cQ.x[18] = 6.6129732907647177e0; 	cQ.w[18] = 3.7493875811399589e-4;
						cQ.x[19] = 7.3507569364194338e0; 	cQ.w[19] = 1.7941616616270753e-4;
						cQ.x[20] = 8.1281459091173179e0; 	cQ.w[20] = 8.2523935092697431e-5;
						cQ.x[21] = 8.9452622755461916e0; 	cQ.w[21] = 3.6480646490280587e-5;
						cQ.x[22] = 9.8022349498040578e0; 	cQ.w[22] = 1.5497209992155702e-5;
						cQ.x[23] = 1.0699199809346890e1; 	cQ.w[23] = 6.3254880868062279e-6;
						cQ.x[24] = 1.1636299818232177e1; 	cQ.w[24] = 2.4804025381841617e-6;
						cQ.x[25] = 1.2613685157937996e1; 	cQ.w[25] = 9.3427092817735856e-7;
						cQ.x[26] = 1.3631513366058145e1; 	cQ.w[26] = 3.3796972810877498e-7;
						cQ.x[27] = 1.4689949483195887e1; 	cQ.w[27] = 1.1739931842906987e-7;
						cQ.x[28] = 1.5789166208402366e1; 	cQ.w[28] = 3.9152666460117425e-8;
						cQ.x[29] = 1.6929344063530748e1; 	cQ.w[29] = 1.2533919800530000e-8;
						cQ.x[30] = 1.8110671566903898e1; 	cQ.w[30] = 3.8508855814103435e-9;
						cQ.x[31] = 1.9333345416721943e1; 	cQ.w[31] = 1.1352650067802054e-9;
						cQ.x[32] = 2.0597570684666645e1; 	cQ.w[32] = 3.2107610902767967e-10;
						cQ.x[33] = 2.1903561020192286e1; 	cQ.w[33] = 8.7096417230648332e-11;
						cQ.x[34] = 2.3251538866027899e1; 	cQ.w[34] = 2.2655716245890225e-11;
						cQ.x[35] = 2.4641735685453456e1; 	cQ.w[35] = 5.6498933794313194e-12;
						cQ.x[36] = 2.6074392201953201e1; 	cQ.w[36] = 1.3504651749637866e-12;
						cQ.x[37] = 2.7549758651893053e1; 	cQ.w[37] = 3.0931346162094001e-13;
						cQ.x[38] = 2.9068095050916087e1; 	cQ.w[38] = 6.7869387110402112e-14;
						cQ.x[39] = 3.0629671474800940e1; 	cQ.w[39] = 1.4262371047812445e-14;
						cQ.x[40] = 3.2234768355582876e1; 	cQ.w[40] = 2.8696621621608344e-15;
						cQ.x[41] = 3.3883676793796597e1; 	cQ.w[41] = 5.5266882811584044e-16;
						cQ.x[42] = 3.5576698887764130e1; 	cQ.w[42] = 1.0185057527135605e-16;
						cQ.x[43] = 3.7314148080920709e1; 	cQ.w[43] = 1.7955213298120343e-17;
						cQ.x[44] = 3.9096349528247114e1; 	cQ.w[44] = 3.0269510335107605e-18;
						cQ.x[45] = 4.0923640482958862e1; 	cQ.w[45] = 4.8782259058116175e-19;
						cQ.x[46] = 4.2796370704691817e1; 	cQ.w[46] = 7.5129177544606858e-20;
						cQ.x[47] = 4.4714902890520744e1; 	cQ.w[47] = 1.1053212066329750e-20;
						cQ.x[48] = 4.6679613130253012e1; 	cQ.w[48] = 1.5528826688231492e-21;
						cQ.x[49] = 4.8690891387554891e1; 	cQ.w[49] = 2.0825248019667178e-22;
						cQ.x[50] = 5.0749142008593753e1; 	cQ.w[50] = 2.6648208774549398e-23;
						cQ.x[51] = 5.2854784260017013e1; 	cQ.w[51] = 3.2523000556521267e-24;
						cQ.x[52] = 5.5008252898239286e1; 	cQ.w[52] = 3.7841616175987121e-25;
						cQ.x[53] = 5.7209998772174173e1; 	cQ.w[53] = 4.1957539145536696e-26;
						cQ.x[54] = 5.9460489461728161e1; 	cQ.w[54] = 4.4310753939163665e-27;
						cQ.x[55] = 6.1760209954572964e1; 	cQ.w[55] = 4.4550960738061500e-28;
						cQ.x[56] = 6.4109663363931410e1; 	cQ.w[56] = 4.2622210868369363e-29;
						cQ.x[57] = 6.6509371690352848e1; 	cQ.w[57] = 3.8781073546455843e-30;
						cQ.x[58] = 6.8959876630719803e1; 	cQ.w[58] = 3.3540852993298416e-31;
						cQ.x[59] = 7.1461740438021006e1; 	cQ.w[59] = 2.7558487347346862e-32;
						cQ.x[60] = 7.4015546835750457e1; 	cQ.w[60] = 2.1498637246741931e-33;
						cQ.x[61] = 7.6621901991151567e1; 	cQ.w[61] = 1.5913963136894115e-34;
						cQ.x[62] = 7.9281435551924119e1; 	cQ.w[62] = 1.1170819746689115e-35;
						cQ.x[63] = 8.1994801751454565e1; 	cQ.w[63] = 7.4310093282819445e-37;
						cQ.x[64] = 8.4762680588122938e1; 	cQ.w[64] = 4.6813625515950693e-38;
						cQ.x[65] = 8.7585779084788683e1; 	cQ.w[65] = 2.7909522621573202e-39;
						cQ.x[66] = 9.0464832635170634e1; 	cQ.w[66] = 1.5735113668960615e-40;
						cQ.x[67] = 9.3400606444521639e1; 	cQ.w[67] = 8.3828847876304648e-42;
						cQ.x[68] = 9.6393897072765977e1; 	cQ.w[68] = 4.2167560748751703e-43;
						cQ.x[69] = 9.9445534089129095e1; 	cQ.w[69] = 2.0010861658834713e-44;
						cQ.x[70] = 1.0255638184825758e2; 	cQ.w[70] = 8.9512017753845699e-46;
						cQ.x[71] = 1.0572734139891817e2; 	cQ.w[71] = 3.7708148208752655e-47;
						cQ.x[72] = 1.0895935253759580e2; 	cQ.w[72] = 1.4945851640195755e-48;
						cQ.x[73] = 1.1225339602070313e2; 	cQ.w[73] = 5.5681807795587594e-50;
						cQ.x[74] = 1.1561049595069244e2; 	cQ.w[74] = 1.9479151853826422e-51;
						cQ.x[75] = 1.1903172235315340e2; 	cQ.w[75] = 6.3918683903792039e-53;
						cQ.x[76] = 1.2251819396402152e2; 	cQ.w[76] = 1.9651813087431218e-54;
						cQ.x[77] = 1.2607108124835155e2; 	cQ.w[77] = 5.6544197387738892e-56;
						cQ.x[78] = 1.2969160967477476e2; 	cQ.w[78] = 1.5207362927895532e-57;
						cQ.x[79] = 1.3338106327281588e2; 	cQ.w[79] = 3.8180775856184579e-59;
						cQ.x[80] = 1.3714078850376025e2; 	cQ.w[80] = 8.9367191262601894e-61;
						cQ.x[81] = 1.4097219847981490e2; 	cQ.w[81] = 1.9473442897375391e-62;
						cQ.x[82] = 1.4487677757099489e2; 	cQ.w[82] = 3.9445422891662265e-64;
						cQ.x[83] = 1.4885608644460317e2; 	cQ.w[83] = 7.4159150213626314e-66;
						cQ.x[84] = 1.5291176758849769e2; 	cQ.w[84] = 1.2919234936181112e-67;
						cQ.x[85] = 1.5704555137672418e2; 	cQ.w[85] = 2.0819221394737772e-69;
						cQ.x[86] = 1.6125926274474058e2; 	cQ.w[86] = 3.0978378698037742e-71;
						cQ.x[87] = 1.6555482855162421e2; 	cQ.w[87] = 4.2480171859192708e-73;
						cQ.x[88] = 1.6993428571864366e2; 	cQ.w[88] = 5.3575541752687071e-75;
						cQ.x[89] = 1.7439979024777800e2; 	cQ.w[89] = 6.2010829537374423e-77;
						cQ.x[90] = 1.7895362724065057e2; 	cQ.w[90] = 6.5720579636670588e-79;
						cQ.x[91] = 1.8359822205850647e2; 	cQ.w[91] = 6.3623825304174234e-81;
						cQ.x[92] = 1.8833615278804584e2; 	cQ.w[92] = 5.6118813553226041e-83;
						cQ.x[93] = 1.9317016420706524e2; 	cQ.w[93] = 4.4976052330512204e-85;
						cQ.x[94] = 1.9810318347914946e2; 	cQ.w[94] = 3.2656789300649514e-87;
						cQ.x[95] = 2.0313833784961355e2; 	cQ.w[95] = 2.1415773771514337e-89;
						cQ.x[96] = 2.0827897466747490e2; 	cQ.w[96] = 1.2642009038210272e-91;
						cQ.x[97] = 2.1352868412296766e2; 	cQ.w[97] = 6.6937523730388343e-94;
						cQ.x[98] = 2.1889132517029559e2; 	cQ.w[98] = 3.1668515385882014e-96;
						cQ.x[99] = 2.2437105520529297e2; 	cQ.w[99] = 1.3331988108895188e-98;
						cQ.x[100] = 2.2997236419317798e2; 	cQ.w[100] = 4.9720167074631327e-101;
						cQ.x[101] = 2.3570011410033001e2; 	cQ.w[101] = 1.6347801791277616e-103;
						cQ.x[102] = 2.4155958468639006e2; 	cQ.w[102] = 4.7134964221305137e-106;
						cQ.x[103] = 2.4755652697313934e2; 	cQ.w[103] = 1.1851002424541479e-108;
						cQ.x[104] = 2.5369722604409365e2; 	cQ.w[104] = 2.5820543464647026e-111;
						cQ.x[105] = 2.5998857527081546e2; 	cQ.w[105] = 4.8417253805291566e-114;
						cQ.x[106] = 2.6643816464709542e2; 	cQ.w[106] = 7.7550688581883223e-117;
						cQ.x[107] = 2.7305438669552003e2; 	cQ.w[107] = 1.0522110136923271e-119;
						cQ.x[108] = 2.7984656447262174e2; 	cQ.w[108] = 1.1982172530332628e-122;
						cQ.x[109] = 2.8682510765704206e2; 	cQ.w[109] = 1.1333834975480497e-125;
						cQ.x[110] = 2.9400170473751538e2; 	cQ.w[110] = 8.8091250671314200e-129;
						cQ.x[111] = 3.0138956219582173e2; 	cQ.w[111] = 5.1031266854191002e-132;
						cQ.x[112] = 3.0900370572897370e2; 	cQ.w[112] = 3.7889105113753785e-134;
						cQ.x[113] = 3.1686136465414833e2; 	cQ.w[113] = 1.6614242788060531e-135;
						cQ.x[114] = 3.2498246980378542e2; 	cQ.w[114] = 3.3159273335592529e-133;
						cQ.x[115] = 3.3339030932831591e2; 	cQ.w[115] = 6.6730040001944718e-134;
						cQ.x[116] = 3.4211240916012430e2; 	cQ.w[116] = 1.2949575356743615e-133;
						cQ.x[117] = 3.5118174138517492e2; 	cQ.w[117] = 9.1532484444758395e-133;
						cQ.x[118] = 3.6063842559961562e2; 	cQ.w[118] = 1.9549729552902828e-134;
						cQ.x[119] = 3.7053219762567725e2; 	cQ.w[119] = 8.0781662197334676e-122;
						cQ.x[120] = 3.8092612308026272e2; 	cQ.w[120] = 3.7898513479231859e-137;
						cQ.x[121] = 3.9190243416377710e2; 	cQ.w[121] = 1.0948592196247296e-132;
						cQ.x[122] = 4.0357221976913419e2; 	cQ.w[122] = 3.6495237116054534e-50;
						cQ.x[123] = 4.1609268503849895e2; 	cQ.w[123] = 7.9553417261007054e-54;
						cQ.x[124] = 4.2970092632479113e2; 	cQ.w[124] = 5.8160707456163802e-67;
						cQ.x[125] = 4.4478945491488902e2; 	cQ.w[125] = 1.6308626485034593e-79;
						cQ.x[126] = 4.6211398112631172e2; 	cQ.w[126] = 8.1797036981379553e-83;
						cQ.x[127] = 4.8363457770593729e2; 	cQ.w[127] = 4.4300469125827045e-82;
						break;
				}
			}
 
			/******************************************************************************/
			//    Input, int KIND, the rule.
			//    1, Legendre,             (a,b)       1.0
			//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
			//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
			//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
			//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
			//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
			//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
			//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

			double r8_epsilon()
			//  Purpose:
			//
			//    R8_EPSILON returns the R8 roundoff unit.
			//
			//  Author:
			//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
			//    C++ version by John Burkardt.
			//
			//  Licensing:
			//
			//    This code is distributed under the GNU LGPL license. 
			//  Reference:
			//
			//  Discussion:
			//
			//    The roundoff unit is a number R which is a power of 2 with the
			//    property that, to the precision of the computer's arithmetic,
			//      1 < 1 + R
			//    but
			//      1 = ( 1 + R / 2 )
			//
			//  Parameters:
			//
			//    Output, double R8_EPSILON, the R8 round-off unit.
			{
				const double value = 2.220446049250313E-016;

				return value;
			}

			double r8_huge()
			//  Purpose:
			//
			//    R8_HUGE returns a "huge" R8.
			//
			//  Author:
			//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
			//    C++ version by John Burkardt.
			//
			//  Licensing:
			//
			//    This code is distributed under the GNU LGPL license. 
			//
			//  Discussion:
			//
			//    The value returned by this function is NOT required to be the
			//    maximum representable R8.  This value varies from machine to machine,
			//    from compiler to compiler, and may cause problems when being printed.
			//    We simply want a "very large" but non-infinite number.
			//
			//  Parameters:
			//
			//    Output, double R8_HUGE, a "huge" R8 value.
			{
				double value;

				value = 1.0E+30;

				return value;
			}

			double r8_sign(double x)
			//  Purpose:
			//
			//    R8_SIGN returns the sign of an R8.
			//
			//  Author:
			//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
			//    C++ version by John Burkardt.
			//
			//  Licensing:
			//
			//    This code is distributed under the GNU LGPL license. 
			//
			//  Parameters:
			//
			//    Input, double X, the number whose sign is desired.
			//    Output, double R8_SIGN, the sign of X.
			{
				double value;

				if (x < 0.0)
				{
					value = -1.0;
				}
				else
				{
					value = 1.0;
				}
				return value;
			}

			void cgqf(int nt, int kind, double alpha, double beta, double a, double b,
				double t[], double wts[])
			//  Purpose:
			//
			//    CGQF computes knots and weights of a Gauss quadrature formula.
			//
			//  Author:
			//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
			//    C++ version by John Burkardt.
			//
			//  Licensing:
			//
			//    This code is distributed under the GNU LGPL license. 
			//
			//  Discussion:
			//
			//    The user may specify the interval (A,B).
			//    Only simple knots are produced.
			//    Use routine EIQFS to evaluate this quadrature formula.
			//
			//  Reference:
			//
			//    Sylvan Elhay, Jaroslav Kautsky,
			//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
			//    Interpolatory Quadrature,
			//    ACM Transactions on Mathematical Software,
			//    Volume 13, Number 4, December 1987, pages 399-415.
			//
			//  Parameters:
			//
			//    Input, int NT, the number of knots.
			//
			//    Input, int KIND, the rule.
			//    1, Legendre,             (a,b)       1.0
			//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
			//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
			//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
			//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
			//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
			//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
			//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
			//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
			//
			//    Input, double ALPHA, the value of Alpha, if needed.
			//    Input, double BETA, the value of Beta, if needed.
			//    Input, double A, B, the interval endpoints, or other parameters.
			//
			//    Output, double T[NT], the knots.
			//    Output, double WTS[NT], the weights.
			{
				int i;
				int *mlt;
				int *ndx;
				//  Compute the Gauss quadrature formula for default values of A and B.
				cdgqf(nt, kind, alpha, beta, t, wts);

				//  Prepare to scale the quadrature formula to other weight function with valid A and B.
				mlt = new int[nt];
				for (i = 0; i < nt; i++)
				{
					mlt[i] = 1;
				}
				ndx = new int[nt];
				for (i = 0; i < nt; i++)
				{
					ndx[i] = i + 1;
				}
				scqf(nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b);

				delete[] mlt;
				delete[] ndx;

				return;
			}

			void cdgqf(int nt, int kind, double alpha, double beta, double t[],
				double wts[])
			//  Purpose:
			//
			//    CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
			//
			//  Author:
			//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
			//    C++ version by John Burkardt.
			//
			//  Licensing:
			//
			//    This code is distributed under the GNU LGPL license. 
			//
			//  Discussion:
			//
			//    This routine computes all the knots and weights of a Gauss quadrature
			//    formula with a classical weight function with default values for A and B,
			//    and only simple knots.
			//    There are no moments checks and no printing is done.
			//    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
			//
			//  Reference:
			//
			//    Sylvan Elhay, Jaroslav Kautsky,
			//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
			//    Interpolatory Quadrature,
			//    ACM Transactions on Mathematical Software,
			//    Volume 13, Number 4, December 1987, pages 399-415.
			//
			//  Parameters:
			//
			//    Input, int NT, the number of knots.
			//
			//    Input, int KIND, the rule.
			//    1, Legendre,             (a,b)       1.0
			//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
			//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
			//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
			//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
			//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
			//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
			//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
			//
			//    Input, double ALPHA, the value of Alpha, if needed.
			//    Input, double BETA, the value of Beta, if needed.
			//    Output, double T[NT], the knots.
			//    Output, double WTS[NT], the weights.
			{
				double *aj;
				double *bj;
				double zemu;

				parchk(kind, 2 * nt, alpha, beta);
				//  Get the Jacobi matrix and zero-th moment.
				aj = new double[nt];
				bj = new double[nt];

				zemu = class_matrix(kind, nt, alpha, beta, aj, bj);
				//  Compute the knots and weights.
				sgqf(nt, aj, bj, zemu, t, wts);

				delete[] aj;
				delete[] bj;

				return;
			}

			double class_matrix(int kind, int m, double alpha, double beta, double aj[],
				double bj[])
			//  Purpose:
			//
			//    CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
			//
			//  Author:
			//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
			//    C++ version by John Burkardt.
			//
			//  Licensing:
			//
			//    This code is distributed under the GNU LGPL license. 
			//
			//  Discussion:
			//
			//    This routine computes the diagonal AJ and sub-diagonal BJ
			//    elements of the order M tridiagonal symmetric Jacobi matrix
			//    associated with the polynomials orthogonal with respect to
			//    the weight function specified by KIND.
			//
			//    For weight functions 1-7, M elements are defined in BJ even
			//    though only M-1 are needed.  For weight function 8, BJ(M) is
			//    set to zero.
			//
			//    The zero-th moment of the weight function is returned in ZEMU.
			//
			//  Reference:
			//
			//    Sylvan Elhay, Jaroslav Kautsky,
			//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
			//    Interpolatory Quadrature,
			//    ACM Transactions on Mathematical Software,
			//    Volume 13, Number 4, December 1987, pages 399-415.
			//
			//  Parameters:
			//
			//    Input, int KIND, the rule.
			//    1, Legendre,             (a,b)       1.0
			//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
			//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
			//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
			//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
			//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
			//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
			//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
			//
			//    Input, int M, the order of the Jacobi matrix.
			//    Input, double ALPHA, the value of Alpha, if needed.
			//    Input, double BETA, the value of Beta, if needed.
			//    Output, double AJ[M], BJ[M], the diagonal and subdiagonal
			//    of the Jacobi matrix.
			//    Output, double CLASS_MATRIX, the zero-th moment.
			{
				double a2b2;
				double ab;
				double aba;
				double abi;
				double abj;
				double abti;
				double apone;
				int i;
				double pi = 3.14159265358979323846264338327950;
				double temp;
				double temp2;
				double zemu;

				temp = r8_epsilon();

				parchk(kind, 2 * m - 1, alpha, beta);

				temp2 = 0.5;

				if (500.0 * temp < fabs(pow(tgamma(temp2), 2) - pi))
				{
					//cout << "\n";
					//cout << "CLASS_MATRIX - Fatal error!\n";
					//cout << "  Gamma function does not match machine parameters.\n";
					//exit(1);
				}

				if (kind == 1)
				{
					ab = 0.0;

					zemu = 2.0 / (ab + 1.0);

					for (i = 0; i < m; i++)
					{
						aj[i] = 0.0;
					}

					for (i = 1; i <= m; i++)
					{
						abi = i + ab * (i % 2);
						abj = 2 * i + ab;
						bj[i - 1] = sqrt(abi * abi / (abj * abj - 1.0));
					}
				}
				else if (kind == 2)
				{
					zemu = pi;

					for (i = 0; i < m; i++)
					{
						aj[i] = 0.0;
					}

					bj[0] = sqrt(0.5);
					for (i = 1; i < m; i++)
					{
						bj[i] = 0.5;
					}
				}
				else if (kind == 3)
				{
					ab = alpha * 2.0;
					zemu = pow(2.0, ab + 1.0) * pow(tgamma(alpha + 1.0), 2)
						/ tgamma(ab + 2.0);

					for (i = 0; i < m; i++)
					{
						aj[i] = 0.0;
					}

					bj[0] = sqrt(1.0 / (2.0 * alpha + 3.0));
					for (i = 2; i <= m; i++)
					{
						bj[i - 1] = sqrt(i * (i + ab) / (4.0 * pow(i + alpha, 2) - 1.0));
					}
				}
				else if (kind == 4)
				{
					ab = alpha + beta;
					abi = 2.0 + ab;
					zemu = pow(2.0, ab + 1.0) * tgamma(alpha + 1.0)
						* tgamma(beta + 1.0) / tgamma(abi);
					aj[0] = (beta - alpha) / abi;
					bj[0] = sqrt(4.0 * (1.0 + alpha) * (1.0 + beta)
						/ ((abi + 1.0) * abi * abi));
					a2b2 = beta * beta - alpha * alpha;

					for (i = 2; i <= m; i++)
					{
						abi = 2.0 * i + ab;
						aj[i - 1] = a2b2 / ((abi - 2.0) * abi);
						abi = abi * abi;
						bj[i - 1] = sqrt(4.0 * i * (i + alpha) * (i + beta) * (i + ab)
							/ ((abi - 1.0) * abi));
					}
				}
				else if (kind == 5)
				{
					zemu = tgamma(alpha + 1.0);

					for (i = 1; i <= m; i++)
					{
						aj[i - 1] = 2.0 * i - 1.0 + alpha;
						bj[i - 1] = sqrt(i * (i + alpha));
					}
				}
				else if (kind == 6)
				{
					zemu = tgamma((alpha + 1.0) / 2.0);

					for (i = 0; i < m; i++)
					{
						aj[i] = 0.0;
					}

					for (i = 1; i <= m; i++)
					{
						bj[i - 1] = sqrt((i + alpha * (i % 2)) / 2.0);
					}
				}
				else if (kind == 7)
				{
					ab = alpha;
					zemu = 2.0 / (ab + 1.0);

					for (i = 0; i < m; i++)
					{
						aj[i] = 0.0;
					}

					for (i = 1; i <= m; i++)
					{
						abi = i + ab * (i % 2);
						abj = 2 * i + ab;
						bj[i - 1] = sqrt(abi * abi / (abj * abj - 1.0));
					}
				}
				else if (kind == 8)
				{
					ab = alpha + beta;
					zemu = tgamma(alpha + 1.0) * tgamma(-(ab + 1.0))
						/ tgamma(-beta);
					apone = alpha + 1.0;
					aba = ab * apone;
					aj[0] = -apone / (ab + 2.0);
					bj[0] = -aj[0] * (beta + 1.0) / (ab + 2.0) / (ab + 3.0);
					for (i = 2; i <= m; i++)
					{
						abti = ab + 2.0 * i;
						aj[i - 1] = aba + 2.0 * (ab + i) * (i - 1);
						aj[i - 1] = -aj[i - 1] / abti / (abti - 2.0);
					}

					for (i = 2; i <= m - 1; i++)
					{
						abti = ab + 2.0 * i;
						bj[i - 1] = i * (alpha + i) / (abti - 1.0) * (beta + i)
							/ (abti * abti) * (ab + i) / (abti + 1.0);
					}
					bj[m - 1] = 0.0;
					for (i = 0; i < m; i++)
					{
						bj[i] = sqrt(bj[i]);
					}
				}

				return zemu;
			}

			void imtqlx(int n, double d[], double e[], double z[])
			//  Purpose:
			//
			//    IMTQLX diagonalizes a symmetric tridiagonal matrix.
			//
			//  Author:
			//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
			//    C++ version by John Burkardt.
			//
			//  Licensing:
			//
			//    This code is distributed under the GNU LGPL license. 
			//
			//  Discussion:
			//
			//    This routine is a slightly modified version of the EISPACK routine to 
			//    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
			//
			//    The authors thank the authors of EISPACK for permission to use this
			//    routine. 
			//
			//    It has been modified to produce the product Q' * Z, where Z is an input 
			//    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
			//    The changes consist (essentialy) of applying the orthogonal transformations
			//    directly to Z as they are generated.
			//
			//  Reference:
			//
			//    Sylvan Elhay, Jaroslav Kautsky,
			//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
			//    Interpolatory Quadrature,
			//    ACM Transactions on Mathematical Software,
			//    Volume 13, Number 4, December 1987, pages 399-415.
			//
			//    Roger Martin, James Wilkinson,
			//    The Implicit QL Algorithm,
			//    Numerische Mathematik,
			//    Volume 12, Number 5, December 1968, pages 377-383.
			//
			//  Parameters:
			//
			//    Input, int N, the order of the matrix.
			//
			//    Input/output, double D(N), the diagonal entries of the matrix.
			//    On output, the information in D has been overwritten.
			//
			//    Input/output, double E(N), the subdiagonal entries of the 
			//    matrix, in entries E(1) through E(N-1).  On output, the information in
			//    E has been overwritten.
			//
			//    Input/output, double Z(N).  On input, a vector.  On output,
			//    the value of Q' * Z, where Q is the matrix that diagonalizes the
			//    input symmetric tridiagonal matrix.
			{
				double b;
				double c;
				double f;
				double g;
				int i;
				int ii;
				int itn = 30;
				int j;
				int k;
				int l;
				int m;
				int mml;
				double p;
				double prec;
				double r;
				double s;

				prec = r8_epsilon();

				if (n == 1)
				{
					return;
				}

				e[n - 1] = 0.0;

				for (l = 1; l <= n; l++)
				{
					j = 0;
					for (; ; )
					{
						for (m = l; m <= n; m++)
						{
							if (m == n)
							{
								break;
							}

							if (fabs(e[m - 1]) <= prec * (fabs(d[m - 1]) + fabs(d[m])))
							{
								break;
							}
						}
						p = d[l - 1];
						if (m == l)
						{
							break;
						}
						if (itn <= j)
						{
							//cout << "\n";
							//cout << "IMTQLX - Fatal error!\n";
							//cout << "  Iteration limit exceeded\n";
							//exit(1);
						}

						j = j + 1;
						g = (d[l] - p) / (2.0 * e[l - 1]);
						r = sqrt(g * g + 1.0);
						g = d[m - 1] - p + e[l - 1] / (g + fabs(r) * r8_sign(g));
						s = 1.0;
						c = 1.0;
						p = 0.0;
						mml = m - l;

						for (ii = 1; ii <= mml; ii++)
						{
							i = m - ii;
							f = s * e[i - 1];
							b = c * e[i - 1];

							if (fabs(g) <= fabs(f))
							{
								c = g / f;
								r = sqrt(c * c + 1.0);
								e[i] = f * r;
								s = 1.0 / r;
								c = c * s;
							}
							else
							{
								s = f / g;
								r = sqrt(s * s + 1.0);
								e[i] = g * r;
								c = 1.0 / r;
								s = s * c;
							}
							g = d[i] - p;
							r = (d[i - 1] - g) * s + 2.0 * c * b;
							p = s * r;
							d[i] = g + p;
							g = c * r - b;
							f = z[i];
							z[i] = s * z[i - 1] + c * f;
							z[i - 1] = c * z[i - 1] - s * f;
						}
						d[l - 1] = d[l - 1] - p;
						e[l - 1] = g;
						e[m - 1] = 0.0;
					}
				}
				//
				//  Sorting.
				//
				for (ii = 2; ii <= m; ii++)
				{
					i = ii - 1;
					k = i;
					p = d[i - 1];

					for (j = ii; j <= n; j++)
					{
						if (d[j - 1] < p)
						{
							k = j;
							p = d[j - 1];
						}
					}

					if (k != i)
					{
						d[k - 1] = d[i - 1];
						d[i - 1] = p;
						p = z[i - 1];
						z[i - 1] = z[k - 1];
						z[k - 1] = p;
					}
				}
				return;
			}

			void parchk(int kind, int m, double alpha, double beta)
			//  Purpose:
			//
			//    PARCHK checks parameters ALPHA and BETA for classical weight functions. 
			//
			//  Author:
			//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
			//    C++ version by John Burkardt.
			//
			//  Licensing:
			//
			//    This code is distributed under the GNU LGPL license. 
			//
			//  Reference:
			//
			//    Sylvan Elhay, Jaroslav Kautsky,
			//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
			//    Interpolatory Quadrature,
			//    ACM Transactions on Mathematical Software,
			//    Volume 13, Number 4, December 1987, pages 399-415.
			//
			//  Parameters:
			//
			//    Input, int KIND, the rule.
			//    1, Legendre,             (a,b)       1.0
			//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
			//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
			//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
			//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
			//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
			//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
			//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
			//
			//    Input, int M, the order of the highest moment to
			//    be calculated.  This value is only needed when KIND = 8.
			//
			//    Input, double ALPHA, BETA, the parameters, if required
			//    by the value of KIND.
			{
				double tmp;

				if (kind <= 0)
				{
					//cout << "\n";
					//cout << "PARCHK - Fatal error!\n";
					//cout << "  KIND <= 0.\n";
					//exit(1);
				}

				//  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
				if (3 <= kind && alpha <= -1.0)
				{
					//cout << "\n";
					//cout << "PARCHK - Fatal error!\n";
					//cout << "  3 <= KIND and ALPHA <= -1.\n";
					//exit(1);
				}

				//  Check BETA for Jacobi.
				if (kind == 4 && beta <= -1.0)
				{
					//cout << "\n";
					//cout << "PARCHK - Fatal error!\n";
					//cout << "  KIND == 4 and BETA <= -1.0.\n";
					//exit(1);
				}

				//  Check ALPHA and BETA for rational.
				if (kind == 8)
				{
					tmp = alpha + beta + m + 1.0;
					if (0.0 <= tmp || tmp <= beta)
					{
						//cout << "\n";
						//cout << "PARCHK - Fatal error!\n";
						//cout << "  KIND == 8 but condition on ALPHA and BETA fails.\n";
						//exit(1);
					}
				}
				return;
			}

			void scqf(int nt, double t[], int mlt[], double wts[], int nwts, int ndx[],
				double swts[], double st[], int kind, double alpha, double beta, double a, double b)
			//  Purpose:
			//
			//    SCQF scales a quadrature formula to a nonstandard interval.
			//
			//  Author:
			//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
			//    C++ version by John Burkardt.
			//
			//  Licensing:
			//
			//    This code is distributed under the GNU LGPL license. 
			//
			//  Discussion:
			//
			//    The arrays WTS and SWTS may coincide.
			//    The arrays T and ST may coincide.
			//
			//  Reference:
			//
			//    Sylvan Elhay, Jaroslav Kautsky,
			//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
			//    Interpolatory Quadrature,
			//    ACM Transactions on Mathematical Software,
			//    Volume 13, Number 4, December 1987, pages 399-415.
			//
			//  Parameters:
			//
			//    Input, int NT, the number of knots.
			//    Input, double T[NT], the original knots.
			//    Input, int MLT[NT], the multiplicity of the knots.
			//    Input, double WTS[NWTS], the weights.
			//    Input, int NWTS, the number of weights.
			//    Input, int NDX[NT], used to index the array WTS.  
			//    For more details see the comments in CAWIQ.
			//
			//    Output, double SWTS[NWTS], the scaled weights.
			//    Output, double ST[NT], the scaled knots.
			//
			//    Input, int KIND, the rule.
			//    1, Legendre,             (a,b)       1.0
			//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
			//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
			//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
			//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
			//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
			//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
			//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
			//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
			//
			//    Input, double ALPHA, the value of Alpha, if needed.
			//    Input, double BETA, the value of Beta, if needed.
			//    Input, double A, B, the interval endpoints.
			{
				double al;
				double be;
				int i;
				int k;
				int l;
				double p;
				double shft;
				double slp;
				double temp;
				double tmp;

				temp = r8_epsilon();

				parchk(kind, 1, alpha, beta);

				if (kind == 1)
				{
					al = 0.0;
					be = 0.0;
					if (fabs(b - a) <= temp)
					{
						//cout << "\n";
						//cout << "SCQF - Fatal error!\n";
						//cout << "  |B - A| too small.\n";
						//exit(1);
					}
					shft = (a + b) / 2.0;
					slp = (b - a) / 2.0;
				}
				else if (kind == 2)
				{
					al = -0.5;
					be = -0.5;
					if (fabs(b - a) <= temp)
					{
						//cout << "\n";
						//cout << "SCQF - Fatal error!\n";
						//cout << "  |B - A| too small.\n";
						//exit(1);
					}
					shft = (a + b) / 2.0;
					slp = (b - a) / 2.0;
				}
				else if (kind == 3)
				{
					al = alpha;
					be = alpha;
					if (fabs(b - a) <= temp)
					{
						//cout << "\n";
						//cout << "SCQF - Fatal error!\n";
						//cout << "  |B - A| too small.\n";
						//exit(1);
					}
					shft = (a + b) / 2.0;
					slp = (b - a) / 2.0;
				}
				else if (kind == 4)
				{
					al = alpha;
					be = beta;

					if (fabs(b - a) <= temp)
					{
						//cout << "\n";
						//cout << "SCQF - Fatal error!\n";
						//cout << "  |B - A| too small.\n";
						//exit(1);
					}
					shft = (a + b) / 2.0;
					slp = (b - a) / 2.0;
				}
				else if (kind == 5)
				{
					if (b <= 0.0)
					{
						//cout << "\n";
						//cout << "SCQF - Fatal error!\n";
						//cout << "  B <= 0\n";
						//exit(1);
					}
					shft = a;
					slp = 1.0 / b;
					al = alpha;
					be = 0.0;
				}
				else if (kind == 6)
				{
					if (b <= 0.0)
					{
						//cout << "\n";
						//cout << "SCQF - Fatal error!\n";
						//cout << "  B <= 0.\n";
						//exit(1);
					}
					shft = a;
					slp = 1.0 / sqrt(b);
					al = alpha;
					be = 0.0;
				}
				else if (kind == 7)
				{
					al = alpha;
					be = 0.0;
					if (fabs(b - a) <= temp)
					{
						//cout << "\n";
						//cout << "SCQF - Fatal error!\n";
						//cout << "  |B - A| too small.\n";
						//exit(1);
					}
					shft = (a + b) / 2.0;
					slp = (b - a) / 2.0;
				}
				else if (kind == 8)
				{
					if (a + b <= 0.0)
					{
						//cout << "\n";
						//cout << "SCQF - Fatal error!\n";
						//cout << "  A + B <= 0.\n";
						//exit(1);
					}
					shft = a;
					slp = a + b;
					al = alpha;
					be = beta;
				}
				else if (kind == 9)
				{
					al = 0.5;
					be = 0.5;
					if (fabs(b - a) <= temp)
					{
						//cout << "\n";
						//cout << "SCQF - Fatal error!\n";
						//cout << "  |B - A| too small.\n";
						//exit(1);
					}
					shft = (a + b) / 2.0;
					slp = (b - a) / 2.0;
				}

				p = pow(slp, al + be + 1.0);

				for (k = 0; k < nt; k++)
				{
					st[k] = shft + slp * t[k];
					l = abs(ndx[k]);

					if (l != 0)
					{
						tmp = p;
						for (i = l - 1; i <= l - 1 + mlt[k] - 1; i++)
						{
							swts[i] = wts[i] * tmp;
							tmp = tmp * slp;
						}
					}
				}
				return;
			}

			void sgqf(int nt, double aj[], double bj[], double zemu, double t[], double wts[])
			//  Purpose:
			//
			//    SGQF computes knots and weights of a Gauss Quadrature formula.
			//
			//  Author:
			//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
			//    C++ version by John Burkardt.
			//
			//  Licensing:
			//
			//    This code is distributed under the GNU LGPL license. 
			//
			//  Discussion:
			//
			//    This routine computes all the knots and weights of a Gauss quadrature
			//    formula with simple knots from the Jacobi matrix and the zero-th
			//    moment of the weight function, using the Golub-Welsch technique.
			//    
			//  Reference:
			//
			//    Sylvan Elhay, Jaroslav Kautsky,
			//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
			//    Interpolatory Quadrature,
			//    ACM Transactions on Mathematical Software,
			//    Volume 13, Number 4, December 1987, pages 399-415.
			//
			//  Parameters:
			//
			//    Input, int NT, the number of knots.
			//
			//    Input, double AJ[NT], the diagonal of the Jacobi matrix.
			//
			//    Input/output, double BJ[NT], the subdiagonal of the Jacobi 
			//    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
			//
			//    Input, double ZEMU, the zero-th moment of the weight function.
			//
			//    Output, double T[NT], the knots.
			//
			//    Output, double WTS[NT], the weights.
			{
				int i;

				//  Exit if the zero-th moment is not positive.
				if (zemu <= 0.0)
				{
					//cout << "\n";
					//cout << "SGQF - Fatal error!\n";
					//cout << "  ZEMU <= 0.\n";
					//exit(1);
				}

				//  Set up vectors for IMTQLX.
				for (i = 0; i < nt; i++)
				{
					t[i] = aj[i];
				}
				wts[0] = sqrt(zemu);
				for (i = 1; i < nt; i++)
				{
					wts[i] = 0.0;
				}

				//  Diagonalize the Jacobi matrix.
				imtqlx(nt, t, bj, wts);

				for (i = 0; i < nt; i++)
				{
					wts[i] = wts[i] * wts[i];
				}

				return;
			}

			Q1<double, e_host> cQ;
	};

} // namespace mt

#endif