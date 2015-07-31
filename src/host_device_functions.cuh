/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef HOST_DEVICE_FUNCTIONS_H
#define HOST_DEVICE_FUNCTIONS_H

#include <type_traits>
#include <algorithm>

#include "math.cuh"
#include "types.hpp"
#include "traits.cuh"

#include <thrust/complex.h>
#include <thrust/swap.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/for_each.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/tuple.h>
#include <thrust/complex.h>

namespace multem
{
	// Input: E_0(keV), Output: lambda (electron wave)
	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	T get_lambda(const T &E_0)
	{
		T emass = 510.99906;
		T hc = 12.3984244;
		T lambda = hc/sqrt(E_0*(2*emass + E_0));
		return lambda;
	}

	// Input: E_0(keV), Output: sigma (Interaction parameter)
	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	T get_sigma(const T &E_0)
	{
		T emass = 510.99906;
		T c_Pi = 3.141592653589793238463;
		T x = (emass + E_0)/(2*emass + E_0);
		T sigma = 2*c_Pi*x/(get_lambda(E_0)*E_0);
		return sigma;
	}

	// Input: E_0(keV), Output: gamma(relativistic factor)
	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	T get_gamma(const T &E_0)
	{
		T emass = 510.99906;
		T gamma = 1 + E_0/emass;
		return gamma;
	}

	// Input: E_0(keV), Output: gamma*lambda/c_Potf
	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	T get_Vr_factor(const T &E_0, const T &theta)
	{
		T c_Potf = 47.877645145863056;
		T fPot = get_gamma(E_0)*get_lambda(E_0)/(c_Potf*cos(theta));
		return fPot;
	}

	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	T get_Scherzer_focus(const T &E_0, const T &Cs3)
	{
		T lambda = get_lambda(E_0);
		return sqrt(Cs3*lambda);
	}

	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	T get_Scherzer_aperture(const T &E_0, const T &Cs3)
	{
		T lambda = get_lambda(E_0);
		return pow(6.0/(Cs3*pow(lambda, 3)), 0.25);
	}

	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	void get_Scherzer_conditions(const T &E_0, const T &Cs3, T &defocus, T &aperture)
	{
		defocus = get_Scherzer_focus(E_0, Cs3);
		aperture = get_Scherzer_aperture(E_0, Cs3);
	}

	template<eDevice dev>
	void device_synchronize()
	{
		if(dev == e_Host)
		{

		}
		else if(dev == e_Device)
		{
			cudaDeviceSynchronize();
		}
	}

	namespace host_device_detail
	{
		template<class TFn, class T>
		inline
		T Root_Finder(TFn fn, T x0, T xe, const T Tol=1e-8, const int itMax=200)
		{
			int it = 0;
			T x, fx, fxe = fn(xe);

			do
			{
				x = 0.5*(x0 + xe);
				fx = fn(x);

				if(fx*fxe<0)			
				{
					x0 = x;
				}
				else
				{
					xe = x;
					fxe = fx;
				}
				it++;
			}while((fabs(fx)>Tol)&&(it < itMax));
 
			return x;
		}

		template<class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void int_fx_x0_xe(TFn fn, const Value_type<TrPP_Coef> &x0, const Value_type<TrPP_Coef> &xe, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rq1, Value_type<TrPP_Coef> &y)
		{
			Value_type<TrPP_Coef> a = 0.5*(xe-x0); 
			Value_type<TrPP_Coef> b = 0.5*(xe+x0);
			Value_type<TrPP_Coef> xi, yi;
			y = 0;

			for(auto i = 0; i< rq1.size; i++)
			{
				xi = a*rq1.x[i] + b;
				fn(xi, rcoef, yi);
				y += a*rq1.w[i]*yi;
			}
		}

		template<class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void int_fx_x0_xe(TFn fn, const Value_type<TrPP_Coef> &x0, const Value_type<TrPP_Coef> &xe, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rq1, Value_type<TrPP_Coef> &y1, Value_type<TrPP_Coef> &y2)
		{
			Value_type<TrPP_Coef> a = 0.5*(xe-x0); 
			Value_type<TrPP_Coef> b = 0.5*(xe+x0);
			Value_type<TrPP_Coef> xi, y1i, y2i;
			y1 = y2 = 0;

			for(auto i = 0; i< rq1.size; i++)
			{
				xi = a*rq1.x[i] + b;
				fn(xi, rcoef, y1i, y2i);
				y1 += a*rq1.w[i]*y1i;
				y2 += a*rq1.w[i]*y2i;
			}
		}

		template<class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void int_fx_x0_pInfty(TFn fn, const Value_type<TrPP_Coef> &x0, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, TrQ1&rq1, Value_type<TrPP_Coef> &y)
		{
			using T = Value_type<TrPP_Coef>;
			T xi, yi;
			y = 0;

			for(auto i = 0; i< rq1.size; i++)
			{
				xi = rq1.x[i] + x0;
				fn(xi, rcoef, yi);
				y += rq1.w[i]*yi;
			}
		}

		template<class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void int_fx_x0_pInfty(TFn fn, const Value_type<TrPP_Coef> &x0, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rq1, Value_type<TrPP_Coef> &y1, Value_type<TrPP_Coef> &y2)
		{
			using T = Value_type<TrPP_Coef>;
			T xi, y1i, y2i;
			y1 = y2 = 0;

			for(auto i = 0; i< rq1.size; i++)
			{
				xi = rq1.x[i] + x0;
				fn(xi, rcoef, y1i, y2i);
				y1 += rq1.w[i]*y1i;
				y2 += rq1.w[i]*y2i;
			}
		}

		template<class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void int_Vz_z0_ze(TFn fn, const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &R, const TrPP_Coef &rcoef, const TrQ1 &rq1, Value_type<TrPP_Coef> &y)
		{
			using T = Value_type<TrPP_Coef>;
			bool split = (z0<0)&&(0<ze);
			T a = (split)?-0.5*z0:0.5*(ze-z0); 
			T b = (split)?0.5*z0:0.5*(ze+z0);
			T zi, ri, yi, R2 = R*R;
			y = 0;

			for(auto i = 0; i< rq1.size; i++)
			{
				zi = a*rq1.x[i] + b;
				ri = sqrt(zi*zi + R2);
				fn(ri, rcoef, yi);
				y += a*rq1.w[i]*yi;
			}

			if(split)
			{
				a = b = 0.5*ze;
				for(auto i = 0; i< rq1.size; i++)
				{
					zi = a*rq1.x[i] + b;
					ri = sqrt(zi*zi + R2);
					fn(ri, rcoef, yi);
					y += a*rq1.w[i]*yi;
				}
			}
		}

		template<class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void int_Vz_dVz_z0_ze(TFn fn, const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &R, const TrPP_Coef &rcoef, const TrQ1 &rq1, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			using T = Value_type<TrPP_Coef>;
			bool split = (z0<0)&&(0<ze);
			T a = (split)?-0.5*z0:0.5*(ze-z0); 
			T b = (split)?0.5*z0:0.5*(ze+z0);
			T zi, ri, yi, dyi, R2 = R*R;
			y = dy = 0;

			for(auto i = 0; i< rq1.size; i++)
			{
				zi = a*rq1.x[i] + b;
				ri = sqrt(zi*zi + R2);
				fn(ri, rcoef, yi, dyi);
				y += a*rq1.w[i]*yi; 
				dy += R*a*rq1.w[i]*dyi/ri;
			}

			if(split)
			{
				a = b = 0.5*ze;
				for(auto i = 0; i< rq1.size; i++)
				{
					zi = a*rq1.x[i] + b;
					ri = sqrt(zi*zi + R2);
					fn(ri, rcoef, yi, dyi);
					y += a*rq1.w[i]*yi; 
					dy += R*a*rq1.w[i]*dyi/ri;
				}
			}
		}

		template<class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void int_VR(TFn fn, const Value_type<TrPP_Coef> &R, const TrPP_Coef &rcoef, const TrQ1 &rq1, Value_type<TrPP_Coef> &y)
		{
			using T = Value_type<TrPP_Coef>;
			T zi, ri, yi, R2 = R*R;
			y = 0;
			for(auto i = 0; i< rq1.size; i++)
			{
				zi = rq1.x[i];
				ri = sqrt(zi*zi + R2);
				fn(ri, rcoef, yi);
				y += rq1.w[i]*yi;
			}
			y *= 2;
		}

		template<class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void int_VR_dVR(TFn fn, const Value_type<TrPP_Coef> &R, const TrPP_Coef &rcoef, const TrQ1 &rq1, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			using T = Value_type<TrPP_Coef>;
			T zi, ri, yi, dyi, R2 = R*R;
			y = dy = 0;
			for(auto i = 0; i< rq1.size; i++)
			{
				zi = rq1.x[i];
				ri = sqrt(zi*zi + R2);
				fn(ri, rcoef, yi, dyi);
				y += rq1.w[i]*yi; 
				dy += R*rq1.w[i]*dyi/ri;
			}
			y *= 2;
			dy *= 2;
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void add_Exponential_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, bool reset=true)
		{
			if(reset)
			{
				y = 0;
			}
			for(auto i = n_0; i< n_e; i++)
			{
				y += rcoef.cl[i]*exp(-rcoef.cnl[i]*x);
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void add_Exponential_dExponential_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy, bool reset=true)
		{
			using T = Value_type<TrPP_Coef>;
			T yt;
			if(reset)
			{
				y = dy = 0;
			}
			for(auto i = n_0; i< n_e; i++)
			{
				y += yt = rcoef.cl[i]*exp(-rcoef.cnl[i]*x);
				dy += -rcoef.cnl[i]*yt;
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void add_Gaussian_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, bool reset=true)
		{
			using T = Value_type<TrPP_Coef>;
			T x2 = x*x;
			if(reset)
			{
				y = 0;
			}
			for(auto i = n_0; i< n_e; i++)
			{
				y += rcoef.cl[i]*exp(-rcoef.cnl[i]*x2);
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void add_Gaussian_dGaussian_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy, bool reset=true)
		{
			using T = Value_type<TrPP_Coef>;
			T yt, x2 = x*x;
			if(reset)
			{
				y = dy = 0;
			}
			for(auto i = n_0; i< n_e; i++)
			{
				y += yt = rcoef.cl[i]*exp(-rcoef.cnl[i]*x2);
				dy += -2*rcoef.cnl[i]*x*yt;
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void add_Lorentzian_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, bool reset=true)
		{
			using T = Value_type<TrPP_Coef>;
			T x2 = x*x;
			if(reset)
			{
				y = 0;
			}
			for(auto i = n_0; i< n_e; i++)
			{
				y += rcoef.cl[i]/(rcoef.cnl[i] + x2);
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void add_Lorentzian_dLorentzian_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy, bool reset=true)
		{
			Value_type<TrPP_Coef> t, yt, x2 = x*x;
			if(reset)
			{
				y = dy = 0;
			}
			for(auto i = n_0; i< n_e; i++)
			{
				t = 1/(rcoef.cnl[i] + x2);
				y += yt = rcoef.cl[i]*t;
				dy += -2*x*yt*t;
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void add_Yukawa_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, bool reset=true)
		{
			using T = Value_type<TrPP_Coef>;
			T ix = 1/x;
			if(reset)
			{
				y = 0;
			}
			for(auto i = n_0; i< n_e; i++)
			{
				y += rcoef.cl[i]*exp(-rcoef.cnl[i]*x)*ix;
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void add_Yukawa_dYukawa_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy, bool reset=true)
		{
			using T = Value_type<TrPP_Coef>;
			T yt, ix = 1/x;
			if(reset)
			{
				y = dy = 0;
			}
			for(auto i = n_0; i< n_e; i++)
			{
				y += yt = rcoef.cl[i]*exp(-rcoef.cnl[i]*x)*ix;
				dy += -(rcoef.cnl[i]+ ix)*yt;
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Pr_Gaussian_feg(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, bool reset=true)
		{
			using T = Value_type<TrPP_Coef>;
			T x2 = x*x;
			if(reset)
			{
				y = 0;
			}
			for(auto i = n_0; i< n_e; i++)
			{
				y += (2*x2 - 3/rcoef.cnl[i])*rcoef.cl[i]*exp(-rcoef.cnl[i]*x2);
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Pr_dPr_Gaussian_feg(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy, bool reset=true)
		{
			using T = Value_type<TrPP_Coef>;
			T yt, x2 = x*x;
			if(reset)
			{
				y = dy = 0;
			}
			for(auto i = n_0; i< n_e; i++)
			{
				yt = rcoef.cl[i]*exp(-rcoef.cnl[i]*x2);
				y += (2*x2 - 3/rcoef.cnl[i])*yt;
				dy += -2*x*(2*rcoef.cnl[i]*x2 - 5)*yt;
			}
		}

		/***************************************************************************/
		/***************************************************************************/

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void feg_Doyle_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gaussian_Fn(x, 0, 4, rcoef, y);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void feg_dfeg_Doyle_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gaussian_dGaussian_Fn(x, 0, 4, rcoef, y, dy);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void feg_Peng_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gaussian_Fn(x, 0, 5, rcoef, y);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void feg_dfeg_Peng_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gaussian_dGaussian_Fn(x, 0, 5, rcoef, y, dy);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void feg_Peng_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gaussian_Fn(x, 0, 5, rcoef, y);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void feg_dfeg_Peng_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gaussian_dGaussian_Fn(x, 0, 5, rcoef, y, dy);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void feg_Kirkland_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Lorentzian_Fn(x, 0, 3, rcoef, y);
			add_Gaussian_Fn(x, 3, 6, rcoef, y, false);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void feg_dfeg_Kirkland_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Lorentzian_dLorentzian_Fn(x, 0, 3, rcoef, y, dy);
			add_Gaussian_dGaussian_Fn(x, 3, 6, rcoef, y, dy, false);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void feg_Weickenmeier_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			using T = Value_type<TrPP_Coef>;
			T x2 = x*x;
	
			y = 0;
			if(x != 0)
			{
				for(auto i = 0; i< 6; i++)
				{
					y += rcoef.cl[i]*(1-exp(-rcoef.cnl[i]*x2));
				}
				y = y/x2;
			}
			else
			{
				for(auto i = 0; i< 6; i++)
				{
					y += rcoef.cl[i]*rcoef.cnl[i];
				}
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void feg_dfeg_Weickenmeier_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			using T = Value_type<TrPP_Coef>;
			T t, x2 = x*x;
			y = dy = 0;
			if(x != 0)
			{
				for(auto i = 0; i< 6; i++)
				{
					t = exp(-rcoef.cnl[i]*x2);
					y += rcoef.cl[i]*(1-t);
					dy += rcoef.cl[i]*(1-(1+rcoef.cnl[i]*x2)*t);
				}
				y = y/x2;
				dy = -2*dy/(x*x2);
			}
			else
			{
				for(auto i = 0; i< 6; i++)
				{
					y += rcoef.cl[i]*rcoef.cnl[i];
				}
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void feg_Lobato_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			using T = Value_type<TrPP_Coef>;
			T t, x2 = x*x;
			y = 0;
			for(auto i = 0; i< 5; i++)
			{
				t = 1/(1 + rcoef.cnl[i]*x2);
				y += rcoef.cl[i]*t*(t + 1);
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void feg_dfeg_Lobato_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			using T = Value_type<TrPP_Coef>;
			T t, x2 = x*x;
			y = dy = 0;
			for(auto i = 0; i< 5; i++)
			{
				t = 1/(1 + rcoef.cnl[i]*x2);
				y += rcoef.cl[i]*t*(t + 1);
				dy += -2*x*rcoef.cl[i]*rcoef.cnl[i]*t*t*(2*t + 1);
			}
		}

		/***************************************************************************/
		/***************************************************************************/

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void fxg_Doyle_0_4(const Value_type<TrPP_Coef> &x, const int &Z, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			feg_Doyle_0_4(x, rcoef, y);
			y = Z - x*x*y;
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void fxg_dfxg_Doyle_0_4(const Value_type<TrPP_Coef> &x, const int &Z, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			feg_dfeg_Doyle_0_4(x, rcoef, y, dy);
			dy = -x*(2*y + x*dy);
			y = Z - x*x*y;
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void fxg_Peng_0_4(const Value_type<TrPP_Coef> &x, const int &Z, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			feg_Peng_0_4(x, rcoef, y);
			y = Z - x*x*y;
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void fxg_dfxg_Peng_0_4(const Value_type<TrPP_Coef> &x, const int &Z, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			feg_dfeg_Peng_0_4(x, rcoef, y, dy);
			dy = -x*(2*y + x*dy);
			y = Z - x*x*y;
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void fxg_Peng_0_12(const Value_type<TrPP_Coef> &x, const int &Z, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			feg_Peng_0_12(x, rcoef, y);
			y = Z - x*x*y;
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void fxg_dfxg_Peng_0_12(const Value_type<TrPP_Coef> &x, const int &Z, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			feg_dfeg_Peng_0_12(x, rcoef, y, dy);
			dy = -x*(2*y + x*dy);
			y = Z - x*x*y;
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void fxg_Kirkland_0_12(const Value_type<TrPP_Coef> &x, const int &Z, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			feg_Kirkland_0_12(x, rcoef, y);
			y = Z - x*x*y;
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void fxg_dfxg_Kirkland_0_12(const Value_type<TrPP_Coef> &x, const int &Z, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			feg_dfeg_Kirkland_0_12(x, rcoef, y, dy);
			dy = -x*(2*y + x*dy);
			y = Z - x*x*y;
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void fxg_Weickenmeier_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gaussian_Fn(x, 0, 6, rcoef, y);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void fxg_dfxg_Weickenmeier_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gaussian_dGaussian_Fn(x, 0, 6, rcoef, y, dy);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void fxg_Lobato_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			using T = Value_type<TrPP_Coef>;
			T t, x2 = x*x;
			y = 0;
			for(auto i = 0; i< 5; i++)
			{
				t = 1/(1+rcoef.cnl[i]*x2);
				y += rcoef.cl[i]*t*t;
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void fxg_dfxg_Lobato_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			using T = Value_type<TrPP_Coef>;
			T t, yt, x2 = x*x;
			y = dy = 0;
			for(auto i = 0; i< 5; i++)
			{
				t = 1/(1+rcoef.cnl[i]*x2);
				y += yt = rcoef.cl[i]*t*t;
				dy += rcoef.cnl[i]*yt*t;
			}
			dy = -4*x*dy;
		}

		/***************************************************************************/
		/***************************************************************************/

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Pr_Doyle_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			Pr_Gaussian_feg(x, 0, 4, rcoef, y);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Pr_dPr_Doyle_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			Pr_dPr_Gaussian_feg(x, 0, 4, rcoef, y, dy);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Pr_Peng_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			Pr_Gaussian_feg(x, 0, 5, rcoef, y);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Pr_dPr_Peng_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			Pr_dPr_Gaussian_feg(x, 0, 5, rcoef, y, dy);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Pr_Peng_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gaussian_Fn(x, 0, 5, rcoef, y);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Pr_dPr_Peng_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			Pr_dPr_Gaussian_feg(x, 0, 5, rcoef, y, dy);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Pr_Kirkland_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Yukawa_Fn(x, 0, 3, rcoef, y);
			Pr_Gaussian_feg(x, 3, 6, rcoef, y, false);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Pr_dPr_Kirkland_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Yukawa_dYukawa_Fn(x, 0, 3, rcoef, y, dy);
			Pr_dPr_Gaussian_feg(x, 3, 6, rcoef, y, dy, false);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Pr_Weickenmeier_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gaussian_Fn(x, 0, 6, rcoef, y);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Pr_dPr_Weickenmeier_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gaussian_dGaussian_Fn(x, 0, 6, rcoef, y, dy);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Pr_Lobato_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Exponential_Fn(x, 0, 5, rcoef, y);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Pr_dPr_Lobato_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Exponential_dExponential_Fn(x, 0, 5, rcoef, y, dy);
		}

		/***************************************************************************/
		/***************************************************************************/

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_Doyle_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gaussian_Fn(x, 0, 4, rcoef, y);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_dVr_Doyle_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gaussian_dGaussian_Fn(x, 0, 4, rcoef, y, dy);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_Peng_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gaussian_Fn(x, 0, 5, rcoef, y);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_dVr_Peng_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gaussian_dGaussian_Fn(x, 0, 5, rcoef, y, dy);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_Peng_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gaussian_Fn(x, 0, 5, rcoef, y);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_dVr_Peng_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gaussian_dGaussian_Fn(x, 0, 5, rcoef, y, dy);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_Kirkland_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Yukawa_Fn(x, 0, 3, rcoef, y);
			add_Gaussian_Fn(x, 3, 6, rcoef, y, false);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_dVr_Kirkland_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Yukawa_dYukawa_Fn(x, 0, 3, rcoef, y, dy);
			add_Gaussian_dGaussian_Fn(x, 3, 6, rcoef, y, dy, false);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_Weickenmeier_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			using T = Value_type<TrPP_Coef>;
			T ix = 1/x;
			y = 0;

			for(auto i = 0; i< 6; i++)
			{
				y += rcoef.cl[i]*erfc(rcoef.cnl[i]*x)*ix;
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_dVr_Weickenmeier_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			using T = Value_type<TrPP_Coef>;
			const T c_Pii2 = 1.772453850905516027298;
			T yt, ix = 1/x, x2 = x*x;
			y = dy = 0;

			for(auto i = 0; i< 6; i++)
			{
				y += yt = rcoef.cl[i]*erfc(rcoef.cnl[i]*x)*ix;
				dy += (-2*rcoef.cl[i]*rcoef.cnl[i]*exp(-rcoef.cnl[i]*rcoef.cnl[i]*x2)/c_Pii2-yt)*ix;
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_Lobato_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			using T = Value_type<TrPP_Coef>;
			T ix = 1/x;
			y = 0;
			for(auto i = 0; i< 5; i++)
			{
				y += rcoef.cl[i]*exp(-rcoef.cnl[i]*x)*ix*(2/rcoef.cnl[i] + x);
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_dVr_Lobato_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			using T = Value_type<TrPP_Coef>;
			T yt, icnl, ix = 1/x;
			y = dy = 0;
			for(auto i = 0; i< 5; i++)
			{
				yt = rcoef.cl[i]*exp(-rcoef.cnl[i]*x)*ix; 
				icnl = 1/rcoef.cnl[i];
				y += yt*(2*icnl + x);
				dy += -yt*(2*icnl*ix + 2 + rcoef.cnl[i]*x);
			}
		}

		/***************************************************************************/
		/***************************************************************************/

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void VR_Doyle_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gaussian_Fn(x, 0, 4, rcoef, y);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void VR_dVR_Doyle_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gaussian_dGaussian_Fn(x, 0, 4, rcoef, y, dy);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void VR_Peng_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gaussian_Fn(x, 0, 5, rcoef, y);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void VR_dVR_Peng_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gaussian_dGaussian_Fn(x, 0, 5, rcoef, y, dy);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void VR_Peng_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			 add_Gaussian_Fn(x, 0, 5, rcoef, y);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void VR_dVR_Peng_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gaussian_dGaussian_Fn(x, 0, 5, rcoef, y, dy);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void VR_Kirkland_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			y = 0;
			for(auto i = 0; i< 3; i++)
			{
				y += rcoef.cl[i]*bessel_k0(rcoef.cnl[i]*x);
			}

			add_Gaussian_Fn(x, 3, 6, rcoef, y, false);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void VR_dVR_Kirkland_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			y = dy = 0;
			for(auto i = 0; i< 3; i++)
			{
				y += rcoef.cl[i]*bessel_k0(rcoef.cnl[i]*x);
				dy += -rcoef.cl[i]*rcoef.cnl[i]*bessel_k1(rcoef.cnl[i]*x);
			}

			add_Gaussian_dGaussian_Fn(x, 3, 6, rcoef, y, dy, false);
		}

		template<class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void VR_Weickenmeier_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y)
		{
			int_VR(Vr_Weickenmeier_0_12<TrPP_Coef>, x, rcoef, rqz, y);
		}

		template<class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void VR_dVR_Weickenmeier_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			int_VR_dVR(Vr_dVr_Weickenmeier_0_12<TrPP_Coef>, x, rcoef, rqz, y, dy);
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void VR_Lobato_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			y = 0;
			for(auto i = 0; i< 5; i++)
			{
				y += rcoef.cl[i]*(2*bessel_k0(rcoef.cnl[i]*x)/rcoef.cnl[i] + x*bessel_k1(rcoef.cnl[i]*x));
			}
		}

		template<class TrPP_Coef>
		DEVICE_CALLABLE FORCEINLINE 
		void VR_dVR_Lobato_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			using T = Value_type<TrPP_Coef>;
			y = dy = 0;
			for(auto i = 0; i< 5; i++)
			{
				T k0 = bessel_k0(rcoef.cnl[i]*x); 
				T k1 = bessel_k1(rcoef.cnl[i]*x);
				y += rcoef.cl[i]*(2*k0/rcoef.cnl[i] + x*k1);
				dy += -rcoef.cl[i]*(rcoef.cnl[i]*x*k0 + 2*k1);
			}
		}

		/***************************************************************************/
		/***************************************************************************/

		template<class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void Vz_Doyle_0_4(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y)
		{
			int_Vz_z0_ze(Vr_Doyle_0_4<TrPP_Coef>, z0, ze, x, rcoef, rqz, y);
		}

		template<class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void Vz_dVz_Doyle_0_4(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			int_Vz_dVz_z0_ze(Vr_dVr_Doyle_0_4<TrPP_Coef>, z0, ze, x, rcoef, rqz, y, dy);
		}

		template<class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void Vz_Peng_0_4(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y)
		{
			int_Vz_z0_ze(Vr_Peng_0_4<TrPP_Coef>, z0, ze, x, rcoef, rqz, y);
		}

		template<class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void Vz_dVz_Peng_0_4(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			int_Vz_dVz_z0_ze(Vr_dVr_Peng_0_4<TrPP_Coef>, z0, ze, x, rcoef, rqz, y, dy);
		}

		template<class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void Vz_Peng_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y)
		{
			int_Vz_z0_ze(Vr_Peng_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y);
		}

		template<class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void Vz_dVz_Peng_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			int_Vz_dVz_z0_ze(Vr_dVr_Peng_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y, dy);
		}

		template<class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void Vz_Kirkland_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y)
		{
			int_Vz_z0_ze(Vr_Kirkland_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y);
		}

		template<class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void Vz_dVz_Kirkland_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			int_Vz_dVz_z0_ze(Vr_dVr_Kirkland_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y, dy);
		}

		template<class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void Vz_Weickenmeier_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y)
		{
			int_Vz_z0_ze(Vr_Weickenmeier_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y);
		}

		template<class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void Vz_dVz_Weickenmeier_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			int_Vz_dVz_z0_ze(Vr_dVr_Weickenmeier_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y, dy);
		}

		template<class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void Vz_Lobato_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y)
		{
			int_Vz_z0_ze(Vr_Lobato_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y);
		}

		template<class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCEINLINE 
		void Vz_dVz_Lobato_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			int_Vz_dVz_z0_ze(Vr_dVr_Lobato_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y, dy);
		}

		/***************************************************************************/
		/***************************************************************************/

		template<class T> 
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_dVrir_Doyle_0_4(const T &r, T *cl, T *cnl, T &Vr, T &dVrir)
		{
			T r2 = r*r;

			T Vr0 = cl[0]*exp(-cnl[0]*r2); 
			T Vr1 = cl[1]*exp(-cnl[1]*r2); 
			T Vr2 = cl[2]*exp(-cnl[2]*r2); 
			T Vr3 = cl[3]*exp(-cnl[3]*r2);

			Vr = Vr0 + Vr1 + Vr2 + Vr3;
			dVrir = -2*(cnl[0]*Vr0 + cnl[1]*Vr1 + cnl[2]*Vr2 + cnl[3]*Vr3);
		}

		template<class T> 
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_dVrir_Peng_0_4_12(const T &r, T *cl, T *cnl, T &Vr, T &dVrir)
		{
			T r2 = r*r;

			T Vr0 = cl[0]*exp(-cnl[0]*r2); 
			T Vr1 = cl[1]*exp(-cnl[1]*r2); 
			T Vr2 = cl[2]*exp(-cnl[2]*r2); 
			T Vr3 = cl[3]*exp(-cnl[3]*r2);
			T Vr4 = cl[4]*exp(-cnl[4]*r2);

			Vr = Vr0 + Vr1 + Vr2 + Vr3 + Vr4;
			dVrir = -2*(cnl[0]*Vr0 + cnl[1]*Vr1 + cnl[2]*Vr2 + cnl[3]*Vr3 + cnl[4]*Vr4);
		}

		template<class T> 
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_dVrir_Kirkland_0_12(const T &r, T *cl, T *cnl, T &Vr, T &dVrir)
		{
			T ir = 1/r; 
			T r2 = r*r;

			T Vr0 = cl[0]*exp(-cnl[0]*r)*ir; 
			T Vr1 = cl[1]*exp(-cnl[1]*r)*ir; 
			T Vr2 = cl[2]*exp(-cnl[2]*r)*ir; 
			T Vr3 = cl[3]*exp(-cnl[3]*r2); 
			T Vr4 = cl[4]*exp(-cnl[4]*r2); 
			T Vr5 = cl[5]*exp(-cnl[5]*r2);

			Vr = Vr0 + Vr1 + Vr2 + Vr3 + Vr4 + Vr5;
			dVrir = -(Vr0*(cnl[0]+ir) + Vr1*(cnl[1]+ir) + Vr2*(cnl[2]+ir) + 2*r*(cnl[3]*Vr3 + cnl[4]*Vr4 + cnl[5]*Vr5))/r;
		}

		template<class T> 
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_dVrir_Weickenmeier_0_12(const T &r, T *cl, T *cnl, T &Vr, T &dVrir)
		{
			T r2 = r*r;
			T c_Pii2 = 1.772453850905516027298;

			T Vr0 = cl[0]*erfc(cnl[0]*r); 
			T Vr1 = cl[1]*erfc(cnl[1]*r); 
			T Vr2 = cl[2]*erfc(cnl[2]*r); 
			T Vr3 = cl[3]*erfc(cnl[3]*r); 
			T Vr4 = cl[4]*erfc(cnl[4]*r); 
			T Vr5 = cl[5]*erfc(cnl[5]*r);

			Vr = (Vr0 + Vr1 + Vr2 + Vr3 + Vr4 + Vr5)/r;
			dVrir = 2*(cl[0]*cnl[0]*exp(-cnl[0]*cnl[0]*r2) + cl[1]*cnl[1]*exp(-cnl[1]*cnl[1]*r2) + cl[2]*cnl[2]*exp(-cnl[2]*cnl[2]*r2)+ 
			cl[3]*cnl[3]*exp(-cnl[3]*cnl[3]*r2) + cl[4]*cnl[4]*exp(-cnl[4]*cnl[4]*r2) + cl[5]*cnl[5]*exp(-cnl[5]*cnl[5]*r2))/c_Pii2;
			dVrir = -(dVrir + Vr)/r2;
		}

		template<class T> 
		DEVICE_CALLABLE FORCEINLINE 
		void Vr_dVrir_Lobato_0_12(const T &r, T *cl, T *cnl, T &Vr, T &dVrir)
		{
			T cnl0r = cnl[0]*r; 
			T cnl1r = cnl[1]*r; 
			T cnl2r = cnl[2]*r; 
			T cnl3r = cnl[3]*r; 
			T cnl4r = cnl[4]*r;

			T Vr0 = cl[0]*exp(-cnl0r); 
			T Vr1 = cl[1]*exp(-cnl1r); 
			T Vr2 = cl[2]*exp(-cnl2r); 
			T Vr3 = cl[3]*exp(-cnl3r); 
			T Vr4 = cl[4]*exp(-cnl4r);

			Vr = Vr0*(2/cnl0r+1) + Vr1*(2/cnl1r+1) + Vr2*(2/cnl2r+1) + Vr3*(2/cnl3r+1)+ Vr4*(2/cnl4r+1);
			dVrir = -(Vr + Vr0*(cnl0r+1) + Vr1*(cnl1r+1) + Vr2*(cnl2r+1) + Vr3*(cnl3r+1)+ Vr4*(cnl4r+1))/(r*r);
		}

		template<class T> 
		DEVICE_CALLABLE FORCEINLINE 
		int unrolledBinarySearch_c_nR(const T &x0, const T *x)
		{
			int i0 = 0, ie = c_nR-1;
			int im = (i0 + ie)>>1; // divide by 2
			if(x0 < x[im]) ie = im; else i0 = im; //64
			im = (i0 + ie)>>1; 	// divide by 2
			if(x0 < x[im]) ie = im; else i0 = im; //32
			im = (i0 + ie)>>1; 	// divide by 2
			if(x0 < x[im]) ie = im; else i0 = im; //16
			im = (i0 + ie)>>1; 	// divide by 2
			if(x0 < x[im]) ie = im; else i0 = im; //8
			im = (i0 + ie)>>1; 	// divide by 2
			if(x0 < x[im]) ie = im; else i0 = im; //4
			im = (i0 + ie)>>1; 	// divide by 2
			if(x0 < x[im]) ie = im; else i0 = im; //2
			im = (i0 + ie)>>1; 	// divide by 2
			if(x0 < x[im]) ie = im; else i0 = im; //1
	
			return i0;
		}

		// Get Local interpolation coefficients
		template<class T>
		DEVICE_CALLABLE FORCEINLINE 
		void cubic_poly_coef(const int &iR, Atom_Vp<T> &atom_Vp)
		{
			T dR2 = 1.0/(atom_Vp.R2[iR+1]-atom_Vp.R2[iR]);
			T V = atom_Vp.c0[iR]; 
			T Vn = atom_Vp.c0[iR+1];
			T dV = atom_Vp.c1[iR]; 
			T dVn = atom_Vp.c1[iR+1];
			T m = (Vn-V)*dR2; 
			T n = dV+dVn;
			atom_Vp.c0[iR] = V-atom_Vp.c0[c_nR-1];
			atom_Vp.c2[iR] = (3.0*m-n-dV)*dR2;
			atom_Vp.c3[iR] = (n-2.0*m)*dR2*dR2;
		}

		// Cubic polynomial evaluation
		template<class T> 
		DEVICE_CALLABLE FORCEINLINE 
		T eval_cubic_poly(int ix, int iy, T R2, const Grid<T> &grid, const Atom_Vp<T> &atom_Vp, int &ixy)
		{
			R2 = max(R2, atom_Vp.R_min2);

			ix -= static_cast<int>(floor(grid.Rx(ix)/grid.lx))*grid.nx;
			iy -= static_cast<int>(floor(grid.Ry(iy)/grid.ly))*grid.ny;

			ix = grid.iRx_shift(ix);
			iy = grid.iRy_shift(iy);
			ixy = grid.ind_col(ix, iy);

			ix = unrolledBinarySearch_c_nR<T>(R2, atom_Vp.R2);

			T dx = R2 - atom_Vp.R2[ix]; 
			T dx2 = dx*dx;
			return atom_Vp.occ*(atom_Vp.c0[ix] + atom_Vp.c1[ix]*dx + atom_Vp.c2[ix]*dx2 + atom_Vp.c3[ix]*dx2*dx);
		}

		template<class TGrid, class TVector>
		DEVICE_CALLABLE FORCEINLINE 
		void fft2_shift(const int &ix, const int &iy, const TGrid &grid, TVector &M_io)
		{
			int ixy = grid.ind_col(ix, iy); 
			int ixy_shift = grid.ind_col(grid.nxh+ix, grid.nyh+iy);
			thrust::swap(M_io[ixy], M_io[ixy_shift]);

			ixy = grid.ind_col(ix, grid.nyh+iy); 
			ixy_shift = grid.ind_col(grid.nxh+ix, iy);
			thrust::swap(M_io[ixy], M_io[ixy_shift]);
		}

		template<class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCEINLINE 
		void bandwidth_limit(const int &ix, const int &iy, const TGrid &grid, const Value_type<TGrid> &g2_min, const Value_type<TGrid> &g2_max, const Value_type<TVector_c> &w, TVector_c &M_io)
		{
			using value_type_r = Value_type<TGrid>;

			int ixy = grid.ind_col(ix, iy); 
			value_type_r g2 = grid.g2_shift(ix, iy);

			if((g2_min <= g2)&&(g2 <= g2_max))
			{
				M_io[ixy] *= w;
			}
			else
			{
 				M_io[ixy] = 0;
			}
		}

		template<class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCEINLINE 
		void phase_multiplication(const int &ix, const int &iy, const TGrid &grid, 
		const TVector_c &exp_x_i, const TVector_c &exp_y_i, TVector_c &psi_i, TVector_c &psi_o)
		{
			int ixy = grid.ind_col(ix, iy);
			psi_o[ixy] = psi_i[ixy]*exp_y_i[ix]*exp_y_i[iy];
		}

		template<class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCEINLINE 
		void propagator_multiplication(const int &ix, const int &iy, const TGrid &grid, 
		const TVector_c &prop_x_i, const TVector_c &prop_y_i, TVector_c &psi_i, TVector_c &psi_o)
		{
			using value_type_r = Value_type<TGrid>;
			using value_type_c = Value_type<TVector_c>;

			int ixy = grid.ind_col(ix, iy);
			value_type_r g2 = grid.g2_shift(ix, iy);

			if((!grid.bwl)||(g2 < grid.gl2_max))
			{
				psi_o[ixy] = static_cast<value_type_c>(grid.inxy)*psi_i[ixy]*prop_x_i[ix]*prop_y_i[iy];
			}
			else
			{
 				psi_o[ixy] = 0;
			}
		}

		template<class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCEINLINE 
		void probe(const int &ix, const int &iy, const TGrid &grid, const Lens<Value_type<TGrid>> &lens, 
		const Value_type<TGrid> &x, const Value_type<TGrid> &y, TVector_c &fPsi_o)
		{
			using value_type_r = Value_type<TGrid>;

			int ixy = grid.ind_col(ix, iy);
			value_type_r gx = grid.gx_shift(ix);
			value_type_r gy = grid.gy_shift(iy);
			value_type_r g2 = gx*gx + gy*gy;

			if((lens.g2_min <= g2)&&(g2 < lens.g2_max))
			{
				value_type_r chi = x*gx + y*gy + g2*(lens.cCs5*g2*g2+lens.cCs3*g2+lens.cf);
				if(nonZero(lens.m)||nonZero(lens.cmfa2)||nonZero(lens.cmfa3))
				{
					value_type_r g = sqrt(g2);
					value_type_r phi = atan2(gy, gx);
					chi += lens.m*phi + lens.cmfa2*g2*sin(2*(phi-lens.afa2)) + lens.cmfa3*g*g2*sin(3*(phi-lens.afa3)); 			
				}	
				fPsi_o[ixy] = thrust::euler(chi); 
			}
			else
			{
 				fPsi_o[ixy] = 0;
			}
		}

		template<class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCEINLINE 
		void apply_CTF(const int &ix, const int &iy, const TGrid &grid, const Lens<Value_type<TGrid>> &lens,
		const Value_type<TGrid> &gxu, const Value_type<TGrid> &gyu, TVector_c &fPsi_i, TVector_c &fPsi_o)
		{
			using value_type_r = Value_type<TGrid>;

			int ixy = grid.ind_col(ix, iy);
			value_type_r gx = grid.gx_shift(ix);
			value_type_r gy = grid.gy_shift(iy);
			value_type_r g2 = gx*gx + gy*gy;

			if((lens.g2_min <= g2)&&(g2 < lens.g2_max))
			{
				g2 = (gx-gxu)*(gx-gxu) + (gy-gyu)*(gy-gyu);
				value_type_r chi = g2*(lens.cCs5*g2*g2+lens.cCs3*g2+lens.cf);
				if(nonZero(lens.cmfa2)||nonZero(lens.cmfa3))
				{
					value_type_r g = sqrt(g2);
					value_type_r phi = atan2(gy, gx);
					chi += lens.cmfa2*g2*sin(2*(phi-lens.afa2)) + lens.cmfa3*g*g2*sin(3*(phi-lens.afa3)); 			
				}
				fPsi_o[ixy] = fPsi_i[ixy]*thrust::euler(chi);
			}
			else
			{
 				fPsi_o[ixy] = 0;
			}
		}

		template<class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCEINLINE 
		void apply_PCTF(const int &ix, const int &iy, const TGrid &grid, const Lens<Value_type<TGrid>> &lens,
		TVector_c &fPsi_i, TVector_c &fPsi_o)
		{
			using value_type_r = Value_type<TGrid>;
			const value_type_r c_Pi = 3.141592653589793238463;

			int ixy = grid.ind_col(ix, iy);
			value_type_r g2 = grid.g2_shift(ix, iy);

			if((lens.g2_min <= g2)&&(g2 < lens.g2_max))
			{			
				value_type_r chi = g2*(lens.cCs3*g2+lens.cf);
				value_type_r c = c_Pi*lens.beta*lens.sf;
				value_type_r u = 1.0 + c*c*g2;

				c = c_Pi*lens.sf*lens.lambda*g2;
				value_type_r sie = 0.25*c*c;
				c = c_Pi*lens.beta*(lens.Cs3*lens.lambda2*g2-lens.f);
				value_type_r tie = c*c*g2;
				value_type_r sti = exp(-(sie+tie)/u);

				fPsi_o[ixy] = fPsi_i[ixy]*thrust::polar(sti, chi);
			}
			else
			{
				fPsi_o[ixy] = 0;
			}
		}

		template<class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCEINLINE 
		void kernel_xyz(const int &ix, const int &iy, const TGrid &grid, const EELS<Value_type<TGrid>> &eels, TVector_c &k_x, TVector_c &k_y, TVector_c &k_z)
		{
			using value_type_r = Value_type<TGrid>;
			using value_type_c = Value_type<TVector_c>;

			int ixy = grid.ind_col(ix, iy);
			value_type_r gx = grid.gx_shift(ix);
			value_type_r gy = grid.gy_shift(iy);
			value_type_r g2 = gx*gx + gy*gy;
				
			if(g2 < eels.gc2)
			{
				value_type_c pos = thrust::euler(eels.x*gx + eels.y*gy);
				value_type_r lorentz = eels.factor/(g2 + eels.ge2);
				k_x[ixy] = value_type_c(gx*lorentz, 0)*pos;
				k_y[ixy] = value_type_c(gy*lorentz, 0)*pos;
				k_z[ixy] = value_type_c(eels.ge*lorentz, 0)*pos;
			}
			else
			{
				k_x[ixy] = 0;
				k_y[ixy] = 0;
				k_z[ixy] = 0;
			}
		}

		template<class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCEINLINE 
		void kernel_x(const int &ix, const int &iy, const TGrid &grid, const EELS<Value_type<TGrid>> &eels, TVector_c &k_x)
		{
			using value_type_r = Value_type<TGrid>;
			using value_type_c = Value_type<TVector_c>;

			int ixy = grid.ind_col(ix, iy);
			value_type_r gx = grid.gx_shift(ix);
			value_type_r gy = grid.gy_shift(iy);
			value_type_r g2 = gx*gx + gy*gy;
				
			if(g2 < eels.gc2)
			{
				value_type_c pos = thrust::euler(eels.x*gx + eels.y*gy);
				value_type_r lorentz = eels.factor/(g2 + eels.ge2);
				k_x[ixy] = value_type_c(gx*lorentz, 0)*pos;
			}
			else
			{
				k_x[ixy] = 0;
			}
		}

		template<class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCEINLINE 
		void kernel_y(const int &ix, const int &iy, const TGrid &grid, const EELS<Value_type<TGrid>> &eels, TVector_c &k_y)
		{
			using value_type_r = Value_type<TGrid>;
			using value_type_c = Value_type<TVector_c>;

			int ixy = grid.ind_col(ix, iy);
			value_type_r gx = grid.gx_shift(ix);
			value_type_r gy = grid.gy_shift(iy);
			value_type_r g2 = gx*gx + gy*gy;
				
			if(g2 < eels.gc2)
			{
				value_type_c pos = thrust::euler(eels.x*gx + eels.y*gy);
				value_type_r lorentz = eels.factor/(g2 + eels.ge2);
				k_y[ixy] = value_type_c(gy*lorentz, 0)*pos;
			}
			else
			{
				k_y[ixy] = 0;
			}
		}

		template<class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCEINLINE 
		void kernel_z(const int &ix, const int &iy, const TGrid &grid, const EELS<Value_type<TGrid>> &eels, TVector_c &k_z)
		{
			using value_type_r = Value_type<TGrid>;
			using value_type_c = Value_type<TVector_c>;

			int ixy = grid.ind_col(ix, iy);
			value_type_r gx = grid.gx_shift(ix);
			value_type_r gy = grid.gy_shift(iy);
			value_type_r g2 = gx*gx + gy*gy;
				
			if(g2 < eels.gc2)
			{
				value_type_c pos = thrust::euler(eels.x*gx + eels.y*gy);
				value_type_r lorentz = eels.factor/(g2 + eels.ge2);
				k_z[ixy] = value_type_c(eels.ge*lorentz, 0)*pos;
			}
			else
			{
				k_z[ixy] = 0;
			}
		}

		template<class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCEINLINE 
		void kernel_mn1(const int &ix, const int &iy, const TGrid &grid, const EELS<Value_type<TGrid>> &eels, TVector_c &k_mn1)
		{
			using value_type_r = Value_type<TGrid>;
			using value_type_c = Value_type<TVector_c>;

			int ixy = grid.ind_col(ix, iy);
			value_type_r gx = grid.gx_shift(ix);
			value_type_r gy = grid.gy_shift(iy);
			value_type_r g2 = gx*gx + gy*gy;
				
			if(g2 < eels.gc2)
			{
				const value_type_r c_i2i2 = 0.70710678118654746; 

				value_type_c pos = thrust::euler(eels.x*gx + eels.y*gy);
				value_type_r lorentz = c_i2i2*eels.factor/(g2 + eels.ge2);
				value_type_c k_x(gx*lorentz, 0);
				value_type_c k_y(0, gy*lorentz);
				k_mn1[ixy] = (k_x - k_y)*pos;
			}
			else
			{
				k_mn1[ixy] = 0;
			}
		}

		template<class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCEINLINE 
		void kernel_mp1(const int &ix, const int &iy, const TGrid &grid, const EELS<Value_type<TGrid>> &eels, TVector_c &k_mp1)
		{
			using value_type_r = Value_type<TGrid>;
			using value_type_c = Value_type<TVector_c>;

			int ixy = grid.ind_col(ix, iy);
			value_type_r gx = grid.gx_shift(ix);
			value_type_r gy = grid.gy_shift(iy);
			value_type_r g2 = gx*gx + gy*gy;
				
			if(g2 < eels.gc2)
			{
				const value_type_r c_i2i2 = 0.70710678118654746; 

				value_type_c pos = thrust::euler(eels.x*gx + eels.y*gy);
				value_type_r lorentz = c_i2i2*eels.factor/(g2 + eels.ge2);
				value_type_c k_x(gx*lorentz, 0);
				value_type_c k_y(0, gy*lorentz);
				k_mp1[ixy] = (k_x + k_y)*pos;
			}
			else
			{
				k_mp1[ixy] = 0;
			}
		}
	} // host_device_detail

	// Electron scattering factors calculation (feg)
	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	void feg(const ePotential_Type &potential_type, const T &g, const rPP_Coef<T> &feg, T &y)
	{
		switch(potential_type)
		{
			case ePT_Doyle_0_4:
				// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
				host_device_detail::feg_Doyle_0_4(g, feg, y);
				break;
			case ePT_Peng_0_4:
				// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
				host_device_detail::feg_Peng_0_4(g, feg, y);
				break;
			case ePT_Peng_0_12:
				// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
				host_device_detail::feg_Peng_0_12(g, feg, y);
				break;
			case ePT_Kirkland_0_12:
				// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
				host_device_detail::feg_Kirkland_0_12(g, feg, y);
				break;
			case ePT_Weickenmeier_0_12:
				// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
				host_device_detail::feg_Weickenmeier_0_12(g, feg, y);
				break;
			case ePT_Lobato_0_12:
				// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
				host_device_detail::feg_Lobato_0_12(g, feg, y);
				break;
		}
	}

	// Electron scattering factor(feg, dfeg) where dfg is the first derivative along g
	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	void feg_dfeg(const ePotential_Type &potential_type, const T &g, const rPP_Coef<T> &feg, T &y, T &dy)
	{
		switch(potential_type)
		{
			case ePT_Doyle_0_4:
				// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
				host_device_detail::feg_dfeg_Doyle_0_4(g, feg, y, dy);
				break;
			case ePT_Peng_0_4:
				// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
				host_device_detail::feg_dfeg_Peng_0_4(g, feg, y, dy);
				break;
			case ePT_Peng_0_12:
				// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
				host_device_detail::feg_dfeg_Peng_0_12(g, feg, y, dy);
				break;
			case ePT_Kirkland_0_12:
				// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
				host_device_detail::feg_dfeg_Kirkland_0_12(g, feg, y, dy);
				break;
			case ePT_Weickenmeier_0_12:
				// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
				host_device_detail::feg_dfeg_Weickenmeier_0_12(g, feg, y, dy);
				break;
			case ePT_Lobato_0_12:
				// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
				host_device_detail::feg_dfeg_Lobato_0_12(g, feg, y, dy);
				break;
		}
	}

	// Electron scattering factor(fg)
	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	void fxg(const ePotential_Type &potential_type, const T &g, const int &Z, const rPP_Coef<T> &fxg, T &y)
	{
		switch(potential_type)
		{
			case ePT_Doyle_0_4:
				// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
				host_device_detail::fxg_Doyle_0_4(g, Z, fxg, y);
				break;
			case ePT_Peng_0_4:
				// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
				host_device_detail::fxg_Peng_0_4(g, Z, fxg, y);
				break;
			case ePT_Peng_0_12:
				// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
				host_device_detail::fxg_Peng_0_12(g, Z, fxg, y);
				break;
			case ePT_Kirkland_0_12:
				// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
				host_device_detail::fxg_Kirkland_0_12(g, Z, fxg, y);
				break;
			case ePT_Weickenmeier_0_12:
				// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
				host_device_detail::fxg_Weickenmeier_0_12(g, fxg, y);
				break;
			case ePT_Lobato_0_12:
				// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
				host_device_detail::fxg_Lobato_0_12(g, fxg, y);
				break;
		}
	}

	// Electron scattering factor(fg, dfg) where dfg is the first derivative along g
	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	void fxg_dfxg(const ePotential_Type &potential_type, const T &g, const int &Z, const rPP_Coef<T> &fxg, T &y, T &dy)
	{
		switch(potential_type)
		{
			case ePT_Doyle_0_4:
				// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
				host_device_detail::fxg_dfxg_Doyle_0_4(g, Z, fxg, y, dy);
				break;
			case ePT_Peng_0_4:
				// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
				host_device_detail::fxg_dfxg_Peng_0_4(g, Z, fxg, y, dy);
				break;
			case ePT_Peng_0_12:
				// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
				host_device_detail::fxg_dfxg_Peng_0_12(g, Z, fxg, y, dy);
				break;
			case ePT_Kirkland_0_12:
				// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
				host_device_detail::fxg_dfxg_Kirkland_0_12(g, Z, fxg, y, dy);
				break;
			case ePT_Weickenmeier_0_12:
				// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
				host_device_detail::fxg_dfxg_Weickenmeier_0_12(g, fxg, y, dy);
				break;
			case ePT_Lobato_0_12:
				// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
				host_device_detail::fxg_dfxg_Lobato_0_12(g, fxg, y, dy);
				break;
		}
	}

	// Electron density (Pr)
	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	void Pr(const ePotential_Type &potential_type, const T &r, const rPP_Coef<T> &Pr, T &y)
	{
		switch(potential_type)
		{
			case ePT_Doyle_0_4:
				// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
				host_device_detail::Pr_Doyle_0_4(r, Pr, y);
				break;
			case ePT_Peng_0_4:
				// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
				host_device_detail::Pr_Peng_0_4(r, Pr, y);
				break;
			case ePT_Peng_0_12:
				// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
				host_device_detail::Pr_Peng_0_12(r, Pr, y);
				break;
			case ePT_Kirkland_0_12:
				// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
				host_device_detail::Pr_Kirkland_0_12(r, Pr, y);
				break;
			case ePT_Weickenmeier_0_12:
				// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
				host_device_detail::Pr_Weickenmeier_0_12(r, Pr, y);
				break;
			case ePT_Lobato_0_12:
				// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
				host_device_detail::Pr_Lobato_0_12(r, Pr, y);
				break;
		}
	}

	// Electron density (Pr, dPr) where dPr is the first derivative along r
	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	void Pr_dPr(const ePotential_Type &potential_type, const T &r, const rPP_Coef<T> &Pr, T &y, T &dy)
	{
		switch(potential_type)
		{
			case ePT_Doyle_0_4:
				// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
				host_device_detail::Pr_dPr_Doyle_0_4(r, Pr, y, dy);
				break;
			case ePT_Peng_0_4:
				// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
				host_device_detail::Pr_dPr_Peng_0_4(r, Pr, y, dy);
				break;
			case ePT_Peng_0_12:
				// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
				host_device_detail::Pr_dPr_Peng_0_12(r, Pr, y, dy);
				break;
			case ePT_Kirkland_0_12:
				// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
				host_device_detail::Pr_dPr_Kirkland_0_12(r, Pr, y, dy);
				break;
			case ePT_Weickenmeier_0_12:
				// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
				host_device_detail::Pr_dPr_Weickenmeier_0_12(r, Pr, y, dy);
				break;
			case ePT_Lobato_0_12:
				// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
				host_device_detail::Pr_dPr_Lobato_0_12(r, Pr, y, dy);
				break;
		}
	}

	// Potential calculation(Vr)
	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	void Vr(const ePotential_Type &potential_type, const T &r, const rPP_Coef<T> &Vr, T &y)
	{
		switch(potential_type)
		{
			case ePT_Doyle_0_4:
				// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
				host_device_detail::Vr_Doyle_0_4(r, Vr, y);
				break;
			case ePT_Peng_0_4:
				// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
				host_device_detail::Vr_Peng_0_4(r, Vr, y);
				break;
			case ePT_Peng_0_12:
				// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
				host_device_detail::Vr_Peng_0_12(r, Vr, y);
				break;
			case ePT_Kirkland_0_12:
				// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
				host_device_detail::Vr_Kirkland_0_12(r, Vr, y);
				break;
			case ePT_Weickenmeier_0_12:
				// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
				host_device_detail::Vr_Weickenmeier_0_12(r, Vr, y);
				break;
			case ePT_Lobato_0_12:
				// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
				host_device_detail::Vr_Lobato_0_12(r, Vr, y);
				break;
		}
	}

	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	void Vr_dVr(const ePotential_Type &potential_type, const T &r, const rPP_Coef<T> &Vr, T &y, T &dy)
	{
		switch(potential_type)
		{
			case ePT_Doyle_0_4:
				// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
				host_device_detail::Vr_dVr_Doyle_0_4(r, Vr, y, dy);
				break;
			case ePT_Peng_0_4:
				// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
				host_device_detail::Vr_dVr_Peng_0_4(r, Vr, y, dy);
				break;
			case ePT_Peng_0_12:
				// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
				host_device_detail::Vr_dVr_Peng_0_12(r, Vr, y, dy);
				break;
			case ePT_Kirkland_0_12:
				// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
				host_device_detail::Vr_dVr_Kirkland_0_12(r, Vr, y, dy);
				break;
			case ePT_Weickenmeier_0_12:
				// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
				host_device_detail::Vr_dVr_Weickenmeier_0_12(r, Vr, y, dy);
				break;
			case ePT_Lobato_0_12:
				// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
				host_device_detail::Vr_dVr_Lobato_0_12(r, Vr, y, dy);
				break;
		}
	}

	// Projected potential (VR)
	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	void VR(const ePotential_Type &potential_type, const T &R, const rPP_Coef<T> &VR, const rQ1<T> &Qz_0_I, T &y)
	{
		switch(potential_type)
		{
			case ePT_Doyle_0_4:
				// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
				host_device_detail::VR_Doyle_0_4(R, VR, y);
				break;
			case ePT_Peng_0_4:
				// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
				host_device_detail::VR_Peng_0_4(R, VR, y);
				break;
			case ePT_Peng_0_12:
				// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
				host_device_detail::VR_Peng_0_12(R, VR, y);
				break;
			case ePT_Kirkland_0_12:
				// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
				host_device_detail::VR_Kirkland_0_12(R, VR, y);
				break;
			case ePT_Weickenmeier_0_12:
				// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
				host_device_detail::VR_Weickenmeier_0_12(R, VR, Qz_0_I, y);
				break;
			case ePT_Lobato_0_12:
				// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
				host_device_detail::VR_Lobato_0_12(R, VR, y);
				break;
		}
	}

	// Projected potential (VR, dVR) where dVr is the first derivative along R
	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	void VR_dVR(const ePotential_Type &potential_type, const T &R, const rPP_Coef<T> &VR, const rQ1<T> &Qz_0_I, T &y, T &dy)
	{
		switch(potential_type)
		{
			case ePT_Doyle_0_4:
				// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
				host_device_detail::VR_dVR_Doyle_0_4(R, VR, y, dy);
				break;
			case ePT_Peng_0_4:
				// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
				host_device_detail::VR_dVR_Peng_0_4(R, VR, y, dy);
				break;
			case ePT_Peng_0_12:
				// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
				host_device_detail::VR_dVR_Peng_0_12(R, VR, y, dy);
				break;
			case ePT_Kirkland_0_12:
				// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
				host_device_detail::VR_dVR_Kirkland_0_12(R, VR, y, dy);
				break;
			case ePT_Weickenmeier_0_12:
				// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
				host_device_detail::VR_dVR_Weickenmeier_0_12(R, VR, Qz_0_I, y, dy);
				break;
			case ePT_Lobato_0_12:
				// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
				host_device_detail::VR_dVR_Lobato_0_12(R, VR, y, dy);
				break;
		}
	}

	// Projected potential (Vz)[z0, ze]
	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	void Vz(const ePotential_Type &potential_type, const T &z0, const T &ze, const T &R, const rPP_Coef<T> &Vr, const rQ1<T> &Qz_a_b, T &y)
	{
		switch(potential_type)
		{
			case ePT_Doyle_0_4:
				// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
				host_device_detail::Vz_Doyle_0_4(z0, ze, R, Vr, Qz_a_b, y);
				break;
			case ePT_Peng_0_4:
				// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
				host_device_detail::Vz_Peng_0_4(z0, ze, R, Vr, Qz_a_b, y);
				break;
			case ePT_Peng_0_12:
				// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
				host_device_detail::Vz_Peng_0_12(z0, ze, R, Vr, Qz_a_b, y);
				break;
			case ePT_Kirkland_0_12:
				// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
				host_device_detail::Vz_Kirkland_0_12(z0, ze, R, Vr, Qz_a_b, y);
				break;
			case ePT_Weickenmeier_0_12:
				// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
				host_device_detail::Vz_Weickenmeier_0_12(z0, ze, R, Vr, Qz_a_b, y);
				break;
			case ePT_Lobato_0_12:
				// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
				host_device_detail::Vz_Lobato_0_12(z0, ze, R, Vr, Qz_a_b, y);
				break;
		}
	}

	// Projected potential (Vz, dVz)[z0, ze] where dVr is the first derivative along R
	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	void Vz_dVz(const ePotential_Type &potential_type, const T &z0, const T &ze, const T &R, const rPP_Coef<T> &Vr, const rQ1<T> &Qz_a_b, T &y, T &dy)
	{
		switch(potential_type)
		{
			case ePT_Doyle_0_4:
				// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
				host_device_detail::Vz_dVz_Doyle_0_4(z0, ze, R, Vr, Qz_a_b, y, dy);
				break;
			case ePT_Peng_0_4:
				// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
				host_device_detail::Vz_dVz_Peng_0_4(z0, ze, R, Vr, Qz_a_b, y, dy);
				break;
			case ePT_Peng_0_12:
				// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
				host_device_detail::Vz_dVz_Peng_0_12(z0, ze, R, Vr, Qz_a_b, y, dy);
				break;
			case ePT_Kirkland_0_12:
				// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
				host_device_detail::Vz_dVz_Kirkland_0_12(z0, ze, R, Vr, Qz_a_b, y, dy);
				break;
			case ePT_Weickenmeier_0_12:
				// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
				host_device_detail::Vz_dVz_Weickenmeier_0_12(z0, ze, R, Vr, Qz_a_b, y, dy);
				break;
			case ePT_Lobato_0_12:
				// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
				host_device_detail::Vz_dVz_Lobato_0_12(z0, ze, R, Vr, Qz_a_b, y, dy);
				break;
		}
	}
	
	template<ePotential_Type potential_type, class T>
	DEVICE_CALLABLE FORCEINLINE
	void Vr_dVrir(const T &r, T *cl, T *cnl, T f, T &Vr, T &dVrir)
	{
		switch(potential_type)
		{
			case ePT_Doyle_0_4:
				// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
				host_device_detail::Vr_dVrir_Doyle_0_4<T>(r, cl, cnl, Vr, dVrir);
				break;
			case ePT_Peng_0_4:
				// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
				host_device_detail::Vr_dVrir_Peng_0_4_12<T>(r, cl, cnl, Vr, dVrir);
				break;
			case ePT_Peng_0_12:
				// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
				host_device_detail::Vr_dVrir_Peng_0_4_12<T>(r, cl, cnl, Vr, dVrir);
				break;
			case ePT_Kirkland_0_12:
				// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
				host_device_detail::Vr_dVrir_Kirkland_0_12<T>(r, cl, cnl, Vr, dVrir);
				break;
			case ePT_Weickenmeier_0_12:
				// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
				host_device_detail::Vr_dVrir_Weickenmeier_0_12<T>(r, cl, cnl, Vr, dVrir);
				break;
			case ePT_Lobato_0_12:
				// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
				host_device_detail::Vr_dVrir_Lobato_0_12<T>(r, cl, cnl, Vr, dVrir);
				break;
		}

		Vr *= f;
		dVrir *= f;
	}

	/***************************************************************************/
	/***************************************************************************/

	namespace functor
	{
		template<class TGrid>
		struct phase_components
		{
			using T = Value_type<TGrid>;

			const TGrid grid;
			const T w;
			const T gxu;
			const T gyu;
			phase_components(TGrid grid_i, T w_i, T gxu_i, T gyu_i): grid(grid_i), w(w_i), gxu(gxu_i), gyu(gyu_i){}

			template<class Ttuple>
			__host__ __device__
			void operator()(Ttuple t)
			{
				int ix = thrust::get<0>(t);
				if(ix < grid.nx)
				{
					T Rx = grid.Rx_shift(ix);
					thrust::get<1>(t) = thrust::euler(w*Rx*gxu);
				}

				int iy = ix;
				if(iy < grid.ny)
				{
					T Ry = grid.Ry_shift(iy);
					thrust::get<2>(t) = thrust::euler(w*Ry*gyu);
				}
			}
		};

		template<class TGrid>
		struct propagator_components
		{
			using T = Value_type<TGrid>;

			const TGrid grid;
			const T w;
			const T gxu;
			const T gyu;
			propagator_components(TGrid grid_i, T w_i, T gxu_i, T gyu_i): grid(grid_i), w(w_i), gxu(gxu_i), gyu(gyu_i){}

			template<class Ttuple>
			__host__ __device__
			void operator()(Ttuple t)
			{
				int ix = thrust::get<0>(t);
				if(ix < grid.nx)
				{
					T gx = grid.gx_shift(ix) + gxu;
					thrust::get<1>(t) = thrust::euler(w*gx*gx);
				}

				int iy = ix;
				if(iy < grid.ny)
				{
					T gy = grid.gy_shift(iy) + gyu;
					thrust::get<2>(t) = thrust::euler(w*gy*gy);
				}
			}
		};

		template<class T>
		struct transmission_funtion
		{
			const eElec_Spec_Int_Model elec_spec_int_model;

			const T w;
			transmission_funtion(T w_i, eElec_Spec_Int_Model ElecSpecIntModel_i): w(w_i), elec_spec_int_model(ElecSpecIntModel_i){}

			template<class U>
			__host__ __device__
			complex<T> operator()(const U &x) const 
			{ 
				T theta = w*static_cast<T>(x);
				return (elec_spec_int_model == eESIM_Weak_Phase_Object)?complex<T>(1.0, theta):thrust::euler(theta);
			}
		};

		template<class T>
		struct scale
		{
			const T w;
			scale(T w_i = T()): w(w_i){}

			__host__ __device__
			T operator()(const T &x) const { return w*x; }
		};

		template<class T>
		struct square
		{
			template<class U>
			__host__ __device__
			T operator()(const U &x) const { return thrust::norm(x); }
		};

		template<class T>
		struct square_scale
		{
			const T w;
			square_scale(T w_i = T()): w(w_i){}

			template<class U>
			__host__ __device__
			T operator()(const U &x) const { return w*thrust::norm(x); }
		};

		template<class T>
		struct add_scale
		{
			const T w;
			add_scale(T w_i = T()): w(w_i){}

			__host__ __device__
			T operator()(const T &lhs, const T &rhs) const{ return w*lhs + rhs; }
		};

		template<class T>
		struct add_scale_i
		{
			const T w;
			add_scale_i(T w_i = T()): w(w_i){}

			__host__ __device__
			T operator()(const T &lhs, const T &rhs) const{ return w*(lhs + rhs); }
		};

		template<class T>
		struct add_square
		{
			template<class U>
			__host__ __device__
			T operator()(const U &lhs, const T &rhs) const { return thrust::norm(lhs) + rhs; }
		};

		template<class T>
		struct add_square_i
		{
			template<class U>
			__host__ __device__
			T operator()(const U &lhs, const U &rhs) const { return thrust::norm(lhs) + thrust::norm(rhs); }
		};

		template<class T>
		struct add_square_scale
		{
			const T w;
			add_square_scale(T w_i = T()): w(w_i){}

			template<class U>
			__host__ __device__
			T operator()(const U &lhs, const T &rhs) const { return w*thrust::norm(lhs) + rhs; }
		};

		template<class T>
		struct add_square_scale_i
		{
			const T w;
			add_square_scale_i(T w_i = T()): w(w_i){}

			template<class U>
			__host__ __device__
			T operator()(const U &lhs, const U &rhs) const { return w*(thrust::norm(lhs) + thrust::norm(rhs)); }
		};

	} //namespace functor

	template<class TVector>
	void scale(TVector &M_io, const Value_type<TVector> &w_i)
	{
		using value_type = Value_type<TVector>;
		thrust::transform(M_io.begin(), M_io.end(), M_io.begin(), functor::scale<value_type>(w_i));
	}

	template<class TVector>
	void fill(TVector &M_io, const Value_type<TVector> &value_i)
	{
		thrust::fill(M_io.begin(), M_io.end(), value_i);
	}


	template<class TVector_1, class TVector_2>
	void assign(TVector_1 &M_i, TVector_2 &M_o)
	{
		M_o.assign(M_i.begin(), M_i.end());
	}

	template<class TVector_1, class TVector_2>
	void assign_square(TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), functor::square<value_type>());
	}

	template<class TVector>
	void assign_scale(Value_type<TVector> w_i, TVector &M_i, TVector &M_o)
	{
		using value_type = Value_type<TVector>;
		thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), functor::scale<value_type>(w_i));
	}

	template<class TVector_1, class TVector_2>
	void assign_square_scale(Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), functor::square_scale<value_type>(w_i));
	}


	template<class TVector_1, class TVector_2>
	void add(const TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), M_o.begin(), thrust::plus<value_type>());
	}

	template<class TVector_1, class TVector_2>
	void add(const TVector_1 &M1_i, const TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_1>;
		thrust::transform(M1_i.begin(), M1_i.end(), M2_i.begin(), M_o.begin(), thrust::plus<value_type>());
	}


	template<class TVector_1, class TVector_2>
	void add_square(const TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_1>;
		thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), M_o.begin(), functor::add_square<value_type>());
	}

	template<class TVector_1, class TVector_2>
	void add_square(const TVector_1 &M1_i, const TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M1_i.begin(), M1_i.end(), M2_i.begin(), M_o.begin(), functor::add_square_i<value_type>());
	}

	template<class TVector_1, class TVector_2>
	void add_scale(Value_type<TVector_1> w_i, const TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_1>;
		thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), M_o.begin(), functor::add_scale<value_type>(w_i));
	}

	template<class TVector_1, class TVector_2>
	void add_scale(Value_type<TVector_1> w_i, const TVector_1 &M1_i, const TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_1>;

		thrust::transform(M1_i.begin(), M1_i.end(), M2_i.begin(), M_o.begin(), functor::add_scale_i<value_type>());
	}

	template<class TVector_1, class TVector_2>
	void add_square_scale(Value_type<TVector_2> w_i, const TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), M_o.begin(), functor::add_square_scale<value_type>(w_i));
	}

	template<class TVector_1, class TVector_2>
	void add_square_scale(Value_type<TVector_1> w_i, const TVector_1 &M1_i, const TVector_1 &M2_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		thrust::transform(M1_i.begin(), M1_i.end(), M2_i.begin(), M_o.begin(), functor::add_square_scale_i<value_type>(w_i));
	}

	template<class TVector>
	void multiply(const TVector &M_i, TVector &M_io)
	{
		using value_type = Value_type<TVector>;
		thrust::transform(M_i.begin(), M_i.end(), M_io.begin(), M_io.begin(), thrust::multiplies<value_type>());
	}

	template<class TVector>
	void multiply(const TVector &M1_i, const TVector &M2_i, TVector &M_o)
	{
		using value_type = Value_type<TVector>;
		thrust::transform(M1_i.begin(), M1_i.end(), M2_i.begin(), M_o.begin(), thrust::multiplies<value_type>());
	}

	template<class TGrid, class TVector>
	Value_type<TVector> sum(TGrid &grid, TVector &M_i)
	{
		Value_type<TVector> sum_o = thrust::reduce(M_i.begin(), M_i.end());
		return sum_o;
	}

	template<class TGrid, class TVector>
	Value_type<TGrid> sum_square(TGrid &grid, TVector &M_i)
	{
		using value_type_r = Value_type<TGrid>;
		value_type_r sum_o = thrust::transform_reduce(M_i.begin(), M_i.end(), 
		functor::square<value_type_r>(), static_cast<value_type_r>(0), thrust::plus<value_type_r>());
		return sum_o;
	}

	template<class TGrid, class TVector>
	void phase_components(TGrid &grid, const Value_type<TGrid> &gxu, const Value_type<TGrid> &gyu, TVector &V_x, TVector &V_y)
	{
		thrust::counting_iterator<int> first(0);
		thrust::counting_iterator<int> last = first + grid.nx_ny_max();

		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(first, V_x.begin(), V_y.begin())), 
						 thrust::make_zip_iterator(thrust::make_tuple(last, V_x.end(), V_y.end())), 
						 functor::phase_components<TGrid>(grid, c_2Pi, gxu, gyu));
	}

	template<class TGrid, class TVector>
	void propagator_components(TGrid &grid, const Value_type<TGrid> &gxu, const Value_type<TGrid> &gyu, const Value_type<TGrid> &w, TVector &V_x, TVector &V_y)
	{
		thrust::counting_iterator<int> first(0);
		thrust::counting_iterator<int> last = first + grid.nx_ny_max();

		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(first, V_x.begin(), V_y.begin())), 
						 thrust::make_zip_iterator(thrust::make_tuple(last, V_x.end(), V_y.end())), 
						 functor::propagator_components<TGrid>(grid, w, gxu, gyu));
	}

	template<class TGrid, class TFFT2, class TVector>
	void propagate(TGrid &grid, TFFT2 &fft2, eSpace space, TVector &prop_x_i, TVector &prop_y_i, TVector &psi_i, TVector &psi_o)
	{
		fft2.forward(psi_i, psi_o); 

		propagator_multiplication(grid, prop_x_i, prop_y_i, psi_o, psi_o);

		if(space == eS_Real)
		{
			fft2.inverse(psi_o);
		}
	}

	template<class TGrid, class TFFT2, class TVector_1, class TVector_2>
	void transmission_funtion(TGrid &grid, TFFT2 &fft2, eElec_Spec_Int_Model elec_spec_int_model, const Value_type<TGrid> w, TVector_1 &V0_i, TVector_2 &Trans_o)
	{	
		using value_type_r = Value_type<TGrid>;

		thrust::transform(V0_i.begin(), V0_i.end(), Trans_o.begin(), functor::transmission_funtion<value_type_r>(w, elec_spec_int_model));

		if(grid.bwl)
		{
			fft2.forward(Trans_o);
			bandwidth_limit(grid, 0, grid.gl_max, grid.inxy, Trans_o);
			fft2.inverse(Trans_o);
		}
	}

	template<class TGrid, class TVector_i, class TVector_o>
	void to_host_shift(TGrid &grid, TVector_i &M_i, TVector_o &M_o)
	{
		copy_to_host(grid, M_i, M_o);
		fft2_shift(grid, M_o);
	}

} //namespace multem

#endif