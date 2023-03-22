/*
 * This file is part of MULTEM.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
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
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef CGPU_FCNS_H
#define CGPU_FCNS_H

#include <type_traits>
#include <algorithm>

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "lin_alg_def.cuh"

#include <thrust/complex.h>
#include <thrust/swap.h>
#include <thrust/extrema.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/for_each.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/tuple.h>
#include <thrust/execution_policy.h>

namespace mt
{
	// E. J. Kirkland - Advanced computing in electron microscopy page: 10-13
	// Input: E_0(keV), Output: lambda (electron wave)
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T get_lambda(const T &E_0)
	{
		T emass = 510.99906;
		T hc = 12.3984244;
		T lambda = hc/sqrt(E_0*(2*emass + E_0));
		return lambda;
	}

	// E. J. Kirkland - Advanced computing in electron microscopy page: 10-13
	// Input: E_0(keV), Output: sigma (Interaction parameter)
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T get_sigma(const T &E_0)
	{
		const T c_Pi = 3.141592653589793238463;
		T emass = 510.99906;
		T x = (emass + E_0)/(2*emass + E_0);
		T sigma = 2*c_Pi*x/(get_lambda(E_0)*E_0);
		return sigma;
	}

	// E. J. Kirkland - Advanced computing in electron microscopy page: 10-13
	// Input: E_0(keV), Output: gamma(relativistic factor)
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T get_gamma(const T &E_0)
	{
		T emass = 510.99906;
		T gamma = 1 + E_0/emass;
		return gamma;
	}

	// Input: theta(mrad), E_0(keV), Output: A^-1
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T rad_2_rAngs(const T &E_0, const T &theta)
	{
		return sin(theta)/get_lambda(E_0);
	}

	// Marc De Graef - Introduction to Conventional Transmission Electron Microscopy page: 608
	// hwhm: Half width at half maximum
	// sigma: Standard deviation
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T hwhm_2_sigma(const T &v)
	{
		T c_hwhm_2_sigma = 0.84932180028801907; // hwhm to sigma 1/(sqrt(2*log(2)))
		return v*c_hwhm_2_sigma;
	}

	// Marc De Graef - Introduction to Conventional Transmission Electron Microscopy page: 608
	// fwhm: Full width at half maximum
	// sigma: Standard deviation
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T fwhm_2_sigma(const T &v)
	{
		T c_fwhm_2_sigma = 0.42466090014400953; // fwhm to sigma 1/(2*sqrt(2*log(2)))
		return v*c_fwhm_2_sigma;
	}

	// Marc De Graef - Introduction to Conventional Transmission Electron Microscopy page: 608
	// iehwgd: e^-1 half-width value of the Gaussian distribution
	// sigma: Standard deviation
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T iehwgd_2_sigma(const T &v)
	{
		T c_iehwgd_2_sigma = 0.70710678118654746; // iehwgd to sigma 1/sqrt(2)
		return v*c_iehwgd_2_sigma;
	}

	// Marc De Graef - Introduction to Conventional Transmission Electron Microscopy page: 608
	// sigma: Standard deviation
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T rad_2_sigma(const T &E_0, const T &theta)
	{
		T q0 = sin(theta)/get_lambda(E_0);
		return iehwgd_2_sigma(q0);
	}

	// E. J. Kirkland - Advanced computing in electron microscopy page: 10-13
	// Input: E_0(keV), Output: gamma*lambda/c_Potf
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T get_Vr_factor(const T &E_0, const T &theta)
	{
		T c_Potf = 47.877645145863056;
		T fPot = get_gamma(E_0)*get_lambda(E_0)/(c_Potf*cos(theta));
		return fPot;
	}

	// E. J. Kirkland - Advanced computing in electron microscopy page: 33
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T get_Scherzer_defocus(const T &E_0, const T &c_30)
	{
		T lambda = get_lambda(E_0);
		T n = 1.0;
		return -copysign(sqrt((2*n-0.5)*fabs(c_30)*lambda), c_30);
	}

	// E. J. Kirkland - Advanced computing in electron microscopy page: 33
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T get_Scherzer_aperture(const T &E_0, const T &c_30)
	{
		T lambda = get_lambda(E_0);
		T n = 1.0;
		return pow(4*(2*n-0.5)*lambda/fabs(c_30), 0.25);
	}

	// E. J. Kirkland - Advanced computing in electron microscopy page: 33
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	void get_Scherzer_conditions(const T &E_0, const T &c_30, T &defocus, T &aperture)
	{
		defocus = get_Scherzer_defocus(E_0, c_30);
		aperture = get_Scherzer_aperture(E_0, c_30);
	}

	// new size
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	int get_new_size(const int &n, T factor)
	{
		return max(static_cast<int>(ceil(n*factor)), 1);
	}

	template <class T>
	FORCE_INLINE 
	Vector<T, e_host> get_rotation_matrix(const T &theta, const r3d<T> &u0)
	{
		Vector<T, e_host> Rm(9);
		T alpha = 1-cos(theta);
		alpha = (isZero(alpha)?0:alpha);
		T beta = sin(theta);
		beta = (isZero(beta)?0:beta);
		Rm[0] = 1.0 + alpha*(u0.x*u0.x-1);
		Rm[1] = u0.y*u0.x*alpha + u0.z*beta;
		Rm[2] = u0.z*u0.x*alpha - u0.y*beta;

		Rm[3] = u0.x*u0.y*alpha - u0.z*beta;
		Rm[4] = 1.0 + alpha*(u0.y*u0.y-1);
		Rm[5] = u0.z*u0.y*alpha + u0.x*beta;

		Rm[6] = u0.x*u0.z*alpha + u0.y*beta;
		Rm[7] = u0.y*u0.z*alpha - u0.x*beta;
		Rm[8] = 1.0 + alpha*(u0.z*u0.z-1);
		return Rm;
	}

	// distance from point to line
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T get_dist_from_p2l(const T &a, const T &b, const T &c, const T &x0, const T &y0)
	{
		return fabs(a*x0+b*y0+c)/sqrt(a*a+b*b);
	}

	// calculate intersection points
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	void get_int_pt_from_p2l(const T &a, const T &b, const T &c, const T &x0, const T &y0, const T &x, const T &y)
	{
		x = (b*(b*x0-a*y0)-a*c)/(a*a+b*b);
		y = -(a*(b*x0-a*y0)-b*c)/(a*a+b*b);
	}

	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T m_bound(const T &x_0, const T &x, const T &x_e)
	{
		return max(x_0, min(x, x_e));
	}

	template <eDevice dev>
	void synchronize_every(int i_sync, int n_sync)
	{

	}

#ifdef __CUDACC__
	template <>
	void synchronize_every<e_device>(int i_sync, int n_sync)
	{
		if(i_sync % n_sync == 0)
		{
			cudaDeviceSynchronize();
		}
	}
#endif

	namespace host_device_detail
	{
		// Kahan summation algorithm
		// https:// en.wikipedia.org/wiki/Kahan_summation_algorithm
		template <class T>
		DEVICE_CALLABLE FORCE_INLINE 
		void kh_sum(T &sum_v, T v, T &error)
		{
			v = v - error;
			T t = sum_v + v;
			error = (t-sum_v)-v;
			sum_v = t;
		}

		template <class TFn, class T>
		inline
		T Root_Finder(TFn fn, T x0, T xe, const T Tol = 1e-8, const int itMax = 200)
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
			}while((fabs(fx)>Tol) && (it < itMax));
 
			return x;
		}

		template <class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void int_fx_x0_xe(TFn fn, const Value_type<TrPP_Coef> &x0, const Value_type<TrPP_Coef> &xe, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rq1, Value_type<TrPP_Coef> &y)
		{
			Value_type<TrPP_Coef> a = 0.5*(xe-x0); 
			Value_type<TrPP_Coef> b = 0.5*(xe+x0);
			Value_type<TrPP_Coef> xi, yi;
			y = 0;

			for(auto i = 0; i< rq1.m_size; i++)
			{
				xi = a*rq1.x[i] + b;
				fn(xi, rcoef, yi);
				y += a*rq1.w[i]*yi;
			}
		}

		template <class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void int_fx_x0_xe(TFn fn, const Value_type<TrPP_Coef> &x0, const Value_type<TrPP_Coef> &xe, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rq1, Value_type<TrPP_Coef> &y1, Value_type<TrPP_Coef> &y2)
		{
			Value_type<TrPP_Coef> a = 0.5*(xe-x0); 
			Value_type<TrPP_Coef> b = 0.5*(xe+x0);
			Value_type<TrPP_Coef> xi, y1i, y2i;
			y1 = y2 = 0;

			for(auto i = 0; i< rq1.m_size; i++)
			{
				xi = a*rq1.x[i] + b;
				fn(xi, rcoef, y1i, y2i);
				y1 += a*rq1.w[i]*y1i;
				y2 += a*rq1.w[i]*y2i;
			}
		}

		template <class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void int_fx_x0_pInfty(TFn fn, const Value_type<TrPP_Coef> &x0, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, TrQ1&rq1, Value_type<TrPP_Coef> &y)
		{
			using T = Value_type<TrPP_Coef>;
			T xi, yi;
			y = 0;

			for(auto i = 0; i< rq1.m_size; i++)
			{
				xi = rq1.x[i] + x0;
				fn(xi, rcoef, yi);
				y += rq1.w[i]*yi;
			}
		}

		template <class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void int_fx_x0_pInfty(TFn fn, const Value_type<TrPP_Coef> &x0, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rq1, Value_type<TrPP_Coef> &y1, Value_type<TrPP_Coef> &y2)
		{
			using T = Value_type<TrPP_Coef>;
			T xi, y1i, y2i;
			y1 = y2 = 0;

			for(auto i = 0; i< rq1.m_size; i++)
			{
				xi = rq1.x[i] + x0;
				fn(xi, rcoef, y1i, y2i);
				y1 += rq1.w[i]*y1i;
				y2 += rq1.w[i]*y2i;
			}
		}

		template <class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void int_Vz_z0_ze(TFn fn, const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &R, const TrPP_Coef &rcoef, const TrQ1 &rq1, Value_type<TrPP_Coef> &y)
		{
			using T = Value_type<TrPP_Coef>;
			bool split = (z0<0) && (0<ze);
			T a = (split)?-0.5*z0:0.5*(ze-z0); 
			T b = (split)?0.5*z0:0.5*(ze+z0);
			T zi, ri, yi, R2 = R*R;
			y = 0;

			for(auto i = 0; i< rq1.m_size; i++)
			{
				zi = a*rq1.x[i] + b;
				ri = sqrt(zi*zi + R2);
				fn(ri, rcoef, yi);
				y += a*rq1.w[i]*yi;
			}

			if(split)
			{
				a = b = 0.5*ze;
				for(auto i = 0; i< rq1.m_size; i++)
				{
					zi = a*rq1.x[i] + b;
					ri = sqrt(zi*zi + R2);
					fn(ri, rcoef, yi);
					y += a*rq1.w[i]*yi;
				}
			}
		}

		template <class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void int_Vz_dVz_z0_ze(TFn fn, const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &R, const TrPP_Coef &rcoef, const TrQ1 &rq1, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			using T = Value_type<TrPP_Coef>;
			bool split = (z0<0) && (0<ze);
			T a = (split)?-0.5*z0:0.5*(ze-z0); 
			T b = (split)?0.5*z0:0.5*(ze+z0);
			T zi, ri, yi, dyi, R2 = R*R;
			y = dy = 0;

			for(auto i = 0; i< rq1.m_size; i++)
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
				for(auto i = 0; i< rq1.m_size; i++)
				{
					zi = a*rq1.x[i] + b;
					ri = sqrt(zi*zi + R2);
					fn(ri, rcoef, yi, dyi);
					y += a*rq1.w[i]*yi; 
					dy += R*a*rq1.w[i]*dyi/ri;
				}
			}
		}

		template <class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void int_VR(TFn fn, const Value_type<TrPP_Coef> &R, const TrPP_Coef &rcoef, const TrQ1 &rq1, Value_type<TrPP_Coef> &y)
		{
			using T = Value_type<TrPP_Coef>;
			T zi, ri, yi, R2 = R*R;
			y = 0;
			for(auto i = 0; i< rq1.m_size; i++)
			{
				zi = rq1.x[i];
				ri = sqrt(zi*zi + R2);
				fn(ri, rcoef, yi);
				y += rq1.w[i]*yi;
			}
			y *= 2;
		}

		template <class TFn, class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void int_VR_dVR(TFn fn, const Value_type<TrPP_Coef> &R, const TrPP_Coef &rcoef, const TrQ1 &rq1, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			using T = Value_type<TrPP_Coef>;
			T zi, ri, yi, dyi, R2 = R*R;
			y = dy = 0;
			for(auto i = 0; i< rq1.m_size; i++)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void add_Exponential_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, bool reset = true)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void add_Exponential_dExponential_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy, bool reset = true)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void add_Gauss_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, bool reset = true)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void add_Gauss_dGauss_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy, bool reset = true)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void add_Lorentzian_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, bool reset = true)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void add_Lorentzian_dLorentzian_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy, bool reset = true)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void add_Yukawa_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, bool reset = true)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void add_Yukawa_dYukawa_Fn(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy, bool reset = true)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_Gauss_feg(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, bool reset = true)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_dPr_Gauss_feg(const Value_type<TrPP_Coef> &x, const int &n_0, const int &n_e, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy, bool reset = true)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void feg_Doyle_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gauss_Fn(x, 0, 4, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void feg_dfeg_Doyle_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gauss_dGauss_Fn(x, 0, 4, rcoef, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void feg_Peng_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gauss_Fn(x, 0, 5, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void feg_dfeg_Peng_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gauss_dGauss_Fn(x, 0, 5, rcoef, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void feg_Peng_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gauss_Fn(x, 0, 5, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void feg_dfeg_Peng_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gauss_dGauss_Fn(x, 0, 5, rcoef, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void feg_Kirkland_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Lorentzian_Fn(x, 0, 3, rcoef, y);
			add_Gauss_Fn(x, 3, 6, rcoef, y, false);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void feg_dfeg_Kirkland_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Lorentzian_dLorentzian_Fn(x, 0, 3, rcoef, y, dy);
			add_Gauss_dGauss_Fn(x, 3, 6, rcoef, y, dy, false);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void feg_Weickenmeier_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			using T = Value_type<TrPP_Coef>;
			T x2 = x*x;
	
			y = 0;
			if(nonZero(x))
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void feg_dfeg_Weickenmeier_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			using T = Value_type<TrPP_Coef>;
			T t, x2 = x*x;
			y = dy = 0;
			if(nonZero(x))
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void feg_Lobato_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void feg_dfeg_Lobato_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void feg_Peng_ion_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gauss_Fn(x, 0, 5, rcoef, y);
			y += rcoef.cl[5]/(rcoef.cnl[5] + x*x);

		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void feg_dfeg_Peng_ion_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			using T = Value_type<TrPP_Coef>;
			add_Gauss_dGauss_Fn(x, 0, 5, rcoef, y, dy);
			T yt, t = 1/(rcoef.cnl[5] + x*x);
			y += yt = rcoef.cl[5]*t;
			dy += -2*x*yt*t;
		}

		/***************************************************************************/
		/***************************************************************************/

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void fxg_Doyle_neutral_0_4(const int &Z, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			feg_Doyle_neutral_0_4(x, rcoef, y);
			y = Z - x*x*y;
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void fxg_dfxg_Doyle_neutral_0_4(const int &Z, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			feg_dfeg_Doyle_neutral_0_4(x, rcoef, y, dy);
			dy = -x*(2*y + x*dy);
			y = Z - x*x*y;
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void fxg_Peng_neutral_0_4(const int &Z, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			feg_Peng_neutral_0_4(x, rcoef, y);
			y = Z - x*x*y;
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void fxg_dfxg_Peng_neutral_0_4(const int &Z, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			feg_dfeg_Peng_neutral_0_4(x, rcoef, y, dy);
			dy = -x*(2*y + x*dy);
			y = Z - x*x*y;
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void fxg_Peng_neutral_0_12(const int &Z, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			feg_Peng_neutral_0_12(x, rcoef, y);
			y = Z - x*x*y;
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void fxg_dfxg_Peng_neutral_0_12(const int &Z, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			feg_dfeg_Peng_neutral_0_12(x, rcoef, y, dy);
			dy = -x*(2*y + x*dy);
			y = Z - x*x*y;
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void fxg_Kirkland_neutral_0_12(const int &Z, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			feg_Kirkland_neutral_0_12(x, rcoef, y);
			y = Z - x*x*y;
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void fxg_dfxg_Kirkland_neutral_0_12(const int &Z, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			feg_dfeg_Kirkland_neutral_0_12(x, rcoef, y, dy);
			dy = -x*(2*y + x*dy);
			y = Z - x*x*y;
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void fxg_Weickenmeier_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gauss_Fn(x, 0, 6, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void fxg_dfxg_Weickenmeier_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gauss_dGauss_Fn(x, 0, 6, rcoef, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void fxg_Lobato_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void fxg_dfxg_Lobato_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void fxg_Peng_ion_0_4(const int &Z, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			feg_Peng_neutral_0_4(x, rcoef, y);
			y = Z - x*x*y;
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void fxg_dfxg_Peng_ion_0_4(const int &Z, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			feg_dfeg_Peng_neutral_0_4(x, rcoef, y, dy);
			dy = -x*(2*y + x*dy);
			y = Z - x*x*y;
		}
		/***************************************************************************/
		/***************************************************************************/

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_Doyle_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			Pr_Gauss_feg(x, 0, 4, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_dPr_Doyle_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			Pr_dPr_Gauss_feg(x, 0, 4, rcoef, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_Peng_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			Pr_Gauss_feg(x, 0, 5, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_dPr_Peng_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			Pr_dPr_Gauss_feg(x, 0, 5, rcoef, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_Peng_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gauss_Fn(x, 0, 5, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_dPr_Peng_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			Pr_dPr_Gauss_feg(x, 0, 5, rcoef, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_Kirkland_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Yukawa_Fn(x, 0, 3, rcoef, y);
			Pr_Gauss_feg(x, 3, 6, rcoef, y, false);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_dPr_Kirkland_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Yukawa_dYukawa_Fn(x, 0, 3, rcoef, y, dy);
			Pr_dPr_Gauss_feg(x, 3, 6, rcoef, y, dy, false);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_Weickenmeier_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gauss_Fn(x, 0, 6, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_dPr_Weickenmeier_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gauss_dGauss_Fn(x, 0, 6, rcoef, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_Lobato_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Exponential_Fn(x, 0, 5, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_dPr_Lobato_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Exponential_dExponential_Fn(x, 0, 5, rcoef, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_Peng_ion_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			Pr_Gauss_feg(x, 0, 5, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Pr_dPr_Peng_ion_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			Pr_dPr_Gauss_feg(x, 0, 5, rcoef, y, dy);
		}

		/***************************************************************************/
		/***************************************************************************/

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_Doyle_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gauss_Fn(x, 0, 4, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_dVr_Doyle_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gauss_dGauss_Fn(x, 0, 4, rcoef, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_Peng_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gauss_Fn(x, 0, 5, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_dVr_Peng_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gauss_dGauss_Fn(x, 0, 5, rcoef, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_Peng_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gauss_Fn(x, 0, 5, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_dVr_Peng_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gauss_dGauss_Fn(x, 0, 5, rcoef, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_Kirkland_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Yukawa_Fn(x, 0, 3, rcoef, y);
			add_Gauss_Fn(x, 3, 6, rcoef, y, false);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_dVr_Kirkland_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Yukawa_dYukawa_Fn(x, 0, 3, rcoef, y, dy);
			add_Gauss_dGauss_Fn(x, 3, 6, rcoef, y, dy, false);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_Weickenmeier_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			using T = Value_type<TrPP_Coef>;
			T ix = 1/x;
			y = 0;

			for(auto i = 0; i< 6; i++)
			{
				y += rcoef.cl[i]*erfc(rcoef.cnl[i]*x)*ix;
			}
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_dVr_Weickenmeier_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_Lobato_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			using T = Value_type<TrPP_Coef>;
			T ix = 1/x;
			y = 0;
			for(auto i = 0; i< 5; i++)
			{
				y += rcoef.cl[i]*exp(-rcoef.cnl[i]*x)*ix*(2/rcoef.cnl[i] + x);
			}
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_dVr_Lobato_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_Peng_ion_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gauss_Fn(x, 0, 5, rcoef, y);
			y += rcoef.cl[5]*exp(-rcoef.cnl[5]*x)/x;
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_dVr_Peng_ion_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			using T = Value_type<TrPP_Coef>;
			add_Gauss_dGauss_Fn(x, 0, 5, rcoef, y, dy);
			T yt, ix = 1/x;
			y += yt = rcoef.cl[5]*exp(-rcoef.cnl[5]*x)*ix;
			dy += -(rcoef.cnl[5]+ ix)*yt;
		}

		/***************************************************************************/
		/***************************************************************************/

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void VR_Doyle_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gauss_Fn(x, 0, 4, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void VR_dVR_Doyle_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gauss_dGauss_Fn(x, 0, 4, rcoef, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void VR_Peng_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gauss_Fn(x, 0, 5, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void VR_dVR_Peng_neutral_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gauss_dGauss_Fn(x, 0, 5, rcoef, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void VR_Peng_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			 add_Gauss_Fn(x, 0, 5, rcoef, y);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void VR_dVR_Peng_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gauss_dGauss_Fn(x, 0, 5, rcoef, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void VR_Kirkland_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			y = 0;
			for(auto i = 0; i< 3; i++)
			{
				y += rcoef.cl[i]*bessel_k0(rcoef.cnl[i]*x);
			}

			add_Gauss_Fn(x, 3, 6, rcoef, y, false);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void VR_dVR_Kirkland_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			y = dy = 0;
			for(auto i = 0; i< 3; i++)
			{
				y += rcoef.cl[i]*bessel_k0(rcoef.cnl[i]*x);
				dy += -rcoef.cl[i]*rcoef.cnl[i]*bessel_k1(rcoef.cnl[i]*x);
			}

			add_Gauss_dGauss_Fn(x, 3, 6, rcoef, y, dy, false);
		}

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void VR_Weickenmeier_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y)
		{
			int_VR(Vr_Weickenmeier_neutral_0_12<TrPP_Coef>, x, rcoef, rqz, y);
		}

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void VR_dVR_Weickenmeier_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			int_VR_dVR(Vr_dVr_Weickenmeier_neutral_0_12<TrPP_Coef>, x, rcoef, rqz, y, dy);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void VR_Lobato_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			y = 0;
			for(auto i = 0; i< 5; i++)
			{
				y += rcoef.cl[i]*(2*bessel_k0(rcoef.cnl[i]*x)/rcoef.cnl[i] + x*bessel_k1(rcoef.cnl[i]*x));
			}
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void VR_dVR_Lobato_neutral_0_12(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
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

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void VR_Peng_ion_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y)
		{
			add_Gauss_Fn(x, 0, 5, rcoef, y);
			y += rcoef.cl[5]*bessel_k0(rcoef.cnl[5]*x);
		}

		template <class TrPP_Coef>
		DEVICE_CALLABLE FORCE_INLINE 
		void VR_dVR_Peng_ion_0_4(const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			add_Gauss_dGauss_Fn(x, 0, 5, rcoef, y, dy);
			y += rcoef.cl[5]*bessel_k0(rcoef.cnl[5]*x);
			dy += -rcoef.cl[5]*rcoef.cnl[5]*bessel_k1(rcoef.cnl[5]*x);
		}

		/***************************************************************************/
		/***************************************************************************/

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vz_Doyle_neutral_0_4(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y)
		{
			int_Vz_z0_ze(Vr_Doyle_neutral_0_4<TrPP_Coef>, z0, ze, x, rcoef, rqz, y);
		}

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vz_dVz_Doyle_neutral_0_4(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			int_Vz_dVz_z0_ze(Vr_dVr_Doyle_neutral_0_4<TrPP_Coef>, z0, ze, x, rcoef, rqz, y, dy);
		}

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vz_Peng_neutral_0_4(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y)
		{
			int_Vz_z0_ze(Vr_Peng_neutral_0_4<TrPP_Coef>, z0, ze, x, rcoef, rqz, y);
		}

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vz_dVz_Peng_neutral_0_4(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			int_Vz_dVz_z0_ze(Vr_dVr_Peng_neutral_0_4<TrPP_Coef>, z0, ze, x, rcoef, rqz, y, dy);
		}

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vz_Peng_neutral_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y)
		{
			int_Vz_z0_ze(Vr_Peng_neutral_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y);
		}

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vz_dVz_Peng_neutral_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			int_Vz_dVz_z0_ze(Vr_dVr_Peng_neutral_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y, dy);
		}

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vz_Kirkland_neutral_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y)
		{
			int_Vz_z0_ze(Vr_Kirkland_neutral_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y);
		}

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vz_dVz_Kirkland_neutral_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			int_Vz_dVz_z0_ze(Vr_dVr_Kirkland_neutral_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y, dy);
		}

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vz_Weickenmeier_neutral_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y)
		{
			int_Vz_z0_ze(Vr_Weickenmeier_neutral_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y);
		}

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vz_dVz_Weickenmeier_neutral_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			int_Vz_dVz_z0_ze(Vr_dVr_Weickenmeier_neutral_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y, dy);
		}

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vz_Lobato_neutral_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y)
		{
			int_Vz_z0_ze(Vr_Lobato_neutral_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y);
		}

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vz_dVz_Lobato_neutral_0_12(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			int_Vz_dVz_z0_ze(Vr_dVr_Lobato_neutral_0_12<TrPP_Coef>, z0, ze, x, rcoef, rqz, y, dy);
		}

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vz_Peng_ion_0_4(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y)
		{
			int_Vz_z0_ze(Vr_Peng_ion_0_4<TrPP_Coef>, z0, ze, x, rcoef, rqz, y);
		}

		template <class TrPP_Coef, class TrQ1>
		DEVICE_CALLABLE FORCE_INLINE 
		void Vz_dVz_Peng_ion_0_4(const Value_type<TrPP_Coef> &z0, const Value_type<TrPP_Coef> &ze, const Value_type<TrPP_Coef> &x, const TrPP_Coef &rcoef, const TrQ1 &rqz, Value_type<TrPP_Coef> &y, Value_type<TrPP_Coef> &dy)
		{
			int_Vz_dVz_z0_ze(Vr_dVr_Peng_ion_0_4<TrPP_Coef>, z0, ze, x, rcoef, rqz, y, dy);
		}
		/***************************************************************************/
		/***************************************************************************/

		template <class T> 
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_dVrir_Doyle_neutral_0_4(const T &r, T *cl, T *cnl, T &Vr, T &dVrir)
		{
			T r2 = r*r;

			T Vr0 = cl[0]*exp(-cnl[0]*r2); 
			T Vr1 = cl[1]*exp(-cnl[1]*r2); 
			T Vr2 = cl[2]*exp(-cnl[2]*r2); 
			T Vr3 = cl[3]*exp(-cnl[3]*r2);

			Vr = Vr0 + Vr1 + Vr2 + Vr3;
			dVrir = -2*(cnl[0]*Vr0 + cnl[1]*Vr1 + cnl[2]*Vr2 + cnl[3]*Vr3);
		}

		template <class T> 
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_dVrir_Peng_neutral_0_4_12(const T &r, T *cl, T *cnl, T &Vr, T &dVrir)
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

		template <class T> 
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_dVrir_Kirkland_neutral_0_12(const T &r, T *cl, T *cnl, T &Vr, T &dVrir)
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

		template <class T> 
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_dVrir_Weickenmeier_neutral_0_12(const T &r, T *cl, T *cnl, T &Vr, T &dVrir)
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

		template <class T> 
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_dVrir_Lobato_neutral_0_12(const T &r, T *cl, T *cnl, T &Vr, T &dVrir)
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

		template <class T> 
		DEVICE_CALLABLE FORCE_INLINE 
		void Vr_dVrir_Peng_ion_0_4(const T &r, T *cl, T *cnl, T &Vr, T &dVrir)
		{
			T ir = 1/r; 
			T r2 = r*r;

			T Vr0 = cl[0]*exp(-cnl[0]*r2); 
			T Vr1 = cl[1]*exp(-cnl[1]*r2); 
			T Vr2 = cl[2]*exp(-cnl[2]*r2); 
			T Vr3 = cl[3]*exp(-cnl[3]*r2);
			T Vr4 = cl[4]*exp(-cnl[4]*r2);
			T Vr5 = cl[5]*exp(-cnl[5]*r)*ir;

			Vr = Vr0 + Vr1 + Vr2 + Vr3 + Vr4 + Vr5;
			dVrir = -2*(cnl[0]*Vr0 + cnl[1]*Vr1 + cnl[2]*Vr2 + cnl[3]*Vr3 + cnl[4]*Vr4)-Vr5*(cnl[5]+ir)/r;
		}

		template <class T> 
		DEVICE_CALLABLE FORCE_INLINE 
		int unrolledBinarySearch_c_nR(const T &x0, const T *x)
		{
			int i0 = 0, ie = c_nR-1;
			int im = (i0 + ie)>>1; // divide by 2
			if(x0 < x[im]) ie = im; else i0 = im; // 64
			im = (i0 + ie)>>1; 	// divide by 2
			if(x0 < x[im]) ie = im; else i0 = im; // 32
			im = (i0 + ie)>>1; 	// divide by 2
			if(x0 < x[im]) ie = im; else i0 = im; // 16
			im = (i0 + ie)>>1; 	// divide by 2
			if(x0 < x[im]) ie = im; else i0 = im; // 8
			im = (i0 + ie)>>1; 	// divide by 2
			if(x0 < x[im]) ie = im; else i0 = im; // 4
			im = (i0 + ie)>>1; 	// divide by 2
			if(x0 < x[im]) ie = im; else i0 = im; // 2
			im = (i0 + ie)>>1; 	// divide by 2
			if(x0 < x[im]) ie = im; else i0 = im; // 1
	
			return i0;
		}

		// cosine tapering
		template <class T>
		DEVICE_CALLABLE FORCE_INLINE 
		T tapering(const T &x_tap, const T &alpha, const T &x)
		{
			return (x_tap<x)?cos(alpha*(x-x_tap)):1;
		}

		// cosine tapering
		template <class T>
		DEVICE_CALLABLE FORCE_INLINE 
		void apply_tapering(const T &x_tap, const T &alpha, const T &x, T &y, T &dy)
		{
			if(x_tap<x)
			{
				T tap, dtap;
				sincos(alpha*(x-x_tap), &dtap, &tap);
				dy = dy*tap - alpha*y*dtap;
				y *= tap;
			}
		}

		// Get Local interpolation coefficients
		template <class TAtom>
		DEVICE_CALLABLE FORCE_INLINE 
		void cubic_poly_coef(const int &iR, TAtom &atom)
		{
			auto idR = 1.0/(atom.R2[iR+1]-atom.R2[iR]);
			auto V = atom.c0[iR]; 
			auto Vn = atom.c0[iR+1];
			auto dV = atom.c1[iR]; 
			auto dVn = atom.c1[iR+1];
			auto m = (Vn-V)*idR; 
			auto n = dV+dVn;
			atom.c2[iR] = (3.0*m-n-dV)*idR;
			atom.c3[iR] = (n-2.0*m)*idR*idR;
		}

		// Cubic polynomial evaluation
		template <class T, class TAtom>
		DEVICE_CALLABLE FORCE_INLINE 
		T eval_cubic_poly(const T &R2, const TAtom &atom)
		{
			const int ix = unrolledBinarySearch_c_nR<T>(R2, atom.R2);

			const T dx = R2 - atom.R2[ix]; 
			return (((atom.c3[ix]*dx + atom.c2[ix])*dx + atom.c1[ix])*dx + atom.c0[ix]);
		}

		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		void fft1_shift(const int &ix, const TGrid &grid_1d, TVector &M_io)
		{
			int ix_shift = grid_1d.iRx_shift(ix);
			thrust::swap(M_io[ix], M_io[ix_shift]);
		}

		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		void fft2_sft_bc(const int &ix, const int &iy, const TGrid &grid_2d, TVector &M_io)
		{
			int ixy = grid_2d.ind_col(ix, iy); 
			int ixy_shift = grid_2d.ind_col(ix, grid_2d.nyh+iy);
			thrust::swap(M_io[ixy], M_io[ixy_shift]);
		}

		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		void fft2_shift(const int &ix, const int &iy, const TGrid &grid_2d, TVector &M_io)
		{
			int ixy = grid_2d.ind_col(ix, iy); 
			int ixy_shift = grid_2d.ind_col(grid_2d.nxh+ix, grid_2d.nyh+iy);
			thrust::swap(M_io[ixy], M_io[ixy_shift]);

			ixy = grid_2d.ind_col(ix, grid_2d.nyh+iy); 
			ixy_shift = grid_2d.ind_col(grid_2d.nxh+ix, iy);
			thrust::swap(M_io[ixy], M_io[ixy_shift]);
		}
 
 		/***************************************************************************/
  		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		void assign_shift_2d(const int &ix, const int &iy, const TGrid &grid_2d, 
		TVector &M_i, TVector &M_o)
		{
			int ixy = grid_2d.ind_col(ix, iy); 
			int ixy_shift = grid_2d.ind_col(grid_2d.nxh+ix, grid_2d.nyh+iy);
			M_o[ixy] = M_i[ixy_shift];
			M_o[ixy_shift] = M_i[ixy];

			ixy = grid_2d.ind_col(ix, grid_2d.nyh+iy); 
			ixy_shift = grid_2d.ind_col(grid_2d.nxh+ix, iy);
			M_o[ixy] = M_i[ixy_shift];
			M_o[ixy_shift] = M_i[ixy];
		}
 
 		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		void add_scale_shift_2d(const int &ix, const int &iy, const TGrid &grid_2d, 
		const Value_type<TVector> &w, TVector &M_i, TVector &M_o)
		{
			int ixy = grid_2d.ind_col(ix, iy); 
			int ixy_shift = grid_2d.ind_col(grid_2d.nxh+ix, grid_2d.nyh+iy);

			M_o[ixy] += w*M_i[ixy_shift];
 			M_o[ixy_shift] += w*M_i[ixy];

			/***************************************************************************/
			ixy = grid_2d.ind_col(ix, grid_2d.nyh+iy); 
			ixy_shift = grid_2d.ind_col(grid_2d.nxh+ix, iy);

			M_o[ixy] += w*M_i[ixy_shift];
 			M_o[ixy_shift] += w*M_i[ixy];
		}

		template <class TGrid, class TVector_1, class TVector_2>
		DEVICE_CALLABLE FORCE_INLINE 
		void add_scale_square_shift_2d(const int &ix, const int &iy, const TGrid &grid_2d, 
		const Value_type<TVector_2> &w, TVector_1 &M_i, TVector_2 &M_o)
		{
			int ixy = grid_2d.ind_col(ix, iy); 
			int ixy_shift = grid_2d.ind_col(grid_2d.nxh+ix, grid_2d.nyh+iy);

			M_o[ixy] += w*::norm(M_i[ixy_shift]);
 			M_o[ixy_shift] += w*::norm(M_i[ixy]);

			/***************************************************************************/
			ixy = grid_2d.ind_col(ix, grid_2d.nyh+iy); 
			ixy_shift = grid_2d.ind_col(grid_2d.nxh+ix, iy);

			M_o[ixy] += w*::norm(M_i[ixy_shift]);
 			M_o[ixy_shift] += w*::norm(M_i[ixy]);
		}

		/***************************************************************************/
 		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		void assign_crop(const int &ix, const int &iy, const TGrid &grid_2d, TVector &M_i, Range_2d &range, TVector &M_o)
		{
			if(range.chk_bound(ix, iy))
			{
				 M_o[range.ind_col_o(ix, iy)] = M_i[grid_2d.ind_col(ix, iy)];
			}
		}

 		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		void assign_crop_shift_2d(const int &ix, const int &iy, const TGrid &grid_2d, TVector &M_i, Range_2d &range, TVector &M_o)
		{
 			int ix_i = ix;
			int iy_i = iy;

			int ix_s = grid_2d.nxh+ix;
			int iy_s = grid_2d.nyh+iy;

			if(range.chk_bound(ix_i, iy_i))
			{
				 M_o[range.ind_col_o(ix_i, iy_i)] = M_i[grid_2d.ind_col(ix_s, iy_s)];
			}

 			if(range.chk_bound(ix_s, iy_s))
			{
				 M_o[range.ind_col_o(ix_s, iy_s)] = M_i[grid_2d.ind_col(ix_i, iy_i)];
			}

			/***************************************************************************/
 			ix_i = ix;
			iy_i = grid_2d.nyh+iy;

			ix_s = grid_2d.nxh+ix;
			iy_s = iy;

			if(range.chk_bound(ix_i, iy_i))
			{
				 M_o[range.ind_col_o(ix_i, iy_i)] = M_i[grid_2d.ind_col(ix_s, iy_s)];
			}

 			if(range.chk_bound(ix_s, iy_s))
			{
				 M_o[range.ind_col_o(ix_s, iy_s)] = M_i[grid_2d.ind_col(ix_i, iy_i)];
			}
		}

 		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		void add_scale_crop_shift_2d(const int &ix, const int &iy, const TGrid &grid_2d, 
		const Value_type<TVector> &w, TVector &M_i, Range_2d &range, TVector &M_o)
		{
 			int ix_i = ix;
			int iy_i = iy;

			int ix_s = grid_2d.nxh+ix;
			int iy_s = grid_2d.nyh+iy;

			if(range.chk_bound(ix_i, iy_i))
			{
				 M_o[range.ind_col_o(ix_i, iy_i)] += w*M_i[grid_2d.ind_col(ix_s, iy_s)];
			}

 			if(range.chk_bound(ix_s, iy_s))
			{
				 M_o[range.ind_col_o(ix_s, iy_s)] += w*M_i[grid_2d.ind_col(ix_i, iy_i)];
			}

			/***************************************************************************/
 			ix_i = ix;
			iy_i = grid_2d.nyh+iy;

			ix_s = grid_2d.nxh+ix;
			iy_s = iy;

			if(range.chk_bound(ix_i, iy_i))
			{
				 M_o[range.ind_col_o(ix_i, iy_i)] += w*M_i[grid_2d.ind_col(ix_s, iy_s)];
			}

 			if(range.chk_bound(ix_s, iy_s))
			{
				 M_o[range.ind_col_o(ix_s, iy_s)] += w*M_i[grid_2d.ind_col(ix_i, iy_i)];
			}
		}

		template <class TGrid, class TVector_1, class TVector_2>
		DEVICE_CALLABLE FORCE_INLINE 
		void add_scale_square_crop_shift_2d(const int &ix, const int &iy, const TGrid &grid_2d, 
		const Value_type<TVector_2> &w, TVector_1 &M_i, Range_2d &range, TVector_2 &M_o)
		{
 			int ix_i = ix;
			int iy_i = iy;

			int ix_s = grid_2d.nxh+ix;
			int iy_s = grid_2d.nyh+iy;

			if(range.chk_bound(ix_i, iy_i))
			{
				 M_o[range.ind_col_o(ix_i, iy_i)] += w*::norm(M_i[grid_2d.ind_col(ix_s, iy_s)]);
			}

 			if(range.chk_bound(ix_s, iy_s))
			{
				 M_o[range.ind_col_o(ix_s, iy_s)] += w*::norm(M_i[grid_2d.ind_col(ix_i, iy_i)]);
			}

			/***************************************************************************/
 			ix_i = ix;
			iy_i = grid_2d.nyh+iy;

			ix_s = grid_2d.nxh+ix;
			iy_s = iy;

			if(range.chk_bound(ix_i, iy_i))
			{
				 M_o[range.ind_col_o(ix_i, iy_i)] += w*::norm(M_i[grid_2d.ind_col(ix_s, iy_s)]);
			}

 			if(range.chk_bound(ix_s, iy_s))
			{
				 M_o[range.ind_col_o(ix_s, iy_s)] += w*::norm(M_i[grid_2d.ind_col(ix_i, iy_i)]);
			}
		}

		/***************************************************************************/
		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		void sum_over_Det(const int &ix, const int &iy, const TGrid &grid_2d, 
		const Value_type<TGrid> &g2_min, const Value_type<TGrid> &g2_max, const TVector &M_i, Value_type<TGrid> &sum)
		{
			auto g2 = grid_2d.g2_shift(ix, iy);
			if((g2_min <= g2) && (g2 < g2_max))
			{
				int ixy = grid_2d.ind_col(ix, iy); 
				sum += M_i[ixy];
			}
		}

		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		void sum_square_over_Det(const int &ix, const int &iy, const TGrid &grid_2d, 
		const Value_type<TGrid> &g2_min, const Value_type<TGrid> &g2_max, const TVector &M_i, Value_type<TGrid> &sum)
		{
			auto g2 = grid_2d.g2_shift(ix, iy);
			if((g2_min <= g2) && (g2 < g2_max))
			{
				int ixy = grid_2d.ind_col(ix, iy);
				sum += norm(M_i[ixy]);
			}
		}

		template <class TGrid, class TVector_1, class TVector_2>
		DEVICE_CALLABLE FORCE_INLINE 
		void sum_square_over_Det(const int &ix, const int &iy, const TGrid &grid_2d, 
		const TVector_1 &S_i, const TVector_2 &M_i, Value_type<TGrid> &sum)
		{
			const int ixy = grid_2d.ind_col(ix, iy);
			sum += S_i[ixy] * norm(M_i[ixy]);
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void bandwidth_limit(const int &ix, const int &iy, const TGrid &grid_2d, TVector_c &M_io)
		{
			const int ixy = grid_2d.ind_col(ix, iy); 
			M_io[ixy] *= grid_2d.bwl_factor_shift(ix, iy)/grid_2d.nxy_r();
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void hard_aperture(const int &ix, const int &iy, const TGrid &grid_2d, const Value_type<TGrid> &g2_max, const Value_type<TGrid> &w, TVector_c &M_io)
		{
			using T = Value_type<TVector_c>;
			const int ixy = grid_2d.ind_col(ix, iy); 
			const auto g2 = grid_2d.g2_shift(ix, iy);

			M_io[ixy] = ((g2 < g2_max))?(T(w)*M_io[ixy]):0;
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void propagate(const int &ix, const int &iy, const TGrid &grid_2d, const Value_type<TGrid> &w, 
		const Value_type<TGrid> &gx_0, const Value_type<TGrid> &gy_0, TVector_c &psi_i, TVector_c &psi_o)
		{
			const int ixy = grid_2d.ind_col(ix, iy);
			const auto m = grid_2d.bwl_factor_shift(ix, iy)/grid_2d.nxy_r();
			const auto theta = w*grid_2d.g2_shift(ix, iy, gx_0, gy_0);

			psi_o[ixy] = polar(m, theta)*psi_i[ixy];
		}

		/********************* phase shifts real space **********************/
		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void exp_r_factor_1d(const int &ix, const TGrid &grid_1d, 
		const Value_type<TGrid> gx, TVector_c &psi_i, TVector_c &psi_o)
		{
			const auto Rx = grid_1d.Rx_shift(ix)-grid_1d.Rx_c();
			psi_o[ix] = psi_i[ix]*euler(gx*Rx)/grid_1d.nx_r();
		}

		template <class TGrid, class TVector_r, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void exp_r_factor_2d_bc(const int &ix, const int &iy, const TGrid &grid_2d, 
		const Value_type<TGrid> alpha, TVector_r &gy, TVector_c &psi_i, TVector_c &psi_o)
		{
			const int ixy = grid_2d.ind_col(ix, iy);
			const auto Ry = grid_2d.Ry_shift(iy)-grid_2d.Ry_c();
			psi_o[ixy] = psi_i[ixy]*euler(alpha*gy[ix]*Ry)/grid_2d.ny_r();
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void exp_r_factor_2d(const int &ix, const int &iy, const TGrid &grid_2d, 
		const Value_type<TGrid> &gx, const Value_type<TGrid> &gy, TVector_c &psi_i, TVector_c &psi_o)
		{
			const int ixy = grid_2d.ind_col(ix, iy);
			const auto Rx = grid_2d.Rx_shift(ix)-grid_2d.Rx_c();
			const auto Ry = grid_2d.Ry_shift(iy)-grid_2d.Ry_c();
			psi_o[ixy] = psi_i[ixy]*euler(gx*Rx + gy*Ry)/grid_2d.nxy_r();
		}

		template <class TGrid, class TVector_r, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void mul_exp_r_factor_2d(const int &ix, const int &iy, const TGrid &grid_2d, 
		TVector_r &gx, TVector_r &gy, TVector_c &psi_i, TVector_c &psi_o)
		{
			using T_c = Value_type<TVector_c>;

			const int ixy = grid_2d.ind_col(ix, iy);
			const auto Rx = grid_2d.Rx_shift(ix)-grid_2d.Rx_c();
			const auto Ry = grid_2d.Ry_shift(iy)-grid_2d.Ry_c();

			T_c exp_sup = 0;
			for(auto it=0; it<gx.size(); it++)
			{
				exp_sup += euler(gx[it]*Rx + gy[it]*Ry);
			}
			psi_o[ixy] = psi_i[ixy]*exp_sup/grid_2d.nxy_r();
		}

		/********************* phase shifts Fourier space **********************/
		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void exp_g_factor_1d(const int &ix, const TGrid &grid_1d, 
		const Value_type<TGrid> Rx, TVector_c &psi_i, TVector_c &psi_o)
		{
			const auto gx = grid_1d.gx_shift(ix);
			psi_o[ix] = psi_i[ix]*euler(Rx*gx)/grid_1d.nx_r();
		}

		template <class TGrid, class TVector_r, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void exp_g_factor_2d_bc(const int &ix, const int &iy, const TGrid &grid_2d, 
		const Value_type<TGrid> alpha, TVector_r &Ry, TVector_c &psi_i, TVector_c &psi_o)
		{
			const int ixy = grid_2d.ind_col(ix, iy);
			const auto gy = grid_2d.gy_shift(iy);
			psi_o[ixy] = psi_i[ixy]*euler(alpha*Ry[ix]*gy)/grid_2d.ny_r();
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void exp_g_factor_2d(const int &ix, const int &iy, const TGrid &grid_2d, 
		const Value_type<TGrid> &Rx, const Value_type<TGrid> &Ry, TVector_c &psi_i, TVector_c &psi_o)
		{
			const int ixy = grid_2d.ind_col(ix, iy);
			const auto gx = grid_2d.gx_shift(ix);
			const auto gy = grid_2d.gy_shift(iy);
			psi_o[ixy] = psi_i[ixy]*euler(Rx*gx + Ry*gy)/grid_2d.nxy_r();
		}

		template <class TGrid, class TVector_r, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void mul_exp_g_factor_2d(const int &ix, const int &iy, const TGrid &grid_2d, 
		TVector_r &Rx, TVector_r &Ry, TVector_c &psi_i, TVector_c &psi_o)
		{
			using T_c = Value_type<TVector_c>;

			const int ixy = grid_2d.ind_col(ix, iy);
			const auto gx = grid_2d.gx_shift(ix);
			const auto gy = grid_2d.gy_shift(iy);

			T_c exp_sup = 0;
			for(auto it=0; it<Rx.size(); it++)
			{
				exp_sup += euler(Rx[it]*gx + Ry[it]*gy);
			}
			psi_o[ixy] = psi_i[ixy]*exp_sup/grid_2d.nxy_r();
		}

		template <class T>
		DEVICE_CALLABLE FORCE_INLINE 
		complex<T> exp_i_chi(const int &ix, const int &iy, const Grid_2d<T> &grid_2d, const Lens<T> &lens, 
		const T &x, const T &y, const T &gxu, const T &gyu)
		{
			auto gx = grid_2d.gx_shift(ix)+gxu;
			auto gy = grid_2d.gy_shift(iy)+gyu;
			auto g2 = gx*gx + gy*gy;

			complex<T> v = 0;

			if((lens.g2_min <= g2) && (g2 < lens.g2_max))
			{
				auto g4 = g2*g2;
				auto g6 = g4*g2;
				auto chi = x*gx + y*gy + lens.eval_c_10(g2) + lens.eval_c_30(g4) + lens.eval_c_50(g6);
				if(lens.is_phi_required())
				{
					auto g = sqrt(g2);
					auto g3 = g2*g;
					auto g5 = g4*g;
					auto phi = atan2(gy, gx);
					chi += lens.eval_m(phi) + lens.eval_c_12(g2, phi);
					chi += lens.eval_c_21_c_23(g3, phi) + lens.eval_c_32_c_34(g4, phi);
					chi += lens.eval_c_41_c_43_c_45(g5, phi) + lens.eval_c_52_c_54_c_56(g6, phi); 
				}	
				v = euler(chi); 
			}

			return v;
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void probe(const int &ix, const int &iy, const TGrid &grid_2d, const Lens<Value_type<TGrid>> &lens, 
		const Value_type<TGrid> &x, const Value_type<TGrid> &y, const Value_type<TGrid> &gxu, 
		const Value_type<TGrid> &gyu, TVector_c &fPsi_o)
		{
			auto v = exp_i_chi(ix, iy, grid_2d, lens, x, y, gxu, gyu);

			int ixy = grid_2d.ind_col(ix, iy);
			fPsi_o[ixy] = v;
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void apply_CTF(const int &ix, const int &iy, const TGrid &grid_2d, const Lens<Value_type<TGrid>> &lens, 
		const Value_type<TGrid> &gxu, const Value_type<TGrid> &gyu, TVector_c &fPsi_i, TVector_c &fPsi_o)
		{
			using T = Value_type<TGrid>;

			T x = 0;
			T y = 0;
			auto v = exp_i_chi(ix, iy, grid_2d, lens, x, y, gxu, gyu);

			int ixy = grid_2d.ind_col(ix, iy);
			fPsi_o[ixy] = v*fPsi_i[ixy];
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void apply_PCTF(const int &ix, const int &iy, const TGrid &grid_2d, const Lens<Value_type<TGrid>> &lens, 
		TVector_c &fPsi_i, TVector_c &fPsi_o)
		{
			using T_r = Value_type<TGrid>;

			const T_r c_Pi = 3.141592653589793238463;

			int ixy = grid_2d.ind_col(ix, iy);
			T_r g2 = grid_2d.g2_shift(ix, iy);

			if((lens.g2_min <= g2) && (g2 < lens.g2_max))
			{			
				T_r chi = g2*(lens.c_c_30*g2 + lens.c_c_10);
				T_r c = c_Pi*lens.si_theta_c*lens.ti_iehwgd;
				T_r u = 1.0 + c*c*g2;

				c = c_Pi*lens.ti_iehwgd*lens.lambda*g2;
				T_r temp_inc = 0.25*c*c;

				c = c_Pi*lens.si_theta_c*(lens.c_30*lens.lambda2*g2 + lens.c_10);
				T_r spa_inc = c*c*g2;

				T_r st_inc = exp(-(spa_inc+temp_inc)/u)/sqrt(u);

				fPsi_o[ixy] = fPsi_i[ixy]*polar(st_inc, chi);
			}
			else
			{
				fPsi_o[ixy] = 0;
			}
		}

		template <class TGrid>
		DEVICE_CALLABLE FORCE_INLINE 
		void Lorentz_factor(const int &ix, const int &iy, const TGrid &grid_2d, const Value_type<TGrid> &gc2, const Value_type<TGrid> &ge2, Value_type<TGrid> &sum)
		{
			using T_r = Value_type<TGrid>;

			T_r g2 = grid_2d.g2_shift(ix, iy);
			if(g2 < gc2)
			{
				sum += 1.0/(g2 + ge2);
			}
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void kernel_xyz(const int &ix, const int &iy, const TGrid &grid_2d, const EELS<Value_type<TGrid>> &eels, TVector_c &k_x, TVector_c &k_y, TVector_c &k_z)
		{
			using T_r = Value_type<TGrid>;
			using T_c = Value_type<TVector_c>;

			int ixy = grid_2d.ind_col(ix, iy);
			T_r gx = grid_2d.gx_shift(ix);
			T_r gy = grid_2d.gy_shift(iy);
			T_r g2 = gx*gx + gy*gy;
				
			if(g2 < eels.gc2)
			{
				T_c pos = euler(eels.x*gx + eels.y*gy);
				T_r lorentz = eels.factor/(g2 + eels.ge2);
				k_x[ixy] = T_c(gx*lorentz, 0)*pos;
				k_y[ixy] = T_c(gy*lorentz, 0)*pos;
				k_z[ixy] = T_c(eels.ge*lorentz, 0)*pos;
			}
			else
			{
				k_x[ixy] = 0;
				k_y[ixy] = 0;
				k_z[ixy] = 0;
			}
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void kernel_x(const int &ix, const int &iy, const TGrid &grid_2d, const EELS<Value_type<TGrid>> &eels, TVector_c &k_x)
		{
			using T_r = Value_type<TGrid>;
			using T_c = Value_type<TVector_c>;

			int ixy = grid_2d.ind_col(ix, iy);
			T_r gx = grid_2d.gx_shift(ix);
			T_r gy = grid_2d.gy_shift(iy);
			T_r g2 = gx*gx + gy*gy;
				
			if(g2 < eels.gc2)
			{
				T_c pos = euler(eels.x*gx + eels.y*gy);
				T_r lorentz = eels.factor/(g2 + eels.ge2);
				k_x[ixy] = T_c(gx*lorentz, 0)*pos;
			}
			else
			{
				k_x[ixy] = 0;
			}
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void kernel_y(const int &ix, const int &iy, const TGrid &grid_2d, const EELS<Value_type<TGrid>> &eels, TVector_c &k_y)
		{
			using T_r = Value_type<TGrid>;
			using T_c = Value_type<TVector_c>;

			int ixy = grid_2d.ind_col(ix, iy);
			T_r gx = grid_2d.gx_shift(ix);
			T_r gy = grid_2d.gy_shift(iy);
			T_r g2 = gx*gx + gy*gy;
				
			if(g2 < eels.gc2)
			{
				T_c pos = euler(eels.x*gx + eels.y*gy);
				T_r lorentz = eels.factor/(g2 + eels.ge2);
				k_y[ixy] = T_c(gy*lorentz, 0)*pos;
			}
			else
			{
				k_y[ixy] = 0;
			}
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void kernel_z(const int &ix, const int &iy, const TGrid &grid_2d, const EELS<Value_type<TGrid>> &eels, TVector_c &k_z)
		{
			using T_r = Value_type<TGrid>;
			using T_c = Value_type<TVector_c>;

			int ixy = grid_2d.ind_col(ix, iy);
			T_r gx = grid_2d.gx_shift(ix);
			T_r gy = grid_2d.gy_shift(iy);
			T_r g2 = gx*gx + gy*gy;
				
			if(g2 < eels.gc2)
			{
				T_c pos = euler(eels.x*gx + eels.y*gy);
				T_r lorentz = eels.factor/(g2 + eels.ge2);
				k_z[ixy] = T_c(eels.ge*lorentz, 0)*pos;
			}
			else
			{
				k_z[ixy] = 0;
			}
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void kernel_mn1(const int &ix, const int &iy, const TGrid &grid_2d, const EELS<Value_type<TGrid>> &eels, TVector_c &k_mn1)
		{
			using T_r = Value_type<TGrid>;
			using T_c = Value_type<TVector_c>;

			int ixy = grid_2d.ind_col(ix, iy);
			T_r gx = grid_2d.gx_shift(ix);
			T_r gy = grid_2d.gy_shift(iy);
			T_r g2 = gx*gx + gy*gy;
				
			if(g2 < eels.gc2)
			{
				const T_r c_i2i2 = 0.70710678118654746; 

				T_c pos = euler(eels.x*gx + eels.y*gy);
				T_r lorentz = c_i2i2*eels.factor/(g2 + eels.ge2);
				T_c k_x(gx*lorentz, 0);
				T_c k_y(0, gy*lorentz);
				k_mn1[ixy] = (k_x - k_y)*pos;
			}
			else
			{
				k_mn1[ixy] = 0;
			}
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void kernel_mp1(const int &ix, const int &iy, const TGrid &grid_2d, const EELS<Value_type<TGrid>> &eels, TVector_c &k_mp1)
		{
			using T_r = Value_type<TGrid>;
			using T_c = Value_type<TVector_c>;

			int ixy = grid_2d.ind_col(ix, iy);
			T_r gx = grid_2d.gx_shift(ix);
			T_r gy = grid_2d.gy_shift(iy);
			T_r g2 = gx*gx + gy*gy;
				
			if(g2 < eels.gc2)
			{
				const T_r c_i2i2 = 0.70710678118654746; 

				T_c pos = euler(eels.x*gx + eels.y*gy);
				T_r lorentz = c_i2i2*eels.factor/(g2 + eels.ge2);
				T_c k_x(gx*lorentz, 0);
				T_c k_y(0, gy*lorentz);
				k_mp1[ixy] = (k_x + k_y)*pos;
			}
			else
			{
				k_mp1[ixy] = 0;
			}
		}
		
		template <class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		void trs(const int &ix, const int &iy, const int &ncols, const int &nrows, TVector &M_i, TVector &M_o)
		{
			int ixy = ix*nrows+iy;
			int ixy_t = iy*ncols+ix;
			M_o[ixy_t] = M_i[ixy];
		}

		/****************** Gaussian convolution ********************/
		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void gauss_cv_1d(const int &ix, const TGrid &grid_1d, 
		const Value_type<TGrid> &alpha, TVector_c &M_io)
		{
			auto fg = exp(-alpha*grid_1d.g2_shift(ix));
			M_io[ix] *= fg/grid_1d.nx_r();
		};

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void gauss_cv_2d(const int &ix, const int &iy, const TGrid &grid_2d, 
		const Value_type<TGrid> &alpha, TVector_c &M_io)
		{
			auto fg = exp(-alpha*grid_2d.g2_shift(ix, iy));
			M_io[grid_2d.ind_col(ix, iy)] *= fg/grid_2d.nxy_r();
		};

		/****************** Gaussian deconvolution ********************/
		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void gauss_dcv_1d(const int &ix, const TGrid &grid_1d, 
		const Value_type<TGrid> &alpha, const Value_type<TGrid> &PSNR, TVector_c &M_io)
		{
			auto fg = exp(-alpha*grid_1d.g2_shift(ix));
			M_io[ix] *= fg/((fg*fg+PSNR)*grid_1d.nx_r());
		};

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void gauss_dcv_2d(const int &ix, const int &iy, const TGrid &grid_2d, 
		const Value_type<TGrid> &alpha, const Value_type<TGrid> &PSNR, TVector_c &M_io)
		{
			auto fg = exp(-alpha*grid_2d.g2_shift(ix, iy));
			M_io[grid_2d.ind_col(ix, iy)] *= fg/((fg*fg+PSNR)*grid_2d.nxy_r());
		};

		/***************** vector col/row x matrix *****************/
		template <class TGrid, class TVector, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void vector_row_x_matrix(const int &ix, const int &iy, const TGrid &grid_2d, 
		TVector &Vr, TVector_c &M_io)
		{
			M_io[grid_2d.ind_row(ix, iy)] *= Vr[ix];
		}

		template <class TGrid, class TVector, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void vector_col_x_matrix(const int &ix, const int &iy, const TGrid &grid_2d, 
		TVector &Vc, TVector_c &M_io)
		{
			M_io[grid_2d.ind_col(ix, iy)] *= Vc[iy];
		}

		/*************** phase correlation functions **************/
		template <class TGrid, class TVector, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void pcf_1d_pp(const int &ix, const int &ix_s, const TGrid &grid_1d, 
		const Butterworth_1d<Value_type<TGrid>> &bw_1d, TVector &M_i, TVector_c &M_o)
		{
			using T = Value_type<TGrid>;

			int ix_n = (ix+1<grid_1d.nx)?(ix+1):ix;
			int ix_o = ix + ix_s;

			T R2 = grid_1d.R2(ix, bw_1d.x_c);
			T fh = bw_1d.eval_norm(R2);
			M_o[ix_o] = (M_i[ix_n]-M_i[ix])*fh;
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void pcf_1d_gaussian(const int &ix, const TGrid &grid_1d, const Gauss_1d<Value_type<TGrid>> &gs_1d, 
		TVector_c &M_1, TVector_c &M_2, TVector_c &pcf)
		{
			using T = Value_type<TGrid>;

			auto z = conj(M_1[ix])*M_2[ix];
			T g2 = grid_1d.g2_shift(ix);
			T m = gs_1d(g2);
			T theta = thrust::arg<T>(z);
			pcf[ix] = (ix == 0)?1:thrust::polar(m, theta);
		}

		template <class TGrid, class TVector, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void pcf_2d_bc_pp(const int &ix, const int &iy, const int &iy_s, const TGrid &grid_2d_i, 
		const TGrid &grid_2d_o, TVector &M_i, TVector &fh, TVector_c &M_o)
		{
			int ixy_i = grid_2d_i.ind_col(ix, iy);
			int iy_n = (iy+1<grid_2d_i.ny)?(iy+1):iy;
			int ixy_in = grid_2d_i.ind_col(ix, iy_n);
			int ixy_o = grid_2d_o.ind_col(ix, iy+iy_s);
			M_o[ixy_o] = (M_i[ixy_in]-M_i[ixy_i])*fh[iy];
		}

		template <class TGrid, class TVector, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void pcf_2d_bc_gaussian(const int &ix, const int &iy, const TGrid &grid_2d, 
		TVector_c &M_1, TVector_c &M_2, TVector &fg, TVector_c &pcf)
		{
			using T = Value_type<TGrid>;

			int ixy = grid_2d.ind_col(ix, iy);
			auto z = conj(M_1[ixy])*M_2[ixy];
			T m = fg[iy];
			T theta = thrust::arg<T>(z);
			pcf[ixy] = (iy == 0)?1:thrust::polar(m, theta);
		}

		template <class TGrid, class TVector, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void pcf_2d_pp(const int &ix, const int &iy, const TGrid &grid_2d, 
		const Butterworth_2d<Value_type<TGrid>> &bw_2d, TVector &M_i, TVector_c &M_o)
		{
			using T = Value_type<TGrid>;

			int ixy = grid_2d.ind_col(ix, iy);
			int iy_n = (iy+1<grid_2d.ny)?(iy+1):iy;
			int ixy_n = grid_2d.ind_col(ix, iy_n);

			T R2 = grid_2d.R2(ix, iy, bw_2d.x_c, bw_2d.y_c);
			T fh = bw_2d.eval_norm(R2);
			M_o[ixy] = (M_i[ixy_n]-M_i[ixy])*fh;
		}

		template <class TGrid, class TVector_c>
		DEVICE_CALLABLE FORCE_INLINE 
		void pcf_2d_gaussian(const int &ix, const int &iy, const TGrid &grid_2d, 
		const Gauss_2d<Value_type<TGrid>> &gs_2d, TVector_c &M_1, TVector_c &M_2, TVector_c &pcf)
		{
			using T = Value_type<TGrid>;

			int ixy = grid_2d.ind_col(ix, iy);
			auto z = conj(M_1[ixy])*M_2[ixy];
			T g2 = grid_2d.g2_shift(ix, iy);
			T m = gs_2d(g2);
			T theta = thrust::arg<T>(z);
			pcf[ixy] = (ixy == 0)?1:thrust::polar(m, theta);
		}

		template <class T>
		DEVICE_CALLABLE FORCE_INLINE 
		r2d<T> af_iscale(const r2d<T> &p, const T &fxy)
		{
			return p/fxy;
		}

		template <class T>
		DEVICE_CALLABLE FORCE_INLINE 
		r2d<T> af_irotate(const T &theta, const r2d<T> &p0, r2d<T> p)
		{
			T sin_t, cos_t;
			sincos(theta, &sin_t, &cos_t);
			p -= p0;
			p = r2d<T>(cos_t*p.x+sin_t*p.y, -sin_t*p.x+cos_t*p.y);
			p += p0;
			return p;
		}

		template <class T>
		DEVICE_CALLABLE FORCE_INLINE 
		r2d<T> af_irot_sca_sft(const T &theta, const r2d<T> &p0, const T &fx, const T &fy, const r2d<T> &ps, r2d<T> p)
		{
			T sin_t, cos_t;
			sincos(theta, &sin_t, &cos_t);
			p.x = (p.x-ps.x)/fx;
 p.y = (p.y-ps.y)/fy;
			p -= p0;
			p = r2d<T>(cos_t*p.x+sin_t*p.y, -sin_t*p.x+cos_t*p.y);
			p += p0;
			return p;
		}

		template <class T>
		DEVICE_CALLABLE FORCE_INLINE 
		r2d<T> af_shx_scy(const r2d<T> &p, const T &a, const T &b)
		{
			return r2d<T>(p.x+a*p.y, b*p.y);
		}

		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		Value_type<TGrid> interp_bl_2d(const r2d<Value_type<TGrid>> &p, const TGrid &grid_2d, TVector &Im)
		{
			using T = Value_type<TGrid>;

			const int ix = min(grid_2d.nx-2, max(0, grid_2d.floor_dRx(p.x)));
			const int iy = min(grid_2d.ny-2, max(0, grid_2d.floor_dRy(p.y)));

			const T f11 = Im[grid_2d.ind_col(ix, iy)];
			const T f12 = Im[grid_2d.ind_col(ix, iy+1)];
			const T f21 = Im[grid_2d.ind_col(ix+1, iy)];
			const T f22 = Im[grid_2d.ind_col(ix+1, iy+1)];

			const T x1 = grid_2d.Rx(ix);
			const T x2 = grid_2d.Rx(ix+1);
			const T y1 = grid_2d.Ry(iy);
			const T y2 = grid_2d.Ry(iy+1);
	
			const T dx1 = (p.x-x1)/(x2-x1);
			const T dx2 = (x2-p.x)/(x2-x1);
			const T dy1 = (p.y-y1)/(y2-y1);
			const T dy2 = (y2-p.y)/(y2-y1);
			T f = dx2*(f11*dy2 + f12*dy1)+dx1*(f21*dy2 + f22*dy1);

			return f;
		};

		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		void sc_2d(const int &ix, const int &iy, const TGrid &grid_2d_i, TVector &M_i, 
		const Value_type<TGrid> &fxy, const TGrid &grid_2d_o, TVector &M_o)
		{
			using T = Value_type<TGrid>;

			r2d<T> p(grid_2d_o.Rx(ix), grid_2d_o.Ry(iy));
			p = af_iscale(p, fxy);

			M_o[grid_2d_o.ind_col(ix, iy)] = interp_bl_2d(p, grid_2d_i, M_i);
		};

		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		void rot_2d(const int &ix, const int &iy, const TGrid &grid_2d_i, TVector &M_i, 
		const Value_type<TGrid> &theta, const r2d<Value_type<TGrid>> &p0, const Value_type<TGrid> &bg, 
		const TGrid &grid_2d_o, TVector &M_o)
		{
			using T = Value_type<TGrid>;

			r2d<T> p(grid_2d_o.Rx(ix), grid_2d_o.Ry(iy));
			p = af_irotate(theta, p0, p);

			M_o[grid_2d_o.ind_col(ix, iy)] = (grid_2d_i.ckb_bound(p))?interp_bl_2d(p, grid_2d_i, M_i):bg;
		};

		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		void rot_sca_sft_2d(const int &ix, const int &iy, const TGrid &grid_2d_i, TVector &M_i, 
		const Value_type<TGrid> &theta, const r2d<Value_type<TGrid>> &p0, const Value_type<TGrid> &fx, const Value_type<TGrid> &fy, 
		const r2d<Value_type<TGrid>> &ps, const Value_type<TGrid> &bg, const TGrid &grid_2d_o, TVector &M_o)
		{
			using T = Value_type<TGrid>;

			r2d<T> p(grid_2d_o.Rx(ix), grid_2d_o.Ry(iy));
			p = af_irot_sca_sft(theta, p0, fx, fy, ps, p);

			M_o[grid_2d_o.ind_col(ix, iy)] = (grid_2d_i.ckb_bound(p))?interp_bl_2d(p, grid_2d_i, M_i):bg;
		};

		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		void shx_scy(const int &ix, const int &iy, const TGrid &grid_2d_i, TVector &M_i, 
		const Value_type<TGrid> &fx, const Value_type<TGrid> &fy, const Value_type<TGrid> &bg, 
		const TGrid &grid_2d_o, TVector &M_o)
		{
			using T = Value_type<TGrid>;

			r2d<T> p(grid_2d_o.Rx(ix), grid_2d_o.Ry(iy));
			p = af_shx_scy(p, fx, fy);

			M_o[grid_2d_o.ind_col(ix, iy)] = (grid_2d_i.ckb_bound(p))?interp_bl_2d(p, grid_2d_i, M_i):bg;
		};

		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		Value_type<TGrid> max_pos_1d(const TGrid &grid_1d, TVector &M)
		{
			int ix_max = thrust::max_element(M.begin(), M.end())-M.begin();

			return grid_1d.Rx(ix_max);
		};

		template <class TGrid, class TVector>
		DEVICE_CALLABLE FORCE_INLINE 
		r2d<Value_type<TGrid>> max_pos_2d(const TGrid &grid_2d, TVector &M)
		{
			using T = Value_type<TGrid>;

			// get maximum index position
			int ixy_max = thrust::max_element(M.begin(), M.end())-M.begin();
			int ix_max, iy_max;
			grid_2d.col_row(ixy_max, ix_max, iy_max);

			return r2d<T>(grid_2d.Rx(ix_max), grid_2d.Ry(iy_max));
		};

		// template <class TGrid, class TVector>
		// DEVICE_CALLABLE FORCE_INLINE 
		// void dilate(const int &ix_i, const int &iy_i, TVector &Im_i, TVector &Im_o)
		// {
		// 	int ix_0 = max(ix_i+nk0, 0);
		// 	int ixe = min(ix_i+nke, nx_i);

		// 	int iy_0 = max(iy_i+nk0, 0);
		// 	int iye = min(iy_i+nke, ny_i);

		// 	for (auto ix = ix_0; ix < ixe; ix++)
		// 	{
		// 		for (auto iy = iy_0; iy < iye; iy++)
		// 		{
		// 			if(Im_i[ix*ny_i+iy]>0.5)
		// 			{	
		// 				Im_o[ix_i*ny_i+iy_i] = 1;
		// 				return;
		// 			}
		// 		}
		// 	}
		// 	Im_o[ix_i*ny_i+iy_i] = 0;
		// }
	} // host_device_detail

	// Electron scattering factors calculation (feg)
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	void feg(const ePotential_Type &potential_type, const int &charge, const T &g, const rPP_Coef<T> &c_feg, T &y)
	{
		if(charge == 0)
		{
			switch(potential_type)
			{
				case ePT_Doyle_0_4:
					// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
					host_device_detail::feg_Doyle_neutral_0_4(g, c_feg, y);
					break;
				case ePT_Peng_0_4:
					// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
					host_device_detail::feg_Peng_neutral_0_4(g, c_feg, y);
					break;
				case ePT_Peng_0_12:
					// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
					host_device_detail::feg_Peng_neutral_0_12(g, c_feg, y);
					break;
				case ePT_Kirkland_0_12:
					// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
					host_device_detail::feg_Kirkland_neutral_0_12(g, c_feg, y);
					break;
				case ePT_Weickenmeier_0_12:
					// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
					host_device_detail::feg_Weickenmeier_neutral_0_12(g, c_feg, y);
					break;
				case ePT_Lobato_0_12:
					// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
					host_device_detail::feg_Lobato_neutral_0_12(g, c_feg, y);
					break;
			}
		}
		else
		{
			host_device_detail::feg_Peng_ion_0_4(g, c_feg, y);
		}
	}

	// Electron scattering factor(c_feg, dfeg) where dfg is the first derivative along g
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	void feg_dfeg(const ePotential_Type &potential_type, const int &charge, const T &g, const rPP_Coef<T> &c_feg, T &y, T &dy)
	{
		if(charge == 0)
		{
			switch(potential_type)
			{
				case ePT_Doyle_0_4:
					// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
					host_device_detail::feg_dfeg_Doyle_neutral_0_4(g, c_feg, y, dy);
					break;
				case ePT_Peng_0_4:
					// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
					host_device_detail::feg_dfeg_Peng_neutral_0_4(g, c_feg, y, dy);
					break;
				case ePT_Peng_0_12:
					// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
					host_device_detail::feg_dfeg_Peng_neutral_0_12(g, c_feg, y, dy);
					break;
				case ePT_Kirkland_0_12:
					// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
					host_device_detail::feg_dfeg_Kirkland_neutral_0_12(g, c_feg, y, dy);
					break;
				case ePT_Weickenmeier_0_12:
					// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
					host_device_detail::feg_dfeg_Weickenmeier_neutral_0_12(g, c_feg, y, dy);
					break;
				case ePT_Lobato_0_12:
					// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
					host_device_detail::feg_dfeg_Lobato_neutral_0_12(g, c_feg, y, dy);
					break;
			}
		}
		else
		{
			host_device_detail::feg_dfeg_Peng_ion_0_4(g, c_feg, y, dy);
		}
	}

	// Electron scattering factor(fg)
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	void fxg(const ePotential_Type &potential_type, const int &charge, const int &Z, const T &g, const rPP_Coef<T> &c_fxg, T &y)
	{
		if(charge == 0)
		{
			switch(potential_type)
			{
				case ePT_Doyle_0_4:
					// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
					host_device_detail::fxg_Doyle_neutral_0_4(Z, g, c_fxg, y);
					break;
				case ePT_Peng_0_4:
					// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
					host_device_detail::fxg_Peng_neutral_0_4(Z, g, c_fxg, y);
					break;
				case ePT_Peng_0_12:
					// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
					host_device_detail::fxg_Peng_neutral_0_12(Z, g, c_fxg, y);
					break;
				case ePT_Kirkland_0_12:
					// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
					host_device_detail::fxg_Kirkland_neutral_0_12(Z, g, c_fxg, y);
					break;
				case ePT_Weickenmeier_0_12:
					// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
					host_device_detail::fxg_Weickenmeier_neutral_0_12(g, c_fxg, y);
					break;
				case ePT_Lobato_0_12:
					// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
					host_device_detail::fxg_Lobato_neutral_0_12(g, c_fxg, y);
					break;
			}
		}
		else
		{
			host_device_detail::fxg_Peng_ion_0_4(Z-charge, g, c_fxg, y);
		}
	}

	// Electron scattering factor(fg, dfg) where dfg is the first derivative along g
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	void fxg_dfxg(const ePotential_Type &potential_type, const int &charge, const int &Z, const T &g, const rPP_Coef<T> &c_fxg, T &y, T &dy)
	{
		if(charge == 0)
		{
			switch(potential_type)
			{
				case ePT_Doyle_0_4:
					// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
					host_device_detail::fxg_dfxg_Doyle_neutral_0_4(Z, g, c_fxg, y, dy);
					break;
				case ePT_Peng_0_4:
					// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
					host_device_detail::fxg_dfxg_Peng_neutral_0_4(Z, g, c_fxg, y, dy);
					break;
				case ePT_Peng_0_12:
					// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
					host_device_detail::fxg_dfxg_Peng_neutral_0_12(Z, g, c_fxg, y, dy);
					break;
				case ePT_Kirkland_0_12:
					// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
					host_device_detail::fxg_dfxg_Kirkland_neutral_0_12(Z, g, c_fxg, y, dy);
					break;
				case ePT_Weickenmeier_0_12:
					// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
					host_device_detail::fxg_dfxg_Weickenmeier_neutral_0_12(g, c_fxg, y, dy);
					break;
				case ePT_Lobato_0_12:
					// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
					host_device_detail::fxg_dfxg_Lobato_neutral_0_12(g, c_fxg, y, dy);
					break;
			}
		}
		else
		{
			host_device_detail::fxg_dfxg_Peng_ion_0_4(Z-charge, g, c_fxg, y, dy);
		}
	}

	// Electron density (Pr)
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	void Pr(const ePotential_Type &potential_type, const int &charge, const T &r, const rPP_Coef<T> &c_Pr, T &y)
	{
		if(charge == 0)
		{
			switch(potential_type)
			{
				case ePT_Doyle_0_4:
					// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
					host_device_detail::Pr_Doyle_neutral_0_4(r, c_Pr, y);
					break;
				case ePT_Peng_0_4:
					// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
					host_device_detail::Pr_Peng_neutral_0_4(r, c_Pr, y);
					break;
				case ePT_Peng_0_12:
					// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
					host_device_detail::Pr_Peng_neutral_0_12(r, c_Pr, y);
					break;
				case ePT_Kirkland_0_12:
					// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
					host_device_detail::Pr_Kirkland_neutral_0_12(r, c_Pr, y);
					break;
				case ePT_Weickenmeier_0_12:
					// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
					host_device_detail::Pr_Weickenmeier_neutral_0_12(r, c_Pr, y);
					break;
				case ePT_Lobato_0_12:
					// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
					host_device_detail::Pr_Lobato_neutral_0_12(r, c_Pr, y);
					break;
			}
		}
		else
		{
			host_device_detail::Pr_Peng_ion_0_4(r, c_Pr, y);
		}
	}

	// Electron density (c_Pr, dPr) where dPr is the first derivative along r
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	void Pr_dPr(const ePotential_Type &potential_type, const int &charge, const T &r, const rPP_Coef<T> &c_Pr, T &y, T &dy)
	{
		if(charge == 0)
		{
			switch(potential_type)
			{
				case ePT_Doyle_0_4:
					// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
					host_device_detail::Pr_dPr_Doyle_neutral_0_4(r, c_Pr, y, dy);
					break;
				case ePT_Peng_0_4:
					// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
					host_device_detail::Pr_dPr_Peng_neutral_0_4(r, c_Pr, y, dy);
					break;
				case ePT_Peng_0_12:
					// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
					host_device_detail::Pr_dPr_Peng_neutral_0_12(r, c_Pr, y, dy);
					break;
				case ePT_Kirkland_0_12:
					// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
					host_device_detail::Pr_dPr_Kirkland_neutral_0_12(r, c_Pr, y, dy);
					break;
				case ePT_Weickenmeier_0_12:
					// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
					host_device_detail::Pr_dPr_Weickenmeier_neutral_0_12(r, c_Pr, y, dy);
					break;
				case ePT_Lobato_0_12:
					// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
					host_device_detail::Pr_dPr_Lobato_neutral_0_12(r, c_Pr, y, dy);
					break;
			}
		}
		else
		{
			host_device_detail::Pr_dPr_Peng_ion_0_4(r, c_Pr, y, dy);
		}
	}

	// Projected_Potential calculation(Vr)
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	void Vr(const ePotential_Type &potential_type, const int &charge, const T &r, const rPP_Coef<T> &c_Vr, T &y)
	{
		if(charge == 0)
		{
			switch(potential_type)
			{
				case ePT_Doyle_0_4:
					// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
					host_device_detail::Vr_Doyle_neutral_0_4(r, c_Vr, y);
					break;
				case ePT_Peng_0_4:
					// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
					host_device_detail::Vr_Peng_neutral_0_4(r, c_Vr, y);
					break;
				case ePT_Peng_0_12:
					// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
					host_device_detail::Vr_Peng_neutral_0_12(r, c_Vr, y);
					break;
				case ePT_Kirkland_0_12:
					// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
					host_device_detail::Vr_Kirkland_neutral_0_12(r, c_Vr, y);
					break;
				case ePT_Weickenmeier_0_12:
					// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
					host_device_detail::Vr_Weickenmeier_neutral_0_12(r, c_Vr, y);
					break;
				case ePT_Lobato_0_12:
					// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
					host_device_detail::Vr_Lobato_neutral_0_12(r, c_Vr, y);
					break;
			}
		}
		else
		{
			host_device_detail::Vr_Peng_ion_0_4(r, c_Vr, y);
		}
	}

	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	void Vr_dVr(const ePotential_Type &potential_type, const int &charge, const T &r, const rPP_Coef<T> &c_Vr, T &y, T &dy)
	{
		if(charge == 0)
		{
			switch(potential_type)
			{
				case ePT_Doyle_0_4:
					// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
					host_device_detail::Vr_dVr_Doyle_neutral_0_4(r, c_Vr, y, dy);
					break;
				case ePT_Peng_0_4:
					// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
					host_device_detail::Vr_dVr_Peng_neutral_0_4(r, c_Vr, y, dy);
					break;
				case ePT_Peng_0_12:
					// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
					host_device_detail::Vr_dVr_Peng_neutral_0_12(r, c_Vr, y, dy);
					break;
				case ePT_Kirkland_0_12:
					// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
					host_device_detail::Vr_dVr_Kirkland_neutral_0_12(r, c_Vr, y, dy);
					break;
				case ePT_Weickenmeier_0_12:
					// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
					host_device_detail::Vr_dVr_Weickenmeier_neutral_0_12(r, c_Vr, y, dy);
					break;
				case ePT_Lobato_0_12:
					// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
					host_device_detail::Vr_dVr_Lobato_neutral_0_12(r, c_Vr, y, dy);
					break;
			}
		}
		else
		{
			host_device_detail::Vr_dVr_Peng_ion_0_4(r, c_Vr, y, dy);
		}
	}

	// Projected potential (VR)
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	void VR(const ePotential_Type &potential_type, const int &charge, const T &R, const rPP_Coef<T> &c_VR, const rQ1<T> &Qz_0_I, T &y)
	{
		if(charge == 0)
		{
			switch(potential_type)
			{
				case ePT_Doyle_0_4:
					// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
					host_device_detail::VR_Doyle_neutral_0_4(R, c_VR, y);
					break;
				case ePT_Peng_0_4:
					// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
					host_device_detail::VR_Peng_neutral_0_4(R, c_VR, y);
					break;
				case ePT_Peng_0_12:
					// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
					host_device_detail::VR_Peng_neutral_0_12(R, c_VR, y);
					break;
				case ePT_Kirkland_0_12:
					// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
					host_device_detail::VR_Kirkland_neutral_0_12(R, c_VR, y);
					break;
				case ePT_Weickenmeier_0_12:
					// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
					host_device_detail::VR_Weickenmeier_neutral_0_12(R, c_VR, Qz_0_I, y);
					break;
				case ePT_Lobato_0_12:
					// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
					host_device_detail::VR_Lobato_neutral_0_12(R, c_VR, y);
					break;
			}
		}
		else
		{
			host_device_detail::VR_Peng_ion_0_4(R, c_VR, y);
		}
	}

	// Projected potential (c_VR, dVR) where dVr is the first derivative along R
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	void VR_dVR(const ePotential_Type &potential_type, const int &charge, const T &R, const rPP_Coef<T> &c_VR, const rQ1<T> &Qz_0_I, T &y, T &dy)
	{
		if(charge == 0)
		{
			switch(potential_type)
			{
				case ePT_Doyle_0_4:
					// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
					host_device_detail::VR_dVR_Doyle_neutral_0_4(R, c_VR, y, dy);
					break;
				case ePT_Peng_0_4:
					// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
					host_device_detail::VR_dVR_Peng_neutral_0_4(R, c_VR, y, dy);
					break;
				case ePT_Peng_0_12:
					// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
					host_device_detail::VR_dVR_Peng_neutral_0_12(R, c_VR, y, dy);
					break;
				case ePT_Kirkland_0_12:
					// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
					host_device_detail::VR_dVR_Kirkland_neutral_0_12(R, c_VR, y, dy);
					break;
				case ePT_Weickenmeier_0_12:
					// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
					host_device_detail::VR_dVR_Weickenmeier_neutral_0_12(R, c_VR, Qz_0_I, y, dy);
					break;
				case ePT_Lobato_0_12:
					// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
					host_device_detail::VR_dVR_Lobato_neutral_0_12(R, c_VR, y, dy);
					break;
			}
		}
		else
		{
			host_device_detail::VR_dVR_Peng_ion_0_4(R, c_VR, y, dy);
		}
	}

	// Projected potential (Vz)[z0, ze]
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	void Vz(const ePotential_Type &potential_type, const int &charge, const T &z0, const T &ze, const T &R, const rPP_Coef<T> &c_Vr, const rQ1<T> &Qz_a_b, T &y)
	{
		if(charge == 0)
		{
			switch(potential_type)
			{
				case ePT_Doyle_0_4:
					// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
					host_device_detail::Vz_Doyle_neutral_0_4(z0, ze, R, c_Vr, Qz_a_b, y);
					break;
				case ePT_Peng_0_4:
					// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
					host_device_detail::Vz_Peng_neutral_0_4(z0, ze, R, c_Vr, Qz_a_b, y);
					break;
				case ePT_Peng_0_12:
					// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
					host_device_detail::Vz_Peng_neutral_0_12(z0, ze, R, c_Vr, Qz_a_b, y);
					break;
				case ePT_Kirkland_0_12:
					// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
					host_device_detail::Vz_Kirkland_neutral_0_12(z0, ze, R, c_Vr, Qz_a_b, y);
					break;
				case ePT_Weickenmeier_0_12:
					// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
					host_device_detail::Vz_Weickenmeier_neutral_0_12(z0, ze, R, c_Vr, Qz_a_b, y);
					break;
				case ePT_Lobato_0_12:
					// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
					host_device_detail::Vz_Lobato_neutral_0_12(z0, ze, R, c_Vr, Qz_a_b, y);
					break;
			}
		}
		else
		{
			host_device_detail::Vz_Peng_ion_0_4(z0, ze, R, c_Vr, Qz_a_b, y);
		}
	}

	// Projected potential (Vz, dVz)[z0, ze] where dVr is the first derivative along R
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	void Vz_dVz(const ePotential_Type &potential_type, const int &charge, const T &z0, const T &ze, const T &R, const rPP_Coef<T> &c_Vr, const rQ1<T> &Qz_a_b, T &y, T &dy)
	{
		if(charge == 0)
		{
			switch(potential_type)
			{
				case ePT_Doyle_0_4:
					// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
					host_device_detail::Vz_dVz_Doyle_neutral_0_4(z0, ze, R, c_Vr, Qz_a_b, y, dy);
					break;
				case ePT_Peng_0_4:
					// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
					host_device_detail::Vz_dVz_Peng_neutral_0_4(z0, ze, R, c_Vr, Qz_a_b, y, dy);
					break;
				case ePT_Peng_0_12:
					// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
					host_device_detail::Vz_dVz_Peng_neutral_0_12(z0, ze, R, c_Vr, Qz_a_b, y, dy);
					break;
				case ePT_Kirkland_0_12:
					// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
					host_device_detail::Vz_dVz_Kirkland_neutral_0_12(z0, ze, R, c_Vr, Qz_a_b, y, dy);
					break;
				case ePT_Weickenmeier_0_12:
					// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
					host_device_detail::Vz_dVz_Weickenmeier_neutral_0_12(z0, ze, R, c_Vr, Qz_a_b, y, dy);
					break;
				case ePT_Lobato_0_12:
					// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
					host_device_detail::Vz_dVz_Lobato_neutral_0_12(z0, ze, R, c_Vr, Qz_a_b, y, dy);
					break;
			}
		}
		else
		{
			host_device_detail::Vz_dVz_Peng_ion_0_4(z0, ze, R, c_Vr, Qz_a_b, y, dy);
		}
	}

	template <ePotential_Type potential_type, int charge, class T>
	DEVICE_CALLABLE FORCE_INLINE
	void Vr_dVrir(const T &r, T *cl, T *cnl, T f, T &Vr, T &dVrir)
	{
		if(charge == 0)
		{
			switch(potential_type)
			{
				case ePT_Doyle_0_4:
					// 1: Doyle and Turner parameterization - 4 Gaussians - [0, 4]
					host_device_detail::Vr_dVrir_Doyle_neutral_0_4<T>(r, cl, cnl, Vr, dVrir);
					break;
				case ePT_Peng_0_4:
					// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
					host_device_detail::Vr_dVrir_Peng_neutral_0_4_12<T>(r, cl, cnl, Vr, dVrir);
					break;
				case ePT_Peng_0_12:
					// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
					host_device_detail::Vr_dVrir_Peng_neutral_0_4_12<T>(r, cl, cnl, Vr, dVrir);
					break;
				case ePT_Kirkland_0_12:
					// 4: Kirkland parameterization - 3 Lorentzian + 3 Gaussians - [0, 12]				
					host_device_detail::Vr_dVrir_Kirkland_neutral_0_12<T>(r, cl, cnl, Vr, dVrir);
					break;
				case ePT_Weickenmeier_0_12:
					// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
					host_device_detail::Vr_dVrir_Weickenmeier_neutral_0_12<T>(r, cl, cnl, Vr, dVrir);
					break;
				case ePT_Lobato_0_12:
					// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
					host_device_detail::Vr_dVrir_Lobato_neutral_0_12<T>(r, cl, cnl, Vr, dVrir);
					break;
			}
		}
		else
		{
			host_device_detail::Vr_dVrir_Peng_ion_0_4<T>(r, cl, cnl, Vr, dVrir);
		}
		Vr *= f;
		dVrir *= f;
	}

	/***************************************************************************/
	/***************************************************************************/

	namespace functor
	{
		template <class T>
		struct transmission_function
		{
			const eElec_Spec_Int_Model elec_spec_int_model;

			const T w;
			transmission_function(T w_i, eElec_Spec_Int_Model ElecSpecIntModel_i): w(w_i), elec_spec_int_model(ElecSpecIntModel_i){}

			template <class U>
			DEVICE_CALLABLE
			complex<T> operator()(const U &x) const 
			{ 
				T theta = w*x;
				return (elec_spec_int_model == eESIM_Weak_Phase_Object)?complex<T>(1.0, theta):euler(theta);
			}
		};

		template <class T>
		struct closest_element
		{
			const T m_val;

			closest_element(T val): m_val(val) {}

			bool operator()(const T &a, const T &b)
			{	
				const T da = 1e-4;
				return fabs(a-m_val-da)<fabs(b-m_val); 
			};
		};

		template <class T>
		struct assign_real
		{
			template <class U>
			DEVICE_CALLABLE
			T operator()(const U &x) const { return x.real(); }
		};

		template <class T>
		struct assign_max_real
		{
			const T vm;
			assign_max_real(T vm_i = T()): vm(vm_i){}

			template <class U>
			DEVICE_CALLABLE
			T operator()(const U &x) const { return ::fmax(vm, x.real()); }
		};

		template <class T>
		struct scale
		{
			const T w;
			scale(T w_i = T()): w(w_i){}

			template <class U>
			DEVICE_CALLABLE
			T operator()(const U &x) const { return w*x; }
		};

		template <class T>
		DEVICE_CALLABLE
		T norm(const thrust::device_reference< thrust::complex<T> >&x) 
		{
		thrust::complex<T> xx = (thrust::complex<T>)x;
			return norm(xx);    
		}

		template <class T>
		struct square
		{
			template <class U>
			DEVICE_CALLABLE
			T operator()(const U &x) const { return norm(x); }
		};

		template <class T>
		struct square_scale
		{
			const T w;
			square_scale(T w_i = T()): w(w_i){}

			template <class U>
			DEVICE_CALLABLE
			T operator()(const U &x) const { return w*norm(x); }
		};

		template <class T>
		struct add
		{
			add(){}

			template <class U, class V>
			DEVICE_CALLABLE
			T operator()(const U &lhs, const V &rhs) const { return lhs + rhs; }
		};

		template <class T>
		struct multiply
		{
			multiply(){}

			template <class U, class V>
			DEVICE_CALLABLE
			T operator()(const U &lhs, const V &rhs) const { return lhs*rhs; }
		};

		template <class T>
		struct add_scale
		{
			const T w;
			add_scale(T w_i = T()): w(w_i){}

			template <class U, class V>
			DEVICE_CALLABLE
			T operator()(const U &lhs, const V &rhs) const{ return w*lhs + rhs; }
		};

		template <class T>
		struct add_scale_i
		{
			const T w1;
			const T w2;
			add_scale_i(T w1_i = T(), T w2_i = T()): w1(w1_i), w2(w2_i){}

			template <class U, class V>
			DEVICE_CALLABLE
			T operator()(const U &lhs, const V &rhs) const{ return w1*lhs + w2*rhs; }
		};

		template <class T>
		struct add_square
		{
			template <class U, class V>
			DEVICE_CALLABLE
			T operator()(const U &lhs, const V &rhs) const { return norm(lhs) + rhs; }
		};

		template <class T>
		struct add_square_i
		{
			template <class U, class V>
			DEVICE_CALLABLE
			T operator()(const U &lhs, const V &rhs) const { return norm(lhs) + norm(rhs); }
		};

		template <class T>
		struct add_scale_square
		{
			const T w;
			add_scale_square(T w_i = T()): w(w_i){}

			template <class U, class V>
			DEVICE_CALLABLE
			T operator()(const U &lhs, const V &rhs) const { return w*norm(lhs) + rhs; }
		};

		template <class T>
		struct add_scale_square_i
		{
			const T w1;
			const T w2;
			add_scale_square_i(T w1_i = T(), T w2_i = T()): w1(w1_i), w2(w2_i){}

			template <class U, class V>
			DEVICE_CALLABLE
			T operator()(const U &lhs, const V &rhs) const { return w1*::norm(lhs) + w2*::norm(rhs); }
		};

		template <class T, class Tr>
		struct square_dif
		{
			const T x_mean;
			square_dif(T x_mean_i = T()): x_mean(x_mean_i){}

			template <class U>
			DEVICE_CALLABLE
			Tr operator()(const U &x) const { return ::norm(x-x_mean); }
		};


		template <class T>
		struct binarize
		{
			const T thr;
			binarize(T thr_i = T()): thr(thr_i){}

			template <class U>
			DEVICE_CALLABLE
			T operator()(const U &x) const { return (x<thr)?0:1; }
		};

		template <class T>
		struct thresholding
		{
			const T thr;
			thresholding(T thr_i = T()): thr(thr_i){}

			template <class U>
			DEVICE_CALLABLE
			T operator()(const U &x) const { return (x<thr)?thr:x; }
		};

		template <class T>
		struct anscombe_forward 
		{
			const T xs;
			anscombe_forward(): xs(3.0/8.0){}

			template <class U>
			DEVICE_CALLABLE
			T operator()(const U &x) const { return (x<0)?0:2*sqrt(x+xs); }
		};

		template <class T>
		struct anscombe_inverse 
		{
			const T a;
			const T b;
			const T c;
			const T d;
			const T e;
			anscombe_inverse(): a(1.0/4.0), b(sqrt(3.0/2.0)/4), 
				c(-11.0/8.0), d(5*sqrt(3.0/2.0)/8), e(-1.0/8.0){}\

			template <class U>
			DEVICE_CALLABLE
			T operator()(const U &x) const 
			{ 
				if(isZero(x))
				{
					return ::fmax(0, a*x+e);
				}
				else
				{
					T ix = 1.0/x;
					return ::fmax(0, a*x*x+e+ix*(b+ix*(c+ix*d)));				
				}
			}
		};
	} // namespace functor

	template <class TVector>
	Value_type<TVector> min_element(TVector &x)
	{
		return *thrust::min_element(x.begin(), x.end());
	}

	template <class TVector, class ...TArg>
	Value_type<TVector> min_element(TVector &x, TArg &...arg)
	{
		return ::fmin(min_element(x), min_element(arg...));
	}

	template <class TVector>
	Value_type<TVector> max_element(TVector &x)
	{
		return *thrust::max_element(x.begin(), x.end());
	}

	template <class TVector, class ...TArg>
	Value_type<TVector> max_element(TVector &x, TArg &...arg)
	{
		return ::fmax(max_element(x), max_element(arg...));
	}
	
	template <class TVector>
	void minmax_element(TVector &x, Value_type<TVector> &x_min, Value_type<TVector> &x_max)
	{
		auto x_min_max = thrust::minmax_element(x.begin(), x.end());
		x_min = *(x_min_max.first);
		x_max = *(x_min_max.second);
	}

	template <class TVector>
	Value_type<TVector> mean(TVector &M_i)
	{
		return thrust::reduce(M_i.begin(), M_i.end())/M_i.size();
	}

	template <class TVector>
	void mean_var(TVector &M_i, Value_type<TVector> &x_mean, Value_type_r<TVector> &x_var)
	{
		using T = Value_type<TVector>;
		using T_r = Value_type_r<TVector>;

		x_mean = mean(M_i);

		x_var = thrust::transform_reduce(M_i.begin(), M_i.end(), 
		functor::square_dif<T, T_r>(x_mean), T_r(0), functor::add<T_r>());

		x_var = x_var/M_i.size();
	}

	template <class TVector>
	void mean_std(TVector &M_i, Value_type<TVector> &x_mean, Value_type_r<TVector> &x_std)
	{
		mean_var(M_i, x_mean, x_std);
		x_std = sqrt(x_std);
	}

	template <class TVector>
	Value_type_r<TVector> variance(TVector &M_i)
	{
		using T = Value_type<TVector>;
		using T_r = Value_type_r<TVector>;

		T x_mean;
		T_r x_var;
		mean_var(M_i, x_mean, x_var);
		return x_var;
	}

	template<class TVector>
	void scale(Value_type<TVector> f, TVector &x)
	{
		using value_type = Value_type<TVector>;
		thrust::transform(x.begin(), x.end(), x.begin(), functor::scale<value_type>(f));
	}

	template <class TVector_c, class TVector_r>
	enable_if_complex_vector_and_real_vector<TVector_c, TVector_r, void>
	assign_real(TVector_c &M_i, TVector_r &M_o)
	{
		using value_type = Value_type<TVector_r>;

		thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), functor::assign_real<value_type>());
	}

	template <class TVector_c, class TVector_r>
	enable_if_complex_vector_and_real_vector<TVector_c, TVector_r, void>
	assign_real(TVector_c &M_i, TVector_r &M_o, Value_type<TVector_r> M_v)
	{
		using value_type = Value_type<TVector_r>;

		if(M_v>0)
		{
			thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), functor::assign_max_real<value_type>(0));
		}
		else
		{
			thrust::transform(M_i.begin(), M_i.end(), M_o.begin(), functor::assign_real<value_type>());
		}
	}

} // namespace mt

#endif
