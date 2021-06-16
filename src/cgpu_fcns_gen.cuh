/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef CGPU_FCNS_GEN_H
	#define CGPU_FCNS_GEN_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include <cmath>
	#include <array>
	#include <vector>
	#include <string>
	#include <cctype>
	#include <locale>

	#include "macros.cuh"
	#include "const_enum.cuh"
	#include "math.cuh"
	#include "type_traits_gen.cuh"

	/********************************* simple functions ************************************/
	/********** static member function are not supported for the cuda compiler *************/
	/***************************************************************************************/
	//template <class T>
	//CGPU_EXEC
	//T fill_nzero(const T& r) // for R_2d and R_3d compatibility 
	//{
	//	return r;
	//}


	/* is equal */
	namespace mt
	{
		template <class T>
		CGPU_EXEC_INL
		dt_bool fcn_is_equal(const T& a, const T& b)
		{
			return a == b;
		}

		template <>
		CGPU_EXEC_INL
		dt_bool fcn_is_equal<dt_int32>(const dt_int32& a, const dt_int32& b)
		{
			return a == b;
		}

		template <>
		CGPU_EXEC_INL
		dt_bool fcn_is_equal<dt_float32>(const dt_float32 &a, const dt_float32 &b)
		{
			const dt_float32 eps_abs = 1e-5f;
			const dt_float32 eps_rel = 1e-4f;

			// check if the numbers are really close -- needed when comparing numbers near zero.
			dt_float32 diff = fabs(a - b);
			if (diff <= eps_abs)
				return true;
 
			// otherwise fall back to Knuth's algorithm
			return diff <= ((fabs(a)<fabs(b)?fabs(b):fabs(a))*eps_rel);
		}

		template <>
		CGPU_EXEC_INL
		dt_bool fcn_is_equal<dt_float64>(const dt_float64 &a, const dt_float64 &b)
		{
			const dt_float64 eps_abs = 1e-13;
			const dt_float64 eps_rel = 1e-8;

			// check if the numbers are really close -- needed when comparing numbers near zero.
			dt_float64 diff = fabs(a - b);
			if (diff <= eps_abs)
				return true;
 
			// otherwise fall back to Knuth's algorithm
			return diff <= ((fabs(a)<fabs(b)?fabs(b):fabs(a))*eps_rel);
		}
	}

	/* is zero */
	namespace mt
	{
		template <class T>
		CGPU_EXEC_INL
		dt_bool fcn_is_zero(const T& x)
		{
			return fcn_is_equal(x, T(0));
		}

		template <class T, class... TArg>
		CGPU_EXEC_INL
		dt_bool fcn_is_zero(const T& x, const TArg&... arg)
		{
			return fcn_is_zero(x) && fcn_is_zero(arg...);
		}

		template <class T>
		CGPU_EXEC_INL
		dt_bool fcn_is_nzero(const T& x)
		{
			return !fcn_is_equal(x, T(0));
		}

		template <class T, class... TArg>
		CGPU_EXEC_INL
		dt_bool fcn_is_nzero(const T& x, const TArg&... arg)
		{
			return fcn_is_nzero(x) && fcn_is_nzero(arg...);
		}
	}

	/* is one */
	namespace mt
	{
		template <class T>
		CGPU_EXEC_INL
		dt_bool fcn_is_one(const T& x)
		{
			return fcn_is_equal(x, T(1));
		}

		template <class T, class... TArg>
		CGPU_EXEC_INL
		dt_bool fcn_is_one(const T& x, const TArg&... arg)
		{
			return fcn_is_one(x) && fcn_is_one(arg...);
		}
	}

	/* if else zero */
	namespace mt
	{
		template <class T>
		CGPU_EXEC
		T fcn_set_zero_eps(const T& x)
		{
			return (fcn_is_zero(x))?T(0):x;
		}		
				
		template <class T>
		CGPU_EXEC
		T fcn_ifzero(const T& x, const T& x_t, const T& x_f)
		{
			return (fcn_is_zero(x))?x_t:x_f;
		}					
		
		template <class T>
		CGPU_EXEC
		T fcn_ifnzero(const T& x, const T& x_t, const T& x_f)
		{
			return (fcn_is_nzero(x))?x_t:x_f;
		}		
	}

	/* div/min/max */
	namespace mt
	{
		template <class T, class U>
		CGPU_EXEC
		T fcn_div(const T& x, const U& y)
		{
			return (fcn_is_zero(y))?T(0): x/static_cast<T>(y);
		}

		template <class T>
		CGPU_EXEC
		T fcn_min(const T& x, const T& y)
		{
			return (x<y)?x:y;
		}	
	
		template <class T>
		CGPU_EXEC
		T fcn_max(const T& x, const T& y)
		{
			return (x>y)?x:y;
		}
	}

	/*  check bounds */
	namespace mt
	{
		template <class T>
		CGPU_EXEC
		dt_bool fcn_chk_bound(const T& x, const T& x_min, const T& x_max)
		{
			return (x_min <= x) && (x < x_max);
		}

		// template <class T, class... TArg>
		// CGPU_EXEC
		// dt_bool fcn_chk_bound(const T& x, const T& x_min, const T& x_max, const TArg&... arg)
		// {
		// 	return fcn_chk_bound(x, x_min, x_max) && fcn_chk_bound(x, arg...);
		// }

		template <class T>
		CGPU_EXEC
		dt_bool fcn_chk_bound_eps(const T& x, const T& x_min, const T& x_max)
		{
			const T eps = epsilon_abs<T>();
			return (x_min - eps <= x) && (x < x_max + eps);
		}

		template <class T, class... TArg>
		CGPU_EXEC
		dt_bool fcn_chk_bound_eps(const T& x, const T& x_min, const T& x_max, const TArg&... arg)
		{
			return fcn_chk_bound_eps(x, x_min, x_max) && fcn_chk_bound_eps(arg...);
		}
	}

	/*  set bounds */
	namespace mt
	{
		template <class T>
		CGPU_EXEC
		enable_if_float<T, T>
		fcn_set_bound(const T& x, const T& x_min, const T& x_max)
		{
			return ::fmax(x_min, ::fmin(x_max, x));
		}

		template <class T>
		CGPU_EXEC
		enable_if_int<T, T>
		fcn_set_bound(const T& x, const T& x_min, const T& x_max)
		{
			return (x<=x_min)?x_min:(x>=x_max)?x_max:x;
		}
	}

	/* locate index: floor/ceil/pbc */
	namespace mt
	{
		// ! apply periodic boundary conditions to index
		template <class T, class ST>
		CGPU_EXEC
		ST fcn_ind_pbc(const ST& ix, const ST& nx) 
		{
			return ix - static_cast<ST>(::floor(T(ix)/T(nx)))*nx;
		}

		// ! cast floor
		template <class ST, class T>
		CGPU_EXEC
		ST fcn_cfloor(const T& x) 
		{
			return static_cast<ST>(::floor(x));
		}

		// ! cast ceil
		template <class ST, class T>
		CGPU_EXEC
		ST fcn_cceil(const T& x) 
		{
			return static_cast<ST>(::ceil(x));
		}

		// ! bound cast floor
		template <class T, class U>
		CGPU_EXEC
		U fcn_bcfloor(const T& x, const U& ix_min, const U& ix_max) 
		{
			return fcn_set_bound(fcn_cfloor<U>(x), ix_min, ix_max);
		}

		// ! bound cast ceil
		template <class T, class U>
		CGPU_EXEC
		U fcn_bcceil(const T& x, const U& ix_min, const U& ix_max) 
		{
			return fcn_set_bound(fcn_cceil<U>(x), ix_min, ix_max);
		}

		// get index and distance
		template <class T, class ST>
		CGPU_EXEC
		void fcn_get_idx_0_idx_n(const T& x, const T& x_max, const T& dx, const dt_bool& pbc, const ST& nx, ST& ix_0, ST& ix_n)
		{
			ix_0 = fcn_cfloor<ST>((x - x_max)/dx);
			auto ix_e = fcn_cceil<ST>((x + x_max)/dx);

			if (!pbc)
			{
				ix_0 = fcn_set_bound(ix_0, ST(0), nx);
				ix_e = fcn_set_bound(ix_e, ST(0), nx);
			}

			ix_n = (ix_0 == ix_e)?ST(0):(ix_e - ix_0 + 1);
		}
	}

	/* index search */
	namespace mt
	{
		template <class T, class TFcn>
		CGPU_EXEC
		dt_int32 fcn_r_2_ir_b_by_fcn(TFcn fcn, const T& x, dt_int32 ix_min, dt_int32 ix_max)
		{
			do 
			{
				const dt_int32 ix_m = (ix_min + ix_max) >> 1;		// divide by 2

				if (x <= fcn(ix_m))
					ix_max = ix_m;
				else
					ix_min = ix_m;
			} while ((ix_max - ix_min) > 1);

			if (x >= fcn(ix_max))
				return ix_max;
			else
				return ix_min;
		}

		template <class T>
		CGPU_EXEC
		dt_int32 fcn_r_2_ir_by_vctr(T* xv, const T& x, dt_int32 ix_min, dt_int32 ix_max)
		{
			do {
				const dt_int32 ix_m = (ix_min + ix_max) >> 1;		// divide by 2

				if (x <= xv[ix_m])
					ix_max = ix_m;
				else
					ix_min = ix_m;
			} while ((ix_max - ix_min) > 1);

			if (x >= xv[ix_max])
				return ix_max;
			else
				return ix_min;
		}

		template <class T>
		CGPU_EXEC
		dt_int32 fcn_r_2_ir_b_by_vctr(T* xv, const T& x, dt_int32 ix_min, dt_int32 ix_max)
		{
			if (x > xv[ix_max])
				return -2;
			else if (x < x[ix_min])
				return -1;

			return fcn_r_2_ir_by_vctr(xv, x, ix_min, ix_max);
		}
	}

	/* find root */
	namespace mt
	{
		template <class T, class TFcn>
		CGPU_EXEC
		T fcn_fd_root(T x_0, T x_e, const T& Tol, const dt_int32& it_max, TFcn fcn)
		{
			dt_int32 it = 0;
			T x, fx, fx_e = fcn(x_e);

			do
			{
				x = 0.5*(x_0 + x_e);
				fx = fcn(x);

				if (fx*fx_e<0)			
				{
					x_0 = x;
				}
				else
				{
					x_e = x;
					fx_e = fx;
				}
				it++;
			} while ((::fabs(fx)>Tol) && (it < it_max));
 
			return x;
		}
	}

	/* point 2 line 2d*/
	namespace mt
	{
		// ax + by + c = 0
		// https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line

		// distance from point to line
		template <class T>
		CGPU_EXEC 
		T fcn_pt_2_ln_dist(const T& a, const T& b, const T& c, const T& x0, const T& y0)
		{
			return ::fabs(a*x0 + b*y0 + c)/::sqrt(a*a + b*b);
		}

		// calculate intersection from point to line
		template <class T>
		CGPU_EXEC 
		void fcn_pt_2_ln_intxn_2d(const T& a, const T& b, const T& c, const T& x0, const T& y0, const T& x, const T& y)
		{
			auto m = ::sqrt(a*a + b*b);

			x = (b*(b*x0-a*y0)-a*c)/m;
			y = -(a*(b*x0-a*y0)-b*c)/m;
		}
	}

	/* length of curve */
	namespace mt
	{
		template <class T>
		CGPU_EXEC 
		T fcn_get_len_crv(dt_int32 ix_0, dt_int32 ix_e, Ctpr<T> x, Ctpr<T> y)
		{
			T ln = 0;
			for(auto ik = ix_0; ik < ix_e-1; ik++)
			{
				const T dx = x[ik + 1] - x[ik];
				const T dy = y[ik + 1] - y[ik];
				ln = ln + ::sqrt(dx*dx + dy*dy);
			}

			return ln;
		}

		template <class T>
		CGPU_EXEC 
		T fcn_get_len_crv(dt_int32 ix_0, dt_int32 ix_e, Ctpr<T> x, Ctpr<T> y, const T& ln_max, dt_int32& iln)
		{
			T ln = 0;

			if (ix_0 < ix_e)
			{
				iln = ix_e;
				for(auto ik = ix_0; ik < ix_e-1; ik++)
				{
					const T dx = x[ik + 1] - x[ik];
					const T dy = y[ik + 1] - y[ik];
					ln = ln + ::sqrt(dx*dx + dy*dy);
					if ((ln_max > 0) && (ln >= ln_max))
					{
		 				iln = ik;
		 				break;
					}
				}
			}
			else
			{
				iln = ix_e;
				for(auto ik = ix_0; ik > ix_e; ik--)
				{
					const T dx = x[ik] - x[ik-1];
					const T dy = y[ik] - y[ik-1];
					ln = ln + ::sqrt(dx*dx + dy*dy);
					if ((ln_max > 0) && (ln >= ln_max))
					{
		 				iln = ik;
		 				break;
					}
				}
			}
		}

		// get index to maximum distance
		template <class T>
		CGPU_EXEC 
		T fcn_get_max_crv_2_ln_dist(const T& x_1, const T& y_1, const T& x_2, const T& y_2, dt_int32 ix_0, dt_int32 ix_e, Tpr<T> x, Tpr<T> y, dt_int32& ind_max)
		{
			const T m = ::sqrt(::square(x_2 - x_1) + ::square(y_2 - y_1));
			const T a = (y_2 - y_1)/m;
			const T b = -(x_2 - x_1)/m;
			const T c = (x_2*y_1 - y_2*x_1)/m;

			T d_max = 0;
			ind_max = 0;
			for(auto ik = ix_0; ik < ix_e; ik++)
			{
				const T d = ::fabs(a*x[ik] + b*y[ik] + c);
				if (d > d_max)
				{
					d_max = d;
					ind_max = ik;
				}
			}

			return d_max;
		}
	}
	
	/* Kahan summation algorithm*/
	namespace mt
	{
		// https:// en.wikipedia.org/wiki/Kahan_summation_algorithm
		template <class T>
		CGPU_EXEC 
		void fcn_kh_sum(T& sum, T val, T& error)
		{
			val = val - error;
			T t = sum + val;
			error = (t-sum) - val;
			sum = t;
		}
	}

	/* tapering/cosine */
	namespace mt
	{
		// ! calculate fermi low-pass filter alpha parameter
		template <class T>
		CGPU_EXEC
		T fcn_coef_tap(const T& r_tap, const T& r_max)
		{ 
			// y = cos(coef_tap*(x-x_tap))
			const auto c_i2pi = T(1.570796326794896619231);			// pi/2
			return (fcn_is_equal(r_tap, r_max))?T(0):c_i2pi/(r_max - r_tap);
		}			
		
		// ! calculate fermi low-pass filter alpha parameter
		template <class T>
		CGPU_EXEC
		T fcn_fermi_lpf_alpha(const T& r_cut, const T& dr_0=T(0.25), const T& fr_0=T(1e-02))
		{ 
			// y = 1/(1+exp(alpha*(x^2-x_c^2)))
			return log(T(1)/fr_0-T(1))/(::square(r_cut+dr_0)-::square(r_cut));
		}		

		// ! fermi low-pass filter value
		template <class T>
		CGPU_EXEC
		T fcn_fermi_lpf(const T& alpha, const T& gl2_max, const T& g2)
		{ 
			// y = 1/(1+exp(alpha*(x^2-x_c^2)))
			return T(1)/(T(1) + ::exp(alpha*(g2 - gl2_max)));
		}

		// ! cosine tapering
		template <class T>
		CGPU_EXEC
		T fcn_cos_tap(const T& x_tap, const T& coef_tap, const T& x)
		{
			return (x<x_tap)?T(1):(::cos(coef_tap*(x-x_tap)));
		}

		// ! apply cosine tapering to y and dy
		template <class T>
		CGPU_EXEC
		void fcn_apply_cos_tap(const T& x_tap, const T& coef_tap, const T& x, T& y, T& dy)
		{
			if (x_tap<x)
			{
				T tap, dtap;
				::sincos(coef_tap*(x-x_tap), &dtap, &tap);
				dy = dy*tap - coef_tap*y*dtap;
				y *= tap;
			}
		}
	}

	/* electron microscopy fcns */
	namespace mt
	{
		// ! exp(-r^2/(2*sigma_r^2)) -> FFT ->exp(2*pi^2*sigma_r^2*g^2)
		template <class T>
		CGPU_EXEC
		T fcn_sigma_r_2_sigma_g(const T& sigma_g)
		{
			return T(1)/(T(6.283185307179586476925)*sigma_g);
		}

		// ! plane wave exponential factor
		template <class T>
		CGPU_EXEC
		T fcn_n_2pi_sft(const T& x, const T& x_c) 
		{ 
			return -T(6.283185307179586476925)*(x-x_c);
		}

		// in: E_0(keV), Output: lambda (electron wave): Angs
		// E. J. Kirkland - Advanced computing in electron microscopy page: 10-13
		template <class T>
		CGPU_EXEC 
		T fcn_lambda(const T& E_0)
		{
			const auto emass = T(510.99906);
			const auto hc = T(12.3984244);

			return hc/::sqrt(E_0*(T(2)*emass + E_0));
		}

		// in: E_0(keV), Output: gamma(relativistic factor): unitless
		// E. J. Kirkland - Advanced computing in electron microscopy page: 10-13
		template <class T>
		CGPU_EXEC 
		T fcn_gamma(const T& E_0)
		{
			const auto emass = T(510.99906);

			return T(1) + E_0/emass;
		}

		// in: E_0(keV), Output: sigma (Interaction parameter): radians/(kV - Angs)
		// E. J. Kirkland - Advanced computing in electron microscopy page: 79
		template <class T>
		CGPU_EXEC 
		T fcn_e_interact_parm_kva(const T& E_0)
		{
			const auto c_2pi = T(6.283185307179586476925);
			const auto emass = T(510.99906);
			T x = (emass + E_0)/(T(2)*emass + E_0);

			return c_2pi*x/(fcn_lambda(E_0)*E_0);
		}

		// in: E_0(keV), Theta, Output: sigma = gamma*lambda/c_pot_factor: radians/(V - Angs), where c_pot_factor = 47.877645145863056
		// E. J. Kirkland - Advanced computing in electron microscopy page: 10-13, 79, 120

		template <class T>
		CGPU_EXEC 
		T fcn_transf_exp_factor(const T& E_0, const T& theta)
		{
			return T(1e-3)*fcn_e_interact_parm_kva(E_0)/::cos(theta);
		}

		// reciprocal Angs to radians
		// in: theta(mrad), E_0(keV), Output: Angs^-1
		template <class T>
		CGPU_EXEC 
		T fcn_rangs_2_rad(const T& E_0, const T& rAngs)
		{
			return ::asin(rAngs*fcn_lambda(E_0));
		}

		// in: theta(mrad), E_0(keV), Output: Angs^-1
		template <class T>
		CGPU_EXEC 
		T fcn_rad_2_rangs(const T& E_0, const T& theta)
		{
			return ::sin(theta)/fcn_lambda(E_0);
		}

		// Marc De Graef - Introduction to Conventional Transmission Electron Microscopy page: 608
		// hwhm: Half width at half maximum
		// sigma: standard deviation
		template <class T>
		CGPU_EXEC 
		T fcn_hwhm_2_sigma(const T& v)
		{
			const auto c_hwhm_2_sigma = T(0.84932180028801907);		// hwhm to sigma 1/(sqrt(2*log(2)))

			return v*c_hwhm_2_sigma;
		}

		// Marc De Graef - Introduction to Conventional Transmission Electron Microscopy page: 608
		// fwhm: Full width at half maximum
		// sigma: standard deviation
		template <class T>
		CGPU_EXEC 
		T fcn_fwhm_2_sigma(const T& v)
		{
			const auto c_fwhm_2_sigma = T(0.42466090014400953);		// fwhm to sigma 1/(2*sqrt(2*log(2)))

			return v*c_fwhm_2_sigma;
		}

		// Marc De Graef - Introduction to Conventional Transmission Electron Microscopy page: 608
		// iehwgd: e^-1 half-width value of the Gaussian distribution
		// sigma: standard deviation
		template <class T>
		CGPU_EXEC 
		T fcn_iehwgd_2_sigma(const T& v)
		{
			const auto c_iehwgd_2_sigma = T(0.70710678118654746);		// iehwgd to sigma 1/sqrt(2)

			return v*c_iehwgd_2_sigma;
		}

		// Marc De Graef - Introduction to Conventional Transmission Electron Microscopy page: 608
		// sigma: standard deviation
		template <class T>
		CGPU_EXEC 
		T fcn_rad_2_sigma(const T& E_0, const T& theta)
		{
			const auto q0 = fcn_rad_2_rangs(E_0, theta);
			return fcn_iehwgd_2_sigma(q0);
		}

		// E. J. Kirkland - Advanced computing in electron microscopy page: 33
		template <class T>
		CGPU_EXEC 
		T fcn_scherzer_defocus(const T& E_0, const T& c_30, T n=1.0)
		{
			n = ::fmax(n, 1);
			const auto lambda = fcn_lambda(E_0);

			return -::copysign(::sqrt((T(2)*n-T(0.5))*::fabs(c_30)*lambda), c_30);
		}

		// E. J. Kirkland - Advanced computing in electron microscopy page: 33
		template <class T>
		CGPU_EXEC 
		T fcn_scherzer_aperture(const T& E_0, const T& c_30, T n=1.0)
		{
			n = ::fmax(n, 1);
			const auto lambda = fcn_lambda(E_0);
			return ::pow(T(4)*(T(2)*n-T(0.5))*lambda/::fabs(c_30), T(0.25));
		}

		// E. J. Kirkland - Advanced computing in electron microscopy page: 33
		template <class T>
		CGPU_EXEC 
		void fcn_scherzer_conditions(const T& E_0, const T& c_30, T& defocus, T& aperture)
		{
			defocus = fcn_scherzer_defocus(E_0, c_30);

			aperture = fcn_scherzer_aperture(E_0, c_30);
		}
	}

	namespace mt
	{
		// ceil - scale size
		template <class T>
		CGPU_EXEC 
		dt_int32 fcn_sc_size_c(const dt_int32& n, T factor)
		{
			return max(dt_int32(::ceil(n*factor)), 1);
		}

		// ! next power of 2
		CGPU_EXEC
		dt_uint64 fcn_next_pow2(dt_uint32 x)
		{
			--x;
			x |= x >> 1;
			x |= x >> 2;
			x |= x >> 4;
			x |= x >> 8;
			x |= x >> 16;

			return ++x;
		}

		// from bytes to mb
		template <class T>
		CGPU_EXEC
		dt_float64 fcn_size_mb(const dt_int64& n)
		{
			return static_cast<dt_float64>(n*sizeof(T)/dt_float64(c_bytes_2_mb));
		}

		// select background option
		template <class T>
		T fcn_select_bg(const eFil_Sel_Typ& bg_opt, const T& v_min, const T& v_max, const T& v_mean, T bg =0)
		{
			T bg_r = 0;

			switch (bg_opt) 
			{
				case efst_min:
					bg_r = v_min;
				break;
				case efst_max:
					bg_r = v_max;
				break;
				case efst_mean:
					bg_r = v_mean;
				break;
				case efst_min_mean:
					bg_r = (v_mean + v_min)/T(2);
				break;
				case efst_max_mean:
					bg_r = (v_mean + v_max)/T(2);
				break;
				case efst_user_def:
					bg_r = bg;
				break;
			}

			return bg_r;
		}

		template <class T>
		dt_init_list<T> fcn_read_col(T* v, const dt_int64& n_r, const dt_int64& n_c, const dt_int64& idx, dt_int64 icol=0)
		{
			const auto ip = icol*n_r + idx;

			switch(n_c)
			{
				case 1:
				{
					return {T(v[ip + 0*n_r])};
				}
				case 2:
				{
					return {T(v[ip + 0*n_r]), T(v[ip + 1*n_r])};
				}
				case 3:
				{
					return {T(v[ip + 0*n_r]), T(v[ip + 1*n_r]), T(v[ip + 2*n_r])};
				}
				case 4:
				{
					return {T(v[ip + 0*n_r]), T(v[ip + 1*n_r]), T(v[ip + 2*n_r]), T(v[ip + 3*n_r])};
				}
				case 5:
				{
					return {T(v[ip + 0*n_r]), T(v[ip + 1*n_r]), T(v[ip + 2*n_r]), T(v[ip + 3*n_r]), T(v[ip + 4*n_r])};
				}
				case 6:
				{
					return {T(v[ip + 0*n_r]), T(v[ip + 1*n_r]), T(v[ip + 2*n_r]), T(v[ip + 3*n_r]), T(v[ip + 4*n_r]), T(v[ip + 5*n_r])};
				}
				case 7:
				{
					return {T(v[ip + 0*n_r]), T(v[ip + 1*n_r]), T(v[ip + 2*n_r]), T(v[ip + 3*n_r]), T(v[ip + 4*n_r]), T(v[ip + 5*n_r]), T(v[ip + 6*n_r])};
				}
				case 8:
				{
					return {T(v[ip + 0*n_r]), T(v[ip + 1*n_r]), T(v[ip + 2*n_r]), T(v[ip + 3*n_r]), T(v[ip + 4*n_r]), T(v[ip + 5*n_r]), T(v[ip + 6*n_r]), T(v[ip + 7*n_r])};
				}
				case 9:
				{
					return {T(v[ip + 0*n_r]), T(v[ip + 1*n_r]), T(v[ip + 2*n_r]), T(v[ip + 3*n_r]), T(v[ip + 4*n_r]), T(v[ip + 5*n_r]), T(v[ip + 6*n_r]), T(v[ip + 7*n_r]), T(v[ip + 8*n_r])};
				}
				case 10:
				{
					return {T(v[ip + 0*n_r]), T(v[ip + 1*n_r]), T(v[ip + 2*n_r]), T(v[ip + 3*n_r]), T(v[ip + 4*n_r]), T(v[ip + 5*n_r]), T(v[ip + 6*n_r]), T(v[ip + 7*n_r]), T(v[ip + 8*n_r]), T(v[ip + 9*n_r])};
				}
			}
		}

	}

	/* strings */
	namespace mt
	{
		static inline
		std::vector<std::string> fcn_str_split(const std::string& s, char delimiter)
		{
			std::vector<std::string> string_elems;
			std::istringstream tokenStream(s);
			std::string token;

			while (std::getline(tokenStream, token, delimiter))
			{
				string_elems.push_back(token);
			}
			return string_elems;
		}

		// to lowercase (in place)
		static inline
		void fcn_str_2_lower_ip(std::string& str)
		{
			std::for_each(str.begin(), str.end(), [](auto &v){ v = tolower(v); });
		}

		// to lowercase
		static inline
		std::string fcn_str_2_lower(std::string str)
		{
			fcn_str_2_lower_ip(str);
			return str;
		}

		// to uppercase (in place)
		static inline
		void fcn_str_2_upper_ip(std::string& str)
		{
			std::for_each(str.begin(), str.end(), [](auto &v){ v = toupper(v); });
		}

		// to uppercase
		static inline
		std::string fcn_str_2_upper(std::string str)
		{
			fcn_str_2_upper_ip(str);
			return str;
		}

		// trim from start (in place)
		static inline
		void fcn_str_ltrim_ip(std::string &s) {
			s.erase(s.begin(), std::find_if (s.begin(), s.end(), [](dt_int32 ch) 
			{
				return !std::isspace(ch);
			}));
		}

		// trim from start
		static inline
		std::string fcn_str_ltrim(std::string s) 
		{
			fcn_str_ltrim_ip(s);
			return s;
		}

		// trim from end (in place)
		static inline
		void fcn_str_rtrim_ip(std::string &s) 
		{
			s.erase(std::find_if(s.rbegin(), s.rend(), [](dt_int32 ch) 
			{
				return !std::isspace(ch);
			}).base(), s.end());
		}

		// trim from end
		static inline
		std::string fcn_rtrim(std::string s) 
		{
			fcn_str_rtrim_ip(s);
			return s;
		}

		// trim from both ends (in place)
		static inline
		void fcn_str_trim_ip_ip(std::string &s) 
		{
			fcn_str_ltrim_ip(s);
			fcn_str_rtrim_ip(s);
		}

		// trim from both ends
		static inline
		std::string fcn_str_trim(std::string &s) 
		{
			fcn_str_trim_ip_ip(s);
		}

		// Compare string
		static inline
		dt_bool fcn_str_cmp(const std::string &s1, const std::string &s2) 
		{
			return (s1.find(s2) != std::string::npos);
		}

		void fcn_str_rem_char_not_of(std::string &str, std::string chars) 
		{
			str.erase(remove_if (str.begin(), str.end(), [chars](const auto& v){ return std::all_of(chars.begin(), chars.end(), [v](const auto& a){return a != v; }); }), str.end() );
		}

		void fcn_str_rem_char_of(std::string &str, std::string chars) 
		{
			str.erase(remove_if (str.begin(), str.end(), [chars](const auto& v){ return std::any_of(chars.begin(), chars.end(), [v](const auto& a){return a == v; }); }), str.end() );
		}

	}

#endif