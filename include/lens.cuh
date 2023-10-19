/*
 * This file is part of Multem.
 * Copyright 2022 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef LENS_H
	#define LENS_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include "const_enum_mt.cuh"
	#include "math_mt.h"
	#include "fcns_cgpu_gen.h"
	#include "grid.h"

	namespace mt
	{
		/********************************* Lens Class **********************************
		 * This class represents an optical lens in electron microscopy, capturing 
		 * various aberration coefficients and other lens properties.
		 * 
		 * Template Argument:
		 * - T: The data type for all numerical values (usually float or double).
		 * 
		 * Members:
		 * - Aberration Coefficients: Various types of aberrations are considered, each
		 *   represented by its own coefficient (c_10, c_12, ..., c_56) and its azimuthal
		 *   angle (phi_12, phi_21, ..., phi_56).
		 * - Aperture Angles: Inner and outer aperture angles of the lens.
		 * - Temporal and Spatial Incidence: Parameters for Gaussian and exponential 
		 *   components of temporal and spatial incidence.
		 * - Defocus: Type and plane of defocus.
		 * - Wavelength: Wavelength of the electron beam.
		 * - Precomputed Coefficients: Coefficients that are precomputed for efficiency.
		 * - Source Sampling: Information related to the sampling of the source.
		 ******************************************************************************/

		template <class T>
		class Lens
		{
		public:
			using value_type = T;

			dt_int32 m;						// Momentum of the vortex

			T c_10;							// Defocus (\AA)
			T c_12;							// 2-fold astigmatism (\AA)
			T phi_12;						// Azimuthal angle of 2-fold astigmatism (rad)

			T c_21;							// Axial coma (\AA)
			T phi_21;						// Azimuthal angle of axial coma (rad)
			T c_23;							// 3-fold astigmatism (\AA)
			T phi_23;						// Azimuthal angle of 3-fold astigmatism (rad)

			T c_30;							// 3rd order spherical aberration (\AA)
			T c_32;							// Axial star aberration (\AA)
			T phi_32;						// Azimuthal angle of axial star aberration (rad)
			T c_34;							// 4-fold astigmatism (\AA)
			T phi_34;						// Azimuthal angle of 4-fold astigmatism (rad)

			T c_41;							// 4th order axial coma (\AA)
			T phi_41;						// Azimuthal angle of 4th order axial coma (rad)
			T c_43;							// 3-lobe aberration (\AA)
			T phi_43;						// Azimuthal angle of 3-lobe aberration (rad)
			T c_45;							// 5-fold astigmatism (\AA)
			T phi_45;						// Azimuthal angle of 5-fold astigmatism (rad)

			T c_50;							// 5th order spherical aberration (\AA)
			T c_52;							// 5th order axial star aberration (\AA)
			T phi_52;						// Azimuthal angle of 5th order axial star aberration (rad)
			T c_54;							// 5th order rosette aberration (\AA)
			T phi_54;						// Azimuthal angle of 5th order rosette aberration (rad)
			T c_56;							// 6-fold astigmatism (\AA)
			T phi_56;						// Azimuthal angle of 6-fold astigmatism (rad)

			T inner_aper_ang;				// Inner aperture (rad);
			T outer_aper_ang;				// Outer aperture (rad);

			T tp_inc_a;						// Height proportion of a normalized Gaussian [0, 1]
			T tp_inc_sigma;					// standard deviation of the defocus spread function for the Gaussian component:
			T tp_inc_beta;					// standard deviation of the defocus spread function for the exponential component:
			dt_int32 tp_inc_npts;			// number of integration points of the defocus spread function

			T tp_inc_iehwgd;				// e^-1 half-width value of the Gaussian distribution

			T spt_inc_a;					// Height proportion of a normalized Gaussian [0, 1]
			T spt_inc_sigma;				// standard deviation of the source spread function for the Gaussian component: For parallel ilumination(\AA^-1);otherwise (\AA)
			T spt_inc_beta;					// standard deviation of the source spread function for the exponential component: For parallel ilumination(\AA^-1);otherwise (\AA)
			dt_int32 spt_inc_rad_npts;		// number of radial integration points
			dt_int32 spt_inc_azm_npts;		// number of azimuth integration points

			T spt_inc_iehwgd;				// e^-1 half-width value of the Gaussian distribution
			T spt_inc_theta_c;				// divergence semi-angle (rad)

			eZero_Def_Typ zero_def_typ;// Defocus type: ezdt_first = 1, ezdt_middle = 2, ezdt_last = 3, ezdt_user_def = 4
			T zero_def_plane;				// plane

			T lambda;						// wavelength(Angstrom)

			T c_c_10;						// -pi*c_10*lambda
			T c_c_12;						// -pi*c_12*lambda

			T c_c_21;						// -2*pi*c_21*lambda^2/3
			T c_c_23;						// -2*pi*c_23*lambda^2/3

			T c_c_30;						// -pi*c_30*lambda^3/2
			T c_c_32;						// -pi*c_32*lambda^3/2
			T c_c_34;						// -pi*c_34*lambda^3/2

			T c_c_41;						// -2*pi*c_41*lambda^4/5 
			T c_c_43;						// -2*pi*c_43*lambda^4/5
			T c_c_45;						// -2*pi*c_45*lambda^4/5

			T c_c_50;						// -pi*c_50*lambda^5/3
			T c_c_52;						// -pi*c_52*lambda^5/3
			T c_c_54;						// -pi*c_54*lambda^5/3
			T c_c_56;						// -pi*c_56*lambda^5/3

			T g2_inner;						// inner_aper_ang/lambda
			T g2_outer;						// outer_aper_ang/lambda

			dt_int32 ngxs;					// number of source sampling points x
			dt_int32 ngys;					// number of source sampling points y
			T dgxs;							// source sampling m_size;
			T dgys;							// source sampling m_size;
			T g2_maxs;						// q maximum square;

			/************************************* constructors ************************************/
			Lens(): m(0), c_10(0), c_12(0), phi_12(0), c_21(0), phi_21(0), c_23(0), phi_23(0), 
				c_30(0), c_32(0), phi_32(0), c_34(0), phi_34(0), 
				c_41(0), phi_41(0), c_43(0), phi_43(0), c_45(0), phi_45(0), 
				c_50(0), c_52(0), phi_52(0), c_54(0), phi_54(0), c_56(0), phi_56(0), 
				inner_aper_ang(0), outer_aper_ang(0), tp_inc_a(1.0), tp_inc_sigma(0), tp_inc_beta(0), tp_inc_npts(0), tp_inc_iehwgd(0), 
				spt_inc_a(1.0), spt_inc_sigma(0), spt_inc_beta(0), spt_inc_rad_npts(0), spt_inc_azm_npts(0), spt_inc_iehwgd(0), spt_inc_theta_c(0), 
				zero_def_typ(ezdt_last), zero_def_plane(0), lambda(0), g2_inner(0), g2_outer(0), ngxs(0), ngys(0), dgxs(0), dgys(0), g2_maxs(0), 
				c_c_10(0), c_c_12(0), c_c_21(0), c_c_23(0), c_c_30(0), c_c_32(0), c_c_34(0), c_c_41(0), 
				c_c_43(0), c_c_45(0), c_c_50(0), c_c_52(0), c_c_54(0), c_c_56(0) {}

			/* copy constructor */
			CGPU_EXEC
			Lens(const Lens<T>& lens): Lens()
			{
				*this = lens;
			}

			/* converting constructor */
			template <class U>
			CGPU_EXEC
			Lens(const Lens<U>& lens): Lens()
			{
				*this = lens;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			Lens<T>& operator=(const Lens<T>& lens)
			{
				if (this != &lens)
				{
					m = lens.m;

					c_10 = lens.c_10;
					c_12 = lens.c_12;
					phi_12 = lens.phi_12;

					c_21 = lens.c_21;
					phi_21 = lens.phi_21;
					c_23 = lens.c_23;
					phi_23 = lens.phi_23;

					c_30 = lens.c_30;
					c_32 = lens.c_32;
					phi_32 = lens.phi_32;
					c_34 = lens.c_34;
					phi_34 = lens.phi_34;

					c_41 = lens.c_41;
					phi_41 = lens.phi_41;
					c_43 = lens.c_43;
					phi_43 = lens.phi_43;
					c_45 = lens.c_45;
					phi_45 = lens.phi_45;

					c_50 = lens.c_50;
					c_52 = lens.c_52;
					phi_52 = lens.phi_52;
					c_54 = lens.c_54;
					phi_54 = lens.phi_54;
					c_56 = lens.c_56;
					phi_56 = lens.phi_56;

					inner_aper_ang = lens.inner_aper_ang;
					outer_aper_ang = lens.outer_aper_ang;

					tp_inc_a = lens.tp_inc_a;
					tp_inc_sigma = lens.tp_inc_sigma;
					tp_inc_beta = lens.tp_inc_beta;
					tp_inc_npts = lens.tp_inc_npts;
					tp_inc_iehwgd = lens.tp_inc_iehwgd;

					spt_inc_a = lens.spt_inc_a;
					spt_inc_sigma = lens.spt_inc_sigma;
					spt_inc_beta = lens.spt_inc_beta;
					spt_inc_rad_npts = lens.spt_inc_rad_npts;
					spt_inc_azm_npts = lens.spt_inc_azm_npts;

					spt_inc_iehwgd = lens.spt_inc_iehwgd;
					spt_inc_theta_c = lens.spt_inc_theta_c;

					zero_def_typ = lens.zero_def_typ;
					zero_def_plane = lens.zero_def_plane;

					lambda = lens.lambda;

					c_c_10 = lens.c_c_10;
					c_c_12 = lens.c_c_12;

					c_c_21 = lens.c_c_21;
					c_c_23 = lens.c_c_23;

					c_c_30 = lens.c_c_30;
					c_c_32 = lens.c_c_32;
					c_c_34 = lens.c_c_34;

					c_c_41 = lens.c_c_41;
					c_c_43 = lens.c_c_43;
					c_c_45 = lens.c_c_45;

					c_c_50 = lens.c_c_50;
					c_c_52 = lens.c_c_52;
					c_c_54 = lens.c_c_54;
					c_c_56 = lens.c_c_56;

					g2_inner = lens.g2_inner;
					g2_outer = lens.g2_outer;

					ngxs = lens.ngxs;
					ngys = lens.ngys;

					dgxs = lens.dgxs;
					dgys = lens.dgys;
					g2_maxs = lens.g2_outers;
				}

				return *this;
			}

			/* converting assignment operator */
			template <class U>
			CGPU_EXEC
			Lens<T>& operator=(const Lens<U>& lens)
			{
				m = lens.m;

				c_10 = T(lens.c_10);
				c_12 = T(lens.c_12);
				phi_12 = T(lens.phi_12);

				c_21 = T(lens.c_21);
				phi_21 = T(lens.phi_21);
				c_23 = T(lens.c_23);
				phi_23 = T(lens.phi_23);

				c_30 = T(lens.c_30);
				c_32 = T(lens.c_32);
				phi_32 = T(lens.phi_32);
				c_34 = T(lens.c_34);
				phi_34 = T(lens.phi_34);

				c_41 = T(lens.c_41);
				phi_41 = T(lens.phi_41);
				c_43 = T(lens.c_43);
				phi_43 = T(lens.phi_43);
				c_45 = T(lens.c_45);
				phi_45 = T(lens.phi_45);

				c_50 = T(lens.c_50);
				c_52 = T(lens.c_52);
				phi_52 = T(lens.phi_52);
				c_54 = T(lens.c_54);
				phi_54 = T(lens.phi_54);
				c_56 = T(lens.c_56);
				phi_56 = T(lens.phi_56);

				inner_aper_ang = T(lens.inner_aper_ang);
				outer_aper_ang = T(lens.outer_aper_ang);

				tp_inc_a = T(lens.tp_inc_a);
				tp_inc_sigma = T(lens.tp_inc_sigma);
				tp_inc_beta = T(lens.tp_inc_beta);
				tp_inc_npts = lens.tp_inc_npts;
				tp_inc_iehwgd = lens.tp_inc_iehwgd;

				spt_inc_a = T(lens.spt_inc_a);
				spt_inc_sigma = T(lens.spt_inc_sigma);
				spt_inc_beta = T(lens.spt_inc_beta);
				spt_inc_rad_npts = lens.spt_inc_rad_npts;
				spt_inc_azm_npts = lens.spt_inc_azm_npts;

				spt_inc_iehwgd = T(lens.spt_inc_iehwgd);
				spt_inc_theta_c = T(lens.spt_inc_theta_c);

				zero_def_typ = lens.zero_def_typ;
				zero_def_plane = T(lens.zero_def_plane);

				lambda = T(lens.lambda);

				c_c_10 = T(lens.c_c_10);
				c_c_12 = T(lens.c_c_12);

				c_c_21 = T(lens.c_c_21);
				c_c_23 = T(lens.c_c_23);

				c_c_30 = T(lens.c_c_30);
				c_c_32 = T(lens.c_c_32);
				c_c_34 = T(lens.c_c_34);

				c_c_41 = T(lens.c_c_41);
				c_c_43 = T(lens.c_c_43);
				c_c_45 = T(lens.c_c_45);

				c_c_50 = T(lens.c_c_50);
				c_c_52 = T(lens.c_c_52);
				c_c_54 = T(lens.c_c_54);
				c_c_56 = T(lens.c_c_56);

				g2_inner = T(lens.g2_inner);
				g2_outer = T(lens.g2_outer);

				ngxs = lens.ngxs;
				ngys = lens.ngys;

				dgxs = T(lens.dgxs);
				dgys = T(lens.dgys);
				g2_maxs = T(lens.g2_outers);

				return *this;
			}

			template <class U>
			void assign(const Lens<U>& lens)
			{
				*this = lens;
			}

			/***************************************************************************************/
			void set_dep_var(const T& E_0, const Grid_2d<T>& grid)
			{
				lambda = fcn_lambda(E_0);

				c_c_10 = fcn_is_zero(c_10)?T(0):-c_pi<T>*c_10*lambda;
				c_c_12 = fcn_is_zero(c_12)?T(0):-c_pi<T>*c_12*lambda;

				c_c_21 = fcn_is_zero(c_21)?T(0):-T(2)*c_pi<T>*c_21*pow(lambda, 2)/T(3);
				c_c_23 = fcn_is_zero(c_23)?T(0):-T(2)*c_pi<T>*c_23*pow(lambda, 2)/T(3);

				c_c_30 = fcn_is_zero(c_30)?T(0):-c_pi<T>*c_30*pow(lambda, 3)/T(2);
				c_c_32 = fcn_is_zero(c_32)?T(0):-c_pi<T>*c_32*pow(lambda, 3)/T(2);
				c_c_34 = fcn_is_zero(c_34)?T(0):-c_pi<T>*c_34*pow(lambda, 3)/T(2);

				c_c_41 = fcn_is_zero(c_41)?T(0):-T(2)*c_pi<T>*c_41*pow(lambda, 4)/T(5);
				c_c_43 = fcn_is_zero(c_43)?T(0):-T(2)*c_pi<T>*c_43*pow(lambda, 4)/T(5);
				c_c_45 = fcn_is_zero(c_45)?T(0):-T(2)*c_pi<T>*c_45*pow(lambda, 4)/T(5);

				c_c_50 = fcn_is_zero(c_50)?T(0):-c_pi<T>*c_50*pow(lambda, 5)/T(3);
				c_c_52 = fcn_is_zero(c_52)?T(0):-c_pi<T>*c_52*pow(lambda, 5)/T(3);
				c_c_54 = fcn_is_zero(c_54)?T(0):-c_pi<T>*c_54*pow(lambda, 5)/T(3);
				c_c_56 = fcn_is_zero(c_56)?T(0):-c_pi<T>*c_56*pow(lambda, 5)/T(3);

				g2_inner = (inner_aper_ang < epsilon_abs<T>())?T(0):pow(rad_2_rangs(inner_aper_ang), 2);
				g2_outer = (outer_aper_ang < epsilon_abs<T>())?T(2)*grid.g2_max():pow(rad_2_rangs(outer_aper_ang), 2);

				tp_inc_a = max(T(0), min(T(1), tp_inc_a));
				set_tp_inc_sigma(tp_inc_sigma);
				tp_inc_beta = max(T(0), tp_inc_beta);
				tp_inc_npts = max(1, tp_inc_npts);

				spt_inc_a = max(T(0), min(T(1), spt_inc_a));
				set_spt_inc_sigma(spt_inc_sigma);
				spt_inc_beta = max(T(0), spt_inc_beta);
				spt_inc_rad_npts = max(1, spt_inc_rad_npts);
				spt_inc_azm_npts = max(1, spt_inc_azm_npts);

				T gmaxs = T(3.5)*spt_inc_sigma;
				g2_maxs = gmaxs*gmaxs;
				T dgs = gmaxs/static_cast<T>(spt_inc_rad_npts);

				dt_int32 n = (dgs<grid.dgx)?static_cast<dt_int32>(::floor(grid.dgx/dgs)+1):1;
				ngxs = static_cast<dt_int32>(::floor(n*gmaxs/grid.dgx)) + 1;
				dgxs = gmaxs/ngxs;

				n = (dgs<grid.dgy)?static_cast<dt_int32>(::floor(grid.dgy/dgs)+1):1;
				ngys = static_cast<dt_int32>(::floor(n*gmaxs/grid.dgy)) + 1;
				dgys = gmaxs/ngys;
			}

			CGPU_EXEC 
			T rad_2_rangs(const T& theta) const 
			{
				return sin(theta)/lambda;
			}

			void set_tp_inc_sigma(const T& ti_sigma)
			{
				tp_inc_sigma = max(T(0), ti_sigma);
				tp_inc_iehwgd = c_2i2<T>*tp_inc_sigma;
			}

			void set_spt_inc_sigma(const T& si_sigma)
			{
				spt_inc_sigma = max(T(0), si_sigma);
				spt_inc_iehwgd = c_2i2<T>*spt_inc_sigma;
				spt_inc_theta_c = asin(spt_inc_iehwgd*lambda);
			}

			void set_defocus(const T& c_10)
			{
				this->c_10 = c_10;
				c_c_10 = (fcn_is_zero(c_10))?T(0):-c_pi<T>*c_10*lambda;
			}

			T get_zero_def_plane(const T& z_min, const T& z_max) const 
			{
				switch(zero_def_typ)
				{
					case ezdt_first:
					{
						return z_min;
					}
					case ezdt_middle:
					{
						return T(0.5)*(z_min + z_max);
					}
					case ezdt_last:
					{
						return z_max;
					}
					default:
					{
						return zero_def_plane;
					}
				}
			};

			dt_bool is_zdt_first() const 
			{
				return mt::is_zdt_first(zero_def_typ);
			}

			dt_bool is_zdt_middle() const 
			{
				return mt::is_zdt_first(zero_def_typ);
			}

			dt_bool is_zdt_last() const 
			{
				return mt::is_zdt_first(zero_def_typ);
			}

			dt_bool is_zdt_user_def() const 
			{
				return mt::is_zdt_first(zero_def_typ);
			}

			T gxs(const dt_int32& ix) const 
			{ 
				return static_cast<T>(ix)*dgxs;
			}
			T gys(const dt_int32& iy) const 
			{ 
				return static_cast<T>(iy)*dgys;
			}
			T g2s(const dt_int32& ix, const dt_int32& iy) const 
			{ 
				T gxi = gxs(ix);
				T gyi = gys(iy);

				return gxi*gxi + gyi*gyi;
			}

			/***************************************************************************************/
			CGPU_EXEC_INL
			T tp_inc_sigma_2() const
			{
				return ::square(tp_inc_sigma);
			}

			CGPU_EXEC_INL
			T tp_inc_beta_2() const
			{
				return ::square(tp_inc_beta);
			}

			CGPU_EXEC_INL
			T spt_inc_sigma_2() const
			{
				return ::square(spt_inc_sigma);
			}

			CGPU_EXEC_INL
			T spt_inc_beta_2() const
			{
				return ::square(spt_inc_beta);
			}

			/***************************************************************************************/
			CGPU_EXEC_INL
			dt_bool is_phi_required() const
			{
				auto bb = fcn_is_nzero(m)||fcn_is_nzero(c_12)||fcn_is_nzero(c_21)||fcn_is_nzero(c_23)||fcn_is_nzero(c_32)||fcn_is_nzero(c_34);
				bb = bb||fcn_is_nzero(c_41)||fcn_is_nzero(c_43)||fcn_is_nzero(c_45)||fcn_is_nzero(c_52)||fcn_is_nzero(c_54)||fcn_is_nzero(c_56);

				return bb;
			}

			CGPU_EXEC_INL
			T eval_m(const T& phi) const
			{
				return T(m)*phi;
			}

			/***************************************************************************************/
			CGPU_EXEC_INL
			T eval_c_10(const T& g2) const
			{
				return c_c_10*g2;
			}

			CGPU_EXEC_INL
			T eval_c_12(const T& g2, const T& phi) const
			{
				return c_c_12*g2*sin(T(2)*(phi-phi_12));
			}

			/***************************************************************************************/
			CGPU_EXEC_INL
			T eval_c_21(const T& g3, const T& phi) const
			{
				return c_c_21*g3*sin(phi-phi_21);
			}

			CGPU_EXEC_INL
			T eval_c_23(const T& g3, const T& phi) const
			{
				return c_c_23*g3*sin(T(3)*(phi-phi_23));
			}

			CGPU_EXEC_INL
			T eval_c_21_c_23(const T& g3, const T& phi) const
			{
				return eval_c_21(g3, phi) + eval_c_23(g3, phi);
			}

			/***************************************************************************************/
			CGPU_EXEC_INL
			T eval_c_30(const T& g4) const
			{
				return c_c_30*g4;
			}

			CGPU_EXEC_INL
			T eval_c_32(const T& g4, const T& phi) const
			{
				return c_c_32*g4*sin(T(2)*(phi-phi_32));
			}

			CGPU_EXEC_INL
			T eval_c_34(const T& g4, const T& phi) const
			{
				return c_c_34*g4*sin(T(4)*(phi-phi_34));
			}

			CGPU_EXEC_INL
			T eval_c_32_c_34(const T& g4, const T& phi) const
			{
				return eval_c_32(g4, phi) + eval_c_34(g4, phi);
			}

			/***************************************************************************************/
			CGPU_EXEC_INL
			T eval_c_41(const T& g5, const T& phi) const
			{
				return c_c_41*g5*sin(phi-phi_41);
			}

			CGPU_EXEC_INL
			T eval_c_43(const T& g5, const T& phi) const
			{
				return c_c_43*g5*sin(T(3)*(phi-phi_43));
			}

			CGPU_EXEC_INL
			T eval_c_45(const T& g5, const T& phi) const
			{
				return c_c_45*g5*sin(T(5)*(phi-phi_45));
			}

			CGPU_EXEC_INL
			T eval_c_41_c_43_c_45(const T& g5, const T& phi) const
			{
				return eval_c_41(g5, phi) + eval_c_43(g5, phi) + eval_c_45(g5, phi);
			}

			/***************************************************************************************/
			CGPU_EXEC_INL
			T eval_c_50(const T& g6) const
			{
				return c_c_50*g6;
			}

			CGPU_EXEC_INL
			T eval_c_52(const T& g6, const T& phi) const
			{
				return c_c_52*g6*sin(T(2)*(phi-phi_52));
			}

			CGPU_EXEC_INL
			T eval_c_54(const T& g6, const T& phi) const
			{
				return c_c_54*g6*sin(T(4)*(phi-phi_54));
			}

			CGPU_EXEC_INL
			T eval_c_56(const T& g6, const T& phi) const
			{
				return c_c_56*g6*sin(T(6)*(phi-phi_56));
			}

			CGPU_EXEC_INL
			T eval_c_52_c_54_c_56(const T& g6, const T& phi) const
			{
				return eval_c_52(g6, phi) + eval_c_54(g6, phi) + eval_c_56(g6, phi);
			}
		};
	}

#endif