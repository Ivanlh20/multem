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

#ifndef FCNS_ELEM_H
	#define FCNS_ELEM_H

	#ifdef _MSC_VER
		#pragma once
	#endif

	#include "const_enum.cuh"
	#include "math.cuh"
	#include "r_2d.cuh"
	#include "r_3d.cuh"

	/***************************************************************************************/
	/********************************* base class definitions ******************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T, eFcn_typ Fcn_typ> class Fcn_Elem;
	}

	/***************************************************************************************/
	/********************************** forward definitions ********************************/
	/***************************************************************************************/
	namespace mt
	{
		template<class T>
		using Fcn_Cos_Tap = Fcn_Elem<T, efcn_cos_tap>;

		template<class T>
		using Fcn_Gauss = Fcn_Elem<T, efcn_gauss>;

		template<class T>
		using Fcn_Exp = Fcn_Elem<T, efcn_exp>;

		template<class T>
		using Fcn_Fermi = Fcn_Elem<T, efcn_fermi>;

		template<class T>
		using Fcn_Butwth = Fcn_Elem<T, efcn_butwth>;

		template<class T>
		using Fcn_Hann = Fcn_Elem<T, efcn_hann>;
	}

	/***************************************************************************************/
	/******************************** cosine tapering function *****************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T>
		class Fcn_Elem<T, efcn_cos_tap>
		{
		public:
			using value_type = T;

			T r_tap;
			T r_max;
			T coef_tap;

			/************************************* constructors ************************************/
			CGPU_EXEC
			Fcn_Elem(): r_tap(0), r_max(0), coef_tap(0){}

			Fcn_Elem(const T& r_tap, const T& r_max)
			{
				set_in_data(r_tap, r_max);
			}

			/* copy constructor */
			CGPU_EXEC
			Fcn_Elem(const Fcn_Elem<T, efcn_cos_tap>& parm)
			{
				*this = parm;
			}

			/* converting constructor */
			template <class U>
			CGPU_EXEC
			Fcn_Elem(const Fcn_Elem<U, efcn_cos_tap>& parm)
			{
				*this = parm;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			Fcn_Elem<T, efcn_cos_tap>& operator=(const Fcn_Elem<T, efcn_cos_tap>& parm)
			{
				if (this != &parm)
				{
					r_tap = parm.r_tap;
					r_max = parm.r_max;
					coef_tap = parm.coef_tap;
				}

				return *this;
			}
			
			/* converting assignment operator */
			template <class U>
			CGPU_EXEC
			Fcn_Elem<T, efcn_cos_tap>& operator=(const Fcn_Elem<U, efcn_cos_tap>& parm)
			{
				r_tap = T(parm.r_tap);
				r_max = T(parm.r_max);
				coef_tap = T(parm.coef_tap);

				return *this;
			}

			template <class U>
			CGPU_EXEC
			void assign(const Fcn_Elem<U, efcn_cos_tap>& parm)
			{
				*this = parm;
			}

			/***************************************************************************************/
			void set_in_data(const T& r_tap, const T& r_max)
			{
				this->r_tap = r_tap;
				this->r_max = r_max;

				coef_tap = fcn_coef_tap(r_tap, r_max);
			}

			CGPU_EXEC
			void clear()
			{
				r_tap = T(0);
				r_max = T(0);
				coef_tap = T(0);
			}

			/***************************************************************************************/
			CGPU_EXEC
			T eval_r(const T& r, const T& r_tap, const T& r_max) const
			{ 
				const T coef_tap = fcn_coef_tap(r_tap, r_max);
				return fcn_cos_tap(r_tap, coef_tap, r);
			}

			CGPU_EXEC
			T eval_r(const T& r) const
			{ 
				return fcn_cos_tap(r_tap, coef_tap, r);
			}

			CGPU_EXEC
			R_2d<T> eval_r(const R_2d<T>& r) const
			{ 
				return {eval_r(r.x), eval_r(r.y)};
			}

			CGPU_EXEC
			R_3d<T> eval_r(const R_3d<T>& r) const
			{ 
				return {eval_r(r.x), eval_r(r.y), eval_r(r.z)};
			}

			/***************************************************************************************/
			CGPU_EXEC
			T operator()(const T& r, const T& r_tap, const T& r_max) const
			{ 
				return eval_r(r, r_tap, r_max);
			}

			CGPU_EXEC
			T operator()(const T& r) const
			{ 
				return eval_r(r);
			}

			CGPU_EXEC
			R_2d<T> operator()(const R_2d<T>& r) const
			{ 
				return eval_r(r);
			}

			CGPU_EXEC
			R_3d<T> operator()(const R_3d<T>& r) const
			{ 
				return eval_r(r);
			}
		};
	}

	/***************************************************************************************/
	/********************************** gaussian function **********************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T>
		class Fcn_Elem<T, efcn_gauss>
		{
		public:
			using value_type = T;

			T a_1;
			T b_1;

			/************************************* constructors ************************************/
			CGPU_EXEC
			Fcn_Elem(): a_1(0), b_1(0) {}

			Fcn_Elem(const T& sigma)
			{
				set_in_data(T(1), sigma);
			}

			Fcn_Elem(const T& a, const T& sigma)
			{
				set_in_data(a, sigma);
			}

			/* copy constructor */
			CGPU_EXEC
			Fcn_Elem(const Fcn_Elem<T, efcn_gauss>& parm)
			{
				*this = parm;
			}

			/* copy assignment operator */
			CGPU_EXEC
			Fcn_Elem<T, efcn_gauss>& operator=(const Fcn_Elem<T, efcn_gauss>& parm)
			{
				if (this != &parm)
				{
					a_1 = parm.a_1;
					b_1 = parm.b_1;
				}

				return *this;
			}

			/* converting assignment operator */
			template <class U>
			CGPU_EXEC
			Fcn_Elem<T, efcn_gauss>& operator=(const Fcn_Elem<U, efcn_gauss>& parm)
			{
				a_1 = T(parm.a_1);
				b_1 = T(parm.b_1);

				return *this;
			}

			template <class U>
			CGPU_EXEC
			void assign(const Fcn_Elem<U, efcn_gauss>& parm)
			{
				*this = parm;
			}

			/***************************************************************************************/
			void set_in_data(const T& a, const T& sigma)
			{
				a_1 = a;
				b_1 = T(1)/(T(2)*sigma*sigma);
			}

			CGPU_EXEC
			void clear()
			{
				a_1 = T(0);
				b_1 = T(0);
			}

			/***************************************************************************************/
			CGPU_EXEC
			T eval_r2(const T& r2, const T& sigma) const
			{ 
				return a_1*exp(-r2/(T(2)*sigma*sigma));
			}

			CGPU_EXEC
			T eval_r2(const T& r2) const
			{ 
				return a_1*exp(-b_1*r2);
			}

			CGPU_EXEC
			R_2d<T> eval_r2(const R_2d<T>& r2) const
			{ 
				return {eval_r2(r2.x), eval_r2(r2.y)};
			}

			CGPU_EXEC
			R_3d<T> eval_r2(const R_3d<T>& r2) const
			{ 
				return {eval_r2(r2.x), eval_r2(r2.y), eval_r2(r2.z)};
			}

			/***************************************************************************************/
			CGPU_EXEC
			T eval_r(const T& r, const T& sigma) const
			{ 
				return eval_r2(::square(r), sigma);
			}

			CGPU_EXEC
			T eval_r(const T& r) const
			{ 
				return eval_r2(::square(r));
			}

			CGPU_EXEC
			R_2d<T> eval_r(const R_2d<T>& r) const
			{ 
				return eval_r2(square(r));
			}

			CGPU_EXEC
			R_3d<T> eval_r(const R_3d<T>& r) const
			{ 
				return eval_r2(square(r));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T operator()(const T& r2, const T& sigma) const
			{ 
				return eval_r2(r2, sigma);
			}

			CGPU_EXEC
			T operator()(const T& r2) const
			{ 
				return eval_r2(r2);
			}

			CGPU_EXEC
			R_2d<T> operator()(const R_2d<T>& r2) const
			{ 
				return eval_r2(r2);
			}

			CGPU_EXEC
			R_3d<T> operator()(const R_3d<T>& r2) const
			{ 
				return eval_r2(r2);
			}
		};
	}

	/***************************************************************************************/
	/********************************* exponential function ********************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T>
		class Fcn_Elem<T, efcn_exp>
		{
		public:
			using value_type = T;

			T a_1;
			T b_1;

			/************************************* constructors ************************************/
			CGPU_EXEC
			Fcn_Elem(): a_1(0), b_1(0) {}

			Fcn_Elem(const T& a, const T& beta)
			{
				set_in_data(a, beta);
			}

			/* copy constructor */
			CGPU_EXEC
			Fcn_Elem(const Fcn_Elem<T, efcn_exp>& parm)
			{
				*this = parm;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			Fcn_Elem<T, efcn_exp>& operator=(const Fcn_Elem<T, efcn_exp>& parm)
			{
				if (this != &parm)
				{
					a_1 = parm.a_1;
					b_1 = parm.b_1;
				}

				return *this;
			}

			/***************************************************************************************/
			void set_in_data(const T& a, const T& beta)
			{
				a_1 = a;
				b_1 = T(1)/(beta);
			}

			CGPU_EXEC
			void clear()
			{
				a_1 = T(0);
				b_1 = T(0);
			}

			/***************************************************************************************/
			CGPU_EXEC
			T eval_r(const T& r, const T& beta) const
			{ 
				return a_1*exp(-r/beta);
			}

			CGPU_EXEC
			T eval_r(const T& r) const
			{ 
				return a_1*exp(-b_1*r);
			}

			CGPU_EXEC
			R_2d<T> eval_r(const R_2d<T>& r) const
			{ 
				return {eval_r(r.x), eval_r(r.y)};
			}

			CGPU_EXEC
			R_3d<T> eval_r(const R_3d<T>& r) const
			{ 
				return {eval_r(r.x), eval_r(r.y), eval_r(r.z)};
			}

			/***************************************************************************************/
			CGPU_EXEC
			T eval_r2(const T& r2, const T& beta) const
			{ 
				return eval_r(::sqrt(r2), beta);
			}

			CGPU_EXEC
			T eval_r2(const T& r2) const
			{ 
				return eval_r(::sqrt(r2));
			}

			CGPU_EXEC
			R_2d<T> eval_r2(const R_2d<T>& r2) const
			{ 
				return eval_r(::sqrt(r2));
			}

			CGPU_EXEC
			R_3d<T> eval_r2(const R_3d<T>& r2) const
			{ 
				return eval_r(::sqrt(r2));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T operator()(const T& r, const T& beta) const
			{ 
				return eval_r(r, beta);
			}

			CGPU_EXEC
			T operator()(const T& r) const
			{ 
				return eval_r(r);
			}

			CGPU_EXEC
			R_2d<T> operator()(const R_2d<T>& r) const
			{ 
				return eval_r(r);
			}

			CGPU_EXEC
			R_3d<T> operator()(const R_3d<T>& r) const
			{ 
				return eval_r(r);
			}
		};
	}
	
	/***************************************************************************************/
	/************************************ fermi function ***********************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T>
		class Fcn_Elem<T, efcn_fermi>
		{
		public:
			using value_type = T;

			T a_1;
			T b_1;
			T b_2;

			/************************************* constructors ************************************/
			CGPU_EXEC
			Fcn_Elem(): a_1(0), b_1(0), b_2(0) {}

			Fcn_Elem(const T& a, const T& alpha, const T& radius)
			{
				set_in_data(a, alpha, radius);
			}

			/* copy constructor */
			CGPU_EXEC
			Fcn_Elem(const Fcn_Elem<T, efcn_fermi>& parm)
			{
				*this = parm;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			Fcn_Elem<T, efcn_fermi>& operator=(const Fcn_Elem<T, efcn_fermi>& parm)
			{
				if (this != &parm)
				{
					a_1 = parm.a_1;
					b_1 = parm.b_1;
					b_2 = parm.b_2;
				}

				return *this;
			}

			/***************************************************************************************/
			void set_in_data(const T& a, const T& alpha, const T& radius)
			{
				a_1 = a;
				b_1 = alpha;
				b_2 = ::square(radius);
			}

			CGPU_EXEC
			void clear()
			{
				a_1 = T(0);
				b_1 = T(0);
				b_2 = T(0);
			}

			/***************************************************************************************/
			CGPU_EXEC
			T eval_r2(const T& r2, const T& radius) const
			{ 
				return a_1/(T(1) + exp(b_1*(r2-::square(radius))));
			}

			CGPU_EXEC
			T eval_r2(const T& r2) const
			{ 
				return a_1/(T(1) + exp(b_1*(r2-b_2)));
			}

			CGPU_EXEC
			R_2d<T> eval_r2(const R_2d<T>& r2) const
			{ 
				return {eval_r2(r2.x), eval_r2(r2.y)};
			}

			CGPU_EXEC
			R_3d<T> eval_r2(const R_3d<T>& r2) const
			{ 
				return {eval_r2(r2.x), eval_r2(r2.y), eval_r2(r2.z)};
			}

			/***************************************************************************************/
			CGPU_EXEC
			T eval_r(const T& r, const T& radius) const
			{ 
				return eval_r2(::square(r), radius);
			}

			CGPU_EXEC
			T eval_r(const T& r) const
			{ 
				return eval_r2(::square(r));
			}

			CGPU_EXEC
			R_2d<T> eval_r(const R_2d<T>& r) const
			{ 
				return eval_r2(square(r));
			}

			CGPU_EXEC
			R_3d<T> eval_r(const R_3d<T>& r) const
			{ 
				return eval_r2(square(r));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T operator()(const T& r2, const T& radius) const
			{ 
				return eval_r2(r2, radius);
			}

			CGPU_EXEC
			T operator()(const T& r2) const
			{ 
				return eval_r2(r2);
			}

			CGPU_EXEC
			R_2d<T> operator()(const R_2d<T>& r2) const
			{ 
				return eval_r2(r2);
			}

			CGPU_EXEC
			R_3d<T> operator()(const R_3d<T>& r2) const
			{ 
				return eval_r2(r2);
			}
		};
	}
	
	/***************************************************************************************/
	/********************************** butterworth function *******************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T>
		class Fcn_Elem<T, efcn_butwth>
		{
		public:
			using value_type = T;

			T a_1;
			dt_int32 b_1;
			T b_2;

			/************************************* constructors ************************************/
			CGPU_EXEC
			Fcn_Elem(): a_1(0), b_1(0), b_2(0) {}

			Fcn_Elem(const T& a, const dt_int32& n, const T& radius)
			{
				set_in_data(a, n, radius);
			}

			/* copy constructor */
			CGPU_EXEC
			Fcn_Elem(const Fcn_Elem<T, efcn_butwth>& parm)
			{
				*this = parm;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			Fcn_Elem<T, efcn_butwth>& operator=(const Fcn_Elem<T, efcn_butwth>& parm)
			{
				if (this != &parm)
				{
					a_1 = parm.a_1;
					b_1 = parm.b_1;
					b_2 = parm.b_2;
				}

				return *this;
			}

			/***************************************************************************************/
			void set_in_data(const T& a, const dt_int32& n, const T& radius)
			{
				a_1 = a;
				b_1 = n;
				b_2 = ::square(radius);
			}

			CGPU_EXEC
			void clear()
			{
				a_1 = T(0);
				b_1 = dt_int32(0);
				b_2 = T(0);
			}

			/***************************************************************************************/
			CGPU_EXEC
			T eval_r2(const T& r2, const dt_int32& n, const T& radius) const
			{ 
				return a_1/(T(1)+pow(r2/::square(radius), n));
			}

			CGPU_EXEC
			T eval_r2(const T& r2) const
			{ 
				return a_1/(T(1)+pow(r2/b_2, b_1));
			}

			CGPU_EXEC
			R_2d<T> eval_r2(const R_2d<T>& r2) const
			{ 
				return {eval_r2(r2.x), eval_r2(r2.y)};
			}

			CGPU_EXEC
			R_3d<T> eval_r2(const R_3d<T>& r2) const
			{ 
				return {eval_r2(r2.x), eval_r2(r2.y), eval_r2(r2.z)};
			}

			/***************************************************************************************/
			CGPU_EXEC
			T eval_r(const T& r, const dt_int32& n, const T& radius) const
			{ 
				return eval_r2(::square(r), n, radius);
			}

			CGPU_EXEC
			T eval_r(const T& r) const
			{ 
				return eval_r2(::square(r));
			}

			CGPU_EXEC
			R_2d<T> eval_r(const R_2d<T>& r) const
			{ 
				return eval_r2(square(r));
			}

			CGPU_EXEC
			R_3d<T> eval_r(const R_3d<T>& r) const
			{ 
				return eval_r2(square(r));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T operator()(const T& r2, const dt_int32& n, const T& radius) const
			{ 
				return eval_r2(r2, n, radius);
			}

			CGPU_EXEC
			T operator()(const T& r2) const
			{ 
				return eval_r2(r2);
			}

			CGPU_EXEC
			R_2d<T> operator()(const R_2d<T>& r2) const
			{ 
				return eval_r2(r2);
			}

			CGPU_EXEC
			R_3d<T> operator()(const R_3d<T>& r2) const
			{ 
				return eval_r2(r2);
			}
		};
	}
	
	/***************************************************************************************/
	/************************************ hann function ************************************/
	/***************************************************************************************/
	namespace mt
	{
		// https:// en.wikipedia.org/wiki/Hann_function
		template <class T>
		class Fcn_Elem<T, efcn_hann>
		{
		public:
			using value_type = T;

			T a_1;
			T b_1;

			/************************************* constructors ************************************/
			CGPU_EXEC
			Fcn_Elem(): a_1(0), b_1(0) {}

			Fcn_Elem(const T& a, const T& l)
			{
				set_in_data(a, l);
			}

			/* copy constructor */
			CGPU_EXEC
			Fcn_Elem(const Fcn_Elem<T, efcn_hann>& parm)
			{
				*this = parm;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			Fcn_Elem<T, efcn_hann>& operator=(const Fcn_Elem<T, efcn_hann>& parm)
			{
				if (this != &parm)
				{
					a_1 = parm.a_1;
					b_1 = parm.b_1;
				}

				return *this;
			}

			/***************************************************************************************/
			void set_in_data(const T& a, const T& l)
			{
				a_1 = a;
				b_1 = c_pi<T>/l;
			}

			CGPU_EXEC
			void clear()
			{
				a_1 = T(0);
				b_1 = T(0);
			}

			/***************************************************************************************/
			CGPU_EXEC
			T eval_r(const T& r, const T& l) const
			{ 
				return a_1*::square(cos(c_pi<T>*r/l)); // r<=l/2, otherwise 0
			}

			CGPU_EXEC
			T eval_r(const T& r) const
			{ 
				return a_1*::square(cos(b_1*r)); // r<=l/2, otherwise 0
			}

			CGPU_EXEC
			R_2d<T> eval_r(const R_2d<T>& r) const
			{ 
				return {eval_r(r.x), eval_r(r.y)};
			}

			CGPU_EXEC
			R_3d<T> eval_r(const R_3d<T>& r) const
			{ 
				return {eval_r(r.x), eval_r(r.y), eval_r(r.z)};
			}

			/***************************************************************************************/
			CGPU_EXEC
			T eval_r2(const T& r2, const T& l) const
			{ 
				return eval_r(::sqrt(r2), l);
			}

			CGPU_EXEC
			T eval_r2(const T& r2) const
			{ 
				return eval_r(::sqrt(r2));
			}

			CGPU_EXEC
			R_2d<T> eval_r2(const R_2d<T>& r2) const
			{ 
				return eval_r(::sqrt(r2));
			}

			CGPU_EXEC
			R_3d<T> eval_r2(const R_3d<T>& r2) const
			{ 
				return eval_r(::sqrt(r2));
			}

			/***************************************************************************************/
			CGPU_EXEC
			T operator()(const T& r, const T& l) const
			{ 
				return eval_r(r, l);
			}

			CGPU_EXEC
			T operator()(const T& r) const
			{ 
				return eval_r(r);
			}

			CGPU_EXEC
			R_2d<T> operator()(const R_2d<T>& r) const
			{ 
				return eval_r(r);
			}

			CGPU_EXEC
			R_3d<T> operator()(const R_3d<T>& r) const
			{ 
				return eval_r(r);
			}
		};
	}

#endif 