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
 * GNU General Public License for more details.Gauss_wd_
 *
 * You should have received a copy of the GNU General Public License
 * along with Multem. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef FCNS_WD_H
	#define FCNS_WD_H

	#ifdef _MSC_VER
		#pragma once
	#endif

	#include "const_enum.cuh"
	#include "cgpu_fcns_gen.cuh"
	#include "cgpu_stream.cuh"
	#include "r_2d.cuh"
	#include "r_3d.cuh"
	#include "fcns_elem.cuh"
	#include "cgpu_detail.cuh"

	/***************************************************************************************/
	/********************************** base class fcn window ******************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T, eDim Dim, eFcn_typ Fcn_typ>
		class Wdb_fcn_xd: public Fcn_Elem<T, Fcn_typ>
		{
		public:
			using value_type = T;

			R_xd<T, Dim> r;
			T r_wd_2;
			T r_max_2;
			T sft;
			T sc;

			/************************************* constructors ************************************/
			CGPU_EXEC
			Wdb_fcn_xd(): Fcn_Elem<T, Fcn_typ>(), r(), r_wd_2(0), r_max_2(0), sft(0), sc(1) {}

			/* copy constructor */
			CGPU_EXEC
			Wdb_fcn_xd(const Wdb_fcn_xd<T, Dim, Fcn_typ>& wd)
			{
				*this = wd;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC 
			Wdb_fcn_xd<T, Dim, Fcn_typ>& operator=(const Wdb_fcn_xd<T, Dim, Fcn_typ>& wd)
			{
				if (this != &wd)
				{
					Fcn_Elem<T, Fcn_typ>::operator=(wd);

					r = wd.r;
					r_wd_2 = wd.r_wd_2;
					r_max_2 = wd.r_max_2;
					sft = wd.sft;
					sc = wd.sc;
				}

				return *this;
			}

			CGPU_EXEC
			void assign(const Wdb_fcn_xd<T, Dim, Fcn_typ>& wd)
			{
				*this = wd;
			}

			/***************************************************************************************/
			CGPU_EXEC
			void clear()
			{
				Fcn_Elem<T, Fcn_typ>::clear();

				r = T(0);
				r_wd_2 = T(0);
				r_max_2 = T(0);
				sft = T(0);
				sc = T(1);
			}

			CGPU_EXEC
			T eval_r2(const T& r2) const
			{ 
				return (r2<r_wd_2)?T(1):(r2<r_max_2)?::fmax(T(0), (Fcn_Elem<T, Fcn_typ>::eval_r2(r2) - sft)/sc):T(0);
			}

			CGPU_EXEC
			T eval_r2(const R_2d<T>& r2) const
			{ 
				return eval_r2(r2.x)*eval_r2(r2.y);
			}

			CGPU_EXEC
			T eval_r2(const R_3d<T>& r2) const
			{ 
				return eval_r2(r2.x)*eval_r2(r2.y)*eval_r2(r2.z);
			}

			void set_in_data(const R_xd<T, Dim>& r, const T& r_wd, const T& r_max)
			{
				this->r = r;
				r_wd_2 = ::square(r_wd);
				r_max_2 = ::square(r_max);
				sft = Fcn_Elem<T, Fcn_typ>::eval_r2(r_max_2);
				sc = Fcn_Elem<T, Fcn_typ>::eval_r2(r_wd_2) - sft;
			}
		};
	}

	/***************************************************************************************/
	/********************************** forward definitions ********************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T, eDim Dim, eFcn_typ Fcn_typ> class Wd_fcn_xd;

		/* Gauss */
		template<class T, eDim Dim>
		using Wd_Gauss_xd = Wd_fcn_xd<T, Dim, efcn_gauss>;

		template<class T>
		using Wd_Gauss_1d = Wd_fcn_xd<T, edim_1, efcn_gauss>;

		template<class T>
		using Wd_Gauss_2d = Wd_fcn_xd<T, edim_2, efcn_gauss>;

		template<class T>
		using Wd_Gauss_3d = Wd_fcn_xd<T, edim_3, efcn_gauss>;

		/* Exp */
		template<class T, eDim Dim>
		using Wd_Exp_xd = Wd_fcn_xd<T, Dim, efcn_exp>;

		template<class T>
		using Wd_Exp_1d = Wd_fcn_xd<T, edim_1, efcn_exp>;

		template<class T>
		using Wd_Exp_2d = Wd_fcn_xd<T, edim_2, efcn_exp>;

		template<class T>
		using Wd_Exp_3d = Wd_fcn_xd<T, edim_3, efcn_exp>;

		/* Fermi */
		template<class T, eDim Dim>
		using Wd_Fermi_xd = Wd_fcn_xd<T, Dim, efcn_fermi>;

		template<class T>
		using Wd_Fermi_1d = Wd_fcn_xd<T, edim_1, efcn_fermi>;

		template<class T>
		using Wd_Fermi_2d = Wd_fcn_xd<T, edim_2, efcn_fermi>;

		template<class T>
		using Wd_Fermi_3d = Wd_fcn_xd<T, edim_3, efcn_fermi>;

		/* Butwth */
		template<class T, eDim Dim>
		using Wd_Butwth_xd = Wd_fcn_xd<T, Dim, efcn_butwth>;

		template<class T>
		using Wd_Butwth_1d = Wd_fcn_xd<T, edim_1, efcn_butwth>;

		template<class T>
		using Wd_Butwth_2d = Wd_fcn_xd<T, edim_2, efcn_butwth>;

		template<class T>
		using Wd_Butwth_3d = Wd_fcn_xd<T, edim_3, efcn_butwth>;

		/* Hann */
		template<class T, eDim Dim>
		using Wd_Hann_xd = Wd_fcn_xd<T, Dim, efcn_hann>;

		template<class T>
		using Wd_Hann_1d = Wd_fcn_xd<T, edim_1, efcn_hann>;

		template<class T>
		using Wd_Hann_2d = Wd_fcn_xd<T, edim_2, efcn_hann>;

		template<class T>
		using Wd_Hann_3d = Wd_fcn_xd<T, edim_3, efcn_hann>;
	}

	/***************************************************************************************/
	/********************************** gaussian window ************************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T, eDim Dim>
		class Wd_fcn_xd<T, Dim, efcn_gauss>: public Wdb_fcn_xd<T, Dim, efcn_gauss>
		{
		public:
			using value_type = T;

			/************************************* constructors ************************************/
			CGPU_EXEC
			Wd_fcn_xd(): Wdb_fcn_xd<T, Dim, efcn_gauss>() {}

			Wd_fcn_xd(const R_xd<T, Dim>& r, const T& sigma, const T& r_wd, const T& r_max) 
			{
				set_in_data(r, sigma, r_wd, r_max);
			}

			/* copy constructor */
			CGPU_EXEC
			Wd_fcn_xd(const Wd_fcn_xd<T, Dim, efcn_gauss>& wd): Wdb_fcn_xd<T, Dim, efcn_gauss>(wd){}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC 
			Wdb_fcn_xd<T, Dim, efcn_gauss>& operator=(const Wdb_fcn_xd<T, Dim, efcn_gauss>& wd)
			{
				Wdb_fcn_xd<T, Dim, efcn_gauss>::operator=(wd);

				return *this;
			}

			/***************************************************************************************/
			void set_in_data(const R_xd<T, Dim>& r, const T& sigma, const T& r_wd, const T& r_max)
			{
				Fcn_Elem<T, efcn_gauss>::set_in_data(T(1), sigma);
				Wdb_fcn_xd<T, Dim, efcn_gauss>::set_in_data(r, r_wd, r_max);
			}

		};
	}

	/***************************************************************************************/
	/********************************* exponential window **********************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T, eDim Dim>
		class Wd_fcn_xd<T, Dim, efcn_exp>: public Wdb_fcn_xd<T, Dim, efcn_exp>
		{
		public:
			using value_type = T;

			/************************************* constructors ************************************/
			CGPU_EXEC
			Wd_fcn_xd(): Wdb_fcn_xd<T, Dim, efcn_exp>() {}

			Wd_fcn_xd(const R_xd<T, Dim>& r, const T& beta, const T& r_wd, const T& r_max) 
			{
				set_in_data(r, beta, r_wd, r_max);
			}

			/* copy constructor */
			CGPU_EXEC
			Wd_fcn_xd(const Wd_fcn_xd<T, Dim, efcn_exp>& wd): Wdb_fcn_xd<T, Dim, efcn_exp>(wd){}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC 
			Wdb_fcn_xd<T, Dim, efcn_exp>& operator=(const Wdb_fcn_xd<T, Dim, efcn_exp>& wd)
			{
				Wdb_fcn_xd<T, Dim, efcn_exp>::operator=(wd);

				return *this;
			}

			/***************************************************************************************/
			void set_in_data(const R_xd<T, Dim>& r, const T& beta, const T& r_wd, const T& r_max)
			{
				Fcn_Elem<T, efcn_exp>::set_in_data(T(1), beta);
				Wdb_fcn_xd<T, Dim, efcn_exp>::set_in_data(r, r_wd, r_max);
			}

		};
	}

	/***************************************************************************************/
	/************************************ fermi window *************************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T, eDim Dim>
		class Wd_fcn_xd<T, Dim, efcn_fermi>: public Wdb_fcn_xd<T, Dim, efcn_fermi>
		{
		public:
			using value_type = T;

			/************************************* constructors ************************************/
			CGPU_EXEC
			Wd_fcn_xd(): Wdb_fcn_xd<T, Dim, efcn_fermi>() {}

			Wd_fcn_xd(const R_xd<T, Dim>& r, const T& alpha, const T& r_wd, const T& r_max) 
			{
				set_in_data(r, alpha, r_wd, r_max);
			}

			/* copy constructor */
			CGPU_EXEC
			Wd_fcn_xd(const Wd_fcn_xd<T, Dim, efcn_fermi>& wd): Wdb_fcn_xd<T, Dim, efcn_fermi>(wd){}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC 
			Wdb_fcn_xd<T, Dim, efcn_fermi>& operator=(const Wdb_fcn_xd<T, Dim, efcn_fermi>& wd)
			{
				Wdb_fcn_xd<T, Dim, efcn_fermi>::operator=(wd);

				return *this;
			}

			/***************************************************************************************/
			void set_in_data(const R_xd<T, Dim>& r, const T& alpha, const T& r_wd, const T& r_max)
			{
				Fcn_Elem<T, efcn_fermi>::set_in_data(T(1), alpha, r_wd);
				Wdb_fcn_xd<T, Dim, efcn_fermi>::set_in_data(r, r_wd, r_max);
			}

		};
	}

	/***************************************************************************************/
	/********************************* butterworth window **********************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T, eDim Dim>
		class Wd_fcn_xd<T, Dim, efcn_butwth>: public Wdb_fcn_xd<T, Dim, efcn_butwth>
		{
		public:
			using value_type = T;

			/************************************* constructors ************************************/
			CGPU_EXEC
			Wd_fcn_xd(): Wdb_fcn_xd<T, Dim, efcn_butwth>() {}

			Wd_fcn_xd(const R_xd<T, Dim>& r, const dt_int32& n, const T& r_wd, const T& r_max) 
			{
				set_in_data(r, n, r_wd, r_max);
			}

			/* copy constructor */
			CGPU_EXEC
			Wd_fcn_xd(const Wd_fcn_xd<T, Dim, efcn_butwth>& wd): Wdb_fcn_xd<T, Dim, efcn_butwth>(wd){}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC 
			Wdb_fcn_xd<T, Dim, efcn_butwth>& operator=(const Wdb_fcn_xd<T, Dim, efcn_butwth>& wd)
			{
				Wdb_fcn_xd<T, Dim, efcn_butwth>::operator=(wd);

				return *this;
			}

			/***************************************************************************************/
			void set_in_data(const R_xd<T, Dim>& r, const dt_int32& n, const T& r_wd, const T& r_max)
			{
				Fcn_Elem<T, efcn_butwth>::set_in_data(T(1), n, r_wd);
				Wdb_fcn_xd<T, Dim, efcn_butwth>::set_in_data(r, r_wd, r_max);
			}

		};
	}

	/***************************************************************************************/
	/************************************* hann window *************************************/
	/***************************************************************************************/
	namespace mt
	{
		template <class T, eDim Dim>
		class Wd_fcn_xd<T, Dim, efcn_hann>: public Wdb_fcn_xd<T, Dim, efcn_hann>
		{
		public:
			using value_type = T;

			/************************************* constructors ************************************/
			CGPU_EXEC
			Wd_fcn_xd(): Wdb_fcn_xd<T, Dim, efcn_hann>() {}

			Wd_fcn_xd(const R_xd<T, Dim>& r, const T& l, const T& r_wd, const T& r_max) 
			{
				set_in_data(r, l, r_wd, r_max);
			}

			/* copy constructor */
			CGPU_EXEC
			Wd_fcn_xd(const Wd_fcn_xd<T, Dim, efcn_hann>& wd): Wdb_fcn_xd<T, Dim, efcn_hann>(wd){}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC 
			Wdb_fcn_xd<T, Dim, efcn_hann>& operator=(const Wdb_fcn_xd<T, Dim, efcn_hann>& wd)
			{
				Wdb_fcn_xd<T, Dim, efcn_hann>::operator=(wd);

				return *this;
			}

			/***************************************************************************************/
			void set_in_data(const R_xd<T, Dim>& r, const T& l, const T& r_wd, const T& r_max)
			{
				Fcn_Elem<T, efcn_hann>::set_in_data(T(1), l);
				Wdb_fcn_xd<T, Dim, efcn_hann>::set_in_data(r, r_wd, r_max);
			}

		};
	}

#endif 