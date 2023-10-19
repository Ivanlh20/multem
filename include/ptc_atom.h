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

#pragma once

#include "macros.h"
#include "math_mt.h"
#include "r_3d.h"
#include "pvctr.h"
#include "vctr_cpu.h"
#include "fcns_cpu.h"
#include "fcns_cgpu_gen.h"
#include "ptc_r_3d.h"

/* template definition */
namespace mt
{
	template <class T> class Ptc_s_Atom;

	template <class T> class Ptc_Atom;
}

/* template Ptc_s_Atom */
namespace mt
{
	template <class T>
	class Ptc_s_Atom: public Ptc_3d_0<T>
	{
	public:
		dt_int32 Z;
		dt_float32 sigma;
		dt_float32 occ;
		dt_int32 tag;
		dt_int32 charge;

		Ptc_s_Atom():Z(0), Ptc_s_3d_0<T>(), sigma(0), occ(0), tag(0), charge(0) {};

		Ptc_s_Atom(const dt_int32& Z, const T& x, const T& y, const T& z, 
		dt_float32 sigma=c_dflt_rms3d, dt_float32 occ=c_dflt_occ, dt_int32 tag=c_dflt_tag, dt_int32 charge=c_dflt_charge);
			
		Ptc_s_Atom(const dt_int32& Z, const R_3d<T>& r, 
		dt_float32 sigma=c_dflt_rms3d, dt_float32 occ=c_dflt_occ, dt_int32 tag=c_dflt_tag, dt_int32 charge=c_dflt_charge);

		/* constructor by pointer */
		template <class U>
		Ptc_s_Atom(U* v, const dt_int64& n_r, const dt_int64& n_c, const dt_int64& idx, dt_int64 icol=0);

		/* copy constructor */
		Ptc_s_Atom(const Ptc_s_Atom<T>& ptc_s);

		/* converting constructor */
		template <class U> 
		Ptc_s_Atom(const Ptc_s_Atom<U>& ptc_s);

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		Ptc_s_Atom<T>& operator=(const Ptc_s_Atom<T>& ptc_s);

		/* converting assignment operator */
		template <class U> 
		Ptc_s_Atom<T>& operator=(const Ptc_s_Atom<U>& ptc_s);

		template <class U> 
		void assign(const Ptc_s_Atom<U>& ptc_s);
	};
}

/* template Ptc_Atom */
namespace mt
{
	template <class T>
	class Ptc_Atom: public Ptc_R_3d<T>
	{
		public:
			using value_type = T;
			using size_type = dt_int64;

			using Ptc_s = Ptc_s_Atom<T>;	

			mutable Vctr_cpu<dt_int32> Z;			// atomic number
			mutable Vctr_cpu<dt_float32> sigma;		// 3d root fcn_mean squared displacement (rmsd)
			mutable Vctr_cpu<dt_float32> occ;		// occupancy
			mutable Vctr_cpu<dt_int32> tag;			// tag
			mutable Vctr_cpu<dt_int32> charge;		// charge

			dt_int32 cols_used;						// number of used columns

			R_2d<dt_int32> Z_lim;
			R_2d<dt_float32> sigma_lim;
			R_2d<dt_float32> occ_lim;
			R_2d<dt_int32> tag_lim;
			R_2d<dt_int32> charge_lim;

			/************************************* constructors ************************************/
			Ptc_Atom();

			template <class U>
			Ptc_Atom(const pVctr_cpu_64<U>& ptc, const R_3d<U>& bs, dt_bool pbc_xy = false, dt_bool b_statistic = true);

			/* copy constructor */
			Ptc_Atom(const Ptc_Atom<T>& ptc);

			/* converting constructor */
			template <class U>
			Ptc_Atom(const Ptc_Atom<U>& ptc);

			/******************************** assignment operators *********************************/
			/* assignment operator */
			Ptc_Atom<T>& operator=(const Ptc_Atom<T>& ptc);

			/* converting assignment operator */
			template <class U> 
			Ptc_Atom<T>& operator=(const Ptc_Atom<U>& ptc);

			template <class U>
			void assign(const Ptc_Atom<U>& ptc);

			/***************************************************************************************/
			virtual size_type cols() const;		

			void clear();

			void resize(size_type new_size);

			void reserve(size_type new_size);

			void shrink_to_fit();

			void push_back(const Ptc_s& ptc_s);

			dt_float32 get_sigma(const size_type& ia) const;

			dt_float32 get_occ(const size_type& ia) const;

			dt_int32 get_tag(const size_type& ia) const;

			dt_int32 get_charge(const size_type& ia) const;

			Ptc_s get(const size_type& ia) const;

			void set(const size_type& ia, const Ptc_s& ptc_s);

			template <class U>
			void set_ptc(const Ptc_Atom<U>& ptc, dt_bool pbc_xy = false, dt_bool b_statistic = true);

			template <class U>
			void set_ptc(const pVctr_cpu_64<U>& ptc, const R_3d<U>& bs, dt_bool pbc_xy = false, dt_bool b_statistic = true);

			/* copy data to pointer */
			template <class U>
			dt_int32 cpy_to_ptr(U *ptc, size_type n_ptc, dt_int32 is_0=0, dt_int32 is_e=8);

			// sort by x
			void sort_by_x();

			// sort by y
			void sort_by_y();

			// sort by z
			void sort_by_z();

			// sort by idx
			void sort_by_idx(const dt_int32 idx = 3);

			/***************************************************************************************/
			virtual void get_statistic();

			// max z value within a tag
			void minmax_z_by_region(const T& tag_v, T& z_min, T& z_max);

			void remove_ptc_out_z_range(const T& z_0, const T& z_e);
	};
}

/* fcns */
namespace mt
{
	template <class T>
	void remove_ptc_out_z_range(const T& z_0, const T& z_e, Ptc_Atom<T>& ptc);
}

#include "../src/ptc_atom.inl"