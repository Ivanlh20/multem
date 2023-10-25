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

#include <type_traits>
#include <vector>
#include <algorithm>
#include <thread>
#include <mutex>
#include <utility>
#include <functional> 

#include "const_enum.h"
#include "type_traits_gen.h"
#include "math_mt.h"
#include "r_2d.h"
#include "r_3d.h"
#include "fcns_cgpu_gen.h"
#include "vctr_cpu.h"	
#include "grid_2d.h"	

namespace mt
{	
	/********************************* 2d array of vectors *********************************/
	template <class TVctr>
	struct AV_2d
	{
		using value_type = Value_type<TVctr>;

		AV_2d():m_dim_0(0), m_dim_1(0), m_size(0), m_size_a(0), m_d_size(0) {}

		AV_2d(dt_int32 rows, dt_int32 cols, dt_int32 d_size, dt_int32 size_a=-1)
		{
				set_size(rows, cols, d_size, size_a);
		}

		void clear()
		{
			for(auto ik = 0; ik < m_size_a; ik++)
			{
				m_data[ik].clear();
			}
			m_data.clear();
		}

		void shrink_to_fit()
		{
			for(auto ik = 0; ik < m_size_a; ik++)
			{
				m_data[ik].shrink_to_fit();
			}
			m_data.shrink_to_fit();
		}

		void clear_shrink_to_fit()
		{
			clear();
			shrink_to_fit();
		}

		void set_size(dt_int32 dim_0, dt_int32 dim_1, dt_int32 d_size, dt_int32 size_a=-1)
		{
			m_dim_0 = dim_0;
			m_dim_1 = dim_1;
			m_size = m_dim_0*m_dim_1;
			m_size_a = (size_a<0)?m_size:min(m_size, size_a);
			m_d_size = max(0, d_size);

			m_data.resize(m_size);

			if (m_d_size>0)
			{
				for(auto ik = 0; ik < m_size_a; ik++)
				{	
					m_data[ik].resize(m_d_size);
				}
			}
		}

		dt_bool sub_2_ind(const dt_int32 ix_0, const dt_int32 ix_1) const 
		{ 
			return (ix_0 + m_dim_0*ix_1);
		}

		dt_bool exists(const dt_int32 ix_0, const dt_int32 ix_1) const 
		{ 
			return (sub_2_ind(ix_0, ix_1) < m_size_a);
		}

		template <class TAV_2d>
		void assign(TAV_2d &av_2d)
		{
			set_size(av_2d.dim_0(), av_2d.dim_1(), av_2d.d_size(), av_2d.size_a());

			if (m_d_size>0)
			{
				for(auto ik = 0; ik < m_size_a; ik++)
				{
					thrust::copy(av_2d[ik].begin(), av_2d[ik].end(), m_data[ik].begin());
				}
			}
		}

 		template <class TAV_2d>
		void cpy_allow_data(TAV_2d &av_2d)
		{
 			if (m_d_size>0)
			{
				for(auto ik = 0; ik < m_size_a; ik++)
				{
					thrust::copy(m_data[ik].begin(), m_data[ik].end(), av_2d[ik].begin());
				}
			}
		}

 		template <class TAV_2d>
		AV_2d<TVctr>& operator=(const TAV_2d &av_2d)
		{
			assign(av_2d);
			return *this;
		}

		CGPU_EXEC_INL
		TVctr& operator[](const dt_int32 idx) { return m_data[idx]; }

		CGPU_EXEC_INL
		const TVctr& operator[](const dt_int32 idx) const { return m_data[idx]; }

		TVctr& operator()(const dt_int32 ix_0, const dt_int32 ix_1) 
		{ 
			return m_data[sub_2_ind(ix_0, ix_1)];
		}

		const TVctr& operator()(const dt_int32 ix_0, const dt_int32 ix_1) const
		{
			return m_data[sub_2_ind(ix_0, ix_1)];
		}

		value_type& operator()(const dt_int32 ix_0, const dt_int32 ix_1, const dt_int32 id) 
		{ 
			return m_data[sub_2_ind(ix_0, ix_1)][id];
		}

		const value_type& operator()(const dt_int32 ix_0, const dt_int32 ix_1, const dt_int32 id) const 
		{ 
			return m_data[sub_2_ind(ix_0, ix_1)][id];
		}

		dt_int32 dim_0() const
		{
			return m_dim_0;
		}

		dt_int32 dim_1() const
		{
			return m_dim_1;
		}

		dt_int32 size() const
		{
			return m_size;
		}

		dt_int32 size_a() const
		{
			return m_size_a;
		}

 		dt_int32 d_size() const
		{
			return m_d_size;
		}

		void fill(value_type value)
		{
			if (m_d_size>0)
			{
				for(auto ik = 0; ik < m_size_a; ik++)
				{
					thrust::fill(m_data[ik].begin(), m_data[ik].end(), value);
				}
			}
		}

		dt_int32 m_dim_0;
		dt_int32 m_dim_1;
		dt_int32 m_size;
		dt_int32 m_size_a;
		dt_int32 m_d_size;
		Vctr_cpu<TVctr> m_data;
	};

	/***************************** atoms for simulated annealing ***************************/
	template <class T>
	class Atom_Data_Sa{
		public:
			using value_type = T;
			using size_type = dt_uint64;

			size_type size() const
			{
				return Z.size();
			}

			dt_bool empty() const
			{
				return size() == 0;
			}

			template <class TAtom_SA> 
			void assign(TAtom_SA &atom_sa)
			{ 
				Z.assign(atom_sa.Z.begin(), atom_sa.Z.end());

				r_min.assign(atom_sa.r_min.begin(), atom_sa.r_min.end());
				r_max.assign(atom_sa.r_max.begin(), atom_sa.r_max.end());
				r_0.assign(atom_sa.r_0.begin(), atom_sa.r_0.end());
				r_d.assign(atom_sa.r_d.begin(), atom_sa.r_d.end());

				r.assign(atom_sa.r.begin(), atom_sa.r.end());
				r_n.assign(atom_sa.r_n.begin(), atom_sa.r_n.end());
				r_opt.assign(atom_sa.r_opt.begin(), atom_sa.r_opt.end());

				chi2.assign(atom_sa.chi2.begin(), atom_sa.chi2.end());
				chi2_n.assign(atom_sa.chi2_n.begin(), atom_sa.chi2_n.end());
				chi2_opt.assign(atom_sa.chi2_opt.begin(), atom_sa.chi2_opt.end());

				df.assign(atom_sa.df.begin(), atom_sa.df.end());
			}

			// resize number of atoms
			void resize(const size_type& new_size, const value_type& value = value_type())
			{
				Z.resize(new_size, value);

				r_min.resize(new_size, value);
				r_max.resize(new_size, value);
				r_0.resize(new_size, value);
				r_d.resize(new_size, value);

				r.resize(new_size, value);
				r_n.resize(new_size, value);
				r_opt.resize(new_size, value);

				chi2.resize(new_size, value);
				chi2_n.resize(new_size, value);
				chi2_opt.resize(new_size, value);

				df.resize(new_size, value);

			}

			// set atoms
			void set_ptc(const size_type& natoms_i, dt_float64 *atoms_i, dt_float64 *atoms_min_i, dt_float64 *atoms_max_i)
			{
				resize(natoms_i);

				for(auto iatoms = 0; iatoms < size(); iatoms++)
				{
					Z[iatoms] = static_cast<dt_int32>(atoms_i[0*natoms_i + iatoms]);	// atomic number
					r[iatoms].x = atoms_i[1*natoms_i + iatoms];							// x-position
					r[iatoms].y = atoms_i[2*natoms_i + iatoms];							// y-position
					r[iatoms].z = atoms_i[3*natoms_i + iatoms];							// z-position

					r_min[iatoms].x = atoms_min_i[0*natoms_i + iatoms];					// x-position
					r_min[iatoms].y = atoms_min_i[1*natoms_i + iatoms];					// y-position
					r_min[iatoms].z = atoms_min_i[2*natoms_i + iatoms];					// z-position

					r_max[iatoms].x = atoms_max_i[0*natoms_i + iatoms];					// x-position
					r_max[iatoms].y = atoms_max_i[1*natoms_i + iatoms];					// y-position
					r_max[iatoms].z = atoms_max_i[2*natoms_i + iatoms];					// z-position

					r_0[iatoms] = r_min[iatoms];
					r_d[iatoms] = r_max[iatoms]-r_min[iatoms];

					df[iatoms] = 1;
				}
			}

			// set atoms
			void set_ptc(const size_type& natoms_i, dt_float64 *atoms_i, R_3d<T> d_i)
			{
				resize(natoms_i);

				for(auto iatoms = 0; iatoms < size(); iatoms++)
				{
					Z[iatoms] = static_cast<dt_int32>(atoms_i[0*natoms_i + iatoms]);	// atomic number

					r[iatoms].x = atoms_i[0*natoms_i + iatoms];							// x-position
					r[iatoms].y = atoms_i[1*natoms_i + iatoms];							// y-position
					r[iatoms].z = atoms_i[2*natoms_i + iatoms];							// z-position

					r_min[iatoms] = r - d_i;
					r_max[iatoms] = r + d_i;

					r_0[iatoms] = r_min[iatoms];
					r_d[iatoms] = r_max[iatoms]-r_min[iatoms];

					df[iatoms] = 1;
				}
			}

			void set_range(dt_int32 Z_i, R_3d<T> r_min_i, R_3d<T> r_max_i)
			{
				for(auto iatoms = 0; iatoms < size(); iatoms++)
				{
					Z[iatoms] = Z_i;

					r_min[iatoms] = r_min_i;
					r_max[iatoms] = r_max_i;
					r_0[iatoms] = r_min[iatoms];
					r_d[iatoms] = r_max[iatoms]-r_min[iatoms];
					df[iatoms] = 1;
				}
			}

			inline
			T norm_2(const dt_int32& iatoms, const R_3d<T>& r)
			{
				auto rd = r_n[iatoms]-r;
				return mt::norm_2(rd);
			}

			Vctr<dt_int32, edev_cpu> Z;

			Vctr<R_3d<T>, edev_cpu> r_min;
			Vctr<R_3d<T>, edev_cpu> r_max;
			Vctr<R_3d<T>, edev_cpu> r_0;
			Vctr<R_3d<T>, edev_cpu> r_d;

			Vctr<R_3d<T>, edev_cpu> r;
			Vctr<R_3d<T>, edev_cpu> r_n;
			Vctr<R_3d<T>, edev_cpu> r_opt;

			Vctr<T, edev_cpu> chi2;
			Vctr<T, edev_cpu> chi2_n;
			Vctr<T, edev_cpu> chi2_opt;

			Vctr<T, edev_cpu> df;
	};

	template <class T>
	struct Atom_Sa
	{
		public:
			using value_type = T;

			T x;
			T y;
			T R2_max;
			T* r2;
			T* c3;
			T* c2;
			T* c1;
			T* c0;

			dt_int32 ix_0;
			dt_int32 ix_n;
			dt_int32 iy_0;
			dt_int32 iy_n;

			dt_int32 *iv;
			T* v;

			Atom_Sa(): x(0), y(0), R2_max(0), r2(nullptr), c3(nullptr), c2(nullptr), c1(nullptr), 
			c0(nullptr), ix_0(1), ix_n(0), iy_0(0), iy_n(0), iv(nullptr), v(nullptr) {}

			inline
			void set_ix0_ixn(const Grid_2d<T>& grid_2d, const T& r_max)
			{
				fcn_get_idx_0_idx_n(x, r_max, grid_2d.drx, grid_2d.pbc_x, grid_2d.nx, ix_0, ix_n);
			}

			inline
			void set_iy0_iyn(const Grid_2d<T>& grid_2d, const T& r_max)
			{
				fcn_get_idx_0_idx_n(y, r_max, grid_2d.dry, grid_2d.pbc_y, grid_2d.ny, iy_0, iy_n);
			}

		#ifdef __CUDACC__
			inline
			D_Grid_Blk get_eval_cubic_poly_gridBT()
			{
				D_Grid_Blk d_grid_blk;
				d_grid_blk.grid = dim3((iy_n+c_thr_2d_x-1)/c_thr_2d_x, (ix_n+c_thr_2d_y-1)/c_thr_2d_y);
				d_grid_blk.blk = dim3(c_thr_2d_x, c_thr_2d_y);
				return d_grid_blk;
			}
		#endif
	};

	/************************* affine parameters for shx and scy ***************************/
	template <class T>
	struct Atp_1
	{
		R_2d<T> f;
		R_2d<T> tr;
		T chi2;

		Atp_1():f(T(0), T(1)), tr(T(0), T(0)), chi2(0) {}

		Atp_1(const R_2d<T>& f_i, const R_2d<T>& tr_i, const T& chi2_i): f(f_i), tr(tr_i), chi2(chi2_i) {}

		template <class TAtp_1> 
		Atp_1<T>& operator=(TAtp_1 &atp_1)
		{
			f = atp_1.f;
			tr = atp_1.tr;
			chi2 = atp_1.chi2;

			return *this;
		}
	};

	/*************************** affine parameters for rotation ****************************/
	template <class T>
	struct Atp_2
	{
		T theta;
		R_2d<T> tr;
		T chi2;

		Atp_2():theta(T(0)), tr(T(0), T(0)), chi2(0) {}

		Atp_2(const T& theta_i, const R_2d<T>& tr_i, const T& chi2_i): theta(theta_i), tr(tr_i), chi2(chi2_i) {}

		template <class TAtp_2> 
		Atp_2<T>& operator=(TAtp_2 &atp_2)
		{
			theta = atp_2.theta;
			tr = atp_2.tr;
			chi2 = atp_2.chi2;

			return *this;
		}
	};
}