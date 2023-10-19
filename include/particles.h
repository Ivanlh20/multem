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

#include <numeric>
#include <algorithm>

#include "macros.h"
#include "math_mt.h"
#include "type_traits_gen.h"
#include "mx_2x2.h"
#include "mx_3x3.h"
#include "cgpu_vctr.cuh"
#include "fcns_cpu.h"
#ifdef __CUDACC__
	#include "fcns_gpu.h"
#endif
#include "grid.h"
#include "range.cuh"
#include "fcns_elem.cuh"
	
#include <thrust/sort.h>

/***************************************************************************************/
namespace mt
{
	namespace ptc_detail
	{
		template <class PTC_i, class PTC_o>
		void set_ptc_pbc_xy(const PTC_i& ptc_i, const dt_bool& pbc_xy, const dt_bool& b_statistic, PTC_o& ptc_o)
		{
			using T = typename PTC_o::value_type;

			ptc_o.bs = ptc_i.bs;

			if (ptc_i.size()==0) 
				return;

			if (!pbc_xy)
			{
				ptc_o = ptc_i;
			}
			else
			{
				auto bs_e = ptc_o.bs - c_dflt_pos_ee;

				for(dt_int64 ia = 0; ia < ptc_i.size(); ia++)
				{
					auto bb_x = fcn_chk_bound(ptc_i.x[ia], T(0), bs_e.x);
					auto bb_y = fcn_chk_bound(ptc_i.y[ia], T(0), bs_e.y);

					if (bb_x && bb_y)
					{
						ptc_o.push_back(ptc_i.get(ia));
					}
				}

				ptc_o.shrink_to_fit();
			}

			if (b_statistic)
			{
				ptc_o.get_statistic();
			}
		}

		template <class U, class R_D, class PTC_o>
		void set_ptc_pbc_xy(const pVctr_cpu_64<U>& p_ptc, const dt_int64& icol, const R_D& bs, const dt_bool& pbc_xy, const dt_bool& b_statistic, PTC_o& ptc_o)
		{
			using T = typename PTC_o::value_type;
			using Ptc_s = typename PTC_o::Ptc_s;

			ptc_o.bs = bs;

			if (p_ptc.empty()) 
				return;

			auto s_0 = p_ptc.s0();
			auto s_1 = p_ptc.s1();

			if (!pbc_xy)
			{
				for(dt_int64 ia = 0; ia < s_0; ia++)
				{
					ptc_o.push_back(Ptc_s(p_ptc.data(), s_0, s_1, ia, icol));
				}
			}
			else
			{
				auto bs_e = ptc_o.bs - c_dflt_pos_ee;

				for(dt_int64 ia = 0; ia < s_0; ia++)
				{
					auto ptc_s = Ptc_s(p_ptc.data(), s_0, s_1, ia, icol);
					auto bb_x = fcn_chk_bound(ptc_s.x, T(0), bs_e.x);
					auto bb_y = fcn_chk_bound(ptc_s.y, T(0), bs_e.y);

					if (bb_x && bb_y)
					{
						ptc_o.push_back(ptc_s);
					}
				}
			}

			ptc_o.shrink_to_fit();

			if (b_statistic)
			{
				ptc_o.get_statistic();
			}
		}

		/***************************************************************************************/
		/*********************************** macros 2d/3d **************************************/
		/***************************************************************************************/
		#define CASE_SORT(N) case N: { thrust::sort(first, last, mt::cgpu_fctr::less_soa<N>()); } break;

		#define ZIPITER_0_2D this->x, this->y
		#define ZIPITER_0_3D this->x, this->y, this->z

		#define ZIPITER_0(DIM) ZIPITER_0_##DIM##D
		#define ZIPITER_1(DIM) ZIPITER_0(DIM), this->c_1
		#define ZIPITER_2(DIM) ZIPITER_1(DIM), this->c_2
		#define ZIPITER_3(DIM) ZIPITER_2(DIM), this->c_3
		#define ZIPITER_4(DIM) ZIPITER_3(DIM), this->c_4
		#define ZIPITER_5(DIM) ZIPITER_4(DIM), this->c_5
		#define ZIPITER_6(DIM) ZIPITER_5(DIM), this->c_6
		#define ZIPITER_7(DIM) ZIPITER_6(DIM), this->c_7
		#define ZIPITER_8(DIM) ZIPITER_7(DIM), this->c_8
		#define ZIPITER_9(DIM) ZIPITER_8(DIM), this->c_9

		#define SWITCH_0_2D(macro) macro(0) macro(1)
		#define SWITCH_1_2D(macro) SWITCH_0_2D(macro) macro(2)
		#define SWITCH_2_2D(macro) SWITCH_1_2D(macro) macro(3)
		#define SWITCH_3_2D(macro) SWITCH_2_2D(macro) macro(4)
		#define SWITCH_4_2D(macro) SWITCH_3_2D(macro) macro(5)
		#define SWITCH_5_2D(macro) SWITCH_4_2D(macro) macro(6)
		#define SWITCH_6_2D(macro) SWITCH_5_2D(macro) macro(7)
		#define SWITCH_7_2D(macro) SWITCH_6_2D(macro) macro(8)
		#define SWITCH_8_2D(macro) SWITCH_7_2D(macro) macro(9)
		#define SWITCH_9_2D(macro) SWITCH_8_2D(macro) macro(10)

		#define SWITCH_0_3D(macro) macro(0) macro(1) macro(2)
		#define SWITCH_1_3D(macro) SWITCH_0_3D(macro) macro(3)
		#define SWITCH_2_3D(macro) SWITCH_1_3D(macro) macro(4)
		#define SWITCH_3_3D(macro) SWITCH_2_3D(macro) macro(5)
		#define SWITCH_4_3D(macro) SWITCH_3_3D(macro) macro(6)
		#define SWITCH_5_3D(macro) SWITCH_4_3D(macro) macro(7)
		#define SWITCH_6_3D(macro) SWITCH_5_3D(macro) macro(8)
		#define SWITCH_7_3D(macro) SWITCH_6_3D(macro) macro(9)
		#define SWITCH_8_3D(macro) SWITCH_7_3D(macro) macro(10)
		#define SWITCH_9_3D(macro) SWITCH_8_3D(macro) macro(11)

		#define FCN_SORT_BY_IDX(DIM, N)									\
		virtual void sort_by_idx(const dt_int32& idx)					\
		{																\
			auto first = fcn_mkzipiter_begin(ZIPITER_##N(DIM));			\
			auto last = fcn_mkzipiter_end(ZIPITER_##N(DIM));			\
																		\
			switch(idx)													\
			{															\
				SWITCH_##N##_##DIM##D(CASE_SORT)						\
			}															\
		}

		#define CD_PASTE(pre, x, y) pre ## _ ## x ## d ## _ ## y
		#define CD_EVAL_PASTE(pre, x, y) CD_PASTE(pre, x, y)
		#define CD_EVAL_2_PASTE(pre, x, y) CD_EVAL_PASTE(pre, x, y)

		#define CD_PTC_S_N(DIM, N) CD_EVAL_PASTE(Ptc_s, DIM, N)
		#define CD_PTC_S_NB(DIM, N) CD_EVAL_2_PASTE(Ptc_s, DIM, DEC(N, 1))

		#define CD_PTC_N(DIM, N) CD_EVAL_PASTE(Ptc, DIM, N)
		#define CD_PTC_NB(DIM, N) CD_EVAL_2_PASTE(Ptc, DIM, DEC(N, 1))

		#define N_DIM(N, DIM) EVAL_2_PASTE(INC, DEC(DIM, 1), N)

		/***************************************************************************************/
		#define C_PTC_S_COEF_DIM_N(DIM, N)																																\
		template <class T>																																				\
		class CD_PTC_S_N(DIM, N): public CD_PTC_S_NB(DIM, N)<T>																											\
		{																																								\
		public:																																							\
			T c_##N;																																					\
																																										\
			CD_PTC_S_N(DIM, N)(): CD_PTC_S_NB(DIM, N)<T>(), c_##N(0) {}																									\
																																										\
			/*	constructor by initializer list */ 																														\
			template <class U>																																			\
			CD_PTC_S_N(DIM, N)(const dt_init_list<U>& list): CD_PTC_S_NB(DIM, N)<T>(list)																				\
			{																																							\
				auto ptr = list.begin();																																\
																																										\
				c_##N = T(ptr[N_DIM(N, DIM)]);																															\
			}																																							\
																																										\
			/* constructor by base class */ 																															\
			template <class U>																																			\
			CD_PTC_S_N(DIM, N)(const CD_PTC_S_NB(DIM, N)<U>& base, const T& c_##N): CD_PTC_S_NB(DIM, N)<T>(base)														\
			{																																							\
				this->c_##N = T(c_##N);																																	\
			}																																							\
																																										\
			/* constructor by pointer */ 																																\
			template <class U>																																			\
			CD_PTC_S_N(DIM, N)(U* v, const dt_int64& n_r, const dt_int64& n_c, const dt_int64& idx, dt_int64 icol=0): CD_PTC_S_NB(DIM, N)<T>(v, n_r, n_c, idx, icol)	\
			{																																							\
				const auto ip = icol*n_r + idx;																															\
				c_##N = (n_c>N_DIM(N, DIM))?T(v[ip + N_DIM(N, DIM)*n_r]):T(1);																							\
			}																																							\
																																										\
			/* copy constructor */ 																																		\
			CD_PTC_S_N(DIM, N)(const CD_PTC_S_N(DIM, N)<T>& ptc_s)																										\
			{																																							\
				*this = ptc_s;																																			\
			}																																							\
																																										\
			/* converting constructor */ 																																\
			template <class U> 																																			\
			CD_PTC_S_N(DIM, N)(const CD_PTC_S_N(DIM, N)<U>& ptc_s)																										\
			{																																							\
				*this = ptc_s;																																			\
			}																																							\
																																										\
			/******************************** assignment operators *********************************/																	\
			/* copy assignment operator */																																\
			CD_PTC_S_N(DIM, N)<T>& operator=(const CD_PTC_S_N(DIM, N)<T>& ptc_s)																						\
			{																																							\
				if (this != &ptc_s)																																		\
				{																																						\
					CD_PTC_S_NB(DIM, N)<T>::operator=(ptc_s);																											\
																																										\
					c_##N = ptc_s.c_##N;																																\
				}																																						\
																																										\
				return *this;																																			\
			}																																							\
																																										\
			/* converting assignment operator */ 																														\
			template <class U> 																																			\
			CD_PTC_S_N(DIM, N)<T>& operator=(const CD_PTC_S_N(DIM, N)<U>& ptc_s)																						\
			{																																							\
				assign(ptc_s);																																			\
																																										\
				return *this;																																			\
			}																																							\
																																										\
			template <class U> 																																			\
			void assign(const CD_PTC_S_N(DIM, N)<U>& ptc_s)																												\
			{ 																																							\
				if ((void*)this != (void*)&ptc_s)																														\
				{																																						\
					CD_PTC_S_NB(DIM, N)<T>::assign(ptc_s);																												\
																																										\
					c_##N = T(ptc_s.c_##N);																																\
				}																																						\
			}																																							\
		}

		/***************************************************************************************/
		#define C_PTC_COEF_DIM_N(DIM, N)																																\
		template <class T>																																				\
		class CD_PTC_N(DIM, N): public CD_PTC_NB(DIM, N)<T>																												\
		{																																								\
			public:																																						\
				using value_type = T;																																	\
				using size_type = dt_int64;																																\
																																										\
				using Ptc_s = CD_PTC_S_N(DIM, N)<T>;																													\
																																										\
				mutable Vctr_cpu<T> c_##N;																																\
																																										\
				R_2d<T> c_##N##_lim;																																	\
																																										\
				/************************************* constructors ************************************/																\
				CD_PTC_N(DIM, N)(): CD_PTC_NB(DIM, N)<T>(), c_##N##_lim(){}																								\
																																										\
				template <class U>																																		\
				CD_PTC_N(DIM, N)(const pVctr_cpu_64<U>& ptc, const dt_int64& icol, const R_##DIM##d<U>& bs, dt_bool pbc_xy = false, dt_bool b_statistic = true)			\
				{																																						\
					set_ptc(ptc, icol, bs, pbc_xy, b_statistic);																										\
				}																																						\
																																										\
				/* copy constructor */																																	\
				CD_PTC_N(DIM, N)(const CD_PTC_N(DIM, N)<T>& ptc)																										\
				{																																						\
					*this = ptc;																																		\
				}																																						\
																																										\
				/* converting constructor */																															\
				template <class U>																																		\
				CD_PTC_N(DIM, N)(const CD_PTC_N(DIM, N)<U>& ptc)																										\
				{																																						\
					*this = ptc;																																		\
				}																																						\
																																										\
				/******************************** assignment operators *********************************/																\
				/* copy assignment operator */																															\
				CD_PTC_N(DIM, N)<T>& operator=(const CD_PTC_N(DIM, N)<T>& ptc)																							\
				{																																						\
					assign(ptc);																																		\
																																										\
					return *this;																																		\
				}																																						\
																																										\
				/* converting assignment operator */																													\
				template <class U> 																																		\
				CD_PTC_N(DIM, N)<T>& operator=(const CD_PTC_N(DIM, N)<U>& ptc)																							\
				{																																						\
					assign(ptc);																																		\
																																										\
					return *this;																																		\
				}																																						\
																																										\
				template <class U>																																		\
				void assign(const CD_PTC_N(DIM, N)<U>& ptc)																												\
				{																																						\
					if ((void*)this != (void*)&ptc)																														\
					{																																					\
						CD_PTC_NB(DIM, N)<T>::assign(ptc);																												\
																																										\
						c_##N = ptc.c_##N;																																\
																																										\
						c_##N##_lim = ptc.c_##N##_lim;																													\
					}																																					\
				}																																						\
																																										\
				/***************************************************************************************/																\
				virtual size_type cols() const																															\
				{																																						\
					return INC(N, DIM);																																	\
				}																																						\
																																										\
				void clear()																																			\
				{																																						\
					CD_PTC_NB(DIM, N)<T>::clear();																														\
																																										\
					c_##N.clear();																																		\
																																										\
					c_##N##_lim = 0;																																	\
				}																																						\
																																										\
				void resize(size_type new_size)																															\
				{																																						\
					new_size = max(size_type(0), new_size);																												\
																																										\
					CD_PTC_NB(DIM, N)<T>::resize(new_size);																												\
																																										\
					c_##N.resize(new_size);																																\
				}																																						\
																																										\
				void reserve(size_type new_size)																														\
				{																																						\
					new_size = max(size_type(0), new_size);																												\
																																										\
					CD_PTC_NB(DIM, N)<T>::reserve(new_size);																											\
																																										\
					c_##N.reserve(new_size);																															\
				}																																						\
																																										\
				void shrink_to_fit()																																	\
				{																																						\
					CD_PTC_NB(DIM, N)<T>::shrink_to_fit();																												\
																																										\
					c_##N.shrink_to_fit();																																\
				}																																						\
																																										\
				void push_back(const Ptc_s& ptc_s)																														\
				{																																						\
					CD_PTC_NB(DIM, N)<T>::push_back(ptc_s);																												\
																																										\
					c_##N.push_back(ptc_s.c_##N);																														\
				}																																						\
																																										\
				Ptc_s get(const size_type& ia) const																													\
				{																																						\
					return {CD_PTC_NB(DIM, N)<T>::get(ia), c_##N[ia]};																									\
				}																																						\
																																										\
				void set(const size_type& ia, const Ptc_s& ptc_s)																										\
				{																																						\
					CD_PTC_NB(DIM, N)<T>::set(ia, ptc_s);																												\
																																										\
					c_##N[ia] = ptc_s.c_##N;																															\
				}																																						\
																																										\
				template <class U>																																		\
				void set_ptc(const CD_PTC_N(DIM, N)<U>& ptc, dt_bool pbc_xy = false, dt_bool b_statistic = true)														\
				{																																						\
					clear();																																			\
					reserve(ptc.size());																																\
					mt::ptc_detail::set_ptc_pbc_xy(ptc, pbc_xy, b_statistic, *this);																					\
				}																																						\
																																										\
				template <class U>																																		\
				void set_ptc(const pVctr_cpu_64<U>& ptc, const dt_int64& icol, const R_##DIM##d<U>& bs, dt_bool pbc_xy = false, dt_bool b_statistic = true)				\
				{																																						\
					clear();																																			\
					reserve(ptc.size());																																\
					mt::ptc_detail::set_ptc_pbc_xy(ptc, icol, bs, pbc_xy, b_statistic, *this);																			\
				}																																						\
																																										\
				/* copy data to pointer */																																\
				template <class U>																																		\
				dt_int32 cpy_to_ptr(U *ptc, size_type n_ptc, dt_int32 is_0=0, dt_int32 is_e=INC(N, DIM)) const															\
				{																																						\
					if (is_0>is_e)																																		\
					{																																					\
						std::swap(is_0, is_e);																															\
					}																																					\
																																										\
					auto n_data = min(n_ptc, this->size());																												\
					auto is = CD_PTC_NB(DIM, N)<T>::cpy_to_ptr(ptc, n_ptc, is_0, is_e);																					\
																																										\
					if (fcn_chk_bound(N_DIM(N, DIM), is_0, is_e))																										\
						memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), c_##N.data(), n_data);																						\
																																										\
					return is;																																			\
				}																																						\
																																										\
				/* get statistic */																																		\
				virtual void get_statistic()																															\
				{																																						\
					if (this->empty())																																	\
					{																																					\
						return;																																			\
					}																																					\
																																										\
					CD_PTC_NB(DIM, N)<T>::get_statistic();																												\
																																										\
					fcn_minmax_element(c_##N, c_##N##_lim.x, c_##N##_lim.y);																							\
				}																																						\
																																										\
				/* sort by idx */																																		\
				FCN_SORT_BY_IDX(DIM, N);																																\
		}		
	};
}
	
/***************************************************************************************/
/******************************** position particle ************************************/
/***************************************************************************************/

/* derived class Ptc_s_2d_0 */
namespace mt
{
	template <class T>
	using Ptc_s_2d_0 = Ptc_s_xd_0<T, edim_2>;

	template <class T>
	using Ptc_R_2d = Ptc_R_xd<T, edim_2>;
}

/* template specialization 2d*/
namespace mt
{
	template <class T>
	class Ptc_s_xd_0<T, edim_2>: public Ptc_s_xd_0<T, edim_1>
	{
	public:
		T y;

		Ptc_s_xd_0(): Ptc_s_xd_0<T, edim_1>(), y(0) {}

		Ptc_s_xd_0(const T& x, const T& y): Ptc_s_xd_0<T, edim_1>(x), y(y){}

		template <class U>
		Ptc_s_xd_0(const U& x, const U& y): Ptc_s_xd_0<T, edim_1>(x), y(y){}

		/* constructor by initializer list */
		template <class U>
		Ptc_s_xd_0(const dt_init_list<U>& list): Ptc_s_xd_0<T, edim_1>(list)
		{
			auto ptr = list.begin();

			y = T(ptr[1]);
		}

		/* constructor by pointer */
		template <class U>
		Ptc_s_xd_0(U* v, const dt_int64& n_r, const dt_int64& n_c, const dt_int64& idx, dt_int64 icol=0): Ptc_s_xd_0<T, edim_1>(v, n_r, n_c, idx, icol)
		{
			const auto ip = icol*n_r + idx;

			y = (n_c>1)?T(v[ip + 1*n_r]):T(0);
		}

		/* copy constructor */
		CGPU_EXEC
		Ptc_s_xd_0(const Ptc_s_xd_0<T, edim_2>& ptc)
		{
			*this = ptc;
		}

		/* converting constructor */
		template <class U> 
		CGPU_EXEC
		Ptc_s_xd_0(const Ptc_s_xd_0<U, edim_2>& ptc)
		{
			*this = ptc;
		}

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		Ptc_s_xd_0<T, edim_2>& operator=(const Ptc_s_xd_0<T, edim_2>& ptc)
		{
			if (this != &ptc)
			{
				Ptc_s_xd_0<T, edim_1>::operator=(ptc);

				y = ptc.y;
			}
			
			return *this;
		}

		/* converting assignment operator */
		template <class U> 
		CGPU_EXEC
		Ptc_s_xd_0<T, edim_2>& operator=(const Ptc_s_xd_0<U, edim_2>& ptc)
		{
			Ptc_s_xd_0<T, edim_1>::operator=(ptc);

			y = T(ptc.y);
			
			return *this;
		}

		template <class U> 
		CGPU_EXEC
		void assign(const Ptc_s_xd_0<U, edim_2>& ptc)
		{ 
			*this = ptc;
		}

		CGPU_EXEC
		void set_pos(const R_2d<T>& r)
		{
			this->x = r.x;
			y = r.y;
		}

		CGPU_EXEC
		R_2d<T> get_pos() const
		{
			return {this->x, y};
		}
	};
}

/* template specialization 2d */
namespace mt
{
	template <class T>
	class Ptc_R_xd<T, edim_2>
	{
		public:
			using value_type = T;
			using size_type = dt_int64;

			using Ptc_s = Ptc_s_2d_0<T>;

			R_2d<T> bs;				// box size

			mutable Vctr_cpu<T> x;
			mutable Vctr_cpu<T> y;

			R_2d<T> x_lim;
			R_2d<T> y_lim;

			R_2d<T> r_mean;			// mean position
			R_2d<T> r_std;			// standard deviation
			R_2d<T> sz;				// size

			/************************************* constructors ************************************/
			Ptc_R_xd(): bs(), x_lim(), y_lim(), r_mean(), r_std(), sz() {}

			template <class U>
			Ptc_R_xd(const pVctr_cpu_64<U>& ptc, const dt_int64& icol, const R_2d<U>& bs, dt_bool pbc_xy = false, dt_bool b_statistic = true)
			{
				set_ptc(ptc, icol, bs, pbc_xy, b_statistic);
			}

			/* copy constructor */
			Ptc_R_xd(const Ptc_R_xd<T, edim_2>& ptc)
			{
				*this = ptc;
			}

			/* converting constructor */
			template <class U>
			Ptc_R_xd(const Ptc_R_xd<U, edim_2>& ptc)
			{
				*this = ptc;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			Ptc_R_xd<T, edim_2>& operator=(const Ptc_R_xd<T, edim_2>& ptc)
			{
				assign(ptc);

				return *this;
			}

			/* converting assignment operator */
			template <class U> 
			Ptc_R_xd<T, edim_2>& operator=(const Ptc_R_xd<U, edim_2>& ptc)
			{
				assign(ptc);

				return *this;
			}

			template <class U>
			void assign(const Ptc_R_xd<U, edim_2>& ptc)
			{
				if ((void*)this != (void*)&ptc)
				{
					bs = ptc.bs;

					x = ptc.x;
					y = ptc.y;

					x_lim = ptc.x_lim;
					y_lim = ptc.y_lim;

					r_mean = ptc.r_mean;
					r_std = ptc.r_std;
					sz = ptc.sz;
				}
			}

			/***************************************************************************************/
			dt_shape_st<size_type> shape() const
			{
				return {size(), cols(), 1, 1};
			}

			size_type size() const
			{
				return x.size();
			}

			dt_int32 size_32() const
			{
				return static_cast<dt_int32>(x.size());
			}	

			virtual size_type cols() const
			{
				return 2;
			}

			dt_bool empty() const
			{
				return size() == 0;
			}

			void clear()
			{
				bs = 0;

				x.clear();
				y.clear();

				x_lim = 0;
				y_lim = 0;

				r_mean = 0;
				r_std = 0;
				sz = 0;
			}

			void resize(size_type new_size)
			{
				new_size = max(size_type(0), new_size);

				x.resize(new_size);
				y.resize(new_size);
			}

			void reserve(size_type new_size)
			{
				new_size = max(size_type(0), new_size);

				x.reserve(new_size);
				y.reserve(new_size);
			}

			void shrink_to_fit()
			{
				x.shrink_to_fit();
				y.shrink_to_fit();
			}

			void push_back(const Ptc_s& ptc_s)
			{
				x.push_back(ptc_s.x);
				y.push_back(ptc_s.y);
			}

			template <class U>
			void set_bs(const R_2d<U>& bs)
			{
				this->bs = bs;
			}

			Ptc_s get(const size_type& ia) const
			{
				return {x[ia], y[ia]};
			}

			void set(const size_type& ia, const Ptc_s& ptc_s)
			{
				x[ia] = ptc_s.x;
				y[ia] = ptc_s.y;
			}

			R_2d<T> get_pos(const size_type& ia) const
			{
				return {x[ia], y[ia]};
			}

			void set_pos(const size_type& ia, const R_2d<T>& r)
			{
				x[ia] = r.x;
				y[ia] = r.y;
			}

			template <class U>
			void set_ptc(const Ptc_R_xd<U, edim_2>& ptc, dt_bool pbc_xy = false, dt_bool b_statistic = true)
			{
				mt::ptc_detail::set_ptc_pbc_xy(ptc, pbc_xy, b_statistic, *this);
			}

			template <class U>
			void set_ptc(const pVctr_cpu_64<U>& ptc, const dt_int64& icol, const R_2d<U>& bs, dt_bool pbc_xy = false, dt_bool b_statistic = true)
			{
				mt::ptc_detail::set_ptc_pbc_xy(ptc, icol, bs, pbc_xy, b_statistic, *this);
			}
		
			/* copy data to pointer */
			template <class U>
			dt_int32 cpy_to_ptr(U *ptc, size_type n_ptc, dt_int32 is_0=0, dt_int32 is_e=2) const
			{
				if (is_0>is_e)
				{
					std::swap(is_0, is_e);
				}

				auto n_data = min(n_ptc, size());
				dt_int32 is = 0;

				if (fcn_chk_bound(0, is_0, is_e))
					memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), x.data(), n_data);				// x-position

				if (fcn_chk_bound(1, is_0, is_e))
					memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), y.data(), n_data);				// y-position

				return is;
			}

			// sort by x
			void sort_by_x()
			{
				sort_by_idx(0);
			}

			// sort by y
			void sort_by_y()
			{
				sort_by_idx(1);
			}

			// sort by idx
			virtual void sort_by_idx(const dt_int32& idx)
			{
				auto first = fcn_mkzipiter_begin(this->x, this->y);
				auto last = fcn_mkzipiter_end(this->x, this->y);

				switch(idx)
				{
					case 0: 
					{ 
						thrust::sort(first, last, mt::cgpu_fctr::less_soa<0>()); 
					} 
					break;
					case 1:
					{ 
						thrust::sort(first, last, mt::cgpu_fctr::less_soa<1>()); 
					} 
					break;
				}
			}

			/***************************************************************************************/
			T norm_2_pbc_xy(const size_type& ia, const R_2d<T>& r_0) const
			{
				auto r = get_pos(ia) - r_0;

				r.x = ::fabs(r.x);
				r.y = ::fabs(r.y);

				r.x = ::fmin(r.x, ::fabs(r.x-bs.x));
				r.y = ::fmin(r.y, ::fabs(r.y-bs.y));

				return mt::norm_2(r);
			}

			T norm_2_pbc(const size_type& ia, const R_2d<T>& r_0) const
			{
				auto r = get_pos(ia) - r_0;

				r.x = fabs(r.x);
				r.y = fabs(r.y);

				r.x = ::fmin(r.x, fabs(r.x-bs.x));
				r.y = ::fmin(r.y, fabs(r.y-bs.y));

				return mt::norm_2(r);
			}

			T norm_2(const size_type& ia, const R_2d<T>& r_0) const
			{
				return mt::norm_2(get_pos(ia) - r_0);
			}

			T norm_2(const size_type& ia, const T& x, const T& y) const
			{
				return mt::norm_2(get_pos(ia) - R_2d<T>(x, y));
			}

			T norm_2(const size_type& ia_0, const size_type& ia_e) const
			{
				return mt::norm_2(get_pos(ia_0) - get_pos(ia_e));
			}

			/***************************************************************************************/
			T norm_pbc_xy(const size_type& ia, const R_2d<T>& r_0) const
			{
				return ::sqrt(this->norm_2_pbc_xy(ia, r_0));
			}				
				
			T norm_pbc(const size_type& ia, const R_2d<T>& r_0) const
			{
				return ::sqrt(this->norm_2_pbc(ia, r_0));
			}

			T norm(const size_type& ia, const R_2d<T>& r_0) const
			{
				return ::sqrt(this->norm_2(ia, r_0));
			}

			T norm(const size_type& ia, const T& x, const T& y) const
			{
				return ::sqrt(this->norm_2(ia, x, y));
			}

			T norm(const size_type& ia_0, const size_type& ia_e) const
			{
				return ::sqrt(this->norm_2(ia_0, ia_e));
			}

			/***************************************************************************************/
			virtual void get_statistic()
			{
				fcn_ptc_pos_statistic(*this);
			}

			void shift(const R_2d<T>& r_sft)
			{
				mt::fcn_ptc_pos_shift(r_sft, *this);
			}

			void recenter(const R_2d<T>& bs)
			{
				mt::fcn_ptc_pos_recenter(bs, *this);
			}

			void recenter()
			{
				mt::fcn_ptc_pos_recenter(bs, *this);
			}

			void apply_ltf(const Mx_2x2<T>& mx, const R_2d<T>& p)
			{
				fcn_ptc_pos_apply_ltf(mx, p, *this); 
			}

			void rotate(const T& theta, const R_2d<T>& p)
			{
				mt::fcn_ptc_pos_rotate(theta, p, *this);
			}
	};
}

/* fcns 2d */
namespace mt
{
	template <class T>
	void fcn_ptc_pos_statistic(Ptc_R_2d<T>& ptc)
	{
		if (ptc.empty())
		{
			return;
		}

		fcn_minmax_element(ptc.x, ptc.x_lim.x, ptc.x_lim.y);
		fcn_minmax_element(ptc.y, ptc.y_lim.x, ptc.y_lim.y);
				
		fcn_mean_std(ptc.x, ptc.r_mean.x, ptc.r_std.x);
		fcn_mean_std(ptc.y, ptc.r_mean.y, ptc.r_std.y);

		ptc.sz = R_2d<T>(ptc.x_lim.y - ptc.x_lim.x, ptc.y_lim.y - ptc.y_lim.x);
	}

	template <class T>
	void fcn_ptc_pos_shift(const R_2d<T>& r_sft, Ptc_R_2d<T>& ptc)
	{
		for(auto ia = 0; ia < ptc.size(); ia++)
		{
			ptc.x[ia] += r_sft.x;
			ptc.y[ia] += r_sft.y;
		}

			fcn_ptc_pos_statistic(ptc);
	}

	template <class T>
	void fcn_ptc_pos_recenter(const R_2d<T>& bs, Ptc_R_2d<T>& ptc)
	{
		const R_2d<T> r_sft = (bs-ptc.sz)/T(2) - R_2d<T>(ptc.x_lim.x, ptc.y_lim.x);

		fcn_ptc_pos_shift(r_sft, ptc);
	}

	template <class T>
	void fcn_ptc_pos_apply_ltf(const Mx_2x2<T>& mx, const R_2d<T>& p, Ptc_R_2d<T>& ptc)
	{
		for(dt_int64 ia = 0; ia < ptc.size(); ia++)
		{
			auto r = mx*ptc.get_pos(ia) + p;

			ptc.x[ia] = r.x;
			ptc.y[ia] = r.y;
		}

		fcn_ptc_pos_statistic(ptc);
	}
				
	template <class T>
	void fcn_ptc_pos_rotate(const T& theta, const R_2d<T>& p, Ptc_R_2d<T>& ptc)
	{
		const auto Rm = fcn_rot_mx_2d(theta);
		const auto p_sft = p - Rm*p;

		fcn_ptc_pos_apply_ltf(Rm, p_sft, ptc);
	}
}



/* derived class Ptc_s_3d_0 */
namespace mt
{
	template <class T>
	using Ptc_s_3d_0 = Ptc_s_xd_0<T, edim_3>;

	template <class T>
	using Ptc_R_3d = Ptc_R_xd<T, edim_3>;
}

/* template specialization 3d*/
namespace mt
{
	template <class T>
	class Ptc_s_xd_0<T, edim_3>: public Ptc_s_xd_0<T, edim_2>
	{
	public:
		T z;

		Ptc_s_xd_0(): Ptc_s_xd_0<T, edim_2>(), z(0) {}

		Ptc_s_xd_0(const T& x, const T& y, const T& z): Ptc_s_xd_0<T, edim_2>(x, y), z(z){}

		template <class U>
		Ptc_s_xd_0(const U& x, const U& y, const U& z): Ptc_s_xd_0<T, edim_2>(x, y), z(z){}

		/* constructor by initializer list */
		template <class U>
		Ptc_s_xd_0(const dt_init_list<U>& list): Ptc_s_xd_0<T, edim_2>(list)
		{
			auto ptr = list.begin();

			z = T(ptr[2]);
		}

		/* constructor by pointer */
		template <class U>
		Ptc_s_xd_0(U* v, const dt_int64& n_r, const dt_int64& n_c, const dt_int64& idx, dt_int64 icol=0): Ptc_s_xd_0<T, edim_2>(v, n_r, n_c, idx, icol)
		{
			const auto ip = icol*n_r + idx;

			z = (n_c>1)?T(v[ip + 2*n_r]):T(0);
		}

		/* copy constructor */
		CGPU_EXEC
		Ptc_s_xd_0(const Ptc_s_xd_0<T, edim_3>& ptc)
		{
			*this = ptc;
		}

		/* converting constructor */
		template <class U> 
		CGPU_EXEC
		Ptc_s_xd_0(const Ptc_s_xd_0<U, edim_3>& ptc)
		{
			*this = ptc;
		}

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		CGPU_EXEC
		Ptc_s_xd_0<T, edim_3>& operator=(const Ptc_s_xd_0<T, edim_3>& ptc)
		{
			if (this != &ptc)
			{
				Ptc_s_xd_0<T, edim_2>::operator=(ptc);

				z = ptc.z;
			}
			
			return *this;
		}

		/* converting assignment operator */
		template <class U> 
		CGPU_EXEC
		Ptc_s_xd_0<T, edim_3>& operator=(const Ptc_s_xd_0<U, edim_3>& ptc)
		{
			Ptc_s_xd_0<T, edim_2>::operator=(ptc);

			z = T(ptc.z);
			
			return *this;
		}

		template <class U> 
		CGPU_EXEC
		void assign(const Ptc_s_xd_0<U, edim_3>& ptc)
		{ 
			*this = ptc;
		}

		CGPU_EXEC
		void set_pos(const R_3d<T>& r)
		{
			this->x = r.x;
			this->y = r.y;
			z = r.z;
		}

		CGPU_EXEC
		R_3d<T> get_pos() const
		{
			return {this->x, this->y, z};
		}
	};
}

/* template specialization 3d */
namespace mt
{
	template <class T>
	class Ptc_R_xd<T, edim_3>
	{
		public:
			using value_type = T;
			using size_type = dt_int64;

			using Ptc_s = Ptc_s_3d_0<T>;

			R_3d<T> bs;				// box size

			mutable Vctr_cpu<T> x;
			mutable Vctr_cpu<T> y;
			mutable Vctr_cpu<T> z;

			R_2d<T> x_lim;			// x limits
			R_2d<T> y_lim;			// y limits
			R_2d<T> z_lim;			// z limits

			R_3d<T> r_mean;			// mean position
			R_3d<T> r_std;			// standard deviation
			R_3d<T> sz;				// size

			/************************************* constructors ************************************/
			Ptc_R_xd(): bs(), x_lim(), y_lim(), z_lim(), r_mean(), r_std(), sz() {}

			template <class U>
			Ptc_R_xd(const pVctr_cpu_64<U>& ptc, const dt_int64& icol, const R_3d<U>& bs, dt_bool pbc_xy = false, dt_bool b_statistic = true)
			{
				set_ptc(ptc, icol, bs, pbc_xy, b_statistic);
			}

			/* copy constructor */
			Ptc_R_xd(const Ptc_R_xd<T, edim_3>& ptc)
			{
				*this = ptc;
			}

			/* converting constructor */
			template <class U>
			Ptc_R_xd(const Ptc_R_xd<U, edim_3>& ptc)
			{
				*this = ptc;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			Ptc_R_xd<T, edim_3>& operator=(const Ptc_R_xd<T, edim_3>& ptc)
			{
				assign(ptc);

				return *this;
			}

			/* converting assignment operator */
			template <class U> 
			Ptc_R_xd<T, edim_3>& operator=(const Ptc_R_xd<U, edim_3>& ptc)
			{
				assign(ptc);

				return *this;
			}

			template <class U>
			void assign(const Ptc_R_xd<U, edim_3>& ptc)
			{
				if ((void*)this != (void*)&ptc)
				{
					bs = ptc.bs;

					x = ptc.x;
					y = ptc.y;
					z = ptc.z;

					x_lim = ptc.x_lim;
					y_lim = ptc.y_lim;
					z_lim = ptc.z_lim;

					r_mean = ptc.r_mean;
					r_std = ptc.r_std;
					sz = ptc.sz;
				}
			}

			/***************************************************************************************/
			dt_shape_st<size_type> shape() const
			{
				return {size(), cols(), 1, 1};
			}

			size_type size() const
			{
				return x.size();
			}
	
			dt_int32 size_32() const
			{
				return static_cast<dt_int32>(x.size());
			}

			virtual size_type cols() const
			{
				return 3;
			}

			dt_bool empty() const
			{
				return size() == 0;
			}

			void clear()
			{
				bs = 0;

				x.clear();
				y.clear();
				z.clear();

				x_lim = 0;
				y_lim = 0;
				z_lim = 0;

				r_mean = 0;
				r_std = 0;
				sz = 0;
			}

			void resize(size_type new_size)
			{
				new_size = max(size_type(0), new_size);

				x.resize(new_size);
				y.resize(new_size);
				z.resize(new_size);
			}

			void reserve(size_type new_size)
			{
				new_size = max(size_type(0), new_size);

				x.reserve(new_size);
				y.reserve(new_size);
				z.reserve(new_size);
			}

			void shrink_to_fit()
			{
				x.shrink_to_fit();
				y.shrink_to_fit();
				z.shrink_to_fit();
			}

			void push_back(const Ptc_s& ptc_s)
			{
				x.push_back(ptc_s.x);
				y.push_back(ptc_s.y);
				z.push_back(ptc_s.z);
			}

			template <class U>
			void set_bs(const R_3d<U>& bs)
			{
				this->bs = bs;
			}

			Ptc_s get(const size_type& ia) const
			{
				return {x[ia], y[ia], z[ia]};
			}

			void set(const size_type& ia, const Ptc_s& ptc_s)
			{
				x[ia] = ptc_s.x;
				y[ia] = ptc_s.y;
				z[ia] = ptc_s.z;
			}

			R_3d<T> get_pos(const size_type& ia) const
			{
				return {x[ia], y[ia], z[ia]};
			}

			void set_pos(const size_type& ia, const R_3d<T>& r)
			{
				x[ia] = r.x;
				y[ia] = r.y;
				z[ia] = r.z;
			}

			template <class U>
			void set_ptc(const Ptc_R_xd<U, edim_3>& ptc, dt_bool pbc_xy = false, dt_bool b_statistic = true)
			{
				clear();
				reserve(ptc.size());

				mt::ptc_detail::set_ptc_pbc_xy(ptc, pbc_xy, b_statistic, *this);
			}

			template <class U>
			void set_ptc(const pVctr_cpu_64<U>& ptc, const dt_int64& icol, const R_3d<U>& bs, dt_bool pbc_xy = false, dt_bool b_statistic = true)
			{
				clear();
				reserve(ptc.size());

				mt::ptc_detail::set_ptc_pbc_xy(ptc, icol, bs, pbc_xy, b_statistic, *this);
			}
		
			/* copy data to pointer */
			template <class U>
			dt_int32 cpy_to_ptr(U *ptc, size_type n_ptc, dt_int32 is_0=0, dt_int32 is_e=3) const
			{
				if (is_0>is_e)
				{
					std::swap(is_0, is_e);
				}

				auto n_data = min(n_ptc, size());
				dt_int32 is = 0;

				if (fcn_chk_bound(0, is_0, is_e))
					memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), x.data(), n_data);				// x-position

				if (fcn_chk_bound(1, is_0, is_e))
					memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), y.data(), n_data);				// y-position

				if (fcn_chk_bound(2, is_0, is_e))
					memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), z.data(), n_data);				// z-position

				return is;
			}

			// sort by x
			void sort_by_x()
			{
				sort_by_idx(0);
			}

			// sort by y
			void sort_by_y()
			{
				sort_by_idx(1);
			}

			// sort by z
			void sort_by_z()
			{
				sort_by_idx(2);
			}

			// sort by idx
			virtual void sort_by_idx(const dt_int32& idx)
			{
				auto first = fcn_mkzipiter_begin(this->x, this->y, this->z);
				auto last = fcn_mkzipiter_end(this->x, this->y, this->z);

				switch(idx)
				{
					case 0: 
					{ 
						thrust::sort(first, last, mt::cgpu_fctr::less_soa<0>()); 
					} 
					break;
					case 1:
					{ 
						thrust::sort(first, last, mt::cgpu_fctr::less_soa<1>()); 
					} 
					break;
					case 2:
					{ 
						thrust::sort(first, last, mt::cgpu_fctr::less_soa<2>()); 
					} 
					break;
				}
			}

			/***************************************************************************************/
			T norm_2_pbc_xy(const size_type& ia, const R_3d<T>& r_0) const
			{
				auto r = get_pos(ia) - r_0;

				r.x = fabs(r.x);
				r.y = fabs(r.y);

				r.x = ::fmin(r.x, fabs(r.x-bs.x));
				r.y = ::fmin(r.y, fabs(r.y-bs.y));

				return mt::norm_2(r);
			}				
				
			T norm_2_pbc(const size_type& ia, const R_3d<T>& r_0) const
			{
				auto r = get_pos(ia) - r_0;

				r.x = fabs(r.x);
				r.y = fabs(r.y);
				r.z = fabs(r.z);

				r.x = ::fmin(r.x, fabs(r.x-bs.x));
				r.y = ::fmin(r.y, fabs(r.y-bs.y));
				r.z = ::fmin(r.z, fabs(r.y-bs.z));

				return mt::norm_2(r);
			}

			T norm_2(const size_type& ia, const R_3d<T>& r_0) const
			{
				return mt::norm_2(get_pos(ia) - r_0);
			}

			T norm_2(const size_type& ia, const T& x, const T& y, const T& z) const
			{
				return mt::norm_2(get_pos(ia) - R_3d<T>(x, y, z));
			}

			T norm_2(const size_type& ia_0, const size_type& ia_e) const
			{
				return mt::norm_2(get_pos(ia_0) - get_pos(ia_e));
			}

			/***************************************************************************************/
			T norm_pbc_xy(const size_type& ia, const R_3d<T>& r_0) const
			{
				return ::sqrt(this->norm_2_pbc_xy(ia, r_0));
			}				
				
			T norm_pbc(const size_type& ia, const R_3d<T>& r_0) const
			{
				return ::sqrt(this->norm_2_pbc(ia, r_0));
			}

			T norm(const size_type& ia, const R_3d<T>& r_0) const
			{
				return ::sqrt(this->norm_2(ia, r_0));
			}

			T norm(const size_type& ia, const T& x, const T& y, const T& z) const
			{
				return ::sqrt(this->norm_2(ia, x, y, z));
			}

			T norm(const size_type& ia_0, const size_type& ia_e) const
			{
				return ::sqrt(this->norm_2(ia_0, ia_e));
			}

			/***************************************************************************************/
			virtual void get_statistic()
			{
				fcn_ptc_pos_statistic(*this);
			}

			void shift(const R_3d<T>& r_sft)
			{
				mt::fcn_ptc_pos_shift(r_sft, *this);
			}

			void recenter(const R_3d<T>& bs)
			{
				mt::fcn_ptc_pos_recenter(bs, *this);
			}

			void recenter()
			{
				mt::fcn_ptc_pos_recenter(bs, *this);
			}

			void recenter_xy(const R_2d<T>& bs)
			{
				fcn_ptc_pos_recenter_xy(bs, *this);
			}

			void recenter_xy()
			{
				cn_ptc_pos_recenter_xy({bs.x, bs.y}, *this);
			}

			void apply_ltf(const Mx_2x2<T>& mx, const R_3d<T>& p)
			{
				fcn_ptc_pos_apply_ltf(mx, p, *this); 
			}

			void rotate(const T& theta, const R_3d<T>& u_0, const R_3d<T>& p)
			{
				mt::fcn_ptc_pos_rotate(theta, u_0, p, *this);
			}
	};
}

/* fcns 3d */
namespace mt
{
	template <class T>
	void fcn_ptc_pos_statistic(Ptc_R_3d<T>& ptc)
	{
		if (ptc.empty())
		{
			return;
		}

		fcn_minmax_element(ptc.x, ptc.x_lim.x, ptc.x_lim.y);
		fcn_minmax_element(ptc.y, ptc.y_lim.x, ptc.y_lim.y);
		fcn_minmax_element(ptc.z, ptc.z_lim.x, ptc.z_lim.y);
				
		fcn_mean_std(ptc.x, ptc.r_mean.x, ptc.r_std.x);
		fcn_mean_std(ptc.y, ptc.r_mean.y, ptc.r_std.y);
		fcn_mean_std(ptc.z, ptc.r_mean.z, ptc.r_std.z);

		ptc.sz = R_3d<T>(ptc.x_lim.y - ptc.x_lim.x, ptc.y_lim.y - ptc.y_lim.x, ptc.z_lim.y - ptc.z_lim.x);
	}
		
	template <class T>
	void fcn_ptc_pos_shift(const R_3d<T>& r_sft, Ptc_R_3d<T>& ptc)
	{
		for(auto ia = 0; ia < ptc.size(); ia++)
		{
			ptc.x[ia] += r_sft.x;
			ptc.y[ia] += r_sft.y;
			ptc.z[ia] += r_sft.z;
		}

		fcn_ptc_pos_statistic(ptc);
	}

	template <class T>
	void fcn_ptc_pos_recenter(const R_3d<T>& bs, Ptc_R_3d<T>& ptc)
	{
		const R_3d<T> r_sft = (bs-ptc.sz)/T(2) - R_3d<T>(ptc.x_lim.x, ptc.y_lim.x, ptc.z_lim.x);

		fcn_ptc_pos_shift(r_sft, ptc);
	}

	template <class T>
	void fcn_ptc_pos_recenter_xy(const R_2d<T>& bs, Ptc_R_3d<T>& ptc)
	{
		const R_2d<T> r_sft = (bs-R_2d<T>(ptc.sz.x, ptc.sz.y))/T(2) - R_2d<T>(ptc.x_lim.x, ptc.y_lim.x);

		fcn_ptc_pos_shift({r_sft.x, r_sft.y, T(0)}, ptc);
	}

	template <class T>
	void fcn_ptc_pos_apply_ltf(const Mx_3x3<T>& mx, const R_3d<T>& p, Ptc_R_3d<T>& ptc)
	{
		for(dt_int64 ia = 0; ia < ptc.size(); ia++)
		{
			auto r = mx*ptc.get_pos(ia) + p;

			ptc.x[ia] = r.x;
			ptc.y[ia] = r.y;
		}

		fcn_ptc_pos_statistic(ptc);
	}
				
	template <class T>
	void fcn_ptc_pos_rotate(const T& theta, const R_3d<T>& u_0, const R_3d<T>& p, Ptc_R_3d<T>& ptc)
	{
		const auto Rm = fcn_rot_mx_3d(theta, u_0);
		const auto p_sft = p - Rm*p;

		fcn_ptc_pos_apply_ltf(Rm, p_sft, ptc);
	}
}


/********************************* atomic particles ************************************/
namespace mt
{
	/******************************* forward declarations **********************************/
	template <class T>
	class Ptc_Atom;

	template <class T>
	void remove_ptc_out_z_range(const T& z_0, const T& z_e, Ptc_Atom<T>& atoms);
}

namespace mt
{
	template <class T>
	class Ptc_s_Atom: public Ptc_s_3d_0<T>
	{
	public:
		dt_int32 Z;
		dt_float32 sigma;
		dt_float32 occ;
		dt_int32 tag;
		dt_int32 charge;

		Ptc_s_Atom():Z(0), Ptc_s_3d_0<T>(), sigma(0), occ(0), tag(0), charge(0) {};

		Ptc_s_Atom(const dt_int32& Z, const T& x, const T& y, const T& z, 
		dt_float32 sigma=c_dflt_rms3d, dt_float32 occ=c_dflt_occ, dt_int32 tag=c_dflt_tag, dt_int32 charge=c_dflt_charge):
		Z(Z), Ptc_s_3d_0<T>(x, y, z), sigma(sigma), occ(occ), tag(tag), charge(charge) {};
			
		Ptc_s_Atom(const dt_int32& Z, const R_3d<T>& r, 
		dt_float32 sigma=c_dflt_rms3d, dt_float32 occ=c_dflt_occ, dt_int32 tag=c_dflt_tag, dt_int32 charge=c_dflt_charge):
		Z(Z), Ptc_s_3d_0<T>(r.x, r.y, r.z), sigma(sigma), occ(occ), tag(tag), charge(charge) {};

		/* constructor by pointer */
		template <class U>
		Ptc_s_Atom(U* v, const dt_int64& n_r, const dt_int64& n_c, const dt_int64& idx, dt_int64 icol=0)
		{
			const auto ip = icol*n_r + idx;

			Z = dt_int32(v[ip + 0*n_r]);								// atomic number
			this->x = (n_c>1)?T(v[ip + 1*n_r]):T(0);					// x-position
			this->y = (n_c>2)?T(v[ip + 2*n_r]):T(0);					// y-position
			this->z = (n_c>3)?T(v[ip + 3*n_r]):T(0);					// z-position

			sigma = (n_c>4)?dt_float32(v[ip + 4*n_r]):c_dflt_rms3d;		// standard deviation
			occ = (n_c>5)?dt_float32(v[ip + 5*n_r]):c_dflt_occ;			// occupancy
			tag = (n_c>6)?dt_int32(v[ip + 6*n_r]):c_dflt_tag;		// tag
			charge = (n_c>7)?dt_int32(v[ip + 7*n_r]):c_dflt_charge;		// charge
		}

		/* copy constructor */
		Ptc_s_Atom(const Ptc_s_Atom<T>& ptc_s)
		{
			*this = ptc_s;
		}

		/* converting constructor */
		template <class U> 
		Ptc_s_Atom(const Ptc_s_Atom<U>& ptc_s)
		{
			*this = ptc_s;
		}

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		Ptc_s_Atom<T>& operator=(const Ptc_s_Atom<T>& ptc_s)
		{
			if (this != &ptc_s)
			{
				Z = ptc_s.Z;
				this->x = ptc_s.x;
				this->y = ptc_s.y;
				this->z = ptc_s.z;
				sigma = ptc_s.sigma;
				occ = ptc_s.occ;
				tag = ptc_s.tag;
				charge = ptc_s.charge;
			}
			
			return *this;
		}

		/* converting assignment operator */
		template <class U> 
		Ptc_s_Atom<T>& operator=(const Ptc_s_Atom<U>& ptc_s)
		{
			assign(ptc_s);
			
			return *this;
		}

		template <class U> 
		void assign(const Ptc_s_Atom<U>& ptc_s)
		{ 
			if ((void*)this != (void*)&ptc_s)
			{
				Z = ptc_s.Z;
				this->x = T(ptc_s.x);
				this->y = T(ptc_s.y);
				this->z = T(ptc_s.z);
				sigma = ptc_s.sigma;
				occ = ptc_s.occ;
				tag = ptc_s.tag;
				charge = ptc_s.charge;
			}
		}
	};
}

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
			mutable Vctr_cpu<dt_int32> tag;		// tag
			mutable Vctr_cpu<dt_int32> charge;		// charge

			dt_int32 cols_used;						// number of used columns

			R_2d<dt_int32> Z_lim;
			R_2d<dt_float32> sigma_lim;
			R_2d<dt_float32> occ_lim;
			R_2d<dt_int32> tag_lim;
			R_2d<dt_int32> charge_lim;

			/************************************* constructors ************************************/
			Ptc_Atom(): Ptc_R_3d<T>(), cols_used(4), Z_lim(), sigma_lim(), occ_lim(), tag_lim(), charge_lim() {}

			template <class U>
			Ptc_Atom(const pVctr_cpu_64<U>& ptc, const R_3d<U>& bs, dt_bool pbc_xy = false, dt_bool b_statistic = true)
			{
				set_ptc(ptc, bs, pbc_xy, b_statistic);
			}

			/* copy constructor */
			Ptc_Atom(const Ptc_Atom<T>& ptc)
			{
				*this = ptc;
			}

			/* converting constructor */
			template <class U>
			Ptc_Atom(const Ptc_Atom<U>& ptc)
			{
				*this = ptc;
			}

			/******************************** assignment operators *********************************/
			/* assignment operator */
			Ptc_Atom<T>& operator=(const Ptc_Atom<T>& ptc)
			{
				assign(ptc);

				return *this;
			}

			/* converting assignment operator */
			template <class U> 
			Ptc_Atom<T>& operator=(const Ptc_Atom<U>& ptc)
			{
				assign(ptc);

				return *this;
			}

			template <class U>
			void assign(const Ptc_Atom<U>& ptc)
			{
				if ((void*)this != (void*)&ptc)
				{
					cols_used = ptc.cols_used;

					Z = ptc.Z;

					Ptc_R_3d<T>::assign(ptc);

					sigma = ptc.sigma;
					occ = ptc.occ;
					tag = ptc.tag;
					charge = ptc.charge;

					Z_lim = ptc.Z_lim;
					sigma_lim = ptc.sigma_lim;
					occ_lim = ptc.occ_lim;
					tag_lim = ptc.tag_lim;
					charge_lim = ptc.charge_lim;
				}
			}

			/***************************************************************************************/
			virtual size_type cols() const
			{
				return 8;	
			}		

			void clear()
			{
				cols_used = 4;

				Z.clear();

				Ptc_R_3d<T>::clear();

				sigma.clear();
				occ.clear();
				tag.clear();
				charge.clear();

				Z_lim = 0;
				sigma_lim = 0;
				occ_lim = 0;
				tag_lim = 0;
				charge_lim = 0;
			}

			void resize(size_type new_size)
			{
				new_size = max(size_type(0), new_size);

				Z.resize(new_size);

				Ptc_R_3d<T>::resize(new_size);

				if (cols_used>4)
					sigma.resize(new_size);

				if (cols_used>5)
					occ.resize(new_size);

				if (cols_used>6)
					tag.resize(new_size);

				if (cols_used>7)
					charge.resize(new_size);
			}

			void reserve(size_type new_size)
			{
				new_size = max(size_type(0), new_size);

				Z.reserve(new_size);

				Ptc_R_3d<T>::reserve(new_size);

				if (cols_used>4)
					sigma.reserve(new_size);

				if (cols_used>5)
					occ.reserve(new_size);

				if (cols_used>6)
					tag.reserve(new_size);

				if (cols_used>7)
					charge.reserve(new_size);
			}

			void shrink_to_fit()
			{
				Z.shrink_to_fit();

				Ptc_R_3d<T>::shrink_to_fit();		

				sigma.shrink_to_fit();
				occ.shrink_to_fit();
				tag.shrink_to_fit();
				charge.shrink_to_fit();
			}

			void push_back(const Ptc_s& ptc_s)
			{
				Z.push_back(ptc_s.Z);							// atomic number

				Ptc_R_3d<T>::push_back(ptc_s);					// xyz

				if (cols_used>4)
					sigma.push_back(ptc_s.sigma);				// standard deviation

				if (cols_used>5)
					occ.push_back(ptc_s.occ);					// occupancy

				if (cols_used>6)
					tag.push_back(abs(ptc_s.tag));				// tag

				if (cols_used>7)
					charge.push_back(ptc_s.charge);				// charge
			}

			dt_float32 get_sigma(const size_type& ia) const
			{
				return (cols_used>4)?sigma[ia]:c_dflt_rms3d;		// standard deviation
			}

			dt_float32 get_occ(const size_type& ia) const
			{
				return (cols_used>5)?occ[ia]:c_dflt_occ;			// occupancy
			}

			dt_int32 get_tag(const size_type& ia) const
			{
				return (cols_used>6)?tag[ia]:c_dflt_tag;			// tag
			}

			dt_int32 get_charge(const size_type& ia) const
			{
				return (cols_used>7)?charge[ia]:c_dflt_charge;		// charge
			}

			Ptc_s get(const size_type& ia) const
			{
				return {Z[ia], this->x[ia], this->y[ia], this->z[ia], get_sigma(ia), get_occ(ia), get_tag(ia), get_charge(ia)};
			}

			void set(const size_type& ia, const Ptc_s& ptc_s)
			{
				Z[ia] = ptc_s.Z;					// atomic number

				Ptc_R_3d<T>::set(ia, ptc_s);		// xyz

				if (cols_used>4)					// standard deviation
					sigma[ia] = ptc_s.sigma;

				if (cols_used>5)					// occupancy
					occ[ia] = ptc_s.occ;

				if (cols_used>6)					// tag
					tag[ia] = ptc_s.tag;

				if (cols_used>7)					// charge
					charge[ia] = ptc_s.charge;
			}

			template <class U>
			void set_ptc(const Ptc_Atom<U>& ptc, dt_bool pbc_xy = false, dt_bool b_statistic = true)
			{
				clear();
				cols_used = ptc.cols_used;
				reserve(ptc.size());

				mt::ptc_detail::set_ptc_pbc_xy(ptc, pbc_xy, b_statistic, *this);
			}

			template <class U>
			void set_ptc(const pVctr_cpu_64<U>& ptc, const R_3d<U>& bs, dt_bool pbc_xy = false, dt_bool b_statistic = true)
			{
				clear();
				cols_used = ptc.s1_32();
				reserve(ptc.size());

				mt::ptc_detail::set_ptc_pbc_xy(ptc, 0, bs, pbc_xy, b_statistic, *this);
			}

			/* copy data to pointer */
			template <class U>
			dt_int32 cpy_to_ptr(U *ptc, size_type n_ptc, dt_int32 is_0=0, dt_int32 is_e=8)
			{
				is_e = min(is_e, cols_used);

				if (is_0>is_e)
				{
					std::swap(is_0, is_e);
				}

				auto n_data = min(n_ptc, this->size());
				dt_int32 is = 0;

				if (fcn_chk_bound(0, is_0, is_e))
					memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), Z.data(), n_data);				// atomic number

				if (fcn_chk_bound(1, is_0, is_e))
					memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), this->x.data(), n_data);		// x-position

				if (fcn_chk_bound(2, is_0, is_e))
					memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), this->y.data(), n_data);		// y-position

				if (fcn_chk_bound(3, is_0, is_e))
					memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), this->z.data(), n_data);		// z-position

				if (fcn_chk_bound(4, is_0, is_e))
					memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), sigma.data(), n_data);			// standard deviation

				if (fcn_chk_bound(5, is_0, is_e))
					memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), occ.data(), n_data);			// occupancy

				if (fcn_chk_bound(6, is_0, is_e))
					memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), tag.data(), n_data);		// tag

				if (fcn_chk_bound(7, is_0, is_e))
					memcpy_cpu_cpu(&(ptc[(is++)*n_ptc]), charge.data(), n_data);		// charge

				return is;
			}

			// sort by x
			void sort_by_x()
			{
				sort_by_idx(1);
			}

			// sort by y
			void sort_by_y()
			{
				sort_by_idx(2);
			}

			// sort by z
			void sort_by_z()
			{
				sort_by_idx(3);
			}

			// sort by idx
			void sort_by_idx(const dt_int32 idx = 3)
			{
				if (cols_used==4)
				{
					auto first = fcn_mkzipiter_begin(Z, this->x, this->y, this->z);
					auto last = fcn_mkzipiter_end(Z, this->x, this->y, this->z);

					switch(idx)	
					{
						SWITCH_1_3D(CASE_SORT)
					}
				}
				else if (cols_used==5)
				{
					auto first = fcn_mkzipiter_begin(Z, this->x, this->y, this->z, sigma);
					auto last = fcn_mkzipiter_end(Z, this->x, this->y, this->z, sigma);

					switch(idx)	
					{
						SWITCH_2_3D(CASE_SORT)
					}
				}
				else if (cols_used==6)
				{
					auto first = fcn_mkzipiter_begin(Z, this->x, this->y, this->z, sigma, occ);
					auto last = fcn_mkzipiter_end(Z, this->x, this->y, this->z, sigma, occ);

					switch(idx)	
					{
						SWITCH_3_3D(CASE_SORT)
					}
				}
				else if (cols_used==7)
				{
					auto first = fcn_mkzipiter_begin(Z, this->x, this->y, this->z, sigma, occ, tag);
					auto last = fcn_mkzipiter_end(Z, this->x, this->y, this->z, sigma, occ, tag);

					switch(idx)	
					{
						SWITCH_4_3D(CASE_SORT)
					}
				}
				else if (cols_used==8)
				{
					auto first = fcn_mkzipiter_begin(Z, this->x, this->y, this->z, sigma, occ, tag, charge);
					auto last = fcn_mkzipiter_end(Z, this->x, this->y, this->z, sigma, occ, tag, charge);

					switch(idx)	
					{
						SWITCH_5_3D(CASE_SORT)
					}
				}
			}


			/***************************************************************************************/
			virtual void get_statistic()
			{
				if (this->empty())
				{
					return;
				}

				mt::fcn_minmax_element(Z, Z_lim.x, Z_lim.y);

				Ptc_R_3d<T>::get_statistic();	

				sigma_lim.x = c_dflt_rms3d;
				sigma_lim.y = c_dflt_rms3d;
				if (cols_used>4)
					mt::fcn_minmax_element(sigma, sigma_lim.x, sigma_lim.y);

				occ_lim.x = c_dflt_occ;
				occ_lim.y = c_dflt_occ;
				if (cols_used>5)
					mt::fcn_minmax_element(occ, occ_lim.x, occ_lim.y);

				tag_lim.x = c_dflt_tag;
				tag_lim.y = c_dflt_tag;
				if (cols_used>6)
					mt::fcn_minmax_element(tag, tag_lim.x, tag_lim.y);

				charge_lim.x = c_dflt_charge;
				charge_lim.y = c_dflt_charge;
				if (cols_used>7)
					mt::fcn_minmax_element(tag, charge_lim.x, charge_lim.y);

				this->bs.z = ::fmax(this->sz.z, this->bs.z);
			}

			// max z value within a tag
			void minmax_z_by_region(const T& tag_v, T& z_min, T& z_max)
			{
				z_min = 1e20;
				z_max = -1e20;
				for(auto iz = 0; iz < this->size(); iz++)
				{
					if (tag[iz]==tag_v)
					{
						z_max = ::fmax(z_max, this->z[iz]);
						z_min = ::fmin(z_min, this->z[iz]);
					}
				} 
			}

			void remove_ptc_out_z_range(const T& z_0, const T& z_e)
			{
				mt::remove_ptc_out_z_range(z_0, z_e, *this);
			}
	};
}

namespace mt
{
	template <class T>
	void remove_ptc_out_z_range(const T& z_0, const T& z_e, Ptc_Atom<T>& ptc)
	{
		dt_int32 ia_z = 0;
		for(dt_int64 ia = 0; ia < ptc.size(); ia++)
		{
			if (fcn_chk_bound(ptc.z[ia], z_0, z_e))
			{
				ptc.Z[ia_z] = ptc.Z[ia];
				ptc.x[ia_z] = ptc.x[ia];
				ptc.y[ia_z] = ptc.y[ia];
				ptc.z[ia_z] = ptc.z[ia];

				if (ptc.cols_used>4)
					ptc.sigma[ia_z] = ptc.sigma[ia];

				if (ptc.cols_used>5)
					ptc.occ[ia_z] = ptc.occ[ia];

				if (ptc.cols_used>6)
					ptc.tag[ia_z] = ptc.tag[ia];

				if (ptc.cols_used>7)
					ptc.charge[ia_z] = ptc.charge[ia];

				ia_z++;
			}
		}

		ptc.resize(ia_z);
		ptc.shrink_to_fit();
	}
 
}

/* forward declarations */
namespace mt
{
	/* 1d */
	template <class T>
	void fcn_ptc_pos_statistic(Ptc_R_1d<T>& ptc);

	template <class T>
	void fcn_ptc_pos_shift(const R_1d<T>& r_sft, Ptc_R_1d<T>& ptc);

	template <class T>
	void fcn_ptc_pos_recenter(const R_1d<T>& bs, Ptc_R_1d<T>& ptc);

	template <class T>
	void fcn_ptc_pos_apply_ltf(const T& mx, const R_1d<T>& p, Ptc_R_1d<T>& ptc);

	/* 2d */
	template <class T>
	void fcn_ptc_pos_statistic(Ptc_R_2d<T>& ptc);

	template <class T>
	void fcn_ptc_pos_shift(const R_2d<T>& r_sft, Ptc_R_2d<T>& ptc);

	template <class T>
	void fcn_ptc_pos_recenter(const R_2d<T>& bs, Ptc_R_2d<T>& ptc);

	template <class T>
	void fcn_ptc_pos_apply_ltf(const Mx_2x2<T>& mx, const R_2d<T>& p, Ptc_R_2d<T>& ptc);
		
	template <class T>
	void fcn_ptc_pos_rotate(const T& theta, const R_2d<T>& p, Ptc_R_2d<T>& ptc);

	/* 3d */
	template <class T>
	void fcn_ptc_pos_statistic(Ptc_R_3d<T>& ptc);

	template <class T>
	void fcn_ptc_pos_shift(const R_3d<T>& bs, Ptc_R_3d<T>& ptc);

	template <class T>
	void fcn_ptc_pos_recenter(const R_3d<T>& bs, Ptc_R_3d<T>& ptc);

	template <class T>
	void fcn_ptc_pos_recenter_xy(const R_2d<T>& bs, Ptc_R_3d<T>& ptc);

	template <class T>
	void fcn_ptc_pos_apply_ltf(const Mx_3x3<T>& mx, const R_3d<T>& p, Ptc_R_3d<T>& ptc);
		
	template <class T>
	void fcn_ptc_pos_rotate(const T& theta, const R_3d<T>& u_0, const R_3d<T>& p, Ptc_R_3d<T>& ptc);
}

/* derived class 2d */
namespace mt
{
	/***************************************************************************************/
	/************************************ pos 2d/c1 ****************************************/
	/***************************************************************************************/
	C_PTC_S_COEF_DIM_N(2, 1);	// Ptc_s_2d_1
	C_PTC_COEF_DIM_N(2, 1);		// Ptc_2d_1

	/***************************************************************************************/
	/*********************************** pos 2d/c1/c2 **************************************/
	/***************************************************************************************/
	C_PTC_S_COEF_DIM_N(2, 2);	// Ptc_s_2d_2
	C_PTC_COEF_DIM_N(2, 2);		// Ptc_2d_2

	/***************************************************************************************/
	/*********************************** pos 2d/c1/c2/c3 ***********************************/
	/***************************************************************************************/
	C_PTC_S_COEF_DIM_N(2, 3);	// Ptc_s_2d_3
	C_PTC_COEF_DIM_N(2, 3);		// Ptc_2d_3

	/***************************************************************************************/
	/********************************* pos 2d/c1/c2/c3/c4 **********************************/
	/***************************************************************************************/
	C_PTC_S_COEF_DIM_N(2, 4);	// Ptc_s_2d_4
	C_PTC_COEF_DIM_N(2, 4);		// Ptc_2d_4

	/***************************************************************************************/
	/******************************** pos 2d/c1/c2/c3/c4/c5 ********************************/
	/***************************************************************************************/
	C_PTC_S_COEF_DIM_N(2, 5);	// Ptc_s_2d_5
	C_PTC_COEF_DIM_N(2, 5);		// Ptc_2d_5

}

/* derived class 3d */
namespace mt
{
	/***************************************************************************************/
	/************************************** pos 3d/c1 **************************************/
	/***************************************************************************************/
	C_PTC_S_COEF_DIM_N(3, 1);	// Ptc_s_3d_1
	C_PTC_COEF_DIM_N(3, 1);		// Ptc_3d_1

	/***************************************************************************************/
	/************************************* pos 3d/c1/c2 ************************************/
	/***************************************************************************************/
	C_PTC_S_COEF_DIM_N(3, 2);	// Ptc_s_3d_2
	C_PTC_COEF_DIM_N(3, 2);		// Ptc_3d_2

	/***************************************************************************************/
	/*********************************** pos 3d/c1/c2/c3 ***********************************/
	/***************************************************************************************/
	C_PTC_S_COEF_DIM_N(3, 3);	// Ptc_s_3d_3
	C_PTC_COEF_DIM_N(3, 3);		// Ptc_3d_3

	/***************************************************************************************/
	/********************************** pos 3d/c1/c2/c3/c4 *********************************/
	/***************************************************************************************/
	C_PTC_S_COEF_DIM_N(3, 4);	// Ptc_s_3d_4
	C_PTC_COEF_DIM_N(3, 4);		// Ptc_3d_4

	/***************************************************************************************/
	/******************************** pos 3d/c1/c2/c3/c4/c5 ********************************/
	/***************************************************************************************/
	C_PTC_S_COEF_DIM_N(3, 5);	// Ptc_s_3d_5
	C_PTC_COEF_DIM_N(3, 5);		// Ptc_3d_5
}
