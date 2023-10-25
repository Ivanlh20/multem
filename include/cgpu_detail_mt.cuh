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

#ifndef CGPU_DETAIL_MT_H
	#define CGPU_DETAIL_MT_H

	#ifdef _MSC_VER
		#pragma once
	#endif

	#include "const_enum_mt.cuh"
	#include "type_traits_mt.cuh"
	#include "lens.cuh"
	#include "energy_loss.cuh"
	#include "cgpu_detail.cuh"

	/* pointer to atomic functions */
	namespace mt
	{
		template <class T>
		using pFcn_clnl_3_1 = T (*)(const T&, Ctpr<T>, Ctpr<T>);

		template <class T>
		using pFcn_clnl_3_2 = void (*)(const T&, Ctpr<T>, Ctpr<T>, T&, T&);

		/***************************************************************************************/
		template <class T>
		using pFcn_clnl_4_1 = T (*)(const T&, const T&, Ctpr<T>, Ctpr<T>);

		template <class T>
		using pFcn_clnl_4_2 = void (*)(const T&, const T&, Ctpr<T>, Ctpr<T>, T&, T&);

		/***************************************************************************************/
		template <class T>
		using pFcn_clnl_6_1 = T (*)(const T&, Ctpr<T>, Ctpr<T>, const dt_int32&, Ctpr<T>, Ctpr<T>);

		template <class T>
		using pFcn_clnl_6_2 = void (*)(const T&, Ctpr<T>, Ctpr<T>, const dt_int32&, Ctpr<T>, Ctpr<T>, T&, T&);
		
		/***************************************************************************************/
		template <class T>
		using pFcn_clnl_8_1 = T (*)(const T&, const T&, const T&, Ctpr<T>, Ctpr<T>, const dt_int32&, Ctpr<T>, Ctpr<T>);

		template <class T>
		using pFcn_clnl_8_2 = void (*)(const T&, const T&, const T&, Ctpr<T>, Ctpr<T>, const dt_int32&, Ctpr<T>, Ctpr<T>, T&, T&);
	}

	/* atomic functions */
	namespace mt
	{
		namespace cgpu_detail_mt
		{
			template <class T>
			CGPU_EXEC_INL 
			T fcn_spt_exp_v(const T& x, const dt_int32& k_0, const dt_int32& k_e, Ctpr<T> cl, Ctpr<T> cnl)
			{
				T y = 0;

				for(auto ik = k_0; ik <k_e; ik++)
				{
					y += cl[ik]*exp(-cnl[ik]*x);
				}

				return y;
			}

			template <class T, dt_bool init = true>
			CGPU_EXEC_INL 
			void fcn_spt_exp_vd(const T& x, const dt_int32& k_0, const dt_int32& k_e, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				if (init)
				{
					y = dy = 0;
				}

				for(auto ik = k_0; ik <k_e; ik++)
				{
					const auto yt = cl[ik]*exp(-cnl[ik]*x);
					y += cl[ik]*exp(-cnl[ik]*x);
					dy += -cnl[ik]*yt;
				}
			}

			template <class T>
			CGPU_EXEC_INL 
			T fcn_spt_gauss_v(const T& x, const dt_int32& k_0, const dt_int32& k_e, Ctpr<T> cl, Ctpr<T> cnl)
			{
				T y = 0;

				const T x2 = x*x;

				for(auto ik = k_0; ik <k_e; ik++)
				{
					y += cl[ik]*exp(-cnl[ik]*x2);
				}

				return y;
			}

			template <class T, dt_bool init = true>
			CGPU_EXEC_INL 
			void fcn_spt_gauss_vd(const T& x, const dt_int32& k_0, const dt_int32& k_e, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				if (init)
				{
					y = dy = 0;
				}

				const T x2 = x*x;

				for(auto ik = k_0; ik <k_e; ik++)
				{
					const auto yt = cl[ik]*exp(-cnl[ik]*x2);
					y += yt;
					dy += -2*cnl[ik]*x*yt;
				}
			}

			template <class T>
			CGPU_EXEC_INL 
			T fcn_spt_lorentzian_v(const T& x, const dt_int32& k_0, const dt_int32& k_e, Ctpr<T> cl, Ctpr<T> cnl)
			{
				T y = 0;

				const T x2 = x*x;

				for(auto ik = k_0; ik <k_e; ik++)
				{
					y += cl[ik]/(cnl[ik] + x2);
				}

				return y;
			}

			template <class T, dt_bool init = true>
			CGPU_EXEC_INL 
			void fcn_spt_lorentzian_vd(const T& x, const dt_int32& k_0, const dt_int32& k_e, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{				
				if (init)
				{
					y = dy = 0;
				}

				const T x2 = x*x;

				for(auto ik = k_0; ik <k_e; ik++)
				{
					const auto t = T(1)/(cnl[ik] + x2);
					const auto yt = cl[ik]*t;
					y += cl[ik]*t;
					dy += T(-2)*x*yt*t;
				}
			}

			template <class T>
			CGPU_EXEC_INL 
			T fcn_spt_yukawa_v(const T& x, const dt_int32& k_0, const dt_int32& k_e, Ctpr<T> cl, Ctpr<T> cnl)
			{
				T y = 0;

				const T ix = T(1)/x;

				for(auto ik = k_0; ik <k_e; ik++)
				{
					y += cl[ik]*exp(-cnl[ik]*x)*ix;
				}

				return y;
			}

			template <class T, dt_bool init = true>
			CGPU_EXEC_INL 
			void fcn_spt_yukawa_vd(const T& x, const dt_int32& k_0, const dt_int32& k_e, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				if (init)
				{
					y = dy = 0;
				}

				const T ix = T(1)/x;

				for(auto ik = k_0; ik <k_e; ik++)
				{
					const auto yt = cl[ik]*exp(-cnl[ik]*x)*ix;
					y += yt;
					dy += -(cnl[ik]+ ix)*yt;
				}
			}

			template <class T>
			CGPU_EXEC_INL 
			T fcn_spt_pr_gauss_feg_v(const T& x, const dt_int32& k_0, const dt_int32& k_e, Ctpr<T> cl, Ctpr<T> cnl)
			{
				T y = 0;

				const T x2 = x*x;

				for(auto ik = k_0; ik <k_e; ik++)
				{
					y += (T(2)*x2 - T(3)/cnl[ik])*cl[ik]*exp(-cnl[ik]*x2);
				}

				return y;
			}

			template <class T, dt_bool init = true>
			CGPU_EXEC_INL 
			void fcn_spt_pr_gauss_feg_vd(const T& x, const dt_int32& k_0, const dt_int32& k_e, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				if (init)
				{
					y = dy = 0;
				}

				const T x2 = x*x;

				for(auto ik = k_0; ik <k_e; ik++)
				{
					const auto yt = cl[ik]*exp(-cnl[ik]*x2);
					y += (T(2)*x2 - T(3)/cnl[ik])*yt;
					dy += -T(2)*x*(T(2)*cnl[ik]*x2 - T(5))*yt;
				}
			}

			/***************************************************************************************/
			/***************************************************************************************/
			template <class T>
			CGPU_EXEC_INL 
			T fcn_int_fx_x0_xe(const T& x, const T& x_0, const T& x_e, Ctpr<T> cl, Ctpr<T> cnl, const dt_int32& n_q, Ctpr<T> qx, Ctpr<T> qw, pFcn_clnl_3_1<T> fcn)
			{
				const T a = 0.5*(x_e-x_0);
				const T b = 0.5*(x_e+x_0);
				T s = 0;

				for(auto ik = 0; ik <n_q; ik++)
				{
					const auto xi = a*qx[ik] + b;
					const auto yi = fcn(xi, cl, cnl);
					s += a*qw[ik]*yi;
				}

				return s;
			}

			template <class T>
			CGPU_EXEC_INL 
			void fcn_int_fx_x0_xe(const T& x, const T& x_0, const T& x_e, Ctpr<T> cl, Ctpr<T> cnl, const dt_int32& n_q, Ctpr<T> qx, Ctpr<T> qw, T& s1, T& s2, pFcn_clnl_3_2<T> fcn)
			{
				const T a = 0.5*(x_e-x_0);
				const T b = 0.5*(x_e+x_0);
				s1 = s2 = 0;

				for(auto ik = 0; ik <n_q; ik++)
				{
					const auto xi = a*qx[ik] + b;
					T y1i, y2i;
					fcn(xi, cl, cnl, y1i, y2i);
					s1 += a*qw[ik]*y1i;
					s2 += a*qw[ik]*y2i;
				}
			}

			template <class T>
			CGPU_EXEC_INL 
			T fcn_int_fx_x0_pinfty(const T& x, const T& x_0, Ctpr<T> cl, Ctpr<T> cnl, const dt_int32& n_q, Ctpr<T> qx, Ctpr<T> qw, pFcn_clnl_3_1<T> fcn)
			{
				T s = 0;

				for(auto ik = 0; ik <n_q; ik++)
				{
					const auto xi = qx[ik] + x_0;
					const auto yi = fcn(xi, cl, cnl);
					s += qw[ik]*yi;
				}

				return s;
			}

			template <class T>
			CGPU_EXEC_INL 
			void fcn_int_fx_x0_pinfty(const T& x, const T& x_0, Ctpr<T> cl, Ctpr<T> cnl, const dt_int32& n_q, Ctpr<T> qx, Ctpr<T> qw, T& s1, T& s2, pFcn_clnl_3_2<T> fcn)
			{
				s1 = s2 = 0;

				for(auto ik = 0; ik <n_q; ik++)
				{
					const auto xi = qx[ik] + x_0;
					T y1i, y2i;
					fcn(xi, cl, cnl, y1i, y2i);
					s1 += qw[ik]*y1i;
					s2 += qw[ik]*y2i;
				}
			}

			template <class T>
			CGPU_EXEC_INL 
			T fcn_int_vz_z0_ze(const T& r, const T& z_0, const T& z_e, Ctpr<T> cl, Ctpr<T> cnl, const dt_int32& n_q, Ctpr<T> qx, Ctpr<T> qw, pFcn_clnl_3_1<T> fcn)
			{
				const dt_bool split = (z_0<0) && (0<z_e);
				T a = (split)?-0.5*z_0:0.5*(z_e-z_0);
				T b = (split)?0.5*z_0:0.5*(z_e+z_0);
				const T r2 = r*r;
				T s = 0;

				for(auto ik = 0; ik <n_q; ik++)
				{
					const auto zi = a*qx[ik] + b;
					const auto ri = ::sqrt(zi*zi + r2);
					const auto fi = fcn(ri, cl, cnl);
					s += a*qw[ik]*fi;
				}

				if (split)
				{
					a = b = 0.5*z_e;
					for(auto ik = 0; ik <n_q; ik++)
					{
						const auto zi = a*qx[ik] + b;
						const auto ri = ::sqrt(zi*zi + r2);
						const auto fi = fcn(ri, cl, cnl);
						s += a*qw[ik]*fi;
					}
				}

				return s;
			}

			template <class T>
			CGPU_EXEC_INL 
			void fcn_int_vz_dvz_z0_ze( const T& r, const T& z_0, const T& z_e, Ctpr<T> cl, Ctpr<T> cnl, const dt_int32& n_q, Ctpr<T> qx, Ctpr<T> qw, T& s, T& ds, pFcn_clnl_3_2<T> fcn)
			{
				dt_bool split = (z_0<0) && (0<z_e);
				T a = (split)?-0.5*z_0:0.5*(z_e-z_0);
				T b = (split)?0.5*z_0:0.5*(z_e+z_0);
				const T r2 = r*r;
				s = ds = 0;

				for(auto ik = 0; ik <n_q; ik++)
				{
					const auto zi = a*qx[ik] + b;
					const auto ri = ::sqrt(zi*zi + r2);
					T fi, dfi;
					fcn(ri, cl, cnl, fi, dfi);
					s += a*qw[ik]*fi;
					ds += r*a*qw[ik]*dfi/ri;
				}

				if (split)
				{
					a = b = 0.5*z_e;
					for(auto ik = 0; ik <n_q; ik++)
					{
						const auto zi = a*qx[ik] + b;
						const auto ri = ::sqrt(zi*zi + r2);
						T fi, dfi;
						fcn(ri, cl, cnl, fi, dfi);
						s += a*qw[ik]*fi;
						ds += r*a*qw[ik]*dfi/ri;
					}
				}
			}

			template <class T>
			CGPU_EXEC_INL 
			T fcn_int_vz_ninfty_pinfty(const T& r, Ctpr<T> cl, Ctpr<T> cnl, const dt_int32& n_q, Ctpr<T> qx, Ctpr<T> qw, pFcn_clnl_3_1<T> fcn)
			{
				const T r2 = r*r;
				T s = 0;

				for(auto ik = 0; ik <n_q; ik++)
				{
					const auto zi = qx[ik];
					const auto ri = ::sqrt(zi*zi + r2);
					const auto fi = fcn(ri, cl, cnl);
					s += qw[ik]*fi;
				}

				return T(2)*s;
			}

			template <class T>
			CGPU_EXEC_INL 
			void fcn_int_vz_dvz_ninfty_pinfty(const T& r, Ctpr<T> cl, Ctpr<T> cnl, const dt_int32& n_q, Ctpr<T> qx, Ctpr<T> qw, T& s, T& ds, pFcn_clnl_3_2<T> fcn)
			{
				const T r2 = r*r;
				s = ds = 0;

				for(auto ik = 0; ik <n_q; ik++)
				{
					const auto zi = qx[ik];
					const auto ri = ::sqrt(zi*zi + r2);
					T fi, dfi;
					fcn(ri, cl, cnl, fi, dfi);
					s += qw[ik]*fi;
					ds += r*qw[ik]*dfi/ri;
				}
				s *= T(2);
				ds *= T(2);
			}

			/***************************************************************************************/
			/*************************************** macros ****************************************/
			/***************************************************************************************/
		#ifdef __CUDACC__
			#define pFCN_TEMPLATE_AF(cnam, TpFcn)															\
			template<eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T, eDev Dev>										\
			class pFcn_##cnam																				\
			{																								\
			public:																							\
				using pFcn = TpFcn<T>;																		\
																											\
				pFcn_##cnam()																				\
				{																							\
					if (Dev==edev_cpu)																		\
						pfcn = fcn_##cnam<atomic_pot_parm_typ, T>;													\
					else																					\
						cudaMemcpyFromSymbol(&pfcn, pgpu_fcn_##cnam<atomic_pot_parm_typ, T>, sizeof(pFcn));		\
				}																							\
																											\
				operator pFcn() 																			\
				{																							\
					return pfcn;																			\
				}																							\
			private:																						\
				pFcn pfcn;																					\
			}

			#define pFCN_TEMPLATE_AF_SPEC(cnam, TpFcn, ppt)													\
			template <class T, eDev Dev>																	\
			class pFcn_##cnam<ppt, T, Dev>																	\
			{																								\
			public:																							\
				using pFcn = TpFcn<T>;																		\
																											\
				pFcn_##cnam()																				\
				{																							\
					if (Dev==edev_cpu)																		\
						pfcn = fcn_##cnam<ppt, T>;															\
					else																					\
						cudaMemcpyFromSymbol(&pfcn, pgpu_fcn_##cnam##_spec<ppt, T>, sizeof(pFcn));			\
				}																							\
																											\
				operator pFcn() 																			\
				{																							\
					return pfcn;																			\
				}																							\
			private:																						\
				pFcn pfcn;																					\
			}
		#else
			#define pFCN_TEMPLATE_AF(cnam, TpFcn)															\
			template<eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T, eDev Dev>										\
			class pFcn_##cnam																				\
			{																								\
			public:																							\
				using pFcn = TpFcn<T>;																		\
																											\
				pFcn_##cnam()																				\
				{																							\
					pfcn = fcn_##cnam<atomic_pot_parm_typ, T>;														\
				}																							\
																											\
				operator pFcn() 																			\
				{																							\
					return pfcn;																			\
				}																							\
			private:																						\
				pFcn pfcn;																					\
			}

			#define pFCN_TEMPLATE_AF_SPEC(cnam, TpFcn, ppt)													\
			template <class T, eDev Dev>																	\
			class pFcn_##cnam<ppt, T, Dev>																	\
			{																								\
			public:																							\
				using pFcn = TpFcn<T>;																		\
																											\
				pFcn_##cnam()																				\
				{																							\
					pfcn = fcn_##cnam<ppt, T>;																\
				}																							\
																											\
				operator pFcn() 																			\
				{																							\
					return pfcn;																			\
				}																							\
			private:																						\
				pFcn pfcn;																					\
			}
		#endif

			/***************************************************************************************/
			/*************************************** feg *******************************************/
			/***************************************************************************************/
			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_doyle_0_4<atomic_pot_parm_typ, T>
			fcn_feg(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_gauss_v(x, 0, 4, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_doyle_0_4<atomic_pot_parm_typ, void>
			fcn_feg_dfeg(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_gauss_vd(x, 0, 4, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_4<atomic_pot_parm_typ, T>
			fcn_feg(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_gauss_v(x, 0, 5, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_4<atomic_pot_parm_typ, void>
			fcn_feg_dfeg(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_gauss_vd(x, 0, 5, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_12<atomic_pot_parm_typ, T>
			fcn_feg(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_gauss_v(x, 0, 5, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_12<atomic_pot_parm_typ, void>
			fcn_feg_dfeg(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_gauss_vd(x, 0, 5, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_kirkland_0_12<atomic_pot_parm_typ, T>
			fcn_feg(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_lorentzian_v(x, 0, 3, cl, cnl) + fcn_spt_gauss_v(x, 3, 6, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_kirkland_0_12<atomic_pot_parm_typ, void>
			fcn_feg_dfeg(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_lorentzian_vd(x, 0, 3, cl, cnl, y, dy);
				fcn_spt_gauss_vd<T, false>(x, 3, 6, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_weickenmeier_0_12<atomic_pot_parm_typ, T>
			fcn_feg(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				const T x2 = x*x;

				T y = 0;

				if (fcn_is_zero(x))
				{
					for(auto ik = 0; ik <6; ik++)
					{
						y += cl[ik]*cnl[ik];
					}

				}
				else
				{
					for(auto ik = 0; ik <6; ik++)
					{
						y += cl[ik]*(T(1)-exp(-cnl[ik]*x2));
					}
					y = y/x2;
				}

				return y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_weickenmeier_0_12<atomic_pot_parm_typ, void>
			fcn_feg_dfeg(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				const T x2 = x*x;

				y = dy = 0;

				if (fcn_is_zero(x))
				{
					for(auto ik = 0; ik <6; ik++)
					{
						y += cl[ik]*cnl[ik];
					}
				}
				else
				{
					for(auto ik = 0; ik <6; ik++)
					{
						const auto t = exp(-cnl[ik]*x2);
						y += cl[ik]*(T(1)-t);
						dy += cl[ik]*(T(1)-(T(1)+cnl[ik]*x2)*t);
					}
					y = y/x2;
					dy = T(-2)*dy/(x*x2);
				}
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_lobato_0_12<atomic_pot_parm_typ, T>
			fcn_feg(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				const T x2 = x*x;

				T y = 0;

				for(auto ik = 0; ik <5; ik++)
				{
					const auto t = T(1)/(T(1) + cnl[ik]*x2);
					y += cl[ik]*t*(t + T(1));
				}

				return y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_lobato_0_12<atomic_pot_parm_typ, void>
			fcn_feg_dfeg(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				const T x2 = x*x;

				y = dy = 0;

				for(auto ik = 0; ik <5; ik++)
				{
					const auto t = T(1)/(T(1) + cnl[ik]*x2);
					y += cl[ik]*t*(t + T(1));
					dy += T(-2)*x*cl[ik]*cnl[ik]*t*t*(T(2)*t + T(1));
				}
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_ion_0_4<atomic_pot_parm_typ, T>
			fcn_feg(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				auto y = fcn_spt_gauss_v(x, 0, 5, cl, cnl);
				y += cl[5]/(cnl[5] + x*x);

				return y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_ion_0_4<atomic_pot_parm_typ, void>
			fcn_feg_dfeg(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_gauss_vd(x, 0, 5, cl, cnl, y, dy);
				const auto t = T(1)/(cnl[5] + x*x);
				const auto yt = cl[5]*t;
				y += yt;
				dy += T(-2)*x*yt*t;
			}

			/***************************************************************************************/
			template <class T>
			using pFcn_feg1 = pFcn_clnl_3_1<T>;

		#ifdef __CUDACC__
			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_feg1<T> pgpu_fcn_feg = fcn_feg<atomic_pot_parm_typ, T>;
		#endif

			pFCN_TEMPLATE_AF(feg, pFcn_feg1);

			/***************************************************************************************/
			template <class T>
			using pFcn_feg2 = pFcn_clnl_3_2<T>;

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_feg2<T> pgpu_fcn_feg_dfeg = fcn_feg_dfeg<atomic_pot_parm_typ, T>;

			pFCN_TEMPLATE_AF(feg_dfeg , pFcn_feg2);

			/***************************************************************************************/
			/***************************************** fxg *****************************************/
			/***************************************************************************************/
			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_doyle_0_4<atomic_pot_parm_typ, T>
			fcn_fxg(const T& x, const T& Z, Ctpr<T> cl, Ctpr<T> cnl)
			{
				auto y = fcn_feg<atomic_pot_parm_typ, T>(x, cl, cnl);
				y = Z - x*x*y;

				return y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_doyle_0_4<atomic_pot_parm_typ, void>
			fcn_fxg_dfxg(const T& x, const T& Z, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_feg_dfeg<atomic_pot_parm_typ, T>(x, cl, cnl, y, dy);
				dy = -x*(T(2)*y + x*dy);
				y = Z - x*x*y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_4<atomic_pot_parm_typ, T>
			fcn_fxg(const T& x, const T& Z, Ctpr<T> cl, Ctpr<T> cnl)
			{
				auto y = fcn_feg<atomic_pot_parm_typ, T>(x, cl, cnl);
				y = Z - x*x*y;

				return y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_4<atomic_pot_parm_typ, void>
			fcn_fxg_dfxg(const T& x, const T& Z, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_feg_dfeg<atomic_pot_parm_typ, T>(x, cl, cnl, y, dy);
				dy = -x*(T(2)*y + x*dy);
				y = Z - x*x*y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_12<atomic_pot_parm_typ, T>
			fcn_fxg(const T& x, const T& Z, Ctpr<T> cl, Ctpr<T> cnl)
			{
				auto y = fcn_feg<atomic_pot_parm_typ, T>(x, cl, cnl);
				y = Z - x*x*y;

				return y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_12<atomic_pot_parm_typ, void>
			fcn_fxg_dfxg(const T& x, const T& Z, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_feg_dfeg<atomic_pot_parm_typ, T>(x, cl, cnl, y, dy);
				dy = -x*(T(2)*y + x*dy);
				y = Z - x*x*y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_kirkland_0_12<atomic_pot_parm_typ, T>
			fcn_fxg(const T& x, const T& Z, Ctpr<T> cl, Ctpr<T> cnl)
			{
				auto y = fcn_feg<atomic_pot_parm_typ, T>(x, cl, cnl);
				y = Z - x*x*y;

				return y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_kirkland_0_12<atomic_pot_parm_typ, void>
			fcn_fxg_dfxg(const T& x, const T& Z, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_feg_dfeg<atomic_pot_parm_typ, T>(x, cl, cnl, y, dy);
				dy = -x*(T(2)*y + x*dy);
				y = Z - x*x*y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_weickenmeier_0_12<atomic_pot_parm_typ, T>
			fcn_fxg(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_gauss_v(x, 0, 6, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_weickenmeier_0_12<atomic_pot_parm_typ, void>
			fcn_fxg_dfxg(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_gauss_vd(x, 0, 6, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_lobato_0_12<atomic_pot_parm_typ, T>
			fcn_fxg(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				const T x2 = x*x;

				T y = 0;

				for(auto ik = 0; ik <5; ik++)
				{
					const auto t = T(1)/(T(1)+cnl[ik]*x2);
					y += cl[ik]*t*t;
				}

				return y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_lobato_0_12<atomic_pot_parm_typ, void>
			fcn_fxg_dfxg(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				const T x2 = x*x;

				y = dy = 0;

				for(auto ik = 0; ik <5; ik++)
				{
					const auto t = T(1)/(T(1)+cnl[ik]*x2);
					const auto yt = cl[ik]*t*t;
					y += yt;
					dy += cnl[ik]*yt*t;
				}
				dy = T(-4)*x*dy;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_ion_0_4<atomic_pot_parm_typ, T>
			fcn_fxg(const T& x, const T& Z, Ctpr<T> cl, Ctpr<T> cnl)
			{
				auto y = fcn_feg<atomic_pot_parm_typ, T>(x, cl, cnl);
				y = Z - x*x*y;

				return y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_ion_0_4<atomic_pot_parm_typ, void>
			fcn_fxg_dfxg(const T& x, const T& Z, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_feg_dfeg<atomic_pot_parm_typ, T>(x, cl, cnl, y, dy);
				dy = -x*(T(2)*y + x*dy);
				y = Z - x*x*y;
			}

			/***************************************************************************************/
			template <class T>
			using pFcn_fxg1 = pFcn_clnl_4_1<T>;

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_fxg1<T> pgpu_fcn_fxg = fcn_fxg<atomic_pot_parm_typ, T>;

			pFCN_TEMPLATE_AF(fxg, pFcn_fxg1);

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_clnl_3_1<T> pgpu_fcn_fxg_spec = fcn_fxg<atomic_pot_parm_typ, T>;

			pFCN_TEMPLATE_AF_SPEC(fxg, pFcn_clnl_3_1, eappt_weickenmeier_0_12);

			pFCN_TEMPLATE_AF_SPEC(fxg, pFcn_clnl_3_1, eappt_lobato_0_12);

			/***************************************************************************************/
			template <class T>
			using pFcn_fxg2 = pFcn_clnl_4_2<T>;

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_fxg2<T> pgpu_fcn_fxg_dfxg = fcn_fxg_dfxg<atomic_pot_parm_typ, T>;

			pFCN_TEMPLATE_AF(fxg_dfxg, pFcn_fxg2);

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_clnl_3_2<T> pgpu_fcn_fxg_dfxg_spec = fcn_fxg_dfxg<atomic_pot_parm_typ, T>;

			pFCN_TEMPLATE_AF_SPEC(fxg_dfxg, pFcn_clnl_3_2, eappt_weickenmeier_0_12);

			pFCN_TEMPLATE_AF_SPEC(fxg_dfxg, pFcn_clnl_3_2, eappt_lobato_0_12);

			/***************************************************************************************/
			/****************************************** pr *****************************************/
			/***************************************************************************************/
			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_doyle_0_4<atomic_pot_parm_typ, T>
			fcn_pr(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_pr_gauss_feg_v(x, 0, 4, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_doyle_0_4<atomic_pot_parm_typ, void>
			fcn_pr_dpr(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_pr_gauss_feg_vd(x, 0, 4, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_4<atomic_pot_parm_typ, T>
			fcn_pr(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_pr_gauss_feg_v(x, 0, 5, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_4<atomic_pot_parm_typ, void>
			fcn_pr_dpr(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_pr_gauss_feg_vd(x, 0, 5, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_12<atomic_pot_parm_typ, T>
			fcn_pr(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_gauss_v(x, 0, 5, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_12<atomic_pot_parm_typ, void>
			fcn_pr_dpr(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_pr_gauss_feg_vd(x, 0, 5, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_kirkland_0_12<atomic_pot_parm_typ, T>
			fcn_pr(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_yukawa_v(x, 0, 3, cl, cnl) + fcn_spt_pr_gauss_feg_v(x, 3, 6, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_kirkland_0_12<atomic_pot_parm_typ, void>
			fcn_pr_dpr(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_yukawa_vd(x, 0, 3, cl, cnl, y, dy);
				fcn_spt_pr_gauss_feg_vd<T, false>(x, 3, 6, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_weickenmeier_0_12<atomic_pot_parm_typ, T>
			fcn_pr(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_gauss_v(x, 0, 6, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_weickenmeier_0_12<atomic_pot_parm_typ, void>
			fcn_pr_dpr(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_gauss_vd(x, 0, 6, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_lobato_0_12<atomic_pot_parm_typ, T>
			fcn_pr(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_exp_v(x, 0, 5, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_lobato_0_12<atomic_pot_parm_typ, void>
			fcn_pr_dpr(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_exp_vd(x, 0, 5, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_ion_0_4<atomic_pot_parm_typ, T>
			fcn_pr(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_pr_gauss_feg_v(x, 0, 5, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_ion_0_4<atomic_pot_parm_typ, void>
			fcn_pr_dpr(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_pr_gauss_feg_vd(x, 0, 5, cl, cnl, y, dy);
			}

			/***************************************************************************************/
			template <class T>
			using pFcn_pr1 = pFcn_clnl_3_1<T>;

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_pr1<T> pgpu_fcn_pr = fcn_pr<atomic_pot_parm_typ, T>;

			pFCN_TEMPLATE_AF(pr, pFcn_pr1);

			/***************************************************************************************/
			template <class T>
			using pFcn_pr2 = pFcn_clnl_3_2<T>;

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_pr2<T> pgpu_fcn_pr_dpr = fcn_pr_dpr<atomic_pot_parm_typ, T>;

			pFCN_TEMPLATE_AF(pr_dpr , pFcn_pr2);

			/***************************************************************************************/
			/***************************************** vr ******************************************/
			/***************************************************************************************/
			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_doyle_0_4<atomic_pot_parm_typ, T>
			fcn_vr(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_gauss_v(x, 0, 4, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_doyle_0_4<atomic_pot_parm_typ, void>
			fcn_vr_dvr(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_gauss_vd(x, 0, 4, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_4<atomic_pot_parm_typ, T>
			fcn_vr(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_gauss_v(x, 0, 5, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_4<atomic_pot_parm_typ, void>
			fcn_vr_dvr(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_gauss_vd(x, 0, 5, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_12<atomic_pot_parm_typ, T>
			fcn_vr(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_gauss_v(x, 0, 5, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_12<atomic_pot_parm_typ, void>
			fcn_vr_dvr(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_gauss_vd(x, 0, 5, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_kirkland_0_12<atomic_pot_parm_typ, T>
			fcn_vr(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_yukawa_v(x, 0, 3, cl, cnl) + fcn_spt_gauss_v(x, 3, 6, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_kirkland_0_12<atomic_pot_parm_typ, void>
			fcn_vr_dvr(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_yukawa_vd(x, 0, 3, cl, cnl, y, dy);
				fcn_spt_gauss_vd<T, false>(x, 3, 6, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_weickenmeier_0_12<atomic_pot_parm_typ, T>
			fcn_vr(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				const T ix = T(1)/x;

				T y = 0;

				for(auto ik = 0; ik <6; ik++)
				{
					y += cl[ik]*erfc(cnl[ik]*x)*ix;
				}

				return y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_weickenmeier_0_12<atomic_pot_parm_typ, void>
			fcn_vr_dvr(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				const T c_pii2 = 1.772453850905516027298;
				const T ix = T(1)/x;
				const T x2 = x*x;

				y = dy = 0;

				for(auto ik = 0; ik <6; ik++)
				{
					const auto yt = cl[ik]*erfc(cnl[ik]*x)*ix;
					y += yt;
					dy += (T(-2)*cl[ik]*cnl[ik]*exp(-cnl[ik]*cnl[ik]*x2)/c_pii2-yt)*ix;
				}
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_lobato_0_12<atomic_pot_parm_typ, T>
			fcn_vr(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				const T ix = T(1)/x;

				T y = 0;

				for(auto ik = 0; ik <5; ik++)
				{
					y += cl[ik]*exp(-cnl[ik]*x)*ix*(T(2)/cnl[ik] + x);
				}

				return y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_lobato_0_12<atomic_pot_parm_typ, void>
			fcn_vr_dvr(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				const T ix = T(1)/x;

				y = dy = 0;

				for(auto ik = 0; ik <5; ik++)
				{
					const auto yt = cl[ik]*exp(-cnl[ik]*x)*ix;
					const auto icnl = T(1)/cnl[ik];
					y += yt*(T(2)*icnl + x);
					dy += -yt*(T(2)*icnl*ix + T(2) + cnl[ik]*x);
				}
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_ion_0_4<atomic_pot_parm_typ, T>
			fcn_vr(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				auto y = fcn_spt_gauss_v(x, 0, 5, cl, cnl);
				y += cl[5]*exp(-cnl[5]*x)/x;

				return y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_ion_0_4<atomic_pot_parm_typ, void>
			fcn_vr_dvr(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_gauss_vd(x, 0, 5, cl, cnl, y, dy);
				const auto ix = T(1)/x;
				const auto yt = cl[5]*exp(-cnl[5]*x)*ix;
				y += yt;
				dy += -(cnl[5]+ ix)*yt;
			}

			 /***************************************************************************************/
			template <class T>
			using pFcn_vr1 = pFcn_clnl_3_1<T>;

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_vr1<T> pgpu_fcn_vr = fcn_vr<atomic_pot_parm_typ, T>;

			pFCN_TEMPLATE_AF(vr, pFcn_vr1);

			 /***************************************************************************************/
			template <class T>
			using pFcn_vr2 = pFcn_clnl_3_2<T>;

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_vr2<T> pgpu_fcn_vr_dvr = fcn_vr_dvr<atomic_pot_parm_typ, T>;

			pFCN_TEMPLATE_AF(vr_dvr , pFcn_vr2);

			/***************************************************************************************/
			/***************************************** vz ******************************************/
			/***************************************************************************************/
			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL
			T fcn_vz(const T& x, const T& z_0, const T& z_e, Ctpr<T> cl, Ctpr<T> cnl, const dt_int32& n_q, Ctpr<T> qx, Ctpr<T> qw)
			{
				return fcn_int_vz_z0_ze(x, z_0, z_e, cl, cnl, n_q, qx, qw, fcn_vr<atomic_pot_parm_typ, T>);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL
			void fcn_vz_dvz(const T& x, const T& z_0, const T& z_e, Ctpr<T> cl, Ctpr<T> cnl, const dt_int32& n_q, Ctpr<T> qx, Ctpr<T> qw, T& y, T& dy)
			{
				fcn_int_vz_dvz_z0_ze(x, z_0, z_e, cl, cnl, n_q, qx, qw, y, dy, fcn_vr_dvr<atomic_pot_parm_typ, T>);
			}

			 /***************************************************************************************/
			template <class T>
			using pFcn_vz1 = pFcn_clnl_8_1<T>;

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_vz1<T> pgpu_fcn_vz = fcn_vz<atomic_pot_parm_typ, T>;

			pFCN_TEMPLATE_AF(vz, pFcn_vz1);

			 /***************************************************************************************/
			template <class T>
			using pFcn_vz2 = pFcn_clnl_8_2<T>;

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_vz2<T> pgpu_fcn_vz_dvz = fcn_vz_dvz<atomic_pot_parm_typ, T>;

			pFCN_TEMPLATE_AF(vz_dvz , pFcn_vz2);

			 /***************************************************************************************/
			/****************************************** vzp *****************************************/
			 /***************************************************************************************/
			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_doyle_0_4<atomic_pot_parm_typ, T>
			fcn_vzp(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_gauss_v(x, 0, 4, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_doyle_0_4<atomic_pot_parm_typ, void>
			fcn_vzp_dvzp(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_gauss_vd(x, 0, 4, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_4<atomic_pot_parm_typ, T>
			fcn_vzp(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				return fcn_spt_gauss_v(x, 0, 5, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_4<atomic_pot_parm_typ, void>
			fcn_vzp_dvzp(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_gauss_vd(x, 0, 5, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_12<atomic_pot_parm_typ, T>
			fcn_vzp(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				 return fcn_spt_gauss_v(x, 0, 5, cl, cnl);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_0_12<atomic_pot_parm_typ, void>
			fcn_vzp_dvzp(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_gauss_vd(x, 0, 5, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_kirkland_0_12<atomic_pot_parm_typ, T>
			fcn_vzp(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				T y = 0;

				for(auto ik = 0; ik <3; ik++)
				{
					y += cl[ik]*bessel_k0(cnl[ik]*x);
				}

				y += fcn_spt_gauss_v(x, 3, 6, cl, cnl);

				return y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_kirkland_0_12<atomic_pot_parm_typ, void>
			fcn_vzp_dvzp(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				y = dy = 0;

				for(auto ik = 0; ik <3; ik++)
				{
					y += cl[ik]*bessel_k0(cnl[ik]*x);
					dy += -cl[ik]*cnl[ik]*bessel_k1(cnl[ik]*x);
				}

				fcn_spt_gauss_vd<T, false>(x, 3, 6, cl, cnl, y, dy);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_weickenmeier_0_12<atomic_pot_parm_typ, T>
			fcn_vzp(const T& x, Ctpr<T> cl, Ctpr<T> cnl, const dt_int32& n_q, Ctpr<T> qx, Ctpr<T> qw)
			{
				return fcn_int_vz_ninfty_pinfty(x, cl, cnl, n_q, qx, qw, fcn_vr<atomic_pot_parm_typ, T>);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_weickenmeier_0_12<atomic_pot_parm_typ, void>
			fcn_vzp_dvzp(const T& x, Ctpr<T> cl, Ctpr<T> cnl, const dt_int32& n_q, Ctpr<T> qx, Ctpr<T> qw, T& y, T& dy)
			{
				fcn_int_vz_dvz_ninfty_pinfty(x, cl, cnl, n_q, qx, qw, y, dy, fcn_vr_dvr<atomic_pot_parm_typ, T>);
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_lobato_0_12<atomic_pot_parm_typ, T>
			fcn_vzp(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				T y = 0;

				for(auto ik = 0; ik <5; ik++)
				{
					y += cl[ik]*(T(2)*bessel_k0(cnl[ik]*x)/cnl[ik] + x*bessel_k1(cnl[ik]*x));
				}

				return y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_lobato_0_12<atomic_pot_parm_typ, void>
			fcn_vzp_dvzp(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				y = dy = 0;

				for(auto ik = 0; ik <5; ik++)
				{
					const auto k0 = bessel_k0(cnl[ik]*x);
					const auto k1 = bessel_k1(cnl[ik]*x);
					y += cl[ik]*(T(2)*k0/cnl[ik] + x*k1);
					dy += -cl[ik]*(cnl[ik]*x*k0 + T(2)*k1);
				}
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_ion_0_4<atomic_pot_parm_typ, T>
			fcn_vzp(const T& x, Ctpr<T> cl, Ctpr<T> cnl)
			{
				auto y = fcn_spt_gauss_v(x, 0, 5, cl, cnl);
				y += cl[5]*bessel_k0(cnl[5]*x);

				return y;
			}

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			CGPU_EXEC_INL 
			enable_if_eappt_peng_ion_0_4<atomic_pot_parm_typ, void>
			fcn_vzp_dvzp(const T& x, Ctpr<T> cl, Ctpr<T> cnl, T& y, T& dy)
			{
				fcn_spt_gauss_vd(x, 0, 5, cl, cnl, y, dy);
				y += cl[5]*bessel_k0(cnl[5]*x);
				dy += -cl[5]*cnl[5]*bessel_k1(cnl[5]*x);
			}

			/***************************************************************************************/
			template <class T>
			using pFcn_vzp1 = pFcn_clnl_3_1<T>;

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_vzp1<T> pgpu_fcn_vzp = fcn_vzp<atomic_pot_parm_typ, T>;

			pFCN_TEMPLATE_AF(vzp, pFcn_vzp1);

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_clnl_6_1<T> pgpu_fcn_vzp_spec = fcn_vzp<atomic_pot_parm_typ, T>;

			pFCN_TEMPLATE_AF_SPEC(vzp, pFcn_clnl_6_1, eappt_weickenmeier_0_12);

			/***************************************************************************************/
			template <class T>
			using pFcn_vzp2 = pFcn_clnl_3_2<T>;

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_vzp2<T> pgpu_fcn_vzp_dvzp = fcn_vzp_dvzp<atomic_pot_parm_typ, T>;

			pFCN_TEMPLATE_AF(vzp_dvzp, pFcn_vzp2);

			template <eAtomic_Pot_Parm_Typ atomic_pot_parm_typ, class T>
			GPU_EXEC pFcn_clnl_6_2<T> pgpu_fcn_vzp_dvzp_spec = fcn_vzp_dvzp<atomic_pot_parm_typ, T>;

			pFCN_TEMPLATE_AF_SPEC(vzp_dvzp, pFcn_clnl_6_2, eappt_weickenmeier_0_12);
		}
	}

	/* detector integration */
	namespace mt
	{
		namespace cgpu_detail_mt
		{
			/* integration over a detector ring */
			template <class T, class U>
			CGPU_EXEC_INL 
			void fcn_int_det_ring(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
			const T& g2_min, const T& g2_max, Ctpr<U> mx_i, KS<U>& sum)
			{
				const auto g2 = grid.g2_sft(ix, iy);
				if (mt::fcn_chk_bound(g2, g2_min, g2_max))
				{
					const auto ixy = grid.sub_2_ind(ix, iy);
					sum += mx_i[ixy];
				}
			}

			/* norm_2 integration over a detector ring */
			template <class T, class U>
			CGPU_EXEC_INL 
			void fcn_int_det_ring_norm_2(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
			const T& g2_min, const T& g2_max, Ctpr<U> mx_i, KS<Value_type_r<U>>& sum)
			{
				const auto g2 = grid.g2_sft(ix, iy);
				if (mt::fcn_chk_bound(g2, g2_min, g2_max))
				{
					const auto ixy = grid.sub_2_ind(ix, iy);
					sum += ::norm_2(mx_i[ixy]);
				}
			}

			/* integration over a detector with sensitivity */
			template <class T, class U>
			CGPU_EXEC_INL 
			void fcn_int_det_sen(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
			Ctpr<U> sen_i, Ctpr<U> mx_i, KS<U>& sum)
			{
				const auto ixy = grid.sub_2_ind(ix, iy);
				sum += sen_i[ixy]*mx_i[ixy];
			}

			/* norm_2 integration over a detector with sensitivity */
			template <class T, class U>
			CGPU_EXEC_INL 
			void fcn_int_det_sen_norm_2(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
			Ctpr<Value_type_r<U>> sen_i, Ctpr<U> mx_i, KS<Value_type_r<U>>& sum)
			{
				const auto ixy = grid.sub_2_ind(ix, iy);
				sum += sen_i[ixy]*::norm_2(mx_i[ixy]);
			}
		}
	}

	/* wave propagation */
	namespace mt
	{
		namespace cgpu_detail_mt
		{
			/* propagate */
			template <class T>
			CGPU_EXEC_INL 
			void fcn_fs_propagate(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
				complex<T>* psi_i, const R_2d<T>& g_0, const T& w_g2, const T& w, complex<T>* psi_o)
			{
				const auto ixy = grid.sub_2_ind(ix, iy);
				const auto theta = w_g2*grid.g2_sft(ix, iy, g_0);

				psi_o[ixy] = w*psi_i[ixy]*euler(theta);
			}

			/* propagate and bandwith limit using a fermi aperture */
			template <class T>
			CGPU_EXEC_INL 
			void fcn_fs_propagate_bw_f(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
				complex<T>* psi_i, const R_2d<T>& g_0, const T& w_g2, const T& g2_cut, const T& alpha, const T& w, complex<T>* psi_o)
			{
				const auto ixy = grid.sub_2_ind(ix, iy);
				const auto theta = w_g2*grid.g2_sft(ix, iy, g_0);
				const auto m = w*fcn_fermi_lpf(alpha, g2_cut, grid.g2_sft(ix, iy)); // bandwith limit does not mater if it includes till illumination
				psi_o[ixy] = m*psi_i[ixy]*euler(theta);
			}

			/* propagate and bandwith limit using a hard aperture */
			template <class T>
			CGPU_EXEC_INL 
			void fcn_fs_propagate_bw_h(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
				complex<T>* psi_i, const R_2d<T>& g_0, const T& w_g2, const T& g2_cut, const T& w, complex<T>* psi_o)
			{
				const auto ixy = grid.sub_2_ind(ix, iy);
				const auto theta = w_g2*grid.g2_sft(ix, iy, g_0);
				const auto m = ((grid.g2_sft(ix, iy) <= g2_cut))?w:T(0);	// bandwith limit does not mater if it includes till illumination
				psi_o[ixy] = m*psi_i[ixy]*euler(theta);
			}
		}
	}

	/* probe - ctf - pctf */
	namespace mt
	{
		namespace cgpu_detail_mt
		{
			/* create probe */
			template <bool bb_phi, class T>
			CGPU_EXEC_INL 
			complex<T> fcn_fs_exp_i_chi(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, const Lens<T>& lens, 
			const R_2d<T>& R, const R_2d<T>& gu)
			{
				const auto gx = grid.gx_sft(ix) + gu.x;
				const auto gy = grid.gy_sft(iy) + gu.y;
				const auto g2 = gx*gx + gy*gy;

				complex<T> v = 0;

				if (fcn_chk_bound(g2, lens.g2_inner, lens.g2_outer))
				{
					const auto g4 = g2*g2;
					const auto g6 = g4*g2;
					const auto chi = R.x*gx + R.y*gy + lens.eval_c_10(g2) + lens.eval_c_30(g4) + lens.eval_c_50(g6);

					if (bb_phi)
					{
						const auto g =::sqrt(g2);
						const auto g3 = g2*g;
						const auto g5 = g4*g;
						const auto phi = atan2(gy, gx);
						chi += lens.eval_m(phi) + lens.eval_c_12(g2, phi);
						chi += lens.eval_c_21_c_23(g3, phi) + lens.eval_c_32_c_34(g4, phi);
						chi += lens.eval_c_41_c_43_c_45(g5, phi) + lens.eval_c_52_c_54_c_56(g6, phi);
					}

					v = euler(chi);
				}

				return v;
			}

			/* create probe */
			template <bool bb_phi, class T>
			CGPU_EXEC_INL 
			void fcn_fs_probe(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, const Lens<T>& lens, 
			const R_2d<T>& r, const R_2d<T>& gu, complex<T>* psi_o)
			{
				const auto v = fcn_fs_exp_i_chi<bb_phi>(ix, iy, grid, lens, r, gu);

				const auto ixy = grid.sub_2_ind(ix, iy);
				psi_o[ixy] = v;
			}

			/* apply coherent transfer function */
			template <bool bb_phi, class T>
			CGPU_EXEC_INL 
			void fcn_fs_apply_ctf(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, complex<T>* psi_i, 
			const Lens<T>& lens, const R_2d<T>& gu, complex<T>* psi_o)
			{
				const auto v = fcn_fs_exp_i_chi<bb_phi>(ix, iy, grid, lens, R_2d<T>(), gu);

				const auto ixy = grid.sub_2_ind(ix, iy);
				psi_o[ixy] = v*psi_i[ixy];
			}

			/* apply partial coherent transfer function */
			// Marc De Graef - Introduction to Conventional Transmission Electron Microscopy page: 611
			// tp_inc_iehwgd: 08, spt_inc_theta_c: 611
			template <class T>
			CGPU_EXEC_INL 
			void fcn_fs_apply_pctf(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, complex<T>* psi_i, 
			const Lens<T>& lens, complex<T>* psi_o)
			{
				const auto ixy = grid.sub_2_ind(ix, iy);
				const auto g2 = grid.g2_sft(ix, iy);

				if (fcn_chk_bound(g2, lens.g2_inner, lens.g2_outer))
				{		
					const auto c_pi = T(3.141592653589793238463);

					const auto chi = lens.eval_c_10(g2) + lens.eval_c_30(g2*g2);
					const auto c_u = c_pi*lens.spt_inc_theta_c*lens.tp_inc_iehwgd;
					const auto u = T(1) + T(2)*c_u*c_u*g2;

					const auto c_tp = c_pi*lens.tp_inc_iehwgd*lens.lambda*g2;
					const auto tp_inc = T(0.5)*c_tp*c_tp*g2*g2;

					const auto c_spt = c_pi*lens.spt_inc_iehwgd*(lens.c_30*lens.lambda_2*g2-lens.c_10);
					const auto spt_inc = c_spt*c_spt*g2;

					const auto st_inc = exp(-(spt_inc + tp_inc)/u); // sqrt(u);

					psi_o[ixy] = psi_i[ixy]*polar(st_inc, chi);
				}
				else
				{
					psi_o[ixy] = T(0);
				}
			}
		}
	}

	/* transmission function */
	namespace mt
	{
		namespace cgpu_detail_mt
		{
			template <eElec_Spec_Interact_Mod esim, class T>
			CGPU_EXEC_INL 
			void fcn_trans_fcn(const dt_int32& ix, const iGrid_1d& igrid, Ctpr<T>* vzp_i, const T& w, complex<T>* tfcn_o)
			{
				const auto ind = igrid.sub_2_ind(ix);

				if (esim == eesim_weak_phase_object)
				{
					tfcn_o[ind] =complex<T>(T(1), w*vzp_i[ind]);
				}
				else
				{
					tfcn_o[ind] =euler(w*vzp_i[ind]);
				}
			}
		}
	}

	/* eels */
	namespace mt
	{
		namespace cgpu_detail_mt
		{
			template <class T>
			CGPU_EXEC_INL 
			void fcn_eels_lorentz_norm_factor(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
				const T& gc2, const T& ge2, KS<T>& sum)
			{
				const auto g2 = grid.g2_sft(ix, iy);
				if (g2 < gc2)
				{
					sum += T(1)/(g2 + ge2);
				}
			}

			template <class T>
			CGPU_EXEC_INL 
			void fcn_eels_w_xyz(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
				const EELS<T>& eels, Tpr<complex<T>> w_x, Tpr<complex<T>> w_y, Tpr<complex<T>> w_z)
			{
				const auto ixy = grid.sub_2_ind(ix, iy);
				const auto gx = grid.gx_sft(ix);
				const auto gy = grid.gy_sft(iy);
				const auto g2 = gx*gx + gy*gy;
				
				if (g2 < eels.gc2)
				{
					const auto pos = euler(eels.x*gx + eels.y*gy);
					const auto lorentz = eels.factor/(g2 + eels.ge2);
					w_x[ixy] = gx*lorentz*pos;
					w_y[ixy] = gy*lorentz*pos;
					w_z[ixy] = eels.ge*lorentz*pos;
				}
				else
				{
					w_x[ixy] = complex<T>(0);
					w_y[ixy] = complex<T>(0);
					w_z[ixy] = complex<T>(0);
				}
			}

			template <class T>
			CGPU_EXEC_INL 
			void fcn_eels_w_x(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, const EELS<T>& eels, Tpr<complex<T>> w_x)
			{
				const auto ixy = grid.sub_2_ind(ix, iy);
				const auto gx = grid.gx_sft(ix);
				const auto gy = grid.gy_sft(iy);
				const auto g2 = gx*gx + gy*gy;
				
				if (g2 < eels.gc2)
				{
					const auto pos = euler(eels.x*gx + eels.y*gy);
					const auto lorentz = eels.factor/(g2 + eels.ge2);
					w_x[ixy] = gx*lorentz*pos;
				}
				else
				{
					w_x[ixy] = complex<T>(0);
				}
			}

			template <class T>
			CGPU_EXEC_INL 
			void fcn_eels_w_y(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, const EELS<T>& eels, Tpr<complex<T>> w_y)
			{
				const auto ixy = grid.sub_2_ind(ix, iy);
				const auto gx = grid.gx_sft(ix);
				const auto gy = grid.gy_sft(iy);
				const auto g2 = gx*gx + gy*gy;
				
				if (g2 < eels.gc2)
				{
					const auto pos = euler(eels.x*gx + eels.y*gy);
					const auto lorentz = eels.factor/(g2 + eels.ge2);
					w_y[ixy] = gy*lorentz*pos;
				}
				else
				{
					w_y[ixy] = complex<T>(0);
				}
			}

			template <class T>
			CGPU_EXEC_INL 
			void fcn_eels_w_z(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, const EELS<T>& eels, Tpr<complex<T>> w_z)
			{
				const auto ixy = grid.sub_2_ind(ix, iy);
				const auto gx = grid.gx_sft(ix);
				const auto gy = grid.gy_sft(iy);
				const auto g2 = gx*gx + gy*gy;
				
				if (g2 < eels.gc2)
				{
					const auto pos = euler(eels.x*gx + eels.y*gy);
					const auto lorentz = eels.factor/(g2 + eels.ge2);
					w_z[ixy] = eels.ge*lorentz*pos;
				}
				else
				{
					w_z[ixy] = complex<T>(0);
				}
			}

			template <class T>
			CGPU_EXEC_INL 
			void fcn_eels_w_mn1(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, const EELS<T>& eels, Tpr<complex<T>> w_mn1)
			{
				const auto ixy = grid.sub_2_ind(ix, iy);
				const auto gx = grid.gx_sft(ix);
				const auto gy = grid.gy_sft(iy);
				const auto g2 = gx*gx + gy*gy;
				
				if (g2 < eels.gc2)
				{
					const auto c_i2i2 = T(0.70710678118654746);

					const auto pos = euler(eels.x*gx + eels.y*gy);
					const auto lorentz = c_i2i2*eels.factor/(g2 + eels.ge2);
					w_mn1[ixy] = complex<T>(gx, -gy)*lorentz*pos;
				}
				else
				{
					w_mn1[ixy] = complex<T>(0);
				}
			}

			template <class T>
			CGPU_EXEC_INL 
			void fcn_eels_w_mp1(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, const EELS<T>& eels, Tpr<complex<T>> w_mp1)
			{
				const auto ixy = grid.sub_2_ind(ix, iy);
				const auto gx = grid.gx_sft(ix);
				const auto gy = grid.gy_sft(iy);
				const auto g2 = gx*gx + gy*gy;
				
				if (g2 < eels.gc2)
				{
					const auto c_i2i2 = T(0.70710678118654746);

					const auto pos = euler(eels.x*gx + eels.y*gy);
					const auto lorentz = c_i2i2*eels.factor/(g2 + eels.ge2);
					w_mp1[ixy] = complex<T>(gx, gy)*lorentz*pos;
				}
				else
				{
					w_mp1[ixy] = complex<T>(0);
				}
			}
		}
	}

#endif
