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

#ifndef CPU_FCNS_MT_H
	#define CPU_FCNS_MT_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include "const_enum_mt.cuh"
	#include "math.cuh"
	#include "types_mt.cuh"
	#include "type_traits_mt.cuh"
	#include "cgpu_vctr.cuh"
	#include "cgpu_fft.cuh"
	#include "cgpu_stream.cuh"
	#include "quad_data.cuh"
	#include "cgpu_detail_mt.cuh"
	#include "cpu_fcns.hpp"

	/* atomic functions */
	namespace mt
	{
		/*********************************** pFcn_clnl_3_x<T> **********************************/
		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_cpu_32<T>& x, const pLNL_Coef_cpu<T>& coef, pVctr_cpu_32<T>& fx, pFcn_clnl_3_1<T> fcn)
		{
			for(auto ix = 0; ix < x.m_size; ix++)
			{
				fx[ix] = fcn(x[ix], coef.cl, coef.cnl);
			}
		}

		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_cpu_32<T>& x, const pLNL_Coef_cpu<T>& coef, pVctr_cpu_32<T>& fx1, pVctr_cpu_32<T>& fx2, pFcn_clnl_3_2<T> fcn)
		{
			for(auto ix = 0; ix < x.m_size; ix++)
			{
				fcn(x[ix], coef.cl, coef.cnl, fx1[ix], fx2[ix]);
			}
		}

		/********************************* pFcn_clnl_4_x<T> ************************************/
		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_cpu_32<T>& x, const T& c, const pLNL_Coef_cpu<T>& coef, pVctr_cpu_32<T>& fx, pFcn_clnl_4_1<T> fcn)
		{
			for(auto ix = 0; ix < x.m_size; ix++)
			{
				fx[ix] = fcn(x[ix], c, coef.cl, coef.cnl);
			}
		}
		
		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_cpu_32<T>& x, const T& c, const pLNL_Coef_cpu<T>& coef, pVctr_cpu_32<T>& fx1, pVctr_cpu_32<T>& fx2, pFcn_clnl_4_2<T> fcn)
		{
			for(auto ix = 0; ix < x.m_size; ix++)
			{
				fcn(x[ix], c, coef.cl, coef.cnl, fx1[ix], fx2[ix]);
			}
		}

		/********************************* pFcn_clnl_6_x<T> ************************************/
		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_cpu_32<T>& x, const pLNL_Coef_cpu<T>& coef, const pQuad_Coef_1d_cpu<T>& quad, pVctr_cpu_32<T>& fx, pFcn_clnl_6_1<T> fcn)
		{
			for(auto ix = 0; ix < x.m_size; ix++)
			{
				fx[ix] = fcn(x[ix], coef.cl, coef.cnl, quad.m_size, quad.x, quad.w);
			}
		}

		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_cpu_32<T>& x, const pLNL_Coef_cpu<T>& coef, const pQuad_Coef_1d_cpu<T>& quad, pVctr_cpu_32<T>& fx1, pVctr_cpu_32<T>& fx2, pFcn_clnl_6_2<T> fcn)
		{
			for(auto ix = 0; ix < x.m_size; ix++)
			{
				fcn(x[ix], coef.cl, coef.cnl, quad.m_size, quad.x, quad.w, fx1[ix], fx2[ix]);
			}
		}

		/********************************* pFcn_clnl_8_x<T> ************************************/
		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_cpu_32<T>& x, const T& z_0, const T& z_e, const pLNL_Coef_cpu<T>& coef, const pQuad_Coef_1d_cpu<T>& quad, pVctr_cpu_32<T>& fx, pFcn_clnl_8_1<T> fcn)
		{
			for(auto ix = 0; ix < x.m_size; ix++)
			{
				fx[ix] = fcn(x[ix], z_0, z_e, coef.cl, coef.cnl, quad.m_size, quad.x, quad.w);
			}
		}
				
		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_cpu_32<T>& x, const T& z_0, const T& z_e, const pLNL_Coef_cpu<T>& coef, const pQuad_Coef_1d_cpu<T>& quad, pVctr_cpu_32<T>& fx1, pVctr_cpu_32<T>& fx2, pFcn_clnl_8_2<T> fcn)
		{
			for(auto ix = 0; ix < x.m_size; ix++)
			{
				fcn(x[ix], z_0, z_e, coef.cl, coef.cnl, quad.m_size, quad.x, quad.w, fx1[ix], fx2[ix]);
			}
		}

		/***************************************************************************************/
		/***************************************************************************************/
		/***************************************************************************************/
		#define SWITCH_FCN_ATOMIC_V(atomic_pot_parm_typ, fcn, ...)										\
		switch(atomic_pot_parm_typ)																	\
		{																						\
			case eappt_doyle_0_4:																\
			{																					\
				return cgpu_detail_mt::fcn<eappt_doyle_0_4, T>(__VA_ARGS__);						\
			}																					\
			case eappt_peng_0_4:																	\
			{																					\
				return cgpu_detail_mt::fcn<eappt_peng_0_4, T>(__VA_ARGS__);						\
			}																					\
			case eappt_peng_0_12:																\
			{																					\
				return cgpu_detail_mt::fcn<eappt_peng_0_12, T>(__VA_ARGS__);						\
			}																					\
			case eappt_kirkland_0_12:															\
			{																					\
				return cgpu_detail_mt::fcn<eappt_kirkland_0_12, T>(__VA_ARGS__);					\
			}																					\
			case eappt_weickenmeier_0_12:														\
			{																					\
				return cgpu_detail_mt::fcn<eappt_weickenmeier_0_12, T>(__VA_ARGS__);				\
			}																					\
			case eappt_lobato_0_12:																\
			{																					\
				return cgpu_detail_mt::fcn<eappt_lobato_0_12, T>(__VA_ARGS__);					\
			}																					\
			case eappt_peng_ion_0_4:																\
			{																					\
				return cgpu_detail_mt::fcn<eappt_peng_ion_0_4, T>(__VA_ARGS__);					\
			}																					\
			default:																			\
			{																					\
				return 0;																		\
			}																					\
		}		

		#define SWITCH_FCN_ATOMIC_VD(atomic_pot_parm_typ, fcn, ...)									\
		switch(atomic_pot_parm_typ)																	\
		{																						\
			case eappt_doyle_0_4:																\
			{																					\
				cgpu_detail_mt::fcn<eappt_doyle_0_4, T>(__VA_ARGS__);							\
			}																					\
			break;																				\
			case eappt_peng_0_4:																	\
			{																					\
				cgpu_detail_mt::fcn<eappt_peng_0_4, T>(__VA_ARGS__);								\
			}																					\
			break;																				\
			case eappt_peng_0_12:																\
			{																					\
				cgpu_detail_mt::fcn<eappt_peng_0_12, T>(__VA_ARGS__);							\
			}																					\
			break;																				\
			case eappt_kirkland_0_12:															\
			{																					\
				cgpu_detail_mt::fcn<eappt_kirkland_0_12, T>(__VA_ARGS__);						\
			}																					\
			break;																				\
			case eappt_weickenmeier_0_12:														\
			{																					\
				cgpu_detail_mt::fcn<eappt_weickenmeier_0_12, T>(__VA_ARGS__);					\
			}																					\
			break;																				\
			case eappt_lobato_0_12:																\
			{																					\
				cgpu_detail_mt::fcn<eappt_lobato_0_12, T>(__VA_ARGS__);							\
			}																					\
			break;																				\
			case eappt_peng_ion_0_4:																\
			{																					\
				cgpu_detail_mt::fcn<eappt_peng_ion_0_4, T>(__VA_ARGS__);							\
			}																					\
			break;																				\
		}

		/************************************** feg ********************************************/
		template <class T>
		CPU_EXEC
		T fcn_feg(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const T& g, const pLNL_Coef_cpu<T>& coef)
		{
			SWITCH_FCN_ATOMIC_V(atomic_pot_parm_typ, fcn_feg, g, coef.cl, coef.cnl);
		}

		template <class T>
		CPU_EXEC
		void fcn_feg_dfeg(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const T& g, const pLNL_Coef_cpu<T>& coef, T& y, T& dy)
		{
			SWITCH_FCN_ATOMIC_VD(atomic_pot_parm_typ, fcn_feg_dfeg, g, coef.cl, coef.cnl, y, dy);
		}

		/*************************************** fxg *******************************************/
		template <class T>
		CPU_EXEC
		T fcn_fxg(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const T& g, const T& Z, const pLNL_Coef_cpu<T>& coef)
		{
			switch(atomic_pot_parm_typ)															
			{																				
				case eappt_doyle_0_4:														
				{																			
					return cgpu_detail_mt::fcn_fxg<eappt_doyle_0_4, T>(g, Z, coef.cl, coef.cnl);
				}																			
				case eappt_peng_0_4:															
				{																			
					return cgpu_detail_mt::fcn_fxg<eappt_peng_0_4, T>(g, Z, coef.cl, coef.cnl);
				}																			
				case eappt_peng_0_12:														
				{																			
					return cgpu_detail_mt::fcn_fxg<eappt_peng_0_12, T>(g, Z, coef.cl, coef.cnl);
				}																			
				case eappt_kirkland_0_12:													
				{																			
					return cgpu_detail_mt::fcn_fxg<eappt_kirkland_0_12, T>(g, Z, coef.cl, coef.cnl);
				}
				case eappt_weickenmeier_0_12:														
				{																			
					return cgpu_detail_mt::fcn_fxg<eappt_weickenmeier_0_12, T>(g, coef.cl, coef.cnl);
				}																			
				case eappt_lobato_0_12:															
				{																			
					return cgpu_detail_mt::fcn_fxg<eappt_lobato_0_12, T>(g, coef.cl, coef.cnl);
				}	
				case eappt_peng_ion_0_4:														
				{																			
					return cgpu_detail_mt::fcn_fxg<eappt_peng_ion_0_4, T>(g, Z, coef.cl, coef.cnl);
				}																			
				default:																	
				{																			
					return 0;																
				}																			
			}
		}

		template <class T>
		CPU_EXEC
		void fcn_fxg_dfxg(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const T& g, const T& Z, const pLNL_Coef_cpu<T>& coef, T& y, T& dy)
		{
			switch(atomic_pot_parm_typ)																
			{																					
				case eappt_doyle_0_4:															
				{																				
					cgpu_detail_mt::fcn_fxg_dfxg<eappt_doyle_0_4, T>(g, Z, coef.cl, coef.cnl, y, dy);
				}																				
				break;																			
				case eappt_peng_0_4:																
				{																				
					cgpu_detail_mt::fcn_fxg_dfxg<eappt_peng_0_4, T>(g, Z, coef.cl, coef.cnl, y, dy);
				}																				
				break;																			
				case eappt_peng_0_12:															
				{																				
					cgpu_detail_mt::fcn_fxg_dfxg<eappt_peng_0_12, T>(g, Z, coef.cl, coef.cnl, y, dy);
				}																				
				break;													
				case eappt_kirkland_0_12:														
				{																				
					cgpu_detail_mt::fcn_fxg_dfxg<eappt_kirkland_0_12, T>(g, Z, coef.cl, coef.cnl, y, dy);
				}																				
				break;	
				case eappt_weickenmeier_0_12:															
				{																				
					cgpu_detail_mt::fcn_fxg_dfxg<eappt_weickenmeier_0_12, T>(g, coef.cl, coef.cnl, y, dy);
				}																				
				break;																			
				case eappt_lobato_0_12:																
				{																				
					cgpu_detail_mt::fcn_fxg_dfxg<eappt_lobato_0_12, T>(g, coef.cl, coef.cnl, y, dy);
				}																				
				break;
				case eappt_peng_ion_0_4:															
				{																				
					cgpu_detail_mt::fcn_fxg_dfxg<eappt_peng_ion_0_4, T>(g, Z, coef.cl, coef.cnl, y, dy);
				}																				
				break;																			
			}
		}

		/*************************************** pr ********************************************/
		template <class T>
		CPU_EXEC
		T fcn_pr(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const T& r, const pLNL_Coef_cpu<T>& coef)
		{
			SWITCH_FCN_ATOMIC_V(atomic_pot_parm_typ, fcn_pr, r, coef.cl, coef.cnl);
		}

		template <class T>
		CPU_EXEC
		void fcn_pr_dpr(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const T& r, const pLNL_Coef_cpu<T>& coef, T& y, T& dy)
		{
			SWITCH_FCN_ATOMIC_VD(atomic_pot_parm_typ, fcn_pr_dpr, r, coef.cl, coef.cnl, y, dy);
		}

		/**************************************** vr *******************************************/
		template <class T>
		CPU_EXEC
		T fcn_vr(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const T& r, const pLNL_Coef_cpu<T>& coef)
		{
			SWITCH_FCN_ATOMIC_V(atomic_pot_parm_typ, fcn_vr, r, coef.cl, coef.cnl);
		}

		template <class T>
		CPU_EXEC
		void fcn_vr_dvr(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const T& r, const pLNL_Coef_cpu<T>& coef, T& y, T& dy)
		{
			SWITCH_FCN_ATOMIC_VD(atomic_pot_parm_typ, fcn_vr_dvr, r, coef.cl, coef.cnl, y, dy);
		}

		/************************************** vz *********************************************/
		template <class T>
		CPU_EXEC
		T fcn_vz(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const T& r, const T& z_0, const T& z_e, const pLNL_Coef_cpu<T>& coef, const pQuad_Coef_1d_cpu<T>& quad)
		{
			SWITCH_FCN_ATOMIC_V(atomic_pot_parm_typ, fcn_vz, r, z_0, z_e, coef.cl, coef.cnl, quad.m_size, quad.x, quad.w);
		}

		template <class T>
		CPU_EXEC
		void fcn_vz_dvz(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const T& r, const T& z_0, const T& z_e, const pLNL_Coef_cpu<T>& coef, const pQuad_Coef_1d_cpu<T>& quad, T& y, T& dy)
		{
			SWITCH_FCN_ATOMIC_VD(atomic_pot_parm_typ, fcn_vz_dvz, r, z_0, z_e, coef.cl, coef.cnl, quad.m_size, quad.x, quad.w, y, dy);
		}

		/************************************* vzp *********************************************/
		template <class T>
		CPU_EXEC
		T fcn_vzp(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const T& r, const pLNL_Coef_cpu<T>& coef, pQuad_Coef_1d_cpu<T>* pquad = nullptr)
		{
			switch(atomic_pot_parm_typ)															
			{																				
				case eappt_doyle_0_4:														
				{																			
					return cgpu_detail_mt::fcn_vzp<eappt_doyle_0_4, T>(r, coef.cl, coef.cnl);
				}																			
				case eappt_peng_0_4:															
				{																			
					return cgpu_detail_mt::fcn_vzp<eappt_peng_0_4, T>(r, coef.cl, coef.cnl);
				}																			
				case eappt_peng_0_12:														
				{																			
					return cgpu_detail_mt::fcn_vzp<eappt_peng_0_12, T>(r, coef.cl, coef.cnl);
				}																			
				case eappt_kirkland_0_12:													
				{																			
					return cgpu_detail_mt::fcn_vzp<eappt_kirkland_0_12, T>(r, coef.cl, coef.cnl);
				}
				case eappt_weickenmeier_0_12:														
				{																			
					return cgpu_detail_mt::fcn_vzp<eappt_weickenmeier_0_12, T>(r, coef.cl, coef.cnl, pquad->m_size, pquad->x, pquad->w);
				}																			
				case eappt_lobato_0_12:															
				{																			
					return cgpu_detail_mt::fcn_vzp<eappt_lobato_0_12, T>(r, coef.cl, coef.cnl);
				}
				case eappt_peng_ion_0_4:														
				{																			
					return cgpu_detail_mt::fcn_vzp<eappt_peng_ion_0_4, T>(r, coef.cl, coef.cnl);
				}																			
				default:																	
				{																			
					return 0;																
				}																			
			}
		}

		template <class T>
		CPU_EXEC
		void fcn_vzp_dvzp(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const T& r, const pLNL_Coef_cpu<T>& coef, T& y, T& dy, pQuad_Coef_1d_cpu<T>* pquad = nullptr)
		{
			switch(atomic_pot_parm_typ)																
			{																					
				case eappt_doyle_0_4:															
				{																				
					cgpu_detail_mt::fcn_vzp_dvzp<eappt_doyle_0_4, T>(r, coef.cl, coef.cnl, y, dy);
				}																				
				break;																			
				case eappt_peng_0_4:																
				{																				
					cgpu_detail_mt::fcn_vzp_dvzp<eappt_peng_0_4, T>(r, coef.cl, coef.cnl, y, dy);
				}																				
				break;																			
				case eappt_peng_0_12:															
				{																				
					cgpu_detail_mt::fcn_vzp_dvzp<eappt_peng_0_12, T>(r, coef.cl, coef.cnl, y, dy);
				}																				
				break;																			
				case eappt_kirkland_0_12:														
				{																				
					cgpu_detail_mt::fcn_vzp_dvzp<eappt_kirkland_0_12, T>(r, coef.cl, coef.cnl, y, dy);
				}																				
				break;	
				case eappt_weickenmeier_0_12:															
				{																				
					cgpu_detail_mt::fcn_vzp_dvzp<eappt_weickenmeier_0_12, T>(r, coef.cl, coef.cnl, pquad->m_size, pquad->x, pquad->w, y, dy);
				}																				
				break;																			
				case eappt_lobato_0_12:																
				{																				
					cgpu_detail_mt::fcn_vzp_dvzp<eappt_lobato_0_12, T>(r, coef.cl, coef.cnl, y, dy);
				}																				
				break;	
				case eappt_peng_ion_0_4:	
				{																				
					cgpu_detail_mt::fcn_vzp_dvzp<eappt_peng_ion_0_4, T>(r, coef.cl, coef.cnl, y, dy);
				}																				
				break;																			
			}
		}
	}

	/* detector integration */
	namespace mt
	{
		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type<TVctr>>
		fcn_int_det_ring(Grid_2d<T>& grid, T g_min, T g_max, TVctr& mx_i, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			T g2_min = ::square(g_min);
			T g2_max = ::square(g_max);

			return fcn_stream_exec_xd_krn_reduce<edim_2>(pstream, grid.nx, grid.ny, std::plus<KS<T>>(), T(0), cgpu_detail_mt::fcn_int_det_ring<T, U>, grid, g2_min, g2_max, mx_i.m_data);
		}

		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type_r<TVctr>>
		fcn_int_det_ring_norm_2(Grid_2d<T>& grid, T g_min, T g_max, TVctr& mx_i, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;
			using Ur = Value_type_r<TVctr>;

			T g2_min = ::square(g_min);
			T g2_max = ::square(g_max);

			return fcn_stream_exec_xd_krn_reduce<edim_2>(pstream, grid.nx, grid.ny, std::plus<KS<Ur>>(), Ur(0), cgpu_detail_mt::fcn_int_det_ring_norm_2<T, U>, grid, g2_min, g2_max, mx_i.m_data);
		}

		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, Value_type<TVctr>>
		fcn_int_det_sen(Grid_2d<T>& grid, TVctr& sen_i, TVctr& mx_i, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			KS<U> sum_total = U(0);

			return fcn_stream_exec_xd_krn_reduce<edim_2>(pstream, grid.nx, grid.ny, std::plus<KS<U>>(), U(0), cgpu_detail_mt::fcn_int_det_sen<T, U>, grid, sen_i.m_data, mx_i.m_data);
		}

		template <class T, class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, Value_type_r<TVctr_2>>
		fcn_int_det_sen_norm_2(Grid_2d<T>& grid, TVctr_1& sen_i, TVctr_2& mx_i, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_2>;
			using Ur = Value_type_r<TVctr_2>;

			return fcn_stream_exec_xd_krn_reduce<edim_2>(pstream, grid.nx, grid.ny, std::plus<KS<Ur>>(), Ur(0), cgpu_detail_mt::fcn_int_det_sen_norm_2<T, U>, grid, sen_i.m_data, mx_i.m_data);
		}
	}

	/* wave propagation */
	namespace mt
	{
		/* propagate */
		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_fs_propagate(Grid_2d<T>& grid, TVctr& psi_i, R_2d<T> g_0, T w_g2, T w, TVctr& psi_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail_mt::fcn_fs_propagate<T, U>, grid, psi_i.m_data, g_0, w_g2, w, psi_o.m_data);
		}

		/* propagate and bandwith limit using a fermi aperture */
		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_fs_propagate_bw_f(Grid_2d<T>& grid, TVctr& psi_i, R_2d<T> g_0, T w_g2, T g2_cut, T alpha, T w, TVctr& psi_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail_mt::fcn_fs_propagate_bw_f<T, U>, grid, psi_i.m_data, g_0, w_g2, g2_cut, alpha, w, psi_o.m_data);
		}

		/* propagate and bandwith limit using a hard aperture */
		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_fs_propagate_bw_h(Grid_2d<T>& grid, TVctr& psi_i, R_2d<T> g_0, T w_g2, T g2_cut, T w, TVctr& psi_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail_mt::fcn_fs_propagate_bw_h<T, U>, grid, psi_i.m_data, g_0, w_g2, g2_cut, w, psi_o.m_data);
		}
	}

	/* probe - ctf - pctf */
	namespace mt
	{
		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_fs_probe(Grid_2d<T>& grid, Lens<T>& lens, R_2d<T> r, R_2d<T> gu, TVctr& psi_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			auto bb_phi_eq = lens.is_phi_required();

			if (bb_phi_eq)
			{
				fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, fcn_fs_probe<true, T, U>, grid, lens, r, gu, psi_o.m_data);
			}
			else
			{
				fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, fcn_fs_probe<false, T, U>, grid, lens, r, gu, psi_o.m_data);
			}

			auto tot_intensity = fcn_sum_norm_2(psi_o, pstream);
			auto sc_fctr = sqrt(T(1)/tot_intensity);

			fcn_scale(sc_fctr, psi_o, pstream);
		}

		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_fs_apply_ctf(Grid_2d<T>& grid, TVctr& psi_i, Lens<T>& lens, R_2d<T> gu, TVctr& psi_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			auto bb_phi_eq = lens.is_phi_required();

			if (bb_phi_eq)
			{
				fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, fcn_fs_apply_ctf<true, T, U>, grid, psi_i, lens, gu, psi_o.m_data);
			}
			else
			{
				fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, fcn_fs_apply_ctf<false, T, U>, grid, psi_i, lens, gu, psi_o.m_data);
			}
		}

		/* apply partial coherent transfer function */
		// Marc De Graef - Introduction to Conventional Transmission Electron Microscopy page: 611
		// tp_inc_iehwgd: 08, spt_inc_theta_c: 611
		template <class T, class TVctr>
		enable_if_vctr_cpu<TVctr, void>
		fcn_fs_apply_pctf(Grid_2d<T>& grid, TVctr& psi_i, Lens<T>& lens, R_2d<T> gu, TVctr& psi_o, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			auto bb_phi_eq = lens.is_phi_required();

			if (bb_phi_eq)
			{
				fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, fcn_fs_apply_pctf<true, T, U>, grid, psi_i, lens, gu, psi_o.m_data);
			}
			else
			{
				fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, fcn_fs_apply_pctf<false, T, U>, grid, psi_i, lens, gu, psi_o.m_data);
			}
		}
	}

	/* transmission function */
	namespace mt
	{
		template <class T, class TVctr_1, class TVctr_2>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_1, TVctr_2, void>
		transmission_fcn(eElec_Spec_Int_Mod esim, TVctr_1& vzp_i, T w, TVctr_2& tfcn_o, Stream_cpu* pstream = nullptr)
		{	
			using U = Value_type<TVctr_2>;

			iGrid_1d igrid(vzp_i.size());

			if (esim==eesim_weak_phase_object)
			{
				fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, cgpu_detail_mt::fcn_trans_fcn<eesim_weak_phase_object, U>, igrid, vzp_i.m_data, w, tfcn_o.m_data);
			}
			else
			{
				fcn_stream_exec_xd_krn<edim_2>(pstream, igrid.nx, cgpu_detail_mt::fcn_trans_fcn<eesim_phase_object, U>, igrid, vzp_i.m_data, w, tfcn_o.m_data);
			}
		}
	}

	/* eels */
	namespace mt
	{
		template <class T>
		T fcn_eels_lorentz_norm_factor_cpu(Grid_2d<T>& grid, EELS<T>& eels, Stream_cpu* pstream = nullptr)
		{
			auto sum_total = fcn_stream_exec_xd_krn_reduce<edim_2>(pstream, grid.nx, grid.ny, std::plus<KS<T>>(), T(0),cgpu_detail_mt::fcn_eels_lorentz_norm_factor<T>, grid, eels.gc2, eels.ge2);

			return sqrt(eels.occ)/sum_total;
		}

		template <class T, class TVctr_c>
		enable_if_vctr_cpu<TVctr_c, void>
		fcn_eels_w_xyz(Grid_2d<T>& grid, EELS<T>& eels, FFT_cpu<T>& fft, TVctr_c& w_x, TVctr_c& w_y, TVctr_c& w_z, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			eels.factor = fcn_eels_lorentz_norm_factor_cpu(grid, eels, pstream);

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail_mt::fcn_eels_w_xyz<T, U>, grid, eels, w_x.m_data, w_y.m_data, w_z.m_data);

			fft.inverse(w_x);
			fft.inverse(w_y);
			fft.inverse(w_z);
		}

		template <class T, class TVctr_c>
		enable_if_vctr_cpu<TVctr_c, void>
		fcn_eels_w_x(Grid_2d<T>& grid, EELS<T>& eels, FFT_cpu<T>& fft, TVctr_c& w_x, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			eels.factor = fcn_eels_lorentz_norm_factor_cpu(grid, eels, pstream);

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail_mt::fcn_eels_w_x<T, U>, grid, eels, w_x.m_data);

			fft.inverse(w_x);
		}

		template <class T, class TVctr_c>
		enable_if_vctr_cpu<TVctr_c, void>
		fcn_eels_w_y(Grid_2d<T>& grid, EELS<T>& eels, FFT_cpu<T>& fft, TVctr_c& w_y, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			eels.factor = fcn_eels_lorentz_norm_factor_cpu(grid, eels, pstream);

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail_mt::fcn_eels_w_y<T, U>, grid, eels, w_y.m_data);

			fft.inverse(w_y);
		}

		template <class T, class TVctr_c>
		enable_if_vctr_cpu<TVctr_c, void>
		fcn_eels_w_z(Grid_2d<T>& grid, EELS<T>& eels, FFT_cpu<T>& fft, TVctr_c& w_z, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			eels.factor = fcn_eels_lorentz_norm_factor_cpu(grid, eels, pstream);

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail_mt::fcn_eels_w_z<T, U>, grid, eels, w_z.m_data);

			fft.inverse(w_z);
		}

		template <class T, class TVctr_c>
		enable_if_vctr_cpu<TVctr_c, void>
		fcn_eels_w_mn1(Grid_2d<T>& grid, EELS<T>& eels, FFT_cpu<T>& fft, TVctr_c& w_mn1, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			eels.factor = fcn_eels_lorentz_norm_factor_cpu(grid, eels, pstream);

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail_mt::fcn_eels_w_mn1<T, U>, grid, eels, w_mn1.m_data);

			fft.inverse(w_mn1);
		}

		template <class T, class TVctr_c>
		enable_if_vctr_cpu<TVctr_c, void>
		fcn_eels_w_mp1(Grid_2d<T>& grid, EELS<T>& eels, FFT_cpu<T>& fft, TVctr_c& w_mp1, Stream_cpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			eels.factor = fcn_eels_lorentz_norm_factor_cpu(grid, eels, pstream);

			fcn_stream_exec_xd_krn<edim_2>(pstream, grid.nx, grid.ny, cgpu_detail_mt::fcn_eels_w_mp1<T, U>, grid, eels, w_mp1.m_data);

			fft.inverse(w_mp1);
		}
	}

	namespace mt
	{
		// /***********************Temporal and Spatial quadratures*********************/
		// template <class T>
		// void obj_lens_temporal_spatial_quadratures(Lens<T>& lens, Quad_Coef_1d<T, edev_cpu>& qt, Quad_Coef_2d<T, edev_cpu>& qs)
		// {
		// 	/*********************Temporal quad_data**********************/
		// 	Quad_Data quad_data;
		// 	quad_data(7, lens.tp_inc_npts, qt);		// 7: int_-infty^infty f(x) x^0 Exp[-x^2] dx
		// 	std::for_each(qt.w.begin(), qt.w.end(), [](T& v) { v = v/c_pii2; });

		// 	/*********************Spatial quad_data**********************/
		// 	qs.reserve((2 * lens.ngxs + 1)*(2 * lens.ngys + 1));
		// 	dt_int32 iqs = 0;
		// 	T sum_w = 0;
		// 	T sum_ee = 0;
		// 	T alpha = 0.5/pow(lens.spt_inc_sigma, 2);

		// 	for(auto ix = -lens.ngxs; ix <= lens.ngxs; ix++)
		// 	{
		// 		for(auto iy = -lens.ngys; iy <= lens.ngys; iy++)
		// 		{
		// 			T gxs = lens.gxs(ix);
		// 			T gys = lens.gys(iy);
		// 			T g2s = gxs*gxs + gys*gys;
		// 			if (g2s < lens.g2_outers)
		// 			{
		// 				qs.x.push_back(gxs);
		// 				qs.y.push_back(gys);
		// 				T v = exp(-alpha*g2s);
		// 				qs.w.push_back(v);
		// 				fcn_kh_sum(sum_w, v, sum_ee);
		// 			}
		// 		}
		// 	}
		// 	qs.resize(iqs);
		// 	std::for_each(qs.w.begin(), qs.w.end(), [sum_w](T& v) { v = v/sum_w; });
		// }
	
		// template <class T>
		// void cond_lens_temporal_spatial_quadratures(Lens<T>& lens, Quad_Coef_1d<dt_float64, edev_cpu>& qt, Quad_Coef_2d<dt_float64, edev_cpu>& qs)
		// {
		// 	/*********************Temporal quad_data**********************/
		// 	dt_bool bb_ti = (lens.tp_inc_npts > 1) && fcn_is_nzero(lens.tp_inc_sigma);
		// 	Quad_Data quad_data;

		// 	if (bb_ti)
		// 	{
		// 		dt_float64 a = (fcn_is_zero(lens.tp_inc_sigma))?0:(lens.tp_inc_a/(sqrt(c_2pi)*lens.tp_inc_sigma));
		// 		dt_float64 b = (fcn_is_zero(lens.tp_inc_sigma))?0:(0.5/lens.tp_inc_sigma_2());
		// 		dt_float64 c = (fcn_is_zero(lens.tp_inc_beta))?0:((1-lens.tp_inc_a)/(2*lens.tp_inc_beta));
		// 		dt_float64 d = (fcn_is_zero(lens.tp_inc_beta))?0:(1.0/lens.tp_inc_beta);

		// 		dt_float64 b_q = 0.5/sqrt(lens.tp_inc_a*lens.tp_inc_sigma_2() + 2*(1-lens.tp_inc_a)*lens.tp_inc_beta_2());

		// 		// 26, generalized hermite, (-inf, inf) |x-a|^alpha*exp(-b*(x-a)^2)
		// 		quad_data(26, lens.tp_inc_npts, qt, 0, 0, 0, b_q);

		// 		for(auto it = 0; it < lens.tp_inc_npts; it++)
		// 		{
		// 			dt_float64 r = abs(qt.x[it]);		// positive
		// 			dt_float64 f = a*exp(-b*r*r + b_q*r*r) + c*exp(-d*r + b_q*r*r);
		// 			qt.w[it] *= f;
		// 		}
		// 	}
		// 	else
		// 	{
		// 		qt.resize(1);
		// 		qt.x[0] = 0;
		// 		qt.w[0] = 1;
		// 	}

		// 	/*********************Spatial quad_data**********************/
		// 	dt_bool bb_si = ((lens.spt_inc_rad_npts > 1) || (lens.spt_inc_azm_npts > 1)) && !fcn_is_zero(lens.spt_inc_sigma, lens.spt_inc_beta);

		// 	if (bb_si)
		// 	{
		// 		dt_float64 a = (fcn_is_zero(lens.spt_inc_sigma))?0:(lens.spt_inc_a/(c_2pi<T>*lens.spt_inc_sigma_2()));
		// 		dt_float64 b = (fcn_is_zero(lens.spt_inc_sigma))?0:(0.5/lens.spt_inc_sigma_2());
		// 		dt_float64 c = (fcn_is_zero(lens.spt_inc_beta))?0:((1-lens.spt_inc_a)/(c_2pi<T>*lens.spt_inc_beta_2()));
		// 		dt_float64 d = (fcn_is_zero(lens.spt_inc_beta))?0:(1.0/lens.spt_inc_beta);

		// 		dt_float64 b_q = 1.0/sqrt((2*lens.spt_inc_a*lens.spt_inc_sigma_2()+ 6*(1-lens.spt_inc_a)*lens.spt_inc_beta_2())/6);

		// 		// radial part
		// 		Quad_Coef_1d<dt_float64, edev_cpu> qr;
		// 		// 25, generalized laguerre, (a, inf) (x-a)^alpha*exp(-b*(x-a))
		// 		quad_data(25, lens.spt_inc_rad_npts, qr, 1, 0, 0, b_q);

		// 		for(auto ir = 0; ir < lens.spt_inc_rad_npts; ir++)
		// 		{
		// 			dt_float64 r = abs(qr.x[ir]);		// positive
		// 			dt_float64 f = a*exp(-b*r*r + b_q*r) + c*exp(-d*r + b_q*r);
		// 			qr.w[ir] *= f;
		// 		}

		// 		// Azimuth part
		// 		Quad_Coef_1d<dt_float64, edev_cpu> qa;
		// 		qa.resize(lens.spt_inc_azm_npts);
		// 		dt_float64 h = c_2pi/lens.spt_inc_azm_npts;
		// 		for(auto ia = 0; ia < lens.spt_inc_azm_npts; ia++)
		// 		{
		// 			qa.x[ia] = ia*h;
		// 			qa.w[ia] = h;
		// 		}

		// 		qs.reserve(lens.spt_inc_rad_npts*lens.spt_inc_azm_npts);
		// 		for(auto ir = 0; ir < lens.spt_inc_rad_npts; ir++)
		// 		{
		// 			for(auto ia = 0; ia < lens.spt_inc_azm_npts; ia++)
		// 			{
		// 				dt_float64 sin_theta, cos_theta;
		// 				sincos(qa.x[ia], &sin_theta, &cos_theta);
		// 				qs.x.push_back(qr.x[ir]*cos_theta);
		// 				qs.y.push_back(qr.x[ir]*sin_theta);
		// 				qs.w.push_back(qr.w[ir]*qa.w[ia]);
		// 			}
		// 		}
		// 	}
		// 	else
		// 	{
		// 		qs.resize(1);
		// 		qs.x[0] = 0;
		// 		qs.y[0] = 0;
		// 		qs.w[0] = 1;
		// 	}
		// }
	}

#endif