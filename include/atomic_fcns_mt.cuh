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

#ifndef ATOMIC_FCNS_MT_H
	#define ATOMIC_FCNS_MT_H

	#ifdef _MSC_VER
		#pragma once
	#endif

	#include "const_enum_mt.cuh"
	#include "math_mt.h"
	#include "type_traits_mt.cuh"
	#include "intrpl_coef.cuh"
	#include "quad_data.cuh"
	#include "cgpu_vctr.cuh"
	#include "atomic_data_mt.cuh"
	#include "detail_cgpu_mt.cuh"
	#include "fcns_mt_cpu.h"

	#ifdef __CUDACC__
		#include "gpu_fcns_mt.cuh"
	#endif

	namespace mt
	{
		#define SWITCH_FCN_EVAL_FCN_COEF_LNL_VCTR(atomic_pot_parm_typ, pfcn, ...)													\
		switch(atomic_pot_parm_typ)																								\
		{																													\
			case eappt_doyle_0_4:																							\
			{																												\
				fcn_eval_fcn_coef_lnl<T>(__VA_ARGS__, detail_cgpu_mt::pfcn<eappt_doyle_0_4, T, Dev>());						\
			}																												\
			break;																											\
			case eappt_peng_0_4:																								\
			{																												\
				fcn_eval_fcn_coef_lnl<T>(__VA_ARGS__, detail_cgpu_mt::pfcn<eappt_peng_0_4, T, Dev>());						\
			}																												\
			break;																											\
			case eappt_peng_0_12:																							\
			{																												\
				fcn_eval_fcn_coef_lnl<T>(__VA_ARGS__, detail_cgpu_mt::pfcn<eappt_peng_0_12, T, Dev>());						\
			}																												\
			break;																											\
			case eappt_kirkland_0_12:																						\
			{																												\
				fcn_eval_fcn_coef_lnl<T>(__VA_ARGS__, detail_cgpu_mt::pfcn<eappt_kirkland_0_12, T, Dev>());					\
			}																												\
			break;																											\
			case eappt_weickenmeier_0_12:																					\
			{																												\
				fcn_eval_fcn_coef_lnl<T>(__VA_ARGS__, detail_cgpu_mt::pfcn<eappt_weickenmeier_0_12, T, Dev>());				\
			}																												\
			break;																											\
			case eappt_lobato_0_12:																							\
			{																												\
				fcn_eval_fcn_coef_lnl<T>(__VA_ARGS__, detail_cgpu_mt::pfcn<eappt_lobato_0_12, T, Dev>());					\
			}																												\
			break;																											\
			case eappt_peng_ion_0_4:																							\
			{																												\
				fcn_eval_fcn_coef_lnl<T>(__VA_ARGS__, detail_cgpu_mt::pfcn<eappt_peng_ion_0_4, T, Dev>());					\
			}																												\
			break;																											\
		}

		/************************************** feg ********************************************/
		template <class T, eDev Dev>
		void fcn_feg(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const pVctr_32<T, Dev>& g, const pLNL_Coef<T, Dev>& coef, pVctr_32<T, Dev>& y)
		{
			SWITCH_FCN_EVAL_FCN_COEF_LNL_VCTR(atomic_pot_parm_typ, pFcn_feg, g, coef, y);
		}

		template <class T, eDev Dev>
		void fcn_feg_dfeg(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const pVctr_32<T, Dev>& g, const pLNL_Coef<T, Dev>& coef, pVctr_32<T, Dev>& y, pVctr_32<T, Dev>& dy)
		{
			SWITCH_FCN_EVAL_FCN_COEF_LNL_VCTR(atomic_pot_parm_typ, pFcn_feg_dfeg, g, coef, y, dy);
		}

		/*************************************** fxg *******************************************/
		template <class T, eDev Dev>
		void fcn_fxg(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const pVctr_32<T, Dev>& g, const T& Z, const pLNL_Coef<T, Dev>& coef, pVctr_32<T, Dev>& y)
		{
			switch(atomic_pot_parm_typ)																					
			{																										
				case eappt_doyle_0_4:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(g, Z, coef, y, detail_cgpu_mt::pFcn_fxg<eappt_doyle_0_4, T, Dev>());		
				}																									
				break;																								
				case eappt_peng_0_4:																					
				{																									
					fcn_eval_fcn_coef_lnl<T>(g, Z, coef, y, detail_cgpu_mt::pFcn_fxg<eappt_peng_0_4, T, Dev>());		
				}																									
				break;																								
				case eappt_peng_0_12:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(g, Z, coef, y, detail_cgpu_mt::pFcn_fxg<eappt_peng_0_12, T, Dev>());		
				}																									
				break;																								
				case eappt_kirkland_0_12:																			
				{																									
					fcn_eval_fcn_coef_lnl<T>(g, Z, coef, y, detail_cgpu_mt::pFcn_fxg<eappt_kirkland_0_12, T, Dev>());	
				}																																																		
				break;	
				case eappt_weickenmeier_0_12:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(g, coef, y, detail_cgpu_mt::pFcn_fxg<eappt_weickenmeier_0_12, T, Dev>());		
				}																									
				break;																								
				case eappt_lobato_0_12:																					
				{																									
					fcn_eval_fcn_coef_lnl<T>(g, coef, y, detail_cgpu_mt::pFcn_fxg<eappt_lobato_0_12, T, Dev>());		
				}																									
				break;	
				case eappt_peng_ion_0_4:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(g, Z, coef, y, detail_cgpu_mt::pFcn_fxg<eappt_peng_ion_0_4, T, Dev>());	
				}																									
				break;																								
			}
		}

		template <class T, eDev Dev>
		void fcn_fxg_dfxg(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const pVctr_32<T, Dev>& g, const T& Z, const pLNL_Coef<T, Dev>& coef, pVctr_32<T, Dev>& y, pVctr_32<T, Dev>& dy)
		{
			switch(atomic_pot_parm_typ)																					
			{																										
				case eappt_doyle_0_4:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(g, Z, coef, y, dy, detail_cgpu_mt::pFcn_fxg_dfxg<eappt_doyle_0_4, T, Dev>());		
				}																									
				break;																								
				case eappt_peng_0_4:																					
				{																									
					fcn_eval_fcn_coef_lnl<T>(g, Z, coef, y, dy, detail_cgpu_mt::pFcn_fxg_dfxg<eappt_peng_0_4, T, Dev>());		
				}																									
				break;																								
				case eappt_peng_0_12:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(g, Z, coef, y, dy, detail_cgpu_mt::pFcn_fxg_dfxg<eappt_peng_0_12, T, Dev>());		
				}																									
				break;																								
				case eappt_kirkland_0_12:																			
				{																									
					fcn_eval_fcn_coef_lnl<T>(g, Z, coef, y, dy, detail_cgpu_mt::pFcn_fxg_dfxg<eappt_kirkland_0_12, T, Dev>());	
				}																																																		
				break;	
				case eappt_weickenmeier_0_12:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(g, coef, y, dy, detail_cgpu_mt::pFcn_fxg_dfxg<eappt_weickenmeier_0_12, T, Dev>());		
				}																									
				break;																								
				case eappt_lobato_0_12:																					
				{																									
					fcn_eval_fcn_coef_lnl<T>(g, coef, y, dy, detail_cgpu_mt::pFcn_fxg_dfxg<eappt_lobato_0_12, T, Dev>());		
				}
				break;	
				case eappt_peng_ion_0_4:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(g, Z, coef, y, dy, detail_cgpu_mt::pFcn_fxg_dfxg<eappt_peng_ion_0_4, T, Dev>());	
				}																									
				break;																								
			}
		}

		/*************************************** pr ********************************************/
		template <class T, eDev Dev>
		void fcn_pr(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const pVctr_32<T, Dev>& r, const pLNL_Coef<T, Dev>& coef, pVctr_32<T, Dev>& y)
		{
			SWITCH_FCN_EVAL_FCN_COEF_LNL_VCTR(atomic_pot_parm_typ, pFcn_pr, r, coef, y);
		}

		template <class T, eDev Dev>
		void fcn_pr_dpr(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const pVctr_32<T, Dev>& r, const pLNL_Coef<T, Dev>& coef, pVctr_32<T, Dev>& y, pVctr_32<T, Dev>& dy)
		{
			SWITCH_FCN_EVAL_FCN_COEF_LNL_VCTR(atomic_pot_parm_typ, pFcn_pr_dpr, r, coef, y, dy);
		}

		/**************************************** vr *******************************************/
		template <class T, eDev Dev>
		void fcn_vr(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const pVctr_32<T, Dev>& r, const pLNL_Coef<T, Dev>& coef, pVctr_32<T, Dev>& y)
		{
			SWITCH_FCN_EVAL_FCN_COEF_LNL_VCTR(atomic_pot_parm_typ, pFcn_vr, r, coef, y);
		}

		template <class T, eDev Dev>
		void fcn_vr_dvr(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const pVctr_32<T, Dev>& r, const pLNL_Coef<T, Dev>& coef, pVctr_32<T, Dev>& y, pVctr_32<T, Dev>& dy)
		{
			SWITCH_FCN_EVAL_FCN_COEF_LNL_VCTR(atomic_pot_parm_typ, pFcn_vr_dvr, r, coef, y, dy);
		}

		/************************************** vz *********************************************/
		template <class T, eDev Dev>
		void fcn_vz(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const pVctr_32<T, Dev>& r, const T& z_0, const T& z_e, const pLNL_Coef<T, Dev>& coef, const pQuad_Coef_1d<T, Dev>& quad, pVctr_32<T, Dev>& y)
		{
			SWITCH_FCN_EVAL_FCN_COEF_LNL_VCTR(atomic_pot_parm_typ, pFcn_vz, r, z_0, z_e, coef, quad, y);
		}

		template <class T, eDev Dev>
		void fcn_vz_dvz(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const pVctr_32<T, Dev>& r, const T& z_0, const T& z_e, const pLNL_Coef<T, Dev>& coef, const pQuad_Coef_1d<T, Dev>& quad, pVctr_32<T, Dev>& y, pVctr_32<T, Dev>& dy)
		{
			SWITCH_FCN_EVAL_FCN_COEF_LNL_VCTR(atomic_pot_parm_typ, pFcn_vz_dvz, r, z_0, z_e, coef, quad, y, dy);
		}

		/************************************* vzp *********************************************/
		template <class T, eDev Dev>
		void fcn_vzp(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const pVctr_32<T, Dev>& r, const pLNL_Coef<T, Dev>& coef, pVctr_32<T, Dev>& y, pQuad_Coef_1d<T, Dev>* pquad = nullptr)
		{
			switch(atomic_pot_parm_typ)																					
			{																										
				case eappt_doyle_0_4:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(r, coef, y, detail_cgpu_mt::pFcn_vzp<eappt_doyle_0_4, T, Dev>());		
				}																									
				break;																								
				case eappt_peng_0_4:																					
				{																									
					fcn_eval_fcn_coef_lnl<T>(r, coef, y, detail_cgpu_mt::pFcn_vzp<eappt_peng_0_4, T, Dev>());		
				}																									
				break;																								
				case eappt_peng_0_12:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(r, coef, y, detail_cgpu_mt::pFcn_vzp<eappt_peng_0_12, T, Dev>());		
				}																									
				break;																								
				case eappt_kirkland_0_12:																			
				{																									
					fcn_eval_fcn_coef_lnl<T>(r, coef, y, detail_cgpu_mt::pFcn_vzp<eappt_kirkland_0_12, T, Dev>());	
				}																																																		
				break;	
				case eappt_weickenmeier_0_12:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(r, coef, *pquad, y, detail_cgpu_mt::pFcn_vzp<eappt_weickenmeier_0_12, T, Dev>());		
				}																									
				break;																								
				case eappt_lobato_0_12:																					
				{																									
					fcn_eval_fcn_coef_lnl<T>(r, coef, y, detail_cgpu_mt::pFcn_vzp<eappt_lobato_0_12, T, Dev>());		
				}																									
				break;	
				case eappt_peng_ion_0_4:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(r, coef, y, detail_cgpu_mt::pFcn_vzp<eappt_peng_ion_0_4, T, Dev>());	
				}																									
				break;																								
			}
		}

		template <class T, eDev Dev>
		void fcn_vzp_dvzp(const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const pVctr_32<T, Dev>& r, const pLNL_Coef<T, Dev>& coef, pVctr_32<T, Dev>& y, pVctr_32<T, Dev>& dy, pQuad_Coef_1d<T, Dev>* pquad = nullptr)
		{
			switch(atomic_pot_parm_typ)																					
			{																										
				case eappt_doyle_0_4:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(r, coef, y, dy, detail_cgpu_mt::pFcn_vzp_dvzp<eappt_doyle_0_4, T, Dev>());
				}																									
				break;																								
				case eappt_peng_0_4:																					
				{																									
					fcn_eval_fcn_coef_lnl<T>(r, coef, y, dy, detail_cgpu_mt::pFcn_vzp_dvzp<eappt_peng_0_4, T, Dev>());	
				}																									
				break;																								
				case eappt_peng_0_12:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(r, coef, y, dy, detail_cgpu_mt::pFcn_vzp_dvzp<eappt_peng_0_12, T, Dev>());
				}																									
				break;																								
				case eappt_kirkland_0_12:																			
				{																									
					fcn_eval_fcn_coef_lnl<T>(r, coef, y, dy, detail_cgpu_mt::pFcn_vzp_dvzp<eappt_kirkland_0_12, T, Dev>());
				}
				break;		
				case eappt_weickenmeier_0_12:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(r, coef, *pquad, y, dy, detail_cgpu_mt::pFcn_vzp_dvzp<eappt_weickenmeier_0_12, T, Dev>());
				}																									
				break;																								
				case eappt_lobato_0_12:																					
				{																									
					fcn_eval_fcn_coef_lnl<T>(r, coef, y, dy, detail_cgpu_mt::pFcn_vzp_dvzp<eappt_lobato_0_12, T, Dev>());		
				}																									
				break;	
				case eappt_peng_ion_0_4:																				
				{																									
					fcn_eval_fcn_coef_lnl<T>(r, coef, y, dy, detail_cgpu_mt::pFcn_vzp_dvzp<eappt_peng_ion_0_4, T, Dev>());	
				}																									
				break;																								
			}
		}

		/***************************************************************************************/
		/********************************** atomic fcns ****************************************/
		/***************************************************************************************/
		template <class T>
		class Atomic_Fcns
		{
		public:
			using value_type = T;

			Atomic_Fcns(): pcoef(nullptr)
			{
				Quad_Data quad_data;
				quad_data(eqt_tanh_sinh_int_n1_p1, c_nqz, quad_a_b);				// 1: eqt_tanh_sinh_int_n1_p1 -> int_-1^1 f(x) dx
				quad_data(eqt_exp_sinh_int_0_pinfty, c_nqz, quad_0_infty);			// 2: eqt_exp_sinh_int_0_pinfty -> int_0^infty f(x) dx

				pquad_a_b = quad_a_b;
				pquad_0_infty = quad_0_infty;
			}

			Atomic_Fcns(Atomic_Coef_cpu<T>& coef): Atomic_Fcns()
			{
				set_atomic_coef(coef);
			}

			void set_atomic_coef(Atomic_Coef_cpu<T>& coef)
			{
				pcoef = &coef;
			}

			/***************************************************************************************/
			T feg(const T& g)
			{
				return fcn_feg<T>(pcoef->atomic_pot_parm_typ, g, pcoef->feg);
			}

			void feg_dfeg(const T& g, T& y, T& dy)
			{
				fcn_feg_dfeg<T>(pcoef->atomic_pot_parm_typ, g, pcoef->feg, y, dy);
			}

			void feg(const pVctr_cpu_32<T>&& g, pVctr_cpu_32<T>&& y)
			{
				fcn_feg<T, edev_cpu>(pcoef->atomic_pot_parm_typ, g, pcoef->feg, y);
			}

			void feg_dfeg(const pVctr_cpu_32<T>&& g, pVctr_cpu_32<T>&& y, pVctr_cpu_32<T>&& dy)
			{
				fcn_feg_dfeg<T, edev_cpu>(pcoef->atomic_pot_parm_typ, g, pcoef->feg, y, dy);
			}

			/***************************************************************************************/
			T fxg(const T& g)
			{
				return fcn_fxg<T>(pcoef->atomic_pot_parm_typ, g, pcoef->Z_diff(), pcoef->fxg);
			}

			void fxg_dfxg(const T& g, T& y, T& dy)
			{
				fcn_fxg_dfxg<T>(pcoef->atomic_pot_parm_typ, g, pcoef->Z_diff(), pcoef->fxg, y, dy);
			}

			void fxg(const pVctr_cpu_32<T>&& g, pVctr_cpu_32<T>&& y)
			{
				fcn_fxg<T, edev_cpu>(pcoef->atomic_pot_parm_typ, g, pcoef->Z_diff(), pcoef->fxg, y);
			}

			void fxg_dfxg(const pVctr_cpu_32<T>&& g, pVctr_cpu_32<T>&& y, pVctr_cpu_32<T>&& dy)
			{
				fcn_fxg_dfxg<T, edev_cpu>(pcoef->atomic_pot_parm_typ, g, pcoef->Z_diff(), pcoef->fxg, y, dy);
			}

			/***************************************************************************************/
			T pr(const T& r)
			{
				return fcn_pr<T>(pcoef->atomic_pot_parm_typ, r, pcoef->pr);
			}

			void pr_dpr(const T& r, T& y, T& dy)
			{
				fcn_pr_dpr<T>(pcoef->atomic_pot_parm_typ, r, pcoef->pr, y, dy);
			}

			void pr(const pVctr_cpu_32<T>&& r, pVctr_cpu_32<T>&& y)
			{
				fcn_pr<T, edev_cpu>(pcoef->atomic_pot_parm_typ, r, pcoef->pr, y);
			}

			void pr_dpr(const pVctr_cpu_32<T>&& r, pVctr_cpu_32<T>&& y, pVctr_cpu_32<T>&& dy)
			{
				fcn_pr_dpr<T, edev_cpu>(pcoef->atomic_pot_parm_typ, r, pcoef->pr, y, dy);
			}
			
			/***************************************************************************************/
			T vr(const T& r)
			{
				return fcn_vr<T>(pcoef->atomic_pot_parm_typ, r, pcoef->vr);
			}

			void vr_dvr(const T& r, T& y, T& dy)
			{
				fcn_vr_dvr<T>(pcoef->atomic_pot_parm_typ, r, pcoef->vr, y, dy);
			}

			void vr(const pVctr_cpu_32<T>&& r, pVctr_cpu_32<T>&& y)
			{
				fcn_vr<T, edev_cpu>(pcoef->atomic_pot_parm_typ, r, pcoef->vr, y);
			}

			void vr_dvr(const pVctr_cpu_32<T>&& r, pVctr_cpu_32<T>&& y, pVctr_cpu_32<T>&& dy)
			{
				fcn_vr_dvr<T, edev_cpu>(pcoef->atomic_pot_parm_typ, r, pcoef->vr, y, dy);
			}

			/***************************************************************************************/
			T vz(const T& z_0, const T& z_e, const T& r)
			{
				return fcn_vz<T>(pcoef->atomic_pot_parm_typ, r, z_0, z_e, pcoef->vr, pquad_a_b);
			}

			void vz_dvz(const T& z_0, const T& z_e, const T& r, T& y, T& dy)
			{
				fcn_vz_dvz<T>(pcoef->atomic_pot_parm_typ, r, z_0, z_e, pcoef->vr, pquad_a_b, y, dy);
			}

			void vz(const T& z_0, const T& z_e, const pVctr_cpu_32<T>&& r, pVctr_cpu_32<T>&& y)
			{
				fcn_vz<T, edev_cpu>(pcoef->atomic_pot_parm_typ, r, z_0, z_e, pcoef->vr, pquad_a_b, y);
			}

			void vz_dvz(const T& z_0, const T& z_e, const pVctr_cpu_32<T>&& r, pVctr_cpu_32<T>&& y, pVctr_cpu_32<T>&& dy)
			{
				fcn_vz_dvz<T, edev_cpu>(pcoef->atomic_pot_parm_typ, r, z_0, z_e, pcoef->vr, pquad_a_b, y, dy);
			}

			/***************************************************************************************/
			T vzp(const T& r)
			{
				return fcn_vzp<T>(pcoef->atomic_pot_parm_typ, r, pcoef->vzp, &pquad_0_infty);
			}

			void vzp_dvzp(const T& r, T& y, T& dy)
			{
				fcn_vzp_dvzp<T>(pcoef->atomic_pot_parm_typ, r, pcoef->vzp, y, dy, &pquad_0_infty);
			}

			void vzp(const pVctr_cpu_32<T>&& r, pVctr_cpu_32<T>&& y)
			{
				fcn_vzp<T, edev_cpu>(pcoef->atomic_pot_parm_typ, r, pcoef->vzp, y, &pquad_0_infty);
			}

			void vzp_dvzp(const pVctr_cpu_32<T>&& r, pVctr_cpu_32<T>&& y, pVctr_cpu_32<T>&& dy)
			{
				fcn_vzp_dvzp<T, edev_cpu>(pcoef->atomic_pot_parm_typ, r, pcoef->vzp, y, dy, &pquad_0_infty);
			}

			/***************************************************************************************/
			T atomic_radius_rms(const dt_int32& dim)
			{
				if (fcn_is_zero(pcoef->vr.cl[0]))
				{
					return T(0);
				}

				auto fcn = [&](const T& r) { return (dim == 3)?vr(r):vzp(r); };

				KS<T> vr_sum = T(0);
				KS<T> vr2_sum = T(0);
				for(auto ik = 0; ik < quad_0_infty.size(); ik++)
				{
					const auto r_ik = quad_0_infty.x[ik];
					const auto vr_ik = quad_0_infty.w[ik]*fcn(r_ik)*pow(r_ik, dim-1);
					vr_sum += vr_ik;
					vr2_sum += vr_ik*r_ik*r_ik;
				}

				return ::sqrt(vr2_sum/vr_sum);
			}

			T atomic_radius_cutoff(const dt_int32& dim, const T& vr_lim)
			{
				if (fcn_is_zero(pcoef->vr.cl[0]) || fabs(vr_lim)<1e-8)
				{
					return T(0);
				}

				auto fcn = [&](const T& r) { return ((dim == 3)?vr(r):vzp(r))-vr_lim; };

				const T r_min = 1e-5;
				const T r_max = 25.0;

				return fcn_fd_root(r_min, r_max, 1e-8, 200, fcn);
			}

		private:
			Atomic_Coef_cpu<T>* pcoef;

			Quad_Coef_1d_cpu<T> quad_a_b;
			Quad_Coef_1d_cpu<T> quad_0_infty;
			pQuad_Coef_1d_cpu<T> pquad_a_b;
			pQuad_Coef_1d_cpu<T> pquad_0_infty;
		};

	}

#endif