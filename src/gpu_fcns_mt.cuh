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

#ifndef GPU_FCNS_MT_H
	#define GPU_FCNS_MT_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include "const_enum_mt.cuh"
	#include "math.cuh"
	#include "types.cuh"
	#include "type_traits_mt.cuh"
	#include "cgpu_vctr.cuh"
	#include "cgpu_fft.cuh"
	#include "cgpu_stream.cuh"
	#include "gpu_detail_mt.cuh"

	/* atomic functions */
	namespace mt
	{
		/*********************************** pFcn_clnl_3_x<T> **********************************/
		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_gpu_32<T>& x, const pLNL_Coef_gpu<T>& coef, pVctr_gpu_32<T>& fx, pFcn_clnl_3_1<T> fcn)
		{
			gpu_detail_mt::fcn_eval_fcn_coef_lnl<<<x.d_grid_size(), x.d_blk_size()>>>(x, coef, fx, fcn);
		}

		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_gpu_32<T>& x, const pLNL_Coef_gpu<T>& coef, pVctr_gpu_32<T>& fx1, pVctr_gpu_32<T>& fx2, pFcn_clnl_3_2<T> fcn)
		{
			gpu_detail_mt::fcn_eval_fcn_coef_lnl<<<x.d_grid_size(), x.d_blk_size()>>>(x, coef, fx1, fx2, fcn);
		}

		/********************************* pFcn_clnl_4_x<T> ************************************/
		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_gpu_32<T>& x, const T& c, const pLNL_Coef_gpu<T>& coef, pVctr_gpu_32<T>& fx, pFcn_clnl_4_1<T> fcn)
		{
			gpu_detail_mt::fcn_eval_fcn_coef_lnl<<<x.d_grid_size(), x.d_blk_size()>>>(x, c, coef, fx, fcn);
		}

		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_gpu_32<T>& x, const T& c, const pLNL_Coef_gpu<T>& coef, pVctr_gpu_32<T>& fx1, pVctr_gpu_32<T>& fx2, pFcn_clnl_4_2<T> fcn)
		{
			gpu_detail_mt::fcn_eval_fcn_coef_lnl<<<x.d_grid_size(), x.d_blk_size()>>>(x, c, coef, fx1, fx2, fcn);
		}
		
		/********************************* pFcn_clnl_6_x<T> ************************************/
		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_gpu_32<T>& x, const pLNL_Coef_gpu<T>& coef, const pQuad_Coef_1d_gpu<T>& quad, pVctr_gpu_32<T>& fx, pFcn_clnl_6_1<T> fcn)
		{
			gpu_detail_mt::fcn_eval_fcn_coef_lnl<<<x.d_grid_size(), x.d_blk_size()>>>(x, coef, quad, fx, fcn);
		}

		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_gpu_32<T>& x, const pLNL_Coef_gpu<T>& coef, const pQuad_Coef_1d_gpu<T>& quad, pVctr_gpu_32<T>& fx1, pVctr_gpu_32<T>& fx2, pFcn_clnl_6_2<T> fcn)
		{
			gpu_detail_mt::fcn_eval_fcn_coef_lnl<<<x.d_grid_size(), x.d_blk_size()>>>(x, coef, quad, fx1, fx2, fcn);
		}

		/********************************* pFcn_clnl_8_x<T> ************************************/
		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_gpu_32<T>& x, const T& z_0, const T& z_e, const pLNL_Coef_gpu<T>& coef, const pQuad_Coef_1d_gpu<T>& quad, pVctr_gpu_32<T>& fx, pFcn_clnl_8_1<T> fcn)
		{
			gpu_detail_mt::fcn_eval_fcn_coef_lnl<<<x.d_grid_size(), x.d_blk_size()>>>(x, z_0, z_e, coef, quad, fx, fcn);
		}

		template <class T>
		void fcn_eval_fcn_coef_lnl(const pVctr_gpu_32<T>& x, const T& z_0, const T& z_e, const pLNL_Coef_gpu<T>& coef, const pQuad_Coef_1d_gpu<T>& quad, pVctr_gpu_32<T>& fx1, pVctr_gpu_32<T>& fx2, pFcn_clnl_8_2<T> fcn)
		{
			gpu_detail_mt::fcn_eval_fcn_coef_lnl<<<x.d_grid_size(), x.d_blk_size()>>>(x, z_0, z_e, coef, quad, fx1, fx2, fcn);
		}
	}

	/* detector integration */
	namespace mt
	{
		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type<TVctr>>
		fcn_int_det_ring(Grid_2d<T>& grid, T g_min, T g_max, TVctr& mx_i, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			T g2_min = ::square(g_min);
			T g2_max = ::square(g_max);

			auto dgb = grid.d_grid_blk(dim3(c_thr_2d_x, c_thr_2d_y));
			TVctr sum_v(dgb.grid_size());

			gpu_detail_mt::fcn_int_det_ring<T, U><<<dgb.grid, dgb.blk, dgb.smems_red()*sizeof(U)>>>(grid, g2_min, g2_max, mx_i.ptr_32(), sum_v.ptr_32());

			return gpu_detail::fcn_sum_rem_cpu(sum_v);
		}

		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type_r<TVctr>>
		fcn_int_det_ring_norm_2(Grid_2d<T>& grid, T g_min, T g_max, TVctr& mx_i, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;
			using Ur = Value_type_r<TVctr>;

			T g2_min = ::square(g_min);
			T g2_max = ::square(g_max);

			auto dgb = grid.d_grid_blk(dim3(c_thr_2d_x, c_thr_2d_y));
			Vctr_gpu<Ur> sum_v(dgb.grid_size());

			gpu_detail_mt::fcn_int_det_ring_norm_2<T, U><<<dgb.grid, dgb.blk, dgb.smems_red()*sizeof(Ur)>>>(grid, g2_min, g2_max, mx_i.ptr_32(), sum_v.ptr_32());

			return gpu_detail::fcn_sum_rem_cpu(sum_v);
		}

		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, Value_type<TVctr>>
		fcn_int_det_sen(Grid_2d<T>& grid, TVctr& sen_i, TVctr& mx_i, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			auto dgb = grid.d_grid_blk(dim3(c_thr_2d_x, c_thr_2d_y));
			TVctr sum_v(dgb.grid_size());

			gpu_detail_mt::fcn_int_det_sen<T, U><<<dgb.grid, dgb.blk, dgb.smems_red()*sizeof(U)>>>(grid, sen_i.ptr_32(), mx_i.ptr_32(), sum_v.ptr_32());

			return gpu_detail::fcn_sum_rem_cpu(sum_v);
		}

		template <class T, class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, Value_type_r<TVctr_2>>
		fcn_int_det_sen_norm_2(Grid_2d<T>& grid, TVctr_1& sen_i, TVctr_2& mx_i, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_2>;
			using Ur = Value_type_r<TVctr_2>;

			auto dgb = grid.d_grid_blk(dim3(c_thr_2d_x, c_thr_2d_y));
			Vctr_gpu<Ur> sum_v(dgb.grid_size());

			gpu_detail_mt::fcn_int_det_sen_norm_2<T, U><<<dgb.grid, dgb.blk, dgb.smems_red()*sizeof(Ur)>>>(grid, sen_i.ptr_32(), mx_i.ptr_32(), sum_v.ptr_32());

			return gpu_detail::fcn_sum_rem_cpu(sum_v);
		}
	}

	/* wave propagation */
	namespace mt
	{
		/* propagate */
		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fs_propagate(Grid_2d<T>& grid, TVctr& psi_i, R_2d<T> g_0, T w_g2, T w, TVctr& psi_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			gpu_detail_mt::fcn_fs_propagate<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, psi_i.ptr_32(), g_0, w_g2, w, psi_o.ptr_32());
		}

		/* propagate and bandwith limit using a fermi aperture */
		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fs_propagate_bw_f(Grid_2d<T>& grid, TVctr& psi_i, R_2d<T> g_0, T w_g2, T g2_cut, T alpha, T w, TVctr& psi_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			gpu_detail_mt::fcn_fs_propagate_bw_f<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, psi_i.ptr_32(), g_0, w_g2, g2_cut, alpha, w, psi_o.ptr_32());
		}

		/* propagate and bandwith limit using a hard aperture */
		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fs_propagate_bw_h(Grid_2d<T>& grid, TVctr& psi_i, R_2d<T> g_0, T w_g2, T g2_cut, T w, TVctr& psi_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			gpu_detail_mt::fcn_fs_propagate_bw_h<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, psi_i.ptr_32(), g_0, w_g2, g2_cut, w, psi_o.ptr_32());
		}
	}
	
	/* probe - ctf - pctf */
	namespace mt
	{
		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fs_probe(Grid_2d<T>& grid, Lens<T>& lens, R_2d<T> r, R_2d<T> gu, TVctr& psi_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			auto bb_phi_eq = lens.is_phi_required();

			if (bb_phi_eq)
			{
				gpu_detail_mt::fcn_fs_probe<true, T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, lens, r, gu, psi_o.ptr_32());
			}
			else
			{
				gpu_detail_mt::fcn_fs_probe<false, T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, lens, r, gu, psi_o.ptr_32());
			}

			auto tot_intensity = fcn_sum_norm_2(psi_o, pstream);
			auto sc_fctr = ::sqrt(T(1)/tot_intensity);

			fcn_scale(sc_fctr, psi_o, pstream);
		}

		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fs_apply_ctf(Grid_2d<T>& grid, TVctr& psi_i, Lens<T>& lens, R_2d<T> gu, TVctr& psi_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			auto bb_phi_eq = lens.is_phi_required();

			if (bb_phi_eq)
			{
				gpu_detail_mt::fcn_fs_apply_ctf<true, T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, psi_i, lens, gu, psi_o.ptr_32());
			}
			else
			{
				gpu_detail_mt::fcn_fs_apply_ctf<false, T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, psi_i, lens, gu, psi_o.ptr_32());
			}
		}

		/* apply partial coherent transfer function */
		// Marc De Graef - Introduction to Conventional Transmission Electron Microscopy page: 611
		// tp_inc_iehwgd: 08, spt_inc_theta_c: 611
		template <class T, class TVctr>
		enable_if_vctr_gpu<TVctr, void>
		fcn_fs_apply_pctf(Grid_2d<T>& grid, TVctr& psi_i, Lens<T>& lens, TVctr& psi_o, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr>;

			auto bb_phi_eq = lens.is_phi_required();

			if (bb_phi_eq)
			{
				gpu_detail_mt::fcn_fs_apply_pctf<true, T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, psi_i, lens, gu, psi_o.ptr_32());
			}
			else
			{
				gpu_detail_mt::fcn_fs_apply_pctf<false, T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, psi_i, lens, gu, psi_o.ptr_32());
			}
		}
	}

	/* transmission function */
	namespace mt
	{
		template <class T, class TVctr_1, class TVctr_2>
		enable_if_vctr_gpu_and_vctr_gpu<TVctr_1, TVctr_2, void>
		transmission_fcn(eElec_Spec_Int_Model esim, TVctr_1& vzp_i, T w, TVctr_2& tfcn_o, Stream_gpu* pstream = nullptr)
		{	
			using U = Value_type<TVctr_2>;

			iGrid_1d igrid(vzp_i.size());

			if (esim==eesim_weak_phase_object)
			{
				gpu_detail_mt::fcn_trans_fcn<eesim_weak_phase_object, U><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, vzp_i.ptr_32(), w, tfcn_o.ptr_32());
			}
			else
			{
				gpu_detail_mt::fcn_trans_fcn<eesim_phase_object, U><<<igrid.d_grid(), igrid.d_blk()>>>(igrid, vzp_i.ptr_32(), w, tfcn_o.ptr_32());
			}
		}
	}
	
	/* eels */
	namespace mt
	{
		template <class T>
		T fcn_eels_lorentz_norm_factor_gpu(Grid_2d<T>& grid, EELS<T>& eels, Stream_gpu* pstream = nullptr)
		{
			auto dgb = grid.d_grid_blk(dim3(c_thr_2d_x, c_thr_2d_y));
			Vctr_gpu<T> sum_v(dgb.grid_size());

			gpu_detail_mt::fcn_eels_lorentz_norm_factor<T><<<dgb.grid, dgb.blk, dgb.smems_red()*sizeof(T)>>>(grid, eels.gc2, eels.ge2, sum_v.ptr_32());

			return ::sqrt(eels.occ)/gpu_detail::fcn_sum_rem_cpu(sum_v);
		}

		template <class T, class TVctr_c>
		enable_if_vctr_gpu<TVctr_c, void>
		fcn_eels_w_xyz(Grid_2d<T>& grid, EELS<T>& eels, FFT_gpu<T>& fft, TVctr_c& w_x, TVctr_c& w_y, TVctr_c& w_z, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			eels.factor = fcn_eels_lorentz_norm_factor_gpu(grid, eels, pstream);

			gpu_detail_mt::fcn_eels_w_xyz<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, eels, w_x.ptr_32(), w_y.ptr_32(), w_z.ptr_32());

			fft.inverse(w_x);
			fft.inverse(w_y);
			fft.inverse(w_z);
		}

		template <class T, class TVctr_c>
		enable_if_vctr_gpu<TVctr_c, void>
		fcn_eels_w_x(Grid_2d<T>& grid, EELS<T>& eels, FFT_gpu<T>& fft, TVctr_c& w_x, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			eels.factor = fcn_eels_lorentz_norm_factor_gpu(grid, eels, pstream);

			gpu_detail_mt::fcn_eels_w_x<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, eels, w_x.ptr_32());

			fft.inverse(w_x);
		}

		template <class T, class TVctr_c>
		enable_if_vctr_gpu<TVctr_c, void>
		fcn_eels_w_y(Grid_2d<T>& grid, EELS<T>& eels, FFT_gpu<T>& fft, TVctr_c& w_y, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			eels.factor = fcn_eels_lorentz_norm_factor_gpu(grid, eels, pstream);

			gpu_detail_mt::fcn_eels_w_y<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, eels, w_y.ptr_32());

			fft.inverse(w_y);
		}

		template <class T, class TVctr_c>
		enable_if_vctr_gpu<TVctr_c, void>
		fcn_eels_w_z(Grid_2d<T>& grid, EELS<T>& eels, FFT_gpu<T>& fft, TVctr_c& w_z, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			eels.factor = fcn_eels_lorentz_norm_factor_gpu(grid, eels, pstream);

			gpu_detail_mt::fcn_eels_w_z<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, eels, w_z.ptr_32());

			fft.inverse(w_z);
		}

		template <class T, class TVctr_c>
		enable_if_vctr_gpu<TVctr_c, void>
		fcn_eels_w_mn1(Grid_2d<T>& grid, EELS<T>& eels, FFT_gpu<T>& fft, TVctr_c& w_mn1, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			eels.factor = fcn_eels_lorentz_norm_factor_gpu(grid, eels, pstream);

			gpu_detail_mt::fcn_eels_w_mn1<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, eels, w_mn1.ptr_32());

			fft.inverse(w_mn1);
		}

		template <class T, class TVctr_c>
		enable_if_vctr_gpu<TVctr_c, void>
		fcn_eels_w_mp1(Grid_2d<T>& grid, EELS<T>& eels, FFT_gpu<T>& fft, TVctr_c& w_mp1, Stream_gpu* pstream = nullptr)
		{
			using U = Value_type<TVctr_c>;

			eels.factor = fcn_eels_lorentz_norm_factor_gpu(grid, eels, pstream);

			gpu_detail_mt::fcn_eels_w_mp1<T, U><<<grid.d_grid(), grid.d_blk()>>>(grid, eels, w_mp1.ptr_32());

			fft.inverse(w_mp1);
		}
	}

#endif