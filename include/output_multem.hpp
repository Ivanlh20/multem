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

#ifndef OUTPUT_MULTEM_H
#define OUTPUT_MULTEM_H

#include <vector>
#include <algorithm>

#include "math_mt.h"
#include "types.cuh"
#include "type_traits_gen.h"
#include "cgpu_stream.cuh"
#include "particles.cuh"
#include "in_classes.cuh"
#include "fcns_cpu.h"
#include "fcns_gpu.h"

namespace mt
{
	inline
		dt_bool is_ot_image_tot_coh(const eEM_Output_Typ &output_type)
	{
		return output_type == eemot_image_tot_coh;
	}

	inline
		dt_bool is_ot_image_tot(const eEM_Output_Typ &output_type)
	{
		return output_type == eemot_image_tot;
	}

	inline
		dt_bool is_ot_m2psi_tot_coh(const eEM_Output_Typ &output_type)
	{
		return output_type == eemot_m2psi_tot_coh;
	}

	inline
		dt_bool is_ot_m2psi_tot(const eEM_Output_Typ &output_type)
	{
		return output_type == eemot_m2psi_tot;
	}

	inline
		dt_bool is_ot_m2psi_tot_psi_coh(const eEM_Output_Typ &output_type)
	{
		return output_type == eemot_m2psi_tot_psi_coh;
	}

	inline
		dt_bool is_ot_psi_coh(const eEM_Output_Typ &output_type)
	{
		return output_type == eemot_psi_coh;
	}

	inline
		dt_bool is_ot_psi_0(const eEM_Output_Typ &output_type)
	{
		return output_type == eemot_psi_0;
	}

	inline
		dt_bool is_ot_V(const eEM_Output_Typ &output_type)
	{
		return output_type == eemot_V;
	}

	inline
		dt_bool is_ot_trans(const eEM_Output_Typ &output_type)
	{
		return output_type == eemot_trans;
	}

	/***************************************************************************************/
	template <class T>
	class Output_Multem:public Multem_In_Parm<T>
	{
	public:
		using T_r = T;
		using T_c = complex<T>;

		using TVctr_hr = host_vector<T>;
		using TVctr_hc = host_vector<complex<T>>;

		using TVctr_dr = device_vector<T>;
		using TVctr_dc = device_vector<complex<T>>;

		Output_Multem(): Multem_In_Parm<T>(), output_type(eemot_m2psi_tot), 
		ndetector(0), nx(0), ny(0), dx(0), dy(0), dr(0), n_thk(0), n_thk_d(0), bb_only_cpu_mem(false) {}

		Output_Multem(Multem_In_Parm<T> *multem_in_parm, dt_bool bb_only_cpu_mem_i=false)
		{
 set_in_data(multem_in_parm, bb_only_cpu_mem_i);
		}

		Output_Multem(const Output_Multem<T>& output_multem)
		{
 assign(output_multem);
		}

		template <class TOutput_Multem>
		void assign(TOutput_Multem &output_multem)
		{
			assign_multem_in_parm(output_multem);

			output_type = output_multem.output_type;
			ndetector = output_multem.ndetector;
			nx = output_multem.nx;
			ny = output_multem.ny;
			dx = output_multem.dx;
			dy = output_multem.dy;
			dr = output_multem.dr;

			x = output_multem.x;
			y = output_multem.y;
			r = output_multem.r;

			n_thk = output_multem.n_thk;

			bb_only_cpu_mem = output_multem.bb_only_cpu_mem;

			image_tot = output_multem.image_tot;

			image_coh = output_multem.image_coh;

			m2psi_tot = output_multem.m2psi_tot;

			m2psi_coh = output_multem.m2psi_coh;

			psi_coh = output_multem.psi_coh;

			V = output_multem.V;

			trans = output_multem.trans;

			psi_0 = output_multem.psi_0;
		}
 
		template <class TOutput_Multem>
		Output_Multem<T>& operator=(const TOutput_Multem &output_multem)
		{
			assign(output_multem);
			return *this;
		}

		template <class TIn_Multislice>
		void set_in_data(TIn_Multislice *multem_in_parm, dt_bool bb_only_cpu_mem_i=false)
		{
			clear();

			bb_only_cpu_mem = bb_only_cpu_mem_i;
			stream.resize(1);

			assign_multem_in_parm(*multem_in_parm);

			// set required number of thickness
			n_thk = this->thick.size();

			// check selected device
			auto bb_is_device = this->system_config.is_gpu() && !bb_only_cpu_mem;

			// get available gpu free memory
			dt_float64 free_gpu_memory_mb = (bb_is_device)?fcn_free_mem<edev_gpu>()-10:0;

			dt_int32 nxy_out = this->output_area.size_s();
			dt_int32 nxy_grid = this->grid_2d.size();

			if (this->system_config.is_gpu())
			{
				psi_zh.resize(nxy_grid);
				m2psi_zh.resize(nxy_grid);
			}

			ndetector = (this->is_STEM_ISTEM_EELS())?1:this->detector.size();
			dt_int32 n_beams = this->number_of_beams();

			set_output_grid();

			set_output_type();

 			auto n_vgpu_a_tr = [&](dt_int32 n_vgpu, dt_int32 nxy)->dt_int32
			{
 				dt_int32 n_vgpu_a = (bb_is_device)?dt_int32(::floor(free_gpu_memory_mb/mt::fcn_size_mb<T_r>(nxy))):0;
				return min(n_vgpu_a, n_vgpu);
			};

			auto n_vgpu_a_tc = [&](dt_int32 n_vgpu, dt_int32 nxy)->dt_int32
			{
 				dt_int32 n_vgpu_a = (bb_is_device)?dt_int32(::floor(free_gpu_memory_mb/mt::fcn_size_mb<T_c>(nxy))):0;
				return min(n_vgpu_a, n_vgpu);
			};

			switch (output_type)
			{
				case eemot_image_tot_coh:
				{
					image_tot.set_size(n_thk, ndetector, nxy_out);
					image_coh.set_size(n_thk, ndetector, nxy_out);

					psi_coh.set_size(n_thk, 1, nxy_grid);

					dt_int32 n_vgpu = n_thk*1;
					dt_int32 n_vgpu_a = n_vgpu_a_tc(n_vgpu, nxy_grid);

					psi_coh_d.set_size(n_thk, 1, nxy_grid, n_vgpu_a);
				}
				break;
				case eemot_image_tot:
				{
					image_tot.set_size(n_thk, ndetector, nxy_out);
				}
				break;
				case eemot_m2psi_tot_coh:
				{
					m2psi_tot.set_size(n_thk, n_beams, nxy_out);
					m2psi_coh.set_size(n_thk, n_beams, nxy_out);
					psi_coh.set_size(n_thk, n_beams, nxy_out);

					dt_int32 n_vgpu = n_thk*n_beams;
					dt_int32 n_vgpu_a = n_vgpu_a_tc(n_vgpu, 2*nxy_out);

					m2psi_tot_d.set_size(n_thk, n_beams, nxy_out, n_vgpu_a);
					m2psi_coh_d.set_size(n_thk, n_beams, nxy_out, n_vgpu_a);
					psi_coh_d.set_size(n_thk, n_beams, nxy_out, n_vgpu_a);
				}
				break;
				case eemot_m2psi_tot:
				{
					m2psi_tot.set_size(n_thk, n_beams, nxy_out);

					dt_int32 n_vgpu = n_thk*n_beams;
					dt_int32 n_vgpu_a = n_vgpu_a_tr(n_vgpu, nxy_out);

					m2psi_tot_d.set_size(n_thk, n_beams, nxy_out, n_vgpu_a);
				}
				break;
				case eemot_m2psi_tot_psi_coh:
				{
					m2psi_tot.set_size(n_thk, n_beams, nxy_out);
					psi_coh.set_size(n_thk, n_beams, nxy_out);

					dt_int32 n_vgpu = n_thk*n_beams;
					dt_int32 n_vgpu_a = n_vgpu_a_tr(n_vgpu, 3*nxy_out);

					m2psi_tot_d.set_size(n_thk, n_beams, nxy_out, n_vgpu_a);
					psi_coh_d.set_size(n_thk, n_beams, nxy_out, n_vgpu_a);
				}
				break;
				case eemot_psi_coh:
				{
					psi_coh.set_size(n_thk, n_beams, nxy_out);

					dt_int32 n_vgpu = n_thk*n_beams;
					dt_int32 n_vgpu_a = n_vgpu_a_tc(n_vgpu, nxy_out);

					psi_coh_d.set_size(n_thk, n_beams, nxy_out, n_vgpu_a);
				}
				break;
				case eemot_V:
				{
					V.set_size(n_thk, 1, nxy_out);
				}
				break;
				case eemot_trans:
				{
					trans.set_size(n_thk, 1, nxy_out);
				}
				break;
				case eemot_psi_0:
				{
					// check the number of beams of this
					psi_0.set_size(n_thk, 1, nxy_out);
				}
				break;
			}
		}

		void joint_data(std::vector<Output_Multem<T>>& output_multem_v)
		{
			dt_int32 n_out_v = output_multem_v.size();
			/***************************************************************************************/
			for(auto ik = 0; ik < image_tot.size(); ik++)
			{
				dt_int32 shift = 0;
 				for(auto iv = 0; iv < n_out_v; iv++)
				{
					auto &image_tot_v = output_multem_v[iv].image_tot[ik];
					thrust::copy(image_tot_v.begin(), image_tot_v.end(), image_tot[ik].begin()+shift);
					shift += image_tot_v.size();
				}
			}
			/***************************************************************************************/
			for(auto ik = 0; ik < image_coh.size(); ik++)
			{
				dt_int32 shift = 0;
 				for(auto iv = 0; iv < n_out_v; iv++)
				{
					auto &image_coh_v = output_multem_v[iv].image_coh[ik];
					thrust::copy(image_coh_v.begin(), image_coh_v.end(), image_coh[ik].begin()+shift);
					shift += image_coh_v.size();
				}
			}
			/***************************************************************************************/
			for(auto ik = 0; ik < m2psi_tot.size(); ik++)
			{
				dt_int32 shift = 0;
 				for(auto iv = 0; iv < n_out_v; iv++)
				{
					auto &m2psi_tot_v = output_multem_v[iv].m2psi_tot[ik];
					thrust::copy(m2psi_tot_v.begin(), m2psi_tot_v.end(), m2psi_tot[ik].begin()+shift);
					shift += m2psi_tot_v.size();
				}
			}
			/***************************************************************************************/
			for(auto ik = 0; ik < m2psi_coh.size(); ik++)
			{
				dt_int32 shift = 0;
 				for(auto iv = 0; iv < n_out_v; iv++)
				{
					auto &m2psi_coh_v = output_multem_v[iv].m2psi_coh[ik];
					thrust::copy(m2psi_coh_v.begin(), m2psi_coh_v.end(), m2psi_coh[ik].begin()+shift);
					shift += m2psi_coh_v.size();
				}
			}
			/***************************************************************************************/
			for(auto ik = 0; ik < psi_coh.size(); ik++)
			{
				dt_int32 shift = 0;
 				for(auto iv = 0; iv < n_out_v; iv++)
				{
					auto &psi_coh_v = output_multem_v[iv].psi_coh[ik];
					thrust::copy(psi_coh_v.begin(), psi_coh_v.end(), psi_coh[ik].begin()+shift);
					shift += psi_coh_v.size();
				}
			}
			/***************************************************************************************/
			for(auto ik = 0; ik < V.size(); ik++)
			{
				dt_int32 shift = 0;
 				for(auto iv = 0; iv < n_out_v; iv++)
				{
					auto &V_v = output_multem_v[iv].V[ik];
					thrust::copy(V_v.begin(), V_v.end(), V[ik].begin()+shift);
					shift += V_v.size();
				}
			}
			/***************************************************************************************/
			for(auto ik = 0; ik < trans.size(); ik++)
			{
				dt_int32 shift = 0;
 				for(auto iv = 0; iv < n_out_v; iv++)
				{
					auto &trans_v = output_multem_v[iv].trans[ik];
					thrust::copy(trans_v.begin(), trans_v.end(), trans[ik].begin()+shift);
					shift += trans_v.size();
				}
			}
			/***************************************************************************************/
			for(auto ik = 0; ik < psi_0.size(); ik++)
			{
				dt_int32 shift = 0;
 				for(auto iv = 0; iv < n_out_v; iv++)
				{
					auto &psi_0_v = output_multem_v[iv].psi_0[ik];
					thrust::copy(psi_0_v.begin(), psi_0_v.end(), psi_0[ik].begin()+shift);
					shift += psi_0_v.size();
				}
			}
		}

		void clear()
		{
			output_type = eemot_m2psi_tot;
			ndetector = 0;
			nx = 0;
			ny = 0;
			dx = 0;
			dy = 0;
			dr = 0;

			n_thk = 0;

			x.clear();
			x.shrink_to_fit();

			y.clear();
			y.shrink_to_fit();

			r.clear();
			r.shrink_to_fit();

			/***************************************************************************************/
			image_tot.clear_shrink_to_fit();
			image_coh.clear_shrink_to_fit();

			m2psi_tot.clear_shrink_to_fit();
			m2psi_coh.clear_shrink_to_fit();
			psi_coh.clear_shrink_to_fit();

			/***************************************************************************************/
			psi_zh.clear();
			psi_zh.shrink_to_fit();

			m2psi_zh.clear();
			m2psi_zh.shrink_to_fit();

 			m2psi_tot_d.clear_shrink_to_fit();
			m2psi_coh_d.clear_shrink_to_fit();
			psi_coh_d.clear_shrink_to_fit();

			/***************************************************************************************/
			V.clear_shrink_to_fit();
			trans.clear_shrink_to_fit();
			psi_0.clear_shrink_to_fit();
		}

		void clean_temporal()
		{
 			psi_zh.clear();
			psi_zh.shrink_to_fit();

			m2psi_zh.clear();
			m2psi_zh.shrink_to_fit();

			switch (output_type)
			{
				case eemot_image_tot_coh:
				{
					psi_coh.clear_shrink_to_fit();

					psi_coh_d.clear_shrink_to_fit();
				}
				break;
				case eemot_image_tot:
				{

				}
				break;
				case eemot_m2psi_tot_coh:
				{
 					psi_coh.clear_shrink_to_fit();

					m2psi_tot_d.clear_shrink_to_fit();
					m2psi_coh_d.clear_shrink_to_fit();
					psi_coh_d.clear_shrink_to_fit();
				}
				break;
				case eemot_m2psi_tot:
				{
					m2psi_tot_d.clear_shrink_to_fit();
				}
				break;
				case eemot_m2psi_tot_psi_coh:
				{
					m2psi_tot_d.clear_shrink_to_fit();
					psi_coh_d.clear_shrink_to_fit();
				}
				break;
				case eemot_psi_coh:
				{
					psi_coh_d.clear_shrink_to_fit();
				}
				break;
			}
		}

		void init()
		{
			image_tot.fill(0);
			image_coh.fill(0);

			m2psi_tot.fill(0);
			m2psi_coh.fill(0);
			psi_coh.fill(0);

			/***************************************************************************************/
			thrust::fill(psi_zh.begin(), psi_zh.end(), T_c(0));
			thrust::fill(m2psi_zh.begin(), m2psi_zh.end(), T_r(0));

			/***************************************************************************************/
 			m2psi_tot_d.fill(0);
			m2psi_coh_d.fill(0);
			psi_coh_d.fill(0);

			/***************************************************************************************/
			V.fill(0);
			trans.fill(0);
			psi_0.fill(0);
		}

		void init_psi_coh()
		{
			psi_coh.fill(0);
			psi_coh_d.fill(0);
		}

		/***************************************************************************************/
		template <class TVctr>
		void add_sc_psi_coh(dt_int32 irow, dt_int32 icol, T_c w, TVctr& phi)
		{
			if (psi_coh_d.exists(irow, icol))
			{
				thrust::transform(thrust::device, phi.begin(), phi.end(), psi_coh_d(irow, icol).begin(), psi_coh_d(irow, icol).begin(), cgpu_fctr::add_scale<T_c>(w));
			}
			else
			{
	 mt::add_sc_to_host(stream, w, phi, psi_coh(irow, icol), &psi_zh);
			}
		}

		template <class TVctr>
		void add_sc_sft_psi_coh(dt_int32 irow, dt_int32 icol, T_c w, TVctr& phi)
		{
			if (psi_coh_d.exists(irow, icol))
			{
				mt::fcn_add_sc_fftsft_2d(this->grid_2d, w, phi, psi_coh_d(irow, icol), &psi_zh);
			}
			else
			{
	 mt::fcn_add_sc_fftsft_2d(this->grid_2d, w, phi, psi_coh(irow, icol), &psi_zh);
			}
		}

 		template <class TVctr>
		void add_sc_crop_sft_psi_coh(dt_int32 irow, dt_int32 icol, T_c w, TVctr& phi)
		{
			if (psi_coh_d.exists(irow, icol))
			{
	 mt::fcn_add_sc_crop_fftsft_2d(this->grid_2d, w, phi, this->output_area, psi_coh_d(irow, icol), &psi_zh);
			}
			else
			{
	 mt::fcn_add_sc_crop_fftsft_2d(this->grid_2d, w, phi, this->output_area, psi_coh(irow, icol), &psi_zh);
			}
		}

 		template <class TVctr>
		void add_sc_crop_sft_m2psi_coh_from_psi(dt_int32 irow, dt_int32 icol, T_r w, TVctr& phi)
		{
			if (m2psi_coh_d.exists(irow, icol))
			{
	 mt::fcn_add_sc_norm_2_crop_fftsft_2d(this->grid_2d, w, phi, this->output_area, m2psi_coh_d(irow, icol), &psi_zh);
			}
			else
			{
	 mt::fcn_add_sc_norm_2_crop_fftsft_2d(this->grid_2d, w, phi, this->output_area, m2psi_coh(irow, icol), &psi_zh);
			}
		}

 		template <class TVctr>
		void add_sc_crop_sft_m2psi_coh_from_m2psi(dt_int32 irow, dt_int32 icol, T_r w, TVctr& m2phi)
		{
			if (m2psi_coh_d.exists(irow, icol))
			{
	 mt::fcn_add_sc_crop_fftsft_2d(this->grid_2d, w, m2phi, this->output_area, m2psi_coh_d(irow, icol), &m2psi_zh);
			}
			else
			{
	 mt::fcn_add_sc_crop_fftsft_2d(this->grid_2d, w, m2phi, this->output_area, m2psi_coh(irow, icol), &m2psi_zh);
			}
		}

 		template <class TVctr>
		void add_sc_crop_sft_m2psi_tot_from_psi(dt_int32 irow, dt_int32 icol, T_r w, TVctr phi)
		{
			if (m2psi_tot_d.exists(irow, icol))
			{
	 mt::fcn_add_sc_norm_2_crop_fftsft_2d(this->grid_2d, w, phi, this->output_area, m2psi_tot_d(irow, icol), &psi_zh);
			}
			else
			{
	 mt::fcn_add_sc_norm_2_crop_fftsft_2d(this->grid_2d, w, phi, this->output_area, m2psi_tot(irow, icol), &psi_zh);
			}
		}

 		template <class TVctr>
		void add_sc_crop_sft_m2psi_tot_from_m2psi(dt_int32 irow, dt_int32 icol, T_r w, TVctr& m2phi)
		{
			if (m2psi_tot_d.exists(irow, icol))
			{
	 mt::fcn_add_sc_crop_fftsft_2d(this->grid_2d, w, m2phi, this->output_area, m2psi_tot_d(irow, icol), &m2psi_zh);
			}
			else
			{
	 mt::fcn_add_sc_crop_fftsft_2d(this->grid_2d, w, m2phi, this->output_area, m2psi_tot(irow, icol), &m2psi_zh);
			}
		}

		/***************************************************************************************/
 		template <class TVctr>
		void set_psi_coh(dt_int32 irow, dt_int32 icol, TVctr& phi)
		{
			if (psi_coh_d.exists(irow, icol))
			{
				thrust::copy(phi.begin(), phi.end(), psi_coh_d(irow, icol).begin());
			}
			else
			{
				thrust::copy(phi.begin(), phi.end(), psi_coh(irow, icol).begin());
			}
		}

 		template <class TVctr>
		void set_sft_psi_coh(dt_int32 irow, dt_int32 icol, TVctr& phi)
		{
			if (psi_coh_d.exists(irow, icol))
			{
	 mt::fcn_fftsft_2d(this->grid_2d, phi, psi_coh_d(irow, icol), &psi_zh);
			}
			else
			{
	 mt::fcn_fftsft_2d(this->grid_2d, phi, psi_coh(irow, icol), &psi_zh);
			}
		}

 		template <class TVctr>
		void set_crop_sft_psi_coh(dt_int32 irow, dt_int32 icol, TVctr& phi)
		{
			if (psi_coh_d.exists(irow, icol))
			{
	 mt::fcn_assign_crop_fftsft_2d(this->grid_2d, phi, this->output_area, psi_coh_d(irow, icol), &psi_zh);
			}
			else
			{
	 mt::fcn_assign_crop_fftsft_2d(this->grid_2d, phi, this->output_area, psi_coh(irow, icol), &psi_zh);
			}
		}

 		template <class TVctr>
		void set_crop_sft_m2psi_coh(dt_int32 irow, dt_int32 icol, TVctr& m2phi)
		{
			if (m2psi_coh_d.exists(irow, icol))
			{
	 mt::fcn_assign_crop_fftsft_2d(this->grid_2d, m2phi, this->output_area, m2psi_coh_d(irow, icol), &m2psi_zh);
			}
			else
			{
	 mt::fcn_assign_crop_fftsft_2d(this->grid_2d, m2phi, this->output_area, m2psi_coh(irow, icol), &m2psi_zh);
			}
		}

 		template <class TVctr>
		void set_crop_sft_m2psi_tot(dt_int32 irow, dt_int32 icol, TVctr& m2phi)
		{
			if (m2psi_tot_d.exists(irow, icol))
			{
	 mt::fcn_assign_crop_fftsft_2d(this->grid_2d, m2phi, this->output_area, m2psi_tot_d(irow, icol), &m2psi_zh);
			}
			else
			{
	 mt::fcn_assign_crop_fftsft_2d(this->grid_2d, m2phi, this->output_area, m2psi_tot(irow, icol), &m2psi_zh);
			}
		}

 		/***************************************************************************************/
 		template <class TVctr>
		void from_psi_coh_2_phi(dt_int32 irow, dt_int32 icol, TVctr& phi)
		{
			if (psi_coh_d.exists(irow, icol))
			{
	 thrust::copy(psi_coh_d(irow, icol).begin(), psi_coh_d(irow, icol).end(), phi.begin());
			}
			else
			{
	 thrust::copy(psi_coh(irow, icol).begin(), psi_coh(irow, icol).end(), phi.begin());
			}
		}

 		template <class TVctr>
		void from_m2psi_coh_2_m2phi(dt_int32 irow, dt_int32 icol, TVctr& m2phi)
		{
			if (m2psi_coh_d.exists(irow, icol))
			{
	 thrust::copy(m2psi_coh_d(irow, icol).begin(), m2psi_coh_d(irow, icol).end(), m2phi.begin());
			}
			else
			{
	 thrust::copy(m2psi_coh(irow, icol).begin(), m2psi_coh(irow, icol).end(), m2phi.begin());
			}		
		}

 		template <class TVctr>
		void from_m2psi_tot_2_m2phi(dt_int32 irow, dt_int32 icol, TVctr& m2phi)
		{
			if (m2psi_tot_d.exists(irow, icol))
			{
				thrust::copy(m2psi_tot_d(irow, icol).begin(), m2psi_tot_d(irow, icol).end(), m2phi.begin());
			}
			else
			{
				thrust::copy(m2psi_tot(irow, icol).begin(), m2psi_tot(irow, icol).end(), m2phi.begin());
			}
		}

		void gather()
		{
			cpy_allow_data(m2psi_tot_d, m2psi_tot);
			cpy_allow_data(m2psi_coh_d, m2psi_coh);
			cpy_allow_data(psi_coh_d, psi_coh);
		}

		dt_int32 size() const { return nx*ny; }

		dt_int32 nxy_pot_g() const { return this->grid_2d.size(); }

		dt_int32 nxy_out_g() const { return this->output_area.size_s(); }

 		dt_int32 nthk() const { return n_thk; }

		/***************************************************************************************/
		inline
		dt_bool is_ot_image_tot_coh() const
		{
			return mt::is_ot_image_tot_coh(output_type);
		}

		inline
		dt_bool is_ot_image_tot() const
		{
			return mt::is_ot_image_tot(output_type);
		}

		inline
		dt_bool is_ot_m2psi_tot_coh() const
		{
			return mt::is_ot_m2psi_tot_coh(output_type);
		}

		inline
		dt_bool is_ot_m2psi_tot() const
		{
			return mt::is_ot_m2psi_tot(output_type);
		}

		inline
		dt_bool is_ot_m2psi_tot_psi_coh() const
		{
			return mt::is_ot_m2psi_tot_psi_coh(output_type);
		}

		inline
		dt_bool is_ot_psi_coh() const
		{
			return mt::is_ot_psi_coh(output_type);
		}

		inline
		dt_bool is_ot_psi_0() const
		{
			return mt::is_ot_psi_0(output_type);
		}

		inline
		dt_bool is_ot_V() const
		{
			return mt::is_ot_V(output_type);
		}

		inline
		dt_bool is_ot_trans() const
		{
			return mt::is_ot_trans(output_type);
		}

		host_vector<dt_float32> extract_data(eAtomic_Vib_Mod_Output fp_ctr, eShow_CData show_data, dt_int32 irow, dt_int32 icol = 0)
		{
			host_vector<dt_float32> data(size());

			switch (output_type)
			{
				case eemot_image_tot_coh:
				{
					data = (fp_ctr == epmo_total)?image_tot(irow, icol):image_coh(irow, icol);
				}
				break;
				case eemot_image_tot:
				{
					data = image_tot(irow, icol);
				}
				break;
				case eemot_m2psi_tot_coh:
				{
					data = (fp_ctr == epmo_total)?m2psi_tot(irow, icol):m2psi_coh(irow, icol);
				}
				break;
				case eemot_m2psi_tot:
				{
					data = m2psi_tot(irow, icol);
				}
				break;
				case eemot_m2psi_tot_psi_coh:
				{
					if (fp_ctr == epmo_total)
					{
						data = m2psi_tot(irow, icol);
					}
					else
					{
						from_complex_to_real(show_data, psi_coh(irow, icol), data);
					}
				}
				break;
				case eemot_psi_coh:
				{
					from_complex_to_real(show_data, psi_coh(irow, icol), data);
				}
				break;
				case eemot_V:
				{
					data = V(irow, icol);
				}
				break;
				case eemot_trans:
				{
					from_complex_to_real(show_data, trans(irow, icol), data);
				}
				break;
				case eemot_psi_0:
				{
					from_complex_to_real(show_data, psi_0(irow, icol), data);
				}
				break;
			}

			return data;
		}

		eEM_Output_Typ output_type;

		dt_int32 ndetector;
		dt_int32 nx;
		dt_int32 ny;
		T_r dx;
		T_r dy;
		T_r dr;

		host_vector<T_r> x;
		host_vector<T_r> y;
		host_vector<T_r> r;

		AV_2d<TVctr_hr> image_tot;
		AV_2d<TVctr_hr> image_coh;

		AV_2d<TVctr_hr> m2psi_tot;
		AV_2d<TVctr_hr> m2psi_coh;
		AV_2d<TVctr_hc> psi_coh;

		AV_2d<TVctr_hr> V;
		AV_2d<TVctr_hc> trans;
		AV_2d<TVctr_hc> psi_0;

		AV_2d<TVctr_hr> image_tot_d;
		AV_2d<TVctr_hr> image_coh_d;

		AV_2d<TVctr_dr> m2psi_tot_d;
		AV_2d<TVctr_dr> m2psi_coh_d;
		AV_2d<TVctr_dc> psi_coh_d;

		dt_bool bb_only_cpu_mem;

		Stream<edev_cpu> stream;
	private:
		Vctr<T_c, edev_cpu> psi_zh;
		Vctr<T_r, edev_cpu> m2psi_zh;

		dt_int32 n_thk;
		dt_int32 n_thk_d;

		template <class TAV_2d_1, class TAV_2d_2>
		void cpy_allow_data(TAV_2d_1 &av_2d_1, TAV_2d_2 &av_2d_2)
		{
			for(auto ik = 0; ik < av_2d_1.size_a(); ik++)
			{
				thrust::copy(av_2d_1[ik].begin(), av_2d_1[ik].end(), av_2d_2[ik].begin());
			}		
		}

		void set_output_grid()
		{
			if (this->is_STEM() || this->is_STEM_ISTEM_EELS())
			{
 				nx = this->scanning.nx;
				ny = this->scanning.ny;

				dx = this->scanning.dR.x;
				dy = this->scanning.dR.y;

				x = this->scanning.x();
				y = this->scanning.y();
				r = this->scanning.mR;
			}
			else
			{
				nx = this->output_area.nx_s();
				ny = this->output_area.ny_s();

				const dt_bool is_RS = this->is_grid_RS();
				dx = (is_RS)?this->grid_2d.drx:this->grid_2d.dgx;
				dy = (is_RS)?this->grid_2d.dry:this->grid_2d.dgy;

				x.resize(nx);
				y.resize(ny);

				for(auto ix = this->output_area.ix_0; ix < this->output_area.ix_e; ix++)
				{
					dt_int32 ix_s = ix-this->output_area.ix_0;
					x[ix_s] = (is_RS)?this->grid_2d.rx(ix):this->grid_2d.gx(ix);
				}

				for(auto iy = this->output_area.iy_0; iy < this->output_area.iy_e; iy++)
				{
					dt_int32 iy_s = iy-this->output_area.iy_0;
					y[iy_s] = (is_RS)?this->grid_2d.ry(iy):this->grid_2d.gy(iy);
				}
			}
		}

		// 1:(image_tot, image_coh);2:(image_tot);3:(m2psi_tot, m2psi_coh);4:(m2psi_tot);
		// 5:(m2psi_tot, psi_coh);6:(psi_coh);7:(psi_0);8:(V);9:(trans)
		void set_output_type()
		{
			if (this->is_STEM() || this->is_STEM_ISTEM_EELS())
			{
				output_type = (this->atomic_vib.coh_contrib)?eemot_image_tot_coh:eemot_image_tot;
			}
			else if (this->is_ISTEM() || this->is_CBED_CBEI() || this->is_PED_HCTEM() ||
				this->is_ED_HRTEM() || this->is_EFTEM())
			{
				output_type = (this->atomic_vib.coh_contrib)?eemot_m2psi_tot_coh:eemot_m2psi_tot;
			}
			else if (this->is_EWFS_EWRS())
			{
				output_type = (this->is_EWFS_EWRS_SC())?eemot_psi_coh:eemot_m2psi_tot_psi_coh;
			}
			else if (this->is_PropFS_PropRS())
			{
				output_type = eemot_psi_coh;
			}
			else if (this->is_IWFS_IWRS())
			{
				output_type = eemot_psi_0;
			}
			else if (this->is_PPFS_PPRS())
			{
				output_type = eemot_V;
			}
			else if (this->is_TFFS_TFRS())
			{
				output_type = eemot_trans;
			}
		}

		template <class TIn_Multislice>
		void assign_multem_in_parm(TIn_Multislice &multem_in_parm)
		{
			this->system_config = multem_in_parm.system_config;

			this->elec_spec_interact_mod = multem_in_parm.elec_spec_interact_mod;
			this->atomic_pot_parm_typ = multem_in_parm.atomic_pot_parm_typ;

			this->operation_mode = multem_in_parm.operation_mode;
			this->slice_storage = multem_in_parm.slice_storage;
			this->reverse_multislice = multem_in_parm.reverse_multislice;
			this->mul_sign = multem_in_parm.mul_sign;
			this->Vrl = multem_in_parm.Vrl;
			this->nR = multem_in_parm.nR;

			this->atomic_vib = multem_in_parm.atomic_vib;

			this->atoms = multem_in_parm.atoms;
			this->is_xtl = multem_in_parm.is_xtl;

			this->rot_in_parm = multem_in_parm.rot_in_parm;

			this->thick_type = multem_in_parm.thick_type;
			this->thick = multem_in_parm.thick;

			this->spec_slic_typ = multem_in_parm.spec_slic_typ;

			this->grid_2d = multem_in_parm.grid_2d;
			this->output_area = multem_in_parm.output_area;

			this->em_sim_typ = multem_in_parm.em_sim_typ;

			this->iw_type = multem_in_parm.iw_type;
			this->iw_psi = multem_in_parm.iw_psi;
			this->beam_pos_2d = multem_in_parm.beam_pos_2d;

			this->E_0 = multem_in_parm.E_0;
			this->theta = multem_in_parm.theta;
			this->phi = multem_in_parm.phi;
			this->nrot = multem_in_parm.nrot;

			this->illum_mod = multem_in_parm.illum_mod;
			this->illum_inc = multem_in_parm.illum_inc;

			this->cond_lens = multem_in_parm.cond_lens;
			this->obj_lens = multem_in_parm.obj_lens;

			this->scanning = multem_in_parm.scanning;
			this->detector = multem_in_parm.detector;

			this->eels_fr = multem_in_parm.eels_fr;

			this->cdl_var_type = multem_in_parm.cdl_var_type;
			this->cdl_var = multem_in_parm.cdl_var;

			this->beam_pos = multem_in_parm.beam_pos;

			this->islice = multem_in_parm.islice;
			this->dp_Shift = multem_in_parm.dp_Shift;
		}
	};

}

#endif