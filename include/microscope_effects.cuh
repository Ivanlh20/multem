/*
 * This file is part of Multem.
 * Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef MICROSCOPE_EFFECTS_H
#define MICROSCOPE_EFFECTS_H

#include "math_mt.h"
#include "types.cuh"
#include "fcns_cpu.h"
#include "fcns_gpu.h"
#include "fcns_gpu.h"
#include "quad_data.cuh"

namespace mt
{
	template <class T, eDev Dev>
	class Microscope_Effects
	{
		public:
			using T_r = T;
			using T_c = complex<T>;

			Microscope_Effects(): multem_in_parm(nullptr), stream(nullptr), fft_2d(nullptr) {}			
			
			void set_in_data(Multem_In_Parm<T_r> *multem_in_parm_i, Stream<Dev> *stream_i, FFT<T_r, Dev> *fft2_i)
			{
				multem_in_parm = multem_in_parm_i;
				stream = stream_i;
				fft_2d = fft2_i;

				psi.resize(multem_in_parm->grid_2d.size());

				if ((multem_in_parm->illum_mod == eim_coherent)||(multem_in_parm->illum_mod == eim_partial_coherent))
				{
					return;
				}

				// Load quadratures
				obj_lens_temporal_spatial_quadratures(multem_in_parm->obj_lens, qt, qs);
			}

			void operator()(Vctr<T_c, Dev>& fpsi, Vctr<T_r, Dev>& m2psi_tot)
			{
				switch(multem_in_parm->illum_mod)
				{
					case eim_coherent:
					{
						CTF_TEM(multem_in_parm->illum_inc, fpsi, m2psi_tot);
					}
					break;
					case eim_partial_coherent:
					{
						PCTF_LI_WPO_TEM(multem_in_parm->illum_inc, fpsi, m2psi_tot);
					}
					break;
					case eim_trans_cross_coef:
					{

					}
					break;
					case eim_full_integration:
					{
						num_int_TEM(multem_in_parm->illum_inc, fpsi, m2psi_tot);
					}
					break;
				}
			}

			template <class TOutput_multislice>
			void operator()(TOutput_multislice &output_multem)
			{
				Vctr<T_c, Dev> psi(multem_in_parm->iw_psi.begin(), multem_in_parm->iw_psi.end());
				mt::fcn_fftsft_2d(*stream, multem_in_parm->grid_2d, psi);
				fft_2d->forward(psi);
				mt::fcn_scale(*stream, multem_in_parm->grid_2d.isize_r(), psi);

				Vctr<T_r, Dev> m2psi_tot(multem_in_parm->grid_2d.size());
				this->operator()(psi, m2psi_tot);
				mt::cpy_to_host(output_multem.stream, m2psi_tot, output_multem.m2psi_tot[0]);
			}

		private:
			void CTF_TEM(const eIllum_Inc &illum_inc, Vctr<T_c, Dev>& fpsi, Vctr<T_r, Dev>& m2psi_tot)
			{
				mt::fcn_apply_ctf(*stream, multem_in_parm->grid_2d, multem_in_parm->obj_lens, 0, 0, fpsi, psi);
				fft_2d->inverse(psi);
				mt::square(*stream, psi, m2psi_tot);
			}

			void PCTF_LI_WPO_TEM(const eIllum_Inc &illum_inc, Vctr<T_c, Dev>& fpsi, Vctr<T_r, Dev>& m2psi_tot)
			{
				T_r tp_inc_sigma = multem_in_parm->obj_lens.tp_inc_sigma;
				T_r spt_inc_sigma = multem_in_parm->obj_lens.spt_inc_sigma;

				switch(illum_inc)
				{
					case eii_temporal:	// Temporal
					{
						multem_in_parm->obj_lens.set_spt_inc_sigma(0);
					}
					break;
					case etst_spatial:	// Spatial
					{
						multem_in_parm->obj_lens.set_tp_inc_sigma(0);
					}
					break;
				}

				mt::fcn_apply_pctf(*stream, multem_in_parm->grid_2d, multem_in_parm->obj_lens, fpsi, psi);
				fft_2d->inverse(psi);
				mt::square(*stream, psi, m2psi_tot);

				multem_in_parm->obj_lens.set_tp_inc_sigma(tp_inc_sigma);
				multem_in_parm->obj_lens.set_spt_inc_sigma(spt_inc_sigma);
			}

			void num_int_TEM(const eIllum_Inc &illum_inc, Vctr<T_c, Dev>& fpsi, Vctr<T_r, Dev>& m2psi_tot)
			{
				T_r c_10_0 = multem_in_parm->obj_lens.c_10;

				fill(*stream, m2psi_tot, 0.0);
				switch(illum_inc)
				{
					case 1:	// Temporal and Spatial
					{
						for(auto i = 0; i<qs.size(); i++)
						{
							for(auto j = 0; j<qt.size(); j++)
							{
								auto c_10 = multem_in_parm->obj_lens.tp_inc_iehwgd*qt.x[j]+c_10_0;
								multem_in_parm->obj_lens.set_defocus(c_10);
								
								mt::fcn_apply_ctf(*stream, multem_in_parm->grid_2d, multem_in_parm->obj_lens, qs.x[i], qs.y[i], fpsi, psi);
								fft_2d->inverse(psi);
								mt::add_scale_norm_2(*stream, qs.w[i]*qt.w[j], psi, m2psi_tot);
							}
						}
					}
					break;
					case 2:	// Temporal
					{
						for(auto j = 0; j<qt.size(); j++)
						{
							auto c_10 = multem_in_parm->obj_lens.tp_inc_iehwgd*qt.x[j]+c_10_0;
							multem_in_parm->obj_lens.set_defocus(c_10);

							mt::fcn_apply_ctf(*stream, multem_in_parm->grid_2d, multem_in_parm->obj_lens, 0.0, 0.0, fpsi, psi);
							fft_2d->inverse(psi);
							mt::add_scale_norm_2(*stream, qt.w[j], psi, m2psi_tot);
						}
					}
					break;
					case 3:	// Spatial
					{
						for(auto i = 0; i<qs.size(); i++)
						{
							mt::fcn_apply_ctf(*stream, multem_in_parm->grid_2d, multem_in_parm->obj_lens, qs.x[i], qs.y[i], fpsi, psi);
							fft_2d->inverse(psi);
							mt::add_scale_norm_2(*stream, qs.w[i], psi, m2psi_tot);
						}
					}
				}

				multem_in_parm->obj_lens.set_defocus(c_10_0);
			}
			
			Multem_In_Parm<T_r> *multem_in_parm;
			Stream<Dev> *stream;
			FFT<T_r, Dev> *fft_2d;

			Vctr<T_c, Dev> psi;

			Quad_Coef_1d<T_r, edev_cpu> qt;
			Quad_Coef_2d<T_r, edev_cpu> qs;
	};

}

#endif