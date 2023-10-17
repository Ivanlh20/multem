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

#ifndef INCIDENT_WAVE_H
	#define INCIDENT_WAVE_H

	#include "math.cuh"
	#include "types.cuh"
	#include "type_traits_gen.cuh"
	#include "cgpu_stream.cuh"
	#include "cgpu_fft.cuh"
	#include "in_classes.cuh"
	#include "output_multem.hpp"
	#include "cpu_fcns.hpp"
	#include "gpu_fcns.cuh"

	namespace mt
	{
		template <class T, eDev Dev>
		class Incident_Wave{
			public:
				using T_r = T;
				using T_c = complex<T>;

				static const eDev device = Dev;

				Incident_Wave(): multem_in_parm(nullptr), stream(nullptr), fft_2d(nullptr) {}

				void set_in_data(Multem_In_Parm<T_r> *multem_in_parm_i, Stream<Dev> *stream_i, FFT<T_r, Dev> *fft2_i)
				{
					multem_in_parm = multem_in_parm_i;
					stream = stream_i;
					fft_2d = fft2_i;

					if (multem_in_parm->is_user_define_wave())
					{
						fpsi_0.assign(multem_in_parm->iw_psi.begin(), multem_in_parm->iw_psi.end());
						mt::fcn_fftsft_2d(*stream, multem_in_parm->grid_2d, fpsi_0);
						fft_2d->forward(fpsi_0);
					}
					else if (multem_in_parm->is_convergent_wave())
					{
						fpsi_0.resize(multem_in_parm->grid_2d.size());
					}
				}

				void operator()(Vctr<T_c, Dev>& psi, R_2d<T_r> gu, Beam_Pos_2d<T_r>& beam_pos_2d, T_r z_init=0)
				{
					switch(multem_in_parm->iw_type)
					{
						case eiwt_plane_wave:
						{
							mt::fill(*stream, psi, T_c(1.0, 0.0));
						}
						break;
						case eiwt_convergent_wave:
						{
							auto f_0 = multem_in_parm->cond_lens.c_10;
							auto f_s = f_0 - (multem_in_parm->cond_lens.zero_def_plane-z_init);
							multem_in_parm->cond_lens.set_defocus(f_s);

							mt::fill(*stream, psi, T_c(0));
							auto R = multem_in_parm->grid_2d.factor_2pi_rv_ctr(beam_pos_2d.p);

							for(auto ib=0; ib<R.size(); ib++)
							{
								mt::fcn_fs_probe(*stream, multem_in_parm->grid_2d, multem_in_parm->cond_lens, R[ib], gu, fpsi_0);
								mt::add(*stream, fpsi_0, psi);
							}
							fft_2d->inverse(psi);

							multem_in_parm->cond_lens.set_defocus(f_0);
						}
						break;
						case eiwt_user_def_Wave:
						{
							// we need to include defocus
							auto f_s = -(multem_in_parm->cond_lens.zero_def_plane-z_init);

							auto R = multem_in_parm->grid_2d.factor_2pi_rv_ctr(beam_pos_2d.p);

							mt::fcn_mul_exp_g_factor_2d(*stream, multem_in_parm->grid_2d, R, fpsi_0, psi);

							fft_2d->inverse(psi);
						}
						break;
					}
				}

				template <class TOutput_multislice>
				void operator()(const eSpace &space, TOutput_multislice &output_multem)
				{
					Vctr<T_c, Dev> psi(multem_in_parm->grid_2d.size());
					this->operator()(psi, multem_in_parm->gu_0(), multem_in_parm->beam_pos_2d);

					if (space == esp_fourier)
					{
						fft_2d->forward(psi);
						mt::fcn_scale(*stream, multem_in_parm->grid_2d.isize_r(), psi);
					}

					mt::fcn_assign_crop_fftsft_2d(multem_in_parm->grid_2d, psi, multem_in_parm->output_area, output_multem.psi_0(0, 0));
				}

			private:
				Multem_In_Parm<T_r> *multem_in_parm;
				Stream<Dev> *stream;
				FFT<T_r, Dev> *fft_2d;

				Vctr<T_c, Dev> fpsi_0;
		};

	}

#endif