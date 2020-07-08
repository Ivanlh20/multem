/*
 * This file is part of MULTEM.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * MULTEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MULTEM is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef INCIDENT_WAVE_H
#define INCIDENT_WAVE_H

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "fft.cuh"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "cpu_fcns.hpp"
#include "gpu_fcns.cuh"

namespace mt
{
	template <class T, eDevice dev>
	class Incident_Wave{
		public:
			using T_r = T;
			using T_c = complex<T>;

			static const eDevice device = dev;

			Incident_Wave(): input_multislice(nullptr), stream(nullptr), fft_2d(nullptr){}

			void set_input_data(Input_Multislice<T_r> *input_multislice_i, Stream<dev> *stream_i, FFT<T_r, dev> *fft2_i)
			{
				input_multislice = input_multislice_i;
				stream = stream_i;
				fft_2d = fft2_i;

				if(input_multislice->is_user_define_wave())
				{
					fpsi_0.assign(input_multislice->iw_psi.begin(), input_multislice->iw_psi.end());
					mt::fft2_shift(*stream, input_multislice->grid_2d, fpsi_0);
					fft_2d->forward(fpsi_0);
				}
				else if(input_multislice->is_convergent_wave())
				{
					fpsi_0.resize(input_multislice->grid_2d.nxy());
				}
			}

			void operator()(Vector<T_c, dev> &psi, T_r gxu, T_r gyu, 
			Vector<T_r, e_host> &x_b, Vector<T_r, e_host> &y_b, T_r z_init=0)
			{
				switch(input_multislice->iw_type)
				{
					case eIWT_Plane_Wave:
					{
						mt::fill(*stream, psi, T_c(1.0, 0.0));
					}
					break;
					case eIWT_Convergent_Wave:
					{
						auto f_0 = input_multislice->cond_lens.c_10;
						auto f_s = f_0 - (input_multislice->cond_lens.zero_defocus_plane-z_init);
						input_multislice->cond_lens.set_defocus(f_s);

						mt::fill(*stream, psi, T_c(0));
						for(auto ib=0; ib<x_b.size(); ib++)
						{
							auto x = input_multislice->grid_2d.exp_factor_Rx(x_b[ib]);
							auto y = input_multislice->grid_2d.exp_factor_Ry(y_b[ib]);

							mt::probe(*stream, input_multislice->grid_2d, input_multislice->cond_lens, x, y, gxu, gyu, fpsi_0);
							mt::add(*stream, fpsi_0, psi);
						}
						fft_2d->inverse(psi);

						input_multislice->cond_lens.set_defocus(f_0);
					}
					break;
					case eIWT_User_Define_Wave:
					{
						// we need to include defocus
						auto f_s = -(input_multislice->cond_lens.zero_defocus_plane-z_init);

						Vector<T_r, e_host> x(x_b.size());
						Vector<T_r, e_host> y(y_b.size());
						for(auto ib=0; ib<x_b.size(); ib++)
						{
							x[ib] = input_multislice->grid_2d.exp_factor_Rx(x_b[ib]);
							y[ib] = input_multislice->grid_2d.exp_factor_Ry(y_b[ib]);
						}

						mt::mul_exp_g_factor_2d(*stream, input_multislice->grid_2d, x, y, fpsi_0, psi);

						fft_2d->inverse(psi);
					}
					break;
				}
			}

			template <class TOutput_multislice>
			void operator()(const eSpace &space, TOutput_multislice &output_multislice)
			{
				Vector<T_c, dev> psi(input_multislice->grid_2d.nxy());
				T_r gxu = input_multislice->gx_0();
				T_r gyu = input_multislice->gy_0();
				Vector<T_r, e_host> &iw_x = input_multislice->iw_x;
				Vector<T_r, e_host> &iw_y = input_multislice->iw_y;

				this->operator()(psi, gxu, gyu, iw_x, iw_y);

				if(space == eS_Reciprocal)
				{
					fft_2d->forward(psi);
					mt::scale(*stream, input_multislice->grid_2d.inxy(), psi);
				}

				mt::copy_to_host(output_multislice.stream, psi, output_multislice.psi_0[0]);
			}

		private:
			Input_Multislice<T_r> *input_multislice;
			Stream<dev> *stream;
			FFT<T_r, dev> *fft_2d;

			Vector<T_c, dev> fpsi_0;
	};

} // namespace mt

#endif