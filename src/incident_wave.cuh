/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
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
#include "fft2.cuh"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "host_functions.hpp"
#include "device_functions.cuh"

namespace multem
{
	template<class T, eDevice dev>
	class Incident_Wave{
		public:
			using value_type_r = T;
			using value_type_c = complex<T>;

			Incident_Wave(): input_multislice(nullptr), stream(nullptr), fft2(nullptr){}

			void set_input_data(Input_Multislice<value_type_r> *input_multislice_i, Stream<dev> *stream_i, FFT2<value_type_r, dev> *fft2_i)
			{
				input_multislice = input_multislice_i;
				stream = stream_i;
				fft2 = fft2_i;

				if(input_multislice->is_user_define_wave())
				{
					fpsi_0.assign(input_multislice->iw_psi.begin(), input_multislice->iw_psi.end());
					multem::fft2_shift(*stream, input_multislice->grid, fpsi_0);
					fft2->forward(fpsi_0);
				}
			}

			void operator()(Vector<value_type_c, dev> &psi, value_type_r z_init=0)
			{
				switch(input_multislice->iw_type)
				{
					case eIWT_Plane_Wave:
					{
						multem::fill(*stream, psi, value_type_c(1.0, 0.0));
					}
					break;
					case eIWT_Convergent_Wave:
					{
						value_type_r x = input_multislice->get_Rx_pos_shift();
						value_type_r y = input_multislice->get_Ry_pos_shift();

						auto f_0 = input_multislice->cond_lens.f;
						auto f_s = f_0 - (input_multislice->cond_lens.zero_defocus_plane-z_init);
						input_multislice->cond_lens.set_defocus(f_s);

						multem::probe(*stream, input_multislice->grid, input_multislice->cond_lens, x, y, psi);
						fft2->inverse(psi);

						input_multislice->cond_lens.set_defocus(f_0);
					}
					break;
					case eIWT_User_Define_Wave:
					{
						value_type_r x = input_multislice->get_Rx_pos_shift();
						value_type_r y = input_multislice->get_Ry_pos_shift();

						multem::phase_factor_2D(*stream, input_multislice->grid, x, y, fpsi_0, psi);
						fft2->inverse(psi);
					}
					break;
				}
			}

			template<class TOutput_multislice>
			void operator()(const eSpace &space, TOutput_multislice &output_multislice)
			{
				Vector<value_type_c, dev> psi(input_multislice->grid.nxy());
				this->operator()(psi);

				if(space == eS_Reciprocal)
				{
					fft2->forward(psi);
					multem::scale(*stream, input_multislice->grid.inxy, psi);
				}

				multem::copy_to_host(output_multislice.stream, psi, output_multislice.psi_0[0]);
				output_multislice.shift();
				output_multislice.clear_temporal_data();
			}

		private:
			Input_Multislice<value_type_r> *input_multislice;
			Stream<dev> *stream;
			FFT2<value_type_r, dev> *fft2;

			Vector<value_type_c, dev> fpsi_0;
	};

} // namespace multem

#endif