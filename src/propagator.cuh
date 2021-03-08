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

#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "fft.cuh"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "cpu_fcns.hpp"
#include "gpu_fcns.cuh"
#include "cgpu_fcns.cuh"

namespace mt
{
	template <class T, eDevice dev>
	class Propagator{
		public:
			using T_r = T;
			using T_c = complex<T>;

			static const eDevice device = dev;

			Propagator(): input_multislice(nullptr), stream(nullptr), fft_2d(nullptr){}

			void set_input_data(Input_Multislice<T_r> *input_multislice_i, Stream<dev> *stream_i, FFT<T_r, dev> *fft2_i)
			{
				input_multislice = input_multislice_i;
				stream = stream_i;
				fft_2d = fft2_i;
			}

			void operator()(const eSpace &space_out, T_r gxu, T_r gyu, 
			T_r z, Vector<T_c, dev> &psi_i, Vector<T_c, dev> &psi_o)
			{
				if(isZero(z))
				{
					if(input_multislice->grid_2d.bwl)
					{
						fft_2d->forward(psi_i, psi_o); 
						mt::bandwidth_limit(*stream, input_multislice->grid_2d, psi_o);
						
						if(space_out == eS_Real)
						{
							fft_2d->inverse(psi_o);
						}
					}
					else
					{
						if(space_out == eS_Reciprocal)
						{
							fft_2d->forward(psi_i, psi_o);
							mt::scale(*stream, input_multislice->grid_2d.inxy(), psi_o);
						}
					}
				}
				else
				{
					fft_2d->forward(psi_i, psi_o); 

					mt::propagate(*stream, input_multislice->grid_2d, input_multislice->get_propagator_factor(z), gxu, gyu, psi_o, psi_o);

					if(space_out == eS_Real)
					{
						fft_2d->inverse(psi_o, psi_o);
					}
				}
			}

			void operator()(const eSpace &space_out, T_r gxu, T_r gyu, 
			T_r z, Vector<T_c, dev> &psi_io)
			{
				this->operator()(space_out, gxu, gyu, z, psi_io, psi_io);
			}

			template <class TOutput_multislice>
			void operator()(const eSpace &space_out, T_r gxu, T_r gyu, 
			T_r z, TOutput_multislice &output_multislice)
			{
				Vector<T_c, dev> psi(input_multislice->iw_psi.begin(), input_multislice->iw_psi.end());
				mt::fft2_shift(*stream, input_multislice->grid_2d, psi);
				this->operator()(space_out, gxu, gyu, z, psi);
				mt::fft2_shift(*stream, input_multislice->grid_2d, psi);
				mt::copy_to_host(output_multislice.stream, psi, output_multislice.psi_coh[0]);
			}

		private:
			Input_Multislice<T_r> *input_multislice;
			Stream<dev> *stream;
			FFT<T_r, dev> *fft_2d;
	};

} // namespace mt

#endif