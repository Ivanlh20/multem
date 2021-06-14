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

#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include "math.cuh"
#include "types.cuh"
#include "type_traits_gen.cuh"
#include "cgpu_stream.cuh"
#include "cgpu_fft.cuh"
#include "in_classes.cuh"
#include "output_multem.hpp"
#include "cpu_fcns.hpp"
#include "gpu_fcns.cuh"
#include "cgpu_fcns.cuh"

namespace mt
{
	template <class T, eDev Dev>
	class Propagator{
		public:
			using T_r = T;
			using T_c = complex<T>;

			static const eDev device = Dev;

			Propagator(): in_multem(nullptr), stream(nullptr), fft_2d(nullptr) {}

			void set_in_data(In_Multem<T_r> *in_multem_i, Stream<Dev> *stream_i, FFT<T_r, Dev> *fft2_i)
			{
				in_multem = in_multem_i;
				stream = stream_i;
				fft_2d = fft2_i;
			}

			void operator()(const eSpace &space_out, T_r gxu, T_r gyu, 
			T_r z, Vctr<T_c, Dev>& psi_i, Vctr<T_c, Dev>& psi_o)
			{
				if (fcn_is_zero(z))
				{
					if (in_multem->grid_2d.bwl)
					{
						fft_2d->forward(psi_i, psi_o);
						mt::fcn_fermi_aperture(*stream, in_multem->grid_2d, psi_o);
						
						if (space_out == eS_Real)
						{
							fft_2d->inverse(psi_o);
						}
					}
					else
					{
						if (space_out == eS_Reciprocal)
						{
							fft_2d->forward(psi_i, psi_o);
							mt::fcn_scale(*stream, in_multem->grid_2d.isize_r(), psi_o);
						}
					}
				}
				else
				{
					fft_2d->forward(psi_i, psi_o);

					mt::fcn_propagate(*stream, in_multem->grid_2d, in_multem->get_propagator_factor(z), gxu, gyu, psi_o, psi_o);

					if (space_out == eS_Real)
					{
						fft_2d->inverse(psi_o, psi_o);
					}
				}
			}

			void operator()(const eSpace &space_out, T_r gxu, T_r gyu, 
			T_r z, Vctr<T_c, Dev>& psi_io)
			{
				this->operator()(space_out, gxu, gyu, z, psi_io, psi_io);
			}

			template <class TOutput_multislice>
			void operator()(const eSpace &space_out, T_r gxu, T_r gyu, 
			T_r z, TOutput_multislice &output_multem)
			{
				Vctr<T_c, Dev> psi(in_multem->iw_psi.begin(), in_multem->iw_psi.end());
				mt::fcn_fftsft_2d(*stream, in_multem->grid_2d, psi);
				this->operator()(space_out, gxu, gyu, z, psi);
				mt::cpy_to_host(output_multem.stream, psi, output_multem.psi_coh[0]);
			}

		private:
			In_Multem<T_r> *in_multem;
			Stream<Dev> *stream;
			FFT<T_r, Dev> *fft_2d;
	};

}

#endif