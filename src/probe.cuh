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
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PROBE_H
#define PROBE_H

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "fft2.cuh"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"

namespace multem
{
	template<class T, eDevice dev>
	class Probe{
		public:
			using value_type_r = T;
			using value_type_c = complex<T>;

			Probe():input_multislice(nullptr), stream(nullptr), fft2(nullptr){}

			void set_input_data(Input_Multislice<value_type_r, dev> *input_multislice_i, Stream<dev> *stream_i, FFT2<value_type_r, dev> *fft2_i)
			{
				input_multislice = input_multislice_i;
				stream = stream_i;
				fft2 = fft2_i;

				probe_0.resize(input_multislice->grid.nxy());
			}

			template<class TOutput_multislice>
			void get(const eSpace &space, TOutput_multislice &output_multislice)
			{
				value_type_r x = input_multislice->get_Rx_pos_shift();
				value_type_r y = input_multislice->get_Ry_pos_shift();

				multem::probe(*stream, input_multislice->grid, input_multislice->lens, x, y, probe_0);

				if(space == eS_Real)
				{
					fft2->inverse(probe_0);
				}

				multem::copy_to_host(output_multislice.stream, input_multislice->grid, probe_0, output_multislice.probe[0]);
				output_multislice.shift();
				output_multislice.clear_temporal_data();
			}

		private:
			Input_Multislice<value_type_r, dev> *input_multislice;
			Stream<dev> *stream;
			FFT2<value_type_r, dev> *fft2;

			Vector<value_type_c, dev> probe_0;
	};

} // namespace multem

#endif