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
#include "types.hpp"
#include "fft2.cuh"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "input_multislice.hpp"

namespace multem
{
	template<class T, eDevice dev>
	class Probe{
		public:
			using value_type_r = T;
			using value_type_c = complex<T>;

			Probe():input_multislice(nullptr){}

			void set_input_data(Input_Multislice<value_type_r, dev> *input_multislice_i)
			{
				input_multislice = input_multislice_i;

				fft2.create_plan(input_multislice->grid.ny, input_multislice->grid.nx, input_multislice->nstream);

				probe_0.resize(input_multislice->grid.nxy());
			}

			template<class TVector_c>
			void get(const eSpace &space, TVector_c &host_probe)
			{
				value_type_r x = input_multislice->get_Rx_pos_shift();
				value_type_r y = input_multislice->get_Ry_pos_shift();

				multem::probe(input_multislice->grid, input_multislice->lens, x, y, probe_0);

				if(space == eS_Real)
				{
					fft2.inverse(probe_0);
				}

				multem::to_host_shift(input_multislice->grid, probe_0, host_probe);
			}

			void cleanup()
			{
				fft2.cleanup();
			}

		private:
			Input_Multislice<value_type_r, dev> *input_multislice;
			FFT2<value_type_r, dev> fft2;

			Vector<value_type_c, dev> probe_0;
	};

} // namespace multem

#endif