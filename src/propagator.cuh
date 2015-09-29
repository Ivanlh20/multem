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

#ifndef PROPAGATOR_H
#define PROPAGATOR_H

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
	class Propagator{
		public:
			using value_type_r = T;
			using value_type_c = complex<T>;

			Propagator():input_multislice(nullptr), stream(nullptr), fft2(nullptr){}

			void set_input_data(Input_Multislice<value_type_r, dev> *input_multislice_i, Stream<dev> *stream_i, FFT2<value_type_r, dev> *fft2_i)
			{
				input_multislice = input_multislice_i;
				stream = stream_i;
				fft2 = fft2_i;

				prop_x.resize(input_multislice->grid.nx);
				prop_y.resize(input_multislice->grid.ny);
			}

			void propagate(const eSpace &space, value_type_r gxu, value_type_r gyu, 
			value_type_r z, Vector<value_type_c, dev> &psi_i, Vector<value_type_c, dev> &psi_o)
			{
				if(isZero(z))
				{
					if(input_multislice->grid.bwl)
					{
						fft2->forward(psi_i, psi_o); 
						multem::bandwidth_limit(*stream, input_multislice->grid, 0, input_multislice->grid.gl_max, input_multislice->grid.inxy, psi_o);
						if(space == eS_Real)
						{
							fft2->inverse(psi_o);
						}
					}
					else
					{
						if(space == eS_Reciprocal)
						{
							fft2->forward(psi_i, psi_o); 
							multem::scale(psi_o, input_multislice->grid.inxy);
						}
					}
				}
				else
				{
					multem::propagator_components(input_multislice->grid, gxu, gyu, input_multislice->lens.prop_factor(z), prop_x, prop_y);
					multem::propagate(*stream, *fft2, input_multislice->grid, space, prop_x, prop_y, psi_i, psi_o);
				}
			}

			void propagate(const eSpace &space, value_type_r gxu, value_type_r gyu, 
			value_type_r z, Vector<value_type_c, dev> &psi_io)
			{
				propagate(space, gxu, gyu, z, psi_io, psi_io);
			}

			template<class TOutput_multislice>
			void propagate(const eSpace &space, value_type_r gxu, value_type_r gyu, 
			value_type_r z, TOutput_multislice &output_multislice)
			{
				propagate(space, gxu, gyu, z, input_multislice->psi_0, input_multislice->psi_0);
				multem::copy_to_host(output_multislice.stream, input_multislice->grid, input_multislice->psi_0, output_multislice.psi_coh[0]);
				output_multislice.shift();
				output_multislice.clear_temporal_data();
			}

		private:
			Input_Multislice<value_type_r, dev> *input_multislice;
			Stream<dev> *stream;
			FFT2<value_type_r, dev> *fft2;

			Vector<value_type_c, dev> prop_x;
			Vector<value_type_c, dev> prop_y;
	};

} // namespace multem

#endif