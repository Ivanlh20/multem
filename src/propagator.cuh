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
#include "types.hpp"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "input_multislice.hpp"

namespace multem
{
	template<class T, eDevice dev>
	class Propagator{
		public:
			using value_type_r = typename T;
			using value_type_c = typename complex<T>;

			Propagator():input_multislice(nullptr), fft2(nullptr){ }

			void set_input_data(Input_Multislice<value_type_r, dev> *input_multislice_io, FFT2<value_type_r, dev> *fft2_i)
			{
				input_multislice = input_multislice_io;
				fft2 = fft2_i;

				prop_x.resize(input_multislice->grid.nx);
				prop_y.resize(input_multislice->grid.ny);
			}

			void propagate(const eSpace &space, value_type_r gxu, value_type_r gyu, 
			value_type_r z, Vector<value_type_c, dev> &psi_i, Vector<value_type_c, dev> &psi_o)
			{
				if(isZero(z))
				{
					fft2->forward(psi_i, psi_o); 

					if(input_multislice->grid.bwl)
					{
						multem::bandwidth_limit(input_multislice->grid, 0, input_multislice->grid.gl_max, input_multislice->grid.inxy, psi_o);
					}
					else
					{
						multem::scale(psi_o, input_multislice->grid.inxy);
					}

					if(space == eS_Real)
					{
						fft2->inverse(psi_o);
					}
				}
				else
				{
					multem::propagator_components(input_multislice->grid, gxu, gyu, input_multislice->lens.prop_factor(z), prop_x, prop_y);
					multem::propagate(input_multislice->grid, *fft2, space, prop_x, prop_y, psi_i, psi_o);
				}
			}

			void propagate(const eSpace &space, value_type_r gxu, value_type_r gyu, 
			value_type_r z, Vector<value_type_c, dev> &psi_io)
			{
				propagate(space, gxu, gyu, z, psi_io, psi_io);
			}

		private:
			Input_Multislice<value_type_r, dev> *input_multislice;
			FFT2<value_type_r, dev> *fft2;

			Vector<value_type_c, dev> prop_x;
			Vector<value_type_c, dev> prop_y;
	};

} // namespace multem

#endif