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

#ifndef ENERGY_LOSS_H
#define ENERGY_LOSS_H

#include "types.cuh"
#include "stream.cuh"
#include "fft2.cuh"
#include "input_multislice.cuh"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"

namespace multem
{
	template<class T, eDevice dev>
	class Energy_Loss
	{
		public:
			using value_type_r = T;
			using value_type_c = complex<T>;

			void set_input_data(Input_Multislice<value_type_r> *input_multislice_i, Stream<dev> *stream_i, FFT2<value_type_r, dev> *fft2_i)
			{
				input_multislice = input_multislice_i;
				stream = stream_i;
				fft2 = fft2_i;

				if(input_multislice->eels_fr.m_selection>2)
				{
					kernel.resize(3);
				}
				else
				{
					kernel.resize(1);
				}

				for(auto ikn =0; ikn<kernel.size(); ikn++)
				{
					kernel[ikn].resize(input_multislice->grid.nxy());
				}

			}

			void set_atom_type(EELS<T> &eels)
			{
				if(eels.m_selection>2)
				{
					multem::kernel_xyz(*stream, input_multislice->grid, eels, *fft2, kernel[0], kernel[1], kernel[2]);
				}
				else if(eels.m_selection == -2)
				{
					multem::kernel_x(*stream, input_multislice->grid, eels, *fft2, kernel[0]);
				}
				else if(eels.m_selection == -1)
				{
					multem::kernel_mn1(*stream, input_multislice->grid, eels, *fft2, kernel[0]);
				}
				else if(eels.m_selection == 0)
				{
					multem::kernel_z(*stream, input_multislice->grid, eels, *fft2, kernel[0]);
				}
				else if(eels.m_selection == 1)
				{
					multem::kernel_mp1(*stream, input_multislice->grid, eels, *fft2, kernel[0]);
				}
				else if(eels.m_selection == 2)
				{
					multem::kernel_y(*stream, input_multislice->grid, eels, *fft2, kernel[0]);
				}
			}

			Vector<Vector<value_type_c, dev>, e_host> kernel;
		private:
			Input_Multislice<value_type_r> *input_multislice;
			Stream<dev> *stream;
			FFT2<value_type_r, dev> *fft2;
	};

} // namespace multem

#endif