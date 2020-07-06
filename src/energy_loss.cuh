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

#ifndef ENERGY_LOSS_H
#define ENERGY_LOSS_H

#include "types.cuh"
#include "stream.cuh"
#include "fft.cuh"
#include "input_multislice.cuh"
#include "cpu_fcns.hpp"
#include "gpu_fcns.cuh"
#include "cgpu_fcns.cuh"

namespace mt
{
	template <class T, eDevice dev>
	class Energy_Loss
	{
		public:
			using T_r = T;
			using T_c = complex<T>;

			void set_input_data(Input_Multislice<T_r> *input_multislice_i, Stream<dev> *stream_i, FFT<T_r, dev> *fft2_i)
			{
				input_multislice = input_multislice_i;
				stream = stream_i;
				fft_2d = fft2_i;

				if(input_multislice->eels_fr.m_selection>2)
				{
					kernel.resize(3);
				}
				else
				{
					kernel.resize(1);
				}

				for(auto ikn = 0; ikn<kernel.size(); ikn++)
				{
					kernel[ikn].resize(input_multislice->grid_2d.nxy());
				}

			}

			void set_atom_type(EELS<T> &eels)
			{
				if(eels.m_selection>2)
				{
					mt::kernel_xyz(*stream, input_multislice->grid_2d, eels, *fft_2d, kernel[0], kernel[1], kernel[2]);
				}
				else if(eels.m_selection == -2)
				{
					mt::kernel_x(*stream, input_multislice->grid_2d, eels, *fft_2d, kernel[0]);
				}
				else if(eels.m_selection == -1)
				{
					mt::kernel_mn1(*stream, input_multislice->grid_2d, eels, *fft_2d, kernel[0]);
				}
				else if(eels.m_selection == 0)
				{
					mt::kernel_z(*stream, input_multislice->grid_2d, eels, *fft_2d, kernel[0]);
				}
				else if(eels.m_selection == 1)
				{
					mt::kernel_mp1(*stream, input_multislice->grid_2d, eels, *fft_2d, kernel[0]);
				}
				else if(eels.m_selection == 2)
				{
					mt::kernel_y(*stream, input_multislice->grid_2d, eels, *fft_2d, kernel[0]);
				}
			}

			Vector<Vector<T_c, dev>, e_host> kernel;
		private:
			Input_Multislice<T_r> *input_multislice;
			Stream<dev> *stream;
			FFT<T_r, dev> *fft_2d;
	};

} // namespace mt

#endif