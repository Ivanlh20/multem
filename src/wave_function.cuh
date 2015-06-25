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

#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H

#include "types.hpp"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "input_multislice.hpp"
#include "transmission.cuh"
#include "propagator.cuh"

namespace multem
{
	template<class T, eDevice dev>
	class Wave_Function: public Transmission<T, dev>
	{
		public:
			using value_type_r = T;
			using value_type_c = complex<T>;

			void set_input_data(Input_Multislice<value_type_r, dev> *input_multislice_io, Stream<value_type_r, dev> *stream_i, FFT2<value_type_r, dev> *fft2_i)
			{
				exp_x.resize(input_multislice_io->grid.nx);
				exp_y.resize(input_multislice_io->grid.ny);

				psi_z.resize(input_multislice_io->grid.nxy());

				prog.set_input_data(input_multislice_io, fft2_i);

				Transmission<T, dev>::set_input_data(input_multislice_io, stream_i, fft2_i);
			}

			void set_plane_wave(Vector<value_type_c, dev> &psi_z)
			{
				multem::fill(psi_z, value_type_c(1.0, 0.0));
			}

			void set_user_input_wave(Vector<value_type_c, dev> &psi_z)
			{
				multem::assign(input_multislice->psi_0, psi_z);
			}

			void set_conv_beam_wave(Vector<value_type_c, dev> &psi_z)
			{
				value_type_r x = input_multislice->get_Rx_pos_shift();
				value_type_r y = input_multislice->get_Ry_pos_shift();

				multem::probe(input_multislice->grid, input_multislice->lens, x, y, psi_z);
				fft2->inverse(psi_z);
			}

			void phase_mul(const value_type_r &gxu, const value_type_r &gyu, Vector<value_type_c, dev> &psi_i, Vector<value_type_c, dev> &psi_o)
			{
				if(input_multislice->dp_Shift || isZero(gxu, gyu))
				{
					if (psi_i.data() != psi_o.data())
					{
						psi_o.assign(psi_i.begin(), psi_i.end());
					}
					return;
				}

				multem::phase_components(input_multislice->grid, gxu, gyu, exp_x, exp_y);
				multem::phase_mul(input_multislice->grid, exp_x, exp_y, psi_i, psi_o);
			}	

			void phase_mul(const value_type_r &gxu, const value_type_r &gyu, Vector<value_type_c, dev> &psi_io)
			{
				phase_mul(gxu, gyu, psi_io, psi_io);
			}

			void psi_0(Vector<value_type_c, dev> &psi_z)
			{
				switch(input_multislice->beam_type)
				{
					case eBT_Plane_Wave:
					{
						set_plane_wave(psi_z);
					}
					break;
					case eBT_Convergent:
					{
						set_conv_beam_wave(psi_z);
					}
					break;
					case eBT_User_Define:
					{
						set_user_input_wave(psi_z);
					}
					break;
				}
			}

			void psi_0()
			{
				psi_0(psi_z);
			}

			void psi(const eSpace &space)
			{
				value_type_r gx_0 = input_multislice->gx_0();
				value_type_r gy_0 = input_multislice->gy_0();

				if(input_multislice->is_multislice())
				{
					for(auto islice=0; islice<slice.size(); islice++)
					{
						transmit(islice, psi_z);
						prog.propagate(eS_Real, gx_0, gy_0, dz(islice), psi_z);
					}
					phase_mul(gx_0, gy_0, psi_z);
					prog.propagate(space, gx_0, gy_0, thickness.z_back_prop[0], psi_z);
				}
				else
				{
					transmit(0, psi_z);
					phase_mul(gx_0, gy_0, psi_z);

					if(space == eS_Reciprocal)
					{
						fft2->forward(psi_z);
						multem::scale(psi_z, input_multislice->grid.inxy);
					}
				}
			}

			void psi(const eSpace &space, int islice_0, int islice_e, Vector<value_type_c, dev> &psi_z)
			{
				value_type_r gx_0 = input_multislice->gx_0();
				value_type_r gy_0 = input_multislice->gy_0();

				if(input_multislice->is_multislice())
				{
					for(auto islice=islice_0; islice<islice_e; islice++)
					{
						transmit(islice, psi_z);
						prog.propagate(eS_Real, gx_0, gy_0, dz(islice), psi_z);
					}
					phase_mul(gx_0, gy_0, psi_z);
					prog.propagate(space, gx_0, gy_0, 0, psi_z);
				}
				else
				{
					transmit(0, psi_z);
					phase_mul(gx_0, gy_0, psi_z);

					if(space == eS_Reciprocal)
					{
						fft2->forward(psi_z);
						multem::scale(psi_z, input_multislice->grid.inxy);
					}
				}
			}

			Vector<value_type_c, dev> psi_z;
			Propagator<value_type_r, dev> prog;
		private:
			Vector<value_type_c, dev> exp_x;
			Vector<value_type_c, dev> exp_y;
	};

} // namespace multem

#endif