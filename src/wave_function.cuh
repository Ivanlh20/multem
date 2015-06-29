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
			void set_input_data(Input_Multislice<value_type_r, dev> *input_multislice_io, Stream<value_type_r, dev> *stream_i, FFT2<value_type_r, dev> *fft2_i)
			{
				exp_x.resize(input_multislice_io->grid.nx);
				exp_y.resize(input_multislice_io->grid.ny);

				psi_z.resize(input_multislice_io->grid.nxy());

				prog.set_input_data(input_multislice_io, fft2_i);

				Transmission::set_input_data(input_multislice_io, stream_i, fft2_i);
			}

			void set_plane_wave()
			{
				multem::fill(psi_z, value_type_c(1.0, 0.0));
			}

			void set_user_input_wave()
			{
				multem::assign(input_multislice->psi_0, psi_z);
			}

			void set_conv_beam_wave()
			{
				value_type_r x = input_multislice->get_Rx_pos_shift();
				value_type_r y = input_multislice->get_Ry_pos_shift();

				multem::probe(input_multislice->grid, input_multislice->lens, x, y, psi_z);
				fft2->inverse(psi_z);
			}

			void phase_mul(const value_type_r &gxu, const value_type_r &gyu, Vector<value_type_c, dev> &psi_io)
			{
				if(input_multislice->dp_Shift || isZero(gxu, gyu))
				{
					return;
				}

				multem::phase_component(input_multislice->grid, gxu, gyu, exp_x, exp_y);
				multem::phase_mul(input_multislice->grid, exp_x, exp_y, psi_io, psi_io);
			}

			void psi_0()
			{
				switch(input_multislice->beam_type)
				{
					case eBT_Plane_Wave:
					{
						set_plane_wave();
					}
					break;
					case eBT_Convergent:
					{
						set_conv_beam_wave();
					}
					break;
					case eBT_User_Define:
					{
						set_user_input_wave();
					}
					break;
				}
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
						prog.propagate(eS_Real, gx_0, gy_0, get_dz(islice), psi_z);
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

			template<class TVector>
			void psi(const eSpace &space, int islice_0, int islice_e, TVector *psi0=nullptr)
			{
				if(psi0!=nullptr)
				{
					multem::assign(*psi0, psi_z);
				}

				value_type_r gx_0 = input_multislice->gx_0();
				value_type_r gy_0 = input_multislice->gy_0();

				if(input_multislice->is_multislice())
				{
					for(auto islice=islice_0; islice<islice_e; islice++)
					{
						transmit(islice, psi_z);
						prog.propagate(eS_Real, gx_0, gy_0, get_dz(islice), psi_z);
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

			void psi(const eSpace &space, const int &islice_e)
			{
				decltype(psi_z) *psi_t = nullptr;
				psi(space, 0, islice_e, psi_t);
			}

			Vector<value_type_c, dev> psi_z;
			Propagator<value_type_r, dev> prog;
		private:
			Vector<value_type_c, dev> exp_x;
			Vector<value_type_c, dev> exp_y;
	};

} // namespace multem

#endif