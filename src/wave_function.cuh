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
#include "fft2.cuh"
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
			using size_type = std::size_t;

			static const eDevice device = dev;

			void set_input_data(Input_Multislice<value_type_r, dev> *input_multislice_i, Stream<value_type_r, dev> *stream_i, FFT2<value_type_r, dev> *fft2_i)
			{
				exp_x.resize(input_multislice_i->grid.nx);
				exp_y.resize(input_multislice_i->grid.ny);

				psi_z.resize(input_multislice_i->grid.nxy());

				if(input_multislice_i->device==e_Device)
				{
					psi_h.resize(input_multislice_i->grid.nxy());
					m2psi_h.resize(input_multislice_i->grid.nxy());
				}

				prog.set_input_data(input_multislice_i, fft2_i);

				Transmission<T, dev>::set_input_data(input_multislice_i, stream_i, fft2_i);
			}

			void set_plane_wave(Vector<value_type_c, dev> &psi_z)
			{
				multem::fill(psi_z, value_type_c(1.0, 0.0));
			}

			void set_user_input_wave(Vector<value_type_c, dev> &psi_z)
			{
				multem::assign(this->input_multislice->psi_0, psi_z);
			}

			void set_conv_beam_wave(Vector<value_type_c, dev> &psi_z)
			{
				value_type_r x = this->input_multislice->get_Rx_pos_shift();
				value_type_r y = this->input_multislice->get_Ry_pos_shift();

				multem::probe(this->input_multislice->grid, this->input_multislice->lens, x, y, psi_z);
				this->fft2->inverse(psi_z);
			}

			void phase_multiplication(const value_type_r &gxu, const value_type_r &gyu, Vector<value_type_c, dev> &psi_i, Vector<value_type_c, dev> &psi_o)
			{
				if(this->input_multislice->dp_Shift || isZero(gxu, gyu))
				{
					if (psi_i.data() != psi_o.data())
					{
						psi_o.assign(psi_i.begin(), psi_i.end());
					}
					return;
				}

				multem::phase_components(this->input_multislice->grid, gxu, gyu, exp_x, exp_y);
				multem::phase_multiplication(this->input_multislice->grid, exp_x, exp_y, psi_i, psi_o);
			}	

			void phase_multiplication(const value_type_r &gxu, const value_type_r &gyu, Vector<value_type_c, dev> &psi_io)
			{
				phase_multiplication(gxu, gyu, psi_io, psi_io);
			}

			void psi_0(Vector<value_type_c, dev> &psi_z)
			{
				switch(this->input_multislice->beam_type)
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

			template<class TVector_c, class TVector_Host_vc>
			void assign_psi(TVector_c &psi, const eSpace &space, const value_type_r &gxu, const value_type_r &gyu, 
			const int &islice, TVector_Host_vc &psi_v)
			{
				int ithk = this->slice.ithk[islice];
				if(0 <= ithk)
				{
					Vector<value_type_c, dev> *psi_zt = &(this->trans_0);
					phase_multiplication(gxu, gyu, psi_z, *psi_zt);
					prog.propagate(space, gxu, gyu, this->thickness.z_back_prop[ithk], *psi_zt);
					multem::copy_to_host(this->input_multislice->grid, *psi_zt, psi_v[ithk], &psi_h);
				}
			}

			template<class TVector_c, class TVector_Host_vc>
			void add_psi(TVector_c &psi, const eSpace &space, const value_type_r &gxu, const value_type_r &gyu, 
			const int &islice, value_type_r w_i, TVector_Host_vc &psi_v)
			{
				int ithk = this->slice.ithk[islice];
				if(0 <= ithk)
				{
					Vector<value_type_c, dev> *psi_zt = &(this->trans_0);
					phase_multiplication(gxu, gyu, psi_z, *psi_zt);
					prog.propagate(space, gxu, gyu, this->thickness.z_back_prop[ithk], *psi_zt);
					multem::add_scale_to_host(this->input_multislice->grid, w_i, *psi_zt, psi_v[ithk], &psi_h);
				}
			}

			template<class TVector_c, class TVector_Host_vr>
			void add_m2psi(TVector_c &psi, const eSpace &space, const value_type_r &gxu, const value_type_r &gyu, 
			const int &islice, value_type_r w_i, TVector_Host_vr &m2psi_v)
			{
				int ithk = this->slice.ithk[islice];
				if(0 <= ithk)
				{
					Vector<value_type_c, dev> *psi_zt = &(this->trans_0);
					phase_multiplication(gxu, gyu, psi_z, *psi_zt);
					prog.propagate(space, gxu, gyu, this->thickness.z_back_prop[ithk], *psi_zt);
					multem::add_square_scale_to_host(this->input_multislice->grid, w_i, *psi_zt, m2psi_v[ithk], &psi_h);
				}
			}

			template<class TVector_c, class TVector_Host_vr, class TVector_Host_vc>
			void add_m2psi_psi(TVector_c &psi, const eSpace &space, const value_type_r &gxu, const value_type_r &gyu, 
			const int &islice, value_type_r w_i, TVector_Host_vr &m2psi_v, TVector_Host_vc &psi_v)
			{
				int ithk = this->slice.ithk[islice];
				if(0 <= ithk)
				{
					Vector<value_type_c, dev> *psi_zt = &(this->trans_0);
					phase_multiplication(gxu, gyu, psi_z, *psi_zt);
					prog.propagate(space, gxu, gyu, this->thickness.z_back_prop[ithk], *psi_zt);
					if(!this->input_multislice->coherent_contribution)
					{
						multem::add_square_scale_to_host(this->input_multislice->grid, w_i, *psi_zt, m2psi_v[ithk], &psi_h);
					}
					else
					{
						multem::add_scale_m2psi_psi(this->input_multislice->grid, w_i, *psi_zt, m2psi_v[ithk], psi_v[ithk], &psi_h);
					}		
				}
			}

			template<class TVector_c>
			void psi_slice(const value_type_r &gxu, const value_type_r &gyu, const int &islice, TVector_c &psi_z)
			{
				this->transmit(islice, psi_z);
				if(this->input_multislice->is_multislice())
				{
					prog.propagate(eS_Real, gxu, gyu, this->dz(islice), psi_z);
				}
			}

			template<class TVector_Host_vc>
			void psi(const eSpace &space, TVector_Host_vc &host_psi_v)
			{
				value_type_r gx_0 = this->input_multislice->gx_0();
				value_type_r gy_0 = this->input_multislice->gy_0();

				for(auto islice=0; islice<this->slice.size(); islice++)
				{
					psi_slice(gx_0, gy_0, islice, psi_z);
					assign_psi(psi_z, space, gx_0, gy_0, islice, host_psi_v);
				}

				for(auto ithk=0; ithk<host_psi_v.size(); ithk++)
				{
					multem::fft2_shift(this->input_multislice->grid, host_psi_v[ithk]);
				}
			}

			//void psi(const eSpace &space, int islice_0, int islice_e, Vector<value_type_c, dev> &psi_z)
			//{
			//	value_type_r gx_0 = this->input_multislice->gx_0();
			//	value_type_r gy_0 = this->input_multislice->gy_0();

			//	if(this->input_multislice->is_multislice())
			//	{
			//		for(auto islice=islice_0; islice<islice_e; islice++)
			//		{
			//			this->transmit(islice, psi_z);
			//			prog.propagate(eS_Real, gx_0, gy_0, dz(islice), psi_z);
			//		}
			//		phase_multiplication(gx_0, gy_0, psi_z);
			//		prog.propagate(space, gx_0, gy_0, 0, psi_z);
			//	}
			//	else
			//	{
			//		this->transmit(0, psi_z);
			//		phase_multiplication(gx_0, gy_0, psi_z);

			//		if(space == eS_Reciprocal)
			//		{
			//			this->fft2->forward(psi_z);
			//			multem::scale(psi_z, this->input_multislice->grid.inxy);
			//		}
			//	}
			//}

			Propagator<value_type_r, dev> prog;

			Vector<value_type_c, dev> psi_z;
			Vector<value_type_c, e_Host> psi_h;
			Vector<value_type_r, e_Host> m2psi_h;
		private:
			Vector<value_type_c, dev> exp_x;
			Vector<value_type_c, dev> exp_y;
	};

} // namespace multem

#endif