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

#include "types.cuh"
#include "fft2.cuh"
#include "input_multislice.hpp"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "microscope_effects.cuh"
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

			void set_input_data(Input_Multislice<value_type_r, dev> *input_multislice_i, Stream<dev> *stream_i, FFT2<value_type_r, dev> *fft2_i)
			{
				exp_x.resize(input_multislice_i->grid.nx);
				exp_y.resize(input_multislice_i->grid.ny);

				psi_z.resize(input_multislice_i->grid.nxy());
				m2psi_z.resize(input_multislice_i->grid.nxy());

				if(input_multislice_i->is_device())
				{
					psi_zh.resize(input_multislice_i->grid.nxy());
					m2psi_zh.resize(input_multislice_i->grid.nxy());
				}

				prog.set_input_data(input_multislice_i, stream_i, fft2_i);

				if(input_multislice_i->is_HRTEM_HCI_EFTEM())
				{
					microscope_effects.set_input_data(input_multislice_i, stream_i, fft2_i);
				}

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

				//value_type_r f0 = this->input_multislice->lens.f;
				//value_type_r f = f0 - (slice.size()>0)?slice.z_int_0[0]:0;
				//this->input_multislice->lens.set_defocus(f);
				multem::probe(*(this->stream), this->input_multislice->grid, this->input_multislice->lens, x, y, psi_z);
				//this->input_multislice->lens.set_defocus(f0);
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
				multem::phase_multiplication(*(this->stream), this->input_multislice->grid, exp_x, exp_y, psi_i, psi_o);
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

			Vector<value_type_c, dev>* get_psi(const int &ithk, const eSpace &space, const value_type_r &gxu, const value_type_r &gyu)
			{
				Vector<value_type_c, dev> *psi_zt = &(this->trans_0);
				phase_multiplication(gxu, gyu, psi_z, *psi_zt);
				prog.propagate(space, gxu, gyu, this->thickness.z_back_prop[ithk], *psi_zt);

				return psi_zt;
			}

			template<class TVector_c, class TOutput_multislice>
			void set_m2psi_tot_psi_coh(TVector_c &psi, const value_type_r &gxu, const value_type_r &gyu, 
			const int &islice, const value_type_r &w_i, TOutput_multislice &output_multislice)
			{
				int ithk = this->slice.ithk[islice];
				if(0 <= ithk)
				{
					Vector<value_type_c, dev> *psi_zt = get_psi(ithk, this->input_multislice->get_simulation_space(), gxu, gyu);
					if(this->input_multislice->is_EWSFS_EWSRS())
					{
						multem::copy_to_host(output_multislice.stream, this->input_multislice->grid, *psi_zt, output_multislice.psi_coh[ithk], &psi_zh);
					}
					else if(this->input_multislice->is_HRTEM_HCI_EFTEM())
					{
						microscope_effects.apply(*psi_zt, m2psi_z);
						multem::add_scale_to_host(output_multislice.stream, this->input_multislice->grid, w_i, m2psi_z, output_multislice.m2psi_tot[ithk], &m2psi_zh);
						if(this->input_multislice->coherent_contribution)
						{
							multem::add_scale_to_host(output_multislice.stream, this->input_multislice->grid, w_i, *psi_zt, output_multislice.psi_coh[ithk], &psi_zh);
						}
					}
					else if(this->input_multislice->is_STEM())
					{
						for(auto iDet=0; iDet<this->input_multislice->det_cir.size(); iDet++)
						{
							value_type_r g_inner = this->input_multislice->det_cir.g_inner(iDet);
							value_type_r g_outer = this->input_multislice->det_cir.g_outer(iDet);
							int iscan = this->input_multislice->iscan;
							output_multislice.image_tot[ithk].image[iDet][iscan] += w_i*sum_square_over_Det(*(this->stream), this->input_multislice->grid, g_inner, g_outer, *psi_zt);
						}

						if(this->input_multislice->coherent_contribution)
						{
							multem::add_scale_to_host(output_multislice.stream, this->input_multislice->grid, w_i, *psi_zt, output_multislice.psi_coh[ithk], &psi_zh);
						}
					}
					else
					{
						if(this->input_multislice->coherent_contribution)
						{
							multem::add_scale_m2psi_psi_to_host(output_multislice.stream, this->input_multislice->grid, w_i, *psi_zt, output_multislice.m2psi_tot[ithk], output_multislice.psi_coh[ithk], &psi_zh);
						}
						else
						{
							multem::assign_square(*psi_zt, m2psi_z);
							multem::add_scale_to_host(output_multislice.stream, this->input_multislice->grid, w_i, m2psi_z, output_multislice.m2psi_tot[ithk], &m2psi_zh);
						}
					}
				}
			}

			template<class TOutput_multislice>
			void set_m2psi_coh(TOutput_multislice &output_multislice)
			{
				if(!this->input_multislice->coherent_contribution || this->input_multislice->is_EWFS_EWRS() || this->input_multislice->is_EWSFS_EWSRS())
				{
					return;
				}

				if(this->input_multislice->is_HRTEM_HCI_EFTEM())
				{
					for(auto ithk=0; ithk < this->thickness.size(); ithk++)
					{
						multem::assign(output_multislice.psi_coh[ithk], psi_z, &psi_zh);
						microscope_effects.apply(psi_z, m2psi_z);
						multem::copy_to_host(output_multislice.stream, this->input_multislice->grid, m2psi_z, output_multislice.m2psi_coh[ithk], &m2psi_zh);
					}
				}
				else if(this->input_multislice->is_STEM())
				{
					for(auto ithk=0; ithk < this->thickness.size(); ithk++)
					{
						multem::assign(output_multislice.psi_coh[ithk], psi_z, &psi_zh);
						for(auto iDet=0; iDet<this->input_multislice->det_cir.size(); iDet++)
						{
							value_type_r g_inner = this->input_multislice->det_cir.g_inner(iDet);
							value_type_r g_outer = this->input_multislice->det_cir.g_outer(iDet);
							int iscan = this->input_multislice->iscan;
							output_multislice.image_coh[ithk].image[iDet][iscan] = sum_square_over_Det(*(this->stream), this->input_multislice->grid, g_inner, g_outer, psi_z);
						}
					}
				}
				else
				{
					for(auto ithk=0; ithk < this->thickness.size(); ithk++)
					{
						multem::assign_square(output_multislice.psi_coh[ithk], output_multislice.m2psi_coh[ithk]);
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

			template<class TOutput_multislice>
			void psi(value_type_r w_i, Vector<value_type_c, dev> &psi_z, TOutput_multislice &output_multislice)
			{
				value_type_r gx_0 = this->input_multislice->gx_0();
				value_type_r gy_0 = this->input_multislice->gy_0();

				for(auto islice=0; islice<this->slice.size(); islice++)
				{
					psi_slice(gx_0, gy_0, islice, psi_z);
					set_m2psi_tot_psi_coh(psi_z, gx_0, gy_0, islice, w_i, output_multislice);
				}
			}

			template<class TOutput_multislice>
			void psi(int islice_0, int islice_e, value_type_r w_i, Vector<value_type_c, dev> &trans, TOutput_multislice &output_multislice)
			{
				int ithk = this->slice.ithk[islice_e];
				if(0 <= ithk)
				{
					value_type_r gx_0 = this->input_multislice->gx_0();
					value_type_r gy_0 = this->input_multislice->gy_0();

					if(this->input_multislice->eels_fr.is_Single_Channelling())
					{
						value_type_r dz = this->dz_m(islice_0, islice_e);
						prog.propagate(eS_Reciprocal, gx_0, gy_0, dz, psi_z);
					}
					else if(this->input_multislice->eels_fr.is_Double_Channelling_FOMS())
					{
						value_type_r dz = this->dz_m(islice_0, islice_e);
						multem::multiply(trans, psi_z);
						prog.propagate(eS_Reciprocal, gx_0, gy_0, dz, psi_z);
					}
					else if(this->input_multislice->eels_fr.is_Double_Channelling_SOMS())
					{
						value_type_r dz = 0.5*this->dz_m(islice_0, islice_e);
						this->prog.propagate(eS_Real, gx_0, gy_0, dz, psi_z);
						multem::multiply(trans, psi_z);
						prog.propagate(eS_Reciprocal, gx_0, gy_0, dz, psi_z);
					}
					else if(this->input_multislice->eels_fr.is_Double_Channelling())
					{
						for(auto islice=islice_0; islice<=islice_e; islice++)
						{
							psi_slice(gx_0, gy_0, islice, psi_z);
						}
						phase_multiplication(gx_0, gy_0, psi_z);
						prog.propagate(eS_Reciprocal, gx_0, gy_0, this->thickness.z_back_prop[ithk], psi_z);
					}

					if(this->input_multislice->is_EELS())
					{
						int iscan = this->input_multislice->iscan;
						output_multislice.image_tot[ithk].image[0][iscan] += w_i*sum_square_over_Det(*(this->stream), this->input_multislice->grid, 0, this->input_multislice->eels_fr.g_collection, psi_z);
					}
					else
					{
						multem::bandwidth_limit(*(this->stream), this->input_multislice->grid, 0, this->input_multislice->eels_fr.g_collection, 1.0, psi_z);
						this->fft2->inverse(psi_z);
						multem::assign_square(psi_z, m2psi_z);
						multem::add_scale_to_host(output_multislice.stream, this->input_multislice->grid, w_i, m2psi_z, output_multislice.m2psi_tot[ithk], &m2psi_zh);
					}
				}
			}

			Propagator<value_type_r, dev> prog;

			Vector<value_type_c, dev> psi_z;
			Vector<value_type_r, dev> m2psi_z;

			Vector<value_type_c, e_host> psi_zh;
			Vector<value_type_r, e_host> m2psi_zh;

			Microscope_Effects<value_type_r, dev> microscope_effects;
		private:
			Vector<value_type_c, dev> exp_x;
			Vector<value_type_c, dev> exp_y;
	};

} // namespace multem

#endif