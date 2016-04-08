/*
 * This file is part of MULTEM.
 * Copyright 2016 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H

#include "types.cuh"
#include "fft2.cuh"
#include "input_multislice.cuh"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "transmission_function.cuh"
#include "incident_wave.cuh"
#include "propagator.cuh"
#include "microscope_effects.cuh"

namespace multem
{
	template<class T, eDevice dev>
	class Wave_Function: public Transmission_Function<T, dev>
	{
		public:
			using value_type_r = T;
			using value_type_c = complex<T>;
			using size_type = std::size_t;

			static const eDevice device = dev;

			void set_input_data(Input_Multislice<value_type_r> *input_multislice_i, Stream<dev> *stream_i, FFT2<value_type_r, dev> *fft2_i)
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

				if(input_multislice_i->is_STEM())
				{
					detector.assign(input_multislice_i->detector);
					if(input_multislice_i->is_detector_matrix())
					{
						for(auto i = 0; i<detector.size(); i++)
						{
							multem::fft2_shift(*stream_i, input_multislice_i->grid, detector.fR[i]);
						}
					}
				}

				incident_wave.set_input_data(input_multislice_i, stream_i, fft2_i);

				prog.set_input_data(input_multislice_i, stream_i, fft2_i);

				if(input_multislice_i->is_ISTEM_CBEI_HRTEM_HCI_EFTEM())
				{
					microscope_effects.set_input_data(input_multislice_i, stream_i, fft2_i);
				} 

				Transmission_Function<T, dev>::set_input_data(input_multislice_i, stream_i, fft2_i);
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

				multem::phase_components(*(this->stream), this->input_multislice->grid, gxu, gyu, exp_x, exp_y);
				multem::phase_multiplication(*(this->stream), this->input_multislice->grid, exp_x, exp_y, psi_i, psi_o);
			}	

			void phase_multiplication(const value_type_r &gxu, const value_type_r &gyu, Vector<value_type_c, dev> &psi_io)
			{
				phase_multiplication(gxu, gyu, psi_io, psi_io);
			}

			Vector<value_type_c, dev>* get_psi(const eSpace &space, const value_type_r &gxu, const value_type_r &gyu, 
			value_type_r z, Vector<value_type_c, dev> &psi_i)
			{
				Vector<value_type_c, dev> *psi_o = &(this->trans_0);
				phase_multiplication(gxu, gyu, psi_i, *psi_o);
				prog.propagate(space, gxu, gyu, z, *psi_o);

				return psi_o;
			}

			value_type_r integrated_intensity_over_det(value_type_r w_i, const int &iDet, Vector<value_type_c, dev> &psi_z)
			{
				value_type_r int_val = 0;
				switch (detector.type)
				{
					case multem::eDT_Circular:
					{
						auto g_inner = detector.g_inner[iDet];
						auto g_outer = detector.g_outer[iDet];
							
						int_val = w_i*multem::sum_square_over_Det(*(this->stream), this->input_multislice->grid, g_inner, g_outer, psi_z);
					}
					break;
					case multem::eDT_Radial:
					{
						int_val = 0;
					}
					break;
					case multem::eDT_Matrix:
					{
						int_val = w_i*multem::sum_square_over_Det(*(this->stream), this->input_multislice->grid, detector.fR[iDet], psi_z);
					}
					break;
				}

				return int_val;
			}

			template<class TOutput_multislice>
			void set_m2psi_tot_psi_coh(Vector<value_type_c, dev> &psi_z_i, const value_type_r &gxu, const value_type_r &gyu, 
			const int &islice, const value_type_r &w_i, TOutput_multislice &output_multislice)
			{
				int ithk = this->slice.ithk[islice];
				if(0 <= ithk)
				{
					auto *psi_z_o = get_psi(this->input_multislice->get_simulation_space(), gxu, gyu, this->thickness.z_back_prop[ithk], psi_z_i);

					if(this->input_multislice->is_STEM())
					{
						for(auto iDet = 0; iDet<detector.size(); iDet++)
						{
							auto iscan = this->input_multislice->iscan;
							output_multislice.image_tot[ithk].image[iDet][iscan] += integrated_intensity_over_det(w_i, iDet, *psi_z_o);
						}

						if(this->input_multislice->coherent_contribution)
						{
							multem::add_scale_to_host(output_multislice.stream, w_i, *psi_z_o, output_multislice.psi_coh[ithk], &psi_zh);
						}
					}
					else if(this->input_multislice->is_EWFS_EWRS_SC())
					{
						multem::copy_to_host(output_multislice.stream, *psi_z_o, output_multislice.psi_coh[ithk], &psi_zh);
					}
					else if(this->input_multislice->is_ISTEM_CBEI_HRTEM_HCI_EFTEM())
					{
						microscope_effects.apply(*psi_z_o, m2psi_z);
						multem::add_scale_to_host(output_multislice.stream, w_i, m2psi_z, output_multislice.m2psi_tot[ithk], &m2psi_zh);

						if(this->input_multislice->coherent_contribution)
						{
							multem::add_scale_to_host(output_multislice.stream, w_i, *psi_z_o, output_multislice.psi_coh[ithk], &psi_zh);
						}
					}
					else
					{
						if(this->input_multislice->coherent_contribution)
						{
							multem::add_scale_m2psi_psi_to_host(output_multislice.stream, w_i, *psi_z_o, output_multislice.m2psi_tot[ithk], output_multislice.psi_coh[ithk], &psi_zh);
						}
						else
						{
							multem::square(*(this->stream), *psi_z_o, m2psi_z);
							multem::add_scale_to_host(output_multislice.stream, w_i, m2psi_z, output_multislice.m2psi_tot[ithk], &m2psi_zh);
						}
					}
				}
			}

			template<class TOutput_multislice>
			void set_m2psi_coh(TOutput_multislice &output_multislice)
			{
				if(!this->input_multislice->coherent_contribution || this->input_multislice->is_EWFS_EWRS())
				{
					return;
				}

				if(this->input_multislice->is_STEM())
				{
					for(auto ithk = 0; ithk < this->thickness.size(); ithk++)
					{
						multem::assign(output_multislice.stream, output_multislice.psi_coh[ithk], psi_z, &psi_zh);
						for(auto iDet = 0; iDet<detector.size(); iDet++)
						{
							auto iscan = this->input_multislice->iscan;
							output_multislice.image_coh[ithk].image[iDet][iscan] += integrated_intensity_over_det(1, iDet, psi_z);
						}
					}
				}
				else if(this->input_multislice->is_ISTEM_CBEI_HRTEM_HCI_EFTEM())
				{
					for(auto ithk = 0; ithk < this->thickness.size(); ithk++)
					{
						multem::assign(output_multislice.stream, output_multislice.psi_coh[ithk], psi_z, &psi_zh);
						microscope_effects.apply(psi_z, m2psi_z);
						multem::copy_to_host(output_multislice.stream, m2psi_z, output_multislice.m2psi_coh[ithk], &m2psi_zh);
					}
				}
				else
				{
					for(auto ithk =0; ithk < this->thickness.size(); ithk++)
					{
						multem::square(output_multislice.stream, output_multislice.psi_coh[ithk], output_multislice.m2psi_coh[ithk]);
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

				for(auto islice = 0; islice<this->slice.size(); islice++)
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
						prog.propagate(eS_Real, gx_0, gy_0, dz, psi_z);
					}
					else if(this->input_multislice->eels_fr.is_Mixed_Channelling())
					{
						value_type_r dz = 0.5*this->dz_m(islice_0, islice_e);
						this->prog.propagate(eS_Real, gx_0, gy_0, dz, psi_z);
						multem::multiply(*(this->stream), trans, psi_z);
						prog.propagate(eS_Real, gx_0, gy_0, dz, psi_z);
					}
					else if(this->input_multislice->eels_fr.is_Double_Channelling())
					{
						for(auto islice = islice_0; islice<= islice_e; islice++)
						{
							psi_slice(gx_0, gy_0, islice, psi_z);
						}
					}

					phase_multiplication(gx_0, gy_0, psi_z);
					prog.propagate(eS_Reciprocal, gx_0, gy_0, this->thickness.z_back_prop[ithk], psi_z);

					if(this->input_multislice->is_EELS())
					{
						int iscan = this->input_multislice->iscan;
						output_multislice.image_tot[ithk].image[0][iscan] += w_i*multem::sum_square_over_Det(*(this->stream), this->input_multislice->grid, 0, this->input_multislice->eels_fr.g_collection, psi_z);
					}
					else
					{
						multem::hard_aperture(*(this->stream), this->input_multislice->grid, this->input_multislice->eels_fr.g_collection, 1.0, psi_z);
						microscope_effects.apply(psi_z, m2psi_z);
						multem::add_scale_to_host(output_multislice.stream, w_i, m2psi_z, output_multislice.m2psi_tot[ithk], &m2psi_zh);
					}
				}
			}

			void set_incident_wave(Vector<value_type_c, dev> &psi)
			{
				auto z_init = this->slice.z_m(0);
				incident_wave(psi, z_init);
			}
			Propagator<value_type_r, dev> prog;

			Vector<value_type_c, dev> psi_z;
			Vector<value_type_r, dev> m2psi_z;

			Vector<value_type_c, e_host> psi_zh;
			Vector<value_type_r, e_host> m2psi_zh;

			Detector<T, dev> detector; 	
			Microscope_Effects<value_type_r, dev> microscope_effects;
		private:
			Vector<value_type_c, dev> exp_x;
			Vector<value_type_c, dev> exp_y;
			Incident_Wave<value_type_r, dev> incident_wave;
	};

} // namespace multem

#endif