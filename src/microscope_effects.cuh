/*
 * This file is part of MULTEM.
 * Copyright 2017 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef MICROSCOPE_EFFECTS_H
#define MICROSCOPE_EFFECTS_H

#include "math.cuh"
#include "types.cuh"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "quadrature.hpp"

namespace mt
{
	template <class T, eDevice dev>
	class Microscope_Effects
	{
		public:
			using T_r = T;
			using T_c = complex<T>;

			Microscope_Effects(): input_multislice(nullptr), stream(nullptr), fft_2d(nullptr){}			
			
			void set_input_data(Input_Multislice<T_r> *input_multislice_i, Stream<dev> *stream_i, FFT<T_r, dev> *fft2_i)
			{
				input_multislice = input_multislice_i;
				stream = stream_i;
				fft_2d = fft2_i;

				psi.resize(input_multislice->grid_2d.nxy());

				if((input_multislice->illumination_model == eIM_Coherent)||(input_multislice->illumination_model == eIM_Partial_Coherent))
				{
					return;
				}

				// Load quadratures
				obj_lens_temporal_spatial_quadratures(input_multislice->obj_lens, qt, qs);
			}

			void operator()(Vector<T_c, dev> &fpsi, Vector<T_r, dev> &m2psi_tot)
			{
				switch(input_multislice->illumination_model)
				{
					case eIM_Coherent:
					{
						CTF_TEM(input_multislice->temporal_spatial_incoh, fpsi, m2psi_tot);
					}
					break;
					case eIM_Partial_Coherent:
					{
						PCTF_LI_WPO_TEM(input_multislice->temporal_spatial_incoh, fpsi, m2psi_tot);
					}
					break;
					case eIM_Trans_Cross_Coef:
					{

					}
					break;
					case eIM_Full_Integration:
					{
						num_int_TEM(input_multislice->temporal_spatial_incoh, fpsi, m2psi_tot);
					}
					break;
				}
			}

			template <class TOutput_multislice>
			void operator()(TOutput_multislice &output_multislice)
			{
				Vector<T_c, dev> psi(input_multislice->iw_psi.begin(), input_multislice->iw_psi.end());
				mt::fft2_shift(*stream, input_multislice->grid_2d, psi);
				fft_2d->forward(psi);
				mt::scale(*stream, input_multislice->grid_2d.inxy(), psi);

				Vector<T_r, dev> m2psi_tot(input_multislice->grid_2d.nxy());
				this->operator()(psi, m2psi_tot);
				mt::copy_to_host(output_multislice.stream, m2psi_tot, output_multislice.m2psi_tot[0]);
			}

		private:
			void CTF_TEM(const eTemporal_Spatial_Incoh &temporal_spatial_incoh, Vector<T_c, dev> &fpsi, Vector<T_r, dev> &m2psi_tot)
			{
				mt::apply_CTF(*stream, input_multislice->grid_2d, input_multislice->obj_lens, 0, 0, fpsi, psi);
				fft_2d->inverse(psi);
				mt::square(*stream, psi, m2psi_tot);
			}

			void PCTF_LI_WPO_TEM(const eTemporal_Spatial_Incoh &temporal_spatial_incoh, Vector<T_c, dev> &fpsi, Vector<T_r, dev> &m2psi_tot)
			{
				T_r dsf_sigma = input_multislice->obj_lens.dsf_sigma;
				T_r ssf_sigma = input_multislice->obj_lens.ssf_sigma;

				switch(temporal_spatial_incoh)
				{
					case eTSI_Temporal:	// Temporal
					{
						input_multislice->obj_lens.set_ssf_sigma(0);
					}
					break;
					case eTSI_Spatial:	// Spatial
					{
						input_multislice->obj_lens.set_dsf_sigma(0);
					}
					break;
				}

				mt::apply_PCTF(*stream, input_multislice->grid_2d, input_multislice->obj_lens, fpsi, psi);
				fft_2d->inverse(psi);
				mt::square(*stream, psi, m2psi_tot);

				input_multislice->obj_lens.set_dsf_sigma(dsf_sigma);
				input_multislice->obj_lens.set_ssf_sigma(ssf_sigma);
			}

			void num_int_TEM(const eTemporal_Spatial_Incoh &temporal_spatial_incoh, Vector<T_c, dev> &fpsi, Vector<T_r, dev> &m2psi_tot)
			{
				T_r c_10_0 = input_multislice->obj_lens.c_10;

				fill(*stream, m2psi_tot, 0.0);
				switch(temporal_spatial_incoh)
				{
					case 1:	// Temporal and Spatial
					{
						for(auto i = 0; i<qs.size(); i++)
						{
							for(auto j = 0; j<qt.size(); j++)
							{
								auto c_10 = input_multislice->obj_lens.dsf_iehwgd*qt.x[j]+c_10_0;
								input_multislice->obj_lens.set_defocus(c_10); 
								
								mt::apply_CTF(*stream, input_multislice->grid_2d, input_multislice->obj_lens, qs.x[i], qs.y[i], fpsi, psi);
								fft_2d->inverse(psi);
								mt::add_scale_square(*stream, qs.w[i]*qt.w[j], psi, m2psi_tot);
							}
						}
					}
					break;
					case 2:	// Temporal
					{
						for(auto j = 0; j<qt.size(); j++)
						{
							auto c_10 = input_multislice->obj_lens.dsf_iehwgd*qt.x[j]+c_10_0;
							input_multislice->obj_lens.set_defocus(c_10); 

							mt::apply_CTF(*stream, input_multislice->grid_2d, input_multislice->obj_lens, 0.0, 0.0, fpsi, psi);
							fft_2d->inverse(psi);
							mt::add_scale_square(*stream, qt.w[j], psi, m2psi_tot);
						}
					}
					break;
					case 3:	// Spatial
					{
						for(auto i = 0; i<qs.size(); i++)
						{
							mt::apply_CTF(*stream, input_multislice->grid_2d, input_multislice->obj_lens, qs.x[i], qs.y[i], fpsi, psi);
							fft_2d->inverse(psi);
							mt::add_scale_square(*stream, qs.w[i], psi, m2psi_tot);
						}
					}
				}

				input_multislice->obj_lens.set_defocus(c_10_0);
			}
			
			Input_Multislice<T_r> *input_multislice;
			Stream<dev> *stream;
			FFT<T_r, dev> *fft_2d;

			Vector<T_c, dev> psi;

			Q1<T_r, e_host> qt;
			Q2<T_r, e_host> qs;
	};

} // namespace mt

#endif