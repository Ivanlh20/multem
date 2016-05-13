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
	template<class T, eDevice dev>
	class Microscope_Effects
	{
		public:
			using value_type_r = T;
			using value_type_c = complex<T>;

			Microscope_Effects(): input_multislice(nullptr), stream(nullptr), fft2(nullptr){}			
			
			void set_input_data(Input_Multislice<value_type_r> *input_multislice_i, Stream<dev> *stream_i, FFT2<value_type_r, dev> *fft2_i)
			{
				input_multislice = input_multislice_i;
				stream = stream_i;
				fft2 = fft2_i;

				psi.resize(input_multislice->grid.nxy());

				if((input_multislice->illumination_model == eIM_Coherent)||(input_multislice->illumination_model == eIM_Partial_Coherent))
				{
					return;
				}

				Stream<e_host> stream_host(input_multislice->nstream);
				/*********************Temporal quadrature**********************/
				Quadrature quadrature;
				quadrature.get(8, input_multislice->obj_lens.nsf, qt); // 8: int_-infty^infty f(x) Exp[-x^2] dx
				mt::scale(stream_host, 1.0/c_Pii2, qt.w);

				/*********************Spatial quadrature**********************/
				qs.resize((2*input_multislice->obj_lens.ngxs+1)*(2*input_multislice->obj_lens.ngys+1));
				int nqs = 0; 
				value_type_r sum_w = 0;
				value_type_r alpha = 0.5/pow(input_multislice->obj_lens.sggs, 2);

				for(auto ix =-input_multislice->obj_lens.ngxs; ix <= input_multislice->obj_lens.ngxs; ix++)
				 {
					 for(auto iy =-input_multislice->obj_lens.ngys; iy <= input_multislice->obj_lens.ngys; iy++)
					 {
						 value_type_r gxs = input_multislice->obj_lens.gxs(ix);
						 value_type_r gys = input_multislice->obj_lens.gys(iy);
						 value_type_r g2s = gxs*gxs + gys*gys;
						if(g2s < input_multislice->obj_lens.g2_maxs)
						{
							qs.x[nqs] = gxs;
							qs.y[nqs] = gys;
							sum_w += qs.w[nqs] = exp(-alpha*g2s);
							nqs++;
						}
					}
				 }
				qs.resize(nqs);
				mt::scale(stream_host, 1.0/sum_w, qs.w);
			}

			void apply(Vector<value_type_c, dev> &fpsi, Vector<value_type_r, dev> &m2psi_tot)
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
						TCC_TEM(input_multislice->temporal_spatial_incoh, fpsi, m2psi_tot);
					}
					break;
				}
			}

			template<class TOutput_multislice>
			void apply(TOutput_multislice &output_multislice)
			{
				Vector<value_type_c, dev> psi(input_multislice->iw_psi.begin(), input_multislice->iw_psi.end());
				mt::fft2_shift(*stream, input_multislice->grid, psi);
				fft2->forward(psi);
				mt::scale(*stream, input_multislice->grid.inxy, psi);

				Vector<value_type_r, dev> m2psi_tot(input_multislice->grid.nxy());
				apply(psi, m2psi_tot);
				mt::copy_to_host(output_multislice.stream, m2psi_tot, output_multislice.m2psi_tot[0]);
				output_multislice.shift();
				output_multislice.clear_temporal_data();
			}

		private:
			void CTF_TEM(const eTemporal_Spatial_Incoh &temporal_spatial_incoh, Vector<value_type_c, dev> &fpsi, Vector<value_type_r, dev> &m2psi_tot)
			{
				mt::apply_CTF(*stream, input_multislice->grid, input_multislice->obj_lens, 0, 0, fpsi, psi);
				fft2->inverse(psi);
				mt::square(*stream, psi, m2psi_tot);
			}

			void PCTF_LI_WPO_TEM(const eTemporal_Spatial_Incoh &temporal_spatial_incoh, Vector<value_type_c, dev> &fpsi, Vector<value_type_r, dev> &m2psi_tot)
			{
				value_type_r sf = input_multislice->obj_lens.sf;
				value_type_r beta = input_multislice->obj_lens.beta;

				switch(temporal_spatial_incoh)
				{
					case eTSI_Temporal:	// Temporal
					{
						input_multislice->obj_lens.beta = 0;
					}
					break;
					case eTSI_Spatial:	// Spatial
					{
						input_multislice->obj_lens.sf = 0;
					}
					break;
				}

				mt::apply_PCTF(*stream, input_multislice->grid, input_multislice->obj_lens, fpsi, psi);
				fft2->inverse(psi);
				mt::square(*stream, psi, m2psi_tot);

				input_multislice->obj_lens.sf = sf;
				input_multislice->obj_lens.beta = beta;
			}

			void TCC_TEM(const eTemporal_Spatial_Incoh &temporal_spatial_incoh, Vector<value_type_c, dev> &fpsi, Vector<value_type_r, dev> &m2psi_tot)
			{
				value_type_r f_0 = input_multislice->obj_lens.f;

				fill(*stream, m2psi_tot, 0.0);
				switch(temporal_spatial_incoh)
				{
					case 1:	// Temporal and Spatial
					{
						for(auto i = 0; i<qs.size(); i++)
						{
							for(auto j = 0; j<qt.size(); j++)
							{
								auto f = input_multislice->obj_lens.sf*qt.x[j]+f_0;
								input_multislice->obj_lens.set_defocus(f); 
								
								mt::apply_CTF(*stream, input_multislice->grid, input_multislice->obj_lens, qs.x[i], qs.y[i], fpsi, psi);
								fft2->inverse(psi);
								mt::add_square_scale(*stream, qs.w[i]*qt.w[j], psi, m2psi_tot);
							}
						}
					}
					break;
					case 2:	// Temporal
					{
						for(auto j = 0; j<qt.size(); j++)
						{
							auto f = input_multislice->obj_lens.sf*qt.x[j]+f_0;
							input_multislice->obj_lens.set_defocus(f); 

							mt::apply_CTF(*stream, input_multislice->grid, input_multislice->obj_lens, 0.0, 0.0, fpsi, psi);
							fft2->inverse(psi);
							mt::add_square_scale(*stream, qt.w[j], psi, m2psi_tot);
						}
					}
					break;
					case 3:	// Spatial
					{
						for(auto i = 0; i<qs.size(); i++)
						{
							mt::apply_CTF(*stream, input_multislice->grid, input_multislice->obj_lens, qs.x[i], qs.y[i], fpsi, psi);
							fft2->inverse(psi);
							mt::add_square_scale(*stream, qs.w[i], psi, m2psi_tot);
						}
					}
				}

				input_multislice->obj_lens.set_defocus(f_0);
			}
			
			Input_Multislice<value_type_r> *input_multislice;
			Stream<dev> *stream;
			FFT2<value_type_r, dev> *fft2;

			Vector<value_type_c, dev> psi;

			Q1<value_type_r, e_host> qt;
			Q2<value_type_r, e_host> qs;
	};

} // namespace mt

#endif