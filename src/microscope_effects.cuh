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

#ifndef MICROSCOPE_EFFECTS_H
#define MICROSCOPE_EFFECTS_H

#include "math.cuh"
#include "types.hpp"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "quadrature.hpp"

namespace multem
{
	template<class T, eDevice dev>
	class Microscope_Effects
	{
		public:
			using value_type_r = typename T;
			using value_type_c = typename complex<T>;

			Microscope_Effects():input_multislice(nullptr), stream(nullptr), fft2(nullptr) {}			
			
			void set_input_data(Input_Multislice<value_type_r, dev> *input_multislice_io, Stream<value_type_r, dev> *stream_i, FFT2<value_type_r, dev> *fft2_i)
			{
				input_multislice = input_multislice_io;
				stream = stream_i;
				fft2 = fft2_i;

				psi.resize(input_multislice->grid.nxy());

				if(input_multislice->microscope_effect == eME_Partial_Coherent)
				{
					return;
				}

				/*********************Temporal quadrature**********************/
				Quadrature quadrature;
				quadrature.get(8, input_multislice->lens.nsf, qt); // 8: int_-infty^infty f(x) Exp[-x^2] dx
				multem::scale(qt.w, 1.0/c_Pii2);

				/*********************Spatial quadrature**********************/
				qs.resize((2*input_multislice->lens.ngxs+1)*(2*input_multislice->lens.ngys+1));
				int nqs = 0; 
				value_type_r sum_w = 0;
				value_type_r alpha = 0.5/pow(input_multislice->lens.sggs, 2);

				for(auto ix=-input_multislice->lens.ngxs; ix <= input_multislice->lens.ngxs; ix++)
				 {
					 for(auto iy=-input_multislice->lens.ngys; iy <= input_multislice->lens.ngys; iy++)
					 {
						 value_type_r gxs = input_multislice->lens.gxs(ix);
						 value_type_r gys = input_multislice->lens.gys(iy);
						 value_type_r g2s = gxs*gxs + gys*gys;
						if(g2s < input_multislice->lens.g2_maxs)
						{
							qs.x[nqs] = gxs;
							qs.y[nqs] = gys;
							sum_w += qs.w[nqs] = exp(-alpha*g2s);
							nqs++;
						}
					}
				 }
				qs.resize(nqs);
				multem::scale(qs.w, 1.0/sum_w);
			}

			void apply(Vector<value_type_c, dev> &fPsi, Vector<value_type_r, dev> &m2psi_tot)
			{
				switch(input_multislice->microscope_effect)
				{
					case eME_Partial_Coherent:
					{
						PC_LI_WPO_TEM(input_multislice->spatial_temporal_effect, fPsi, m2psi_tot);
					}
					break;
					case eME_Transmission_Cross_Coefficient:
					{
						TCC_TEM(input_multislice->spatial_temporal_effect, fPsi, m2psi_tot);
					}
					break;
				}
			}

		private:
			void PC_LI_WPO_TEM(const eSpatial_Temporal_Effect &spatial_temporal_effect, Vector<value_type_c, dev> &fPsi, Vector<value_type_r, dev> &m2psi_tot)
			{
				value_type_r sf = input_multislice->lens.sf;
				value_type_r beta = input_multislice->lens.beta;

				switch(spatial_temporal_effect)
				{
					case eSTE_Temporal:	// Temporal
					{
						input_multislice->lens.beta = 0;
					}
					break;
					case eSTE_Spatial:	// Spatial
					{
						input_multislice->lens.sf = 0;
					}
					break;
				}

				multem::apply_PCTF(input_multislice->grid, input_multislice->lens, fPsi, psi);
				fft2->inverse(psi);
				assign_square(psi, m2psi_tot);

				input_multislice->lens.sf = sf;
				input_multislice->lens.beta = beta;
			}

			void TCC_TEM(const eSpatial_Temporal_Effect &spatial_temporal_effect, Vector<value_type_c, dev> &fPsi, Vector<value_type_r, dev> &m2psi_tot)
			{
				value_type_r f_0 = input_multislice->lens.f;
				value_type_r cf_0 = input_multislice->lens.cf;

				fill(m2psi_tot, 0.0);
				switch(spatial_temporal_effect)
				{
					case 1:	// Temporal and Spatial
					{
						for(auto i=0; i<qs.size(); i++)
						{
							for(auto j=0; j<qt.size(); j++)
							{
								input_multislice->lens.f = input_multislice->lens.sf*qt.x[j]+f_0; 
								input_multislice->lens.cf = c_Pi*input_multislice->lens.lambda*input_multislice->lens.f;
								
								multem::apply_CTF(input_multislice->grid, input_multislice->lens, qs.x[i], qs.y[i], fPsi, psi);
								fft2->inverse(psi);
								multem::add_square_scale(qs.w[i]*qt.w[j], psi, m2psi_tot);
							}
						}
					}
					break;
					case 2:	// Temporal
					{
						for(auto j=0; j<qt.size(); j++)
						{
							input_multislice->lens.f = input_multislice->lens.sf*qt.x[j]+f_0; 
							input_multislice->lens.cf = c_Pi*input_multislice->lens.lambda*input_multislice->lens.f;
								
							multem::apply_CTF(input_multislice->grid, input_multislice->lens, 0.0, 0.0, fPsi, psi);
							fft2->inverse(psi);
							multem::add_square_scale(qt.w[j], psi, m2psi_tot);
						}
					}
					break;
					case 3:	// Spatial
					{
						for(auto i=0; i<qs.size(); i++)
						{
							multem::apply_CTF(input_multislice->grid, input_multislice->lens, qs.x[i], qs.y[i], fPsi, psi);
							fft2->inverse(psi);
							multem::add_square_scale(qs.w[i], psi, m2psi_tot);
						}
					}
				}

				input_multislice->lens.f = f_0;
				input_multislice->lens.cf = cf_0;
			}
			
			Input_Multislice<value_type_r, dev> *input_multislice;
			Stream<value_type_r, dev> *stream;
			FFT2<value_type_r, dev> *fft2;

			Vector<value_type_c, dev> psi;

			Q1<value_type_r, e_Host> qt;
			Q2<value_type_r, e_Host> qs;
	};

} //namespace multem

#endif