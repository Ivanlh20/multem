/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * Multem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version of the License, or
 * (at your option) any later version.
 *
 * Multem is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Multem. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef ENERGY_LOSS_H
	#define ENERGY_LOSS_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include "const_enum_mt.cuh"
	#include "math.cuh"
	#include "type_traits_gen.cuh"
	#include "cgpu_stream.cuh"
	#include "cgpu_fft.cuh"
	#include "cgpu_fcns_gen.cuh"
	// #include "in_classes.cuh"
	// #include "cpu_fcns.hpp"
	// #include "gpu_fcns.cuh"
	// #include "cgpu_fcns.cuh"

	namespace mt
	{
		template <class T>
		class EELS
		{
		public:
			using value_type = T;

			eSpace space;			// real or reciprocal space

			T E_0;					// incident energy
			dt_int32 m_sel;			// selection rule
			T coll_angle;			// collection angle(rad)
			eChan_Typ chan_type;	// channelling type
			T g_coll;				// collection angle(recirpocal Angstroms)

			dt_int32 Z;				// atomic type
			T x;					// x position
			T y;					// y position
			T occ;					// y occupancy
			T E_loss;				// energy loss
			T ge;					// relativistic corrected effective scattering momentum transfer
			T ge2;					// ge square
			T gc;					// ge cut-off at the Bethe Ridge
			T gc2;					// gc square
			T lrtz_factor;			// lorentz factor = sum 1/(g^2+ge^2)

			/************************************* constructors ************************************/
			EELS(): space(esp_real), E_0(0), m_sel(0), coll_angle(0), chan_type(eCT_double_chan), 
			g_coll(0), lrtz_factor(1), Z(0), x(0), y(0), occ(0), E_loss(0), ge(0), ge2(0), gc(0), gc2(0){}

			/* copy constructor */
			EELS(const EELS<T>& eels)
			{
				*this = eels;
			}

			/* converting constructor */
			template <class U> 
			EELS(const EELS<U>& eels)
			{
				*this = eels;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			CGPU_EXEC
			EELS<T>& operator=(const EELS<T>& eels)
			{
				if (this != &eels)
				{
					space = eels.space;

					E_0 = eels.E_0;
					m_sel = eels.m_sel;
					coll_angle = eels.coll_angle;
					chan_type = eels.chan_type;
					g_coll = eels.g_coll;

					Z = eels.Z;
					x = eels.x;
					y = eels.y;
					occ = eels.occ;
					E_loss = eels.E_loss;
					ge = eels.ge;
					ge2 = eels.ge2;
					gc = eels.gc;
					gc2 = eels.gc2;
					lrtz_factor = eels.lrtz_factor;
				}

				return *this;
			}

			/* converting assignment operator */
			template <class U>
			CGPU_EXEC
			EELS<T>& operator=(const EELS<U>& eels)
			{
				space = eels.space;

				E_0 = T(eels.E_0);
				m_sel = eels.m_sel;
				coll_angle = T(eels.coll_angle);
				chan_type = eels.chan_type;
				g_coll = T(eels.g_coll);

				Z = eels.Z;
				x = T(eels.x);
				y = T(eels.y);
				occ = T(eels.occ);
				E_loss = T(eels.E_loss);
				ge = T(eels.ge);
				ge2 = T(eels.ge2);
				gc = T(eels.gc);
				gc2 = T(eels.gc2);
				lrtz_factor = T(eels.lrtz_factor);
			
				return *this;
			}

			template <class U> 
			void assign(const EELS<U>& eels)
			{ 
				*this = eels;
			}

			/***************************************************************************************/
			void set_in_data(const eSpace& space, const T& E_0, dt_int32 m_sel, const T& coll_angle, eChan_Typ chan_type, const dt_int32& Z, const T& E_loss)
			{
				this->space = space;
				this->E_0 = E_0;
				this->m_sel = m_sel;
				set_coll_angle(coll_angle);
				this->chan_type = chan_type;

				this->Z = Z;
				set_e_loss(E_loss);
			}

			void set_coll_angle(const T& coll_angle)
			{
				this->coll_angle = coll_angle;
				g_coll = coll_angle/fcn_lambda(E_0);
			}

			void set_e_loss(const T& E_loss)
			{
				const auto gamma = fcn_gamma(E_0);
				const auto lambda = fcn_lambda(E_0);

				auto theta = theta_efect(E_loss, E_0);
				ge = theta/(lambda*gamma*gamma);
				ge2 = pow(gamma*ge, 2);
				gc = sqrt(T(2)*theta)/lambda;
				gc2 = pow(gc, 2);
			}

			// effective scattering angle
			T theta_efect(const T& E_loss, const T& E_0)
			{
				const auto emass = T(510.99906);
				const auto x = (emass + E_0)/(T(2)*emass + E_0);
				return E_loss*x/E_0;
			}

			dt_bool is_Single_Chan() const
			{
				return chan_type == ect_single_chan;
			}

			dt_bool is_Mixed_Chan() const
			{
				return chan_type == ect_mixed_chan;
			}

			dt_bool is_Double_Chan() const
			{
				return chan_type == eCT_double_chan;
			}

 			dt_bool is_real_space() const
			{
				return space == esp_real;
			}

 			dt_bool is_reciprocal_space() const
			{
				return space == esp_fourier;
			}
		};

		/*
		template <class T, eDev Dev>
		class Energy_Loss
		{
		public:
			using T_r = T;
			using T_c = complex<T>;

			void set_in_data(In_Multem<T_r> *in_multem, Stream<Dev> *stream, FFT<T_r, Dev> *fft2)
			{
				this->in_multem = in_multem;
				this->stream = stream;
				this->fft_2d = fft2;

				if (in_multem->eels_fr.m_sel>2)
				{
					kernel.resize(3);
				}
				else
				{
					kernel.resize(1);
				}

				for(auto ikn = 0; ikn<kernel.size(); ikn++)
				{
					kernel[ikn].resize(in_multem->grid_2d.size());
				}

			}

			void set_atom_type(EELS<T>& eels)
			{
				if (eels.m_sel>2)
				{
					mt::fcn_eels_w_xyz(*stream, in_multem->grid_2d, eels, *fft_2d, kernel[0], kernel[1], kernel[2]);
				}
				else if (eels.m_sel == -2)
				{
					mt::fcn_eels_w_x(*stream, in_multem->grid_2d, eels, *fft_2d, kernel[0]);
				}
				else if (eels.m_sel == -1)
				{
					mt::fcn_eels_w_mn1(*stream, in_multem->grid_2d, eels, *fft_2d, kernel[0]);
				}
				else if (eels.m_sel == 0)
				{
					mt::fcn_eels_w_z(*stream, in_multem->grid_2d, eels, *fft_2d, kernel[0]);
				}
				else if (eels.m_sel == 1)
				{
					mt::fcn_eels_w_mp1(*stream, in_multem->grid_2d, eels, *fft_2d, kernel[0]);
				}
				else if (eels.m_sel == 2)
				{
					mt::fcn_eels_w_y(*stream, in_multem->grid_2d, eels, *fft_2d, kernel[0]);
				}
			}

			Vctr<Vctr<T_c, Dev>, edev_cpu> kernel;
		private:
			In_Multem<T_r> *in_multem;
			Stream<Dev> *stream;
			FFT<T_r, Dev> *fft_2d;
		};
		*/
	}

#endif