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

#ifndef OUTPUT_TOMOGRAPHY_H
#define OUTPUT_TOMOGRAPHY_H

#include "types.cuh"
#include "traits.cuh"
#include "atomic_data_mt.hpp"
#include "input_tomography.cuh"

namespace mt
{
	template <class TVector_r>
	class Output_Tomography: public Input_Tomography<Value_type<TVector_r>>
	{
		public:
			using T_r = Value_type<TVector_r>;
			static const bool is_vector = !is_rmatrix_r<TVector_r>::value;

			template <class TInput_Tomography>
			void set_input_data(TInput_Tomography *input_tomography_i)
			{ 
				set_input_tomography(input_tomography_i);
				if(is_vector)
				{
					// temp.resize(input_tomography_i->atoms.size());
					// chi2.resize(input_tomography_i->atoms.size());

					Z.resize(input_tomography_i->atoms.size());
					x.resize(input_tomography_i->atoms.size());
					y.resize(input_tomography_i->atoms.size());
					z.resize(input_tomography_i->atoms.size());
				}
			}

			TVector_r temp;
			TVector_r chi2;
			TVector_r Z;
			TVector_r x;
			TVector_r y;
			TVector_r z;
		private:
			template <class TInput_Tomography>
			void set_input_tomography(TInput_Tomography *input_tomography)
			{ 
				this->precision = input_tomography->precision;
				this->device = input_tomography->device;
				this->cpu_nthread = input_tomography->cpu_nthread;
				this->gpu_device = input_tomography->gpu_device;
				this->gpu_nstream = input_tomography->gpu_nstream;

	
				this->spec_rot_u0 = input_tomography->spec_rot_u0;
				this->spec_rot_center_p = input_tomography->spec_rot_center_p;

				this->grid_2d = input_tomography->grid_2d;

				temp.clear();
				chi2.clear();
			}
	};

	using Output_Tomography_Matlab = Output_Tomography<rmatrix_r>;

	template <class T>
	using Output_Tomography_Vector = Output_Tomography<Vector<T, e_host>>;

} // namespace mt

#endif