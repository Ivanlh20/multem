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

#ifndef INPUT_SUPERPOSITION_H
#define INPUT_SUPERPOSITION_H

#include <algorithm>
#include <vector>

#include "math.cuh"
#include "types.cuh"
#include "matlab_types.cuh"
#include "memory_info.cuh"

#include "host_device_functions.cuh"
#include "host_functions.hpp"
#include "device_functions.cuh"

namespace multem
{
	bool is_gpu_available();

	template<class T>
	class Input_Superposition
	{
		public:
			using value_type = T;

			ePrecision precision;
			eDevice device; 									// eP_float = 1, eP_double = 2
			int cpu_nthread; 									// Number of threads
			int gpu_device; 									// GPU device
			int gpu_nstream; 									// Number of streams

			int nstream;

			T ff_sigma;
			Grid<T> grid; 										// gridBT information
			Atom_Data_Sp<T> atoms;

			Input_Superposition(): precision(eP_double), device(e_host), cpu_nthread(4), 
				gpu_device(0), gpu_nstream(8), ff_sigma(3), nstream(1){};

			void validate_parameters()
			{
				cpu_nthread = max(1, cpu_nthread);
				gpu_nstream = max(1, gpu_nstream);
				nstream = (is_host())?cpu_nthread:gpu_nstream;

				if(!is_float() && !is_double())
					precision = eP_float;

				if(!is_host() && !is_device())
					device = e_host;

				ff_sigma = max(ff_sigma, T(0));

				set_device();
			}

			T alpha(const int &iatom)
			{
				return 0.5/pow(atoms.sigma[iatom], 2);
			}

			T R_max(const int &iatom)
			{
				return ff_sigma*atoms.sigma[iatom];
			}

			// maximum number of pixels
			int get_nv()
			{
				auto l_x = atoms.l_x + 2*atoms.sigma_max;
				auto l_y = atoms.l_y + 2*atoms.sigma_max;
				return max(grid.nx_dRx(l_x), grid.ny_dRy(l_y));
			}

			bool is_host() const
			{
				return device == multem::e_host;
			}

			bool is_device() const
			{
				return device == multem::e_device;
			}

			void set_device()
			{
				if(is_device())
				{
					if(!is_gpu_available())
					{
						device = multem::e_host;
					} 
					else
					{
						auto ngpu = number_of_gpu_available();
						gpu_device = min(max(0, gpu_device), ngpu-1);
						cudaSetDevice(gpu_device);
					}
				}
				else
				{
					device = multem::e_host;
				}
			}

			bool is_float() const
			{
				return precision == multem::eP_float;
			}

			bool is_double() const
			{
				return precision == multem::eP_double;
			}

			bool is_float_host() const
			{
				return is_float() && is_host();
			}

			bool is_double_host() const
			{
				return is_double() && is_host();
			}

			bool is_float_device() const
			{
				return is_float() && is_device();
			}

			bool is_double_device() const
			{
				return is_double() && is_device();
			}
	};

	template<class TVector_r>
	class Output_Superposition: public Input_Superposition<Value_type<TVector_r>>
	{
		public:
			using value_type_r = Value_type<TVector_r>;
			static const bool is_vector = !is_rmatrix_r<TVector_r>::value;

			template<class TInput_Superposition>
			void set_input_data(TInput_Superposition *input_superposition_i)
			{ 
				set_input_superposition(input_superposition_i);
				stream.resize(this->cpu_nthread);
			}

			Stream<e_host> stream;
			TVector_r Im;
		private:
			template<class TInput_Superposition>
			void set_input_superposition(TInput_Superposition *input_superposition)
			{ 
				this->precision = input_superposition->precision;
				this->device = input_superposition->device;
				this->cpu_nthread = input_superposition->cpu_nthread;
				this->gpu_device = input_superposition->gpu_device;
				this->gpu_nstream = input_superposition->gpu_nstream;

				this->grid = input_superposition->grid;
			}
	};

	using Output_Superposition_Matlab = Output_Superposition<rmatrix_r>;

	template<class T>
	using Output_Superposition_Vector = Output_Superposition<Vector<T, e_host>>;
} // namespace multem

#endif