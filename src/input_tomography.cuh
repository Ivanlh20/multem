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

#ifndef INPUT_TOMOGRAPHY_H
#define INPUT_TOMOGRAPHY_H

#include <algorithm>
#include <vector>

#include "math.cuh"
#include "types.cuh"
#include "memory_info.cuh"
#include "lin_alg_def.cuh"

#include "atomic_data_mt.hpp"
#include "cgpu_fcns.cuh"
#include "cpu_fcns.hpp"
#include "gpu_fcns.cuh"

namespace mt
{
	bool is_gpu_available();

	template <class T>
	class Input_Tomography
	{
		public:
			using value_type = T;

			ePrecision precision;
			eDevice device; 									// eP_float = 1, eP_double = 2
			int cpu_nthread; 									// Number of threads
			int gpu_device; 									// GPU device
			int gpu_nstream; 									// Number of streams

			r3d<T> spec_rot_u0; 										// unitary vector			
			r3d<T> spec_rot_center_p; 										// rotation point

			int nstream;
			Grid_2d<T> grid_2d; 										// grid_bt information

			T r0_min;
			T rTemp;

			eInput_Atoms input_atoms;
			Atom_Data_Sa<T> atoms;

			int Z;
			Vector<T, e_host> r;
			Vector<T, e_host> fr;

			Vector<T, e_host> angle;
			Vector<Vector<T, e_host>, e_host> image;

			Input_Tomography(): precision(eP_double), device(e_host), cpu_nthread(4), 
				gpu_device(0), gpu_nstream(8), spec_rot_u0(1, 0, 0), spec_rot_center_p(1, 0, 0), 
				nstream(0), r0_min(0.75), rTemp(0.8), input_atoms(eIA_no), Z(0){};

			void validate_parameters()
			{
				cpu_nthread = max(1, cpu_nthread);
				gpu_nstream = max(1, gpu_nstream);
				nstream = (is_host())?cpu_nthread:gpu_nstream;

				if(!is_float() && !is_double())
					precision = eP_float;

				if(!is_host() && !is_device())
					device = e_host;

				set_device();

				spec_rot_u0.normalized();

				if(!is_input_atoms())
				{
					T lxyz = ::fmax(grid_2d.lx, grid_2d.ly);
					r3d<T> r_min(0, 0, 0);
					r3d<T> r_max(lxyz, lxyz, lxyz);
					atoms.set_range(Z, r_min, r_max);
				}

				if(r0_min<= 0)
				{
					r0_min = 0.75;
				}

				if((rTemp<= 0)||(rTemp>= 1))
				{
					rTemp = 0.8;
				}

			}

			bool is_input_atoms() const
			{
				return input_atoms == mt::eIA_yes;
			}

			bool is_host() const
			{
				return device == mt::e_host;
			}

			bool is_device() const
			{
				return device == mt::e_device;
			}

			void set_device()
			{
				if(is_device())
				{
					if(!is_gpu_available())
					{
						device = mt::e_host;
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
					device = mt::e_host;
				}
			}

			bool is_float() const
			{
				return precision == mt::eP_float;
			}

			bool is_double() const
			{
				return precision == mt::eP_double;
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

} // namespace mt

#endif