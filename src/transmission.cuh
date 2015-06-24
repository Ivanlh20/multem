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

#ifndef TRANSMISSION_H
#define TRANSMISSION_H

#include "math.cuh"
#include "types.hpp"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "quadrature.hpp"
#include "input_multislice.hpp"
#include "potential.cuh"
#include "memory_info.cuh"

namespace multem
{
	template<class T, eDevice dev>
	class Transmission: public Potential<T, dev>
	{
		public:
			using value_type_c = typename complex<T>;

			Vector<value_type_c, dev> trans_0;

			Transmission():fft2(nullptr) {}

			void set_input_data(Input_Multislice<value_type_r, dev> *input_multislice_io, Stream<value_type_r, dev> *stream_i, FFT2<value_type_r, dev> *fft2_i)
			{
				Potential::set_input_data(input_multislice_io, stream_i);

				fft2 = fft2_i;

				trans_0.resize(input_multislice->grid.nxy());

				if((!input_multislice->fast_cal)||(input_multislice->interaction_model!=eESIM_Multislice))
				{
					memory_slice.clear();
					return;
				}

				int n_slice_sig = (input_multislice->fp_dim.z)?(int)ceil(3.0*atoms.sigma_max/input_multislice->grid.dz):0;
				int n_slice_req = slice.size() + 2*n_slice_sig;

				memory_slice.set_input_data(n_slice_req, input_multislice->grid.nxy());

				if(memory_slice.slice_mem_type==eSMT_Transmission)
				{
					memory_slice.resize_vector(input_multislice->grid.nxy(), trans_v);
				}
				else
				{
					memory_slice.resize_vector(input_multislice->grid.nxy(), Vp_v);
				}
			}

			void move_atoms(const int &iconf)
			{
				Potential::move_atoms(iconf);

				value_type_r fPot = input_multislice->Vr_factor();

				for(auto i=0; i< memory_slice.n_slice_cur(slice.size()); i++)
				{
					if(memory_slice.slice_mem_type==eSMT_Transmission)
					{
						projected_potential(i);
						transmission_fun(input_multislice->grid, input_multislice->interaction_model, fPot, V0, trans_v[i]);
						bandwidth_limit(input_multislice->grid, *fft2, trans_v[i]);
					}
					else
					{
						projected_potential(i, Vp_v[i]);
					}
				}
			}

			void trans(const int &islice, Vector<value_type_c, dev> &trans_0)
			{
				value_type_r fPot = input_multislice->Vr_factor();

				if(islice < memory_slice.n_slice_cur(slice.size()))
				{
					if(memory_slice.slice_mem_type==eSMT_Transmission)
					{
						trans_0.assign(trans_v[islice].begin(), trans_v[islice].end());
					}
					else
					{
						transmission_fun(input_multislice->grid, input_multislice->interaction_model, fPot, Vp_v[islice], trans_0);
						bandwidth_limit(input_multislice->grid, *fft2, trans_0);
					}
				}
				else
				{
					projected_potential(islice);
					transmission_fun(input_multislice->grid, input_multislice->interaction_model, fPot, V0, trans_0);
					bandwidth_limit(input_multislice->grid, *fft2, trans_0); 
				}
			}

			void trans(const int &islice)
			{
				trans(islice, trans_0);
			}

			void transmit(const int &islice, Vector<value_type_c, dev> &psi_io)
			{
				trans(islice);
				multem::multiply(trans_0, psi_io);
			}

		private:
			struct Memory_Slice
			{
				double free_memory;
				double total_memory;
				int n_slice_req;
				int n_slice_Allow;
				eSlice_Memory_Type slice_mem_type;

				Memory_Slice():free_memory(0), total_memory(0), n_slice_req(0), 
					n_slice_Allow(0), slice_mem_type(eSMT_none) {}

				void clear()
				{
					free_memory = total_memory = 0;
					n_slice_req = n_slice_Allow = 0;
					slice_mem_type = eSMT_none;
				}

				void set_input_data(const int &nSlice_req_i, const int &nxy_i)
				{
					n_slice_req = nSlice_req_i;
					memory_info<dev>(total_memory, free_memory);
					free_memory = free_memory - 10.0;

					if(free_memory/sizeMb<value_type_c>(nxy_i) >= n_slice_req)
					{
						slice_mem_type = eSMT_Transmission;
						n_slice_Allow = static_cast<int>(floor(free_memory/sizeMb<value_type_c>(nxy_i)));
					}
					else
					{
						slice_mem_type = eSMT_Potential;
						n_slice_Allow = static_cast<int>(floor(free_memory/sizeMb<value_type_r>(nxy_i)));
					}
					n_slice_Allow = min(n_slice_Allow, n_slice_req);

					if(n_slice_Allow == 0 )
					{
						slice_mem_type = eSMT_none;
					}
				}

				template<class U>
				void resize_vector(const int &nxy_i, U &vector)
				{
					vector.resize(n_slice_Allow);
					//vector.shrink_to_fit(); //this line produce a error --> thrust library
					for(auto i=0; i<n_slice_Allow; i++)
					{
						vector[i].resize(nxy_i);
						vector[i].shrink_to_fit();
					}
				}

				int n_slice_cur(const int &n_slice_i)
				{
					return min(n_slice_Allow, n_slice_i);
				}

			};

			Memory_Slice memory_slice;

		protected:
			Vector<Vector<value_type_c, dev>, e_Host> trans_v;
			Vector<Vector<value_type_r, dev>, e_Host> Vp_v;

			FFT2<value_type_r, dev> *fft2;
	};

} // namespace multem

#endif