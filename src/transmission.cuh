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
#include "fft2.cuh"
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
			using value_type_r = T;
			using value_type_c = complex<T>;
			using size_type = std::size_t;

			static const eDevice device = dev;

			Transmission():fft2(nullptr){}

			void set_input_data(Input_Multislice<value_type_r, dev> *input_multislice_i, Stream<value_type_r, dev> *stream_i, FFT2<value_type_r, dev> *fft2_i)
			{
				Potential<T, dev>::set_input_data(input_multislice_i, stream_i);
				fft2 = fft2_i;

				trans_0.resize(this->input_multislice->grid.nxy());

				if(!this->input_multislice->slice_storage)
				{
					memory_slice.clear();
					return;
				}

				int n_slice_sig = (this->input_multislice->fp_dim.z)?(int)ceil(3.0*this->atoms.sigma_max/this->input_multislice->grid.dz):0;
				int n_slice_req = this->slice.size() + 2*n_slice_sig;

				memory_slice.set_input_data(n_slice_req, this->input_multislice->grid.nxy());

				if(memory_slice.is_potential())
				{
					memory_slice.resize_vector(this->input_multislice->grid.nxy(), Vp_v);
				}
				else if(memory_slice.is_transmission())
				{
					memory_slice.resize_vector(this->input_multislice->grid.nxy(), trans_v);
				}
			}

			void move_atoms(const int &iconf)
			{
				Potential<T, dev>::move_atoms(iconf);

				value_type_r fPot = this->input_multislice->Vr_factor();

				for(auto i=0; i< memory_slice.n_slice_cur(this->slice.size()); i++)
				{
					if(memory_slice.is_potential())
					{
						this->projected_potential(i, Vp_v[i]);
					}
					else if(memory_slice.is_transmission())
					{
						this->projected_potential(i, this->V0);
						multem::transmission_funtion(this->input_multislice->grid, *fft2, this->input_multislice->interaction_model, fPot, this->V0, trans_v[i]);
					}
				}
			}

			void trans(const int &islice, Vector<value_type_c, dev> &trans_0)
			{
				value_type_r fPot = this->input_multislice->Vr_factor();

				if(islice < memory_slice.n_slice_cur(this->slice.size()))
				{
					if(memory_slice.is_potential())
					{
						multem::transmission_funtion(this->input_multislice->grid, *fft2, this->input_multislice->interaction_model, fPot, Vp_v[islice], trans_0);
					}
					else if(memory_slice.is_transmission())
					{
						multem::assign(trans_v[islice], trans_0);
					}
				}
				else
				{
					this->projected_potential(islice, this->V0);
					multem::transmission_funtion(this->input_multislice->grid, *fft2, this->input_multislice->interaction_model, fPot, this->V0, trans_0);
				}
			}

			void trans(const int &islice_0, const int &islice_e, Vector<value_type_c, dev> &trans_0)
			{
				value_type_r fPot = this->input_multislice->Vr_factor();

				this->projected_potential(islice_0, islice_e, this->V0);
				multem::transmission_funtion(this->input_multislice->grid, *fft2, this->input_multislice->interaction_model, fPot, this->V0, trans_0);
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

			Vector<value_type_c, dev> trans_0;
		private:
			struct Memory_Slice
			{
				public:
					int n_slice_req;
					int n_slice_Allow;
					eSlice_Memory_Type slice_mem_type;

					Memory_Slice():n_slice_req(0), n_slice_Allow(0), slice_mem_type(eSMT_none){}

					void clear()
					{
						n_slice_req = n_slice_Allow = 0;
						slice_mem_type = eSMT_none;
					}

					void set_input_data(const int &nSlice_req_i, const int &nxy_i)
					{
						n_slice_req = nSlice_req_i;
						double free_memory = get_free_memory<dev>() - 10;

						if(number_slices<value_type_c>(free_memory, nxy_i) >= n_slice_req)
						{
							slice_mem_type = eSMT_Transmission;
							n_slice_Allow = number_slices<value_type_c>(free_memory, nxy_i);
						}
						else
						{
							slice_mem_type = eSMT_Potential;
							n_slice_Allow = number_slices<value_type_r>(free_memory, nxy_i);
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
						for(auto i=0; i < n_slice_Allow; i++)
						{
							vector[i].resize(nxy_i);
						}
					}

					int n_slice_cur(const int &n_slice_i)
					{
						return min(n_slice_Allow, n_slice_i);
					}

					bool is_transmission() const
					{
						return slice_mem_type == eSMT_Transmission;
					}

					bool is_potential() const
					{
						return slice_mem_type == eSMT_Potential;
					}

				private:
					template<class U>
					int number_slices(const double &memory, const int &nxy)
					{
						return static_cast<int>(floor(memory/sizeMb<U>(nxy)));
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