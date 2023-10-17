/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>transmission_fcn
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

#ifndef TRANSMISSION_FCNS_H
#define TRANSMISSION_FCNS_H

#include "math.cuh"
#include "types.cuh"
#include "type_traits_gen.cuh"
#include "cgpu_stream.cuh"
#include "quad_data.cuh"
#include "cgpu_info.cuh"
#include "in_classes.cuh"
#include "output_multem.hpp"
#include "projected_potential.cuh"

namespace mt
{
	template <class T, eDev Dev>
	class Transmission_Fcn: public Projected_Potential<T, Dev>
	{
		public:
			using T_r = T;
			using T_c = complex<T>;
			using size_type = dt_uint64;

			Transmission_Fcn(): Projected_Potential<T, Dev>(), fft_2d(nullptr) {}

			void set_in_data(Multem_In_Parm<T_r> *multem_in_parm_i, Stream<Dev> *stream_i, FFT<T_r, Dev> *fft2_i)
			{
				Projected_Potential<T, Dev>::set_in_data(multem_in_parm_i, stream_i);
				fft_2d = fft2_i;

				trans_0.resize(this->multem_in_parm->grid_2d.size());

				if (!this->multem_in_parm->slice_storage)
				{
					memory_slice.clear();
					return;
				}

				dt_int32 n_slice_sig = (this->multem_in_parm->atomic_vib.dim_z)?(dt_int32)ceil(3.0*this->atoms.sigma_max/this->multem_in_parm->grid_2d.sli_thick):0;
				dt_int32 n_slice_req = this->slicing.slice.size() + 2*n_slice_sig;

				memory_slice.set_in_data(n_slice_req, this->multem_in_parm->grid_2d.size());

				if (memory_slice.is_potential())
				{
					memory_slice.resize_vector(this->multem_in_parm->grid_2d.size(), Vp_v);
				}
				else if (memory_slice.is_transmission())
				{
					memory_slice.resize_vector(this->multem_in_parm->grid_2d.size(), trans_v);
				}
			}

			void trans(T_r w, Vctr<T_r, Dev>& V0_i, Vctr<T_c, Dev>& Trans_o)
			{	
				mt::transmission_fcn(*(this->stream), this->multem_in_parm->grid_2d, this->multem_in_parm->elec_spec_interact_mod, w, V0_i, Trans_o);

				if (this->multem_in_parm->grid_2d.bwl)
				{
					fft_2d->forward(Trans_o);
					mt::fcn_fermi_aperture(*(this->stream), this->multem_in_parm->grid_2d, Trans_o);
					fft_2d->inverse(Trans_o);
				}
			}

			void trans(const dt_int32& islice, Vctr<T_c, Dev>& trans_0)
			{
				if (islice < memory_slice.n_slice_cur(this->slicing.slice.size()))
				{
					if (memory_slice.is_potential())
					{
						trans(this->multem_in_parm->Vr_factor(), Vp_v[islice], trans_0);
					}
					else if (memory_slice.is_transmission())
					{
						mt::assign(trans_v[islice], trans_0);
					}
				}
				else
				{
					Projected_Potential<T, Dev>::operator()(islice, this->V_0);
					// this->operator()(islice, this->V_0);
					trans(this->multem_in_parm->Vr_factor(), this->V_0, trans_0);
				}
			}

			void trans(const dt_int32& islice_0, const dt_int32& islice_e, Vctr<T_c, Dev>& trans_0)
			{
				Projected_Potential<T, Dev>::operator()(islice_0, islice_e, this->V_0);
				// this->operator()(islice_0, islice_e, this->V_0);
				trans(this->multem_in_parm->Vr_factor(), this->V_0, trans_0);
			}

			template <class TOutput_multislice>
			void trans(const dt_int32& islice, TOutput_multislice &output_multem)
			{
				trans(islice, trans_0);
				mt::cpy_to_host(output_multem.stream, trans_0, output_multem.trans[0]);
			}

			void move_atoms(const dt_int32& fp_iconf)
			{
				Projected_Potential<T, Dev>::move_atoms(fp_iconf);

				// Calculate transmission functions
				for(auto islice = 0; islice< memory_slice.n_slice_cur(this->slicing.slice.size()); islice++)
				{
					if (memory_slice.is_potential())
					{
						Projected_Potential<T, Dev>::operator()(islice, Vp_v[islice]);
					}
					else if (memory_slice.is_transmission())
					{
						Projected_Potential<T, Dev>::operator()(islice, this->V_0);
						trans(this->multem_in_parm->Vr_factor(), this->V_0, trans_v[islice]);
					}
				}
			}

			void transmit(const dt_int32& islice, Vctr<T_c, Dev>& psi_io)
			{
				trans(islice, trans_0);
				mt::ew_mult(*(this->stream), trans_0, psi_io);
			}

			Vctr<T_c, Dev> trans_0;
		private:
			struct Memory_Slice
			{
				public:
					dt_int32 n_slice_req;
					dt_int32 n_slice_Allow;
					eSlice_Memory_Typ slice_mem_type;

					Memory_Slice(): n_slice_req(0), n_slice_Allow(0), slice_mem_type(esmt_none) {}

					void clear()
					{
						n_slice_req = n_slice_Allow = 0;
						slice_mem_type = esmt_none;
					}

					void set_in_data(const dt_int32& nSlice_req_i, const dt_int32& nxy_i)
					{
						n_slice_req = nSlice_req_i;
						dt_float64 free_memory = fcn_free_mem<Dev>() - 10;

						if (number_slices<T_c>(free_memory, nxy_i) >= n_slice_req)
						{
							slice_mem_type = esmt_transmission;
							n_slice_Allow = number_slices<T_c>(free_memory, nxy_i);
						}
						else
						{
							slice_mem_type = esmt_potential;
							n_slice_Allow = number_slices<T_r>(free_memory, nxy_i);
						}
						n_slice_Allow = min(n_slice_Allow, n_slice_req);

						if (n_slice_Allow == 0 )
						{
							slice_mem_type = esmt_none;
						}
					}

					template <class U>
					void resize_vector(const dt_int32& nxy_i, U &vector)
					{
						vector.resize(n_slice_Allow);
						for(auto i = 0; i < n_slice_Allow; i++)
						{
							vector[i].resize(nxy_i);
						}
					}

					dt_int32 n_slice_cur(const dt_int32& n_slice_i)
					{
						return min(n_slice_Allow, n_slice_i);
					}

					dt_bool is_transmission() const
					{
						return slice_mem_type == esmt_transmission;
					}

					dt_bool is_potential() const
					{
						return slice_mem_type == esmt_potential;
					}

				private:
					template <class U>
					dt_int32 number_slices(const dt_float64 &memory, const dt_int32& nxy)
					{
						return static_cast<dt_int32>(::floor(memory/mt::fcn_size_mb<U>(nxy)));
					}
			};

			Memory_Slice memory_slice;

		protected:
			Vctr<Vctr<T_c, Dev>, edev_cpu> trans_v;
			Vctr<Vctr<T_r, Dev>, edev_cpu> Vp_v;

			FFT<T_r, Dev> *fft_2d;
	};

}

#endif