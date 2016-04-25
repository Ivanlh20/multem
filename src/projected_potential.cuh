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

#ifndef PROJECTED_POTENTIAL_H
#define PROJECTED_POTENTIAL_H

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "quadrature.hpp"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "specimen.hpp"

namespace mt
{
	template<class T, eDevice dev>
	class Projected_Potential: public Specimen<T>{
		public:
			using value_type_r = T;
			using size_type = std::size_t;

			static const eDevice device = dev;

			Projected_Potential(): stream(nullptr){}

			void set_input_data(Input_Multislice<value_type_r> *input_multislice_i, Stream<dev> *stream_i)
			{	
				Specimen<T>::set_input_data(input_multislice_i);
				stream = stream_i;

				Quadrature quadrature;
				quadrature.get(0, c_nqz, qz); // 0: int_-1^1 y(x) dx - TanhSinh quadrature

				atom_type.resize(c_nAtomsTypes);
				for(auto i = 0; i<atom_type.size(); i++)
				{
					atom_type[i].assign(this->ptr_atom_type()->at(i));
				}

				int nv = max(this->input_multislice->grid.nx_dRx(this->atoms.l_x_Int), this->input_multislice->grid.ny_dRy(this->atoms.l_y_Int));

				stream_data.resize(stream->size());
				for(auto i = 0; i<stream_data.size(); i++)
				{
					stream_data.iv[i].resize(nv);
					stream_data.v[i].resize(nv);
					stream_data.c0[i].resize(c_nR);
					stream_data.c1[i].resize(c_nR);
					stream_data.c2[i].resize(c_nR);
					stream_data.c3[i].resize(c_nR);
				}

				atom_Vp.resize(stream->size());

				V_0.resize(this->input_multislice->grid.nxy());
			}

			void projected_potential(const value_type_r &z_0, const value_type_r &z_e, const int &iatom_0, const int &iatom_e, Vector<value_type_r, dev> &V)
			{
				mt::fill(*stream, V, 0.0);

				int iatom = iatom_0;
				while (iatom <= iatom_e)
				{
					stream->set_n_act_stream(iatom_e-iatom+1);
					set_atom_Vp(z_0, z_e, iatom, *stream, atom_Vp);
					get_cubic_poly_coef_Vz(*stream, atom_Vp);
					mt::eval_cubic_poly<true>(*stream, this->input_multislice->grid, atom_Vp, V);
					iatom += stream->n_act_stream;
				}

				stream->synchronize();
			}

			void projected_potential(const int &islice_0, const int &islice_e, Vector<value_type_r, dev> &V)
			{
				if((islice_0<0) || (islice_e>=this->slice.size()))
				{
					mt::fill(*stream, V, 0.0);
					return;
				}
				projected_potential(this->slice.z_0[islice_0], this->slice.z_e[islice_e], this->slice.iatom_0[islice_0], this->slice.iatom_e[islice_e], V);
			}

			void projected_potential(const int &islice_0, const int &islice_e)
			{
				projected_potential(islice_0, islice_e, V_0);
			}

			void projected_potential(const int &islice, Vector<value_type_r, dev> &V)
			{
				projected_potential(islice, islice, V);
			}

			void projected_potential(const int &islice)
			{
				projected_potential(islice, islice, V_0);
			}

			template<class TOutput_multislice>
			void projected_potential(const int &islice, TOutput_multislice &output_multislice)
			{
				projected_potential(islice, islice, V_0);
				mt::copy_to_host(output_multislice.stream, V_0, output_multislice.V[0]);
				output_multislice.shift();
				output_multislice.clear_temporal_data();
			}

			Vector<value_type_r, dev> V_0;
			Stream<dev> *stream;
		private:
			struct Stream_Data
			{
				using value_type = T;
				using size_type = std::size_t;

				static const eDevice device = dev;

				size_type size() const
				{
					return iv.size();
				}

				void resize(const size_type &new_size)
				{
					iv.resize(new_size);
					v.resize(new_size);
					c0.resize(new_size);
					c1.resize(new_size);
					c2.resize(new_size);
					c3.resize(new_size);
				}

				Vector<Vector<int, dev>, e_host> iv;
				Vector<Vector<T, dev>, e_host> v;
				Vector<Vector<T, dev>, e_host> c0; 		// zero coefficient
				Vector<Vector<T, dev>, e_host> c1; 		// first coefficient
				Vector<Vector<T, dev>, e_host> c2; 		// second coefficient
				Vector<Vector<T, dev>, e_host> c3; 		// third coefficient
			};

			void set_atom_Vp(const value_type_r &z_0, const value_type_r &z_e, int iatom, Stream<dev> &stream, Vector<Atom_Vp<value_type_r, dev>, e_host> &atom_Vp)
			{
				for(auto istream = 0; istream < stream.n_act_stream; istream++)
				{
					auto iZ = this->atoms.Z[iatom]-1;
					auto charge = this->atoms.charge[iatom];
					int icharge = atom_type[iZ].charge_to_idx(charge);
					auto &coef = atom_type[iZ].coef[icharge];

					atom_Vp[istream].charge = charge;
					atom_Vp[istream].x = this->atoms.x[iatom];
					atom_Vp[istream].y = this->atoms.y[iatom];
					atom_Vp[istream].occ = this->atoms.occ[iatom];
					atom_Vp[istream].R2_min = coef.R2_min();
					atom_Vp[istream].R2_max = coef.R2_max();
					atom_Vp[istream].R2 = raw_pointer_cast(coef.R2.data());
					atom_Vp[istream].set_ix0_ixn(this->input_multislice->grid, coef.R_max);
					atom_Vp[istream].set_iy0_iyn(this->input_multislice->grid, coef.R_max);
					atom_Vp[istream].R2_tap = coef.R2_tap();
					atom_Vp[istream].tap_cf = coef.tap_cf;
					atom_Vp[istream].iv = raw_pointer_cast(stream_data.iv[istream].data());
					atom_Vp[istream].v = raw_pointer_cast(stream_data.v[istream].data());

					if(this->input_multislice->is_subslicing())
					{
						atom_Vp[istream].z0h = 0.5*(z_0 - this->atoms.z[iatom]); 
						atom_Vp[istream].zeh = 0.5*(z_e - this->atoms.z[iatom]);
						atom_Vp[istream].split = (atom_Vp[istream].z0h<0)&&(0<atom_Vp[istream].zeh);
						atom_Vp[istream].cl = raw_pointer_cast(coef.Vr.cl.data());
						atom_Vp[istream].cnl = raw_pointer_cast(coef.Vr.cnl.data());
						atom_Vp[istream].c0 = raw_pointer_cast(stream_data.c0[istream].data());
						atom_Vp[istream].c1 = raw_pointer_cast(stream_data.c1[istream].data());
						atom_Vp[istream].c2 = raw_pointer_cast(stream_data.c2[istream].data());
						atom_Vp[istream].c3 = raw_pointer_cast(stream_data.c3[istream].data());
					}
					else
					{
						atom_Vp[istream].c0 = raw_pointer_cast(coef.ciVR.c0.data());
						atom_Vp[istream].c1 = raw_pointer_cast(coef.ciVR.c1.data());
						atom_Vp[istream].c2 = raw_pointer_cast(coef.ciVR.c2.data());
						atom_Vp[istream].c3 = raw_pointer_cast(coef.ciVR.c3.data());
					}
					iatom++;
				}
			}
			
			void get_cubic_poly_coef_Vz(Stream<dev> &stream, Vector<Atom_Vp<value_type_r, dev>, e_host> &atom_Vp)
			{
				if(this->input_multislice->is_subslicing())
				{
					mt::linear_Vz(stream, this->input_multislice->potential_type, qz, atom_Vp);
					mt::cubic_poly_coef(stream, atom_Vp);
				}
			}

			Q1<value_type_r, dev> qz;
			Vector<Atom_Type<value_type_r, dev>, e_host> atom_type; // Atom types

			Stream_Data stream_data;
			Vector<Atom_Vp<value_type_r, dev>, e_host> atom_Vp;
	};

} // namespace mt
#endif