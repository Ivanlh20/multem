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
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef SUPERPOSITION_H
#define SUPERPOSITION_H

#include <random>
#include <algorithm>

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "r3d.cuh"
#include "stream.cuh"
#include "atomic_data.hpp"
#include "input_output_superposition.cuh"
#include "box_occ.hpp"
#include "cubic_spline.hpp"

namespace multem
{
	template<class T, eDevice dev>
	class Superposition
	{
		public:
			using value_type_r = T;
			using size_type = std::size_t;

			static const eDevice device = dev;

			Superposition(): input_superposition(nullptr){}

			void set_input_data(Input_Superposition<T> *input_superposition_i)
			{	
				input_superposition = input_superposition_i;
				stream.resize(input_superposition->nstream);

				int nv = input_superposition->get_nv();

				stream_data.resize(stream.size());
				for(auto i = 0; i<stream_data.size(); i++)
				{
					stream_data.iv[i].resize(nv);
					stream_data.v[i].resize(nv);
					stream_data.R2[i].resize(c_nR);
					stream_data.c0[i].resize(c_nR);
					stream_data.c1[i].resize(c_nR);
					stream_data.c2[i].resize(c_nR);
					stream_data.c3[i].resize(c_nR);
				}
				atom_sp.resize(stream.size());
				Im.resize(input_superposition->grid.nxy());
			}

			void exec(Vector<T, dev> &Im)
			{
				multem::fill(stream, Im, 0.0);

				int iatom_0 = 0;
				int iatom_e = input_superposition->atoms.size()-1;

				int iatom = iatom_0;
				while (iatom <= iatom_e)
				{
					stream.set_n_act_stream(iatom_e-iatom+1);
					set_atom_sp(iatom, stream, atom_sp);

					multem::linear_Gaussian(stream, atom_sp);
					multem::cubic_poly_coef(stream, atom_sp);
					multem::eval_cubic_poly<false>(stream, input_superposition->grid, atom_sp, Im);
					iatom += stream.n_act_stream;
				}

				stream.synchronize();
			}

			template<class Output_Superposition>
			void run(Output_Superposition &output_superposition)
			{
				exec(Im);
				multem::copy_to_host(output_superposition.stream, Im, output_superposition.Im);
			}

			Input_Superposition<T> *input_superposition;
			Stream<dev> stream;
			Vector<value_type_r, dev> Im;
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
					R2.resize(new_size);
					c0.resize(new_size);
					c1.resize(new_size);
					c2.resize(new_size);
					c3.resize(new_size);
				}

				Vector<Vector<int, dev>, e_host> iv;
				Vector<Vector<T, dev>, e_host> v;
				Vector<Vector<T, dev>, e_host> R2; 		
				Vector<Vector<T, dev>, e_host> c0; 		// zero coefficient
				Vector<Vector<T, dev>, e_host> c1; 		// first coefficient
				Vector<Vector<T, dev>, e_host> c2; 		// second coefficient
				Vector<Vector<T, dev>, e_host> c3; 		// third coefficient
			};

			void set_atom_sp(int iatom, Stream<dev> &stream, Vector<Atom_Sp<T, dev>, e_host> &atom_sp)
			{
				for(auto istream = 0; istream < stream.n_act_stream; istream++)
				{
					atom_sp[istream].x = input_superposition->atoms.x[iatom];
					atom_sp[istream].y = input_superposition->atoms.y[iatom];
					atom_sp[istream].occ = 1;
					atom_sp[istream].a = input_superposition->atoms.a[iatom];
					atom_sp[istream].alpha = input_superposition->alpha(iatom);
					auto R_max = input_superposition->R_max(iatom);
					auto R2_max = R_max*R_max;
					atom_sp[istream].dtR = (exp(-atom_sp[istream].alpha*R2_max)-1)/T(c_nR-1);
					atom_sp[istream].R2_tap = pow(0.85*R_max, 2);
					atom_sp[istream].tap_cf = c_i2Pi/(R2_max-atom_sp[istream].R2_tap);
					atom_sp[istream].R2_max = R2_max;
					atom_sp[istream].R2 = raw_pointer_cast(stream_data.R2[istream].data());
					atom_sp[istream].set_ix0_ixn(input_superposition->grid, R_max);
					atom_sp[istream].set_iy0_iyn(input_superposition->grid, R_max);
					atom_sp[istream].c0 = raw_pointer_cast(stream_data.c0[istream].data());
					atom_sp[istream].c1 = raw_pointer_cast(stream_data.c1[istream].data());
					atom_sp[istream].c2 = raw_pointer_cast(stream_data.c2[istream].data());
					atom_sp[istream].c3 = raw_pointer_cast(stream_data.c3[istream].data());
					atom_sp[istream].iv = raw_pointer_cast(stream_data.iv[istream].data());
					atom_sp[istream].v = raw_pointer_cast(stream_data.v[istream].data());
					iatom++;
				}
			}

			Stream_Data stream_data;
			Vector<Atom_Sp<T, dev>, e_host> atom_sp;
	};

} // namespace multem
#endif