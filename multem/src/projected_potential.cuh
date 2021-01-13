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

#ifndef PROJECTED_POTENTIAL_H
#define PROJECTED_POTENTIAL_H

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "quadrature.hpp"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "cpu_fcns.hpp"
#include "gpu_fcns.cuh"
#include "cgpu_fcns.cuh"
#include "spec.hpp"

namespace mt
{
	template <class T, eDevice dev>
	class Projected_Potential: public Spec<T>{
		public:
			using size_type = std::size_t;

			static const eDevice device = dev;

			Projected_Potential(): stream(nullptr), n_atoms_s(512){}

			void set_input_data(Input_Multislice<T> *input_multislice_i, Stream<dev> *stream_i)
			{	
				Spec<T>::set_input_data(input_multislice_i);
				stream = stream_i;

				Quadrature quadrature;
				quadrature(0, c_nqz, qz); // 0: int_-1^1 y(x) dx - TanhSinh quadrature

				atom_type.resize(c_nAtomsTypes);
				for(auto iatom_type = 0; iatom_type<atom_type.size(); iatom_type++)
				{
					atom_type[iatom_type].assign(Spec<T>::atom_type[iatom_type]);
				}

				n_atoms_s = (device==e_host)?(stream->size()):512;

				int nv = max(this->input_multislice->grid_2d.nx_dRx(this->atoms.l_x_int), this->input_multislice->grid_2d.ny_dRy(this->atoms.l_y_int));

				stream_data.resize(n_atoms_s);

				for(auto i = 0; i<stream_data.size(); i++)
				{
					if(device==e_host)
					{
						stream_data.iv[i].resize(nv);
						stream_data.v[i].resize(nv);
					}

					if(this->input_multislice->is_subslicing())
					{
						stream_data.c0[i].resize(c_nR);
						stream_data.c1[i].resize(c_nR);
						stream_data.c2[i].resize(c_nR);
						stream_data.c3[i].resize(c_nR);
					}
				}

				atom_Vp_h.resize(n_atoms_s);
				atom_Vp.resize(n_atoms_s);

				V_0.resize(this->input_multislice->grid_2d.nxy());
			}

			/************************Host************************/
			template<eDevice devn = dev>
			enable_if_dev_host<devn, void>
			operator()(const T &z_0, const T &z_e, const int &iatom_0, const int &iatom_e, Vector<T, dev> &V)
			{
				auto eval_cubic_poly = [](Stream<e_host> &stream, Grid_2d<T> &grid_2d, 
				Vector<Atom_Vp<T>, e_host> &atom, Vector<T, dev> &M_o)
				{
					if(stream.n_act_stream<= 0)
					{
						return;
					}

					for(auto istream = 0; istream < stream.n_act_stream-1; istream++)
					{
						stream[istream] = std::thread(std::bind(host_detail::eval_cubic_poly<T>, std::ref(stream), std::ref(grid_2d), std::ref(atom[istream]), std::ref(M_o)));
					}

					host_detail::eval_cubic_poly<T>(stream, grid_2d, atom[stream.n_act_stream-1], M_o);

					stream.synchronize();
				};

				mt::fill(*stream, V, T(0));

				int iatoms = iatom_0;
				while (iatoms <= iatom_e)
				{
					stream->set_n_act_stream(iatom_e-iatoms+1);
					set_atom_Vp(z_0, z_e, iatoms, stream->n_act_stream, atom_Vp);
					//get_cubic_poly_coef_Vz(*stream, atom_Vp);
					eval_cubic_poly(*stream, this->input_multislice->grid_2d, atom_Vp, V);
					iatoms += stream->n_act_stream;
				}

				stream->synchronize();
			}

			/***********************Device***********************/
		#ifdef __CUDACC__
			template<eDevice devn = dev>
			enable_if_dev_device<devn, void>
			operator()(const T &z_0, const T &z_e, const int &iatom_0, const int &iatom_e, Vector<T, dev> &V)
			{
				auto get_eval_cubic_poly_gridBT = [](int natoms)->Grid_BT
				{
					Grid_BT grid_bt;
					grid_bt.Blk = dim3(natoms, 1, 1);
					grid_bt.Thr = dim3(c_thrnxny, c_thrnxny, 1);

					return grid_bt;
				};

				mt::fill(*stream, V, T(0));

				int iatoms = iatom_0;
				while (iatoms <= iatom_e)
				{
					int n_atoms = min(n_atoms_s, iatom_e-iatoms+1);
					set_atom_Vp(z_0, z_e, iatoms, n_atoms, atom_Vp);
					//get_cubic_poly_coef_Vz(*stream, atom_Vp_h);

					auto grid_bt = get_eval_cubic_poly_gridBT(n_atoms);
					device_detail::eval_cubic_poly<T><<<grid_bt.Blk, grid_bt.Thr>>>(this->input_multislice->grid_2d, atom_Vp, V);

					iatoms += n_atoms;
				}
			}
		#endif

			void operator()(const int &islice_0, const int &islice_e, Vector<T, dev> &V)
			{
				if((islice_0<0) || (islice_e>=this->slicing.slice.size()))
				{
					mt::fill(*stream, V, 0.0);
					return;
				}
				this->operator()(this->slicing.slice[islice_0].z_0, this->slicing.slice[islice_e].z_e, 
				this->slicing.slice[islice_0].iatom_0, this->slicing.slice[islice_e].iatom_e, V);
			}

			void operator()(const int &islice_0, const int &islice_e)
			{
				this->operator()(islice_0, islice_e, V_0);
			}

			void operator()(const int &islice, Vector<T, dev> &V)
			{
				this->operator()(islice, islice, V);
			}

			void operator()(const int &islice)
			{
				this->operator()(islice, islice, V_0);
			}

			template <class TOutput_multislice>
			void operator()(const int &islice, TOutput_multislice &output_multislice)
			{
				this->operator()(islice, islice, V_0);
				mt::copy_to_host(output_multislice.stream, V_0, output_multislice.V[0]);
			}

			Vector<T, dev> V_0;
			Stream<dev> *stream;
		private:
			int n_atoms_s;

			struct Stream_Data
			{
				using value_type = T;
				using size_type = std::size_t;

				static const eDevice device = e_host;

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

			void set_atom_Vp(const T &z_0, const T &z_e, int iatoms, int n_atoms, Vector<Atom_Vp<T>, dev> &atom_Vp)
			{
				for(auto istream = 0; istream < n_atoms; istream++)
				{
					auto iZ = (this->atoms.Z[iatoms] % 1000)-1;
					auto charge = this->atoms.charge[iatoms];
					int icharge = atom_type[iZ].charge_to_idx(charge);
					auto &coef = atom_type[iZ].coef[icharge];

					atom_Vp_h[istream].charge = charge;
					atom_Vp_h[istream].x = this->atoms.x[iatoms];
					atom_Vp_h[istream].y = this->atoms.y[iatoms];
					atom_Vp_h[istream].occ = this->atoms.occ[iatoms];
					atom_Vp_h[istream].R2_min = coef.R2_min();
					atom_Vp_h[istream].R2_max = coef.R2_max();
					atom_Vp_h[istream].R2 = raw_pointer_cast(coef.R2.data());
					atom_Vp_h[istream].set_ix0_ixn(this->input_multislice->grid_2d, coef.R_max);
					atom_Vp_h[istream].set_iy0_iyn(this->input_multislice->grid_2d, coef.R_max);
					atom_Vp_h[istream].R2_tap = coef.R2_tap();
					atom_Vp_h[istream].tap_cf = coef.tap_cf;

					if(device==e_host)
					{
						atom_Vp_h[istream].iv = raw_pointer_cast(stream_data.iv[istream].data());
						atom_Vp_h[istream].v = raw_pointer_cast(stream_data.v[istream].data());
					}

					if(this->input_multislice->is_subslicing())
					{
						atom_Vp_h[istream].z0h = 0.5*(z_0 - this->atoms.z[iatoms]); 
						atom_Vp_h[istream].zeh = 0.5*(z_e - this->atoms.z[iatoms]);
						atom_Vp_h[istream].split = (atom_Vp_h[istream].z0h<0) && (0<atom_Vp_h[istream].zeh);
						atom_Vp_h[istream].cl = raw_pointer_cast(coef.Vr.cl.data());
						atom_Vp_h[istream].cnl = raw_pointer_cast(coef.Vr.cnl.data());
						atom_Vp_h[istream].c0 = raw_pointer_cast(stream_data.c0[istream].data());
						atom_Vp_h[istream].c1 = raw_pointer_cast(stream_data.c1[istream].data());
						atom_Vp_h[istream].c2 = raw_pointer_cast(stream_data.c2[istream].data());
						atom_Vp_h[istream].c3 = raw_pointer_cast(stream_data.c3[istream].data());
					}
					else
					{
						atom_Vp_h[istream].c0 = raw_pointer_cast(coef.ciVR.c0.data());
						atom_Vp_h[istream].c1 = raw_pointer_cast(coef.ciVR.c1.data());
						atom_Vp_h[istream].c2 = raw_pointer_cast(coef.ciVR.c2.data());
						atom_Vp_h[istream].c3 = raw_pointer_cast(coef.ciVR.c3.data());
					}
					iatoms++;
				}
				thrust::copy(atom_Vp_h.begin(), atom_Vp_h.begin()+n_atoms, atom_Vp.begin());
			}
			
			void get_cubic_poly_coef_Vz(Stream<dev> &stream, Vector<Atom_Vp<T>, e_host> &atom_Vp)
			{
				if(this->input_multislice->is_subslicing())
				{
					mt::linear_Vz(stream, this->input_multislice->potential_type, qz, atom_Vp);
					mt::cubic_poly_coef(stream, atom_Vp);
				}
			}

			Q1<T, dev> qz;
			Vector<Atom_Type<T, dev>, e_host> atom_type; // Atom types

			Stream_Data stream_data;
			Vector<Atom_Vp<T>, e_host> atom_Vp_h;
			Vector<Atom_Vp<T>, dev> atom_Vp;
	};

} // namespace mt
#endif
