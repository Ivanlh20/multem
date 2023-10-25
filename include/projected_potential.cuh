/*
 * This file is part of Multem.
 * Copyright 2022 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef PROJECTED_POTENTIAL_H
#define PROJECTED_POTENTIAL_H

#include "math_mt.h"
#include "types.cuh"
#include "type_traits_gen.h"
#include "cgpu_stream.cuh"
#include "quad_data.cuh"
#include "in_classes.cuh"
#include "output_multem.hpp"
#include "fcns_cpu.h"
#include "fcns_gpu.h"
#include "fcns_gpu.h"
#include "spec.hpp"

namespace mt
 
	template <class T, eDev Dev>
	class Projected_Potential: public Spec<T>{
		public:
			using size_type = dt_uint64;

			static const eDev device = Dev;

			Projected_Potential(): stream(nullptr), n_atoms_s(512) {}

			void set_in_data(Multem_In_Parm<T> *multem_in_parm_i, Stream<Dev> *stream_i)
			{	
				Spec<T>::set_in_data(multem_in_parm_i);
				stream = stream_i;

				Quad_Data quad_data;
				quad_data(0, c_nqz, qz);		// 0: int_-1^1 y(x) dx - tanh_sinh quad_data

				atom_type.resize(c_n_atom_typ);
				for(auto iatom_type = 0; iatom_type<atom_type.size(); iatom_type++)
				{
					atom_type[iatom_type].assign(Spec<T>::atom_type[iatom_type]);
				}

				n_atoms_s = (device==edev_cpu)?(stream->size()):256;

				dt_int32 nv = max(this->multem_in_parm->grid_2d.rx_2_irx_cd(this->atoms.bs_x_int), this->multem_in_parm->grid_2d.ry_2_iry_cd(this->atoms.bs_y_int));

				stream_data.resize(n_atoms_s);

				for(auto i = 0; i<stream_data.size(); i++)
				{
					if (device==edev_cpu)
					{
						stream_data.iv[i].resize(nv);
						stream_data.v[i].resize(nv);
					}

					if (this->multem_in_parm->is_spec_slic_by_dz_sub())
					{
						stream_data.c0[i].resize(c_nR);
						stream_data.c1[i].resize(c_nR);
						stream_data.c2[i].resize(c_nR);
						stream_data.c3[i].resize(c_nR);
					}
				}

				atom_Vp_h.resize(n_atoms_s);
				atom_Vp.resize(n_atoms_s);

				V_0.resize(this->multem_in_parm->grid_2d.size());
			}

			/***************************************** cpu *****************************************/
			template <eDev devn = Dev>
			enable_if_edev_cpu<devn, void>
			operator()(const T& z_0, const T& z_e, const dt_int32& iatom_0, const dt_int32& iatom_e, Vctr<T, Dev>& V)
			{
				auto fcn_eval_poly3 = [](Stream<edev_cpu>& stream, Grid_2d<T>& grid_2d, 
				Vctr<Ptc_pVp<T>, edev_cpu>& atom, Vctr<T, Dev>& mx_o)
				{
					if (stream.n_stream_act<= 0)
					{
						return;
					}

					for(auto istm = 0; istm < stream.n_stream_act-1; istm++)
					{
						stream[istm] = std::thread(std::bind(cpu_detail::fcn_eval_poly3<T>, std::ref(stream), std::ref(grid_2d), std::ref(atom[istm]), std::ref(mx_o)));
					}

					cpu_detail::fcn_eval_poly3<T>(stream, grid_2d, atom[stream.n_stream_act-1], mx_o);

					stream.synchronize();
				};

				mt::fill(*stream, V, T(0));

				dt_int32 iatoms = iatom_0;
				while (iatoms <= iatom_e)
				{
					stream->set_n_stream_act(iatom_e-iatoms+1);
					set_atom_Vp(z_0, z_e, iatoms, stream->n_stream_act, atom_Vp);
					// get_cubic_poly_coef_Vz(*stream, atom_Vp);
					fcn_eval_poly3(*stream, this->multem_in_parm->grid_2d, atom_Vp, V);
					iatoms += stream->n_stream_act;
				}

				stream->synchronize();
			}

			/*************************************** device ****************************************/
		#ifdef __CUDACC__
			template <eDev devn = Dev>
			enable_if_edev_gpu<devn, void>
			operator()(const T& z_0, const T& z_e, const dt_int32& iatom_0, const dt_int32& iatom_e, Vctr<T, Dev>& V)
			{
				auto get_eval_cubic_poly_gridBT = [](dt_int32 natoms)->D_Grid_Blk
				{
					D_Grid_Blk d_grid_blk;
					d_grid_blk.grid = dim3(natoms, 1, 1);
					d_grid_blk.blk = dim3(c_thr_2d_x, c_thr_2d_y, 1);

					return d_grid_blk;
				};

				mt::fill(*stream, V, T(0));

				dt_int32 iatoms = iatom_0;
				while (iatoms <= iatom_e)
				{
					dt_int32 n_atoms = min(n_atoms_s, iatom_e-iatoms+1);
					set_atom_Vp(z_0, z_e, iatoms, n_atoms, atom_Vp);
					// get_cubic_poly_coef_Vz(*stream, atom_Vp_h);

					auto d_grid_blk = get_eval_cubic_poly_gridBT(n_atoms);
					gpu_detail::fcn_eval_poly3<T><<<d_grid_blk.grid, d_grid_blk.blk>>>(this->multem_in_parm->grid_2d, atom_Vp, V);

					iatoms += n_atoms;
				}
			}
		#endif

			void operator()(const dt_int32& islice_0, const dt_int32& islice_e, Vctr<T, Dev>& V)
			{
				if ((islice_0<0) || (islice_e>=this->slicing.slice.size()))
				{
					mt::fill(*stream, V, 0.0);
					return;
				}
				this->operator()(this->slicing.slice[islice_0].z_0, this->slicing.slice[islice_e].z_e, 
				this->slicing.slice[islice_0].iatom_0, this->slicing.slice[islice_e].iatom_e, V);
			}

			void operator()(const dt_int32& islice_0, const dt_int32& islice_e)
			{
				this->operator()(islice_0, islice_e, V_0);
			}

			void operator()(const dt_int32& islice, Vctr<T, Dev>& V)
			{
				this->operator()(islice, islice, V);
			}

			void operator()(const dt_int32& islice)
			{
				this->operator()(islice, islice, V_0);
			}

			template <class TOutput_multislice>
			void operator()(const dt_int32& islice, TOutput_multislice &output_multem)
			{
				this->operator()(islice, islice, V_0);
				mt::cpy_to_host(output_multem.stream, V_0, output_multem.V[0]);
			}

			Vctr<T, Dev> V_0;
			Stream<Dev> *stream;
		private:
			dt_int32 n_atoms_s;

			struct Stream_Data
			{
				using value_type = T;
				using size_type = dt_uint64;

				static const eDev device = edev_cpu;

				size_type size() const
				{
					return iv.size();
				}

				void resize(const size_type& new_size)
				{
					iv.resize(new_size);
					v.resize(new_size);
					c0.resize(new_size);
					c1.resize(new_size);
					c2.resize(new_size);
					c3.resize(new_size);
				}

				Vctr<Vctr<dt_int32, Dev>, edev_cpu> iv;
				Vctr<Vctr<T, Dev>, edev_cpu> v;
				Vctr<Vctr<T, Dev>, edev_cpu> c0;		// zero coefficient
				Vctr<Vctr<T, Dev>, edev_cpu> c1;		// first coefficient
				Vctr<Vctr<T, Dev>, edev_cpu> c2;		// second coefficient
				Vctr<Vctr<T, Dev>, edev_cpu> c3;		// third coefficient
			};

			void set_atom_Vp(const T& z_0, const T& z_e, dt_int32 iatoms, dt_int32 n_atoms, Vctr<Ptc_pVp<T>, Dev>& atom_Vp)
			{
				for(auto istm = 0; istm < n_atoms; istm++)
				{
					auto iZ = this->atoms.Z[iatoms]-1;
					auto charge = this->atoms.charge[iatoms];
					dt_int32 icharge = atom_type[iZ].charge_to_ind(charge);
					auto &coef = atom_type[iZ].coef[icharge];

					atom_Vp_h[istm].charge = charge;
					atom_Vp_h[istm].x = this->atoms.x[iatoms];
					atom_Vp_h[istm].y = this->atoms.y[iatoms];
					atom_Vp_h[istm].occ = this->atoms.occ[iatoms];
					atom_Vp_h[istm].R2_min = coef.R2_min();
					atom_Vp_h[istm].R2_max = coef.R2_max();
					atom_Vp_h[istm].R2 = raw_pointer_cast(coef.R2.data());
					atom_Vp_h[istm].set_ix0_ixn(this->multem_in_parm->grid_2d, coef.R_max);
					atom_Vp_h[istm].set_iy0_iyn(this->multem_in_parm->grid_2d, coef.R_max);
					atom_Vp_h[istm].R2_tap = coef.R2_tap();
					atom_Vp_h[istm].coef_tap = coef.coef_tap;

					if (device==edev_cpu)
					{
						atom_Vp_h[istm].iv = raw_pointer_cast(stream_data.iv[istm].data());
						atom_Vp_h[istm].v = raw_pointer_cast(stream_data.v[istm].data());
					}

					if (this->multem_in_parm->is_spec_slic_by_dz_sub())
					{
						atom_Vp_h[istm].z0h = 0.5*(z_0 - this->atoms.z[iatoms]);
						atom_Vp_h[istm].zeh = 0.5*(z_e - this->atoms.z[iatoms]);
						atom_Vp_h[istm].split = (atom_Vp_h[istm].z0h<0) && (0<atom_Vp_h[istm].zeh);
						atom_Vp_h[istm].cl = raw_pointer_cast(coef.Vr.cl.data());
						atom_Vp_h[istm].cnl = raw_pointer_cast(coef.Vr.cnl.data());
						atom_Vp_h[istm].c0 = raw_pointer_cast(stream_data.c0[istm].data());
						atom_Vp_h[istm].c1 = raw_pointer_cast(stream_data.c1[istm].data());
						atom_Vp_h[istm].c2 = raw_pointer_cast(stream_data.c2[istm].data());
						atom_Vp_h[istm].c3 = raw_pointer_cast(stream_data.c3[istm].data());
					}
					else
					{
						atom_Vp_h[istm].c0 = raw_pointer_cast(coef.ciVR.c0.data());
						atom_Vp_h[istm].c1 = raw_pointer_cast(coef.ciVR.c1.data());
						atom_Vp_h[istm].c2 = raw_pointer_cast(coef.ciVR.c2.data());
						atom_Vp_h[istm].c3 = raw_pointer_cast(coef.ciVR.c3.data());
					}
					iatoms++;
				}
				thrust::copy(atom_Vp_h.begin(), atom_Vp_h.begin()+n_atoms, atom_Vp.begin());
			}
			
			void get_cubic_poly_coef_Vz(Stream<Dev>& stream, Vctr<Ptc_pVp<T>, edev_cpu>& atom_Vp)
			{
				if (this->multem_in_parm->is_spec_slic_by_dz_sub())
				{
					mt::linear_Vz(stream, this->multem_in_parm->atomic_pot_parm_typ, qz, atom_Vp);
					mt::fcn_vd_2_coef_poly3(stream, atom_Vp);
				}
			}

			Quad_Coef_1d<T, Dev> qz;
			Vctr<Atomic_Info_1<T, Dev>, edev_cpu> atom_type;		// Atom types

			Stream_Data stream_data;
			Vctr<Ptc_pVp<T>, edev_cpu> atom_Vp_h;
			Vctr<Ptc_pVp<T>, Dev> atom_Vp;
	};

}
#endif