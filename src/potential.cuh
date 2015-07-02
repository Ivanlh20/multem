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

#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "math.cuh"
#include "types.hpp"
#include "quadrature.hpp"
#include "input_multislice.hpp"
#include "host_functions.hpp"
#include "device_functions.cuh"
#include "host_device_functions.cuh"
#include "specimen.hpp"

namespace multem
{
	template<class T, eDevice dev>
	class Potential: public Specimen<T, dev>{
		public:
			Potential(): stream(nullptr){ }

			void set_input_data(Input_Multislice<value_type_r, dev> *input_multislice_io, Stream<value_type_r, dev> *stream_i)
			{	
				Specimen::set_input_data(input_multislice_io);
				stream = stream_i;

				Quadrature quadrature;
				quadrature.get(0, c_nqz, qz); // 0: int_-1^1 y(x) dx - TanhSinh quadrature

				atom_type.resize(c_nAtomsTypes); 
				for(auto i=0; i<atom_type.size(); i++)
				{
					atom_type[i].assign(Specimen::atom_type[i]);
				}

				int nv = max(input_multislice->grid.nx_dRx(atoms.l_x_Int), input_multislice->grid.ny_dRy(atoms.l_y_Int));

				iv.resize(stream->size());
				v.resize(stream->size());
				ciV0.resize(stream->size());
				for(auto i=0; i<stream->size(); i++)
				{
					iv[i].resize(nv);
					v[i].resize(nv);
					ciV0[i].resize(c_nR);
				}

				atom_Vp.resize(stream->size());

				V0.resize(input_multislice->grid.nxy());
			}

			void projected_potential(const value_type_r &z_0, const value_type_r &z_e, const int &iatom_0, const int &iatom_e, Vector<value_type_r, dev> &V)
			{
				multem::fill(V, 0.0);

				int iatom = iatom_0;
				while (iatom <= iatom_e)
				{
					stream->n_act_stream = min(stream->size(), iatom_e-iatom+1);
					set_atom_Vp(z_0, z_e, iatom, *stream, atom_Vp);
					get_cubic_poly_coef_Vz(*stream, atom_Vp);
					multem::eval_cubic_poly(input_multislice->grid, *stream, atom_Vp, V);
					iatom += stream->n_act_stream;
				}

				stream->synchronize();
			}

			void projected_potential(const int &islice_0, const int &islice_e, Vector<value_type_r, dev> &V)
			{
				if((islice_0<0) || (islice_e>=slice.size()))
				{
					multem::fill(V, 0.0);
					return;
				}
				projected_potential(slice.z_0[islice_0], slice.z_e[islice_e], slice.iatom_0[islice_0], slice.iatom_e[islice_e], V);
			}

			void projected_potential(const int &islice_0, const int &islice_e)
			{
				projected_potential(islice_0, islice_e, V0);
			}

			void projected_potential(const int &islice, Vector<value_type_r, dev> &V)
			{
				projected_potential(islice, islice, V);
			}

			void projected_potential(const int &islice)
			{
				projected_potential(islice, islice, V0);
			}

			Vector<value_type_r, dev> V0;
		private:

			void set_atom_Vp(const value_type_r &z_0, const value_type_r &z_e, int iatom, Stream<value_type_r, dev> &stream, Vector<Atom_Vp<value_type_r>, e_Host> &atom_Vp)
			{
				for(auto istream = 0; istream < stream.n_act_stream; istream++)
				{
					int iZ = atoms.Z[iatom]-1;
					atom_Vp[istream].x = atoms.x[iatom];
					atom_Vp[istream].y = atoms.y[iatom];
					atom_Vp[istream].occ = atoms.occ[iatom];
					atom_Vp[istream].R_min2 = atom_type[iZ].R_min2;
					atom_Vp[istream].R_max2 = atom_type[iZ].R_max2;
					atom_Vp[istream].R2 = raw_pointer_cast(atom_type[iZ].R2.data());
					atom_Vp[istream].set_ix0_ixn(input_multislice->grid, atom_type[iZ].R_max);
					atom_Vp[istream].set_iy0_iyn(input_multislice->grid, atom_type[iZ].R_max);
					atom_Vp[istream].iv = raw_pointer_cast(iv[istream].data());
					atom_Vp[istream].v = raw_pointer_cast(v[istream].data());

					if(input_multislice->is_subslicing())
					{
						atom_Vp[istream].z0h = 0.5*(z_0 - atoms.z[iatom]); 
						atom_Vp[istream].zeh = 0.5*(z_e - atoms.z[iatom]);
						atom_Vp[istream].split = (atom_Vp[istream].z0h<0)&&(0<atom_Vp[istream].zeh);
						atom_Vp[istream].cl = raw_pointer_cast(atom_type[iZ].Vr.cl.data());
						atom_Vp[istream].cnl = raw_pointer_cast(atom_type[iZ].Vr.cnl.data());
						atom_Vp[istream].c0 = raw_pointer_cast(ciV0[istream].c0.data());
						atom_Vp[istream].c1 = raw_pointer_cast(ciV0[istream].c1.data());
						atom_Vp[istream].c2 = raw_pointer_cast(ciV0[istream].c2.data());
						atom_Vp[istream].c3 = raw_pointer_cast(ciV0[istream].c3.data());
					}
					else
					{
						atom_Vp[istream].c0 = raw_pointer_cast(atom_type[iZ].ciVR.c0.data());
						atom_Vp[istream].c1 = raw_pointer_cast(atom_type[iZ].ciVR.c1.data());
						atom_Vp[istream].c2 = raw_pointer_cast(atom_type[iZ].ciVR.c2.data());
						atom_Vp[istream].c3 = raw_pointer_cast(atom_type[iZ].ciVR.c3.data());
					}
					iatom++;
				}
			}
			
			void get_cubic_poly_coef_Vz(Stream<value_type_r, dev> &stream, Vector<Atom_Vp<value_type_r>, e_Host> &atom_Vp)
			{
				if(input_multislice->is_subslicing())
				{
					multem::get_cubic_poly_coef_Vz(input_multislice->potential_type, qz, stream, atom_Vp);
				}
			}

			Vector<Atom_Type<value_type_r, dev>, e_Host> atom_type; // Atom types
			Stream<value_type_r, dev> *stream;

			Q1<value_type_r, dev> qz;
			Vector<Vector<int, dev>, e_Host> iv;
			Vector<Vector<value_type_r, dev>, e_Host> v;
			Vector<CI_Coef<value_type_r, dev>, e_Host> ciV0;

			Vector<Atom_Vp<value_type_r>, e_Host> atom_Vp;
	};

} // namespace multem
#endif