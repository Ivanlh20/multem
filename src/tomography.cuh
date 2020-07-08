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

#ifndef TOMOGRAPHY_H
#define TOMOGRAPHY_H

#include <random>
#include <algorithm>

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "lin_alg_def.cuh"
#include "stream.cuh"
#include "atomic_data.hpp"
#include "input_tomography.cuh"
#include "output_tomography.hpp"
#include "box_occ.hpp"
#include "cubic_spline.hpp"

namespace mt
{
	template <class T, eDevice dev>
	class Tomography
	{
		public:
			using value_type = T;
			using size_type = std::size_t;

			static const eDevice device = dev;

			Tomography(): input_tomography(nullptr), rand_trial(0), chi2_0(0), Temp_max(0), 
			Temp_min(0), rTemp(0), chi2(0), chi2_n(0), chi2_opt(0){}

			void set_input_data(Input_Tomography<T> *input_tomography_i)
			{	
				input_tomography = input_tomography_i;
				atoms.assign(input_tomography->atoms);
				box.set_input_data(input_tomography->r0_min, input_tomography->grid_2d.lx, input_tomography->grid_2d.ly, input_tomography->grid_2d.lx);

				image.resize(input_tomography->image.size());
				assign_input_image();

				chi2_0 = 0;
				T IA_exp = 0;
				for(auto irot = 0; irot< image.size(); irot++)
				{
					chi2_0 += mt::sum_square(input_tomography->grid_2d, image[irot]);
					IA_exp += mt::sum(input_tomography->grid_2d, image[irot]);
				}
				IA_exp = IA_exp*input_tomography->grid_2d.dRx*input_tomography->grid_2d.dRy/(atoms.size()*image.size());

				T IA_sim = 0;
				for(auto ir = 0; ir<input_tomography->r.size(); ir++)
				{
					T &r = input_tomography->r[ir];
					T &fr = input_tomography->fr[ir];
					IA_sim += r*fr;
				}
				IA_sim = IA_sim*(input_tomography->r[1]-input_tomography->r[0]);

				T IA_factor = IA_exp/(c_2Pi*IA_sim);

				for(auto ir = 0; ir<input_tomography->r.size(); ir++)
				{
					T &fr = input_tomography->fr[ir];
					fr = fr*IA_factor;
				}

				/***************************************************************************/
				Atomic_Data atomic_data;
				atomic_data.Load_Data(ePT_Lobato_0_12);

				atom_type_host.resize(c_nAtomsTypes); 
				atom_type.resize(c_nAtomsTypes); 
				mt::Cubic_Spline spline;
				for(auto i = 0; i<c_nAtomsTypes; i++)
				{
					auto Z = i+1;
					atomic_data.To_atom_type_CPU(Z, c_Vrl, c_nR, input_tomography->grid_2d.dR_min(), atom_type_host[i]);
					if(input_tomography->Z == Z)
					{
						mt::assign(input_tomography->r, atom_type_host[i].R);
						mt::square(input_tomography->r, atom_type_host[i].R2);
						spline.set_points(atom_type_host[i].R2, input_tomography->fr);
						spline.get_coeff(atom_type_host[i].ciVR);
						atom_type_host[i].R_min = 0;
						atom_type_host[i].R2_min = pow(atom_type_host[i].R_min, 2);
						atom_type_host[i].R_max = atom_type_host[i].R.back();
						atom_type_host[i].R2_max = pow(atom_type_host[i].R_max, 2);
					}
					atom_type[i].assign(atom_type_host[i]);
				}

				/***************************************************************************/
				stream.resize(input_tomography->nstream);

				int nv = max(input_tomography->grid_2d.nx_dRx(2*input_tomography->grid_2d.lx), input_tomography->grid_2d.ny_dRy(2*input_tomography->grid_2d.ly));
				stream_data.resize(stream.size());
				for(auto i = 0; i<stream_data.size(); i++)
				{
					stream_data.iv[i].resize(nv);
					stream_data.v[i].resize(nv);
				}

				atom_Ip.resize(stream.size());

				rand_trial = 1000;
			}

			void init_point()
			{
				int ntrial = 100;

				chi2 = chi2_0;
				set_atoms_range();
				for(auto itrial = 0; itrial<ntrial; itrial++)
				{
					move_atoms();
					chi2_n = cost_function();
					if(chi2_n < chi2)
					{
						chi2 = chi2_n;
						atoms.r = atoms.r_n;
					}
				}

				chi2_opt = chi2;
				atoms.r_opt = atoms.r;

				Temp_max = 2*chi2_opt;
				Temp_min = 1e-06*chi2_opt;
			}

			template <class Output_Tomography>
			void run(Output_Tomography &output_tomography)
			{
				init_point();

				// T Temp = Temp_max;
				// while (Temp>Temp_min)
				// {
				// 	for(auto itemp = 0; itemp<ntemp; itemp++)
				// 	{
				// 		T p = 0;
				// 		for(auto ip = 0; ip<np; ip++)
				// 		{
				// 			move_atoms();
				// 			chi2_n = cost_function();
				// 			T delta_chi2 = chi2_n-chi2;
				// 			if(delta_chi2<0)
				// 			{
				// 				p = p + 1;
				// 				chi2 = chi2_n;
				// 				atoms.r = atoms.r_n;
				// 				set_atom_range();	

				// 				if(chi2_n < chi2_opt)
				// 				{
				// 					chi2_opt = chi2_n;
				// 					r_opt = r_n;
				// 				}
				// 			}
				// 			else if(exp(-delta_chi2/Temp) > rand.temp())
				// 			{
				// 				p = p + 1;
				// 				chi2 = chi2_n;
				// 				atoms.r = atoms.r_n;
				// 				set_atom_range();			
				// 			}
				// 		}

				// 		for(int iatoms = 0; iatoms<atoms.size(); iatoms++)
				// 		{
				// 			atoms.df[iatoms] = fmin(atoms.df[iatoms]*range_factor(p/np), 1.0);
				// 		}
				// 	}

				// 	Temp = rTemp*Temp;
				// 	chi2 = chi2_opt;
				// 	atoms.r = atoms.r_opt;
				// }

				for(auto iatoms = 0; iatoms<atoms.size(); iatoms++)
				{
					output_tomography.Z[iatoms] = atoms.Z[iatoms];
					auto r = atoms.r[iatoms];
					output_tomography.x[iatoms] = r.x;
					output_tomography.y[iatoms] = r.y;
					output_tomography.z[iatoms] = r.z;
				}
			}

		private:
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
				}

				Vector<Vector<int, dev>, e_host> iv;
				Vector<Vector<T, dev>, e_host> v;
			};

			void set_atom_Ip(int Z_i, const r3d<T> &r_i, Stream<dev> &stream, Vector<Atom_Sa<T>, e_host> &atom_Ip)
			{
				for(auto istream = 0; istream < stream.n_act_stream; istream++)
				{
					int iZ = Z_i-1;
					atom_Ip[istream].x = r_i.x;
					atom_Ip[istream].y = r_i.y;
					atom_Ip[istream].R2_max = atom_type[iZ].R2_max;
					atom_Ip[istream].R2 = raw_pointer_cast(atom_type[iZ].R2.data());
					atom_Ip[istream].set_ix0_ixn(input_tomography->grid_2d, atom_type[iZ].R_max);
					atom_Ip[istream].set_iy0_iyn(input_tomography->grid_2d, atom_type[iZ].R_max);
					atom_Ip[istream].c0 = raw_pointer_cast(atom_type[iZ].ciVR.c0.data());
					atom_Ip[istream].c1 = raw_pointer_cast(atom_type[iZ].ciVR.c1.data());
					atom_Ip[istream].c2 = raw_pointer_cast(atom_type[iZ].ciVR.c2.data());
					atom_Ip[istream].c3 = raw_pointer_cast(atom_type[iZ].ciVR.c3.data());
					atom_Ip[istream].iv = raw_pointer_cast(stream_data.iv[istream].data());
					atom_Ip[istream].v = raw_pointer_cast(stream_data.v[istream].data());

				}
			}

			void set_atom_range(const r3d<T> &r_i, const T &df, const r3d<T> &r_min, const r3d<T> &r_max, r3d<T> &r_0, r3d<T> &r_d)
			{
				r3d<T> dr = (r_max-r_min)*df;
				r_0 = fmax(r_min, r_i-dr);
				r_d = fmin(r_max, r_i+dr)-r_0;
			}

			r3d<T> move_atom(const r3d<T> &r_0, const r3d<T> &r_d)
			{
				for(auto itrial = 0; itrial < rand_trial; itrial++)
				{
					auto r = r_0+r_d*rand();
					if(box.check_r_min(atoms, r))
					{
						return r;
					}
				}
				return r_0;
			}

			T atom_cost_function(const r3d<T> &r_i)
			{
				stream.set_n_act_stream(1);
				T chi2 = 0;
				for(auto irot = 0; irot<input_tomography->angle.size(); irot++)
				{
					auto Rm = get_rotation_matrix(input_tomography->angle[irot], input_tomography->spec_rot_u0);
					auto r = r_i.rotate(Rm, input_tomography->spec_rot_center_p);
					set_atom_Ip(input_tomography->Z, r, stream, atom_Ip);
					chi2 += mt::atom_cost_function(input_tomography->grid_2d, atom_Ip[0], image[irot]);
				}
				return chi2;
			}

			void subtract_atom_in_all_rot(const r3d<T> &r_i)
			{
				stream.set_n_act_stream(1);
				for(auto irot = 0; irot<input_tomography->angle.size(); irot++)
				{
					auto Rm = get_rotation_matrix(input_tomography->angle[irot], input_tomography->spec_rot_u0);
					auto r = r_i.rotate(Rm, input_tomography->spec_rot_center_p);
					set_atom_Ip(input_tomography->Z, r, stream, atom_Ip);
					mt::subtract_atom(stream, input_tomography->grid_2d, atom_Ip, image[irot]);
				}
			}

			void set_atoms_range()
			{
				for(int iatoms = 0; iatoms<atoms.size(); iatoms++)
				{
					set_atom_range(atoms.r[iatoms], atoms.df[iatoms], atoms.r_min[iatoms], atoms.r_max[iatoms], atoms.r_0[iatoms], atoms.r_d[iatoms]);
				}
			}

			void move_atoms()
			{
				box.init();
				for(int iatoms = 0; iatoms<atoms.size(); iatoms++)
				{
					auto r = move_atom(atoms.r_0[iatoms], atoms.r_d[iatoms]);
					box.set_occ(r, iatoms);
					atoms.r_n[iatoms] = r;
				}
			}

			T cost_function()
			{
				stream.set_n_act_stream(1);
				assign_input_image();
				T chi2 = 0;
				for(auto irot = 0; irot<input_tomography->angle.size(); irot++)
				{
					auto Rm = get_rotation_matrix(input_tomography->angle[irot], input_tomography->spec_rot_u0);
					for(int iatoms = 0; iatoms<atoms.size(); iatoms++)
					{
						auto r = atoms.r_n[iatoms].rotate(Rm, input_tomography->spec_rot_center_p);
						// r = r.rotate(Rm, input_tomography->spec_rot_center_p);
						set_atom_Ip(input_tomography->Z, r, stream, atom_Ip);
						mt::subtract_atom(stream, input_tomography->grid_2d, atom_Ip, image[irot]);
					}
					chi2 += mt::sum_square(input_tomography->grid_2d, image[irot]);
				}
				return chi2;
			}

			void assign_input_image()
			{
				for(auto irot = 0; irot< image.size(); irot++)
				{
					mt::assign(input_tomography->image[irot], image[irot]);
				}
			}

			inline T range_factor(const T &p)
			{
				T m = 2.0;
				T p1 = 0.6;
				T p2 = 0.4;

				if(p>p1)
				{
					return 1.0 + m*(p-p1)/p2;
				}
				else if(p<p2)
				{
					return 1.0/(1.0 + m*(p2-p)/p2);
				}
				else
				{
					return 1.0;
				}
			}

			Input_Tomography<T> *input_tomography;

			Stream<dev> stream;

			Box_Occ<T> box;
			Rand_3d<T, e_host> rand;
			int rand_trial;
			T chi2_0;

			T Temp_max;
			T Temp_min;
			T rTemp;

			T chi2;
			T chi2_n;
			T chi2_opt;

			Atom_Data_Sa<T> atoms;
			Vector<Vector<T, dev>, e_host> image;
			Vector<Atom_Sa<T>, e_host> atom_Ip;
			Stream_Data stream_data;

			Vector<Atom_Type<T, e_host>, e_host> atom_type_host; // Atom types
			Vector<Atom_Type<T, dev>, e_host> atom_type;			// Atom types
	};

} // namespace mt
#endif