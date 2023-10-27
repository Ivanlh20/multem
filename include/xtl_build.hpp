/*
 * This file is part of Multem.
 * Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef XTL_BUILD_H
	#define XTL_BUILD_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include <vector>
	#include <cstdlib>
	#include <string>

	#include "types.cuh"
	#include "particles.cuh"
	#include "xtl_build_in_parm.hpp"
	#include "space_group.hpp"

	namespace mt
	{
		// from space group number to xtl system number
		dt_int32 xtl_sgn_2_csn(const dt_int32& sgn)
		{
			if (fcn_chk_bound(sgn, 1, 3))					// triclinic or anorthic
			{
				return 1;
			}
			else if (fcn_chk_bound(sgn, 3, 16))				// monoclinic
			{
				return 2;
			}
			else if (fcn_chk_bound(sgn, 16, 75))			// orthorhombic
			{
				return 3;
			}
			else if (fcn_chk_bound(sgn, 75, 143))			// tetragonal
			{
				return 4;
			}
			else if (fcn_chk_bound(sgn, 143, 168))			// rhombohedral or trigonal
			{
				return 5;
			}
			else if (fcn_chk_bound(sgn, 168, 195))			// hexagonal
			{
				return 6;
			}
			else if (fcn_chk_bound(sgn, 195, 231))			// cubic
			{
				return 7;
			}

			return 1;
		}

		// from xtl system string to xtl system number
		dt_int32 xtl_css_2_csn(std::string cs_str)
		{
			fcn_str_rtrim_ip(cs_str);

			cs_str = mt::fcn_str_2_lower(cs_str);
			//std::transform(cs_str.begin(), cs_str.end(), cs_str.begin(), ::tolower);

			if (fcn_str_cmp(cs_str, "triclinic")||fcn_str_cmp(cs_str, "a"))				// triclinic or anorthic
			{
				return 1;
			}
			else if (fcn_str_cmp(cs_str, "monoclinic")||fcn_str_cmp(cs_str, "m")) 		// monoclinic
			{
				return 2;
			}
			else if (fcn_str_cmp(cs_str, "orthorhombic")||fcn_str_cmp(cs_str, "o")) 	// orthorhombic
			{
				return 3;
			}
			else if (fcn_str_cmp(cs_str, "tetragonal")||fcn_str_cmp(cs_str, "t")) 		// tetragonal
			{
				return 4;
			}
			else if (fcn_str_cmp(cs_str, "rhombohedral")||fcn_str_cmp(cs_str, "r")) 	// rhombohedral or trigonal
			{
				return 5;
			}
			else if (fcn_str_cmp(cs_str, "hexagonal")||fcn_str_cmp(cs_str, "h")) 		// hexagonal
			{
				return 6;
			}
			else if (fcn_str_cmp(cs_str, "cubic")||fcn_str_cmp(cs_str, "c")) 			// cubic
			{
				return 7;
			}

			return 1;
		}

		// from xtl system number to space group range
		std::pair<dt_int32, dt_int32> xtl_csn_2_sgr(const dt_int32& csn)
		{
			switch(csn)				
			{
				case 1:						// triclinic or anorthic
					return {1, 2};
				case 2:						// monoclinic
					return {3, 15};
				case 3:						// orthorhombic
					return {16, 74};
				case 4:						// tetragonal
					return {75, 142};
				case 5:						// rhombohedral or trigonal
					return {143, 167};
				case 6:						// hexagonal
					return {168, 194};
				case 7:						// cubic
					return {195, 230};
			}

			return {0, 0};
		}

		// from xtl system string to space group range
		std::pair<dt_int32, dt_int32> xtl_css_2_sgr(const std::string& cs_str)
		{
			auto csn = xtl_css_2_csn(cs_str);

			return xtl_csn_2_sgr(csn);
		}

		template <class T>
		class Xtl_Build{
			public:
				Xtl_Build(): n_a(0), n_b(0), n_c(0), a(0), b(0), c(0) {};

				Xtl_Build(const Xtl_Build_In_Parm<T>& xtl_build_in_parm)
				{
					set_in_data(xtl_build_in_parm);
				}

				void set_in_data(const Xtl_Build_In_Parm<T>& xtl_build_in_parm)
				{
					a = xtl_build_in_parm.a;
					b = xtl_build_in_parm.b;
					c = xtl_build_in_parm.c;

					alpha = xtl_build_in_parm.alpha;
					beta = xtl_build_in_parm.beta;
					gamma = xtl_build_in_parm.gamma;

					n_a = xtl_build_in_parm.n_a;
					n_b = xtl_build_in_parm.n_b;
					n_c = xtl_build_in_parm.n_c;

					sgn = xtl_build_in_parm.sgn;
					pbc = xtl_build_in_parm.pbc;

					asym_uc = xtl_build_in_parm.asym_uc;
					base = xtl_build_in_parm.base;

					if (xtl_build_in_parm.base.size()==0)
					{
						base = space_group(asym_uc, sgn);
					}
				}
			
				Ptc_Atom<T> operator()() const
				{
					const auto dsm = direct_struct_metric();
					const R_3d<T> box_n = R_3d<T>(T(n_a), T(n_b), T(n_c)) + T(1e-4);

					const dt_int32 n_base = base.size();

					dt_int32 nt_a = n_a;
					dt_int32 nt_b = n_b;
					dt_int32 nt_c = n_c;

					if (!pbc)
					{
						nt_a++;
						nt_b++;
						nt_c++;
					}
					else
					{
						nt_a = max(1, nt_a);
						nt_b = max(1, nt_b);
						nt_c = max(1, nt_c);
					}

					Ptc_Atom<T> xtl;
					xtl.reserve(nt_a*nt_b*nt_c*n_base);
					xtl.cols_used = base.cols_used;

					for(auto ic = 0; ic < nt_c; ic++)
					{
						for(auto ib = 0; ib < nt_b; ib++)
						{
							for(auto ia = 0; ia < nt_a; ia++)
							{
								for(auto ik = 0; ik < n_base; ik++)
								{
									auto atom = base.get(ik);

									atom.x += T(ia);
									atom.y += T(ib);
									atom.z += T(ic);

									if ((atom.x<box_n.x) && (atom.y<box_n.y) && (atom.z<box_n.z))
									{
										atom.set_pos(dsm*atom.get_pos());

										xtl.push_back(atom);
									}
								}
							}
						}
					}
					xtl.shrink_to_fit();

					return xtl;
				}

				Mx_3x3<T> direct_struct_metric() const
				{
					T cos_alpha = cos(alpha);
					T cos_beta = cos(beta);
					T cos_gamma = cos(gamma);
 
					T sin_gamma = sin(gamma);
 
					T F_bga = cos_beta*cos_gamma - cos_alpha;
 
					T vol2 = ::square(a*b*c)*(1-::square(cos_alpha)-::square(cos_beta)-::square(cos_gamma)+2*cos_alpha*cos_beta*cos_gamma);
					T vol = ::sqrt(vol2);
 
					return {a, 0, 0, b*cos_gamma, b*sin_gamma, 0, c*cos_beta, -c*F_bga/sin_gamma, vol/(a*b*sin_gamma)};
				}

			private:

				T a;		// lattice constant a
				T b;		// lattice constant b
				T c;		// lattice constant c

				T alpha;		// angle between b & c
				T beta;		// angle between c & a
				T gamma;		// angle between a & b

				dt_int32 n_a;
				dt_int32 n_b;
				dt_int32 n_c;

				dt_int32 sgn;		// space group number
				dt_bool pbc;

				Ptc_Atom<T> asym_uc;		// normalized positions in the asymmetric unit cell
				Ptc_Atom<T> base;		// normalized positions of atom in the unit cell
				Space_Group<T> space_group;
		};

	}

#endif
