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

#ifndef BOX_OCC_H
#define BOX_OCC_H

#include <vector>
#include <deque>

#include "safe_types.cuh"
#include "math.cuh"
#include "lin_alg_def.cuh"
#include "atom_data_api.h"

namespace mt
{
	template <class T>
	class Box_Occ
	{
		public:
			Box_Occ(): d_min(0), d2_min(0), a_min(0), l_x(0), l_y(0), l_z(0), nx(0), ny(0), nz(0), nxy(0){};

			void init()
			{
				std::fill(occ.begin(), occ.end(), T(-1));
			}

			void set_input_data(T d_min_i, T lx_i, T ly_i, T lz_i)
			{	
				d_min = d_min_i;
				d2_min = pow(d_min, 2);
				a_min = d_min/mt::c_3i2;
				l_x = lx_i;
				l_y = ly_i;
				l_z = lz_i;
				nx = static_cast<int64_t>(ceil(l_x/a_min));
				ny = static_cast<int64_t>(ceil(l_y/a_min));
				nz = static_cast<int64_t>(ceil(l_z/a_min));
				nxy = nx*ny;
				occ.clear();
				occ.resize(nxy*nz, -1);
			}

			std::size_t xyz_2_ind(const int &ix, const int &iy, const int &iz) const
			{
				return (int64_t(iz)*nxy+ int64_t(iy)*nx + int64_t(ix));
			}

			int get_occ(const int &ix, const int &iy, const int &iz) const
			{
				return occ[xyz_2_ind(ix, iy, iz)];
			}

			void set_occ(const int &ix, const int &iy, const int &iz, const int &val)
			{
				occ[xyz_2_ind(ix, iy, iz)] = val;
			}

			void set_occ(const r3d<T> &r, const int &iatoms)
			{
				const int ix = static_cast<int>(floor(r.x/a_min));
				const int iy = static_cast<int>(floor(r.y/a_min));
				const int iz = static_cast<int>(floor(r.z/a_min));
				set_occ(ix, iy, iz, iatoms);
			}

			std::vector<int> range_loop(const int &k, const int &n_k)
			{
				if(k==0)
				{
					return std::vector<int>{k, k+1, k+2};
				}
				else if(k==1)
				{
					return std::vector<int>{k-1, k, k+1, k+2};
				}
				else if(k==n_k-1)
				{
					return std::vector<int>{k-2, k-1, k};
				}
				else if(k==n_k-2)
				{
					return std::vector<int>{k-2, k-1, k, k+1};
				}
				else
				{
					return std::vector<int>{k-2, k-1, k, k+1, k+2};
				}
			}

			std::vector<int> range_loop_pbc(const int &k, const int &n_k)
			{
				if(k==0)
				{
					return std::vector<int>{n_k-2, n_k-1, k, k+1, k+2};
				}
				else if(k==1)
				{
					return std::vector<int>{n_k-1, k-1, k, k+1, k+2};
				}
				else if(k==n_k-1)
				{
					return std::vector<int>{k-2, k-1, k, 0, 1};
				}
				else if(k==n_k-2)
				{
					return std::vector<int>{k-2, k-1, k, k+1, 0};
				}
				else
				{
					return std::vector<int>{k-2, k-1, k, k+1, k+2};
				}
			}	 

			template <class TAtom>
			bool check_r_min(TAtom &atoms, const r3d<T> &r)
			{
				const int ix_i = static_cast<int>(floor(r.x/a_min));
				const int iy_i = static_cast<int>(floor(r.y/a_min));
				const int iz_i = static_cast<int>(floor(r.z/a_min));

				if (get_occ(ix_i, iy_i, iz_i)>-1)
				{
					return false;
				}

				auto ax = range_loop_pbc(ix_i, nx);
				auto ay = range_loop_pbc(iy_i, ny);
				auto az = range_loop(iz_i, nz);

				for(auto iz: az)
				{
					for(auto iy: ay)
					{
						for(auto ix: ax)
						{
							auto iatoms = get_occ(ix, iy, iz);
							if(iatoms>-1)
							{
								auto d2 = atoms.norm_pbc_xy(iatoms, r);
								if(d2<d2_min)
								{
									return false;
								}
							}
						}
					}
				}

				return true;
			}

			T d_min;
			T d2_min;
			T a_min;
			T l_x;
			T l_y;
			T l_z;
			int64_t nx;
			int64_t ny;
			int64_t nz;
			int64_t nxy;
      std::vector<int> occ;
	};

} // namespace mt
#endif

