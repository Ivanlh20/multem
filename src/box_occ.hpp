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

#ifndef BOX_OCC_H
#define BOX_OCC_H

#include "types.cuh"
#include "math.cuh"
#include "r3d.cuh"
#include "atom_data.hpp"

namespace multem
{
	class Box_Occ
	{
		public:
			Box_Occ(): r_min(0), r2_min(0), a_min(0), lx(0), ly(0), lz(0), nx(0), ny(0), nz(0), nxy(0){};

			void init()
			{
				for(auto i = 0; i<occ.size(); i++)
				{
					if(occ[i]>-2)
					{
						occ[i] = -1;
					}
				}
			}

			void set_input_data(double r_min_i, double lx_i, double ly_i, double lz_i)
			{	
				r_min = r_min_i;
				r2_min = pow(r_min, 2);
				a_min = r_min/multem::c_3i2;
				lx = lx_i;
				ly = ly_i;
				lz = lz_i;
				nx = static_cast<int>(ceil(lx/a_min));
				ny = static_cast<int>(ceil(ly/a_min));
				nz = static_cast<int>(ceil(lz/a_min));
				nxy = nx*ny;
				occ.resize(nxy*nz, -1);
			}

			inline
			int64_t xyz_2_ind(const int &ix, const int &iy, const int &iz)
			{
				return int64_t(iz)*nxy+ int64_t(iy)*nx + int64_t(ix);
			}

			template<class T>
			inline void set_occ(const r3d<T> &r, const int &iatom)
			{
				int ix = static_cast<int>(floor(r.x/a_min));
				int iy = static_cast<int>(floor(r.y/a_min));
				int iz = static_cast<int>(floor(r.z/a_min));
				auto ixyz = xyz_2_ind(ix, iy, iz);
				occ[ixyz] = iatom;
			}

			template<class T>
			inline bool check_r_min(const Atom_SA<T> &atoms, const r3d<T> &r_i)
			{
				int ix_i = static_cast<int>(floor(r_i.x/a_min));
				int iy_i = static_cast<int>(floor(r_i.y/a_min));
				int iz_i = static_cast<int>(floor(r_i.z/a_min));
				auto ixyz = xyz_2_ind(ix_i, iy_i, iz_i);

				if (occ[ixyz]>-1)
				{
					return false;
				}

				auto get_limits = [](const int &k, const int &nk, int &k_0, int &k_e)->void
				{
					k_0 = (k-1<0)?0:(k-1);
					k_e = (k+2>nk)?nk:(k+2);
				};

				int ix_0, ix_e;
				get_limits(ix_i, nx, ix_0, ix_e);
				int iy_0, iy_e;
				get_limits(iy_i, ny, iy_0, iy_e);
				int iz_0, iz_e;
				get_limits(iz_i, nz, iz_0, iz_e);

				for(auto iz =iz_0; iz<iz_e; iz++)
				{
					for(auto iy =iy_0; iy<iy_e; iy++)
					{
						for(auto ix =ix_0; ix<ix_e; ix++)
						{
							auto iatom = occ[xyz_2_ind(ix, iy, iz)];
							if((iatom>-1)&&(norm(atoms.r_n[iatom]-r_i)<r2_min))
							{
								return false;
							}
						}
					}
				}

				return true;
			}

		private:
			double r_min;
			double r2_min;
			double a_min;
			double lx;
			double ly;
			double lz;
			int nx;
			int ny;
			int nz;
			int64_t nxy;

			Vector<int, e_host> occ;
	};

} // namespace multem
#endif