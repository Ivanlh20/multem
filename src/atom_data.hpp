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

#ifndef ATOM_DATA_H
#define ATOM_DATA_H

#include <numeric>
#include <algorithm>

#include "math.cuh"
#include "types.hpp"

namespace multem
{
	template<class T>
	class Atom_Data{
		public:
			using value_type = typename T;
			using size_type = std::size_t;

			Atom_Data(): l_x(0), l_y(0), l_z(0), 
				Z_min(0), Z_max(0), x_min(0), x_max(0), 
				y_min(0), y_max(0), z_min(0), z_max(0),
				sigma_min(0), sigma_max(0), occ_min(0), 
				occ_max(0), R_Int_min(0), R_Int_max(0), 
				x_mean(0), y_mean(0), z_mean(0), x_std(0), 
				y_std(0), z_std(0), s_x(0), s_y(0),s_z(0),
				x_Int_min(0), x_Int_max(0), y_Int_min(0),
				y_Int_max(0), z_Int_min(0), z_Int_max(0),
				s_x_Int(0), s_y_Int(0), s_z_Int(0) {}

			size_type size() const
			{
				return Z.size();
			}

			// resize number of atoms
			void resize(const size_type &new_size, const value_type &value = value_type())
			{
				Z.resize(new_size, value);
				Z.shrink_to_fit();

				x.resize(new_size, value);
				x.shrink_to_fit();

				y.resize(new_size, value);
				y.shrink_to_fit();

				z.resize(new_size, value);
				z.shrink_to_fit();

				sigma.resize(new_size, static_cast<float>(value));
				sigma.shrink_to_fit();

				occ.resize(new_size, static_cast<float>(value));
				occ.shrink_to_fit();
			}

			// set atoms
			void set_Atoms(const size_type &natomsM_i, double *atomsM_i, T l_x_i=0, T l_y_i=0, T l_z_i=0, T a_i=1.0, T b_i=1.0, T c_i=1.0)
			{
				resize(natomsM_i);

				l_x = l_x_i;
				l_y = l_y_i;
				l_z = l_z_i;

				for(auto i = 0; i < natomsM_i; i++)
				{		
					Z[i] = static_cast<int>(atomsM_i[0*natomsM_i + i]); 		// Atomic number
					x[i] = a_i*atomsM_i[1*natomsM_i + i]; 						// x-position
					y[i] = b_i*atomsM_i[2*natomsM_i + i]; 						// y-position
					z[i] = c_i*atomsM_i[3*natomsM_i + i]; 						// z-position
					sigma[i] = static_cast<float>(atomsM_i[4*natomsM_i + i]);	// Standard deviation
					occ[i] = static_cast<float>(atomsM_i[5*natomsM_i + i]); 	// Occupancy
				}

				get_Statistic();
			}

			// set atoms
			void set_Atoms(const Atom_Data<T> &atoms, bool PBC_xy_i=false, Vector<Atom_Type<T, Host>, Host> *atom_type=0)
			{
				resize(atoms.size());

				l_x = atoms.l_x;
				l_y = atoms.l_y;

				T dl = 1e-04;
				T lx_b = l_x - dl;
				T ly_b = l_y - dl;

				size_type j = 0;
				for(auto i = 0; i < size(); i++)
				{
					if((!PBC_xy_i)||((atoms.x[i]<lx_b)&&(atoms.y[i]<ly_b)))
					{
						Z[j] = atoms.Z[i];
						x[j] = atoms.x[i];
						y[j] = atoms.y[i];
						z[j] = atoms.z[i];
						sigma[j] = atoms.sigma[i];
						occ[j] = atoms.occ[i];
						j++;
					}
				}

				resize(j);

				get_Statistic(atom_type);
			}
		
			// get statistic
			void get_Statistic(Vector<Atom_Type<T, Host>, Host> *atom_type=0)
			{
				if(Z.empty())
				{
					return;
				}

				bool bAtomTypes = (atom_type == 0)?false:true;

				Z_min = Z_max = Z[0];

				x_min = x_max = x[0];
				y_min = y_max = y[0];
				z_min = z_max = z[0];

				sigma_min = sigma_max = sigma[0];
				occ_min = occ_max = occ[0];

				R_Int_min = R_Int_max = (bAtomTypes)?(*atom_type)[Z[0]-1].R_max:2.0;

				x_mean = y_mean = z_mean = 0.0;
				x_std = y_std = z_std = 0.0;

				Z_unique.resize(c_nAtomsTypes);
				std::fill(Z_unique.begin(), Z_unique.end(), 0);

				for(auto iAtom = 0; iAtom < size(); iAtom++)
				{
					Z_unique[Z[iAtom]]++; 

					Z_min = min(Z[iAtom], Z_min);
					Z_max = max(Z[iAtom], Z_max);

					x_min = min(x[iAtom], x_min);
					x_max = max(x[iAtom], x_max);

					y_min = min(y[iAtom], y_min);
					y_max = max(y[iAtom], y_max);

					z_min = min(z[iAtom], z_min);
					z_max = max(z[iAtom], z_max);

					sigma_min = min(sigma[iAtom], sigma_min);
					sigma_max = max(sigma[iAtom], sigma_max);

					occ_min = min(occ[iAtom], occ_min);
					occ_max = max(occ[iAtom], occ_max);

					if(bAtomTypes)
					{
						R_Int_min = min((*atom_type)[Z[iAtom]-1].R_max, R_Int_min);
						R_Int_max = max((*atom_type)[Z[iAtom]-1].R_max, R_Int_max);
					}

					x_mean += x[iAtom];
					y_mean += y[iAtom];
					z_mean += z[iAtom];

					x_std += x[iAtom]*x[iAtom];
					y_std += y[iAtom]*y[iAtom];
					z_std += z[iAtom]*z[iAtom];
				}

				auto it = std::remove_if(Z_unique.begin(), Z_unique.end(), [](const int &Z){return Z==0;} );
				Z_unique.resize(std::distance(Z_unique.begin(),it));

				T nAtoms = static_cast<T>(size());

				x_mean /= nAtoms;
				y_mean /= nAtoms;
				z_mean /= nAtoms;

				x_std = sqrt(x_std/nAtoms - x_mean*x_mean);
				y_std = sqrt(y_std/nAtoms - y_mean*y_mean);
				z_std = sqrt(z_std/nAtoms - z_mean*z_mean);

				s_x = x_max - x_min;
				s_y = y_max - y_min;
				s_z = z_max - z_min;

				x_Int_min = x_min - R_Int_max;
				x_Int_max = x_max + R_Int_max;

				y_Int_min = y_min - R_Int_max;
				y_Int_max = y_max + R_Int_max;

				z_Int_min = z_min - R_Int_max;
				z_Int_max = z_max + R_Int_max;

				s_x_Int = x_Int_max - x_Int_min;
				s_y_Int = y_Int_max - y_Int_min;
				s_z_Int = z_Int_max - z_Int_min;

				if(isZero(l_x))
				{
					l_x = s_x;
				}

				if(isZero(l_y))
				{
					l_y = s_y;
				}

				if(isZero(l_z))
				{
					l_z = s_z;
				}

				l_x_Int = l_x + 2.0*R_Int_max;
				l_y_Int = l_y + 2.0*R_Int_max;
				l_z_Int = l_z + 2.0*R_Int_max;
			}

			// Sort atoms along z-axis.
			void Sort_by_z()
			{
				Vector<int, Host> idx_sort(size());
				Vector<int, Host> val_i(size());
				Vector<T, Host> val_d(size());

				// Sort atoms along z-axis.
				std::iota(idx_sort.begin(), idx_sort.end(), 0);
				std::sort(idx_sort.begin(), idx_sort.end(), [&](int &i, int &j) {return z[i] < z[j]; });

				std::transform(idx_sort.begin(), idx_sort.end(), val_i.begin(), [&](int &i) {return Z[i]; });
				Z.assign(val_i.begin(), val_i.end());

				std::transform(idx_sort.begin(), idx_sort.end(), val_d.begin(), [&](int &i) {return x[i]; });
				x.assign(val_d.begin(), val_d.end());

				std::transform(idx_sort.begin(), idx_sort.end(), val_d.begin(), [&](int &i) {return y[i]; });
				y.assign(val_d.begin(), val_d.end());

				std::transform(idx_sort.begin(), idx_sort.end(), val_d.begin(), [&](int &i) {return z[i]; });
				z.assign(val_d.begin(), val_d.end());

				std::transform(idx_sort.begin(), idx_sort.end(), val_d.begin(), [&](int &i) {return sigma[i]; });
				sigma.assign(val_d.begin(), val_d.end());

				std::transform(idx_sort.begin(), idx_sort.end(), val_d.begin(), [&](int &i) {return occ[i]; });
				occ.assign(val_d.begin(), val_d.end());
			}

			// get Layers: Require that the atoms to be sorted along z
			void get_z_layer()
			{
				z_layer.clear();

				if(size() == 0 )
				{
					return;
				}

				T zb, zc, zm, zm_c;
				z_layer.reserve(size());
				zb = zm = z[0]; zm_c = 1.0;
				for(auto iAtom = 0; iAtom < size(); iAtom++)
				{
					zc = z[iAtom];
					if(abs(zb-zc)<Epsilon<T>::rel)
					{
						zm += zc;
						zm_c += 1.0;
					}
					else
					{
						z_layer.push_back(static_cast<float>(zm/zm_c));
						zm = zc;
						zm_c = 1.0;
					}
					zb = zc;
				} 
				z_layer.push_back(static_cast<float>(zm/zm_c));
				z_layer.shrink_to_fit();
			}

			// find atoms in slice
			void find_by_z(T z_0, T z_e, int &iz_0, int &iz_e, bool Inc_Borders=true)
			{
				z_0 = (Inc_Borders)?z_0-Epsilon<T>::rel:z_0;
				z_e = (Inc_Borders)?z_e+Epsilon<T>::rel:z_e;

				if((z_0>z_e)||(z.back()<z_0)||(z_e<z.front()))
				{
					iz_0 = 1;
					iz_e = 0;
					return;
				}

				iz_0 = (z_0 <= z.front())?0:std::lower_bound(z.begin(), z.end(), z_0)-z.begin();
				iz_e = (z.back() <= z_e)?z.size()-1:std::lower_bound(z.begin(), z.end(), z_e)-z.begin()-1;

				if((iz_0 > iz_e)||(z[iz_e] < z_0)||(z_e < z[iz_0]))
				{
					iz_0 = 1; 
					iz_e = 0;
				}
			}

			T l_x; 			// Box size-x
			T l_y; 			// Box size-y
			T l_z; 			// Box size-z

			Vector<int, Host> Z;
			Vector<T, Host> x;
			Vector<T, Host> y;
			Vector<T, Host> z;
			Vector<float, Host> sigma;
			Vector<float, Host> occ;

			Vector<int, Host> Z_unique;
			Vector<float, Host> z_layer;

			int Z_min;
			int Z_max;

			T x_min;
			T x_max;

			T y_min;
			T y_max;

			T z_min;
			T z_max;

			float sigma_min;
			float sigma_max;

			float occ_min;
			float occ_max;

			T R_Int_min;
			T R_Int_max;

			T x_mean;
			T y_mean;
			T z_mean;

			T x_std;
			T y_std;
			T z_std;

			T s_x; 			// size-x
			T s_y; 			// size-y
			T s_z; 			// size-z

			T x_Int_min;
			T x_Int_max;

			T y_Int_min;
			T y_Int_max;

			T z_Int_min;
			T z_Int_max;

			T s_x_Int;
			T s_y_Int;
			T s_z_Int;

			T l_x_Int;
			T l_y_Int;
			T l_z_Int;
	};

} // namespace multem

#endif