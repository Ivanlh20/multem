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

#ifndef ATOM_DATA_H
#define ATOM_DATA_H

#include <numeric>
#include <algorithm>

#include "math.cuh"
#include "types.cuh"
#include "lin_alg_def.cuh"

#include <thrust/sort.h>

namespace multem
{
	template<class T>
	class Atom_Data{
		public:
			using value_type = T;
			using size_type = std::size_t;

			Atom_Data(): l_x(0), l_y(0), l_z(0), 
				Z_min(0), Z_max(0), x_min(0), x_max(0), 
				y_min(0), y_max(0), z_min(0), z_max(0), 
				sigma_min(0), sigma_max(0), occ_min(0), 
				occ_max(0), R_Int_min(0), R_Int_max(0), 
				x_mean(0), y_mean(0), z_mean(0), x_std(0), 
				y_std(0), z_std(0), s_x(0), s_y(0), s_z(0), 
				x_Int_min(0), x_Int_max(0), y_Int_min(0), 
				y_Int_max(0), z_Int_min(0), z_Int_max(0), 
				s_x_Int(0), s_y_Int(0), s_z_Int(0), 
				l_x_Int(0), l_y_Int(0), l_z_Int(0){}

			size_type size() const
			{
				return Z.size();
			}

			bool empty() const
			{
				return size() == 0;
			}

			// resize number of atoms
			void resize(const size_type &new_size, const value_type &value = value_type())
			{
				Z.resize(new_size, value);
				Z.shrink_to_fit();

				charge.resize(new_size, value);
				charge.shrink_to_fit();

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


			template<class TOut, class TIn>
			Atom<TOut> read_atom(const int &nr, const int &nc, TIn *atoms, const int &iatom)
			{
				Atom<TOut> atom;
				atom.Z = static_cast<int>(atoms[0*nr + iatom]); 						// Atomic number
				atom.x = atoms[1*nr + iatom]; 										// x-position
				atom.y = atoms[2*nr + iatom]; 										// y-position
				atom.z = atoms[3*nr + iatom]; 										// z-position
				atom.sigma = static_cast<float>((nc>4)?(atoms[4*nr + iatom]):0.085);	// Standard deviation
				atom.occ = static_cast<float>((nc>5)?(atoms[5*nr + iatom]):1.0); 		// Occupancy
				atom.charge = static_cast<int>((nc>6)?(atoms[6*nr + iatom]):0);		// charge

				return atom;
			}

			// set atoms
			void set_Atoms(const size_type &nr_atoms_i, const size_type &nc_atoms_i, double *atoms_i, T l_x_i = 0, T l_y_i = 0, T l_z_i = 0, T a_i = 1.0, T b_i = 1.0, T c_i = 1.0)
			{
				resize(nr_atoms_i);

				l_x = l_x_i;
				l_y = l_y_i;
				l_z = l_z_i;

				for(auto i = 0; i < nr_atoms_i; i++)
				{	
					auto atom = read_atom<T>(nr_atoms_i, nc_atoms_i, atoms_i, i);
					Z[i] = atom.Z; 				// Atomic number
					charge[i] = atom.charge; 	// charge
					x[i] = a_i*atom.x; 			// x-position
					y[i] = b_i*atom.y; 			// y-position
					z[i] = c_i*atom.z; 			// z-position
					sigma[i] = atom.sigma;		// Standard deviation
					occ[i] = atom.occ; 			// Occupancy
				}

				get_Statistic();
			}

			// set atoms
			void set_Atoms(const Atom_Data<T> &atoms, bool PBC_xy_i =false, Vector<Atom_Type<T, e_host>, e_host> *atom_type = 0)
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
						charge[j] = atoms.charge[i];
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
			void get_Statistic(Vector<Atom_Type<T, e_host>, e_host> *atom_type_ptr =nullptr)
			{
				if(Z.empty())
				{
					return;
				}

				bool bAtomTypes = (atom_type_ptr == nullptr)?false:true;

				Z_min = Z_max = Z[0];

				x_min = x_max = x[0];
				y_min = y_max = y[0];
				z_min = z_max = z[0];

				sigma_min = sigma_max = sigma[0];
				occ_min = occ_max = occ[0];

				R_Int_min = R_Int_max = (bAtomTypes)?(*atom_type_ptr)[Z[0]-1].coef[0].R_max:2.5;

				x_mean = y_mean = z_mean = 0.0;
				x_std = y_std = z_std = 0.0;

				for(auto iAtom = 0; iAtom < size(); iAtom++)
				{
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
						R_Int_min = min((*atom_type_ptr)[Z[iAtom]-1].coef[0].R_max, R_Int_min);
						R_Int_max = max((*atom_type_ptr)[Z[iAtom]-1].coef[0].R_max, R_Int_max);
					}

					x_mean += x[iAtom];
					y_mean += y[iAtom];
					z_mean += z[iAtom];

					x_std += x[iAtom]*x[iAtom];
					y_std += y[iAtom]*y[iAtom];
					z_std += z[iAtom]*z[iAtom];
				}

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
				// if(!thrust::is_sorted(z.begin(), z.end()));
				auto first = thrust::make_zip_iterator(thrust::make_tuple(Z.begin(), charge.begin(), x.begin(), y.begin(), z.begin(), sigma.begin(), occ.begin()));
				auto last = thrust::make_zip_iterator(thrust::make_tuple(Z.end(), charge.end(), x.end(), y.end(), z.end(), sigma.end(), occ.end()));

				thrust::sort(first, last, sort_atoms_by_z());
			}

			// get Layers: Require that the atoms to be sorted along z
			void get_z_layer()
			{
				z_layer.clear();

				if(size() == 0)
				{
					return;
				}

				T zb, zc, zm, zm_c;
				z_layer.reserve(size());
				zb = zm = z[0]; zm_c = 1.0;
				for(auto iAtom = 0; iAtom < size(); iAtom++)
				{
					zc = z[iAtom];
					if(fabs(zb-zc) < Epsilon<float>::rel)
					{
						zm += zc;
						zm_c += 1.0;
					}
					else
					{
						z_layer.push_back(zm/zm_c);
						zm = zc;
						zm_c = 1.0;
					}
					zb = zc;
				} 
				z_layer.push_back(zm/zm_c);
				z_layer.shrink_to_fit();
			}

			// get spacing
			template<class TVector>
			T get_spacing(size_type ix, const TVector &x)
			{
				ix = (ix <= 0)?1:min(ix, x.size()-1);
				return (x.size()>1)?x[ix]-x[ix-1]:0.0;
			}

			// find atoms in slice
			int find_by_z(T z_e, bool Inc_Borders =true)
			{
				int iz_e =-1;
				z_e = (Inc_Borders)?z_e+Epsilon<T>::rel:z_e;

				if(z_e<z.front())
				{
					return iz_e;
				}

				iz_e = (z.back() <= z_e)?(z.size()-1):(std::lower_bound(z.begin(), z.end(), z_e)-z.begin()-1);

				if(z[iz_e] < z.front())
				{
					iz_e = -1;
				}

				return iz_e;
			}

			// find atoms in slice
			void find_by_z(T z_0, T z_e, int &iz_0, int &iz_e, bool Inc_Borders =true)
			{
				z_0 = (Inc_Borders)?(z_0-Epsilon<T>::rel):z_0;
				z_e = (Inc_Borders)?(z_e+Epsilon<T>::rel):z_e;

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

			// get z slicing
			void get_z_slice(const ePotential_Slicing &potential_slicing, const T &dz_i, const Atom_Data<T> &atoms, Vector<T, e_host> &z_slice)
			{
				z_slice.clear();

				auto quot = [&](const T &A, const T &B)->int
				{
					return (int)floor(A/B+Epsilon<float>::rel);
				};

				switch(potential_slicing)
				{
					case ePS_Planes:
					{
						if(z_layer.size() == 1)
						{
							z_slice.resize(2);
							z_slice[0] = z_Int_min;
							z_slice[1] = z_Int_max;
							z_slice.shrink_to_fit();
							return;
						}
						T dz_Bot = get_spacing(0, z_layer);
						T dz_Top = get_spacing(z_layer.size() - 1, z_layer);
						T layer_0 = z_layer.front() - 0.5*dz_Bot;
						T layer_e = z_layer.back() + 0.5*dz_Top;
						int nz_Bot = (atoms.z_min<layer_0)?quot(layer_0-atoms.z_min, dz_Bot) + 1:0;
						int nz_Top = (atoms.z_max>layer_e)?quot(atoms.z_max-layer_e, dz_Top) + 1:0;
						z_slice.resize(z_layer.size()+nz_Bot+nz_Top+1);
						int j = 0;
						for(auto i = 0; i <= nz_Bot; i++)
						{
							z_slice[j++] = layer_0-(nz_Bot-i)*dz_Bot;
						}
						for(auto i = 1; i<z_layer.size(); i++)
						{
							T dz = get_spacing(i, z_layer);
							z_slice[j++] = z_layer[i-1] + 0.5*dz;
						}
						for(auto i = 0; i <= nz_Top; i++)
						{
							z_slice[j++] = layer_e+i*dz_Bot;
						}
					}
					break;
					case ePS_dz_Proj:
					{
						int nz = quot(s_z, dz_i )+ 1;
						T layer_0 = z_min-0.5*(nz*dz_i-s_z);
						T layer_e = layer_0 + nz*dz_i;
						int nz_Bot = (atoms.z_min<layer_0)?quot(layer_0-atoms.z_min, dz_i) + 1:0;
						int nz_Top = (atoms.z_max>layer_e)?quot(atoms.z_max-layer_e, dz_i) + 1:0;
						z_slice.resize(nz+nz_Bot+nz_Top+1);
		
						z_slice[0] = layer_0-nz_Bot*dz_i;
						for(auto i = 1; i < z_slice.size(); i++)
						{
							z_slice[i] = z_slice[i-1] + dz_i;
						}
					}
					break;
					case ePS_dz_Sub:
					{
						auto get_dz_b = [](const T &z_0, const T &z_e, const T &dz)->T
						{
							T dz_b = fmod(fabs(z_0-z_e), dz);
							dz_b += (dz_b<((dz>2.0)?0.25:0.50)*dz)?dz:0.0;
							return dz_b;
						};
						/*******************************************************************/
						int nz = quot(s_z, dz_i) + 1;
						T dz_Bot = get_dz_b(atoms.z_Int_min, z_min-0.5*(nz*dz_i-s_z), dz_i);
						T dz_Top = get_dz_b(atoms.z_Int_max, z_max+0.5*(nz*dz_i-s_z), dz_i);

						z_slice.resize(quot(atoms.s_z_Int-dz_Bot-dz_Top, dz_i) + 3);
						/*******************************************************************/
						z_slice[0] = atoms.z_Int_min;
						for(auto i = 1; i<z_slice.size(); i++)
						{
							T dz = (i == 1)?dz_Bot:(i == z_slice.size()-1)?dz_Top:dz_i;
							z_slice[i] = z_slice[i-1] + dz;
						}
					}
					break;
				}
				z_slice.shrink_to_fit();
			}

			bool is_periodic_along_z()
			{
				if(z_layer.empty())
				{
					return false;
				}
				else if(z_layer.size()<2)
				{
					return true;
				}

				Vector<T, e_host> dis(z_layer.size()-1);
				for(auto i = 0; i<dis.size(); i++)
				{
					dis[i] = z_layer[i+1] - z_layer[i];
				}
				T min_dis = *std::min_element(dis.begin(), dis.end());
				T max_dis = *std::max_element(dis.begin(), dis.end());

				T dz = 0.01;
				int nhist = static_cast<int>(ceil((max_dis-min_dis)/dz)+1);
				Vector<int, e_host> hist(nhist, 0);
				for(auto i = 0; i<dis.size(); i++)
				{
					int j = static_cast<int>(floor((dis[i]-min_dis)/dz+0.5));
					if(j<hist.size())
					{
						hist[j]++;
					}
				}

				int counter = 0;
				for(auto i = 0; i<hist.size(); i++)
				{
					if(hist[i]>0)
					{
						counter++;
					}
				}

				if(counter>10)
				{
					return false;
				}
				else
				{
					return true;
				}
			}

			inline
			r3d<T> to_r3d(const int &iatom)
			{
				return r3d<T>(x[iatom], y[iatom], z[iatom]);
			}

			inline
			T norm(const int &iatom, const T &x_i, const T &y_i, const T &z_i)
			{
				auto x_d = x[iatom]-x_i; 
				auto y_d = y[iatom]-y_i; 
				auto z_d = z[iatom]-z_i;
				return x_d*x_d + y_d*y_d + z_d*z_d;
			}

			inline
			T distance(const int &iatom, const T &x_i, const T &y_i, const T &z_i)
			{
				return sqrt(norm(iatom, x_i, y_i, z_i));
			}

			T l_x; 			// Box m_size-x
			T l_y; 			// Box m_size-y
			T l_z; 			// Box m_size-z

			Vector<int, e_host> Z;
			Vector<int, e_host> charge;
			Vector<T, e_host> x;
			Vector<T, e_host> y;
			Vector<T, e_host> z;
			Vector<float, e_host> sigma;
			Vector<float, e_host> occ;

			Vector<T, e_host> z_layer;

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

			T s_x; 			// m_size-x
			T s_y; 			// m_size-y
			T s_z; 			// m_size-z

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

		private:
			struct sort_atoms_by_z
			{
				template<class Ttuple1, class Ttuple2>
				__host__ __device__
				bool operator()(const Ttuple1 &t1, const Ttuple2 &t2)
				{
					return thrust::get<4>(t1) < thrust::get<4>(t2);
				}
			};
	};

} // namespace multem

#endif