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

#ifndef XTL_BUILD_H
#define XTL_BUILD_H

#ifdef _MSC_VER
#pragma once
#endif// _MSC_VER

#include <vector>
#include <cstdlib>

#include "types.cuh"
#include "atomic_data_mt.hpp"

namespace mt
{
	template <class T>
	class Crystal_Spec{
		public:
			Crystal_Spec(): na(0), nb(0), nc(0), a(0), b(0), c(0){};

			void operator()(const int &na_i, const int &nb_i, const int &nc_i, T a_i, T b_i, T c_i, Vector<Atom_Data<T>, e_host> &ulay_i, Atom_Data<T> &Atoms_o)
			{
				na = na_i;
				nb = nb_i;
				nc = nc_i; 

				a = a_i;
				b = b_i;
				c = c_i;

				ulay.resize(ulay_i.size());
				lays.resize(ulay_i.size()); 

				for(auto i = 0; i < ulay.size(); i++)
				{
					ulay[i].set_atoms(ulay_i[i]);
					lays[i].resize(ulay[i].size()*(na+1)*(nb+1));
				}

				auto nAtomslays = stack_layers(na, nb, a, b, c, ulay, lays);
				Atoms_o.resize(nc*nAtomslays + lays[0].size());

				Atoms_o.l_x = static_cast<T>(na)*a;
				Atoms_o.l_y = static_cast<T>(nb)*b;

				std::size_t l = 0;
				for(auto k = 0; k < nc; k++)
				{
					for(auto i = 0; i < lays.size(); i++)
					{
						for(auto j = 0; j < lays[i].size(); j++)
						{
							Atoms_o.Z[l] = lays[i].Z[j];
							Atoms_o.x[l] = lays[i].x[j];
							Atoms_o.y[l] = lays[i].y[j];
							Atoms_o.z[l] = lays[i].z[j] + c*static_cast<T>(k);
							Atoms_o.sigma[l] = lays[i].sigma[j];
							Atoms_o.occ[l] = lays[i].occ[j];
							Atoms_o.region[l] = lays[i].region[j];
							Atoms_o.charge[l] = lays[i].charge[j];
							l++;
						}
					}
				}

				// Last layer
				for(auto j = 0; j < lays[0].size(); j++)
				{
					Atoms_o.Z[l] = lays[0].Z[j];
					Atoms_o.charge[l] = lays[0].charge[j];
					Atoms_o.x[l] = lays[0].x[j];
					Atoms_o.y[l] = lays[0].y[j];
					Atoms_o.z[l] = lays[0].z[j] + c*static_cast<T>(nc);
					Atoms_o.sigma[l] = lays[0].sigma[j];
					Atoms_o.occ[l] = lays[0].occ[j];
					Atoms_o.region[l] = lays[0].region[j];
					Atoms_o.charge[l] = lays[0].charge[j];
					l++;
				}
				Atoms_o.get_statistic();
			}

		private:
			void ulayer_2_layer(const int &na, const int &nb, T a, T b, T c, Atom_Data<T> &ulay, Atom_Data<T> &lay)
			{
				T x, y;

				T xmin = 0.0 - Epsilon<T>::rel; 
				T xmax = na*a + Epsilon<T>::rel;

				T ymin = 0.0 - Epsilon<T>::rel; 
				T ymax = nb*b + Epsilon<T>::rel;

				std::size_t l = 0;
				for(auto j = 0; j <= nb; j++)
				{
					for(auto i = 0; i <= na; i++)
					{
						for(auto k = 0; k < ulay.size(); k++)
						{
							x = (i + ulay.x[k])*a;
							y = (j + ulay.y[k])*b; 			
							if(Check_Bound(x, xmin, xmax, y, ymin, ymax))
							{
								lay.Z[l] = ulay.Z[k];
								lay.x[l] = x;
								lay.y[l] = y;
								lay.z[l] = c*ulay.z[k];
								lay.sigma[l] = ulay.sigma[k];
								lay.occ[l] = ulay.occ[k];
								lay.region[l] = ulay.region[k];
								lay.charge[l] = ulay.charge[k];
								l++;
							}
						}
					}
				}
				lay.resize(l);
			}

			std::size_t stack_layers(const int &na, const int &nb, T a, T b, T c, Vector<Atom_Data<T>, e_host> &ulay, Vector<Atom_Data<T>, e_host> &lays)
			{
				std::size_t nAtoms_lays = 0;
				for(auto i = 0; i < ulay.size(); i++)
				{
					ulayer_2_layer(na, nb, a, b, c, ulay[i], lays[i]);
					nAtoms_lays += lays[i].size();
				}
				return nAtoms_lays;
			}

			int na;
			int nb;
			int nc;

			T a;
			T b;
			T c;

			Vector<Atom_Data<T>, e_host> ulay;
			Vector<Atom_Data<T>, e_host> lays;
	};

} // namespace mt

#endif
