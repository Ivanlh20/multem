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

#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <vector>
#include <cstdlib>

#include "types.hpp"
#include "atom_data.hpp"

namespace multem
{
	template<class T>
	class Crystal{
		public:
			Crystal():na(0), nb(0), nc(0), a(0), b(0), c(0) {};

			void Create3DCrystal(const int &na_i, const int &nb_i, const int &nc_i, T a_i, T b_i, T c_i, Vector<Atom_Data<T>, e_Host> &uLayer_i, Atom_Data<T> &Atoms_o)
			{
				na = na_i;
				nb = nb_i;
				nc = nc_i; 

				a = a_i;
				b = b_i;
				c = c_i;

				uLayer.resize(uLayer_i.size());
				Layers.resize(uLayer_i.size()); 

				for(auto i = 0; i < uLayer.size(); i++)
				{
					uLayer[i].set_Atoms(uLayer_i[i]);
					Layers[i].resize(uLayer[i].size()*(na+1)*(nb+1));
				}

				auto nAtomsLayers = StackLayers(na, nb, a, b, c, uLayer, Layers);
				Atoms_o.resize(nc*nAtomsLayers + Layers[0].size());

				Atoms_o.l_x = static_cast<T>(na)*a;
				Atoms_o.l_y = static_cast<T>(nb)*b;

				std::size_t l = 0;
				for(auto k = 0; k < nc; k++)
				{
					for(auto i = 0; i < Layers.size(); i++)
					{
						for(auto j = 0; j < Layers[i].size(); j++)
						{
							Atoms_o.Z[l] = Layers[i].Z[j];
							Atoms_o.x[l] = Layers[i].x[j];
							Atoms_o.y[l] = Layers[i].y[j];
							Atoms_o.z[l] = Layers[i].z[j] + c*static_cast<T>(k);
							Atoms_o.sigma[l] = Layers[i].sigma[j];
							Atoms_o.occ[l] = Layers[i].occ[j];
							l++;
						}
					}
				}

				// Last layer
				for(auto j = 0; j < Layers[0].size(); j++)
				{
					Atoms_o.Z[l] = Layers[0].Z[j];
					Atoms_o.x[l] = Layers[0].x[j];
					Atoms_o.y[l] = Layers[0].y[j];
					Atoms_o.z[l] = Layers[0].z[j] + c*static_cast<T>(nc);
					Atoms_o.sigma[l] = Layers[0].sigma[j];
					Atoms_o.occ[l] = Layers[0].occ[j];
					l++;
				}
				Atoms_o.get_Statistic();
			}

		private:
			void uLayer_2_Layer(const int &na, const int &nb, T a, T b, T c, Atom_Data<T> &uLayer, Atom_Data<T> &Layer)
			{
				T x, y;

				T xmin = 0.0 - Epsilon<T>::rel; 
				T xmax = static_cast<T>(na)*a + Epsilon<T>::rel;

				T ymin = 0.0 - Epsilon<T>::rel; 
				T ymax = static_cast<T>(nb)*b + Epsilon<T>::rel;

				std::size_t l = 0;
				for(auto j = 0; j <= nb; j++)
				{
					for(auto i = 0; i <= na; i++)
					{
						for(auto k = 0; k < uLayer.size(); k++)
						{
							x = (static_cast<T>(i) + uLayer.x[k])*a;
							y = (static_cast<T>(j) + uLayer.y[k])*b; 			
							if(Check_Bound(x, xmin, xmax, y, ymin, ymax))
							{
								Layer.Z[l] = uLayer.Z[k];
								Layer.x[l] = x;
								Layer.y[l] = y;
								Layer.z[l] = c*uLayer.z[k];
								Layer.sigma[l] = uLayer.sigma[k];
								Layer.occ[l] = uLayer.occ[k];
								l++;
							}
						}
					}
				}
				Layer.resize(l);
			}

			std::size_t StackLayers(const int &na, const int &nb, T a, T b, T c, Vector<Atom_Data<T>, e_Host> &uLayer, Vector<Atom_Data<T>, e_Host> &Layers)
			{
				std::size_t nAtoms_Layers = 0;
				for(auto i = 0; i < uLayer.size(); i++)
				{
					uLayer_2_Layer(na, nb, a, b, c, uLayer[i], Layers[i]);
					nAtoms_Layers += Layers[i].size();
				}
				return nAtoms_Layers;
			}

			int na;
			int nb;
			int nc;

			T a;
			T b;
			T c;

			Vector<Atom_Data<T>, e_Host> uLayer;
			Vector<Atom_Data<T>, e_Host> Layers;
	};

} // namespace multem

#endif
