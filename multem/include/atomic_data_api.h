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

#ifndef ATOMIC_DATA_API_H
#define ATOMIC_DATA_API_H

#include <string>
#include <multem/config.h>
#include "math.cuh"
#include "safe_types.cuh"

namespace mt
{
	class Atomic_Data{
		public:
			DLL_PUBLIC Atomic_Data(ePotential_Type PotPar_i = ePT_none);
			DLL_PUBLIC void Load_Data(ePotential_Type PotPar_i);

			DLL_PUBLIC std::string Z_name(const int &Z) const;
      DLL_PUBLIC std::vector<float> extract_major_edges(const int &Z);
      DLL_PUBLIC std::vector<float> extract_minor_edges(const int &Z);
			DLL_PUBLIC bool Z_bound(const int &Z) const;
			DLL_PUBLIC int mass_number(const int &Z) const;
			DLL_PUBLIC double atomic_mass(const int &Z) const;
			DLL_PUBLIC double nuclear_radius(const int &Z) const;
			DLL_PUBLIC double atomic_radius(const int &Z) const;

			template <class TAtom_Type>
			DLL_PUBLIC void To_atom_type_CPU(int Z_i, double Vrl_i, int nR_i, double R_min_i, TAtom_Type &atom_type);

		private:
      class Data;
      std::unique_ptr<Data> pimpl;
	};

} // namespace mt
#endif
