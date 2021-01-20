/*
 * This file is part of multem-python.
 * Copyright 2021 Diamond Light Source
 * Copyright 2021 Rosalind Franklin Institute
 *
 * Author: James Parkhurst
 *
 * multem-python is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * multem-python is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with multem-python. If not, see <http:// www.gnu.org/licenses/>.
 */
#ifndef MULTEM_PYTHON_ADD_AMORP_LAY
#define MULTEM_PYTHON_ADD_AMORP_LAY

#include <pybind11/pybind11.h>
#include <multem/multem.h>
#include <multem/ext/atom_data.h>

namespace py = pybind11;

namespace mt {

  template <typename T>
  Atom_Data<T> add_amorp_lay(Atom_Data<T> atoms,
                             T l_z,
                             T r_min,
                             int Z,
                             T rms_3d,
                             T rho,
                             mt::eAmorp_Lay_Type lay_pos = mt::eALT_Top,
                             int seed = 300183) {
    auto region = atoms.region_max + 1;

    atoms.amorp_lay_info.resize(1);
    if (lay_pos == mt::eALT_Top) {
      atoms.amorp_lay_info[0].z_0 = atoms.z_min - l_z;
      atoms.amorp_lay_info[0].z_e = atoms.z_min;
      atoms.amorp_lay_info[0].dz = 2.0;
      atoms.amorp_lay_info[0].type = mt::eALT_Top;
      atoms.amorp_lay_info[0].region = region;
    } else {
      atoms.amorp_lay_info[0].z_0 = atoms.z_max;
      atoms.amorp_lay_info[0].z_e = atoms.z_max + l_z;
      atoms.amorp_lay_info[0].dz = 2.0;
      atoms.amorp_lay_info[0].type = mt::eALT_Bottom;
      atoms.amorp_lay_info[0].region = region;
    }

    mt::Amorp_Spec<T> spec;
    spec.create(atoms, r_min, Z, rms_3d, rho, seed);

    return atoms;
  }


  template <typename T>
  Atom_Data<T> add_amorp_lay_from_atom_list(
      const std::vector<Atom<T>> &atoms,
      T l_x,
      T l_y,
      T l_z,
      T r_min,
      int Z,
      T rms_3d,
      T rho,
      mt::eAmorp_Lay_Type lay_pos = mt::eALT_Top,
      int seed = 300183) {
    Atom_Data<T> atom_data;
    atom_data.l_x = l_x;
    atom_data.l_y = l_y;
    py::detail::Helpers<Atom_Data<T>>::set_spec_atoms_internal(atom_data, atoms);
    return add_amorp_lay(atom_data, l_z, r_min, Z, rms_3d, rho, lay_pos, seed);
  }

}  // namespace mt

void export_add_amorp_lay(py::module_ m) {
  m.def("add_amorp_lay", &mt::add_amorp_lay<float>);
  m.def("add_amorp_lay", &mt::add_amorp_lay<double>);
  m.def("add_amorp_lay", &mt::add_amorp_lay_from_atom_list<float>);
  m.def("add_amorp_lay", &mt::add_amorp_lay_from_atom_list<double>);
}

#endif
