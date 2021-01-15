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
#ifndef MULTEM_PYTHON_AMORP_SPEC
#define MULTEM_PYTHON_AMORP_SPEC

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt {

template <typename T>
Atom_Data<T> amorp_spec(T l_z, T r_min, int Z, T rms_3d, T rho, int seed = 300183) {
  mt::Atom_Data<T> atoms;

  atoms.amorp_lay_info.resize(1);
  atoms.amorp_lay_info[0].z_0 = 0;
  atoms.amorp_lay_info[0].z_e = l_z;
  atoms.amorp_lay_info[0].dz = 2.0;
  atoms.amorp_lay_info[0].type = mt::eALT_Bottom;
  atoms.amorp_lay_info[0].region = 0;

  mt::Amorp_Spec<T> spec;
  spec.create(atoms, r_min, Z, rms_3d, rho, seed);

  return atoms;
}

}  // namespace mt

void export_amorp_spec(py::module_ m) {
  m.def("amorp_spec", &mt::amorp_spec<float>);
  m.def("amorp_spec", &mt::amorp_spec<double>);
}

#endif
