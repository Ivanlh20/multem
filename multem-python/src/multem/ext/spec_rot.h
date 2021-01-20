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
#ifndef MULTEM_PYTHON_SPEC_ROT
#define MULTEM_PYTHON_SPEC_ROT

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/ext/atom_data.h>

namespace py = pybind11;

namespace mt {

  template <typename T>
  void spec_rot(Atom_Data<T> &atoms,
                T theta,
                mt::r3d<T> u0,
                mt::eRot_Point_Type &rot_point_type,
                mt::r3d<T> p0) {
    u0.normalized();

    if (rot_point_type == mt::eRPT_geometric_center) {
      p0 = r3d<T>(atoms.x_mean, atoms.y_mean, atoms.z_mean);
    }

    mt::rotate_atoms(atoms, (T)(theta * mt::c_deg_2_rad), u0, p0);
  }
  
  template <typename T>
  Atom_Data<T> spec_rot2(const std::vector<Atom<T>> &spec_atoms,
                 T theta,
                 r3d<T> u0,
                 eRot_Point_Type &rot_point_type,
                 r3d<T> p0) {
    Atom_Data<T> atoms;
    py::detail::Helpers<Atom_Data<T>>::set_spec_atoms_internal(atoms, spec_atoms);
    spec_rot(atoms, theta, u0, rot_point_type, p0);
    return atoms;
  }

}  // namespace mt

void export_spec_rot(py::module_ m) {
  m.def("spec_rot", &mt::spec_rot<float>);
  m.def("spec_rot", &mt::spec_rot<double>);
  m.def("spec_rot", &mt::spec_rot2<float>);
  m.def("spec_rot", &mt::spec_rot2<double>);
}

#endif
