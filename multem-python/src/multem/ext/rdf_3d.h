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
#ifndef MULTEM_PYTHON_PR
#define MULTEM_PYTHON_PR

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt {

template <typename T>
std::tuple<std::vector<T>, std::vector<T>> get_rdf_3d(const mt::Atom_Data<T> &atoms,
                                                      T r_max,
                                                      std::size_t nr) {
  std::vector<T> r(nr);
  std::vector<T> rdf(nr);

  mt::rdf_3d(atoms, r_max, nr, r, rdf);

  return std::make_tuple(r, rdf);
}

}  // namespace mt

void export_rdf_3d(py::module_ m) {
  m.def("rdf_3d", &mt::get_rdf_3d<float>);
  m.def("rdf_3d", &mt::get_rdf_3d<double>);
}

#endif
