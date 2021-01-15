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
#ifndef MULTEM_PYTHON_CRYSTAL_PARAMETERS
#define MULTEM_PYTHON_CRYSTAL_PARAMETERS

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt {

/**
 * A class to hold crystal parameters
 */
template <typename T>
class CrystalParameters {
public:
  typedef std::vector<Atom<T>> Layer;

  int na;
  int nb;
  int nc;
  T a;
  T b;
  T c;
  std::vector<Layer> layers;

  CrystalParameters() : na(0), nb(0), nc(0), a(0), b(0), c(0) {}
};

}  // namespace mt

void export_crystal_parameters(py::module_ m) {
  typedef mt::CrystalParameters<double> Type;

  py::class_<Type>(m, "CrystalParameters")
    .def(py::init<>())
    .def_readwrite("na", &Type::na)
    .def_readwrite("nb", &Type::nb)
    .def_readwrite("nc", &Type::nc)
    .def_readwrite("a", &Type::a)
    .def_readwrite("b", &Type::b)
    .def_readwrite("c", &Type::c)
    .def_readwrite("layers", &Type::layers);
}

#endif
