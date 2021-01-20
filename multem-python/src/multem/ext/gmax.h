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
#ifndef MULTEM_PYTHON_GMAX
#define MULTEM_PYTHON_GMAX

#include <pybind11/pybind11.h>
#include <multem/multem.h>

namespace py = pybind11;

namespace mt {

template <typename T>
T gmax(T lx, int nx, T ly, int ny) {
  return 0.5 * std::min(nx / lx, ny / ly);
}

}  // namespace mt

void export_gmax(py::module_ m) {
  m.def("gmax", &mt::gmax<float>);
  m.def("gmax", &mt::gmax<double>);
}

#endif
