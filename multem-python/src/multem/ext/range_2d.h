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
#ifndef MULTEM_PYTHON_RANGE_2D_H
#define MULTEM_PYTHON_RANGE_2D_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  /**
   * Define wrapper function for the mt::Range_2d class
   */
  template <>
  struct Helpers<mt::Range_2d> {
    /**
     * Get the state
     */
    static py::tuple getstate(const mt::Range_2d &self) {
      return py::make_tuple(
        self.ix_0, self.ix_e, self.iy_0, self.iy_e, self.ixy_0, self.ixy_e);
    }

    /**
     * Set the state
     */
    static mt::Range_2d setstate(py::tuple obj) {
      mt::Range_2d self;
      self.ix_0 = obj[0].cast<int>();
      self.ix_e = obj[1].cast<int>();
      self.iy_0 = obj[2].cast<int>();
      self.iy_e = obj[3].cast<int>();
      self.ixy_0 = obj[4].cast<int>();
      self.ixy_e = obj[5].cast<int>();
      return self;
    }
  };

}}  // namespace pybind11::detail

template <typename T>
void wrap_range_2d(py::module_ m) {
  typedef mt::Range_2d Type;

  // Wrap the mt::Input class
  py::class_<Type>(m, "Range_2d")
    .def(py::init<>())
    .def_readwrite("ix_0", &Type::ix_0)
    .def_readwrite("ix_e", &Type::ix_e)
    .def_readwrite("iy_0", &Type::iy_0)
    .def_readwrite("iy_e", &Type::iy_e)
    .def_readwrite("ixy_0", &Type::ixy_0)
    .def_readwrite("ixy_e", &Type::ixy_e)
    .def(py::pickle(&py::detail::Helpers<Type>::getstate,
                    &py::detail::Helpers<Type>::setstate));
}

void export_range_2d(py::module_ m) {
  wrap_range_2d<double>(m);
}

#endif
