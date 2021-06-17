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
#ifndef MULTEM_PYTHON_SCANNING_H
#define MULTEM_PYTHON_SCANNING_H

#include <pybind11/pybind11.h>
#include <multem/multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  /**
   * Define wrapper function for the mt::Scanning class
   */
  template <typename T>
  struct Helpers<mt::Scanning<T>> {
    /**
     * Get the state
     */
    static py::tuple getstate(const mt::Scanning<T> &self) {
      return py::make_tuple(
        self.type, self.pbc, self.spxs, self.ns, self.x0, self.y0, self.xe, self.ye);
    }

    /**
     * Set the state
     */
    static mt::Scanning<T> setstate(py::tuple obj) {
      mt::Scanning<T> self;
      self.type = obj[0].cast<mt::eScanning_Type>();
      self.pbc = obj[1].cast<bool>();
      self.spxs = obj[2].cast<bool>();
      self.ns = obj[3].cast<int>();
      self.x0 = obj[4].cast<T>();
      self.y0 = obj[5].cast<T>();
      self.xe = obj[6].cast<T>();
      self.ye = obj[7].cast<T>();
      return self;
    }
  };

  py::object Scanning_constructor(const std::string &dtype) {
    MULTEM_ASSERT(dtype == "float" || dtype == "double");
    return (dtype == "float"
      ? py::cast(mt::Scanning<float>())
      : py::cast(mt::Scanning<double>()));
  }

}}  // namespace pybind11::detail

template <typename T>
void wrap_scanning(py::module_ m, const char *name) {
  typedef mt::Scanning<T> Type;

  // Wrap the mt::Input class
  py::class_<Type>(m, name)
    .def(py::init<>())
    .def_readwrite("type", &Type::type)
    .def_readwrite("pbc", &Type::pbc)
    .def_readwrite("spxs", &Type::spxs)
    .def_readwrite("ns", &Type::ns)
    .def_readwrite("x0", &Type::x0)
    .def_readwrite("y0", &Type::y0)
    .def_readwrite("xe", &Type::xe)
    .def_readwrite("ye", &Type::ye)
    .def(py::pickle(&py::detail::Helpers<Type>::getstate,
                    &py::detail::Helpers<Type>::setstate));
}

void export_scanning(py::module_ m) {

  m.def(
      "Scanning", 
      &py::detail::Scanning_constructor, 
      py::arg("dtype") = "double");

  wrap_scanning<float>(m, "Scanning_f");
  wrap_scanning<double>(m, "Scanning_d");
}

#endif
