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
#ifndef MULTEM_PYTHON_SYSTEM_CONFIGURATION_H
#define MULTEM_PYTHON_SYSTEM_CONFIGURATION_H

#include <pybind11/pybind11.h>
#include <multem/multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  /**
   * Define helper function for the multem::SystemConfiguration class
   */
  template <>
  struct Helpers<mt::System_Configuration> {
    /**
     * Get the state
     */
    static py::tuple getstate(const mt::System_Configuration &self) {
      return py::make_tuple(self.device,
                            self.precision,
                            self.cpu_ncores,
                            self.cpu_nthread,
                            self.gpu_device,
                            self.gpu_nstream);
    }

    /**
     * Set the state
     */
    static mt::System_Configuration setstate(py::tuple obj) {
      mt::System_Configuration self;
      self.device = obj[0].cast<mt::eDevice>();
      self.precision = obj[1].cast<mt::ePrecision>();
      self.cpu_ncores = obj[2].cast<std::size_t>();
      self.cpu_nthread = obj[3].cast<std::size_t>();
      self.gpu_device = obj[4].cast<std::size_t>();
      self.gpu_nstream = obj[5].cast<std::size_t>();
      return self;
    }
  };

}}  // namespace pybind11::detail

void export_system_configuration(py::module_ m) {
  typedef mt::System_Configuration Type;
  py::class_<Type>(m, "System_Configuration")
    .def(py::init<>())
    .def_readwrite("precision", &Type::precision)
    .def_readwrite("device", &Type::device)
    .def_readwrite("cpu_ncores", &Type::cpu_ncores)
    .def_readwrite("cpu_nthread", &Type::cpu_nthread)
    .def_readwrite("gpu_device", &Type::gpu_device)
    .def_readwrite("gpu_nstream", &Type::gpu_nstream)
    .def(py::pickle(&py::detail::Helpers<Type>::getstate,
                    &py::detail::Helpers<Type>::setstate));
}

#endif
