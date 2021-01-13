/*
 *  system_configuration.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */

#ifndef MULTEM_PYTHON_SYSTEM_CONFIGURATION_H
#define MULTEM_PYTHON_SYSTEM_CONFIGURATION_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  /**
   * Define helper function for the multem::SystemConfiguration class
   */
  template <>
  struct Helpers <mt::System_Configuration> {
   
    /**
     * Get the state
     */
    static py::tuple getstate(const mt::System_Configuration &self) {
      return py::make_tuple(
        self.device,
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

}}


void export_system_configuration(py::module_ m)
{
  typedef mt::System_Configuration Type;
  py::class_<Type>(m, "System_Configuration")
    .def(py::init<>())
    .def_readwrite("precision", &Type::precision)
    .def_readwrite("device", &Type::device)
    .def_readwrite("cpu_ncores", &Type::cpu_ncores)
    .def_readwrite("cpu_nthread", &Type::cpu_nthread)
    .def_readwrite("gpu_device", &Type::gpu_device)
    .def_readwrite("gpu_nstream", &Type::gpu_nstream)
    .def(py::pickle(
        &py::detail::Helpers<Type>::getstate,
        &py::detail::Helpers<Type>::setstate))
    ;
}

#endif
