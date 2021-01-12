/*
 *  scanning.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */

#ifndef MULTEM_PYTHON_SCANNING_H
#define MULTEM_PYTHON_SCANNING_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {
  
  /**
   * Define wrapper function for the mt::Scanning class
   */
  template <>
  template <typename T>
  struct Helpers <mt::Scanning<T>> {

    /**
     * Get the state
     */
    static py::tuple getstate(const mt::Scanning<T> &self) {
      return py::make_tuple(
          self.type,
          self.pbc,
          self.spxs,
          self.ns,
          self.x0,
          self.y0,
          self.xe,
          self.ye);
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
  
}}


template <typename T>
void wrap_scanning(py::module_ m)
{
  typedef mt::Scanning<T> Type;

  // Wrap the mt::Input class
  py::class_<Type>(m, "Scanning")
    .def(py::init<>())
    .def_readwrite("type", &Type::type)
    .def_readwrite("pbc", &Type::pbc)
    .def_readwrite("spxs", &Type::spxs)
    .def_readwrite("ns", &Type::ns)
    .def_readwrite("x0", &Type::x0)
    .def_readwrite("y0", &Type::y0)
    .def_readwrite("xe", &Type::xe)
    .def_readwrite("ye", &Type::ye)
    .def(py::pickle(
        &py::detail::Helpers<Type>::getstate,
        &py::detail::Helpers<Type>::setstate))
    ;
}

void export_scanning(py::module_ m) {
  wrap_scanning<double>(m);
}

#endif


