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
   * Define wrapper function for the mt::ScanningData class
   */
  template <typename T>
  struct ScanningWrapper : public mt::ScanningData<T> {

    /**
     * Get the state
     */
    static py::tuple getstate(const ScanningWrapper &self) {
      return py::make_tuple(
          self.get_type(),
          self.get_pbc(),
          self.get_spxs(),
          self.get_ns(),
          self.get_x0(),
          self.get_y0(),
          self.get_xe(),
          self.get_ye());
    }

    /**
     * Set the state
     */
    static ScanningWrapper setstate(py::tuple obj) {
      ScanningWrapper self;
      self.set_type(obj[0].cast<mt::eScanning_Type>());
      self.set_pbc(obj[1].cast<bool>());
      self.set_spxs(obj[2].cast<bool>());
      self.set_ns(obj[3].cast<int>());
      self.set_x0(obj[4].cast<T>());
      self.set_y0(obj[5].cast<T>());
      self.set_xe(obj[6].cast<T>());
      self.set_ye(obj[7].cast<T>());
      return self;
    }

  };
  
}}


template <typename Module, typename T>
void wrap_scanning(Module m)
{
  typedef py::detail::ScanningWrapper<T> Type;

  // Wrap the mt::Input class
  py::class_<Type>(m, "Scanning")
    .def(py::init<>())
    .def_property(
        "type",
        &Type::get_type,
        &Type::set_type)
    .def_property(
        "pbc",
        &Type::get_pbc,
        &Type::set_pbc)
    .def_property(
        "spxs",
        &Type::get_spxs,
        &Type::set_spxs)
    .def_property(
        "ns",
        &Type::get_ns,
        &Type::set_ns)
    .def_property(
        "x0",
        &Type::get_x0,
        &Type::set_x0)
    .def_property(
        "y0",
        &Type::get_y0,
        &Type::set_y0)
    .def_property(
        "xe",
        &Type::get_xe,
        &Type::set_xe)
    .def_property(
        "ye",
        &Type::get_ye,
        &Type::set_ye)
    .def(py::pickle(
        &Type::getstate,
        &Type::setstate))
    ;
}

template <typename Module>
void export_scanning(Module m) {
  wrap_scanning<Module, double>(m);
}

#endif


