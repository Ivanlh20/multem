/*
 *  range_2d.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
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
  struct Helpers <mt::Range_2d> {

    /**
     * Get the state
     */
    static py::tuple getstate(const mt::Range_2d &self) {
      return py::make_tuple(
          self.ix_0,
          self.ix_e,
          self.iy_0,
          self.iy_e,
          self.ixy_0,
          self.ixy_e);
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
  

}}


template <typename T>
void wrap_range_2d(py::module_ m)
{
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
    .def(py::pickle(
        &py::detail::Helpers<Type>::getstate,
        &py::detail::Helpers<Type>::setstate))
    ;
}

void export_range_2d(py::module_ m) {
  wrap_range_2d<double>(m);
}

#endif
