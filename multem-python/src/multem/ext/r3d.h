/*
 *  r3d.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */

#ifndef MULTEM_PYTHON_R3D_H
#define MULTEM_PYTHON_R3D_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {
  
  /**
   * Type cast a mt::r3d<T> object to a tuple
   */
  template <> 
  template <typename T>
  class type_caster<mt::r3d<T>> {
  public:
  
    PYBIND11_TYPE_CASTER(mt::r3d<T>, _("mt::r3d<T>"));

    bool load(handle src, bool convert) {
      if (py::isinstance<py::tuple>(src)) {
        py::tuple t = py::cast<py::tuple>(src);
        if (py::len(t) == 3) {
          value.x = py::cast<double>(t[0]);
          value.y = py::cast<double>(t[1]);
          value.z = py::cast<double>(t[2]);
          return true;
        }
      }
      return false;
    }

    static handle cast(mt::r3d<T> src, return_value_policy policy, handle parent) {
      return py::make_tuple(
        src.x, 
        src.y, 
        src.z).release();
    }
  };
 
}}

void export_r3d(py::module_ m) {
}

#endif





