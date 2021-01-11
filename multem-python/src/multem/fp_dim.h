/*
 *  fp_dim.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */

#ifndef MULTEM_PYTHON_FP_DIM_H
#define MULTEM_PYTHON_FP_DIM_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {
  
  /**
   * Type cast a mt::FP_Dim object to a tuple
   */
  template <> 
  class type_caster<mt::FP_Dim> {
  public:
  
    PYBIND11_TYPE_CASTER(mt::FP_Dim, _("mt::FP_Dim"));

    bool load(handle src, bool convert) {
      if (py::isinstance<py::tuple>(src)) {
        py::tuple t = py::cast<py::tuple>(src);
        if (py::len(t) == 3) {
          value.x = py::cast<bool>(t[0]);
          value.y = py::cast<bool>(t[1]);
          value.z = py::cast<bool>(t[2]);
          return true;
        }
      }
      return false;
    }

    static handle cast(mt::FP_Dim src, return_value_policy policy, handle parent) {
      return py::make_tuple(
        src.x, 
        src.y, 
        src.z).release();
    }
  };
 
}}

template <typename Module>
void export_fp_dim(Module m) {
}

#endif






