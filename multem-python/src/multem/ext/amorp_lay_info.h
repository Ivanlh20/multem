/*
 *  amorp_lay_info.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */

#ifndef MULTEM_PYTHON_AMORP_LAY_INFO_H
#define MULTEM_PYTHON_AMORP_LAY_INFO_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {
  
  /**
   * Type cast a mt::Amorp_Lay_Info<T> object to a tuple
   */
  template <> 
  template <typename T>
  class type_caster<mt::Amorp_Lay_Info<T>> {
  public:
  
    PYBIND11_TYPE_CASTER(mt::Amorp_Lay_Info<T>, _("mt::Amorp_Lay_Info<T>"));

    bool load(handle src, bool convert) {
      if (py::isinstance<py::tuple>(src)) {
        py::tuple t = py::cast<py::tuple>(src);
        if (py::len(t) == 3) {
          value.z_0 = py::cast<double>(t[0]);
          value.z_e = py::cast<double>(t[1]);
          value.dz = py::cast<double>(t[2]);
          return true;
        }
      }
      return false;
    }

    static handle cast(mt::Amorp_Lay_Info<T> src, return_value_policy policy, handle parent) {
      return py::make_tuple(
        src.z_0, 
        src.z_e, 
        src.dz).release();
    }
  };
 
}}


void export_amorp_lay_info(py::module_ m) {
}

#endif




