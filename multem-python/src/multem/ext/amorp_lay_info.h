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
#ifndef MULTEM_PYTHON_AMORP_LAY_INFO_H
#define MULTEM_PYTHON_AMORP_LAY_INFO_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <multem/multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  /**
   * Type cast a mt::Amorp_Lay_Info<T> object to a tuple
   */
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

    static handle cast(mt::Amorp_Lay_Info<T> src,
                       return_value_policy policy,
                       handle parent) {
      return py::make_tuple(src.z_0, src.z_e, src.dz).release();
    }
  };

}}  // namespace pybind11::detail

void export_amorp_lay_info(py::module_ m) {}

#endif
