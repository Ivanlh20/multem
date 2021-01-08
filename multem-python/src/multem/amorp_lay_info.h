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
#include <multem.h>
#include <multem/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {
  
  /**
   * Define wrapper function for the mt::Amorp_Lay_Info class
   */
  template <typename T>
  struct Amorp_Lay_Info_Wrapper : public mt::Amorp_Lay_Info<T> {

    /**
     * Get the state
     */
    static py::tuple getstate(const Amorp_Lay_Info_Wrapper &self) {
      return py::make_tuple(
          self.z_0,
          self.z_e,
          self.dz,
          self.region,
          self.type);
    }

    /**
     * Set the state
     */
    static Amorp_Lay_Info_Wrapper setstate(py::tuple obj) {
      Amorp_Lay_Info_Wrapper self;
      self.z_0 = obj[0].cast<T>();
      self.z_e = obj[1].cast<T>();
      self.dz = obj[2].cast<T>();
      self.region = obj[3].cast<int>();
      self.type = obj[4].cast<mt::eAmorp_Lay_Type>();
      return self;
    }

  };
  
}}


template <typename Module, typename T>
void wrap_amorp_lay_info(Module m)
{
  typedef py::detail::Amorp_Lay_Info_Wrapper<T> Type;

  // Wrap the mt::Input class
  py::class_<Type>(m, "Amorp_Lay_Info")
    .def(py::init<>())
    .def_readwrite("z_0", &Type::z_0)
    .def_readwrite("z_e", &Type::z_e)
    .def_readwrite("dz", &Type::dz)
    .def_readwrite("region", &Type::region)
    .def_readwrite("type", &Type::type)
    .def_property_readonly("lz", &Type::lz)
    .def("is_at_top", &Type::is_at_top)
    .def("is_at_bottom", &Type::is_at_bottom)
    .def("is_at_middle", &Type::is_at_middle)
    .def(py::pickle(
        &Type::getstate,
        &Type::setstate))
    ;
}

template <typename Module>
void export_amorp_lay_info(Module m) {
  wrap_amorp_lay_info<Module, float>(m);
  wrap_amorp_lay_info<Module, double>(m);
}

#endif




