/*
 *  detector.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */

#ifndef MULTEM_PYTHON_DETECTOR_H
#define MULTEM_PYTHON_DETECTOR_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {
  
  /**
   * Define wrapper function for the mt::DetectorData class
   */
  template <typename T>
  struct Helpers <mt::DetectorData<T>> {

    /**
     * Get the state
     */
    static py::tuple getstate(const mt::DetectorData<T> &self) {
      return py::make_tuple(
          self.get_type(),
          self.get_fx(),
          self.get_fR(),
          self.get_grid_1d(),
          self.get_grid_2d(),
          self.get_fn(),
          self.get_g_inner(),
          self.get_g_outer(),
          self.get_inner_ang(),
          self.get_outer_ang());

    }

    /**
     * Set the state
     */
    static void setstate(mt::DetectorData<T> &self, py::tuple obj) {
      self.set_type(obj[0].cast<mt::eDetector_Type>());
      self.set_fx(obj[1].cast<std::vector<std::vector<T>>>());
      self.set_fR(obj[2].cast<std::vector<std::vector<T>>>());
      self.set_grid_1d(obj[3].cast<std::vector<mt::Grid_2d<T>>>());
      self.set_grid_2d(obj[4].cast<std::vector<mt::Grid_2d<T>>>());
      self.set_fn(obj[5].cast<std::vector<std::string>>());
      self.set_g_inner(obj[6].cast<std::vector<T>>());
      self.set_g_outer(obj[7].cast<std::vector<T>>());
      self.set_inner_ang(obj[8].cast<std::vector<T>>());
      self.set_outer_ang(obj[9].cast<std::vector<T>>());
    }

  };
  
}}


template <typename Module, typename T>
void wrap_detector(Module m)
{
  typedef mt::DetectorData<T> Type;

  // Wrap the mt::Input class
  py::class_<Type>(m, "Detector")
    .def("size", &Type::size)
    .def("clear", &Type::clear)
    .def("resize", &Type::resize)
    .def("is_detector_circular", &Type::is_detector_circular)
    .def("is_detector_radial", &Type::is_detector_radial)
    .def("is_detector_matrix", &Type::is_detector_matrix)
    .def_property(
        "type",
        &Type::get_type,
        &Type::set_type)
    .def_property(
        "fx",
        &Type::get_fx,
        &Type::set_fx)
    .def_property(
        "fR",
        &Type::get_fR,
        &Type::set_fR)
    .def_property(
        "grid_1d",
        &Type::get_grid_1d,
        &Type::set_grid_1d)
    .def_property(
        "grid_2d",
        &Type::get_grid_2d,
        &Type::set_grid_2d)
    .def_property(
        "fn",
        &Type::get_fn,
        &Type::set_fn)
    .def_property(
        "inner_ang", 
        &Type::get_inner_ang,
        &Type::set_inner_ang)
    .def_property(
        "outer_ang", 
        &Type::get_outer_ang,
        &Type::set_outer_ang)
    ;
}

template <typename Module>
void export_detector(Module m) {
  wrap_detector<Module, double>(m);
}

#endif

