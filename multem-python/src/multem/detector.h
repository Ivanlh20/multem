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
  struct DetectorWrapper : public mt::DetectorData<T> {
    std::vector<T> inner_ang; 
    std::vector<T> outer_ang; 

    /**
     * Get the state
     */
    static py::tuple getstate(const DetectorWrapper &self) {
      return py::make_tuple(
          self.get_type(),
          self.get_fx(),
          self.get_fR(),
          self.get_grid_1d(),
          self.get_grid_2d(),
          self.get_fn(),
          self.get_g_inner(),
          self.get_g_outer(),
          self.inner_ang,
          self.outer_ang);

    }

    /**
     * Set the state
     */
    static DetectorWrapper setstate(py::tuple obj) {
      DetectorWrapper self;
      self.set_type(obj[0].cast<mt::eDetector_Type>());
      self.set_fx(obj[1].cast<std::vector<std::vector<T>>>());
      self.set_fR(obj[2].cast<std::vector<std::vector<T>>>());
      self.set_grid_1d(obj[3].cast<std::vector<mt::Grid_1d<T>>>());
      self.set_grid_2d(obj[4].cast<std::vector<mt::Grid_2d<T>>>());
      self.set_fn(obj[5].cast<std::vector<std::string>>());
      self.set_g_inner(obj[6].cast<std::vector<T>>());
      self.set_g_outer(obj[7].cast<std::vector<T>>());
      self.set_inner_ang(obj[8].cast<std::vector<T>>());
      self.set_outer_ang(obj[9].cast<std::vector<T>>());
      return self;
    }

  };
  
}}


template <typename Module>
void export_detector(Module m)
{
  typedef double T;

  // Wrap the mt::Input class
  py::class_<py::detail::DetectorWrapper<T>>(m, "Detector")
    .def(py::init<>())
    .def("size", &py::detail::DetectorWrapper<T>::size)
    .def("clear", &py::detail::DetectorWrapper<T>::clear)
    .def("resize", &py::detail::DetectorWrapper<T>::resize)
    .def("is_detector_circular", &py::detail::DetectorWrapper<T>::is_detector_circular)
    .def("is_detector_radial", &py::detail::DetectorWrapper<T>::is_detector_radial)
    .def("is_detector_matrix", &py::detail::DetectorWrapper<T>::is_detector_matrix)
    .def_property(
        "type",
        &py::detail::DetectorWrapper<T>::get_type,
        &py::detail::DetectorWrapper<T>::set_type)
    .def_property(
        "fx",
        &py::detail::DetectorWrapper<T>::get_fx,
        &py::detail::DetectorWrapper<T>::set_fx)
    .def_property(
        "fR",
        &py::detail::DetectorWrapper<T>::get_fR,
        &py::detail::DetectorWrapper<T>::set_fR)
    .def_property(
        "grid_1d",
        &py::detail::DetectorWrapper<T>::get_grid_1d,
        &py::detail::DetectorWrapper<T>::set_grid_1d)
    .def_property(
        "grid_2d",
        &py::detail::DetectorWrapper<T>::get_grid_2d,
        &py::detail::DetectorWrapper<T>::set_grid_2d)
    .def_property(
        "fn",
        &py::detail::DetectorWrapper<T>::get_fn,
        &py::detail::DetectorWrapper<T>::set_fn)
    .def_readwrite("inner_ang", &py::detail::DetectorWrapper<T>::inner_ang)
    .def_readwrite("outer_ang", &py::detail::DetectorWrapper<T>::outer_ang)
    ;
}

#endif

