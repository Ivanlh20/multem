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

