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
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  /**
   * Define wrapper function for the mt::Detector class
   */
  template <typename T>
  struct Helpers<mt::Detector<T, mt::e_host>> {
    /**
     * Get the state
     */
    static py::tuple getstate(const mt::Detector<T, mt::e_host> &self) {
      return py::make_tuple(self.type,
                            self.fx,
                            self.fR,
                            self.grid_1d,
                            self.grid_2d,
                            self.fn,
                            self.g_inner,
                            self.g_outer,
                            self.inner_ang,
                            self.outer_ang);
    }

    /**
     * Set the state
     */
    static mt::Detector<T, mt::e_host> setstate(py::tuple obj) {
      mt::Detector<T, mt::e_host> self;
      self.type = obj[0].cast<mt::eDetector_Type>();
      self.fx = obj[1].cast<std::vector<std::vector<T>>>();
      self.fR = obj[2].cast<std::vector<std::vector<T>>>();
      self.grid_1d = obj[3].cast<std::vector<mt::Grid_2d<T>>>();
      self.grid_2d = obj[4].cast<std::vector<mt::Grid_2d<T>>>();
      self.fn = obj[5].cast<std::vector<std::string>>();
      self.g_inner = obj[6].cast<std::vector<T>>();
      self.g_outer = obj[7].cast<std::vector<T>>();
      self.inner_ang = obj[8].cast<std::vector<T>>();
      self.outer_ang = obj[9].cast<std::vector<T>>();
      return self;
    }
  };

}}  // namespace pybind11::detail

template <typename T>
void wrap_detector(py::module_ m) {
  typedef mt::Detector<T, mt::e_host> Type;

  // Wrap the mt::Input class
  py::class_<Type>(m, "Detector")
    .def(py::init<>())
    .def("size", &Type::size)
    .def("clear", &Type::clear)
    .def("resize", &Type::resize)
    .def("is_detector_circular", &Type::is_detector_circular)
    .def("is_detector_radial", &Type::is_detector_radial)
    .def("is_detector_matrix", &Type::is_detector_matrix)
    .def_readwrite("type", &Type::type)
    .def_readwrite("fx", &Type::fx)
    .def_readwrite("fR", &Type::fR)
    .def_readwrite("grid_1d", &Type::grid_1d)
    .def_readwrite("grid_2d", &Type::grid_2d)
    .def_readwrite("fn", &Type::fn)
    .def_readwrite("inner_ang", &Type::inner_ang)
    .def_readwrite("outer_ang", &Type::outer_ang)
    .def(py::pickle(&py::detail::Helpers<Type>::getstate,
                    &py::detail::Helpers<Type>::setstate));
}

void export_detector(py::module_ m) {
  wrap_detector<double>(m);
}

#endif
