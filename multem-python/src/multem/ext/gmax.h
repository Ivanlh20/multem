/*
 *  gmax.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */
#ifndef MULTEM_PYTHON_GMAX
#define MULTEM_PYTHON_GMAX

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt {

  template <typename T>
  T gmax(T lx, int nx, T ly, int ny) {
    return 0.5 * std::min(nx / lx, ny / ly);
  }

}


void export_gmax(py::module_ m) {
  m.def("gmax", &mt::gmax<float>);
  m.def("gmax", &mt::gmax<double>);
}

#endif


