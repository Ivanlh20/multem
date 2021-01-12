/*
 *  multislice.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */

#ifndef MULTEM_PYTHON_MULTISLICE_H
#define MULTEM_PYTHON_MULTISLICE_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  object tem_simulation(
      mt::Input<double> &input, 
      const mt::SystemConfiguration &sys_conf) {
    return py::cast(mt::tem_simulation<double>(input));
  }

}}

void export_multislice(py::module_ &m) {
  m.def("tem_simulation", &py::detail::tem_simulation);
}

#endif


