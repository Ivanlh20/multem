/*
 *  enums.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */

#ifndef MULTEM_PYTHON_ENUMS_H
#define MULTEM_PYTHON_ENUMS_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/serialization.h>

namespace py = pybind11;

template <typename Module>
void export_enums(Module m)
{
  py::enum_<mt::eDevice>(m, "eDevice")
    .value("host", mt::eDevice::e_host)
    .value("device", mt::eDevice::e_device)
    .value("host_device", mt::eDevice::e_host_device)
    .export_values()
    ;

  py::enum_<mt::ePrecision>(m, "ePrecision")
    .value("float", mt::ePrecision::eP_float)
    .value("double", mt::ePrecision::eP_double)
    .export_values()
    ;
}

#endif

