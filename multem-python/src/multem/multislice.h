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

template <typename Module, typename T, mt::eDevice dev>
void wrap_multislice(Module m)
{
  typedef mt::MultisliceData<T,dev> Type;

  // Wrap the mt::Input class
  py::class_<Type>(m, "Multislice")
    .def(py::init<>())
    .def("set_input_data", &Type::set_input_data)
    .def("__call__", &Type::operator())
    ;
}

template <typename Module>
void export_multislice(Module m) {
  wrap_multislice<Module, float, mt::e_host>(m);
  wrap_multislice<Module, float, mt::e_device>(m);
  wrap_multislice<Module, double, mt::e_host>(m);
  wrap_multislice<Module, double, mt::e_device>(m);
}

#endif


