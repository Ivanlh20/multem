#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <multem.h>
#include <multem/system_configuration.h>
#include <multem/input_multislice.h>

namespace py = pybind11;

PYBIND11_MODULE(multem_ext, m)
{
  export_system_configuration(m);
  export_input_multislice(m);
}

