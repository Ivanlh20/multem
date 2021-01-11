#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <multem.h>
#include <multem/enums.h>
#include <multem/system_configuration.h>
#include <multem/grid_2d.h>
#include <multem/amorp_lay_info.h>
#include <multem/atom_data.h>
#include <multem/scanning.h>
#include <multem/detector.h>
#include <multem/input_multislice.h>
#include <multem/output_multislice.h>
#include <multem/multislice.h>

namespace py = pybind11;

PYBIND11_MODULE(multem_ext, m)
{
  export_enums(m);
  export_system_configuration(m);
  export_grid_2d(m);
  export_amorp_lay_info(m);
  export_atom_data(m);
  export_scanning(m);
  export_detector(m);
  export_input_multislice(m);
  export_output_multislice(m);
  export_multislice(m);
}

