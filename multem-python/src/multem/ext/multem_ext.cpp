#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <multem.h>
#include <multem/ext/enums.h>
#include <multem/ext/system_configuration.h>
#include <multem/ext/grid_2d.h>
#include <multem/ext/range_2d.h>
#include <multem/ext/r3d.h>
#include <multem/ext/fp_dim.h>
#include <multem/ext/lens.h>
#include <multem/ext/eels.h>
#include <multem/ext/amorp_lay_info.h>
#include <multem/ext/atom_data.h>
#include <multem/ext/scanning.h>
#include <multem/ext/detector.h>
#include <multem/ext/input_multislice.h>
#include <multem/ext/output_multislice.h>
#include <multem/ext/multislice.h>

namespace py = pybind11;

PYBIND11_MODULE(multem_ext, m)
{
  export_enums(m);
  export_system_configuration(m);
  export_grid_2d(m);
  export_range_2d(m);
  export_r3d(m);
  export_fp_dim(m);
  export_lens(m);
  export_eels(m);
  export_amorp_lay_info(m);
  export_atom_data(m);
  export_scanning(m);
  export_detector(m);
  export_input_multislice(m);
  export_output_multislice(m);
  export_multislice(m);
}

