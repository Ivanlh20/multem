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
#include <multem/ext/add_amorp_lay.h>
#include <multem/ext/amorp_spec.h>
#include <multem/ext/crystal_by_layers.h>
#include <multem/ext/crystal_parameters.h>
#include <multem/ext/feg.h>
#include <multem/ext/gmax.h>
#include <multem/ext/min_spl.h>
#include <multem/ext/misc.h>
/* #include <multem/ext/pr.h> */
/* #include <multem/ext/rdf_3d.h> */
/* #include <multem/ext/spec_rot.h> */
/* #include <multem/ext/vp.h> */

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
  export_add_amorp_lay(m);
  export_amorp_spec(m);
  export_crystal_by_layers(m);
  export_crystal_parameters(m);
  export_feg(m);
  export_gmax(m);
  export_min_spl(m);
  export_misc(m);
  /* export_pr(m); */
  /* export_rdf_3d(m); */
  /* export_spec_rot(m); */
  /* export_vp(m); */
}

