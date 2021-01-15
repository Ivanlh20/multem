/*
 * This file is part of multem-python.
 * Copyright 2021 Diamond Light Source
 * Copyright 2021 Rosalind Franklin Institute
 *
 * Author: James Parkhurst
 *
 * multem-python is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * multem-python is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with multem-python. If not, see <http:// www.gnu.org/licenses/>.
 */
#include <pybind11/pybind11.h>

#include <multem.h>
#include <multem/ext/add_amorp_lay.h>
#include <multem/ext/amorp_lay_info.h>
#include <multem/ext/amorp_spec.h>
#include <multem/ext/atom_data.h>
#include <multem/ext/crystal_by_layers.h>
#include <multem/ext/crystal_parameters.h>
#include <multem/ext/detector.h>
#include <multem/ext/eels.h>
#include <multem/ext/enums.h>
#include <multem/ext/feg.h>
#include <multem/ext/fp_dim.h>
#include <multem/ext/gmax.h>
#include <multem/ext/grid_2d.h>
#include <multem/ext/incident_wave.h>
#include <multem/ext/input_multislice.h>
#include <multem/ext/lens.h>
#include <multem/ext/microscope_abberations.h>
#include <multem/ext/min_spl.h>
#include <multem/ext/misc.h>
#include <multem/ext/output_multislice.h>
#include <multem/ext/projected_potential.h>
#include <multem/ext/propagate.h>
#include <multem/ext/r3d.h>
#include <multem/ext/range_2d.h>
#include <multem/ext/rdf_3d.h>
#include <multem/ext/scanning.h>
#include <multem/ext/spec_planes.h>
#include <multem/ext/spec_rot.h>
#include <multem/ext/spec_slicing.h>
#include <multem/ext/system_configuration.h>
#include <multem/ext/tem_simulation.h>
#include <multem/ext/transmission_function.h>
#include <multem/ext/wave_function.h>

namespace py = pybind11;

PYBIND11_MODULE(multem_ext, m) {
  export_add_amorp_lay(m);
  export_amorp_lay_info(m);
  export_amorp_spec(m);
  export_atom_data(m);
  export_crystal_by_layers(m);
  export_crystal_parameters(m);
  export_detector(m);
  export_eels(m);
  export_enums(m);
  export_feg(m);
  export_fp_dim(m);
  export_gmax(m);
  export_grid_2d(m);
  export_input_multislice(m);
  export_incident_wave(m);
  export_lens(m);
  export_microscope_abberations(m);
  export_min_spl(m);
  export_misc(m);
  export_output_multislice(m);
  export_projected_potential(m);
  export_propagate(m);
  export_range_2d(m);
  export_r3d(m);
  export_rdf_3d(m);
  export_scanning(m);
  export_spec_planes(m);
  export_spec_rot(m);
  export_spec_slicing(m);
  export_system_configuration(m);
  export_tem_simulation(m);
  export_transmission_function(m);
  export_wave_function(m);
}
