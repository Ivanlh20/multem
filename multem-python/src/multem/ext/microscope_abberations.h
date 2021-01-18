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
#ifndef MULTEM_PYTHON_MICROSCOPE_ABBERATIONS_H
#define MULTEM_PYTHON_MICROSCOPE_ABBERATIONS_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/ext/serialization.h>
#include <multem/ext/traits.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  template <typename T>
  void check_input_for_microscope_abberations(mt::Input_Multislice<T> &input) {
    bool pbc_xy = true;

    // Check input system conf
    MULTEM_ASSERT(input.system_conf.precision = type_precision<T>::value);

    // Set the simulation type to IWRS and set some parameters
    input.simulation_type = mt::eTEMST_HRTEM;
    input.atoms.l_z = 0;
    input.grid_2d.dz = 0.25;
    input.grid_2d.bwl = false;
    input.obj_lens.zero_defocus_type = mt::eZDT_Last;
    input.obj_lens.zero_defocus_plane = 0.0;

    // Set grid input data
    input.grid_2d.set_input_data(input.grid_2d.nx,
                                 input.grid_2d.ny,
                                 input.grid_2d.lx,
                                 input.grid_2d.ly,
                                 input.grid_2d.dz,
                                 input.grid_2d.bwl,
                                 pbc_xy);

    // Delete input wave if not user defined
    if (!(input.is_user_define_wave())) {
      input.iw_psi.clear();
    }

    // Set lens input data
    input.cond_lens.set_input_data(input.E_0, input.grid_2d);
    input.obj_lens.set_input_data(input.E_0, input.grid_2d);

    // Validate the parameters
    input.validate_parameters();
  }

  template <typename T>
  object microscope_abberations(mt::Input_Multislice<T> &input) {
    check_input_for_microscope_abberations(input);
    return py::cast(mt::microscope_abberations<T>(input));
  }

}}  // namespace pybind11::detail

void export_microscope_abberations(py::module_ &m) {
  m.def("microscope_abberations", &py::detail::microscope_abberations<float>);
  m.def("microscope_abberations", &py::detail::microscope_abberations<double>);
}

#endif
