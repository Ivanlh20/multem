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
#ifndef MULTEM_PYTHON_PROJECTED_POTENTIAL_H
#define MULTEM_PYTHON_PROJECTED_POTENTIAL_H

#include <pybind11/pybind11.h>
#include <multem/multem.h>
#include <multem/ext/serialization.h>
#include <multem/ext/traits.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  template <typename T>
  void check_input_for_projected_potential(mt::Input_Multislice<T> &input) {
    bool pbc_xy = true;

    // Check input system conf
    MULTEM_ASSERT(input.system_conf.precision = type_precision<T>::value);

    // Set some parameters
    input.simulation_type = mt::eTEMST_PPRS;
    input.pn_coh_contrib = true;
    input.pn_single_conf = true;

    // Get atom statistic or delete atoms
    if (input.is_specimen_required()) {
      input.atoms.get_statistic();
    } else {
      input.atoms.clear();
    }

    // Delete thick if needed
    if (!input.is_specimen_required() || input.is_whole_spec()) {
      input.thick.clear();
    }

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
  object projected_potential(mt::Input_Multislice<T> &input) {
    check_input_for_projected_potential(input);
    return py::cast(mt::projected_potential<T>(input));
  }

}}  // namespace pybind11::detail

void export_projected_potential(py::module_ &m) {
  m.def("projected_potential", &py::detail::projected_potential<float>);
  m.def("projected_potential", &py::detail::projected_potential<double>);
}

#endif
