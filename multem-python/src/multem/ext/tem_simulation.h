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
#ifndef MULTEM_PYTHON_TEM_SIMULATION_H
#define MULTEM_PYTHON_TEM_SIMULATION_H

#include <pybind11/pybind11.h>
#include <multem/multem.h>
#include <multem/ext/serialization.h>
#include <multem/ext/traits.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  template <typename T>
  void check_input_for_tem_simulation(mt::Input_Multislice<T> &input) {
    bool pbc_xy = true;

    // Check input system conf
    MULTEM_ASSERT(input.system_conf.precision = type_precision<T>::value);

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

    // Set scanning grid
    if (input.is_scanning()) {
      input.scanning.set_grid();
    }

    // Check STEM, EELS and EFTEM parameters
    if (input.is_STEM()) {
      auto lambda = input.cond_lens.lambda;
      switch (input.detector.type) {
      case mt::eDT_Circular: {
        auto &detector = input.detector;
        auto inner_ang = detector.inner_ang;
        auto outer_ang = detector.outer_ang;
        std::vector<T> g_inner(inner_ang.size());
        std::vector<T> g_outer(outer_ang.size());
        for (auto i = 0; i < detector.size(); i++) {
          g_inner[i] = sin(inner_ang[i] * mt::c_mrad_2_rad) / lambda;
          g_outer[i] = sin(outer_ang[i] * mt::c_mrad_2_rad) / lambda;
        }
        detector.g_inner = g_inner;
        detector.g_outer = g_outer;
      } break;
      };
    } else if (input.is_EELS()) {
      input.eels_fr.set_input_data(mt::eS_Reciprocal,
                                   input.E_0,
                                   input.eels_fr.E_loss,
                                   input.eels_fr.m_selection,
                                   input.eels_fr.collection_angle,
                                   input.eels_fr.channelling_type,
                                   input.eels_fr.Z);
    } else if (input.is_EFTEM()) {
      input.eels_fr.set_input_data(mt::eS_Real,
                                   input.E_0,
                                   input.eels_fr.E_loss,
                                   input.eels_fr.m_selection,
                                   input.eels_fr.collection_angle,
                                   input.eels_fr.channelling_type,
                                   input.eels_fr.Z);
    }

    // Validate the parameters
    input.validate_parameters();
  }

  template <typename T>
  object tem_simulation(mt::Input_Multislice<T> &input) {
    check_input_for_tem_simulation(input);
    return py::cast(mt::tem_simulation<T>(input));
  }

}}  // namespace pybind11::detail

void export_tem_simulation(py::module_ &m) {
  m.def("tem_simulation", &py::detail::tem_simulation<float>);
  m.def("tem_simulation", &py::detail::tem_simulation<double>);
}

#endif
