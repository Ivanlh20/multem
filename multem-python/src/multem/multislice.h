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

namespace pybind11 { namespace detail {

  template <typename T>
  void check_input_multislice(mt::Input<T> &input) {
    bool pbc_xy = true;

    // Get atom statistic or delete atoms
    if (input.is_specimen_required()) {
      input.get_atoms().get_statistic();
    } else {
      input.get_atoms().clear();
    }

    // Delete thick if needed
    if (!input.is_specimen_required() || input.is_whole_spec()) {
      input.set_thick(std::vector<double>());
    }

    // Set grid input data
    input.get_grid_2d().set_input_data(
        input.get_grid_2d().nx, 
        input.get_grid_2d().ny, 
        input.get_grid_2d().lx, 
        input.get_grid_2d().ly, 
        input.get_grid_2d().dz, 
        input.get_grid_2d().bwl, 
        pbc_xy);

    // Delete input wave if not user defined
    if (!(input.is_user_define_wave())) {
      input.set_iw_psi(std::vector<std::complex<T>>());
    }

    // Set lens input data
    input.get_cond_lens().set_input_data(
        input.get_E_0(), 
        input.get_grid_2d());
    input.get_obj_lens().set_input_data(
        input.get_E_0(), 
        input.get_grid_2d());

    // Set scanning grid
    if (input.is_scanning()) {
      input.get_scanning().set_grid();
    }

    // Check STEM, EELS and EFTEM parameters
    if (input.is_STEM()) {
      T lambda = input.get_cond_lens().lambda;
      switch (input.get_detector().type) {
      case mt::eDT_Circular:
        {
          auto& detector = input.get_detector();
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
        }
        break;
      };
    } else if (input.is_EELS()) {
      input.get_eels_fr().set_input_data(
          mt::eS_Reciprocal, 
          input.get_E_0(), 
          input.get_eels_fr().E_loss, 
          input.get_eels_fr().m_selection, 
          input.get_eels_fr().collection_angle, 
          input.get_eels_fr().channelling_type, 
          input.get_eels_fr().Z);
    } else if (input.is_EFTEM()) {
      input.get_eels_fr().set_input_data(
          mt::eS_Real, 
          input.get_E_0(), 
          input.get_eels_fr().E_loss, 
          input.get_eels_fr().m_selection, 
          input.get_eels_fr().collection_angle, 
          input.get_eels_fr().channelling_type, 
          input.get_eels_fr().Z);
    }

    // Validate the parameters
    input.validate_parameters();
  }

  object tem_simulation(
      mt::Input<double> &input, 
      const mt::System_Configuration &sys_conf) {
    input.set_system_conf(sys_conf);
    check_input_multislice(input);
    return py::cast(mt::tem_simulation<double>(input));
  }

}}

void export_multislice(py::module_ &m) {
  m.def("tem_simulation", &py::detail::tem_simulation);
}

#endif


