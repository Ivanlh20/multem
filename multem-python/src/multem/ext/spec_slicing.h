/*
 *  spec_slicing.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */

#ifndef MULTEM_PYTHON_SPEC_SLICING_H
#define MULTEM_PYTHON_SPEC_SLICING_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  template <typename T>
  void check_input_for_spec_slicing(mt::Input_Multislice<T> &input) {
    bool pbc_xy = true;

    // Set parameters
    input.pn_coh_contrib = false;
    input.pn_single_conf = true;
    input.grid_2d.bwl = false;

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
    input.grid_2d.set_input_data(
        input.grid_2d.nx, 
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
    input.cond_lens.set_input_data(
        input.E_0, 
        input.grid_2d);
    input.obj_lens.set_input_data(
        input.E_0, 
        input.grid_2d);

    // Set scanning grid
    if (input.is_scanning()) {
      input.scanning.set_grid();
    }

    // Validate the parameters
    input.validate_parameters();
  }

  template <typename T>
  object spec_slicing(mt::Input_Multislice<T> &input) {
    check_input_for_spec_slicing(input);
    return py::cast(mt::spec_slicing<T>(input));
  }

}}

void export_spec_slicing(py::module_ &m) {
  m.def("spec_slicing", &py::detail::spec_slicing<float>);
  m.def("spec_slicing", &py::detail::spec_slicing<double>);
}

#endif




