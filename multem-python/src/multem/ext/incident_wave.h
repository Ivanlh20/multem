/*
 *  incident_wave.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */

#ifndef MULTEM_PYTHON_INCIDENT_WAVE_H
#define MULTEM_PYTHON_INCIDENT_WAVE_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  template <typename T>
  void check_input_for_incident_wave(mt::Input_Multislice<T> &input) {
    bool pbc_xy = true;
      
    // Set the simulation type to IWRS and set some parameters
    input.simulation_type = mt::eTEMST_IWRS; 
    input.atoms.l_z = 0;
    input.grid_2d.dz = 0.25;
    input.grid_2d.bwl = false;

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

    // Validate the parameters
    input.validate_parameters();
  }

  template <typename T>
  object incident_wave(
      mt::Input_Multislice<T> &input, 
      const mt::System_Configuration &sys_conf) {
    input.system_conf = sys_conf;
    check_input_for_incident_wave(input);
    return py::cast(mt::incident_wave<T>(input));
  }

}}

void export_incident_wave(py::module_ &m) {
  m.def("incident_wave", &py::detail::incident_wave<float>);
  m.def("incident_wave", &py::detail::incident_wave<double>);
}

#endif



