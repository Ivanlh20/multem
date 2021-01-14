/*
 *  fxg.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */
#ifndef MULTEM_FXG
#define MULTEM_FXG

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt { 

  template <typename T>
  py::object fxg(
      mt::ePotential_Type potential_type,
      int Z, 
      int charge,
      T g) {

    mt::Atom_Type<T, mt::e_host> atom_type;
    mt::Atomic_Data atomic_data(potential_type);
    atomic_data.To_atom_type_CPU(Z, mt::c_Vrl, mt::c_nR, 0.0, atom_type);

    mt::Atom_Cal<T> atomic_fcns_mt;
    atomic_fcns_mt.Set_Atom_Type(potential_type, charge, &atom_type);
    atomic_fcns_mt.fxg_dfxg(g.m_size, g.real, fxg.real, dfxg.real);
  }

}


void export_fxg(py::module_ m) {
  m.def("fxg", &mt::fxg<T>);
  m.def("fxg", &mt::fxg<double>);
}

#endif





