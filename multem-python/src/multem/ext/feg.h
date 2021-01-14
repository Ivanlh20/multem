/*
 *  feg.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */
#ifndef MULTEM_PYTHON_FEG
#define MULTEM_PYTHON_FEG

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt {

  template <typename T>
  Atom_Data<T> feg(
      ePotential_Type potential_type,
      int Z,
      int charge,
      T g) {

    mt::Atom_Data<T> atoms;

    mt::Atom_Type<double, mt::e_host> atom_type;
    mt::Atomic_Data atomic_data(potential_type);
    atomic_data.To_atom_type_CPU(Z, mt::c_Vrl, mt::c_nR, 0.0, atom_type);

    mt::Atom_Cal<double> atomic_fcns_mt;
    atomic_fcns_mt.Set_Atom_Type(potential_type, charge, &atom_type);
    atomic_fcns_mt.feg_dfeg(g.m_size, g.real, feg.real, dfeg.real);

    return atoms;
  }

}


void export_feg(py::module_ m) {
  m.def("feg", &mt::feg<double>);
}

#endif


