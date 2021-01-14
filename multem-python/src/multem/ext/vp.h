/*
 *  vp.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */
#ifndef MULTEM_PYTHON_AMORP_SPEC
#define MULTEM_PYTHON_AMORP_SPEC

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt {

  template <typename T>
  Atom_Data<T> vp(
      Potential_Type potential_type,
      int Z,
      int charge,
      T R) {

    mt::Atom_Type<double, mt::e_host> atom_type;
    mt::Atomic_Data atomic_data(potential_type);
    atomic_data.To_atom_type_CPU(Z, mt::c_Vrl, mt::c_nR, 0.0, atom_type);

    mt::Atom_Cal<double> atomic_fcns_mt;
    atomic_fcns_mt.Set_Atom_Type(potential_type, charge, &atom_type);
    atomic_fcns_mt.VR_dVR(R.m_size, R.real, VR.real, dVR.real);

    return atoms;
  }
  
  template <typename T>
  Atom_Data<T> vr(
      Potential_Type potential_type,
      int Z,
      int charge,
      T r) {

    mt::Atom_Type<double, mt::e_host> atom_type;
    mt::Atomic_Data atomic_data(potential_type);
    atomic_data.To_atom_type_CPU(Z, mt::c_Vrl, mt::c_nR, 0.0, atom_type);

    mt::Atom_Cal<double> atomic_fcns_mt;
    atomic_fcns_mt.Set_Atom_Type(potential_type, charge, &atom_type);
    atomic_fcns_mt.Vr_dVr(r.m_size, r.real, Vr.real, dVr.real);

    return atoms;
  }
  
  template <typename T>
  Atom_Data<T> vr(
      Potential_Type potential_type,
      int Z,
      int charge,
      double z0,
      double ze,
      T R) {

    auto VR = mx_create_matrix<rmatrix_r>(R.rows, R.cols, plhs[0]);
    auto dVR = mx_create_matrix<rmatrix_r>(R.rows, R.cols, plhs[1]);

    mt::Atom_Type<double, mt::e_host> atom_type;
    mt::Atomic_Data atomic_data(potential_type);
    atomic_data.To_atom_type_CPU(Z, mt::c_Vrl, mt::c_nR, 0.0, atom_type);

    mt::Atom_Cal<double> atomic_fcns_mt;
    atomic_fcns_mt.Set_Atom_Type(potential_type, charge, &atom_type);
    atomic_fcns_mt.Vz_dVz(z0, ze, R.m_size, R.real, VR.real, dVR.real);

    return atoms;
  }

}


void export_vp(py::module_ m) {
  m.def("vp", &mt::vp<float>);
  m.def("vp", &mt::vp<double>);

  m.def("vr", &mt::vp<float>);
  m.def("vr", &mt::vp<double>);
  
  m.def("vz", &mt::vz<float>);
  m.def("vz", &mt::vz<double>);
}

#endif


