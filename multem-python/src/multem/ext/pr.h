/*
 *  pr.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */
#ifndef MULTEM_PYTHON_PR
#define MULTEM_PYTHON_PR

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt {
  
  template <typename T>
  T pr(
    ePotential_Type potential_type,
    int Z,
    int charge,
    T r) {

    auto Pr = mx_create_matrix<rmatrix_r>(r.rows, r.cols, plhs[0]);
    auto dPr = mx_create_matrix<rmatrix_r>(r.rows, r.cols, plhs[1]);

    mt::Atom_Type<T, mt::e_host> atom_type;
    mt::Atomic_Data atomic_data(potential_type);
    atomic_data.To_atom_type_CPU(Z, mt::c_Vrl, mt::c_nR, 0.0, atom_type);

    mt::Atom_Cal<T> atomic_fcns_mt;
    atomic_fcns_mt.Set_Atom_Type(potential_type, charge, &atom_type);
    atomic_fcns_mt.Pr_dPr(r.m_size, r.real, Pr.real, dPr.real);
 }

}

void export_pr(py::module_ m) {
  m.def("pr", &mt::pr<float>);
  m.def("pr", &mt::pr<double>);
}

#endif





