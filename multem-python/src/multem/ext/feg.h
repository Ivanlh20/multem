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
  std::tuple<std::vector<T>, std::vector<T>> feg(
      ePotential_Type potential_type,
      int Z,
      int charge,
      std::vector<T> g) {

    std::vector<T> feg(g.size());
    std::vector<T> dfeg(g.size());

    mt::Atom_Type<T, mt::e_host> atom_type;
    mt::Atomic_Data atomic_data(potential_type);
    atomic_data.To_atom_type_CPU(Z, mt::c_Vrl, mt::c_nR, 0.0, atom_type);

    mt::Atom_Cal<T> atomic_fcns_mt;
    atomic_fcns_mt.Set_Atom_Type(potential_type, charge, &atom_type);
    atomic_fcns_mt.feg_dfeg(g.size(), g.data(), feg.data(), dfeg.data());

    return std::make_tuple(feg, dfeg);
  }
  
  template <typename T>
  std::tuple<std::vector<T>, std::vector<T>, std::vector<T>> fxeg_data(int Z, int type) {
    int ng;
    T g[1024], g2[1024], fxg[1024], feg[1024];

    mt::fxeg_Tabulated_Data fxeg_data;
    fxeg_data.ReadTabData(Z, type, 1, ng, g, g2, fxg, feg);
   
    std::vector<T> g_out(ng);
    std::vector<T> fxg_out(ng);
    std::vector<T> feg_out(ng);

    std::copy(g, g + ng, g_out.begin());
    std::copy(fxg, fxg + ng, fxg_out.begin());
    std::copy(feg, feg + ng, feg_out.begin());

    return std::make_tuple(g_out, fxg_out, feg_out);
  }
  
  template <typename T>
  std::tuple<std::vector<T>, std::vector<T>> fxg(
      mt::ePotential_Type potential_type,
      int Z, 
      int charge,
      std::vector<T> g) {
    
    std::vector<T> fxg(g.size());
    std::vector<T> dfxg(g.size());

    mt::Atom_Type<T, mt::e_host> atom_type;
    mt::Atomic_Data atomic_data(potential_type);
    atomic_data.To_atom_type_CPU(Z, mt::c_Vrl, mt::c_nR, 0.0, atom_type);

    mt::Atom_Cal<T> atomic_fcns_mt;
    atomic_fcns_mt.Set_Atom_Type(potential_type, charge, &atom_type);
    atomic_fcns_mt.fxg_dfxg(g.size(), g.data(), fxg.data(), dfxg.data());

    return std::make_tuple(fxg, dfxg);
  }

}


void export_feg(py::module_ m) {
  m.def("feg", &mt::feg<float>);
  m.def("feg", &mt::feg<double>);

  m.def("fxeg_data", &mt::fxeg_data<double>);
  
  m.def("fxg", &mt::fxg<float>);
  m.def("fxg", &mt::fxg<double>);
}

#endif


