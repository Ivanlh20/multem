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
#ifndef MULTEM_PYTHON_FEG
#define MULTEM_PYTHON_FEG

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt {

template <typename T>
std::tuple<std::vector<T>, std::vector<T>> feg(ePotential_Type potential_type,
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
std::tuple<std::vector<T>, std::vector<T>> fxg(mt::ePotential_Type potential_type,
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

template <typename T>
std::tuple<std::vector<T>, std::vector<T>> pr(ePotential_Type potential_type,
                                              int Z,
                                              int charge,
                                              std::vector<T> r) {
  std::vector<T> Pr(r.size());
  std::vector<T> dPr(r.size());

  mt::Atom_Type<T, mt::e_host> atom_type;
  mt::Atomic_Data atomic_data(potential_type);
  atomic_data.To_atom_type_CPU(Z, mt::c_Vrl, mt::c_nR, 0.0, atom_type);

  mt::Atom_Cal<T> atomic_fcns_mt;
  atomic_fcns_mt.Set_Atom_Type(potential_type, charge, &atom_type);
  atomic_fcns_mt.Pr_dPr(r.size(), r.data(), Pr.data(), dPr.data());

  return std::make_tuple(Pr, dPr);
}

template <typename T>
std::tuple<std::vector<T>, std::vector<T>> vp(ePotential_Type potential_type,
                                              int Z,
                                              int charge,
                                              std::vector<T> R) {
  std::vector<T> VR(R.size());
  std::vector<T> dVR(R.size());

  mt::Atom_Type<T, mt::e_host> atom_type;
  mt::Atomic_Data atomic_data(potential_type);
  atomic_data.To_atom_type_CPU(Z, mt::c_Vrl, mt::c_nR, 0.0, atom_type);

  mt::Atom_Cal<T> atomic_fcns_mt;
  atomic_fcns_mt.Set_Atom_Type(potential_type, charge, &atom_type);
  atomic_fcns_mt.VR_dVR(R.size(), R.data(), VR.data(), dVR.data());

  return std::make_tuple(VR, dVR);
}

template <typename T>
std::tuple<std::vector<T>, std::vector<T>> vr(ePotential_Type potential_type,
                                              int Z,
                                              int charge,
                                              std::vector<T> r) {
  std::vector<T> Vr(r.size());
  std::vector<T> dVr(r.size());

  mt::Atom_Type<T, mt::e_host> atom_type;
  mt::Atomic_Data atomic_data(potential_type);
  atomic_data.To_atom_type_CPU(Z, mt::c_Vrl, mt::c_nR, 0.0, atom_type);

  mt::Atom_Cal<T> atomic_fcns_mt;
  atomic_fcns_mt.Set_Atom_Type(potential_type, charge, &atom_type);
  atomic_fcns_mt.Vr_dVr(r.size(), r.data(), Vr.data(), dVr.data());

  return std::make_tuple(Vr, dVr);
}

template <typename T>
std::tuple<std::vector<T>, std::vector<T>> vz(ePotential_Type potential_type,
                                              int Z,
                                              int charge,
                                              double z0,
                                              double ze,
                                              std::vector<T> R) {
  std::vector<T> Vz(R.size());
  std::vector<T> dVz(R.size());

  mt::Atom_Type<T, mt::e_host> atom_type;
  mt::Atomic_Data atomic_data(potential_type);
  atomic_data.To_atom_type_CPU(Z, mt::c_Vrl, mt::c_nR, 0.0, atom_type);

  mt::Atom_Cal<T> atomic_fcns_mt;
  atomic_fcns_mt.Set_Atom_Type(potential_type, charge, &atom_type);
  atomic_fcns_mt.Vz_dVz(z0, ze, R.size(), R.data(), Vz.data(), dVz.data());

  return std::make_tuple(Vz, dVz);
}

}  // namespace mt

void export_feg(py::module_ m) {
  m.def("feg", &mt::feg<float>);
  m.def("feg", &mt::feg<double>);

  m.def("fxeg_data", &mt::fxeg_data<double>);

  m.def("fxg", &mt::fxg<float>);
  m.def("fxg", &mt::fxg<double>);

  m.def("pr", &mt::pr<float>);
  m.def("pr", &mt::pr<double>);

  m.def("vp", &mt::vp<float>);
  m.def("vp", &mt::vp<double>);

  m.def("vr", &mt::vp<float>);
  m.def("vr", &mt::vp<double>);

  m.def("vz", &mt::vz<float>);
  m.def("vz", &mt::vz<double>);
}

#endif
