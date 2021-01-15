/*
 *  amorp_spec.h
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
Atom_Data<T> amorp_spec(T l_z, T r_min, int Z, T rms_3d, T rho, int seed = 300183) {
  mt::Atom_Data<T> atoms;

  atoms.amorp_lay_info.resize(1);
  atoms.amorp_lay_info[0].z_0 = 0;
  atoms.amorp_lay_info[0].z_e = l_z;
  atoms.amorp_lay_info[0].dz = 2.0;
  atoms.amorp_lay_info[0].type = mt::eALT_Bottom;
  atoms.amorp_lay_info[0].region = 0;

  mt::Amorp_Spec<T> spec;
  spec.create(atoms, r_min, Z, rms_3d, rho, seed);

  return atoms;
}

}  // namespace mt

void export_amorp_spec(py::module_ m) {
  m.def("amorp_spec", &mt::amorp_spec<float>);
  m.def("amorp_spec", &mt::amorp_spec<double>);
}

#endif
