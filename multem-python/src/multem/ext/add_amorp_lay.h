/*
 *  add_amorp_lay.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */
#ifndef MULTEM_PYTHON_ADD_AMORP_LAY
#define MULTEM_PYTHON_ADD_AMORP_LAY

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt {

  template <typename T>
  Atom_Data<T> add_amorp_lay(
      Atom_Data<T> atoms, 
      T l_z,
      T r_min,
      int Z,
      T rms_3d,
      T rho,
      mt::eAmorp_Lay_Type lay_pos = mt::eALT_Top,
      int seed = 300183) {

    auto region = atoms.region_max+1;

    atoms.amorp_lay_info.resize(1);
    if(lay_pos == mt::eALT_Top)
    {
      atoms.amorp_lay_info[0].z_0 = atoms.z_min-l_z;
      atoms.amorp_lay_info[0].z_e = atoms.z_min;
      atoms.amorp_lay_info[0].dz = 2.0;
      atoms.amorp_lay_info[0].type = mt::eALT_Top;
      atoms.amorp_lay_info[0].region = region;
    }
    else
    {
      atoms.amorp_lay_info[0].z_0 = atoms.z_max;
      atoms.amorp_lay_info[0].z_e = atoms.z_max+l_z;
      atoms.amorp_lay_info[0].dz = 2.0;
      atoms.amorp_lay_info[0].type = mt::eALT_Bottom;
      atoms.amorp_lay_info[0].region = region;
    }

    mt::Amorp_Spec<T> spec;
    spec.create(atoms, r_min, Z, rms_3d, rho, seed);

    return atoms;
  }

}


void export_add_amorp_lay(py::module_ m) {
  m.def("add_amorp_lay", &mt::add_amorp_lay<float>);
  m.def("add_amorp_lay", &mt::add_amorp_lay<double>);
}

#endif
