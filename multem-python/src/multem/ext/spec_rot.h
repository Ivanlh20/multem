/*
 *  spec_rot.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */
#ifndef MULTEM_PYTHON_SPEC_ROT
#define MULTEM_PYTHON_SPEC_ROT

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt {

  template <typename T>
  void spec_rot(
      Atom_Data<T> &atoms,
      double theta,
      const mt::r3d<T> &u0,
      mt::eRot_Point_Type &rot_point_type,
      const mt::r3d<T> &p0) {

    u0.normalized();
	
    if(rot_point_type == mt::eRPT_geometric_center) {
      p0 = r3d<T>(atoms.x_mean, atoms.y_mean, atoms.z_mean);
    }

    mt::rotate_atoms(atoms, theta * mt::c_deg_2_rad, u0, p0);
  }

}


void export_spec_rot(py::module_ m) {
  m.def("spec_rot", &mt::spec_rot<float>);
  m.def("spec_rot", &mt::spec_rot<double>);
}

#endif


