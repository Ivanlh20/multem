/*
 *  rdf_3d.h
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
std::tuple<std::vector<T>, std::vector<T>> get_rdf_3d(const mt::Atom_Data<T> &atoms,
                                                      T r_max,
                                                      std::size_t nr) {
  std::vector<T> r(nr);
  std::vector<T> rdf(nr);

  mt::rdf_3d(atoms, r_max, nr, r, rdf);

  return std::make_tuple(r, rdf);
}

}  // namespace mt

void export_rdf_3d(py::module_ m) {
  m.def("rdf_3d", &mt::get_rdf_3d<float>);
  m.def("rdf_3d", &mt::get_rdf_3d<double>);
}

#endif
