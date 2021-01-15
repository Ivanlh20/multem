/*
 *  crystal_parameters.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of
 *  which is included in the root directory of this package.
 */
#ifndef MULTEM_PYTHON_CRYSTAL_PARAMETERS
#define MULTEM_PYTHON_CRYSTAL_PARAMETERS

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt {

/**
 * A class to hold crystal parameters
 */
template <typename T>
class CrystalParameters {
public:
  typedef std::vector<Atom<T>> Layer;

  int na;
  int nb;
  int nc;
  T a;
  T b;
  T c;
  std::vector<Layer> layers;

  CrystalParameters() : na(0), nb(0), nc(0), a(0), b(0), c(0) {}
};

}  // namespace mt

void export_crystal_parameters(py::module_ m) {
  typedef mt::CrystalParameters<double> Type;

  py::class_<Type>(m, "CrystalParameters")
    .def(py::init<>())
    .def_readwrite("na", &Type::na)
    .def_readwrite("nb", &Type::nb)
    .def_readwrite("nc", &Type::nc)
    .def_readwrite("a", &Type::a)
    .def_readwrite("b", &Type::b)
    .def_readwrite("c", &Type::c)
    .def_readwrite("layers", &Type::layers);
}

#endif
