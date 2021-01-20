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
#ifndef MULTEM_PYTHON_CRYSTAL_BY_LAYERS
#define MULTEM_PYTHON_CRYSTAL_BY_LAYERS

#include <pybind11/pybind11.h>
#include <multem/multem.h>
#include <multem/ext/crystal_parameters.h>

namespace py = pybind11;

namespace mt {

template <typename T>
Atom_Data<T> crystal_by_layers(const CrystalParameters<T> &params) {
  // Get the layer data
  std::vector<mt::Atom_Data<T>> layers(params.layers.size());
  for (auto i = 0; i < params.layers.size(); ++i) {
    layers[i].resize(params.layers[i].size());
    for (auto j = 0; j < params.layers[i].size(); ++j) {
      layers[i].Z[j] = params.layers[i][j].Z;
      layers[i].x[j] = params.layers[i][j].x;
      layers[i].y[j] = params.layers[i][j].y;
      layers[i].z[j] = params.layers[i][j].z;
      layers[i].sigma[j] = params.layers[i][j].sigma;
      layers[i].occ[j] = params.layers[i][j].occ;
      layers[i].region[j] = abs(params.layers[i][j].region);
      layers[i].charge[j] = params.layers[i][j].charge;
    }
  }

  // Get the atoms from the crystal specification
  mt::Crystal_Spec<T> crystal_spec;
  mt::Atom_Data<T> atoms;
  crystal_spec(
    params.na, params.nb, params.nc, params.a, params.b, params.c, layers, atoms);

  return atoms;
}

}  // namespace mt

void export_crystal_by_layers(py::module_ m) {
  m.def("crystal_by_layers", &mt::crystal_by_layers<float>);
  m.def("crystal_by_layers", &mt::crystal_by_layers<double>);
}

#endif
