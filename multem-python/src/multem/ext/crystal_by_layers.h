/*
 *  crystal_by_layers.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */
#ifndef MULTEM_PYTHON_CRYSTAL_BY_LAYERS
#define MULTEM_PYTHON_CRYSTAL_BY_LAYERS

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt {

  template <typename T>
  Atom_Data<T> crystal_by_layers(const CrystalParameters<T> &params) {
    
    // Get the layer data 
    std::vector<mt::Atom_Data<double>> layers(params.layers.size());
    for(auto i = 0; i < params.layers.size(); ++i) {
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
    mt::Crystal_Spec<double> crystal_spec;
    mt::Atom_Data<double> atoms;
    crystal_spec(
      params.na, 
      params.nb, 
      params.nc, 
      params.a, 
      params.b, 
      params.c, 
      layers, 
      atoms);

    return atoms;
  }

}


void export_crystal_by_layers(py::module_ m) {
  m.def("crystal_by_layers", &mt::crystal_by_layers<float>);
  m.def("crystal_by_layers", &mt::crystal_by_layers<double>);
}

#endif


