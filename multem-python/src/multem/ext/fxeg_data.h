/*
 *  fxeg_data.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */
#ifndef MULTEM_FXEG_DATA
#define MULTEM_FXEG_DATA

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace pybind { namespace detail {

  template <typename T>
  py::object fxeg_data(int Z, int type) {
    int ng;
    T g[1024], g2[1024], fxg[1024], feg[1024];

    mt::fxeg_Tabulated_Data fxeg_data;
    fxeg_data.ReadTabData(Z, type, 1, ng, g, g2, fxg, feg);
   
    std::vector<T> g_out(ng);
    std::vector<T> fxg_out(ng);
    std::vector<T> feg_out(ng);

    std::copy(g, g + ng, g_o.begin());
    std::copy(fxg, fxg + ng, fxg_o.begin());
    std::copy(feg, feg + ng, feg_o.begin());

    return py::make_tuple(g, fxg, feg);
  }

}}


void export_fxeg_data(py::module_ m) {
  m.def("fxeg_data", &mt::fxeg_data<double>);
}

#endif




