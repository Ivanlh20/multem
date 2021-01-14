/*
 *  min_spl.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */
#ifndef MULTEM_PYTHON_MIN_SPL
#define MULTEM_PYTHON_MIN_SPL

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt {
  
  template <typename T>
  std::tuple<int,int> min_spl(T E0, T lx, T ly, T theta) {

    auto g_req = mt::rad_2_rAngs(E0, theta * mt::c_mrad_2_rad);
    auto nx = static_cast<int>(std::round(2*lx*g_req));
    auto ny = static_cast<int>(std::round(2*ly*g_req));

    mt::Prime_Num pm;
    nx = pm(nx, mt::eDST_Closest);
    ny = pm(ny, mt::eDST_Closest);

    return std::make_tuple(nx, ny);
 }

}

void export_min_spl(py::module_ m) {
  m.def("min_spl", &mt::min_spl<double>);
}

#endif




