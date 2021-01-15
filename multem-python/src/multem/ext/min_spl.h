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
#ifndef MULTEM_PYTHON_MIN_SPL
#define MULTEM_PYTHON_MIN_SPL

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace mt {

template <typename T>
std::tuple<int, int> min_spl(T E0, T lx, T ly, T theta) {
  auto g_req = mt::rad_2_rAngs(E0, theta * mt::c_mrad_2_rad);
  auto nx = static_cast<int>(std::round(2 * lx * g_req));
  auto ny = static_cast<int>(std::round(2 * ly * g_req));

  mt::Prime_Num pm;
  nx = pm(nx, mt::eDST_Closest);
  ny = pm(ny, mt::eDST_Closest);

  return std::make_tuple(nx, ny);
}

}  // namespace mt

void export_min_spl(py::module_ m) {
  m.def("min_spl", &mt::min_spl<double>);
}

#endif
