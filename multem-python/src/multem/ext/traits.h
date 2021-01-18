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
#ifndef MULTEM_PYTHON_TRAITS_H
#define MULTEM_PYTHON_TRAITS_H

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  template <typename T>
  struct type_precision {
  };

  template <>
  struct type_precision <float> {
    static const mt::ePrecision value = mt::eP_float;
  };
  
  template <>
  struct type_precision <double> {
    static const mt::ePrecision value = mt::eP_double;
  };

}}

#endif
