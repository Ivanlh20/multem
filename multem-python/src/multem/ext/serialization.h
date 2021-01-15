/*
 *  serialization.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of
 *  which is included in the root directory of this package.
 */

#ifndef MULTEM_PYTHON_SERIALIZATION_H
#define MULTEM_PYTHON_SERIALIZATION_H

#include <pybind11/pybind11.h>
#include <multem.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  /**
   * Define helper functions for a type
   */
  template <typename T>
  struct Helpers {};

}}  // namespace pybind11::detail

#endif
