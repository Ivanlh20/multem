/*
 *  grid_2d.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */

#ifndef MULTEM_PYTHON_GRID_2D_H
#define MULTEM_PYTHON_GRID_2D_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {
  
  /**
   * Define wrapper function for the mt::Grid_2d_ class
   */
  template <>
  template <typename T>
  struct Helpers <mt::Grid_2d<T>> {

    /**
     * Get the state
     */
    static py::tuple getstate(const mt::Grid_2d<T> &self) {
      return py::make_tuple(
          self.nx,
          self.ny,
          self.nxh,
          self.nyh,
          self.lx,
          self.ly,
          self.dz,
          self.bwl,
          self.pbc_xy,
          self.Rx_0,
          self.Ry_0,
          self.dRx,
          self.dRy,
          self.dgx,
          self.dgy,
          self.gl2_max,
          self.alpha);
    }

    /**
     * Set the state
     */
    static mt::Grid_2d<T> setstate(py::tuple obj) {
      mt::Grid_2d<T> self;
      self.nx = obj[0].cast<int>();
      self.ny = obj[1].cast<int>();
      self.nxh = obj[2].cast<int>();
      self.nyh = obj[3].cast<int>();
      self.lx = obj[4].cast<T>();
      self.ly = obj[5].cast<T>();
      self.dz = obj[6].cast<T>();
      self.bwl = obj[7].cast<bool>();
      self.pbc_xy = obj[8].cast<bool>();
      self.Rx_0 = obj[9].cast<T>();
      self.Ry_0 = obj[10].cast<T>();
      self.dRx = obj[11].cast<T>();
      self.dRy = obj[12].cast<T>();
      self.dgx = obj[13].cast<T>();
      self.dgy = obj[14].cast<T>();
      self.gl2_max = obj[15].cast<T>();
      self.alpha = obj[16].cast<T>();
      return self;
    }

  };
  

}}


template <typename T>
void wrap_grid_2d(py::module_ m)
{
  typedef mt::Grid_2d<T> Type;

  // Wrap the mt::Input class
  py::class_<Type>(m, "Grid_2d")
    .def(py::init<>())
    .def_readwrite("nx", &Type::nx)
    .def_readwrite("ny", &Type::ny)
    .def_readwrite("nxh", &Type::nxh)
    .def_readwrite("nyh", &Type::nyh)
    .def_readwrite("lx", &Type::lx)
    .def_readwrite("ly", &Type::ly)
    .def_readwrite("dz", &Type::dz)
    .def_readwrite("bwl", &Type::bwl)
    .def_readwrite("pbc_xy", &Type::pbc_xy)
    .def_readwrite("Rx_0", &Type::Rx_0)
    .def_readwrite("Ry_0", &Type::Ry_0)
    .def_readwrite("dRx", &Type::dRx)
    .def_readwrite("dRy", &Type::dRy)
    .def_readwrite("dgx", &Type::dgx)
    .def_readwrite("dgy", &Type::dgy)
    .def_readwrite("gl2_max", &Type::gl2_max)
    .def_readwrite("alpha", &Type::alpha)
    .def(py::pickle(
        &py::detail::Helpers<Type>::getstate,
        &py::detail::Helpers<Type>::setstate))
    ;
}

void export_grid_2d(py::module_ m) {
  wrap_grid_2d<double>(m);
}

#endif




