/*
 *  eels.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */

#ifndef MULTEM_PYTHON_EELS_H
#define MULTEM_PYTHON_EELS_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {
  
  /**
   * Define wrapper function for the mt::EELS class
   */
  template <>
  template <typename T>
  struct Helpers <mt::EELS<T>> {

    /**
     * Get the state
     */
    static py::tuple getstate(const mt::EELS<T> &self) {
      return py::make_tuple(
        self.space,
        self.E_0,
        self.E_loss,
        self.ge,
        self.ge2,
        self.gc,
        self.gc2,
        self.m_selection,
        self.collection_angle,
        self.channelling_type,
        self.factor,
        self.Z,
        self.x,
        self.y,
        self.occ,
        self.g_collection);
    }

    /**
     * Set the state
     */
    static mt::EELS<T> setstate(py::tuple obj) {
      mt::EELS<T> self;
      self.space = obj[0].cast<mt::eSpace>();
      self.E_0 = obj[1].cast<T>();
      self.E_loss = obj[2].cast<T>();
      self.ge = obj[3].cast<T>();
      self.ge2 = obj[4].cast<T>();
      self.gc = obj[5].cast<T>();
      self.gc2 = obj[6].cast<T>();
      self.m_selection = obj[7].cast<int>();
      self.collection_angle = obj[8].cast<T>();
      self.channelling_type = obj[9].cast<mt::eChannelling_Type>();
      self.factor = obj[10].cast<T>();
      self.Z = obj[11].cast<int>();
      self.x = obj[12].cast<T>();
      self.y = obj[13].cast<T>();
      self.occ = obj[14].cast<T>();
      self.g_collection = obj[15].cast<T>();
      return self;
    }

  };
  

}}


template <typename T>
void wrap_eels(py::module_ m)
{
  typedef mt::EELS<T> Type;

  // Wrap the mt::Input class
  py::class_<Type>(m, "EELS")
    .def(py::init<>())
    .def_readwrite("space", &Type::space)
    .def_readwrite("E_0", &Type::E_0)
    .def_readwrite("E_loss", &Type::E_loss)
    .def_readwrite("ge", &Type::ge)
    .def_readwrite("ge2", &Type::ge2)
    .def_readwrite("gc", &Type::gc)
    .def_readwrite("gc2", &Type::gc2)
    .def_readwrite("m_selection", &Type::m_selection)
    .def_readwrite("collection_angle", &Type::collection_angle)
    .def_readwrite("channelling_type", &Type::channelling_type)
    .def_readwrite("factor", &Type::factor)
    .def_readwrite("Z", &Type::Z)
    .def_readwrite("x", &Type::x)
    .def_readwrite("y", &Type::y)
    .def_readwrite("occ", &Type::occ)
    .def_readwrite("g_collection", &Type::g_collection)
    .def(py::pickle(
        &py::detail::Helpers<Type>::getstate,
        &py::detail::Helpers<Type>::setstate))
    ;
}

void export_eels(py::module_ m) {
  wrap_eels<double>(m);
}

#endif





