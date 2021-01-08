/*
 *  atom_data.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */

#ifndef MULTEM_PYTHON_ATOM_DATA_H
#define MULTEM_PYTHON_ATOM_DATA_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {
  
  /**
   * Define wrapper function for the mt::AtomData class
   */
  template <typename T>
  struct AtomDataWrapper : public mt::AtomData<T> {

    /**
     * Get the state
     */
    static py::tuple getstate(const AtomDataWrapper &self) {
      return py::make_tuple(
          self.get_dz(),
          self.get_l_x(),
          self.get_l_y(),
          self.get_l_z(),
          self.get_ct_na(),
          self.get_ct_nb(),
          self.get_ct_nc(),
          self.get_ct_a(),
          self.get_ct_b(),
          self.get_ct_c(),
          self.get_ct_x0(),
          self.get_ct_y0(),
          self.get_amorphous_parameters(),
          self.get_spec_atoms());
    }

    /**
     * Set the state
     */
    static AtomDataWrapper setstate(py::tuple obj) {
      AtomDataWrapper self;
      self.set_dz(obj[0].cast<T>());
      self.set_l_x(obj[1].cast<T>());
      self.set_l_y(obj[2].cast<T>());
      self.set_l_z(obj[3].cast<T>());
      self.set_ct_na(obj[4].cast<int>());
      self.set_ct_nb(obj[5].cast<int>());
      self.set_ct_nc(obj[6].cast<int>());
      self.set_ct_a(obj[7].cast<T>());
      self.set_ct_b(obj[8].cast<T>());
      self.set_ct_c(obj[9].cast<T>());
      self.set_ct_x0(obj[10].cast<T>());
      self.set_ct_y0(obj[11].cast<T>());
      self.set_amorphous_parameters(obj[12].cast<std::vector<mt::Amorp_Lay_Info<T>>>());
      self.set_spec_atoms(obj[13].cast<std::vector<mt::Atom<T>>>());
      self.get_statistic();
      return self;
    }

  };
  
}}


template <typename Module, typename T>
void wrap_atom_data(Module m)
{
  typedef py::detail::AtomDataWrapper<T> Type;

  // Wrap the mt::Input class
  py::class_<Type>(m, "AtomData")
    .def(py::init<>())
    .def_property(
        "dz",
        &Type::get_dz,
        &Type::set_dz)
    .def_property(
        "l_x",
        &Type::get_l_x,
        &Type::set_l_x)
    .def_property(
        "l_y",
        &Type::get_l_y,
        &Type::set_l_y)
    .def_property(
        "l_z",
        &Type::get_l_z,
        &Type::set_l_z)
    .def_property(
        "ct_na",
        &Type::get_ct_na,
        &Type::set_ct_na)
    .def_property(
        "ct_nb",
        &Type::get_ct_nb,
        &Type::set_ct_nb)
    .def_property(
        "ct_nc",
        &Type::get_ct_nc,
        &Type::set_ct_nc)
    .def_property(
        "ct_a",
        &Type::get_ct_a,
        &Type::set_ct_a)
    .def_property(
        "ct_b",
        &Type::get_ct_b,
        &Type::set_ct_b)
    .def_property(
        "ct_c",
        &Type::get_ct_c,
        &Type::set_ct_c)
    .def_property(
        "ct_x0",
        &Type::get_ct_x0,
        &Type::set_ct_x0)
    .def_property(
        "ct_y0",
        &Type::get_ct_y0,
        &Type::set_ct_y0)
    .def_property(
        "amorphous_parameters",
        &Type::get_amorphous_parameters,
        &Type::set_amorphous_parameters)
    .def_property(
        "spec_atoms",
        &Type::get_spec_atoms,
        &Type::set_spec_atoms)
    .def("get_statistic", &Type::get_statistic)
    .def(py::pickle(
        &Type::getstate,
        &Type::setstate))
    ;
}

template <typename Module>
void export_atom_data(Module m) {
  wrap_atom_data<Module, float>(m);
  wrap_atom_data<Module, double>(m);
}

#endif



