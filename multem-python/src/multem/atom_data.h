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
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <multem.h>
#include <multem/serialization.h>

namespace py = pybind11;

// Make the vector of atoms opaque
PYBIND11_MAKE_OPAQUE(std::vector<mt::Atom<double>>);

namespace pybind11 { namespace detail {
  
  /**
   * Define wrapper function for the AtomList class
   */
  template <>
  template <typename T>
  struct Helpers <std::vector<mt::Atom<T>>> {

    /**
     * Get the state
     */
    static py::list getstate(const std::vector<mt::Atom<T>> &self) {
      py::list result;
      for (auto x : self) {
        result.append(x);
      }
      return result;
    }

    /**
     * Set the state
     */
    static std::vector<mt::Atom<T>> setstate(py::list obj) {
      std::vector<mt::Atom<T>> self;
      for (auto x : obj) {
        self.push_back(x.cast<mt::Atom<T>>());
      }
      return self;
    }

  };
  
  /**
   * Define wrapper function for the mt::AtomData class
   */
  template <>
  template <typename T>
  struct Helpers <mt::AtomData<T>> {

    /**
     * Get the state
     */
    static py::tuple getstate(const mt::AtomData<T> &self) {
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
    static void setstate(mt::AtomData<T>& self, py::tuple obj) {
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
    }

  };
  

  /**
   * Type cast a mt::Atom object to a tuple
   */
  template <> 
  template <typename T>
  class type_caster<mt::Atom<T>> {
  public:
  
    PYBIND11_TYPE_CASTER(mt::Atom<T>, _("mt::Atom<T>"));

    bool load(handle src, bool convert) {
      if (py::isinstance<py::tuple>(src)) {
        py::tuple t = py::cast<py::tuple>(src);
        if (py::len(t) == 8) {
          value.Z = py::cast<int>(t[0]);
          value.x = py::cast<double>(t[1]);
          value.y = py::cast<double>(t[2]);
          value.z = py::cast<double>(t[3]);
          value.sigma = py::cast<double>(t[4]);
          value.occ = py::cast<double>(t[5]);
          value.region = py::cast<int>(t[6]);
          value.charge = py::cast<int>(t[7]);
          return true;
        }
      }
      return false;
    }

    static handle cast(mt::Atom<T> src, return_value_policy policy, handle parent) {
      return py::make_tuple(
        src.Z, 
        src.x, 
        src.y, 
        src.z, 
        src.sigma, 
        src.occ, 
        src.region, 
        src.charge).release();
    }
  };
  
}}

template <typename Module, typename T>
void wrap_atom(Module m)
{
  typedef std::vector<mt::Atom<T>> Type;
  
  // Wrap the vector of atoms
  py::bind_vector<Type>(m, "AtomList")
    .def(py::pickle(
        &py::detail::Helpers<Type>::getstate,
        &py::detail::Helpers<Type>::setstate))
    ;

  // Allow implicit conversion
  py::implicitly_convertible<py::list, Type>();
}

template <typename Module, typename T>
void wrap_atom_data(Module m)
{
  typedef mt::AtomData<T> Type;

  // Wrap the mt::Input class
  py::class_<Type>(m, "AtomData")
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
    ;
}

template <typename Module>
void export_atom_data(Module m) {
  wrap_atom<Module, double>(m);
  wrap_atom_data<Module, double>(m);
}

#endif



