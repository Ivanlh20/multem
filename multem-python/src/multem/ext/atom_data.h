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
#ifndef MULTEM_PYTHON_ATOM_DATA_H
#define MULTEM_PYTHON_ATOM_DATA_H

#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <multem/multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace mt {

template <typename T>
class Atom {
public:
  int Z;
  T x;
  T y;
  T z;
  T sigma;
  T occ;
  int region;
  int charge;

  Atom() : Z(0), x(0), y(0), z(0), sigma(0), occ(0), region(0), charge(0) {}

  Atom(int Z_in, T x_in, T y_in, T z_in, T sigma_in, T occ_in, T region_in, T charge_in)
      : Z(Z_in),
        x(x_in),
        y(y_in),
        z(z_in),
        sigma(sigma_in),
        occ(occ_in),
        region(region_in),
        charge(charge_in) {}
};

}  // namespace mt

// Make the vector of atoms opaque
PYBIND11_MAKE_OPAQUE(std::vector<mt::Atom<float>>);
PYBIND11_MAKE_OPAQUE(std::vector<mt::Atom<double>>);

namespace pybind11 { namespace detail {

  /**
   * Define wrapper function for the AtomList class
   */
  template <typename T>
  struct Helpers<std::vector<mt::Atom<T>>> {
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
   * Define wrapper function for the mt::Atom_Data class
   */
  template <typename T>
  struct Helpers<mt::Atom_Data<T>> {
    /**
     * Get the state
     */
    static py::tuple getstate(const mt::Atom_Data<T> &self) {
      return py::make_tuple(self.dz,
                            self.l_x,
                            self.l_y,
                            self.l_z,
                            self.ct_na,
                            self.ct_nb,
                            self.ct_nc,
                            self.ct_a,
                            self.ct_b,
                            self.ct_c,
                            self.ct_x0,
                            self.ct_y0,
                            self.amorp_lay_info,
                            self.Z,
                            self.x,
                            self.y,
                            self.z,
                            self.sigma,
                            self.occ,
                            self.region,
                            self.charge);
    }

    /**
     * Set the state
     */
    static mt::Atom_Data<T> setstate(py::tuple obj) {
      mt::Atom_Data<T> self;
      self.dz = obj[0].cast<T>();
      self.l_x = obj[1].cast<T>();
      self.l_y = obj[2].cast<T>();
      self.l_z = obj[3].cast<T>();
      self.ct_na = obj[4].cast<int>();
      self.ct_nb = obj[5].cast<int>();
      self.ct_nc = obj[6].cast<int>();
      self.ct_a = obj[7].cast<T>();
      self.ct_b = obj[8].cast<T>();
      self.ct_c = obj[9].cast<T>();
      self.ct_x0 = obj[10].cast<T>();
      self.ct_y0 = obj[11].cast<T>();
      self.amorp_lay_info = obj[12].cast<std::vector<mt::Amorp_Lay_Info<T>>>();
      self.Z = obj[13].cast<std::vector<int>>();
      self.x = obj[14].cast<std::vector<T>>();
      self.y = obj[15].cast<std::vector<T>>();
      self.z = obj[16].cast<std::vector<T>>();
      self.sigma = obj[17].cast<std::vector<T>>();
      self.occ = obj[18].cast<std::vector<T>>();
      self.region = obj[19].cast<std::vector<int>>();
      self.charge = obj[20].cast<std::vector<int>>();
      self.get_statistic();
      return self;
    }

    static std::vector<mt::Atom<T>> get_spec_atoms(const mt::Atom_Data<T> &self) {
      std::vector<mt::Atom<T>> result(self.size());
      for (auto i = 0; i < result.size(); ++i) {
        result[i] = mt::Atom<T>(self.Z[i],
                                self.x[i],
                                self.y[i],
                                self.z[i],
                                self.sigma[i],
                                self.occ[i],
                                self.region[i],
                                self.charge[i]);
      }
      return result;
    }

    template <typename U>
    static void set_spec_atoms_internal(
        mt::Atom_Data<T> &self, const std::vector<mt::Atom<U>> &spec_atoms) {
      self.resize(spec_atoms.size());
      for (auto i = 0; i < spec_atoms.size(); ++i) {
        self.Z[i] = spec_atoms[i].Z;
        self.x[i] = spec_atoms[i].x;
        self.y[i] = spec_atoms[i].y;
        self.z[i] = spec_atoms[i].z;
        self.sigma[i] = spec_atoms[i].sigma;
        self.occ[i] = spec_atoms[i].occ;
        self.region[i] = spec_atoms[i].region;
        self.charge[i] = spec_atoms[i].charge;
      }
    }


    static void set_spec_atoms(mt::Atom_Data<T>& self,
                               const py::object &spec_atoms) {
      try {
        Helpers<mt::Atom_Data<T>>::set_spec_atoms_internal(
            self, 
            spec_atoms.cast<const std::vector<mt::Atom<double>>&>());
      } catch (const std::exception &e) {
        Helpers<mt::Atom_Data<T>>::set_spec_atoms_internal(
            self, 
            spec_atoms.cast<const std::vector<mt::Atom<float>>&>());
      }
    }
  };

  /**
   * Type cast a mt::Atom object to a tuple
   */
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
               src.Z, src.x, src.y, src.z, src.sigma, src.occ, src.region, src.charge)
        .release();
    }
  };
  
  py::object Atom_Data_constructor(const std::string &dtype) {
    MULTEM_ASSERT(dtype == "float" || dtype == "double");
    return (dtype == "float"
      ? py::cast(mt::Atom_Data<float>())
      : py::cast(mt::Atom_Data<double>()));
  }

}}  // namespace pybind11::detail

template <typename T>
void wrap_atom(py::module_ m, const char *name) {
  typedef std::vector<mt::Atom<T>> Type;

  // Register the numpy datatype
  PYBIND11_NUMPY_DTYPE(mt::Atom<T>, Z, x, y, z, sigma, occ, region, charge);

  // Wrap the vector of atoms
  py::bind_vector<Type>(m, name, py::buffer_protocol())
   .def_buffer([](Type &self) -> py::buffer_info {
        return py::buffer_info(
            self.data(),                                  /* Pointer to buffer */
            sizeof(mt::Atom<T>),                          /* Size of one scalar */
            py::format_descriptor<mt::Atom<T>>::format(), /* Python struct-style format descriptor */
            1,                                            /* Number of dimensions */
            { self.size() },                              /* Buffer dimensions */
            { sizeof(mt::Atom<T>) }                       /* Strides (in bytes) for each index */
        );
    })
    .def(py::pickle(&py::detail::Helpers<Type>::getstate,
                    &py::detail::Helpers<Type>::setstate));

  // Allow implicit conversion
  py::implicitly_convertible<py::list, Type>();
  py::implicitly_convertible<py::array_t<mt::Atom<T>>, Type>();


}

template <typename T>
void wrap_atom_data(py::module_ m, const char *name) {
  typedef mt::Atom_Data<T> Type;

  // Wrap the mt::Input class
  py::class_<Type>(m, name)
    .def_readwrite("dz", &Type::dz)
    .def_readwrite("l_x", &Type::l_x)
    .def_readwrite("l_y", &Type::l_y)
    .def_readwrite("l_z", &Type::l_z)
    .def_readwrite("ct_na", &Type::ct_na)
    .def_readwrite("ct_nb", &Type::ct_nb)
    .def_readwrite("ct_nc", &Type::ct_nc)
    .def_readwrite("ct_a", &Type::ct_a)
    .def_readwrite("ct_b", &Type::ct_b)
    .def_readwrite("ct_c", &Type::ct_c)
    .def_readwrite("ct_x0", &Type::ct_x0)
    .def_readwrite("ct_y0", &Type::ct_y0)
    .def_readwrite("amorp_lay_info", &Type::amorp_lay_info)
    .def_readwrite("Z", &Type::Z)
    .def_readwrite("x", &Type::x)
    .def_readwrite("y", &Type::y)
    .def_readwrite("z", &Type::z)
    .def_readwrite("sigma", &Type::sigma)
    .def_readwrite("occ", &Type::occ)
    .def_readwrite("region", &Type::region)
    .def_readwrite("charge", &Type::charge)
    .def_readwrite("Z", &Type::Z)
    .def_readonly("Z_unique", &Type::Z_unique)
    .def_readonly("Z_min", &Type::Z_min)
    .def_readonly("Z_max", &Type::Z_max)
    .def_readonly("x_min", &Type::x_min)
    .def_readonly("x_max", &Type::x_max)
    .def_readonly("y_min", &Type::y_min)
    .def_readonly("y_max", &Type::y_max)
    .def_readonly("z_min", &Type::z_min)
    .def_readonly("z_max", &Type::z_max)
    .def_readonly("sigma_min", &Type::sigma_min)
    .def_readonly("sigma_max", &Type::sigma_max)
    .def_readonly("occ_min", &Type::occ_min)
    .def_readonly("occ_max", &Type::occ_max)
    .def_readonly("region_min", &Type::region_min)
    .def_readonly("region_max", &Type::region_max)
    .def_readonly("R_int_min", &Type::R_int_min)
    .def_readonly("R_int_max", &Type::R_int_max)
    .def_readonly("x_mean", &Type::x_mean)
    .def_readonly("y_mean", &Type::y_mean)
    .def_readonly("z_mean", &Type::z_mean)
    .def_readonly("x_std", &Type::x_std)
    .def_readonly("y_std", &Type::y_std)
    .def_readonly("z_std", &Type::z_std)
    .def_readonly("s_x", &Type::s_x)
    .def_readonly("s_y", &Type::s_y)
    .def_readonly("s_z", &Type::s_z)
    .def_readonly("x_int_min", &Type::x_int_min)
    .def_readonly("x_int_max", &Type::x_int_max)
    .def_readonly("y_int_min", &Type::y_int_min)
    .def_readonly("y_int_max", &Type::y_int_max)
    .def_readonly("z_int_min", &Type::z_int_min)
    .def_readonly("z_int_max", &Type::z_int_max)
    .def_readonly("s_x_int", &Type::s_x_int)
    .def_readonly("s_y_int", &Type::s_y_int)
    .def_readonly("s_z_int", &Type::s_z_int)
    .def_readonly("l_x_int", &Type::l_x_int)
    .def_readonly("l_y_int", &Type::l_y_int)
    .def_readonly("l_z_int", &Type::l_z_int)
    .def_property("spec_atoms",
                  &py::detail::Helpers<Type>::get_spec_atoms,
                  &py::detail::Helpers<Type>::set_spec_atoms)
    .def("get_statistic", static_cast<void (Type::*)()>(&Type::get_statistic))
    .def("sort_by_z", &Type::sort_by_z)
    .def("__len__", &Type::size)
    .def(py::pickle(&py::detail::Helpers<Type>::getstate,
                    &py::detail::Helpers<Type>::setstate));
}

void export_atom_data(py::module_ m) {
  
  m.def(
      "Atom_Data", 
      &py::detail::Atom_Data_constructor, 
      py::arg("dtype") = "double");
  
  wrap_atom<float>(m, "AtomList_f");
  wrap_atom<double>(m, "AtomList_d");
  wrap_atom_data<float>(m, "Atom_Data_f");
  wrap_atom_data<double>(m, "Atom_Data_d");
}

#endif
