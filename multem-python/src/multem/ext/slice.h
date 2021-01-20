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
#ifndef MULTEM_PYTHON_SLICE_H
#define MULTEM_PYTHON_SLICE_H

#include <pybind11/pybind11.h>
#include <multem/multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

// Make the vector of atoms opaque
PYBIND11_MAKE_OPAQUE(std::vector<mt::Slice<float>>);
PYBIND11_MAKE_OPAQUE(std::vector<mt::Slice<double>>);

namespace pybind11 { namespace detail {

  /**
   * Define wrapper function for the mt::Slice<T> class
   */
  template <>
  template <typename T>
  struct Helpers<mt::Slice<T>> {
    /**
     * Get the state
     */
    static py::tuple getstate(const mt::Slice<T> &self) {
      return py::make_tuple(
        self.z_0, self.z_e, self.z_int_0, self.z_int_e, self.iatom_0, self.iatom_e, self.ithk);
    }

    /**
     * Set the state
     */
    static mt::Slice<T> setstate(py::tuple obj) {
      mt::Slice<T> self;
      self.z_0 = obj[0].cast<T>();
      self.z_e = obj[1].cast<T>();
      self.z_int_0 = obj[2].cast<T>();
      self.z_int_e = obj[3].cast<T>();
      self.iatom_0 = obj[4].cast<int>();
      self.iatom_e = obj[5].cast<int>();
      self.ithk = obj[6].cast<int>();
      return self;
    }
  };
  
  /**
   * Define wrapper function for the SliceList class
   */
  template <>
  template <typename T>
  struct Helpers<std::vector<mt::Slice<T>>> {
    /**
     * Get the state
     */
    static py::list getstate(const std::vector<mt::Slice<T>> &self) {
      py::list result;
      for (auto x : self) {
        result.append(x);
      }
      return result;
    }

    /**
     * Set the state
     */
    static std::vector<mt::Slice<T>> setstate(py::list obj) {
      std::vector<mt::Slice<T>> self;
      for (auto x : obj) {
        self.push_back(x.cast<mt::Slice<T>>());
      }
      return self;
    }
  };

}}  // namespace pybind11::detail

template <typename T>
void wrap_slice(py::module_ m, const char *name) {
  typedef mt::Slice<T> Type;

  // Wrap the mt::Input class
  py::class_<Type>(m, name)
    .def(py::init<>())
    .def_readwrite("z_0", &Type::z_0)
    .def_readwrite("z_e", &Type::z_e)
    .def_readwrite("z_int_0", &Type::z_int_0)
    .def_readwrite("z_int_e", &Type::z_int_e)
    .def_readwrite("iatom_0", &Type::iatom_0)
    .def_readwrite("iatom_e", &Type::iatom_e)
    .def_readwrite("ithk", &Type::ithk)
    .def(py::pickle(&py::detail::Helpers<Type>::getstate,
                    &py::detail::Helpers<Type>::setstate));

}

template <typename T>
void wrap_slice_list(py::module_ m, const char *name) {
  typedef std::vector<mt::Slice<T>> Type;
  
  // Register the numpy datatype
  PYBIND11_NUMPY_DTYPE(mt::Slice<T>, z_0, z_e, z_int_0, z_int_e, iatom_0, iatom_e, ithk);

  // Wrap the vector of atoms
  py::bind_vector<Type>(m, name, py::buffer_protocol())
   .def_buffer([](Type &self) -> py::buffer_info {
        return py::buffer_info(
            self.data(),                                    /* Pointer to buffer */
            sizeof(mt::Slice<T>),                           /* Size of one scalar */
            py::format_descriptor<mt::Slice<T>>::format(),  /* Python struct-style format descriptor */
            1,                                              /* Number of dimensions */
            { self.size() },                                /* Buffer dimensions */
            { sizeof(mt::Slice<T>) }                        /* Strides (in bytes) for each index */
        );
    })
    .def(py::pickle(&py::detail::Helpers<Type>::getstate,
                    &py::detail::Helpers<Type>::setstate));

  // Allow implicit conversion
  py::implicitly_convertible<py::list, Type>();
  py::implicitly_convertible<py::array_t<mt::Slice<T>>, Type>();
}

void export_slice(py::module_ m) {
  wrap_slice<float>(m, "Slice_f");
  wrap_slice<double>(m, "Slice_d");
  wrap_slice_list<float>(m, "SliceList_f");
  wrap_slice_list<double>(m, "SliceList_d");
}

#endif

