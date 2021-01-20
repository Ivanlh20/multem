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
#ifndef MULTEM_PYTHON_OUTPUT_MULTISLICE_H
#define MULTEM_PYTHON_OUTPUT_MULTISLICE_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  /**
   * Define helper function for the multem::Det_Int class
   */
  template <>
  template <typename T>
  struct Helpers<mt::Det_Int<T>> {
    /**
     * Get the state
     */
    static py::tuple getstate(const mt::Det_Int<T> &self) {
      return py::cast(self.image);
    }

    /**
     * Set the state
     */
    static mt::Det_Int<T> setstate(py::tuple obj) {
      mt::Det_Int<T> self;
      self.image = obj.cast<std::vector<T>>();
      return self;
    }
  };

  /**
   * Define helper function for the multem::Output_Multislice class
   */
  template <>
  template <typename T>
  struct Helpers<mt::Output_Multislice<T>> {
    /**
     * Get the state
     */
    static py::tuple getstate(const mt::Output_Multislice<T> &self) {
      return py::make_tuple(self.output_type(),
                            self.ndetector(),
                            self.nx(),
                            self.ny(),
                            self.dx(),
                            self.dy(),
                            self.dr(),
                            self.get_x(),
                            self.get_y(),
                            self.get_r(),
                            self.get_image_tot(),
                            self.get_image_coh(),
                            self.get_m2psi_tot(),
                            self.get_m2psi_coh(),
                            self.get_psi_coh(),
                            self.get_V(),
                            self.get_trans(),
                            self.get_psi_0());
    }

    /**
     * Set the state
     */
    static mt::Output_Multislice<T> setstate(py::tuple obj) {
      typedef std::vector<T> vector_type;
      typedef std::vector<std::complex<T>> complex_vector_type;
      mt::Output_Multislice<T> self;
      self.output_type() = obj[0].cast<mt::eTEM_Output_Type>();
      self.ndetector() = obj[1].cast<int>();
      self.nx() = obj[2].cast<int>();
      self.ny() = obj[3].cast<int>();
      self.dx() = obj[4].cast<T>();
      self.dy() = obj[5].cast<T>();
      self.dr() = obj[6].cast<T>();
      self.set_x(obj[7].cast<std::vector<T>>());
      self.set_y(obj[8].cast<std::vector<T>>());
      self.set_r(obj[9].cast<std::vector<T>>());
      self.set_image_tot(obj[10].cast<std::vector<mt::Det_Int<vector_type>>>());
      self.set_image_coh(obj[11].cast<std::vector<mt::Det_Int<vector_type>>>());
      self.set_m2psi_tot(obj[12].cast<std::vector<vector_type>>());
      self.set_m2psi_coh(obj[13].cast<std::vector<vector_type>>());
      self.set_psi_coh(obj[14].cast<std::vector<complex_vector_type>>());
      self.set_V(obj[15].cast<std::vector<vector_type>>());
      self.set_trans(obj[16].cast<std::vector<complex_vector_type>>());
      self.set_psi_0(obj[17].cast<std::vector<complex_vector_type>>());
      return self;
    }

    static mt::eTEM_Output_Type get_output_type(const mt::Output_Multislice<T> &self) {
      return self.output_type();
    }

    static void set_output_type(mt::Output_Multislice<T> &self,
                                mt::eTEM_Output_Type output_type) {
      self.output_type() = output_type;
    }

    static int get_ndetector(const mt::Output_Multislice<T> &self) {
      return self.ndetector();
    }

    static void set_ndetector(mt::Output_Multislice<T> &self, int ndetector) {
      self.ndetector() = ndetector;
    }

    static int get_nx(const mt::Output_Multislice<T> &self) {
      return self.nx();
    }

    static void set_nx(mt::Output_Multislice<T> &self, int nx) {
      self.nx() = nx;
    }

    static int get_ny(const mt::Output_Multislice<T> &self) {
      return self.ny();
    }

    static void set_ny(mt::Output_Multislice<T> &self, int ny) {
      self.ny() = ny;
    }

    static T get_dx(const mt::Output_Multislice<T> &self) {
      return self.dx();
    }

    static void set_dx(mt::Output_Multislice<T> &self, T dx) {
      self.dx() = dx;
    }

    static T get_dy(const mt::Output_Multislice<T> &self) {
      return self.dy();
    }

    static void set_dy(mt::Output_Multislice<T> &self, T dy) {
      self.dy() = dy;
    }

    static T get_dr(const mt::Output_Multislice<T> &self) {
      return self.dr();
    }

    static void set_dr(mt::Output_Multislice<T> &self, T dr) {
      self.dr() = dr;
    }
  };
  
  py::object Output_Multislice_constructor(const std::string &dtype) {
    MULTEM_ASSERT(dtype == "float" || dtype == "double");
    return (dtype == "float"
      ? py::cast(mt::Output_Multislice<float>())
      : py::cast(mt::Output_Multislice<double>()));
  }


}}  // namespace pybind11::detail

template <typename T>
void wrap_det_int(py::module_ m, const char *name) {

  typedef std::vector<T> ValueType;
  typedef mt::Det_Int<ValueType> Type;

  py::class_<Type>(m, name, py::buffer_protocol())
    .def(py::init<>())
    .def(py::init([](py::iterable it) {
        auto self = std::unique_ptr<Type>(new Type());
        self->image.reserve(py::len_hint(it));
        for (py::handle h : it)
           self->image.push_back(h.cast<ValueType>());
        return self.release();
    }))
    .def_readwrite("image", &Type::image)
    .def(py::pickle(&py::detail::Helpers<Type>::getstate,
                    &py::detail::Helpers<Type>::setstate))
    ;

  // Allow implicit conversion
  py::implicitly_convertible<py::list, Type>();
  py::implicitly_convertible<py::array, Type>();
}

template <typename T>
void wrap_output_multislice(py::module_ m, const char *name) {
  typedef mt::Input_Multislice<T> Parent;
  typedef mt::Output_Multislice<T> Type;

  py::class_<Type, Parent>(m, name)
    .def(py::init<>())
    .def_property("output_type",
                  &py::detail::Helpers<Type>::get_output_type,
                  &py::detail::Helpers<Type>::set_output_type)
    .def_property("ndetector",
                  &py::detail::Helpers<Type>::get_ndetector,
                  &py::detail::Helpers<Type>::set_ndetector)
    .def_property(
      "nx", &py::detail::Helpers<Type>::get_nx, &py::detail::Helpers<Type>::set_nx)
    .def_property(
      "ny", &py::detail::Helpers<Type>::get_ny, &py::detail::Helpers<Type>::set_ny)
    .def_property(
      "dx", &py::detail::Helpers<Type>::get_dx, &py::detail::Helpers<Type>::set_dx)
    .def_property(
      "dy", &py::detail::Helpers<Type>::get_dy, &py::detail::Helpers<Type>::set_dy)
    .def_property(
      "dr", &py::detail::Helpers<Type>::get_dr, &py::detail::Helpers<Type>::set_dr)
    .def_property("x", &Type::get_x, &Type::set_x)
    .def_property("y", &Type::get_y, &Type::set_y)
    .def_property("r", &Type::get_r, &Type::set_r)
    .def_property("image_tot", &Type::get_image_tot, &Type::set_image_tot)
    .def_property("image_coh", &Type::get_image_coh, &Type::set_image_coh)
    .def_property("m2psi_tot", &Type::get_m2psi_tot, &Type::set_m2psi_tot)
    .def_property("m2psi_coh", &Type::get_m2psi_coh, &Type::set_m2psi_coh)
    .def_property(
      "psi_coh",
      &Type::get_psi_coh,
      static_cast<void (Type::*)(const std::vector<std::vector<std::complex<T>>> &)>(
        &Type::set_psi_coh))
    .def_property("V", &Type::get_V, &Type::set_V)
    .def_property("trans", &Type::get_trans, &Type::set_trans)
    .def_property("psi_0", &Type::get_psi_0, &Type::set_psi_0)
    .def("gather", &Type::gather)
    .def("clean_temporal", &Type::clean_temporal)
    .def(py::pickle(&py::detail::Helpers<Type>::getstate,
                    &py::detail::Helpers<Type>::setstate));
}

void export_output_multislice(py::module_ m) {

  m.def(
      "Output_Multislice", 
      &py::detail::Output_Multislice_constructor, 
      py::arg("dtype") = "double");

  wrap_det_int<float>(m, "Det_Int_f");
  wrap_det_int<double>(m, "Det_Int_d");
  wrap_output_multislice<float>(m, "Output_Multislice_f");
  wrap_output_multislice<double>(m, "Output_Multislice_d");
}

#endif
