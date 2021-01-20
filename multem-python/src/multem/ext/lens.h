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
#ifndef MULTEM_PYTHON_LENS_H
#define MULTEM_PYTHON_LENS_H

#include <pybind11/pybind11.h>
#include <multem/multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  /**
   * Define wrapper function for the mt::Lens class
   */
  template <>
  template <typename T>
  struct Helpers<mt::Lens<T>> {
    /**
     * Get the state
     */
    static py::tuple getstate(const mt::Lens<T> &self) {
      return py::make_tuple(self.m,
                            self.c_10,
                            self.c_12,
                            self.phi_12,
                            self.c_21,
                            self.phi_21,
                            self.c_23,
                            self.phi_23,
                            self.c_30,
                            self.c_32,
                            self.phi_32,
                            self.c_34,
                            self.phi_34,
                            self.c_41,
                            self.phi_41,
                            self.c_43,
                            self.phi_43,
                            self.c_45,
                            self.phi_45,
                            self.c_50,
                            self.c_52,
                            self.phi_52,
                            self.c_54,
                            self.phi_54,
                            self.c_56,
                            self.phi_56,
                            self.inner_aper_ang,
                            self.outer_aper_ang,
                            self.ti_a,
                            self.ti_sigma,
                            self.ti_beta,
                            self.ti_npts,
                            self.ti_iehwgd,
                            self.si_a,
                            self.si_sigma,
                            self.si_beta,
                            self.si_rad_npts,
                            self.si_azm_npts,
                            self.si_iehwgd,
                            self.si_theta_c,
                            self.zero_defocus_type,
                            self.zero_defocus_plane,
                            self.lambda);
    }

    /**
     * Set the state
     */
    static mt::Lens<T> setstate(py::tuple obj) {
      mt::Lens<T> self;
      self.m = obj[0].cast<int>();
      self.c_10 = obj[1].cast<T>();
      self.c_12 = obj[2].cast<T>();
      self.phi_12 = obj[3].cast<T>();
      self.c_21 = obj[4].cast<T>();
      self.phi_21 = obj[5].cast<T>();
      self.c_23 = obj[6].cast<T>();
      self.phi_23 = obj[7].cast<T>();
      self.c_30 = obj[8].cast<T>();
      self.c_32 = obj[9].cast<T>();
      self.phi_32 = obj[10].cast<T>();
      self.c_34 = obj[11].cast<T>();
      self.phi_34 = obj[12].cast<T>();
      self.c_41 = obj[13].cast<T>();
      self.phi_41 = obj[14].cast<T>();
      self.c_43 = obj[15].cast<T>();
      self.phi_43 = obj[16].cast<T>();
      self.c_45 = obj[17].cast<T>();
      self.phi_45 = obj[18].cast<T>();
      self.c_50 = obj[19].cast<T>();
      self.c_52 = obj[20].cast<T>();
      self.phi_52 = obj[21].cast<T>();
      self.c_54 = obj[22].cast<T>();
      self.phi_54 = obj[23].cast<T>();
      self.c_56 = obj[24].cast<T>();
      self.phi_56 = obj[25].cast<T>();
      self.inner_aper_ang = obj[26].cast<T>();
      self.outer_aper_ang = obj[27].cast<T>();
      self.ti_a = obj[28].cast<T>();
      self.ti_sigma = obj[29].cast<T>();
      self.ti_beta = obj[30].cast<T>();
      self.ti_npts = obj[31].cast<int>();
      self.ti_iehwgd = obj[32].cast<T>();
      self.si_a = obj[33].cast<T>();
      self.si_sigma = obj[34].cast<T>();
      self.si_beta = obj[35].cast<T>();
      self.si_rad_npts = obj[36].cast<int>();
      self.si_azm_npts = obj[37].cast<int>();
      self.si_iehwgd = obj[38].cast<T>();
      self.si_theta_c = obj[39].cast<T>();
      self.zero_defocus_type = obj[40].cast<mt::eZero_Defocus_Type>();
      self.zero_defocus_plane = obj[41].cast<T>();
      self.lambda = obj[42].cast<T>();
      return self;
    }
  };
  
  py::object Lens_constructor(const std::string &dtype) {
    MULTEM_ASSERT(dtype == "float" || dtype == "double");
    return (dtype == "float"
      ? py::cast(mt::Lens<float>())
      : py::cast(mt::Lens<double>()));
  }

}}  // namespace pybind11::detail

template <typename T>
void wrap_lens(py::module_ m, const char *name) {
  typedef mt::Lens<T> Type;

  // Wrap the mt::Input class
  py::class_<Type>(m, name)
    .def(py::init<>())
    .def_readwrite("m", &Type::m)
    .def_readwrite("c_10", &Type::c_10)
    .def_readwrite("c_12", &Type::c_12)
    .def_readwrite("phi_12", &Type::phi_12)
    .def_readwrite("c_21", &Type::c_21)
    .def_readwrite("phi_21", &Type::phi_21)
    .def_readwrite("c_23", &Type::c_23)
    .def_readwrite("phi_23", &Type::phi_23)
    .def_readwrite("c_30", &Type::c_30)
    .def_readwrite("c_32", &Type::c_32)
    .def_readwrite("phi_32", &Type::phi_32)
    .def_readwrite("c_34", &Type::c_34)
    .def_readwrite("phi_34", &Type::phi_34)
    .def_readwrite("c_41", &Type::c_41)
    .def_readwrite("phi_41", &Type::phi_41)
    .def_readwrite("c_43", &Type::c_43)
    .def_readwrite("phi_43", &Type::phi_43)
    .def_readwrite("c_45", &Type::c_45)
    .def_readwrite("phi_45", &Type::phi_45)
    .def_readwrite("c_50", &Type::c_50)
    .def_readwrite("c_52", &Type::c_52)
    .def_readwrite("phi_52", &Type::phi_52)
    .def_readwrite("c_54", &Type::c_54)
    .def_readwrite("phi_54", &Type::phi_54)
    .def_readwrite("c_56", &Type::c_56)
    .def_readwrite("phi_56", &Type::phi_56)
    .def_readwrite("inner_aper_ang", &Type::inner_aper_ang)
    .def_readwrite("outer_aper_ang", &Type::outer_aper_ang)
    .def_readwrite("ti_a", &Type::ti_a)
    .def_readwrite("ti_sigma", &Type::ti_sigma)
    .def_readwrite("ti_beta", &Type::ti_beta)
    .def_readwrite("ti_npts", &Type::ti_npts)
    .def_readwrite("ti_iehwgd", &Type::ti_iehwgd)
    .def_readwrite("si_a", &Type::si_a)
    .def_readwrite("si_sigma", &Type::si_sigma)
    .def_readwrite("si_beta", &Type::si_beta)
    .def_readwrite("si_rad_npts", &Type::si_rad_npts)
    .def_readwrite("si_azm_npts", &Type::si_azm_npts)
    .def_readwrite("si_iehwgd", &Type::si_iehwgd)
    .def_readwrite("si_theta_c", &Type::si_theta_c)
    .def_readwrite("zero_defocus_type", &Type::zero_defocus_type)
    .def_readwrite("zero_defocus_plane", &Type::zero_defocus_plane)
    .def_readwrite("lambda_", &Type::lambda)
    .def(py::pickle(&py::detail::Helpers<Type>::getstate,
                    &py::detail::Helpers<Type>::setstate));
}

void export_lens(py::module_ m) {

  m.def(
      "Lens", 
      &py::detail::Lens_constructor, 
      py::arg("dtype") = "double");

  wrap_lens<float>(m, "Lens_f");
  wrap_lens<double>(m, "Lens_d");
}

#endif
