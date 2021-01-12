/*
 *  output_multislice.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */

#ifndef MULTEM_PYTHON_OUTPUT_MULTISLICE_H
#define MULTEM_PYTHON_OUTPUT_MULTISLICE_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {


  /**
   * A class to hold output data
   */
  /* class Data { */
  /* public: */

  /*   std::vector< Image<double> > image_tot; */
  /*   std::vector< Image<double> > image_coh; */
    
  /*   Image<double> m2psi_tot; */
  /*   Image<double> m2psi_coh; */
  /*   Image<double> V; */
  /*   Image< std::complex<double> > psi_coh; */

  /* }; */
  
  /**
   * Get the py::buffer_info from a std::vector
   * @param self A std::vector
   * @returns A py::buffer_info object
   */
  /* template <typename T> */
  /* py::buffer_info as_buffer_info(std::vector<T> &self) { */
  /*   return py::buffer_info( */
  /*     self.data(), */
  /*     sizeof(T), */
  /*     py::format_descriptor<T>::format(), */
  /*     1, */
  /*     { self.size() }, */
  /*     { sizeof(T) }); */
  /* } */
  
  /**
   * Get the py::buffer_info from a multem::Image
   * @param self A multem::Image
   * @returns A py::buffer_info object
   */
  /* template <typename T> */
  /* py::buffer_info as_buffer_info(multem::Image<T> &self) { */
  /*   typedef typename multem::Image<T>::value_type value_type; */
  /*   return py::buffer_info( */
  /*       self.data.data(), */ 
  /*       sizeof(value_type), */ 
  /*       py::format_descriptor<value_type>::format(), */
  /*       2, */
  /*       { */ 
  /*         self.shape[0], */ 
  /*         self.shape[1] */ 
  /*       }, */
  /*       { */ 
  /*         sizeof(value_type) * self.shape[1], */
  /*         sizeof(value_type) */ 
  /*       }); */
  /* } */


  /**
   * Get the py::array_t from a std::vector
   * @param self A std::vector
   * @returns A py::array_t object
   */
  /* template <typename T> */
  /* py::array_t<T> as_array_t(std::vector<T> &self) { */
  /*   return py::array_t<T>(as_buffer_info(self)); */
  /* } */

  /**
   * Get the py::array_t from a multem::Image
   * @param self A multem::Image
   * @returns A py::array_t object
   */
  /* template <typename T> */
  /* py::array_t<T> as_array_t(multem::Image<T> &self) { */
  /*   return py::array_t<T>(as_buffer_info(self)); */
  /* } */
  
  /**
   * Define helper functions for the mt::DetInt class
   */
  //template <>
  //template <typename T>
  //struct Helpers < mt::DetInt<T> > {

  //  /**
  //   * Create a mt::DetInt from a py::array_t
  //   * @param array The py::array_t object
  //   * @returns The multem::Image object
  //   */
  //  static mt::DetInt<T> init_from_array_t(py::array_t<T> array) {
  //    py::buffer_info buffer = array.request();
  //    /* MULTEM_ASSERT(buffer.ndim == 2); */
  //    /* MULTEM_ASSERT(buffer.shape[0] >= 0); */
  //    /* MULTEM_ASSERT(buffer.shape[1] >= 0); */
  //    return mt::DetInt<T>(
  //      (T *) buffer.ptr, 
  //      typename mt::DetInt<T>::shape_type({
  //        (std::size_t) buffer.shape[0], 
  //        (std::size_t) buffer.shape[1]}));
  //  }

  //  /**
  //   * Create a py::buffer_info object from a mt::DetInt object
  //   * @param self The mt::DetInt object
  //   * @returns The py::buffer_info object
  //   */
  //  static py::buffer_info as_buffer_info(mt::DetInt<T> &self) {
  //    return py::detail::as_buffer_info(self);
  //  }
  //};


  /**
   * Define helper function for the multem::SystemConfiguration class
   */
  template <>
  template <typename T>
  struct Helpers <mt::Output<T>> {

    /**
     * Get the state
     */
    static py::tuple getstate(const mt::Output<T> &self) {
      return py::make_tuple(
        self.get_output_type(),
        self.get_ndetector(),
        self.get_nx(),
        self.get_ny(),
        self.get_dx(),
        self.get_dy(),
        self.get_dr(),
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
        self.get_psi_0(),
        self.get_thk_gpu());
    }

    /**
     * Set the state
     */
    static mt::Output<T> setstate(py::tuple obj) {
      typedef typename mt::Output<T>::vector_type vector_type;
      typedef typename mt::Output<T>::complex_vector_type complex_vector_type;
      mt::Output<T> self;
      self.set_output_type(obj[0].cast<mt::eTEM_Output_Type>());
      self.set_ndetector(obj[1].cast<int>());
      self.set_nx(obj[2].cast<int>());
      self.set_ny(obj[3].cast<int>());
      self.set_dx(obj[4].cast<T>());
      self.set_dy(obj[5].cast<T>());
      self.set_dr(obj[6].cast<T>());
      self.set_x(obj[7].cast<std::vector<T>>());
      self.set_y(obj[10].cast<std::vector<T>>());
      self.set_r(obj[11].cast<std::vector<T>>());
      self.set_image_tot(obj[12].cast<std::vector<mt::DetInt<vector_type>>>());
      self.set_image_coh(obj[13].cast<std::vector<mt::DetInt<vector_type>>>());
      self.set_m2psi_tot(obj[14].cast<std::vector<vector_type>>());
      self.set_m2psi_coh(obj[15].cast<std::vector<vector_type>>());
      self.set_psi_coh(obj[16].cast<std::vector<complex_vector_type>>());
      self.set_V(obj[17].cast<std::vector<vector_type>>());
      self.set_trans(obj[18].cast<std::vector<complex_vector_type>>());
      self.set_psi_0(obj[19].cast<std::vector<complex_vector_type>>());
      self.set_thk_gpu(obj[20].cast<std::vector<bool>>());
      return self;
    }
  };

}}

template <typename T>
void wrap_det_int(const py::handle &m, const char *name)
{
  /* return py::class_< mt::DetInt<T> >(m, name, py::buffer_protocol()) */
  /*     .def(py::init<>()) */
  /*     .def(py::init(&py::detail::Helpers<mt::DetInt<T>>::init_from_array_t)) */
  /*     .def_buffer(&py::detail::Helpers<mt::DetInt<T>>::as_buffer_info) */
  /*     ; */
}

template <typename T>
void wrap_output_multislice(const py::handle &m, const char *name)
{
  typedef mt::Input<T> Parent;
  typedef mt::Output<T> Type;

  py::class_<Type, Parent>(m, name)
    .def(py::init<>())
    .def_property(
        "output_type", 
        &Type::get_output_type,
        &Type::set_output_type)
    .def_property(
        "ndetector", 
        &Type::get_ndetector,
        &Type::set_ndetector)
    .def_property(
        "nx", 
        &Type::get_nx,
        &Type::set_nx)
    .def_property(
        "ny", 
        &Type::get_ny,
        &Type::set_ny)
    .def_property(
        "dx", 
        &Type::get_dx,
        &Type::set_dx)
    .def_property(
        "dy", 
        &Type::get_dy,
        &Type::set_dy)
    .def_property(
        "dr", 
        &Type::get_dr,
        &Type::set_dr)
    .def_property(
        "x", 
        &Type::get_x,
        &Type::set_x)
    .def_property(
        "y", 
        &Type::get_y,
        &Type::set_y)
    .def_property(
        "r", 
        &Type::get_r,
        &Type::set_r)
    .def_property(
        "image_tot", 
        &Type::get_image_tot,
        &Type::set_image_tot)
    .def_property(
        "image_coh", 
        &Type::get_image_coh,
        &Type::set_image_coh)
    .def_property(
        "m2psi_tot", 
        &Type::get_m2psi_tot,
        &Type::set_m2psi_tot)
    .def_property(
        "m2psi_coh", 
        &Type::get_m2psi_coh,
        &Type::set_m2psi_coh)
    .def_property(
        "psi_coh", 
        &Type::get_psi_coh,
        &Type::set_psi_coh)
    .def_property(
        "V", 
        &Type::get_V,
        &Type::set_V)
    .def_property(
        "trans", 
        &Type::get_trans,
        &Type::set_trans)
    .def_property(
        "psi_0", 
        &Type::get_psi_0,
        &Type::set_psi_0)
    .def_property(
        "thk_gpu", 
        &Type::get_thk_gpu,
        &Type::set_thk_gpu)
    /* .def_readonly("data", &Type::data) */
    .def("gather", &Type::gather)
    .def("clean_temporal", &Type::clean_temporal)
    .def(py::pickle(
        &py::detail::Helpers<Type>::getstate,
        &py::detail::Helpers<Type>::setstate))
    ;
}

void export_output_multislice(const py::handle &m) {
  wrap_det_int<double>(m, "DetIntD");
  /* wrap_output_multislice<float>(m, "OutputF"); */
  wrap_output_multislice<double>(m, "OutputD");
}

#endif

