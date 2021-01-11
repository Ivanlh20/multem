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
   * A class to represent an image
   */
  template <typename T>
  class Image {
  public:

    typedef T value_type;
    typedef std::array<std::size_t, 2> shape_type;

    std::vector<value_type> data;
    shape_type shape;

    Image()
      : shape({ 0, 0 }) {}

    /**
     * Construct the image from the pointer
     * @param data_ The data pointer
     * @param shape_ The image shape (Y, X)
     */
    template <typename U>
    Image(const U *data_, shape_type shape_)
      : shape(shape_) {
      std::size_t size = shape[0]*shape[1];
      data.assign(data_, data_ + size);
    }
  };

  /**
   * A class to represent a complex image
   */
  template <typename T>
  class Image< std::complex<T> > {
  public:
    
    typedef std::complex<T> value_type;
    typedef std::array<std::size_t, 2> shape_type;

    std::vector<value_type> data;
    shape_type shape;

    Image()
      : shape({ 0, 0 }) {}

    /**
     * Construct the image from the pointer
     * @param data_ The data pointer
     * @param shape_ The image shape (Y, X)
     */
    template <typename U>
    Image(const U *data_, shape_type shape_)
      : shape(shape_) {
      std::size_t size = shape[0]*shape[1];
      data.resize(size);
      for (auto i = 0; i < size; ++i) {
        data[i] = value_type(data_[i].real(), data_[i].imag());
      }
    }
    
  };

  /**
   * A class to hold output data
   */
  class Data {
  public:

    std::vector< Image<double> > image_tot;
    std::vector< Image<double> > image_coh;
    
    Image<double> m2psi_tot;
    Image<double> m2psi_coh;
    Image<double> V;
    Image< std::complex<double> > psi_coh;

  };


  /**
   * Define helper function for the multem::SystemConfiguration class
   */
  template <typename T>
  struct OutputWrapper : public mt::Output<T> {

    std::vector<Data> data;

    /**
     * Get the state
     */
    static py::tuple getstate(const OutputWrapper &self) {
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
    static OutputWrapper setstate(py::tuple obj) {
      OutputWrapper self;
      /* self.set_output_type(obj[0].cast<mt::eTEM_Output_Type>()); */
      /* self.set_ndetector(obj[1].cast<int>()); */
      /* self.set_nx(obj[2].cast<int>()); */
      /* self.set_ny(obj[3].cast<int>()); */
      /* self.set_dx(obj[4].cast<T>()); */
      /* self.set_dy(obj[5].cast<T>()); */
      /* self.set_dr(obj[6].cast<T>()); */
      /* self.set_x(obj[7].cast<std::vector<T>>()); */
      /* self.set_y(obj[10].cast<std::vector<T>>()); */
      /* self.set_r(obj[11].cast<std::vector<T>>()); */
      /* self.set_image_tot(obj[12].cast<std::vector<mt::DetInt<vector_type>>>()); */
      /* self.set_image_coh(obj[13].cast<std::vector<mt::DetInt<vector_type>>>()); */
      /* self.set_m2psi_tot(obj[14].cast<vector_type>()); */
      /* self.set_m2psi_coh(obj[15].cast<vector_type>()); */
      /* self.set_psi_coh(obj[16].cast<complex_vector_type>()); */
      /* self.set_V(obj[17].cast<vector_type>()); */
      /* self.set_trans(obj[18].cast<complex_vector_type>()); */
      /* self.set_psi_0(obj[19].cast<complex_vector_type>()); */
      /* self.set_thk_gpu(obj[20].cast<std::vector<bool>>()); */
      return self;
    }
  };

}}


template <typename Module, typename T>
void wrap_output_multislice(Module m)
{
  typedef mt::Input<T> Parent;
  typedef py::detail::OutputWrapper<T> Type;

  py::class_<Type, Parent>(m, "Output")
    .def(py::init<>())
    .def_property_readonly("output_type", &Type::get_output_type)
    .def_property_readonly("ndetector", &Type::get_ndetector)
    .def_property_readonly("nx", &Type::get_nx)
    .def_property_readonly("ny", &Type::get_ny)
    .def_property_readonly("dx", &Type::get_dx)
    .def_property_readonly("dy", &Type::get_dy)
    .def_property_readonly("dr", &Type::get_dr)
    .def_property_readonly("x", &Type::get_x)
    .def_property_readonly("y", &Type::get_y)
    .def_property_readonly("r", &Type::get_r)
    .def_property_readonly("image_tot", &Type::get_image_tot)
    .def_property_readonly("image_coh", &Type::get_image_coh)
    .def_property_readonly("m2psi_tot", &Type::get_m2psi_tot)
    .def_property_readonly("m2psi_coh", &Type::get_m2psi_coh)
    .def_property_readonly("psi_coh", &Type::get_psi_coh)
    .def_property_readonly("V", &Type::get_V)
    .def_property_readonly("trans", &Type::get_trans)
    .def_property_readonly("psi_0", &Type::get_psi_0)
    .def_property_readonly("thk_gpu", &Type::get_thk_gpu)
    .def_readonly("data", &Type::data)
    .def("gather", &Type::gather)
    .def("clean_temporal", &Type::clean_temporal)
    .def(py::pickle(
        &Type::getstate,
        &Type::setstate))
    ;
}

template <typename Module>
void export_output_multislice(Module m) {
  wrap_output_multislice<Module, double>(m);
}

#endif

