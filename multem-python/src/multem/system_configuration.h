/*
 *  system_configuration.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */
#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  /**
   * Define helper function for the multem::SystemConfiguration class
   */
  template <>
  struct Helpers <mt::SystemConfiguration> {
   
    /**
     * Get the state
     */
    static py::tuple getstate(const mt::SystemConfiguration &self) {
      return py::make_tuple(
        self.get_device(),
        self.get_precision(),
        self.get_cpu_ncores(),
        self.get_cpu_nthread(),
        self.get_gpu_device(),
        self.get_gpu_nstream());
    }

    /**
     * Set the state
     */
    static mt::SystemConfiguration setstate(py::tuple obj) {
      mt::SystemConfiguration self;
      self.set_device(obj[0].cast<mt::eDevice>());
      self.set_precision(obj[1].cast<mt::ePrecision>());
      self.set_cpu_ncores(obj[2].cast<std::size_t>());
      self.set_cpu_nthread(obj[3].cast<std::size_t>());
      self.set_gpu_device(obj[4].cast<std::size_t>());
      self.set_gpu_nstream(obj[5].cast<std::size_t>());
      return self;
    }
  };

}}


template <typename Module>
void export_system_configuration(Module m)
{
  py::class_<mt::SystemConfiguration>(m, "SystemConfiguration")
    .def(py::init<>())
    .def_property("precision", 
        &mt::SystemConfiguration::get_precision,
        &mt::SystemConfiguration::set_precision)
    .def_property("device", 
        &mt::SystemConfiguration::get_device,
        &mt::SystemConfiguration::set_device)
    .def_property("cpu_ncores", 
        &mt::SystemConfiguration::get_cpu_ncores,
        &mt::SystemConfiguration::set_cpu_ncores)
    .def_property("cpu_nthread", 
        &mt::SystemConfiguration::get_cpu_nthread,
        &mt::SystemConfiguration::set_cpu_nthread)
    .def_property("gpu_device", 
        &mt::SystemConfiguration::get_gpu_device,
        &mt::SystemConfiguration::set_gpu_device)
    .def_property("gpu_nstream", 
        &mt::SystemConfiguration::get_gpu_nstream,
        &mt::SystemConfiguration::set_gpu_device)
    .def(py::pickle(
        &py::detail::Helpers<mt::SystemConfiguration>::getstate,
        &py::detail::Helpers<mt::SystemConfiguration>::setstate))
    ;
}


