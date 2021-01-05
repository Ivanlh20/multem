/*
 * This file is part of MULTEM.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
 * Copyright 2021 Diamond Light Source
 * Copyright 2021 Rosalind Franklin Institute
 *
 * Author: James Parkhurst
 *
 * MULTEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MULTEM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#include <iostream>
#include <multem.h>
#include <types.cuh>

namespace mt {

  struct SystemConfiguration::Impl {
    System_Configuration data;
  };

  SystemConfiguration::SystemConfiguration() 
    : impl_(std::make_unique<Impl>()) {}
  
  SystemConfiguration::SystemConfiguration(const SystemConfiguration& other)
      : impl_(std::make_unique<Impl>(*other.impl_)) {}

  SystemConfiguration::SystemConfiguration(SystemConfiguration&& other) = default;

  SystemConfiguration& SystemConfiguration::operator=(const SystemConfiguration &other) {
    *impl_ = *other.impl_;
    return *this;
  }

  SystemConfiguration& SystemConfiguration::operator=(SystemConfiguration&&) = default;

  SystemConfiguration::~SystemConfiguration() = default;

  void SystemConfiguration::set_precision(ePrecision precision) {
    impl_->data.precision = precision;
  }

  ePrecision SystemConfiguration::get_precision() const {
    return impl_->data.precision;
  }
  
  void SystemConfiguration::set_device(eDevice device) {
    impl_->data.device = device;
  }

  eDevice SystemConfiguration::get_device() const {
    return impl_->data.device;
  }

  void SystemConfiguration::set_cpu_ncores(int cpu_ncores) {
    impl_->data.cpu_ncores = cpu_ncores;
  }

  int SystemConfiguration::get_cpu_ncores() const {
    return impl_->data.cpu_ncores;
  }
  
  void SystemConfiguration::set_cpu_nthread(int cpu_nthread) {
    impl_->data.cpu_nthread = cpu_nthread;
  }

  int SystemConfiguration::get_cpu_nthread() const {
    return impl_->data.cpu_nthread;
  }
  
  void SystemConfiguration::set_gpu_device(int gpu_device) {
    impl_->data.gpu_device = gpu_device;
  }

  int SystemConfiguration::get_gpu_device() const {
    return impl_->data.gpu_device;
  }
  
  void SystemConfiguration::set_gpu_nstream(int gpu_nstream) {
    impl_->data.gpu_nstream = gpu_nstream;
  }

  int SystemConfiguration::get_gpu_nstream() const {
    return impl_->data.gpu_nstream;
  }
  
  void SystemConfiguration::set_nstream(int nstream) {
    impl_->data.nstream = nstream;
  }

  int SystemConfiguration::get_nstream() const {
    return impl_->data.nstream;
  }
  
  void SystemConfiguration::set_active(bool active) {
    impl_->data.active = active;
  }

  bool SystemConfiguration::get_active() const {
    return impl_->data.active;
  }

  bool SystemConfiguration::is_host() const {
    return impl_->data.is_host();
  }

  bool SystemConfiguration::is_device() const {
    return impl_->data.is_device();
  }

  bool SystemConfiguration::is_float() const {
    return impl_->data.is_float();
  }

  bool SystemConfiguration::is_double() const {
    return impl_->data.is_double();
  }

  bool SystemConfiguration::is_float_host() const {
    return impl_->data.is_float_host();
  }

  bool SystemConfiguration::is_double_host() const {
    return impl_->data.is_double_host();
  }

  bool SystemConfiguration::is_float_device() const {
    return impl_->data.is_float_device();
  }

  bool SystemConfiguration::is_double_device() const {
    return impl_->data.is_double_device();
  }

  void test(const SystemConfiguration &a) {
    std::cout << "Hola Mundo" << std::endl;
  }

}
