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

#ifndef MULTEM_H
#define MULTEM_H

#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_DLL
    #ifdef __GNUC__
      #define DLL_PUBLIC __attribute__ ((dllexport))
    #else
      #define DLL_PUBLIC __declspec(dllexport)
    #endif
  #else
    #ifdef __GNUC__
      #define DLL_PUBLIC __attribute__ ((dllimport))
    #else
      #define DLL_PUBLIC __declspec(dllimport)
    #endif
  #endif
#else
  #if __GNUC__ >= 4
    #define DLL_PUBLIC __attribute__ ((visibility ("default")))
  #else
    #define DLL_PUBLIC
  #endif
#endif

#include <array>
#include <cassert>
#include <complex>
#include <string>
#include <vector>
#include <memory>
#include <multem/constants.h>
#include <lin_alg_def.cuh>
#include <safe_types.cuh>
#include <atom_data_api.h>
#include <input_multislice_api.h>
#include <output_multislice_api.h>
#include <amorp_spec_api.h>
#include <xtl_build.hpp>
#include <atom_cal_api.h>
#include <fxeg_data_api.h>

namespace mt {
  

  /****************************************************************************
   * Logging functions
   ***************************************************************************/

  DLL_PUBLIC std::string to_string(const System_Configuration &sys_conf, const std::string &prefix="");

  template <typename T>
  DLL_PUBLIC std::string to_string(const Lens<T> &lens, const std::string &prefix="");

  template <typename T>
  DLL_PUBLIC std::string to_string(const Input_Multislice<T> &input, const std::string &prefix="");

  /****************************************************************************
   * Simulation functions
   ***************************************************************************/

  template <typename T>
  DLL_PUBLIC mt::Output_Multislice<T> tem_simulation(Input_Multislice<T>&);

}  // namespace mt

#endif
