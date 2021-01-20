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

#include <array>
#include <cassert>
#include <complex>
#include <string>
#include <vector>
#include <memory>

#include <multem/amorp_spec.h>
#include <multem/atom_cal.h>
#include <multem/atom_data.h>
#include <multem/atomic_data.h>
#include <multem/box_occ.h>
#include <multem/cgpu_rand.h>
#include <multem/config.h>
#include <multem/constants.h>
#include <multem/error.h>
#include <multem/fxeg_data.h>
#include <multem/input_multislice.h>
#include <multem/lin_alg_def.h>
#include <multem/math.h>
#include <multem/multem.h>
#include <multem/output_multislice.h>
#include <multem/stream.h>
#include <multem/types.h>
#include <multem/xtl_build.h>

namespace mt {
  
	template<class T>
	DLL_PUBLIC void rdf_3d(const Atom_Data<T> &atoms, T r_max, int nr, std::vector<T> &r, std::vector<T> &rdf);

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
  DLL_PUBLIC mt::Output_Multislice<T> incident_wave(Input_Multislice<T>&);
  
  template <typename T>
  DLL_PUBLIC mt::Output_Multislice<T> microscope_aberrations(Input_Multislice<T>&);
  
  template <typename T>
  DLL_PUBLIC mt::Output_Multislice<T> projected_potential(Input_Multislice<T>&);
  
  template <typename T>
  DLL_PUBLIC mt::Output_Multislice<T> propagate(Input_Multislice<T>&);
  
  template <typename T>
  DLL_PUBLIC std::vector<T> spec_planes(Input_Multislice<T>&);
  
  template <typename T>
  DLL_PUBLIC std::tuple<Atom_Data<T>, std::vector<Slice<T>>> spec_slicing(Input_Multislice<T>&);
  
  template <typename T>
  DLL_PUBLIC mt::Output_Multislice<T> tem_simulation(Input_Multislice<T>&);
  
  template <typename T>
  DLL_PUBLIC mt::Output_Multislice<T> transmission_function(Input_Multislice<T>&);
  
  template <typename T>
  DLL_PUBLIC mt::Output_Multislice<T> wave_function(Input_Multislice<T>&);
  

}  // namespace mt

#endif
