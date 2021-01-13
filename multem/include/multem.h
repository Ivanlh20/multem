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
#include <multem/constants.h>
#include <lin_alg_def.cuh>
#include <safe_types.cuh>
#include <atom_data_api.h>
#include <input_multislice_api.h>
#include <output_multislice_api.h>

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

namespace mt {
  
  /****************************************************************************
   * The Scanning interface
   ***************************************************************************/
  
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
  
  /****************************************************************************
   * The Detector interface
   ***************************************************************************/

  /**
   * Interface for Det_Int. This can be removed simply by using std::vector in
   * the Det_Int class rather than thrust::vector
   */
  template <typename T>
  class DetInt {
  public:
		using value_type = typename T::value_type;
		using size_type = std::size_t;

		size_type size() const {
			return image.size();
		}

    std::vector<T> image;
  };

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
