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

#ifndef MULTEM_CONFIG_H
#define MULTEM_CONFIG_H

// The BUILDING_LIBRARY macro is set by cmake when building the shared library.
// In this case we want to set the appropriate DEVICE_CALLABLE and FORCE_INLINE
// flags. However, we don't want to FORCE_INLINE when compiling against the
// library because some functions will have their implementation visible
#ifdef BUILDING_LIBRARY
  #ifndef DEVICE_CALLABLE
  	#ifdef __CUDACC__
  		#define DEVICE_CALLABLE __host__ __device__
  		#define FORCE_INLINE __forceinline__
  	#else
  		#define DEVICE_CALLABLE
  		#define FORCE_INLINE inline
  	#endif
  #endif
#else
  #define DEVICE_CALLABLE
  #define FORCE_INLINE
#endif

#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_LIBRARY
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

#endif
