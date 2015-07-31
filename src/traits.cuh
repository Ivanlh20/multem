/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
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
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TRAITS_H
#define TRAITS_H

#include <type_traits>
#include "math.cuh"
#include "types.hpp"

namespace multem
{
	template <class TVector>
	using Value_type = typename TVector::value_type;

	template <class TVector>
	using Size_type = typename TVector::size_type;

	template<class T>
	struct is_bool: std::integral_constant<bool, std::is_same<T, bool>::value> {};

	template<class T>
	struct is_int: std::integral_constant<bool, std::is_same<T, int>::value || std::is_same<T, unsigned int>::value> {};

	template<class T>
	struct is_float: std::integral_constant<bool, std::is_same<T, float>::value> {};

	template<class T>
	struct is_double: std::integral_constant<bool, std::is_same<T, double>::value> {};

	template<class T>
	struct is_fundamental: std::integral_constant<bool, std::is_fundamental<T>::value || std::is_same<T, complex<float>>::value || std::is_same<T, complex<double>>::value> {};

	template<class T>
	struct is_rmatrix_r: std::integral_constant<bool, std::is_same<T, rmatrix_r>::value> {};

	template<class T>
	struct is_rmatrix_c: std::integral_constant<bool, std::is_same<T, rmatrix_c>::value> {};

	template<class T>
	struct is_enum_bool: std::integral_constant<bool, std::is_enum<T>::value || is_bool<T>::value> {};

	namespace detail_traits
	{
		template<class T>
		struct check_type 
		{
			typedef void type;
		};

		template<class T, class U = void>
		struct has_value_type
		{
			static const bool value = false;
		};

		template<class T>
		struct has_value_type<T, typename check_type<typename T::value_type>::type>
		{
			static const bool value = true;
		};

		template<class T, class Enable = void>
		struct has_device_member: std::integral_constant<bool, false> {};

		template<class T>
		struct has_device_member<T, typename std::enable_if<std::is_class<T>::value>::type>
		{
			struct Fallback { int device; };
			struct Derived: T, Fallback {};

			template<typename C, C> struct ChT; 

			template<typename C> static char (&f(ChT<int Fallback::*, &C::device>*))[1]; 
			template<typename C> static char (&f(...))[2]; 

			static const bool value = sizeof(f<Derived>(0)) == 2;
		}; 
	}

	template<class T, class Enable = void>
	struct is_host_vector: std::integral_constant<bool, false> {};

	template<class T>
	struct is_host_vector<T, typename std::enable_if<detail_traits::has_value_type<T>::value>::type>: std::integral_constant<bool, std::is_same<T, host_vector<Value_type<T>>>::value || std::is_same<T, vector<Value_type<T>>>::value> {};

	template<class T, class Enable = void>
	struct is_device_vector: std::integral_constant<bool, false> {};

	template<class T>
	struct is_device_vector<T, typename std::enable_if<detail_traits::has_value_type<T>::value>::type>: std::integral_constant<bool, std::is_same<T, device_vector<Value_type<T>>>::value> {};

	template<class T, eDevice dev=e_none>
	struct is_Vector: std::integral_constant<bool, (dev == e_none)?(is_host_vector<T>::value || is_device_vector<T>::value):((dev == e_Host)?is_host_vector<T>::value:is_device_vector<T>::value)> {};

	template<class T, eDevice dev, class Enable = void>
	struct is_Host_Device: std::integral_constant<bool, is_Vector<T, dev>::value> {};

	template<class T, eDevice dev>
	struct is_Host_Device<T, dev, typename std::enable_if<detail_traits::has_device_member<T>::value>::type>: std::integral_constant<bool, T::device == dev> {};

	template<class T>
	struct is_Host: std::integral_constant<bool, is_Host_Device<T, e_Host>::value> {};

	template<class T>
	struct is_Device: std::integral_constant<bool, is_Host_Device<T, e_Device>::value> {};

	template <class T, class U>
	using enable_if_host_vector = typename std::enable_if<is_host_vector<T>::value, U>::type;

	template <class T, class U>
	using enable_if_device_vector = typename std::enable_if<is_device_vector<T>::value, U>::type;

	template <class T, class U>
	using enable_if_Host = typename std::enable_if<is_Host<T>::value, U>::type;

	template <class T, class U>
	using enable_if_Device = typename std::enable_if<is_Device<T>::value, U>::type;

	template <class T, class U>
	using enable_if_float = typename std::enable_if<is_float<T>::value, U>::type;

	template <class T, class U>
	using enable_if_double = typename std::enable_if<is_double<T>::value, U>::type;

	template <class T, class U>
	using enable_if_int = typename std::enable_if<is_int<T>::value, U>::type;

	template <class T, class U>
	using enable_if_bool = typename std::enable_if<is_bool<T>::value, U>::type;

	template <class T, class U>
	using enable_if_floating_point = typename std::enable_if<std::is_floating_point<T>::value, U>::type;
	
	template <class T, class U>
	using enable_if_enum_bool = typename std::enable_if<is_enum_bool<T>::value, U>::type;
		
	template <class T, class U>
	using enable_if_pointer = typename std::enable_if<std::is_pointer<T>::value, U>::type;

	template <class T, class U>
	using enable_if_rmatrix_r = typename std::enable_if<is_rmatrix_r<T>::value, U>::type;

	template <class T, class U>
	using enable_if_rmatrix_c = typename std::enable_if<is_rmatrix_c<T>::value, U>::type;

	template <class T, class U>
	using enable_if_rmatrix = typename std::enable_if<is_rmatrix_r<T>::value || is_rmatrix_c<T>::value, U>::type;	

	template <int simulation_type, class U>
	using enable_if_STEM = typename std::enable_if<simulation_type == eST_STEM, U>::type;

	template <int simulation_type, class U>
	using enable_if_ISTEM = typename std::enable_if<simulation_type == eST_ISTEM, U>::type;

	template <int simulation_type, class U>
	using enable_if_CBED = typename std::enable_if<simulation_type == eST_CBED, U>::type;

	template <int simulation_type, class U>
	using enable_if_CBEI = typename std::enable_if<simulation_type == eST_CBEI, U>::type;

	template <int simulation_type, class U>
	using enable_if_ED = typename std::enable_if<simulation_type == eST_ED, U>::type;

	template <int simulation_type, class U>
	using enable_if_HRTEM = typename std::enable_if<simulation_type == eST_HRTEM, U>::type;

	template <int simulation_type, class U>
	using enable_if_PED = typename std::enable_if<simulation_type == eST_PED, U>::type;

	template <int simulation_type, class U>
	using enable_if_HCI = typename std::enable_if<simulation_type == eST_HCI, U>::type;

	template <int simulation_type, class U>
	using enable_if_EWFS = typename std::enable_if<simulation_type == eST_EWFS, U>::type;

	template <int simulation_type, class U>
	using enable_if_EWRS = typename std::enable_if<simulation_type == eST_EWRS, U>::type;

	template <int simulation_type, class U>
	using enable_if_EELS = typename std::enable_if<simulation_type == eST_EELS, U>::type;

	template <int simulation_type, class U>
	using enable_if_EFTEM = typename std::enable_if<simulation_type == eST_EFTEM, U>::type;

	template <int simulation_type, class U>
	using enable_if_ProbeFS = typename std::enable_if<simulation_type == eST_ProbeFS, U>::type;

	template <int simulation_type, class U>
	using enable_if_ProbeRS = typename std::enable_if<simulation_type == eST_ProbeRS, U>::type;
} // namespace multem

#endif