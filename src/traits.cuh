/*
 * This file is part of MULTEM.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef TRAITS_H
#define TRAITS_H

#include <type_traits>
#include "math.cuh"
#include "types.cuh"

namespace mt
{
	template <class T>
	struct is_bool: std::integral_constant<bool, std::is_same<T, bool>::value> {};

	template <class T>
	struct is_int: std::integral_constant<bool, std::is_same<T, int>::value || std::is_same<T, unsigned int>::value> {};

	template <class T>
	struct is_float: std::integral_constant<bool, std::is_same<T, float>::value> {};

	template <class T>
	struct is_double: std::integral_constant<bool, std::is_same<T, double>::value> {};

	template <class T>
	struct is_real: std::integral_constant<bool, is_float<T>::value || is_double<T>::value> {};

	template <class T>
	struct is_cfloat: std::integral_constant<bool, std::is_same<T, complex<float>>::value> {};

	template <class T>
	struct is_cdouble: std::integral_constant<bool, std::is_same<T, complex<double>>::value> {};

	template <class T>
	struct is_complex: std::integral_constant<bool, is_cfloat<T>::value || is_cdouble<T>::value> {};

	template <class T>
	struct is_fundamental: std::integral_constant<bool, std::is_fundamental<T>::value || is_complex<T>::value> {};

	template <class T>
	struct is_enum_bool: std::integral_constant<bool, std::is_enum<T>::value || is_bool<T>::value> {};

	namespace detail_traits
	{
		template <class T>
		struct check_type 
		{
			typedef void type;
		};

		template <class T, class U = void>
		struct has_value_type
		{
			static const bool value = false;
		};

		template <class T>
		struct has_value_type<T, typename check_type<typename T::value_type>::type>
		{
			static const bool value = true;
		};

		template <class TVector, class Enable = void>
		struct get_value_type_r
		{
			using type = typename TVector::value_type;
		};

		template <class TVector>
		struct get_value_type_r<TVector, typename std::enable_if<is_complex<typename TVector::value_type>::value>::type>
		{
			using type = typename TVector::value_type::value_type;
		};

		template <class T, class Enable = void>
		struct has_device_member: std::integral_constant<bool, false> {};

		template <class T>
		struct has_device_member<T, typename std::enable_if<std::is_class<T>::value>::type>
		{
			struct Fallback { int device; };
			struct Derived: T, Fallback {};

			template <typename C, C> struct ChT; 

			template <typename C> static char (&f(ChT<int Fallback::*, &C::device>*))[1]; 
			template <typename C> static char (&f(...))[2]; 

			static const bool value = sizeof(f<Derived>(0)) == 2;
		}; 
	}

	/**************************Output data type**************************/
	template <eData_Type type>
	using Data_Type = typename std::conditional<type == eDT_float, float, typename std::conditional<type == eDT_double, double, \
	typename std::conditional<type == eDT_cfloat, complex<float>, complex<double>>::type>::type >::type;


	template <class TVector>
	using Value_type = typename TVector::value_type;

	template <class TVector>
	using Value_type_r = typename detail_traits::get_value_type_r<TVector>::type;

	template <class TVector>
	using Size_type = typename TVector::size_type;

	template <class T, class Enable = void>
	struct is_host_vector: std::integral_constant<bool, false> {};

	template <class T>
	struct is_host_vector<T, typename std::enable_if<detail_traits::has_value_type<T>::value>::type>: std::integral_constant<bool, std::is_same<T, host_vector<Value_type<T>>>::value || std::is_same<T, vector<Value_type<T>>>::value> {};

	template <class T, class Enable = void>
	struct is_device_vector: std::integral_constant<bool, false> {};

	template <class T>
	struct is_device_vector<T, typename std::enable_if<detail_traits::has_value_type<T>::value>::type>: std::integral_constant<bool, std::is_same<T, device_vector<Value_type<T>>>::value> {};

	template <class T>
	struct is_complex_host_vector: std::integral_constant<bool, is_host_vector<T>::value && is_complex<typename T::value_type>::value> {};

	template <class T>
	struct is_complex_device_vector: std::integral_constant<bool, is_device_vector<T>::value && is_complex<typename T::value_type>::value> {};

	template <class T>
	struct is_host_device_vector: std::integral_constant<bool, is_host_vector<T>::value || is_device_vector<T>::value> {};

	template <class T, eDevice dev =e_host_device>
	struct is_Vector: std::integral_constant<bool, (dev == e_host_device)?(is_host_device_vector<T>::value):((dev == e_host)?is_host_vector<T>::value:is_device_vector<T>::value)> {};

	template <class T, eDevice dev, class Enable = void>
	struct is_host_device: std::integral_constant<bool, is_Vector<T, dev>::value> {};

	template <class T, eDevice dev>
	struct is_host_device<T, dev, typename std::enable_if<detail_traits::has_device_member<T>::value>::type>: std::integral_constant<bool, T::device == dev> {};

	template <class T>
	struct is_host: std::integral_constant<bool, is_host_device<T, e_host>::value> {};

	template <class T>
	struct is_device: std::integral_constant<bool, is_host_device<T, e_device>::value> {};

	template <eDevice dev>
	struct is_dev_host: std::integral_constant<bool, dev == e_host> {};

	template <eDevice dev>
	struct is_dev_device: std::integral_constant<bool, dev == e_device> {};

	template <class T1, class T2>
	struct is_host_vector_and_host_vector: std::integral_constant<bool, is_host_vector<T1>::value && is_host_vector<T2>::value> {};

	template <class T1, class T2>
	struct is_host_vector_and_device_vector: std::integral_constant<bool, is_host_vector<T1>::value && is_device_vector<T2>::value> {};

	template <class T1, class T2>
	struct is_device_vector_and_host_vector: std::integral_constant<bool, is_device_vector<T1>::value && is_host_vector<T2>::value> {};

	template <class T1, class T2>
	struct is_device_vector_and_device_vector: std::integral_constant<bool, is_device_vector<T1>::value && is_device_vector<T2>::value> {};

	template <class T1, class T2>
	struct is_complex_host_vector_and_host_vector: std::integral_constant<bool, is_complex_host_vector<T1>::value && is_host_vector<T2>::value> {};

	template <class T1, class T2>
	struct is_complex_device_vector_and_device_vector: std::integral_constant<bool, is_complex_device_vector<T1>::value && is_device_vector<T2>::value> {};

	template <class T1, class T2>
	struct is_complex_vector_and_real_vector: std::integral_constant<bool, is_complex<typename T1::value_type>::value && is_real<typename T2::value_type>::value> {};

	template <class T, class U>
	using enable_if_real_host_vector = typename std::enable_if<is_host_vector<T>::value && is_real<typename T::value_type>::value, U>::type;

	template <class T, class U>
	using enable_if_complex_host_vector = typename std::enable_if<is_complex_host_vector<T>::value, U>::type;

	template <class T, class U>
	using enable_if_host_vector = typename std::enable_if<is_host_vector<T>::value, U>::type;

	template <class T, class U>
	using enable_if_real_device_vector = typename std::enable_if<is_device_vector<T>::value && is_real<typename T::value_type>::value, U>::type;

	template <class T, class U>
	using enable_if_complex_device_vector = typename std::enable_if<is_complex_device_vector<T>::value, U>::type;

	template <class T, class U>
	using enable_if_device_vector = typename std::enable_if<is_device_vector<T>::value, U>::type;

	template <class T1, class T2, class U>
	using enable_if_host_vector_and_device_vector = typename std::enable_if<is_host_vector_and_device_vector<T1, T2>::value, U>::type;

	template <class T1, class T2, class U>
	using enable_if_device_vector_and_host_vector = typename std::enable_if<is_device_vector_and_host_vector<T1, T2>::value, U>::type;

	template <class T1, class T2, class U>
	using enable_if_host_vector_and_host_vector = typename std::enable_if<is_host_vector_and_host_vector<T1, T2>::value, U>::type;

	template <class T1, class T2, class U>
	using enable_if_device_vector_and_device_vector = typename std::enable_if<is_device_vector_and_device_vector<T1, T2>::value, U>::type;

	template <class T1, class T2, class U>
	using enable_if_complex_host_vector_and_host_vector = typename std::enable_if<is_complex_host_vector_and_host_vector<T1, T2>::value, U>::type;

	template <class T1, class T2, class U>
	using enable_if_complex_device_vector_and_device_vector = typename std::enable_if<is_complex_device_vector_and_device_vector<T1, T2>::value, U>::type;

	template <class T1, class T2, class U>
	using enable_if_complex_vector_and_real_vector = typename std::enable_if<is_complex_vector_and_real_vector<T1, T2>::value, U>::type;

	template <class T, class U>
	using enable_if_host = typename std::enable_if<is_host<T>::value, U>::type;

	template <class T, class U>
	using enable_if_device = typename std::enable_if<is_device<T>::value, U>::type;

	template <eDevice dev, class U>
	using enable_if_dev_host = typename std::enable_if<is_dev_host<dev>::value, U>::type;

	template <eDevice dev, class U>
	using enable_if_dev_device = typename std::enable_if<is_dev_device<dev>::value, U>::type;

	template <class T, class U>
	using enable_if_float = typename std::enable_if<is_float<T>::value, U>::type;

	template <class T, class U>
	using enable_if_double = typename std::enable_if<is_double<T>::value, U>::type;

	template <class T, class U>
 using enable_if_real = typename std::enable_if<is_real<T>::value, U>::type;

	template <class T, class U>
	using enable_if_cfloat = typename std::enable_if<is_cfloat<T>::value, U>::type;

	template <class T, class U>
	using enable_if_cdouble = typename std::enable_if<is_cdouble<T>::value, U>::type;

	template <class T, class U>
	using enable_if_complex = typename std::enable_if<is_complex<T>::value, U>::type;

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

	template <int simulation_type, class U>
	using enable_if_STEM = typename std::enable_if<simulation_type == eTEMST_STEM, U>::type;

	template <int simulation_type, class U>
	using enable_if_ISTEM = typename std::enable_if<simulation_type == eTEMST_ISTEM, U>::type;

	template <int simulation_type, class U>
	using enable_if_CBED = typename std::enable_if<simulation_type == eTEMST_CBED, U>::type;

	template <int simulation_type, class U>
	using enable_if_CBEI = typename std::enable_if<simulation_type == eTEMST_CBEI, U>::type;

	template <int simulation_type, class U>
	using enable_if_ED = typename std::enable_if<simulation_type == eTEMST_ED, U>::type;

	template <int simulation_type, class U>
	using enable_if_HRTEM = typename std::enable_if<simulation_type == eTEMST_HRTEM, U>::type;

	template <int simulation_type, class U>
	using enable_if_PED = typename std::enable_if<simulation_type == eTEMST_PED, U>::type;

	template <int simulation_type, class U>
	using enable_if_HCTEM = typename std::enable_if<simulation_type == eTEMST_HCTEM, U>::type;

	template <int simulation_type, class U>
	using enable_if_EWFS = typename std::enable_if<simulation_type == eTEMST_EWFS, U>::type;

	template <int simulation_type, class U>
	using enable_if_EWRS = typename std::enable_if<simulation_type == eTEMST_EWRS, U>::type;

	template <int simulation_type, class U>
	using enable_if_EELS = typename std::enable_if<simulation_type == eTEMST_EELS, U>::type;

	template <int simulation_type, class U>
	using enable_if_EFTEM = typename std::enable_if<simulation_type == eTEMST_EFTEM, U>::type;

	template <int simulation_type, class U>
	using enable_if_ProbeFS = typename std::enable_if<simulation_type == eTEMST_IWFS, U>::type;

	template <int simulation_type, class U>
	using enable_if_ProbeRS = typename std::enable_if<simulation_type == eTEMST_IWRS, U>::type;


	/*********************************type***********************************/
 template<class TVector>
 enable_if_float<Value_type<TVector>, int>
	matrix_type(TVector &vector){ return 1; }

 template<class TVector>
 enable_if_double<Value_type<TVector>, int>
	matrix_type(TVector &vector){ return 2; }

 template<class TVector>
 enable_if_cfloat<Value_type<TVector>, int>
	matrix_type(TVector &vector){ return 3; }

 template<class TVector>
 enable_if_cdouble<Value_type<TVector>, int>
	matrix_type(TVector &vector){ return 4; }

} // namespace mt

#endif
