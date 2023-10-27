/*
* This file is part of Multem.
* Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
*
* Multem is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version of the License, or
* (at your option) any later version.
*
* Multem is distributed in the hope that it will be useful, 
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Multem. If not, see <http:// www.gnu.org/licenses/>.
*/

#pragma once

#include <type_traits>

#include "macros.h"
#include "const_enum.h"

#ifdef __CUDACC__
	#include <cuda.h>
	#include <math.h>
	#include <thrust/complex.h>
#endif

/* same decay */
namespace mt
{
	template <class T, class U>
	struct is_same_decay: std::integral_constant<bool, std::is_same<typename std::decay<T>::type, typename std::decay<U>::type>::value> {};		
		
	template <class T, class U>
	struct is_diff_decay: std::integral_constant<bool, !is_same_decay<T, U>::value> {};
}	
	
/* enable if */
namespace mt
{
	template <class T, class U, class V=void>
	using enable_if_same_decay = typename std::enable_if<is_same_decay<T, U>::value, V>::type;

	template <class T, class U, class V=void>
	using enable_if_diff_decay = typename std::enable_if<is_diff_decay<T, U>::value, V>::type;
}

/* enum types - verification */
namespace mt
{
	dt_bool is_ebool(const eData_Typ& ed_typ)
	{ 
		return edt_bool == ed_typ;
	}

	dt_bool is_eint8(const eData_Typ& ed_typ)
	{ 
		return edt_int8 == ed_typ;
	}

	dt_bool is_euint8(const eData_Typ& ed_typ)
	{ 
		return edt_uint8 == ed_typ;
	}

	dt_bool is_eint16(const eData_Typ& ed_typ)
	{ 
		return edt_int16 == ed_typ;
	}

	dt_bool is_euint16(const eData_Typ& ed_typ)
	{ 
		return edt_uint16 == ed_typ;
	}

	dt_bool is_eint32(const eData_Typ& ed_typ)
	{ 
		return edt_int32 == ed_typ;
	}

	dt_bool is_euint32(const eData_Typ& ed_typ)
	{ 
		return edt_uint32 == ed_typ;
	}

	dt_bool is_eint64(const eData_Typ& ed_typ)
	{ 
		return edt_int64 == ed_typ;
	}

	dt_bool is_euint64(const eData_Typ& ed_typ)
	{ 
		return edt_uint64 == ed_typ;
	}

	dt_bool is_efloat32(const eData_Typ& ed_typ)
	{ 
		return edt_float32 == ed_typ;
	}

	dt_bool is_efloat64(const eData_Typ& ed_typ)
	{ 
		return edt_float64 == ed_typ;
	}

	dt_bool is_ecint8(const eData_Typ& ed_typ)
	{ 
		return edt_cint8 == ed_typ;
	}

	dt_bool is_ecuint8(const eData_Typ& ed_typ)
	{ 
		return edt_cuint8 == ed_typ;
	}

	dt_bool is_ecint16(const eData_Typ& ed_typ)
	{ 
		return edt_cint16 == ed_typ;
	}

	dt_bool is_ecuint16(const eData_Typ& ed_typ)
	{ 
		return edt_cuint16 == ed_typ;
	}

	dt_bool is_ecint32(const eData_Typ& ed_typ)
	{ 
		return edt_cint32 == ed_typ;
	}

	dt_bool is_ecuint32(const eData_Typ& ed_typ)
	{ 
		return edt_cuint32 == ed_typ;
	}

	dt_bool is_ecint64(const eData_Typ& ed_typ)
	{ 
		return edt_cint64 == ed_typ;
	}

	dt_bool is_ecuint64(const eData_Typ& ed_typ)
	{ 
		return edt_cuint64 == ed_typ;
	}

	dt_bool is_ecfloat32(const eData_Typ& ed_typ)
	{ 
		return edt_cfloat32 == ed_typ;
	}

	dt_bool is_ecfloat64(const eData_Typ& ed_typ)
	{ 
		return edt_cfloat64 == ed_typ;
	}
}

/* real types - verification */
namespace mt
{
	template <class T>
	struct is_bool : std::integral_constant<dt_bool, is_same_decay<T, dt_bool>::value> {};		

	template <class T>
	struct is_enum : std::integral_constant<dt_bool, std::is_enum<typename std::decay<T>::type>::value> {};

	template <class T>
	struct is_enum_bool: std::integral_constant<dt_bool, is_enum<T>::value || is_bool<T>::value> {};

	template <class T>
	struct is_int8 : std::integral_constant<dt_bool, is_same_decay<T, dt_int8>::value> {};

	template <class T>
	struct is_uint8 : std::integral_constant<dt_bool, is_same_decay<T, dt_uint8>::value> {};

	template <class T>
	struct is_int16 : std::integral_constant<dt_bool, is_same_decay<T, dt_int16>::value> {};

	template <class T>
	struct is_uint16 : std::integral_constant<dt_bool, is_same_decay<T, dt_uint16>::value> {};

	template <class T>
	struct is_int32 : std::integral_constant<dt_bool, is_same_decay<T, dt_int32>::value> {};

	template <class T>
	struct is_uint32 : std::integral_constant<dt_bool, is_same_decay<T, dt_uint32>::value> {};

	template <class T>
	struct is_int64 : std::integral_constant<dt_bool, is_same_decay<T, dt_int64>::value> {};

	template <class T>
	struct is_uint64 : std::integral_constant<dt_bool, is_same_decay<T, dt_uint64>::value> {};

	template <class T>
	struct is_float32 : std::integral_constant<dt_bool, is_same_decay<T, dt_float32>::value> {};

	template <class T>
	struct is_float64 : std::integral_constant<dt_bool, is_same_decay<T, dt_float64>::value> {};

	template <typename T> 
	struct is_int: std::integral_constant<dt_bool, is_int8<T>::value || is_uint8<T>::value \
	|| is_int16<T>::value || is_uint16<T>::value || is_int32<T>::value || is_uint32<T>::value \
	|| is_int64<T>::value || is_uint64<T>::value> {};
		
	template <class T>
	struct is_float : std::integral_constant<dt_bool, is_float32<T>::value || is_float64<T>::value> {};

	template <class T>
	struct is_real : std::integral_constant<dt_bool, is_bool<T>::value || is_int8<T>::value \
		|| is_uint8<T>::value || is_int16<T>::value || is_uint16<T>::value || is_int32<T>::value \
		|| is_uint32<T>::value || is_int64<T>::value || is_uint64<T>::value || is_float<T>::value> {};
}

/* std complex types - verification */
namespace mt
{
	template <class T, class U = typename std::decay<T>::type>
	struct is_std_cint8 : std::integral_constant<dt_bool, std::is_same<U, dt_std_cint8>::value> {};

	template <class T, class U = typename std::decay<T>::type>
	struct is_std_cuint8 : std::integral_constant<dt_bool, std::is_same<U, dt_std_cuint8>::value> {};

	template <class T, class U = typename std::decay<T>::type>
	struct is_std_cint16 : std::integral_constant<dt_bool, std::is_same<U, dt_std_cint16>::value> {};

	template <class T, class U = typename std::decay<T>::type>
	struct is_std_cuint16 : std::integral_constant<dt_bool, std::is_same<U, dt_std_cuint16>::value> {};

	template <class T, class U = typename std::decay<T>::type>
	struct is_std_cint32 : std::integral_constant<dt_bool, std::is_same<U, dt_std_cint32>::value> {};

	template <class T, class U = typename std::decay<T>::type>
	struct is_std_cuint32 : std::integral_constant<dt_bool, std::is_same<U, dt_std_cuint32>::value> {};

	template <class T, class U = typename std::decay<T>::type>
	struct is_std_cint64 : std::integral_constant<dt_bool, std::is_same<U, dt_std_cint64>::value> {};

	template <class T, class U = typename std::decay<T>::type>
	struct is_std_cuint64 : std::integral_constant<dt_bool, std::is_same<U, dt_std_cuint64>::value> {};

	template <class T, class U = typename std::decay<T>::type>
	struct is_std_cfloat32 : std::integral_constant<dt_bool, std::is_same<U, dt_std_cfloat32>::value> {};

	template <class T, class U = typename std::decay<T>::type>
	struct is_std_cfloat64 : std::integral_constant<dt_bool, std::is_same<U, dt_std_cfloat64>::value> {};

	template <class T>
	struct is_std_cfloat : std::integral_constant<dt_bool, is_std_cfloat32<T>::value || is_std_cfloat64<T>::value> {};

	template <class T>
	struct is_std_cmplx : std::integral_constant<dt_bool, is_std_cint8<T>::value \
		|| is_std_cuint8<T>::value || is_std_cint16<T>::value || is_std_cuint16<T>::value || is_std_cint32<T>::value \
		|| is_std_cuint32<T>::value || is_std_cint64<T>::value || is_std_cuint64<T>::value || is_std_cfloat<T>::value> {};
}

/* thrust types - verification */
namespace mt
{
#ifdef __CUDACC__
	template <class T>
	struct is_thr_cint8: std::integral_constant<dt_bool, is_same_decay<T, dt_thr_cint8>::value> {};

	template <class T>
	struct is_thr_cuint8: std::integral_constant<dt_bool, is_same_decay<T, dt_thr_cuint8>::value> {};

	template <class T>
	struct is_thr_cint16: std::integral_constant<dt_bool, is_same_decay<T, dt_thr_cint16>::value> {};

	template <class T>
	struct is_thr_cuint16: std::integral_constant<dt_bool, is_same_decay<T, dt_thr_cuint16>::value> {};

	template <class T>
	struct is_thr_cint32: std::integral_constant<dt_bool, is_same_decay<T, dt_thr_cint32>::value> {};

	template <class T>
	struct is_thr_cuint32: std::integral_constant<dt_bool, is_same_decay<T, dt_thr_cuint32>::value> {};

	template <class T>
	struct is_thr_cint64: std::integral_constant<dt_bool, is_same_decay<T, dt_thr_cint64>::value> {};

	template <class T>
	struct is_thr_cuint64: std::integral_constant<dt_bool, is_same_decay<T, dt_thr_cuint64>::value> {};

	template <class T>
	struct is_thr_cfloat32: std::integral_constant<dt_bool, is_same_decay<T, dt_thr_cfloat32>::value> {};

	template <class T>
	struct is_thr_cfloat64: std::integral_constant<dt_bool, is_same_decay<T, dt_thr_cfloat64>::value> {};

	template <class T>
	struct is_thr_cfloat: std::integral_constant<dt_bool, is_thr_cfloat32<T>::value || is_thr_cfloat64<T>::value> {};

	template <class T>
	struct is_thr_cmplx: std::integral_constant<dt_bool, is_thr_cint8<T>::value \
	|| is_thr_cuint8<T>::value || is_thr_cint16<T>::value || is_thr_cuint16<T>::value || is_thr_cint32<T>::value \
	|| is_thr_cuint32<T>::value || is_thr_cint64<T>::value || is_thr_cuint64<T>::value || is_thr_cfloat<T>::value> {};
#endif
}

/* complex types - verification */
namespace mt
{
#ifdef __CUDACC__
	template <class T>
	struct is_cint8: std::integral_constant<dt_bool, is_std_cint8<T>::value || is_thr_cint8<T>::value> {};

	template <class T>
	struct is_cuint8: std::integral_constant<dt_bool, is_std_cuint8<T>::value || is_thr_cuint8<T>::value> {};

	template <class T>
	struct is_cint16: std::integral_constant<dt_bool, is_std_cint16<T>::value || is_thr_cint16<T>::value> {};

	template <class T>
	struct is_cuint16: std::integral_constant<dt_bool, is_std_cuint16<T>::value || is_thr_cuint16<T>::value> {};

	template <class T>
	struct is_cint32: std::integral_constant<dt_bool, is_std_cint32<T>::value || is_thr_cint32<T>::value> {};

	template <class T>
	struct is_cuint32: std::integral_constant<dt_bool, is_std_cuint32<T>::value || is_thr_cuint32<T>::value> {};

	template <class T>
	struct is_cint64: std::integral_constant<dt_bool, is_std_cint64<T>::value || is_thr_cint64<T>::value> {};

	template <class T>
	struct is_cuint64: std::integral_constant<dt_bool, is_std_cuint64<T>::value || is_thr_cuint64<T>::value> {};

	template <class T>
	struct is_cfloat32: std::integral_constant<dt_bool, is_std_cfloat32<T>::value || is_thr_cfloat32<T>::value> {};

	template <class T>
	struct is_cfloat64: std::integral_constant<dt_bool, is_std_cfloat64<T>::value || is_thr_cfloat64<T>::value> {};

	template <class T>
	struct is_cmplx: std::integral_constant<dt_bool, is_std_cmplx<T>::value || is_thr_cmplx<T>::value> {};
#else
	template <class T>
	struct is_cint8: std::integral_constant<dt_bool, is_std_cint8<T>::value> {};

	template <class T>
	struct is_cuint8: std::integral_constant<dt_bool, is_std_cuint8<T>::value> {};

	template <class T>
	struct is_cint16: std::integral_constant<dt_bool, is_std_cint16<T>::value> {};

	template <class T>
	struct is_cuint16: std::integral_constant<dt_bool, is_std_cuint16<T>::value> {};

	template <class T>
	struct is_cint32: std::integral_constant<dt_bool, is_std_cint32<T>::value> {};

	template <class T>
	struct is_cuint32: std::integral_constant<dt_bool, is_std_cuint32<T>::value> {};

	template <class T>
	struct is_cint64: std::integral_constant<dt_bool, is_std_cint64<T>::value> {};

	template <class T>
	struct is_cuint64: std::integral_constant<dt_bool, is_std_cuint64<T>::value> {};

	template <class T>
	struct is_cfloat32: std::integral_constant<dt_bool, is_std_cfloat32<T>::value> {};

	template <class T>
	struct is_cfloat64: std::integral_constant<dt_bool, is_std_cfloat64<T>::value> {};

	template <class T>
	struct is_cmplx: std::integral_constant<dt_bool, is_std_cmplx<T>::value> {};
#endif

	template <class T>
	struct is_cfloat: std::integral_constant<dt_bool, is_cfloat32<T>::value || is_cfloat64<T>::value> {};

	template <class T>
	struct is_number: std::integral_constant<dt_bool, is_real<T>::value || is_cmplx<T>::value> {};
}

/* real types - enable */
namespace mt
{
	template <class T, class U=void>
	using enable_if_bool = typename std::enable_if<is_bool<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_enum = typename std::enable_if<is_enum<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_enum_bool = typename std::enable_if<is_enum_bool<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_int8 = typename std::enable_if<is_int8 <T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_uint8 = typename std::enable_if<is_uint8<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_int16 = typename std::enable_if<is_int16 <T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_uint16 = typename std::enable_if<is_uint16 <T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_int32 = typename std::enable_if<is_int32<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_uint32 = typename std::enable_if<is_uint32<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_int64 = typename std::enable_if<is_int64<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_uint64 = typename std::enable_if<is_uint64<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_float32 = typename std::enable_if<is_float32<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_float64 = typename std::enable_if<is_float64<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_int = typename std::enable_if<is_int<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_float = typename std::enable_if<is_float<T>::value, U>::type;

	template <class T, class U, class V=void>
	using enable_if_float_float = typename std::enable_if<is_float<T>::value && is_float<U>::value, V>::type;

	template <class T, class U=void>
	using enable_if_real = typename std::enable_if<is_real<T>::value, U>::type;

	template <class T, class U, class V=void>
	using enable_if_real_real = typename std::enable_if<is_real<T>::value && is_real<U>::value, V>::type;
}

/* std complex types - enable */
namespace mt
{
	template <class T, class U=void>
	using enable_if_std_cint8 = typename std::enable_if<is_std_cint8 <T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_std_cuint8 = typename std::enable_if<is_std_cuint8<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_std_cint16 = typename std::enable_if<is_std_cint16 <T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_std_cuint16 = typename std::enable_if<is_std_cuint16 <T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_std_cint32 = typename std::enable_if<is_std_cint32<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_std_cuint32 = typename std::enable_if<is_std_cuint32<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_std_cint64 = typename std::enable_if<is_std_cint64<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_std_cuint64 = typename std::enable_if<is_std_cuint64<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_std_cfloat32 = typename std::enable_if<is_std_cfloat32<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_std_cfloat64 = typename std::enable_if<is_std_cfloat64<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_std_cfloat = typename std::enable_if<is_std_cfloat<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_std_cmplx = typename std::enable_if<is_std_cmplx<T>::value, U>::type;
}

/* thrust complex types - enable */
namespace mt
{
#ifdef __CUDACC__
	template <class T, class U=void>
	using enable_if_thr_cint8 = typename std::enable_if<is_thr_cint8 <T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_thr_cuint8 = typename std::enable_if<is_thr_cuint8<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_thr_cint16 = typename std::enable_if<is_thr_cint16 <T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_thr_cuint16 = typename std::enable_if<is_thr_cuint16 <T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_thr_cint32 = typename std::enable_if<is_thr_cint32<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_thr_cuint32 = typename std::enable_if<is_thr_cuint32<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_thr_cint64 = typename std::enable_if<is_thr_cint64<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_thr_cuint64 = typename std::enable_if<is_thr_cuint64<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_thr_cfloat32 = typename std::enable_if<is_thr_cfloat32<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_thr_cfloat64 = typename std::enable_if<is_thr_cfloat64<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_thr_cfloat = typename std::enable_if<is_thr_cfloat<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_thr_cmplx = typename std::enable_if<is_thr_cmplx<T>::value, U>::type;
#endif
}

/* complex types - enable */
namespace mt
{
	template <class T, class U=void>
	using enable_if_cint8 = typename std::enable_if<is_cint8 <T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_cuint8 = typename std::enable_if<is_cuint8<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_cint16 = typename std::enable_if<is_cint16 <T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_cuint16 = typename std::enable_if<is_cuint16 <T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_cint32 = typename std::enable_if<is_cint32<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_cuint32 = typename std::enable_if<is_cuint32<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_cint64 = typename std::enable_if<is_cint64<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_cuint64 = typename std::enable_if<is_cuint64<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_cfloat32 = typename std::enable_if<is_cfloat32<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_cfloat64 = typename std::enable_if<is_cfloat64<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_cfloat = typename std::enable_if<is_cfloat<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_cmplx = typename std::enable_if<is_cmplx<T>::value, U>::type;

	template <class T, class U, class V=void>
	using enable_if_cmplx_cmplx = typename std::enable_if<is_cmplx<T>::value && is_cmplx<U>::value, V>::type;

	template <class T, class U=void>
	using enable_if_number = typename std::enable_if<is_number<T>::value, U>::type;

	template <class T, class U, class V=void>
	using enable_if_number_number = typename std::enable_if<is_number<T>::value && is_number<U>::value, V>::type;
}

/* eDev/eDim - verification */
namespace mt
{
	dt_bool is_edev_cpu(const eDev& edev)
	{
		return edev_cpu == edev;
	}

	dt_bool is_edev_gpu(const eDev& edev)
	{
		return edev_gpu == edev;
	}

	dt_bool is_edim_1(const eDim& edim)
	{
		return edim_1 == edim;
	}

	dt_bool is_edim_2(const eDim& edim)
	{
		return edim_2 == edim;
	}

	dt_bool is_edim_3(const eDim& edim)
	{
		return edim_3 == edim;
	}
}

/* ptr/eDev/eDim - enable */
namespace mt
{
	template <class T, class U=void>
	using enable_if_ptr = typename std::enable_if<std::is_pointer<T>::value, U>::type;

	template <eDev Dev, class U=void>
	using enable_if_edev_cpu = typename std::enable_if<Dev==edev_cpu, U>::type;

	template <eDev Dev, class U=void>
	using enable_if_edev_gpu = typename std::enable_if<Dev==edev_gpu, U>::type;

	template <eDim Dim, class U=void>
	using enable_if_edim_1 = typename std::enable_if<Dim==edim_1, U>::type;

	template <eDim Dim, class U=void>
	using enable_if_edim_2 = typename std::enable_if<Dim==edim_2, U>::type;

	template <eDim Dim, class U=void>
	using enable_if_edim_3 = typename std::enable_if<Dim==edim_3, U>::type;
}

/* types conversion */
namespace mt
{
	template <class T>
	using chg_btw_int32_int64 = typename std::conditional<is_int32<T>::value, dt_int64, dt_int32>::type;

	// ! change to complementary floating point type: this include complex number
#ifdef __CUDACC__
	template <class T>
	using chg_2_compl_float_type = typename std::conditional<is_float32<T>::value, dt_float64, \
		typename std::conditional<is_float64<T>::value, dt_float32, \
		typename std::conditional<is_std_cfloat32<T>::value, std::complex<dt_float64>, \
		typename std::conditional<is_std_cfloat64<T>::value, std::complex<dt_float32>, \
		typename std::conditional<is_thr_cfloat32<T>::value, thrust::complex<dt_float64>, thrust::complex<dt_float32>>::type>::type>::type>::type>::type;
#else
	template <class T>
	using chg_2_compl_float_type = typename std::conditional<is_float32<T>::value, dt_float64, \
		typename std::conditional<is_float64<T>::value, dt_float32, \
		typename std::conditional<is_std_cfloat32<T>::value, std::complex<dt_float64>, std::complex<dt_float32>>::type>::type>::type;
#endif

	// ! large type sel
	template <class T, class U>
	using sel_lg_type = typename std::conditional<is_float32<T>::value && is_float32<U>::value, dt_float32, \
		typename std::conditional<is_float64<T>::value && is_float32<U>::value, dt_float64, \
		typename std::conditional<is_float32<T>::value && is_float64<U>::value, dt_float64, \
		typename std::conditional<is_float<T>::value && is_int<U>::value, T, \
		typename std::conditional<is_int<T>::value && is_float<U>::value, U, dt_float64>::type>::type>::type>::type>::type;	/* types - enable */
}

/* detail_traits */
namespace mt
{
	namespace detail_traits
	{
		template <class T>
		struct check_type 
		{
			typedef void type;
		};

		/************************************ value type****************************************/
		template <class T, class U = void>
		struct sValue_type
		{
			using type = T;
		};

		template <class T>
		struct sValue_type<T, typename check_type<typename T::value_type>::type>
		{
			using type = typename T::value_type;
		};

		template <class T, class U = void>
		struct sValue_type_r
		{
			using type = typename sValue_type<T>::type;
		};

		template <class T>
		struct sValue_type_r<T, typename check_type<typename T::value_type::value_type>::type>
		{
			using type = typename T::value_type::value_type;
		};

		/************************************* size type ***************************************/
		template <class T, class U = void>
		struct sSize_type
		{
			using type = T;
		};

		template <class T>
		struct sSize_type<T, typename check_type<typename T::size_type>::type>
		{
			using type = typename T::size_type;
		};

		template <class T, class U = void>
		struct sSize_type_r
		{
			using type = typename sSize_type<T>::type;
		};

		template <class T>
		struct sSize_type_r<T, typename check_type<typename T::value_type::size_type>::type>
		{
			using type = typename T::value_type::size_type;
		};
	}
}
	
/* vector value type - extraction */
namespace mt
{
	template <class TVctr>
	using Value_type = typename detail_traits::sValue_type<TVctr>::type;

	template <class TVctr>
	using Value_type_r = typename detail_traits::sValue_type_r<TVctr>::type;
}

/* vector size type - extraction */
namespace mt
{
	template <class TVctr>
	using Size_type = typename detail_traits::sSize_type<TVctr>::type;

	template <class TVctr>
	using Size_type_r = typename detail_traits::sSize_type_r<TVctr>::type;
}

/* vector type - verification */
namespace mt
{
	template <class T>
	struct is_vctr_cpu: std::integral_constant<dt_bool, T::device==edev_cpu> {};

	template <class T>
	struct is_vctr_gpu: std::integral_constant<dt_bool, T::device==edev_gpu> {};

	template <class T>
	struct is_vctr: std::integral_constant<dt_bool, is_vctr_cpu<T>::value || is_vctr_gpu<T>::value> {};

	template <class T, class U>
	struct is_vctr_and_vctr: std::integral_constant<dt_bool, is_vctr<T>::value || is_vctr<U>::value> {};

	/***************************************************************************************/
	template <class T>
	struct is_rvctr_cpu: std::integral_constant<dt_bool, is_vctr_cpu<T>::value && is_real<typename T::value_type>::value> {};

	template <class T>
	struct is_vctr_rgpu: std::integral_constant<dt_bool, is_vctr_gpu<T>::value && is_real<typename T::value_type>::value> {};

	template <class T>
	struct is_rvctr: std::integral_constant<dt_bool, is_rvctr_cpu<T>::value || is_vctr_rgpu<T>::value> {};
		
	/***************************************************************************************/
	template <class T>
	struct is_cvctr_cpu: std::integral_constant<dt_bool, is_vctr_cpu<T>::value && is_cmplx<typename T::value_type>::value> {};

	template <class T>
	struct is_cvctr_gpu: std::integral_constant<dt_bool, is_vctr_gpu<T>::value && is_cmplx<typename T::value_type>::value> {};

	template <class T>
	struct is_cvctr: std::integral_constant<dt_bool, is_cvctr_cpu<T>::value || is_cvctr_gpu<T>::value> {};

	/***************************************************************************************/
	template <class T, class U>
	struct is_vctr_cpu_and_vctr_cpu: std::integral_constant<dt_bool, is_vctr_cpu<T>::value && is_vctr_cpu<U>::value> {};

	template <class T, class U>
	struct is_vctr_cpu_and_vctr_gpu: std::integral_constant<dt_bool, is_vctr_cpu<T>::value && is_vctr_gpu<U>::value> {};

	template <class T, class U>
	struct is_vctr_gpu_and_vctr_cpu: std::integral_constant<dt_bool, is_vctr_gpu<T>::value && is_vctr_cpu<U>::value> {};

	template <class T, class U>
	struct is_vctr_gpu_and_vctr_gpu: std::integral_constant<dt_bool, is_vctr_gpu<T>::value && is_vctr_gpu<U>::value> {};

	/***************************************************************************************/
	template <class T, class U>
	struct is_cvctr_cpu_and_vctr_cpu: std::integral_constant<dt_bool, is_cvctr_cpu<T>::value && is_vctr_cpu<U>::value> {};		
		
	template <class T, class U>
	struct is_cvctr_cpu_and_rvctr_cpu: std::integral_constant<dt_bool, is_cvctr_cpu<T>::value && is_rvctr_cpu<U>::value> {};

	template <class T, class U>
	struct is_cvctr_gpu_and_vctr_gpu: std::integral_constant<dt_bool, is_cvctr_gpu<T>::value && is_vctr_gpu<U>::value> {};
		
	template <class T, class U>
	struct is_cvctr_gpu_and_vctr_rgpu: std::integral_constant<dt_bool, is_cvctr_gpu<T>::value && is_vctr_rgpu<U>::value> {};

	template <class T, class U>
	struct is_cvctr_and_rvctr: std::integral_constant<dt_bool, is_cvctr<T>::value && is_rvctr<T>::value> {};
}

/* vector type - enable */
namespace mt
{
	template <class T, class U=void>
	using enable_if_vctr_cpu = typename std::enable_if<is_vctr_cpu<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_vctr_gpu = typename std::enable_if<is_vctr_gpu<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_vctr = typename std::enable_if<is_vctr<T>::value, U>::type;

	template <class T, class U, class V=void>
	using enable_if_vctr_and_vctr = typename std::enable_if<is_vctr_and_vctr<T, U>::value, V>::type;

	/***************************************************************************************/
	template <class T, class U=void>
	using enable_if_rvctr_cpu = typename std::enable_if<is_rvctr_cpu<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_vctr_rgpu = typename std::enable_if<is_vctr_rgpu<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_rvctr = typename std::enable_if<is_rvctr<T>::value, U>::type;

	/***************************************************************************************/
	template <class T, class U=void>
	using enable_if_cvctr_cpu = typename std::enable_if<is_cvctr_cpu<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_cvctr_gpu = typename std::enable_if<is_cvctr_gpu<T>::value, U>::type;

	template <class T, class U=void>
	using enable_if_cvctr = typename std::enable_if<is_cvctr<T>::value, U>::type;

	/***************************************************************************************/
	template <class T, class U, class V=void>
	using enable_if_vctr_cpu_and_vctr_cpu = typename std::enable_if<is_vctr_cpu_and_vctr_cpu<T, U>::value, V>::type;

	template <class T, class U, class V=void>
	using enable_if_vctr_cpu_and_vctr_gpu = typename std::enable_if<is_vctr_cpu_and_vctr_gpu<T, U>::value, V>::type;

	template <class T, class U, class V=void>
	using enable_if_vctr_gpu_and_vctr_cpu = typename std::enable_if<is_vctr_gpu_and_vctr_cpu<T, U>::value, V>::type;

	template <class T, class U, class V=void>
	using enable_if_vctr_gpu_and_vctr_gpu = typename std::enable_if<is_vctr_gpu_and_vctr_gpu<T, U>::value, V>::type;

	/***************************************************************************************/
	template <class T, class U, class V=void>
	using enable_if_cvctr_cpu_and_vctr_cpu = typename std::enable_if<is_cvctr_cpu_and_vctr_cpu<T, U>::value, V>::type;		
		
	template <class T, class U, class V=void>
	using enable_if_cvctr_cpu_and_rvctr_cpu = typename std::enable_if<is_cvctr_cpu_and_rvctr_cpu<T, U>::value, V>::type;

	template <class T, class U, class V=void>
	using enable_if_cvctr_gpu_and_vctr_gpu = typename std::enable_if<is_cvctr_gpu_and_vctr_gpu<T, U>::value, V>::type;

	template <class T, class U, class V=void>
	using enable_if_cvctr_gpu_and_vctr_rgpu = typename std::enable_if<is_cvctr_gpu_and_vctr_rgpu<T, U>::value, V>::type;

	template <class T, class U, class V=void>
	using enable_if_cvctr_and_rvctr = typename std::enable_if<is_cvctr_and_rvctr<T, U>::value, V>::type;
}

/* matrxi type - enable */
namespace mt
{
	template <class TVctr>
	enable_if_float32<Value_type<TVctr>, dt_int32>
	matrix_type(TVctr& vector){ return 1; }

	template <class TVctr>
	enable_if_float64<Value_type<TVctr>, dt_int32>
	matrix_type(TVctr& vector){ return 2; }

	template <class TVctr>
	enable_if_cfloat32<Value_type<TVctr>, dt_int32>
	matrix_type(TVctr& vector){ return 3; }

	template <class TVctr>
	enable_if_cfloat64<Value_type<TVctr>, dt_int32>
	matrix_type(TVctr& vector){ return 4; }
}