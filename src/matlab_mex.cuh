/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef MATLAB_MEX_H
	#define MATLAB_MEX_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include <algorithm>
	#include <cmath>
	#include <string>

	#include "const_enum.cuh"
	#include "type_traits_gen.cuh"
	#include "math.cuh"
	#include "memcpy.cuh"
	#include "r_2d.cuh"
	#include "r_3d.cuh"
	#include "cgpu_vctr.cuh"
	#include "cgpu_info.cuh"

	#include <mex.h>

	template <class T>
	mt::enable_if_real<T, void>
	mex_print(T* ptr, dt_int32 n)
	{
		for(auto ik = 0; ik <n; ik++)
		{
			mexPrintf("ik = %d, value = %8.5f\n", ik, dt_float64(ptr[ik]));
		}
	}

	template <class T>
	mt::enable_if_real<T, void>
	mex_print(const mt::pVctr_cpu_64<T>& vector)
	{
		for(auto ik = 0; ik <vector.size(); ik++)
		{
			mexPrintf("ik = %d, value = %8.5f\n", ik, dt_float64(vector[ik]));
		}
	}

	template <class T>
	mt::enable_if_cmplx<T, void>
	mex_print(const mt::pVctr_cpu_64<T>& vector)
	{
		for(auto ik = 0; ik <vector.size(); ik++)
		{
			mexPrintf("ik = %d, value = (%8.5f, %8.5f)\n", ik, dt_float64(vector[ik].real()), dt_float64(vector[ik].imag()));
		}
	}

	/***************************************************************************************/
	/********************************* get dimensions **************************************/
	/***************************************************************************************/
	CPU_EXEC_INL
	dt_uint32 mex_get_s0(const mxArray* mxA)
	{
		auto n_dim = mxGetNumberOfDimensions(mxA);
		auto dim = mxGetDimensions(mxA);

		return (n_dim > 0)?static_cast<dt_int32>(dim[0]):0;
	}

	CPU_EXEC_INL
	dt_uint32 mex_get_s1(const mxArray* mxA)
	{
		auto n_dim = mxGetNumberOfDimensions(mxA);
		auto dim = mxGetDimensions(mxA);

		return (n_dim>1)?static_cast<dt_int32>(dim[1]):0;
	}

	CPU_EXEC_INL
	dt_uint32 mex_get_s2(const mxArray* mxA)
	{
		auto n_dim = mxGetNumberOfDimensions(mxA);
		auto dim = mxGetDimensions(mxA);

		return (n_dim>2)?static_cast<dt_int32>(dim[2]):0;
	}

	CPU_EXEC_INL
	dt_uint32 mex_get_s3(const mxArray* mxA)
	{
		auto n_dim = mxGetNumberOfDimensions(mxA);
		auto dim = mxGetDimensions(mxA);

		return (n_dim>3)?static_cast<dt_int32>(dim[3]):0;
	}

	CPU_EXEC_INL
	dt_shape mex_get_mx_shape(const mxArray* mxA)
	{
		auto n_dim = min(mwSize(4), mxGetNumberOfDimensions(mxA));
		auto dim = mxGetDimensions(mxA);

		dt_shape shape{1, 1, 1, 1};

		for (auto ik = 0; ik<n_dim; ik++) 
		{
			shape[ik] = dt_uint32(dim[ik]);
		}

		return shape;
	}

	CPU_EXEC_INL
	dt_shape_st<dt_uint64> mex_get_mx_shape_uint64(const mxArray* mxA)
	{
		auto n_dim = min(mwSize(4), mxGetNumberOfDimensions(mxA));
		auto dim = mxGetDimensions(mxA);

		dt_shape_st<dt_uint64> shape{1, 1, 1, 1};

		for (auto ik = 0; ik<n_dim; ik++) 
		{
			shape[ik] = dt_uint64(dim[ik]);
		}

		return shape;
	}

	CPU_EXEC_INL
	dt_uint64 mex_get_size(const mxArray* mxA)
	{
		return dt_uint64(mxGetNumberOfElements(mxA));
	}

	/***************************************************************************************/
	/************************************ field exist **************************************/
	/***************************************************************************************/
	dt_bool mex_field_exist(const mxArray* mxA, const dt_int32& idx, const char *field_name)
	{
		const mxArray* fPtr = mxGetField(mxA, idx, field_name);

		return (fPtr!=nullptr);
	}

	dt_bool mex_field_exist(const mxArray* mxA, const char *field_name)
	{
		return mex_field_exist(mxA, 0, field_name);
	}

	/***************************************************************************************/
	/********************************** get string/number **********************************/
	/***************************************************************************************/
	std::string mex_get_string(const mxArray* mxA)
	{
		return std::string(mxArrayToString(mxA));
	}

	/***************************************************************************************/
	/*********************** pVector: template for reading matlab data *********************/
	/***************************************************************************************/
	template <class T>
	mt::pVctr_cpu_64<T> mex_get_pvctr(const mxArray* mxA) 
	{ 
		return mt::pVctr_cpu_64<T>();
	}

	/***************************************************************************************/
	/************************************ get real data ************************************/
	/***************************************************************************************/
	template <>
	mt::pVctr_cpu_64<dt_int8> mex_get_pvctr<dt_int8>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_int8*>(mxGetInt8s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_uint8> mex_get_pvctr<dt_uint8>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_uint8*>(mxGetUint8s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_int16> mex_get_pvctr<dt_int16>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_int16*>(mxGetInt16s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_uint16> mex_get_pvctr<dt_uint16>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_uint16*>(mxGetUint16s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_int32> mex_get_pvctr<dt_int32>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_int32*>(mxGetInt32s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_uint32> mex_get_pvctr<dt_uint32>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_uint32*>(mxGetUint32s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_int64> mex_get_pvctr<dt_int64>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_int64*>(mxGetInt64s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_uint64> mex_get_pvctr<dt_uint64>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_uint64*>(mxGetUint64s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_float32> mex_get_pvctr<dt_float32>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_float32*>(mxGetSingles(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_float64> mex_get_pvctr<dt_float64>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_float64*>(mxGetDoubles(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	/***************************************************************************************/
	/********************************** get complex data ***********************************/
	/***************************************************************************************/
	template <>
	mt::pVctr_cpu_64<dt_cint8> mex_get_pvctr<dt_cint8>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_cint8*>(mxGetComplexInt8s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_cuint8> mex_get_pvctr<dt_cuint8>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_cuint8*>(mxGetComplexUint8s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_cint16> mex_get_pvctr<dt_cint16>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_cint16*>(mxGetComplexInt16s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_cuint16> mex_get_pvctr<dt_cuint16>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_cuint16*>(mxGetComplexUint16s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_cint32> mex_get_pvctr<dt_cint32>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_cint32*>(mxGetComplexInt32s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_cuint32> mex_get_pvctr<dt_cuint32>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_cuint32*>(mxGetComplexUint32s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_cint64> mex_get_pvctr<dt_cint64>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_cint64*>(mxGetComplexInt64s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_cuint64> mex_get_pvctr<dt_cuint64>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_cuint64*>(mxGetComplexUint64s(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_cfloat32> mex_get_pvctr<dt_cfloat32>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_cfloat32*>(mxGetComplexSingles(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	template <>
	mt::pVctr_cpu_64<dt_cfloat64> mex_get_pvctr<dt_cfloat64>(const mxArray* mxA)
	{
		return {reinterpret_cast<dt_cfloat64*>(mxGetComplexDoubles(mxA)), mex_get_mx_shape_uint64(mxA)};
	}

	/***************************************************************************************/
	/*********************************** read data type ************************************/
	/***************************************************************************************/
	eData_Typ mex_get_data_type(const mxArray* mxA) 
	{
		auto category = mxGetClassID(mxA);

		if (mxIsComplex(mxA))
		{
			switch (category) 
			{
				case mxINT8_CLASS:
				{
					return edt_cint8;
				}
				case mxUINT8_CLASS:		
				{
					return edt_cuint8;
				}
				case mxINT16_CLASS:
				{
					return edt_cint16;
				}
				case mxUINT16_CLASS:
				{
					return edt_cuint16;
				}
				case mxINT32_CLASS:
				{
					return edt_cint32;
				}
				case mxUINT32_CLASS:
				{
					return edt_cuint32;
				}
				case mxINT64_CLASS:
				{
					return edt_cint64;
				}
				case mxUINT64_CLASS:
				{
					return edt_cuint64;
				}
				case mxSINGLE_CLASS:
				{
					return edt_cfloat32;
				} 
				case mxDOUBLE_CLASS:
				{
					return edt_cfloat64;
				}
			}
		}
		else
		{
			switch (category) 
			{
				case mxLOGICAL_CLASS:
				{
					return edt_bool;
				}				
				case mxINT8_CLASS:
				{
					return edt_int8;
				}
				case mxUINT8_CLASS:		
				{
					return edt_uint8;
				}
				case mxINT16_CLASS:
				{
					return edt_int16;
				}
				case mxUINT16_CLASS:
				{
					return edt_uint16;
				}
				case mxINT32_CLASS:
				{
					return edt_int32;
				}
				case mxUINT32_CLASS:
				{
					return edt_uint32;
				}
				case mxINT64_CLASS:
				{
					return edt_int64;
				}
				case mxUINT64_CLASS:
				{
					return edt_uint64;
				}
				case mxSINGLE_CLASS:
				{
					return edt_float32;
				} 
				case mxDOUBLE_CLASS:
				{
					return edt_float64;
				}
			}
		}

		return edt_none;
	}

	eData_Typ mex_get_data_type_from_field(const mxArray* mxA, const dt_int32& idx, const char *field_name)
	{
		if (mex_field_exist(mxA, field_name))
		{
			return mex_get_data_type(mxGetField(mxA, idx, field_name));

		}
		else
		{
			mexPrintf("Field name '%s' do not exist\n", field_name);
			return edt_none;
		}
	}

	eData_Typ mex_get_data_type_from_field(const mxArray* mxA, const char *field_name)
	{
		return mex_get_data_type_from_field(mxA, 0, field_name);
	}
	/***************************************************************************************/
	/********************** Vector: template for reading matlab data ***********************/
	/***************************************************************************************/
	template <class T>
	mt::Vctr_cpu<T> mex_get_vctr(const mxArray* mxA) 
	{
		auto data_type = mex_get_data_type(mxA);

		switch (data_type) 
		{
			case edt_bool:
			{
				return mex_get_pvctr<dt_bool>(mxA);
			}
			case edt_int8:
			{
				return mex_get_pvctr<dt_int8>(mxA);
			}
			case edt_uint8:		
			{
				return mex_get_pvctr<dt_uint8>(mxA);
			}
			case edt_int16:
			{
				return mex_get_pvctr<dt_int16>(mxA);
			}
			case edt_uint16:
			{
				return mex_get_pvctr<dt_uint16>(mxA);
			}
			case edt_int32:
			{
				return mex_get_pvctr<dt_int32>(mxA);
			}
			case edt_uint32:
			{
				return mex_get_pvctr<dt_uint32>(mxA);
			}
			case edt_int64:
			{
				return mex_get_pvctr<dt_int64>(mxA);
			}
			case edt_uint64:
			{
				return mex_get_pvctr<dt_uint64>(mxA);
			}
			case edt_float32:
			{
				return mex_get_pvctr<dt_float32>(mxA);
			} 
			case edt_float64:
			{
				return mex_get_pvctr<dt_float64>(mxA);
			}
		}

		return mt::Vctr_cpu<T>();
	}

	template <class T>
	mt::Vctr_r_2d_cpu<T> mex_get_vctr_r_2d(const mxArray* mxA) 
	{
		auto data_type = mex_get_data_type(mxA);

		switch (data_type) 
		{
			case edt_bool:
			{
				auto ptr = mex_get_pvctr<dt_bool>(mxA);
				return {ptr.data(), ptr.s0()};
			}
			case edt_int8:
			{
				auto ptr = mex_get_pvctr<dt_int8>(mxA);
				return {ptr.data(), ptr.s0()};
			}
			case edt_uint8:		
			{
				auto ptr = mex_get_pvctr<dt_uint8>(mxA);
				return {ptr.data(), ptr.s0()};
			}
			case edt_int16:
			{
				auto ptr = mex_get_pvctr<dt_int16>(mxA);
				return {ptr.data(), ptr.s0()};
			}
			case edt_uint16:
			{
				auto ptr = mex_get_pvctr<dt_uint16>(mxA);
				return {ptr.data(), ptr.s0()};
			}
			case edt_int32:
			{
				auto ptr = mex_get_pvctr<dt_int32>(mxA);
				return {ptr.data(), ptr.s0()};
			}
			case edt_uint32:
			{
				auto ptr = mex_get_pvctr<dt_uint32>(mxA);
				return {ptr.data(), ptr.s0()};
			}
			case edt_int64:
			{
				auto ptr = mex_get_pvctr<dt_int64>(mxA);
				return {ptr.data(), ptr.s0()};
			}
			case edt_uint64:
			{
				auto ptr = mex_get_pvctr<dt_uint64>(mxA);
				return {ptr.data(), ptr.s0()};
			}
			case edt_float32:
			{
				auto ptr = mex_get_pvctr<dt_float32>(mxA);
				return {ptr.data(), ptr.s0()};
			} 
			case edt_float64:
			{
				auto ptr = mex_get_pvctr<dt_float64>(mxA);
				return {ptr.data(), ptr.s0()};
			}
		}

		return mt::Vctr_cpu<T>();
	}

	/***************************************************************************************/
	/********************************** get dimensions *************************************/
	/***************************************************************************************/
	mt::Vctr_cpu<dt_int32> mex_get_vctr_shape(const mxArray* mxA) 
	{
		return mex_get_vctr<dt_int32>(mxA);
	}

	dt_shape mex_get_shape(const mxArray* mxA)
	{
		auto dim = mex_get_vctr<dt_int32>(mxA); // nx, ny, nz
		auto n_dim = min(4, dim.size_32());

		dt_shape shape{1, 1, 1, 1};

		for (auto ik = 0; ik<n_dim; ik++) 
		{
			shape[ik] = dim[ik];
		}

		shape.swap(0, 1);

		return shape;
	}

	/***************************************************************************************/
	/************************************ get r center *************************************/
	/***************************************************************************************/
	template <class T>
	mt::R_3d<T> mex_get_r_3d_r_c(const mxArray* mxA, mt::R_3d<T> r_c_0 = mt::R_3d<T>(0, 0, 0)) 
	{
		const auto vctr = mex_get_vctr<T>(mxA);
		const dt_int32 n_vctr = vctr.size_32();
		mt::R_3d<T> r = r_c_0;

		if (n_vctr<2)
			r = mt::R_3d<T>(vctr[0], vctr[0], vctr[0]);
		else if (n_vctr<3)
			r = mt::R_3d<T>(vctr[0], vctr[1], vctr[1]);
		else if (n_vctr<3)
			r = mt::R_3d<T>(vctr[0], vctr[1], vctr[2]);

		return r;
	}

	/***************************************************************************************/
	/*********************************** get pixel size ************************************/
	/***************************************************************************************/
	template <class T>
	mt::R_3d<T> mex_get_r_3d_pxs(const mxArray* mxA) 
	{
		const auto vctr = mex_get_vctr<T>(mxA);
		const dt_int32 n_vctr = vctr.size_32();
		mt::R_3d<T> r(T(1), T(1), T(1));

		if (n_vctr<1)
			r = mt::R_3d<T>(T(1), T(1), T(1));
		else if (n_vctr<2)
			r = mt::R_3d<T>(vctr[0], vctr[0], vctr[0]);
		else if (n_vctr<3)
			r = mt::R_3d<T>(vctr[0], vctr[1], vctr[1]);
		else if (n_vctr<3)
			r = mt::R_3d<T>(vctr[0], vctr[1], vctr[2]);

		return r;
	}

	/***************************************************************************************/
	/*********************************** get box size **************************************/
	/***************************************************************************************/
	template <class T>
	mt::R_3d<T> mex_get_r_3d_bs(const mxArray* mxA, mt::R_3d<T> f) 
	{
		return mex_get_r_3d_pxs<T>(mxA)*f;
	}	
	
	template <class T, class U>
	mt::R_3d<T> mex_get_r_3d_bs(const mxArray* mxA, mt::Vctr_cpu<U>& f) 
	{
		return mex_get_r_3d_bs<T>(mxA, mt::R_3d<T>(f.m_data, f.size_32()));
	}	

	template <class T, class ST>
	mt::R_3d<T> mex_get_r_3d_bs(const mxArray* mxA, const dt_shape_st<ST>& shape) 
	{
		return mex_get_r_3d_bs<T>(mxA, mt::R_3d<T>(shape[1], shape[0], shape[2]));
	}

	/***************************************************************************************/
	/************************************ get number ***************************************/
	/***************************************************************************************/
	template <class T>
	dt_bool mex_get_bool(const mxArray* mxA)
	{
		auto data = static_cast<dt_int32>(std::round(mxGetScalar(mxA)));
		return data > 0;
	}

	template <class T>
	T mex_get_enum(const mxArray* mxA)
	{
		auto data = static_cast<dt_int32>(std::round(mxGetScalar(mxA)));
		return static_cast<T>(data);
	}

	dt_cfloat64 mex_get_cmplx_num(const mxArray* mxA)
	{
		auto data_type = mex_get_data_type(mxA);

		switch (data_type) 
		{
			case edt_cint8:
			{
				auto r = mxGetComplexInt8s(mxA);
				return {dt_float64(r[0].real), dt_float64(r[0].imag)};
			}
			case edt_cuint8:		
			{
				auto r = mxGetComplexUint8s(mxA);
				return {dt_float64(r[0].real), dt_float64(r[0].imag)};
			}
			case edt_cint16:
			{
				auto r = mxGetComplexInt16s(mxA);
				return {dt_float64(r[0].real), dt_float64(r[0].imag)};
			}
			case edt_cuint16:
			{
				auto r = mxGetComplexUint16s(mxA);
				return {dt_float64(r[0].real), dt_float64(r[0].imag)};
			}
			case edt_cint32:
			{
				auto r = mxGetComplexInt32s(mxA);
				return {dt_float64(r[0].real), dt_float64(r[0].imag)};
			}
			case edt_cuint32:
			{
				auto r = mxGetComplexUint32s(mxA);
				return {dt_float64(r[0].real), dt_float64(r[0].imag)};
			}
			case edt_cint64:
			{
				auto r = mxGetComplexInt64s(mxA);
				return {dt_float64(r[0].real), dt_float64(r[0].imag)};
			}
			case edt_cuint64:
			{
				auto r = mxGetComplexUint64s(mxA);
				return {dt_float64(r[0].real), dt_float64(r[0].imag)};
			}
			case edt_cfloat32:
			{
				auto r = mxGetComplexSingles(mxA);
				return {dt_float64(r[0].real), dt_float64(r[0].imag)};
			} 
			case edt_cfloat64:
			{
				auto r = mxGetComplexDoubles(mxA);
				return {dt_float64(r[0].real), dt_float64(r[0].imag)};
			} 
			default:
			{
				return {mxGetScalar(mxA), dt_float64(0)};
			}
		}

		return {dt_float64(0), dt_float64(0)};
	}

	template <class T>
	T mex_get_num(const mxArray* mxA)
	{
		auto data = static_cast<dt_int32>(std::round(mxGetScalar(mxA)));
		return static_cast<T>(data);
	}

	template <>
	dt_int32 mex_get_num<dt_int32>(const mxArray* mxA)
	{
		return static_cast<dt_int32>(std::round(mxGetScalar(mxA)));
	}

	template <>
	dt_float32 mex_get_num<dt_float32>(const mxArray* mxA)
	{
		return static_cast<dt_float32>(mxGetScalar(mxA));
	}

	template <>
	dt_float64 mex_get_num<dt_float64>(const mxArray* mxA)
	{
		return mxGetScalar(mxA);
	}

	template <>
	dt_cfloat32 mex_get_num<dt_cfloat32>(const mxArray* mxA)
	{
		auto num = mex_get_cmplx_num(mxA);

		return {dt_float32(num.real()), dt_float32(num.imag())};
	}

	template <>
	dt_cfloat64 mex_get_num<dt_cfloat64>(const mxArray* mxA)
	{
		auto num = mex_get_cmplx_num(mxA);

		return {num.real(), num.imag()};
	}


	/***************************************************************************************/
	/******************************** get r_2d/r_3d data ***********************************/
	/***************************************************************************************/
	template <class T>
	mt::R_2d<T> mex_get_r_2d(const mxArray* mxA, mt::R_2d<T> r_0 = mt::R_2d<T>(0, 0))
	{
		if (mxA!=nullptr)
		{
			auto pvctr = mex_get_vctr<T>(mxA);
			return (vctr.size()>=2)?mt::R_2d<T>(vctr.data()):r_0;
		}
		else
		{
			return {T(0), T(0)};
		}
	}

	template <class T>
	mt::R_3d<T> mex_get_r_3d(const mxArray* mxA, mt::R_3d<T> r_0 = mt::R_3d<T>(0, 0, 0))
	{
		if (mxA!=nullptr)
		{
			auto vctr = mex_get_vctr<T>(mxA);
			return (vctr.size()>=3)?mt::R_3d<T>(vctr.data()):r_0;
		}
		else
		{
			return {T(0), T(0), T(0)};
		}
	}

	/***************************************************************************************/
	/********************************** get data by field **********************************/
	/***************************************************************************************/
	template <class T>
	mt::pVctr_cpu_64<T> mex_get_pvctr_from_field(const mxArray* mxA, const dt_int32& idx, const char *field_name)
	{
		if (mex_field_exist(mxA, field_name))
		{
			return mex_get_pvctr<T>(mxGetField(mxA, idx, field_name));
		}
		else
		{
			mexPrintf("Field name '%s' do not exist\n", field_name);
			return mex_get_pvctr<T>(nullptr);
		}
	}

	template <class T>
	mt::pVctr_cpu_64<T> mex_get_pvctr_from_field(const mxArray* mxA, const char *field_name)
	{
		return mex_get_pvctr_from_field<T>(mxA, 0, field_name);
	}

	template <class T>
	T mex_get_num_from_field(const mxArray* mxA, const dt_int32& idx, const char *field_name)
	{
		if (mex_field_exist(mxA, idx, field_name))
		{
			return mex_get_num<T>(mxGetField(mxA, idx, field_name));
		}
		else
		{
			return T(0);
		}
	}

	template <class T>
	T mex_get_num_from_field(const mxArray* mxA, const char *field_name)
	{
		return mex_get_num_from_field<T>(mxA, 0, field_name);
	}

	template <class T>
	mt::R_2d<T> mex_get_r_2d_from_field(const mxArray* mxA, const dt_int32& idx, const char *field_name, mt::R_2d<T> r_d = mt::R_2d<T>(0, 0))
	{
		if (mex_field_exist(mxA, field_name))
		{
			return mex_get_r_2d<T>(mxGetField(mxA, idx, field_name));
		}
		else
		{
			mexPrintf("Field name '%s' do not exist\n", field_name);
			return mex_get_r_2d<T>(nullptr);
		}
	}

	template <class T>
	mt::R_2d<T> mex_get_r_2d_from_field(const mxArray* mxA, const char *field_name, mt::R_2d<T> r_d = mt::R_2d<T>(0, 0))
	{
		return mex_get_r_2d_from_field<T>(mxA, 0, field_name, r_d);
	}

	template <class T>
	mt::R_3d<T> mex_get_r_3d_from_field(const mxArray* mxA, const dt_int32& idx, const char *field_name, mt::R_3d<T> r_d = mt::R_3d<T>(0, 0, 0))
	{
		if (mex_field_exist(mxA, field_name))
		{
			return mex_get_r_3d<T>(mxGetField(mxA, idx, field_name));
		}
		else
		{
			mexPrintf("Field name '%s' do not exist\n", field_name);
			return mex_get_r_3d<T>(nullptr);
		}
	}

	template <class T>
	mt::R_3d<T> mex_get_r_3d_from_field(const mxArray* mxA, const char *field_name, mt::R_3d<T> r_d = mt::R_3d<T>(0, 0, 0))
	{
		return mex_get_r_3d_from_field<T>(mxA, 0, field_name, r_d);
	}

	/***************************************************************************************/
	/*********************** template for creation of matlab data **************************/
	/***************************************************************************************/
	template <class T>
	mt::pVctr_cpu_64<T> mex_create_pVctr(dt_shape_64 shape, mxArray*& mex_data)
	{
		return mt::pVctr_cpu_64<T>(mex_data);
	}
	
	/***************************************************************************************/
	template <>
	mt::pVctr_cpu_64<dt_int8> mex_create_pVctr<dt_int8>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxINT8_CLASS, mxREAL);

		return mex_get_pvctr<dt_int8>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_uint8> mex_create_pVctr<dt_uint8>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxUINT8_CLASS, mxREAL);

		return mex_get_pvctr<dt_uint8>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_int16> mex_create_pVctr<dt_int16>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxINT16_CLASS, mxREAL);

		return mex_get_pvctr<dt_int16>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_uint16> mex_create_pVctr<dt_uint16>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxUINT16_CLASS, mxREAL);

		return mex_get_pvctr<dt_uint16>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_int32> mex_create_pVctr<dt_int32>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxINT32_CLASS, mxREAL);

		return mex_get_pvctr<dt_int32>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_uint32> mex_create_pVctr<dt_uint32>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxUINT32_CLASS, mxREAL);

		return mex_get_pvctr<dt_uint32>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_int64> mex_create_pVctr<dt_int64>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxINT64_CLASS, mxREAL);

		return mex_get_pvctr<dt_int64>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_uint64> mex_create_pVctr<dt_uint64>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxUINT64_CLASS, mxREAL);

		return mex_get_pvctr<dt_uint64>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_float32> mex_create_pVctr<dt_float32>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxSINGLE_CLASS, mxREAL);

		return mex_get_pvctr<dt_float32>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_float64> mex_create_pVctr<dt_float64>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxDOUBLE_CLASS, mxREAL);

		return mex_get_pvctr<dt_float64>(mex_data);
	}

	/***************************************************************************************/
	/******************************* create complex data ***********************************/
	/***************************************************************************************/
	template <>
	mt::pVctr_cpu_64<dt_cint8> mex_create_pVctr<dt_cint8>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxINT8_CLASS, mxCOMPLEX);

		return mex_get_pvctr<dt_cint8>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_cuint8> mex_create_pVctr<dt_cuint8>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxUINT8_CLASS, mxCOMPLEX);

		return mex_get_pvctr<dt_cuint8>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_cint16> mex_create_pVctr<dt_cint16>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxINT16_CLASS, mxCOMPLEX);

		return mex_get_pvctr<dt_cint16>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_cuint16> mex_create_pVctr<dt_cuint16>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxUINT16_CLASS, mxCOMPLEX);

		return mex_get_pvctr<dt_cuint16>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_cint32> mex_create_pVctr<dt_cint32>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxINT32_CLASS, mxCOMPLEX);

		return mex_get_pvctr<dt_cint32>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_cuint32> mex_create_pVctr<dt_cuint32>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxUINT32_CLASS, mxCOMPLEX);

		return mex_get_pvctr<dt_cuint32>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_cint64> mex_create_pVctr<dt_cint64>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxINT64_CLASS, mxCOMPLEX);

		return mex_get_pvctr<dt_cint64>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_cuint64> mex_create_pVctr<dt_cuint64>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxUINT64_CLASS, mxCOMPLEX);

		return mex_get_pvctr<dt_cuint64>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_cfloat32> mex_create_pVctr<dt_cfloat32>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxSINGLE_CLASS, mxCOMPLEX);

		return mex_get_pvctr<dt_cfloat32>(mex_data);
	}

	template <>
	mt::pVctr_cpu_64<dt_cfloat64> mex_create_pVctr<dt_cfloat64>(dt_shape_64 shape, mxArray*& mex_data)
	{
		mex_data = mxCreateNumericArray(shape.dim(), shape.data(), mxDOUBLE_CLASS, mxCOMPLEX);

		return mex_get_pvctr<dt_cfloat64>(mex_data);
	}

	template <class T>
	mt::pVctr_cpu_64<T> mex_create_num(mxArray*& mex_data)
	{
		return mex_create_pVctr<T>({1, 1}, mex_data);
	}

	/***************************************************************************************/
	/**************************** create number/r_2d/r_3d data *****************************/
	/***************************************************************************************/
	template <class T>
	mt::pVctr_cpu_64<T> mex_create_r_2d(mxArray*& mex_data)
	{
		return mex_create_pVctr<T>({2, 1}, mex_data);
	}

	template <class T>
	mt::pVctr_cpu_64<T> mex_create_r_3d(mxArray*& mex_data)
	{
		return mex_create_pVctr<T>({3, 1}, mex_data);
	}

	/***************************************************************************************/
	template <class T, class U>
	mt::enable_if_number<U, mt::pVctr_cpu_64<T>>
	mex_create_set_pVctr(mxArray*& mex_data, const mt::pVctr_cpu_64<U>& vctr)
	{
		auto pdata = mex_create_pVctr<T>(vctr.shape(), mex_data);
		mt::memcpy_cpu_cpu(pdata.data(), vctr.data(), vctr.size());
	
		return pdata;
	}

	template <class T, class U>
	mt::enable_if_r_2d<U, mt::pVctr_cpu_64<T>>
	mex_create_set_pVctr(mxArray*& mex_data, const mt::pVctr_cpu_64<U>& vctr)
	{
		auto pdata = mex_create_pVctr<T>({vctr.size(), 2}, mex_data);
		mt::memcpy_pos_cpu_cpu(pdata.data(), vctr.data(), vctr.size());
	
		return pdata;
	}

	template <class T, class U>
	mt::enable_if_r_3d<U, mt::pVctr_cpu_64<T>>
	mex_create_set_pVctr(mxArray*& mex_data, const mt::pVctr_cpu_64<U>& vctr)
	{
		auto pdata = mex_create_pVctr<T>({vctr.size(), 3}, mex_data);
		mt::memcpy_pos_cpu_cpu(pdata.data(), vctr.data(), vctr.size());
	
		return pdata;
	}	

	template <class T, class U>
	mt::pVctr_cpu_64<T> mex_create_set_pVctr(mxArray*& mex_data, const mt::Vctr_r_2d<U, mt::edev_cpu>& vctr_r_2d)
	{
		auto pdata = mex_create_pVctr<T>({vctr_r_2d.size(), 2}, mex_data);
		mt::memcpy_pos_cpu_cpu(pdata.data(), vctr_r_2d.data(), vctr_r_2d.size());
	
		return pdata;
	}

	template <class T, class U>
	mt::pVctr_cpu_64<T> mex_create_set_pVctr(mxArray*& mex_data, const mt::Vctr_r_3d<U, mt::edev_cpu>& vctr_r_3d)
	{
		auto pdata = mex_create_pVctr<T>({vctr_r_3d.size(), 3}, mex_data);
		mt::memcpy_pos_cpu_cpu(pdata.data(), vctr_r_3d.data(), vctr_r_3d.size());
	
		return pdata;
	}	

	template <class T, class U>
	mt::pVctr_cpu_64<T> mex_create_set_num(mxArray*& mex_data, const U& val)
	{
		auto pval = mex_create_num<T>(mex_data);
		pval[0] = T(val);
	
		return pval;
	}

	/***************************************************************************************/
	template <class T>
	mt::pVctr_cpu_64<T> mex_create_pVctr_field(dt_shape_64 shape, mxArray*& mex_data, const dt_int32& idx, const char *field_name)
	{
		mxArray* mxfield = nullptr;
		auto pfield = mex_create_pVctr<T>(shape, mxfield);
		mxSetField(mex_data, idx, field_name, mxfield);

		return pfield;
	}
	
	template <class T>
	mt::pVctr_cpu_64<T> mex_create_pVctr_field(dt_shape_64 shape, mxArray*& mex_data, const char *field_name)
	{
		return mex_create_pVctr_field<T>(shape, mex_data, 0, field_name);
	}
	
	template <class T>
	mt::pVctr_cpu_64<T> mex_create_num_field(mxArray*& mex_data, const dt_int32& idx, const char *field_name)
	{
		return mex_create_pVctr_field<T>({1, 1}, mex_data, 0, field_name);
	}
	
	template <class T>
	mt::pVctr_cpu_64<T> mex_create_num_field(mxArray*& mex_data, const char *field_name)
	{
		return mex_create_num_field<T>(mex_data, 0, field_name);
	}
	
	/***************************************************************************************/
	template <class T, class U>
	mt::pVctr_cpu_64<T> mex_create_set_pVctr_field(mxArray*& mex_data, const dt_int32& idx, const char *field_name, const mt::pVctr_cpu_64<U>& field_value)
	{
		auto pdata = mex_create_pVctr_field<T>(field_value.shape(), mex_data, idx, field_name);
		mt::memcpy_cpu_cpu(pdata.data(), field_value.data(), field_value.size());

		return pdata;
	}
	
	template <class T, class U>
	mt::pVctr_cpu_64<T> mex_create_set_pVctr_field(mxArray*& mex_data, const char *field_name, const mt::pVctr_cpu_64<U>& field_value)
	{
		return mex_create_set_pVctr_field<T>(mex_data, 0, field_name, field_value);
	}
	
	template <class T, class U>
	mt::pVctr_cpu_64<T> mex_create_set_num_field(mxArray*& mex_data, const dt_int32& idx, const char *field_name, const U& field_value)
	{
		auto pdata = mex_create_num_field<T>(mex_data, idx, field_name);
		pdata[0] = field_value;

		return pdata;
	}
	
	template <class T, class U>
	mt::pVctr_cpu_64<T> mex_create_set_num_field(mxArray*& mex_data, const char *field_name, const U& field_value)
	{
		 return mex_create_set_num_field<T>(mex_data, 0, field_name, field_value);
	}

	/***************************************************************************************/
	/********************************** run mexFunction ************************************/
	/***************************************************************************************/
	#define MEX_RUN_FCN_INT(FCN, idx)											\
	{																			\
		if (mxIsNumeric(prhs[idx]))												\
		{																		\
			auto category = mxGetClassID(prhs[idx]);							\
																				\
			switch (category)													\
			{																	\
				case mxINT8_CLASS:												\
				{																\
					FCN<dt_int8>(nlhs, plhs, nrhs, prhs);						\
				}																\
				break;															\
				case mxUINT8_CLASS:												\
				{																\
					FCN<dt_uint8>(nlhs, plhs, nrhs, prhs);						\
				}																\
				break;															\
				case mxINT16_CLASS:												\
				{																\
					FCN<dt_int16>(nlhs, plhs, nrhs, prhs);						\
				}																\
				break;															\
				case mxUINT16_CLASS:											\
				{																\
					FCN<dt_uint16>(nlhs, plhs, nrhs, prhs);						\
				}																\
				break;															\
				case mxINT32_CLASS:												\
				{																\
					FCN<dt_int32>(nlhs, plhs, nrhs, prhs);						\
				} 																\
				break;															\
				case mxUINT32_CLASS:											\
				{																\
					FCN<dt_uint32>(nlhs, plhs, nrhs, prhs);						\
				} 																\
				break;															\
				case mxINT64_CLASS:												\
				{																\
					FCN<dt_int64>(nlhs, plhs, nrhs, prhs);						\
				}																\
				break;															\
				case mxUINT64_CLASS:											\
				{																\
					FCN<dt_uint64>(nlhs, plhs, nrhs, prhs);						\
				}																\
				break;															\
			}																	\
		}																		\
	}

	#define MEX_RUN_FCN_FLOAT(FCN, idx)											\
	{																			\
		if (mxIsNumeric(prhs[idx]))												\
		{																		\
			auto category = mxGetClassID(prhs[idx]);							\
																				\
			switch (category)													\
			{																	\
				case mxSINGLE_CLASS:											\
				{																\
					FCN<dt_float32>(nlhs, plhs, nrhs, prhs);					\
				} 																\
				break;															\
				case mxDOUBLE_CLASS:											\
				{																\
					FCN<dt_float64>(nlhs, plhs, nrhs, prhs);					\
				}																\
				break;															\
				default:														\
				{																\
					mexPrintf("Error: Input data must be floating point\n");	\
				}																\
				break;															\
			}																	\
		}																		\
	}

	#define MEX_RUN_FCN_FLOAT_OUT(FCN, idx)										\
	{																			\
		if (mxIsNumeric(prhs[idx]))												\
		{																		\
			auto category = mxGetClassID(prhs[idx]);							\
																				\
			switch (category)													\
			{																	\
				case mxSINGLE_CLASS:											\
				{																\
					FCN<dt_float32>(nlhs, plhs, nrhs, prhs);					\
				} 																\
				break;															\
				default:														\
				{																\
					FCN<dt_float64>(nlhs, plhs, nrhs, prhs);					\
				}																\
				break;															\
			}																	\
		}																		\
	}

	#define MEX_RUN_FCN_REAL(FCN, idx)											\
	{																			\
		if (mxIsNumeric(prhs[idx]))												\
		{																		\
			auto category = mxGetClassID(prhs[idx]);							\
																				\
			switch (category)													\
			{																	\
				case mxINT8_CLASS:												\
				{																\
					FCN<dt_int8>(nlhs, plhs, nrhs, prhs);						\
				}																\
				break;															\
				case mxUINT8_CLASS:												\
				{																\
					FCN<dt_uint8>(nlhs, plhs, nrhs, prhs);						\
				}																\
				break;															\
				case mxINT16_CLASS:												\
				{																\
					FCN<dt_int16>(nlhs, plhs, nrhs, prhs);						\
				}																\
				break;															\
				case mxUINT16_CLASS:											\
				{																\
					FCN<dt_uint16>(nlhs, plhs, nrhs, prhs);						\
				}																\
				break;															\
				case mxINT32_CLASS:												\
				{																\
					FCN<dt_int32>(nlhs, plhs, nrhs, prhs);						\
				} 																\
				break;															\
				case mxUINT32_CLASS:											\
				{																\
					FCN<dt_uint32>(nlhs, plhs, nrhs, prhs);						\
				} 																\
				break;															\
				case mxINT64_CLASS:												\
				{																\
					FCN<dt_int64>(nlhs, plhs, nrhs, prhs);						\
				}																\
				break;															\
				case mxUINT64_CLASS:											\
				{																\
					FCN<dt_uint64>(nlhs, plhs, nrhs, prhs);						\
				}																\
				break;															\
				case mxSINGLE_CLASS:											\
				{																\
					FCN<dt_float32>(nlhs, plhs, nrhs, prhs);					\
				} 																\
				break;															\
				case mxDOUBLE_CLASS:											\
				{																\
					FCN<dt_float64>(nlhs, plhs, nrhs, prhs);					\
				}																\
				break;															\
			}																	\
		}																		\
	}


	/***************************************************************************************/
	/**************************** mex system configuration *********************************/
	/***************************************************************************************/
	dt_bool mex_is_system_config(const mxArray* mxA)
	{
		return mxIsStruct(mxA) && mex_field_exist(mxA, "device") && mex_field_exist(mxA, "precision");
	}
	
	mt::System_Config mex_read_system_config(const mxArray* mex_in)
	{
		mt::System_Config system_config;
	
		if (mex_is_system_config(mex_in))
		{
			system_config.device = mex_get_num_from_field<mt::eDev>(mex_in, "device");
			system_config.precision = mex_get_num_from_field<mt::ePrecision>(mex_in, "precision");
			// system_config.cpu_n_proc = 1;
			system_config.cpu_n_proc = mex_get_num_from_field<dt_int32>(mex_in, "cpu_n_proc");
			system_config.cpu_n_thread = mex_get_num_from_field<dt_int32>(mex_in, "cpu_n_thread");
	
			auto gpu_device = mex_get_pvctr_from_field<dt_float64>(mex_in, "gpu_device");
			system_config.gpu_device.assign(gpu_device.begin(), gpu_device.end());
			system_config.gpu_n_stream = 1;
			// system_config.gpu_n_stream = mex_get_num_from_field<dt_int32>(mex_in, "gpu_n_stream");
			system_config.idx_0 = 1;

			system_config.validate_parameters();
		}
		else
		{
			system_config.cpu_n_thread = 1;
			system_config.idx_0 = 0;
		}

		return system_config;
	}

	mt::System_Config mex_read_set_system_config(const mxArray* mex_in)
	{
		auto system_config = mex_read_system_config(mex_in);

		if (system_config.is_gpu())
		{
			system_config.set_gpu();
		}

		return system_config;
	}

	/***************************************************************************************/
	/*************** this option set the output type using the in data ******************/
	/***************************************************************************************/
	#define MEX_RUN_FCN_INT_SYS_CONF_IN_TYP(FCN, idx_0)							\
	{																			\
		auto bb_system_config = mex_is_system_config(prhs[0]);					\
		auto idx = (bb_system_config)?idx_0+1:idx_0;							\
		MEX_RUN_FCN_INT(FCN, idx);												\
	}

	#define MEX_RUN_FCN_FLOAT_SYS_CONF_IN_TYP(FCN, idx_0)						\
	{																			\
		auto bb_system_config = mex_is_system_config(prhs[0]);					\
		auto idx = (bb_system_config)?idx_0+1:idx_0;							\
		MEX_RUN_FCN_FLOAT(FCN, idx);											\
	}

	#define MEX_RUN_FCN_REAL_SYS_CONF_IN_TYP(FCN, idx_0)						\
	{																			\
		auto bb_system_config = mex_is_system_config(prhs[0]);					\
		auto idx = (bb_system_config)?idx_0+1:idx_0;							\
		MEX_RUN_FCN_REAL(FCN, idx);												\
	}

	/***************************************************************************************/
	/************ this option set the output type using system configuration ***************/
	/***************************************************************************************/
	#define MEX_RUN_FCN_DEV(system_config, FCN, T)								\
	if (system_config.is_cpu())													\
	{																			\
		FCN<T, mt::edev_cpu>system_config, (nlhs, plhs, nrhs, prhs);			\
	}																			\
	else if (system_config.is_gpu())											\
	{																			\
		FCN<T, mt::edev_gpu>(system_config, nlhs, plhs, nrhs, prhs);			\
	}

#ifdef __CUDACC__
	#define MEX_RUN_FCN_FLOAT_SYS_CONF(FCN, idx_0)										\
	{																					\
		auto system_config = mex_read_set_system_config(prhs[idx_0]);					\
																						\
		if (system_config.is_float32_cpu())												\
		{																				\
			FCN<dt_float32, mt::edev_cpu>(system_config, nlhs, plhs, nrhs, prhs);		\
		}																				\
		else if (system_config.is_float64_cpu())										\
		{																				\
			FCN<dt_float64, mt::edev_cpu>(system_config, nlhs, plhs, nrhs, prhs);		\
		}																				\
		else if (system_config.is_float32_gpu())										\
		{																				\
			FCN<dt_float32, mt::edev_gpu>(system_config, nlhs, plhs, nrhs, prhs);		\
		}																				\
		else if (system_config.is_float64_gpu())										\
		{																				\
			FCN<dt_float64, mt::edev_gpu>(system_config, nlhs, plhs, nrhs, prhs);		\
		}																				\
	}
#else
	#define MEX_RUN_FCN_FLOAT_SYS_CONF(FCN, idx_0)										\
	{																					\
		auto system_config = mex_read_set_system_config(prhs[idx_0]);					\
																						\
		if (system_config.is_float32_cpu())												\
		{																				\
			FCN<dt_float32, mt::edev_cpu>(system_config, nlhs, plhs, nrhs, prhs);		\
		}																				\
		else if (system_config.is_float64_cpu())										\
		{																				\
			FCN<dt_float64, mt::edev_cpu>(system_config, nlhs, plhs, nrhs, prhs);		\
		}																				\
	}
	#endif

#endif