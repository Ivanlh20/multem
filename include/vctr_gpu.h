/*
* This file is part of Multem.
* Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
*
* Multem is destroy software: you can redistribute it and/or modify
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

#include "vctr_cpu.h"

/* template definition */
namespace mt
{
#ifndef VCTR_DEC
	#define VCTR_DEC
	template <class T, eDev Dev> class Vctr;

	template <class T, eDev Dev, class ST> class pVctr;
#endif
}

/* derived class */
namespace mt
{
	template <class T>
	using Vctr_gpu = Vctr<T, edev_gpu>;

	/***************************************************************************************/
	using Vctr_uint32_gpu = Vctr_gpu<dt_uint32>;

	using Vctr_int32_gpu = Vctr_gpu<dt_int32>;

	using Vctr_uint64_gpu = Vctr_gpu<dt_uint64>;

	using Vctr_int64_gpu = Vctr_gpu<dt_int64>;

	/***************************************************************************************/
	template <class T>
	using Vctr_r_2d_gpu = Vctr<R_2d<T>, edev_gpu>;

	/***************************************************************************************/
	template <class T>
	using Vctr_r_3d_gpu = Vctr<R_3d<T>, edev_gpu>;

	/***************************************************************************************/
	template <class T>
	using Vctr_Mx_2x2_gpu = Vctr<Mx_2x2<T>, edev_gpu>;

	/***************************************************************************************/
	template <class T>
	using Vctr_Mx_3x3_gpu = Vctr<Mx_3x3<T>, edev_gpu>;
}

#ifdef __CUDACC__
/* gpu vector */
namespace mt
{
	template <class T>
	class Vctr<T, edev_gpu>
	{
	public:
		using value_type = T;
		using size_type = dt_int64;
		static const eDev device = edev_gpu;

		mutable T* m_data;
		size_type m_s0;
		size_type m_s1;
		size_type m_s2;
		size_type m_s3;
		size_type m_size;
		size_type m_capacity;

		size_type m_pitch_s1;
		size_type m_pitch_s2;
		size_type m_pitch_s3;

		/************************************* constructors ************************************/
		explicit Vctr();

		Vctr(const dt_init_list_f64& data);

		Vctr(size_type s0);

		Vctr(size_type s0, const T& value);

		explicit Vctr(const dt_shape_st<size_type>& shape);

		explicit Vctr(const dt_shape_st<size_type>& shape, const T& value);

		/* copy constructor */
		Vctr(const Vctr<T, edev_gpu>& vctr);

		/* Move constructor */
		Vctr(Vctr<T, edev_gpu>&& vctr);

		/* converting constructor */
		template <class U>
		Vctr(const Vctr<U, edev_gpu>& vctr);

		template <class U>
		Vctr(U* first, U* last);

		template <class U, class V=T, class = enable_if_r_nd<V>>
		Vctr(U *p, dt_int64 n_p, size_type icol=0);

		template <class U>
		Vctr(const std::vector<U>& vctr);

		// from gpu pVctr to Vctr
		template <class U, class STU>
		Vctr(const pVctr<U, edev_gpu, STU>& pvctr);

		// from cpu pVctr to Vctr
		template <class U, class STU>
		Vctr(const pVctr<U, edev_cpu, STU>& pvctr);

		template <class U>
		Vctr(const Vctr<U, edev_cpu>& vctr);

		template <class U>
		Vctr(const thrust::host_vector<U>& vctr);

		template <class U>
		Vctr(const thrust::device_vector<U>& vctr);

		~Vctr();

		/******************************** assignment operators *********************************/
		/* copy assignment operator */
		Vctr<T, edev_gpu>& operator=(const Vctr<T, edev_gpu>& vctr);

		/* Move assignment operator */
		Vctr<T, edev_gpu>& operator=(Vctr<T, edev_gpu>&& vctr);

		/* converting assignment operator */
		template <class U>
		Vctr<T, edev_gpu>& operator=(const Vctr<U, edev_gpu>& vctr);

		template <class U>
		Vctr<T, edev_gpu>& operator=(const std::vector<U>& vctr);

		template <class U>
		Vctr<T, edev_gpu>& operator=(const Vctr<U, edev_cpu>& vctr);

		template <class U>
		Vctr<T, edev_gpu>& operator=(const thrust::host_vector<U>& vctr);

		template <class U>
		Vctr<T, edev_gpu>& operator=(const thrust::device_vector<U>& vctr);

		template <class U>
		void assign(const Vctr<U, edev_gpu>& vctr, T* pvctr_cpu = nullptr);

		template <class U>
		void assign(U* first, U* last);
			
		// from gpu pVctr to Vctr
		template <class U, class STU>
		void assign(const pVctr<U, edev_gpu, STU>& pvctr);

		template <class U>
		void assign(const thrust::device_ptr<U>& first, const thrust::device_ptr<U>& last, T* pvctr_cpu = nullptr);

		// from cpu pVctr to Vctr
		template <class U, class STU>
		void assign(const pVctr<U, edev_cpu, STU>& pvctr);

		template <class U>
		void assign(const std::vector<U>& vctr, T* pvctr_cpu = nullptr);

		template <class U>
		void assign(const Vctr<U, edev_cpu>& vctr, T* pvctr_cpu = nullptr);

		template <class U>
		void assign(const thrust::host_vector<U>& vctr, T* pvctr_cpu = nullptr);

		template <class U>
		void assign(const thrust::device_vector<U>& vctr, T* pvctr_cpu = nullptr);

		/**************** user define conversion operators *******************/
		pVctr_gpu_32<T> ptr_32() const;

		pVctr_gpu_64<T> ptr_64() const;

		/* user define conversion for pointer Vctr */
		operator pVctr_gpu_32<T>() const;

		operator pVctr_gpu_64<T>() const;

		/* user define conversion */
		operator std::vector<T>() const;

		/* user define conversion in which T is the complemented precision */
		template <class U = chg_2_compl_float_type<T>>
		operator std::vector<U>() const;

		/* user define conversion to output type std::vector<thrust::complex<dt_float32>> */
		template <class U=T, class = enable_if_std_cfloat<U>>
		operator std::vector<thrust::complex<dt_float32>>() const;

		/* user define conversion to output type std::vector<thrust::complex<dt_float64>> */
		template <class U=T, class = enable_if_std_cfloat<U>>
		operator std::vector<thrust::complex<dt_float64>>() const;

		/* user define conversion */
		operator thrust::host_vector<T>() const;

		/* user define conversion in which T is the complemented precision */
		template <class U = chg_2_compl_float_type<T>>
		operator thrust::host_vector<U>() const;

		/* user define conversion to output type thrust::host_vector<std::complex<dt_float32>> */
		template <class U=T, class = enable_if_thr_cfloat<U>>
		operator thrust::host_vector<std::complex<dt_float32>>() const;

		/* user define conversion to output type thrust::host_vector<std::complex<dt_float64>> */
		template <class U=T, class = enable_if_thr_cfloat<U>>
		operator thrust::host_vector<std::complex<dt_float64>>() const;

		/***************************************************************************************/
		/* user define conversion */
		operator thrust::device_vector<T>() const;

		/* user define conversion in which T is the complemented precision */
		template <class U = chg_2_compl_float_type<T>>
		operator thrust::device_vector<U>() const;

		/***************************************************************************************/
		template <class U>
		void cpy_to_cpu_ptr(U* pdata, size_type n_data, T* pvctr_cpu = nullptr);

		template <class U>
		void cpy_to_cpu_ptr(U* first, U* last, T* pvctr_cpu = nullptr);

		template <class U>
		void cpy_to_gpu_ptr(U* pdata, size_type n_data, U* pvctr_cpu = nullptr);

		template <class U>
		void cpy_to_gpu_ptr(U* first, U* last, U* pvctr_cpu = nullptr);

		/***************************************************************************************/
		template <class U, class V=T, class = enable_if_cmplx<V>>
		void cpy_real_to_cpu_ptr(U* pdata, size_type n_data, T* pvctr_cpu = nullptr);

		template <class U, class V=T, class = enable_if_cmplx<V>>
		void cpy_real_to_cpu_ptr(U* first, U* last, T* pvctr_cpu = nullptr);

		template <class U, class V=T, class = enable_if_cmplx<V>>
		void cpy_real_to_gpu_ptr(U* pdata, size_type n_data, U* pvctr_cpu = nullptr);

		template <class U, class V=T, class = enable_if_cmplx<V>>
		void cpy_real_to_gpu_ptr(U* first, U* last, U* pvctr_cpu = nullptr);

		/***************************************************************************************/
		template <class U, class V=T, class = enable_if_cmplx<V>>
		void cpy_imag_to_cpu_ptr(U* pdata, size_type n_data, T* pvctr_cpu = nullptr);

		template <class U, class V=T, class = enable_if_cmplx<V>>
		void cpy_imag_to_cpu_ptr(U* first, U* last, T* pvctr_cpu = nullptr);

		template <class U, class V=T, class = enable_if_cmplx<V>>
		void cpy_imag_to_gpu_ptr(U* pdata, size_type n_data, U* pvctr_cpu = nullptr);

		template <class U, class V=T, class = enable_if_cmplx<V>>
		void cpy_imag_to_gpu_ptr(U* first, U* last, U* pvctr_cpu = nullptr);

		/***************************************************************************************/
		template <class U, class V=T, class = enable_if_r_nd<V>>
		void set_pos_from_cpu_ptr(U *p, size_type n_p, dt_int64 icol=0);

		template <class U, class V=T, class = enable_if_r_nd<V>>
		void cpy_pos_to_cpu_ptr(U *p, size_type icol=0);

		/***************************************************************************************/
		void resize(const dt_shape_st<size_type>& shape);

		void resize(const dt_shape_st<size_type>& shape, const T& value);

		void reserve(const dt_shape_st<size_type>& shape);

		void shrink_to_fit();

		template <class U>
		void push_back(const U& val);

		template <class U>
		void push_back(const Vctr<U, edev_cpu>& vctr);

		void pop_back();

		void fill(T val);

		/***************************************************************************************/
		size_type s0() const;

		size_type s1() const;

		size_type s2() const;

		size_type s3() const;

		dt_int32 s0_32() const;

		dt_int32 s1_32() const;

		dt_int32 s2_32() const;

		dt_int32 s3_32() const;
			
		dt_int64 s0_64() const;

		dt_int64 s1_64() const;

		dt_int64 s2_64() const;

		dt_int64 s3_64() const;

		size_type s0h() const;

		size_type s1h() const;

		size_type s2h() const;

		size_type s3h() const;

		dt_shape_st<size_type> shape() const;

		dt_shape_st<size_type> shape_2d_trs() const;

		size_type shape_size() const;

		size_type pitch_s1() const;

		size_type pitch_s2() const;

		size_type pitch_s3() const;

		size_type size() const;

		dt_int32 size_32() const;

		dt_int64 size_64() const;

		iGrid_1d igrid_1d() const;

		iGrid_2d igrid_2d() const;

		iGrid_3d igrid_3d() const;

		iGrid_1d_64 igrid_1d_64() const;

		iGrid_2d_64 igrid_2d_64() const;

		iGrid_3d_64 igrid_3d_64() const;

		size_type capacity() const;

		dt_bool empty() const;

		dt_bool is_1d() const;

		void clear();

		void clear_shrink_to_fit();

		size_type sub_2_ind(const size_type& ix_0) const;

		size_type sub_2_ind(const size_type& ix_0, const size_type& ix_1) const;

		size_type sub_2_ind(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2) const;

		size_type sub_2_ind(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2, const size_type& ix_3) const;

		// T& operator[](const size_type& iy);

		// const T& operator[](const size_type& iy) const;

		// T& operator()(const size_type& iy);

		// const T& operator()(const size_type& iy) const;

		// T& operator()(const size_type& ix_0, const size_type& ix_1);

		// const T& operator()(const size_type& ix_0, const size_type& ix_1) const;

		// T& operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2);

		// const T& operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2) const;

		// T& operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2, const size_type& ix_3);

		// const T& operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2, const size_type& ix_3) const;

		thrust::device_ptr<T> begin() noexcept;

		const thrust::device_ptr<T> begin() const;

		thrust::device_ptr<T> end() noexcept;

		const thrust::device_ptr<T> end() const;

		T* data();

		const T* data() const;

		template <class U>
		U data_cast();

		template <class U>
		const U data_cast() const;

		// T& front();

		// const T& front() const;

		// T& back();

		// const T& back() const;

		// set shape
		void set_shape(const dt_shape_st<size_type>& shape);

		void trs_shape_2d();

		template <class U, eDev Dev>
		void set_shape(const Vctr<U, Dev>& vctr, dt_bool bb_size = true);
		
		FCNS_DEF_GPU_GRID_BLK_VCTR;

	private:

		void set_picth();

		void set_capacity(size_type size_r);

		void set_shape_cstr(dt_shape_st<size_type>& shape);

		// reallocate and copy memory
		void allocate(dt_shape_st<size_type> shape, dt_bool bb_reserve=false);

		// initialization
		void init();

		// destroy memory on the device
		void destroy();
	};
}

#include "../src/vctr_gpu.inl"

#endif