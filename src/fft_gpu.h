/*
* This file is part of Multem.
* Copyright 2022 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifdef __CUDACC__
	#include <cuda.h>
	#include <cuda_runtime.h>
	#include <cufft.h>
#endif

#include "macros.h"
#include "const_enum.h"
#include "math_mt.h"
#include "vctr_gpu.h"
#include "stream_gpu.h"

/* forward declaration */
namespace mt
{
#ifndef FFT_H
	#define FFT_H
	template <class T, eDev Dev> class FFT;
#endif
}

#ifdef __CUDACC__
/* gpu fourier transform - dt_float32 */
namespace mt
{
	template <>
	class FFT<dt_float32, edev_gpu>
	{
	public:
		using value_type = dt_float32;
		static const eDev device = edev_gpu;

		using TVctr_c = Vctr<complex<value_type>, edev_gpu>;

		FFT();

		FFT(const dt_int32& s0, Stream_gpu* pstream = nullptr);

		FFT(const dt_int32& s0, const dt_int32& s1, Stream_gpu* pstream = nullptr);

		~FFT();

		void cleanup();

		void destroy_plan();

		void create_plan_1d(const dt_int32& s0, dt_int32 n_thread=1);

		void create_plan_1d_batch(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread=1);

		void create_plan_2d(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread=1);

		void create_plan_2d_batch(const dt_int32& s0, const dt_int32& s1, const dt_int32& nz, dt_int32 n_thread=1);

		void set_stream(cudaStream_t &stream);

		void forward(TVctr_c& mx);

		void inverse(TVctr_c& mx);

		void forward(TVctr_c& mx_i, TVctr_c& mx_o);

		void inverse(TVctr_c& mx_i, TVctr_c& mx_o);

	private:
		cufftHandle plan_forward;
		cufftHandle plan_backward;
	};
}

/* gpu fourier transform - dt_float64 */
namespace mt
{
	template <>
	class FFT<dt_float64, edev_gpu>
	{
	public:
		using value_type = dt_float64;
		static const eDev device = edev_gpu;

		using TVctr_c = Vctr<complex<value_type>, edev_gpu>;

		FFT();

		FFT(const dt_int32& s0, Stream_gpu* pstream = nullptr);

		FFT(const dt_int32& s0, const dt_int32& s1, Stream_gpu* pstream = nullptr);

		~FFT();

		void cleanup();

		void destroy_plan();

		void create_plan_1d(const dt_int32& s0, dt_int32 n_thread=1);

		void create_plan_1d_batch(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread=1);

		void create_plan_2d(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread=1);

		void create_plan_2d_batch(const dt_int32& s0, const dt_int32& s1, const dt_int32& s2, dt_int32 n_thread=1);

		void forward(TVctr_c& mx);

		void inverse(TVctr_c& mx);

		void forward(TVctr_c& mx_i, TVctr_c& mx_o);

		void inverse(TVctr_c& mx_i, TVctr_c& mx_o);

	private:
		cufftHandle plan_forward;
		cufftHandle plan_backward;
	};
}

#include "detail/fft_gpu.inl"

#endif