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

#include "fft_gpu.h"

/* gpu fourier transform - dt_float32 */
namespace mt
{
	FFT<dt_float32, edev_gpu>::FFT(): plan_forward(0), plan_backward(0) {}

	FFT<dt_float32, edev_gpu>::FFT(const dt_int32& s0, Stream_gpu* pstream): FFT()
	{ 
		create_plan_1d(s0);
	}

	FFT<dt_float32, edev_gpu>::FFT(const dt_int32& s0, const dt_int32& s1, Stream_gpu* pstream): FFT()
	{ 
		create_plan_2d(s0, s1);
	}

	FFT<dt_float32, edev_gpu>::~FFT()
	{
		destroy_plan();
	}

	void FFT<dt_float32, edev_gpu>::cleanup()
	{
		destroy_plan();
	}

	void FFT<dt_float32, edev_gpu>::destroy_plan()
	{	
		if (plan_backward == 0)
		{
			return;
		}

		cudaDeviceSynchronize();

		cufftDestroy(plan_forward);
		plan_backward = plan_forward = 0;
	}

	void FFT<dt_float32, edev_gpu>::create_plan_1d(const dt_int32& s0, dt_int32 n_thread)
	{
		destroy_plan();

		cufftPlan1d(&plan_forward, s0, CUFFT_C2C, 1);

		plan_backward = plan_forward;
	}

	void FFT<dt_float32, edev_gpu>::create_plan_1d_batch(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread)
	{
		destroy_plan();

		dt_int32 rank = 1;							// Dimensionality of the transform (1, 2, or 3). 
		dt_int32 n[] = {s0};						// 1d transforms of length s0*s1
		dt_int32 how_many = s1;
		dt_int32 idist = s0, odist = s0;			// distance between two successive in elements in the least significant
		dt_int32 istride = 1, ostride = 1;			// distance between two elements in the same column
		dt_int32 *inembed = n, *onembed = n;		// Pointer of size rank that indicates the storage dimensions

		cufftPlanMany(&plan_forward, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, how_many);

		plan_backward = plan_forward;
	}

	void FFT<dt_float32, edev_gpu>::create_plan_2d(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread)
	{
		destroy_plan();

		// https:// docs.nvidia.com/cuda/cufft/index.html
		cufftPlan2d(&plan_forward, s1, s0, CUFFT_C2C);

		plan_backward = plan_forward;
	}

	void FFT<dt_float32, edev_gpu>::create_plan_2d_batch(const dt_int32& s0, const dt_int32& s1, const dt_int32& nz, dt_int32 n_thread)
	{
		destroy_plan();

		dt_int32 rank = 2;							// Dimensionality of the transform (1, 2, or 3). 
		dt_int32 n[] = {s1, s0};					// 2d transforms of length s0*s1
		dt_int32 how_many = nz;
		dt_int32 idist = s0*s1, odist = s0*s1;		// distance between two successive in elements in the least significant
		dt_int32 istride = 1, ostride = 1;			// distance between two elements in the same column
		dt_int32 *inembed = n, *onembed = n;		// Pointer of size rank that indicates the storage dimensions

		cufftPlanMany(&plan_forward, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, how_many);

		plan_backward = plan_forward;
	}

	void FFT<dt_float32, edev_gpu>::set_stream(cudaStream_t &stream)
	{
		cufftSetStream(plan_forward, stream);
	}

	void FFT<dt_float32, edev_gpu>::forward(TVctr_c& mx)
	{
		forward(mx, mx);
	}

	void FFT<dt_float32, edev_gpu>::inverse(TVctr_c& mx)
	{
		inverse(mx, mx);
	}

	void FFT<dt_float32, edev_gpu>::forward(TVctr_c& mx_i, TVctr_c& mx_o)
	{
		auto pdata_i = mx_i.data_cast<cufftComplex*>();
		auto pdata_o = mx_o.data_cast<cufftComplex*>();

		cufftExecC2C(plan_forward, pdata_i, pdata_o, CUFFT_FORWARD);
	}

	void FFT<dt_float32, edev_gpu>::inverse(TVctr_c& mx_i, TVctr_c& mx_o)
	{
		auto pdata_i = mx_i.data_cast<cufftComplex*>();
		auto pdata_o = mx_o.data_cast<cufftComplex*>();

		cufftExecC2C(plan_backward, pdata_i, pdata_o, CUFFT_INVERSE);
	}
}

/* gpu fourier transform - dt_float64 */
namespace mt
{
	FFT<dt_float64, edev_gpu>::FFT(): plan_forward(0), plan_backward(0) {}

	FFT<dt_float64, edev_gpu>::FFT(const dt_int32& s0, Stream_gpu* pstream): FFT()
	{ 
		create_plan_1d(s0);
	}

	FFT<dt_float64, edev_gpu>::FFT(const dt_int32& s0, const dt_int32& s1, Stream_gpu* pstream): FFT()
	{ 
		create_plan_2d(s0, s1);
	}

	FFT<dt_float64, edev_gpu>::~FFT()
	{
		destroy_plan();
	}

	void FFT<dt_float64, edev_gpu>::cleanup()
	{
		destroy_plan();
	}

	void FFT<dt_float64, edev_gpu>::destroy_plan()
	{	
		if (plan_backward == 0)
		{
			return;
		}

		cudaDeviceSynchronize();

		cufftDestroy(plan_forward);
		plan_backward = plan_forward = 0;
	}

	void FFT<dt_float64, edev_gpu>::create_plan_1d(const dt_int32& s0, dt_int32 n_thread)
	{
		destroy_plan();

		cufftPlan1d(&plan_forward, s0, CUFFT_Z2Z, 1);

		plan_backward = plan_forward;
	}

	void FFT<dt_float64, edev_gpu>::create_plan_1d_batch(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread)
	{
		destroy_plan();

		dt_int32 rank = 1;		// Dimensionality of the transform (1, 2, or 3). 
		dt_int32 n[] = {s0};		// 1d transforms of length s0*s1
		dt_int32 how_mas0 = s1;
		dt_int32 idist = s0, odist = s0;		// distance between two successive in elements in the least significant
		dt_int32 istride = 1, ostride = 1;		// distance between two elements in the same column
		dt_int32 *inembed = n, *onembed = n;		// Pointer of size rank that indicates the storage dimensions

		cufftPlanMany(&plan_forward, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, how_mas0);

		plan_backward = plan_forward;
	}

	void FFT<dt_float64, edev_gpu>::create_plan_2d(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread)
	{
		destroy_plan();

		// https:// docs.nvidia.com/cuda/cufft/index.html
		cufftPlan2d(&plan_forward, s1, s0, CUFFT_Z2Z);

		plan_backward = plan_forward;
	}

	void FFT<dt_float64, edev_gpu>::create_plan_2d_batch(const dt_int32& s0, const dt_int32& s1, const dt_int32& s2, dt_int32 n_thread)
	{
		destroy_plan();

		dt_int32 rank = 2;		// Dimensionality of the transform (1, 2, or 3). 
		dt_int32 n[] = {s1, s0};		// 2d transforms of length s0*s1
		dt_int32 how_mas0 = s2;
		dt_int32 idist = s0*s1, odist = s0*s1;		// distance between two successive in elements in the least significant
		dt_int32 istride = 1, ostride = 1;		// distance between two elements in the same column
		dt_int32 *inembed = n, *onembed = n;		// Pointer of size rank that indicates the storage dimensions

		cufftPlanMany(&plan_forward, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, how_mas0);

		plan_backward = plan_forward;
	}

	void FFT<dt_float64, edev_gpu>::forward(TVctr_c& mx)
	{
		forward(mx, mx);
	}

	void FFT<dt_float64, edev_gpu>::inverse(TVctr_c& mx)
	{
		inverse(mx, mx);
	}

	void FFT<dt_float64, edev_gpu>::forward(TVctr_c& mx_i, TVctr_c& mx_o)
	{
		auto pdata_i = mx_i.data_cast<cufftDoubleComplex*>();
		auto pdata_o = mx_o.data_cast<cufftDoubleComplex*>();

		cufftExecZ2Z(plan_forward, pdata_i, pdata_o, CUFFT_FORWARD);
	}

	void FFT<dt_float64, edev_gpu>::inverse(TVctr_c& mx_i, TVctr_c& mx_o)
	{
		auto pdata_i = mx_i.data_cast<cufftDoubleComplex*>();
		auto pdata_o = mx_o.data_cast<cufftDoubleComplex*>();

		cufftExecZ2Z(plan_backward, pdata_i, pdata_o, CUFFT_INVERSE);
	}
}