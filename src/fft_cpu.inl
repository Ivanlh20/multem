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

#include "fft_cpu.h"

/* cpu fourier transform - dt_float32 */
namespace mt
{
	FFT<dt_float32, edev_cpu>::FFT(): plan_forward(nullptr), plan_backward(nullptr)
	{ 
		fftwf_init_threads();
	}

	FFT<dt_float32, edev_cpu>::FFT(const dt_int32& s0, Stream_cpu* pstream): FFT()
	{ 
		const dt_int32 n_thread = (fcn_is_single_thread(pstream))?1:pstream->size();

		create_plan_1d(s0, n_thread);
	}

	FFT<dt_float32, edev_cpu>::FFT(const dt_int32& s0, const dt_int32& s1, Stream_cpu* pstream): FFT()
	{ 
		const dt_int32 n_thread = (fcn_is_single_thread(pstream))?1:pstream->size();

		create_plan_2d(s0, s1, n_thread);
	}

	FFT<dt_float32, edev_cpu>::~FFT()
	{
		destroy_plan();
	}

	void FFT<dt_float32, edev_cpu>::cleanup()
	{
		destroy_plan();

		fftwf_cleanup_threads();
	}

	void FFT<dt_float32, edev_cpu>::destroy_plan()
	{
		if (plan_backward == nullptr)
		{
			return;
		}

		fftwf_destroy_plan(plan_forward);
		fftwf_destroy_plan(plan_backward);

		plan_backward = plan_forward = nullptr;
	}

	void FFT<dt_float32, edev_cpu>::create_plan_1d(const dt_int32& s0, dt_int32 n_thread)
	{
		destroy_plan();

		fftwf_import_wisdom_from_filename("fftwf_1d.wisdom");

		fftwf_plan_with_nthreads(n_thread);

		TVctr_c M(s0);

		auto pdata = M.data_cast<fftwf_complex*>();

		plan_forward = fftwf_plan_dft_1d(s0, pdata, pdata, FFTW_FORWARD, FFTW_MEASURE);
		plan_backward = fftwf_plan_dft_1d(s0, pdata, pdata, FFTW_BACKWARD, FFTW_MEASURE);

		fftwf_export_wisdom_to_filename("fftwf_1d.wisdom");
	}

	void FFT<dt_float32, edev_cpu>::create_plan_1d_batch(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread)
	{
		destroy_plan();

		fftwf_import_wisdom_from_filename("fftwf_1d_batch.wisdom");

		fftwf_plan_with_nthreads(n_thread);

		TVctr_c M(s0*s1);

		auto pdata = M.data_cast<fftwf_complex*>();

		dt_int32 rank = 1;							// Dimensionality of the transform (1, 2, or 3). 
		dt_int32 n[] = {s0};						// 1d transforms of length s0*s1
		dt_int32 how_many = s1;
		dt_int32 idist = s0, odist = s0;			// distance between two successive in elements in the least significant
		dt_int32 istride = 1, ostride = 1;			// distance between two elements in the same column
		dt_int32 *inembed = n, *onembed = n;		// Pointer of size rank that indicates the storage dimensions

		plan_forward = fftwf_plan_many_dft(rank, n, how_many, pdata, inembed, istride, idist, pdata, onembed, ostride, odist, FFTW_FORWARD, FFTW_MEASURE);
		plan_backward = fftwf_plan_many_dft(rank, n, how_many, pdata, inembed, istride, idist, pdata, onembed, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);

		fftwf_export_wisdom_to_filename("fftwf_1d_batch.wisdom");
	}

	void FFT<dt_float32, edev_cpu>::create_plan_2d(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread)
	{
		destroy_plan();

		fftwf_import_wisdom_from_filename("fftwf_2d.wisdom");

		fftwf_plan_with_nthreads(n_thread);

		TVctr_c M(s0*s1);

		auto pdata = M.data_cast<fftwf_complex*>();

		plan_forward = fftwf_plan_dft_2d(s1, s0, pdata, pdata, FFTW_FORWARD, FFTW_MEASURE);
		plan_backward = fftwf_plan_dft_2d(s1, s0, pdata, pdata, FFTW_BACKWARD, FFTW_MEASURE);

		fftwf_export_wisdom_to_filename("fftwf_2d.wisdom");
	}

	void FFT<dt_float32, edev_cpu>::create_plan_2d_batch(const dt_int32& s0, const dt_int32& s1, const dt_int32& s2, dt_int32 n_thread)
	{
		destroy_plan();

		fftwf_import_wisdom_from_filename("fftwf_2d_batch.wisdom");

		fftwf_plan_with_nthreads(n_thread);

		TVctr_c M(s0*s1*s2);

		auto pdata = M.data_cast<fftwf_complex*>();

		dt_int32 rank = 2;								// Dimensionality of the transform (1, 2, or 3). 
		dt_int32 n[] = {s1, s0};						// 2d transforms of length s0*s1
		dt_int32 how_many = s2;
		dt_int32 idist = s0*s1, odist = s0*s1;			// distance between two successive in elements in the least significant
		dt_int32 istride = 1, ostride = 1;				// distance between two elements in the same column
		dt_int32 *inembed = n, *onembed = n;			// Pointer of size rank that indicates the storage dimensions

		plan_forward = fftwf_plan_many_dft(rank, n, how_many, pdata, inembed, istride, idist, pdata, onembed, ostride, odist, FFTW_FORWARD, FFTW_MEASURE);
		plan_backward = fftwf_plan_many_dft(rank, n, how_many, pdata, inembed, istride, idist, pdata, onembed, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);

		fftwf_export_wisdom_to_filename("fftwf_2d_batch.wisdom");
	}

	void FFT<dt_float32, edev_cpu>::forward(TVctr_c& mx)
	{
		auto pdata = mx.data_cast<fftwf_complex*>();

		fftwf_execute_dft(plan_forward, pdata, pdata);
	}

	void FFT<dt_float32, edev_cpu>::inverse(TVctr_c& mx)
	{
		auto pdata = mx.data_cast<fftwf_complex*>();

		fftwf_execute_dft(plan_backward, pdata, pdata);
	}

	void FFT<dt_float32, edev_cpu>::forward(TVctr_c& mx_i, TVctr_c& mx_o)
	{
		if (mx_i.data() != mx_o.data())
		{
			mx_o.assign(mx_i.begin(), mx_i.end());
		}

		forward(mx_o);
	}

	void FFT<dt_float32, edev_cpu>::inverse(TVctr_c& mx_i, TVctr_c& mx_o)
	{
		if (mx_i.data() != mx_o.data())
		{
			mx_o.assign(mx_i.begin(), mx_i.end());
		}

		inverse(mx_o);
	}
}

/* cpu fourier transform - dt_float64 */
namespace mt
{
	FFT<dt_float64, edev_cpu>::FFT(): plan_forward(nullptr), plan_backward(nullptr)
	{ 
		fftw_init_threads();
	}

	FFT<dt_float64, edev_cpu>::FFT(const dt_int32& s0, Stream_cpu* pstream): FFT()
	{ 
		const dt_int32 n_thread = (fcn_is_single_thread(pstream))?1:pstream->size();

		create_plan_1d(s0, n_thread);
	}

	FFT<dt_float64, edev_cpu>::FFT(const dt_int32& s0, const dt_int32& s1, Stream_cpu* pstream): FFT()
	{ 
		const dt_int32 n_thread = (fcn_is_single_thread(pstream))?1:pstream->size();

		create_plan_2d(s0, s1, n_thread);
	}

	FFT<dt_float64, edev_cpu>::~FFT()
	{
		destroy_plan();
	}

	void FFT<dt_float64, edev_cpu>::cleanup()
	{
		destroy_plan();

		fftw_cleanup_threads();
	}

	void FFT<dt_float64, edev_cpu>::destroy_plan()
	{
		if (plan_backward == nullptr)
		{
			return;
		}

		fftw_destroy_plan(plan_forward);
		fftw_destroy_plan(plan_backward);

		plan_backward = plan_forward = nullptr;
	}

	void FFT<dt_float64, edev_cpu>::create_plan_1d(const dt_int32& s0, dt_int32 n_thread)
	{
		destroy_plan();

		fftw_import_wisdom_from_filename("fftw_1d.wisdom");

		fftw_plan_with_nthreads(n_thread);

		TVctr_c M(s0);

		auto pdata = M.data_cast<fftw_complex*>();

		plan_forward = fftw_plan_dft_1d(s0, pdata, pdata, FFTW_FORWARD, FFTW_MEASURE);
		plan_backward = fftw_plan_dft_1d(s0, pdata, pdata, FFTW_BACKWARD, FFTW_MEASURE);

		fftw_export_wisdom_to_filename("fftw_1d.wisdom");
	}

	void FFT<dt_float64, edev_cpu>::create_plan_1d_batch(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread)
	{
		destroy_plan();

		fftw_import_wisdom_from_filename("fftw_1d_batch.wisdom");

		fftw_plan_with_nthreads(n_thread);

		TVctr_c M(s0*s1);

		auto pdata = M.data_cast<fftw_complex*>();

		dt_int32 rank = 1;							// Dimensionality of the transform (1, 2, or 3). 
		dt_int32 n[] = {s0};						// 1d transforms of length s0*s1
		dt_int32 how_many = s1;
		dt_int32 idist = s0, odist = s0;			// distance between two successive in elements in the least significant
		dt_int32 istride = 1, ostride = 1;			// distance between two elements in the same column
		dt_int32 *inembed = n, *onembed = n;		// Pointer of size rank that indicates the storage dimensions

		plan_forward = fftw_plan_many_dft(rank, n, how_many, pdata, inembed, istride, idist, pdata, onembed, ostride, odist, FFTW_FORWARD, FFTW_MEASURE);
		plan_backward = fftw_plan_many_dft(rank, n, how_many, pdata, inembed, istride, idist, pdata, onembed, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);

		fftw_export_wisdom_to_filename("fftw_1d_batch.wisdom");
	}

	void FFT<dt_float64, edev_cpu>::create_plan_2d(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread)
	{
		destroy_plan();

		fftw_import_wisdom_from_filename("fftw_2d.wisdom");

		fftw_plan_with_nthreads(n_thread);

		TVctr_c M(s0*s1);

		auto pdata = M.data_cast<fftw_complex*>();

		plan_forward = fftw_plan_dft_2d(s1, s0, pdata, pdata, FFTW_FORWARD, FFTW_MEASURE);
		plan_backward = fftw_plan_dft_2d(s1, s0, pdata, pdata, FFTW_BACKWARD, FFTW_MEASURE);

		fftw_export_wisdom_to_filename("fftw_2d.wisdom");
	}

	void FFT<dt_float64, edev_cpu>::create_plan_2d_batch(const dt_int32& s0, const dt_int32& s1, const dt_int32& s2, dt_int32 n_thread)
	{
		destroy_plan();

		fftwf_import_wisdom_from_filename("fftw_2d_batch.wisdom");

		fftwf_plan_with_nthreads(n_thread);

		TVctr_c M(s0*s1*s2);

		auto pdata = M.data_cast<fftw_complex*>();

		dt_int32 rank = 2;							// Dimensionality of the transform (1, 2, or 3). 
		dt_int32 n[] = {s1, s0};					// 2d transforms of length s0*s1
		dt_int32 how_many = s2;
		dt_int32 idist = s0*s1, odist = s0*s1;		// distance between two successive in elements in the least significant
		dt_int32 istride = 1, ostride = 1;			// distance between two elements in the same column
		dt_int32 *inembed = n, *onembed = n;		// Pointer of size rank that indicates the storage dimensions

		plan_forward = fftw_plan_many_dft(rank, n, how_many, pdata, inembed, istride, idist, pdata, onembed, ostride, odist, FFTW_FORWARD, FFTW_MEASURE);
		plan_backward = fftw_plan_many_dft(rank, n, how_many, pdata, inembed, istride, idist, pdata, onembed, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);

		fftwf_export_wisdom_to_filename("fftwf_2d_batch.wisdom");
	}

	void FFT<dt_float64, edev_cpu>::forward(TVctr_c& mx)
	{
		auto pdata = mx.data_cast<fftw_complex*>();

		fftw_execute_dft(plan_forward, pdata, pdata);
	}

	void FFT<dt_float64, edev_cpu>::inverse(TVctr_c& mx)
	{
		auto pdata = mx.data_cast<fftw_complex*>();

		fftw_execute_dft(plan_backward, pdata, pdata);
	}

	void FFT<dt_float64, edev_cpu>::forward(TVctr_c& mx_i, TVctr_c& mx_o)
	{
		if (mx_i.data() != mx_o.data())
		{
			mx_o.assign(mx_i.begin(), mx_i.end());
		}

		forward(mx_o);
	}

	void FFT<dt_float64, edev_cpu>::inverse(TVctr_c& mx_i, TVctr_c& mx_o)
	{
		if (mx_i.data() != mx_o.data())
		{
			mx_o.assign(mx_i.begin(), mx_i.end());
		}

		inverse(mx_o);
	}
}