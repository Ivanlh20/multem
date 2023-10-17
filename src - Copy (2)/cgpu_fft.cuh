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

#ifndef CGPU_FFT_H
	#define CGPU_FFT_H

	#include <fftw3.h>

	#ifdef __CUDACC__
		#include <cuda.h>
		#include <cuda_runtime.h>
		#include <cufft.h>
	#endif

	#include "macros.cuh"
	#include "const_enum.cuh"
	#include "math.cuh"
	#include "cgpu_vctr.cuh"
	#include "cgpu_stream.cuh"

	/* cpu fourier transform */
	namespace mt
	{
		template <class T, eDev Dev> class FFT;

		template <>
		class FFT<dt_float32, edev_cpu>
		{
		public:
			using value_type = dt_float32;
			static const eDev device = edev_cpu;

			using TVctr_c = Vctr<complex<value_type>, edev_cpu>;

			FFT(): plan_forward(nullptr), plan_backward(nullptr)
			{ 
				fftwf_init_threads();
			}

			FFT(const dt_int32& s0, Stream_cpu* pstream = nullptr): FFT()
			{ 
				const dt_int32 n_thread = (fcn_is_single_thread(pstream))?1:pstream->size();

				create_plan_1d(s0, n_thread);
			}

			FFT(const dt_int32& s0, const dt_int32& s1, Stream_cpu* pstream = nullptr): FFT()
			{ 
				const dt_int32 n_thread = (fcn_is_single_thread(pstream))?1:pstream->size();

				create_plan_2d(s0, s1, n_thread);
			}

			~FFT()
			{
				destroy_plan();
			}

			void cleanup()
			{
				destroy_plan();

				fftwf_cleanup_threads();
			}

			void destroy_plan()
			{
				if (plan_backward == nullptr)
				{
					return;
				}

				fftwf_destroy_plan(plan_forward);
				fftwf_destroy_plan(plan_backward);

				plan_backward = plan_forward = nullptr;
			}

			void create_plan_1d(const dt_int32& s0, dt_int32 n_thread=1)
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

			void create_plan_1d_batch(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread=1)
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

			void create_plan_2d(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread=1)
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

			void create_plan_2d_batch(const dt_int32& s0, const dt_int32& s1, const dt_int32& s2, dt_int32 n_thread=1)
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

			void forward(TVctr_c& mx)
			{
				auto pdata = mx.data_cast<fftwf_complex*>();

				fftwf_execute_dft(plan_forward, pdata, pdata);
			}

			void inverse(TVctr_c& mx)
			{
				auto pdata = mx.data_cast<fftwf_complex*>();

				fftwf_execute_dft(plan_backward, pdata, pdata);
			}

			void forward(TVctr_c& mx_i, TVctr_c& mx_o)
			{
				if (mx_i.data() != mx_o.data())
				{
					mx_o.assign(mx_i.begin(), mx_i.end());
				}

				forward(mx_o);
			}

			void inverse(TVctr_c& mx_i, TVctr_c& mx_o)
			{
				if (mx_i.data() != mx_o.data())
				{
					mx_o.assign(mx_i.begin(), mx_i.end());
				}

				inverse(mx_o);
			}
		private:
			fftwf_plan plan_forward;
			fftwf_plan plan_backward;
		};

		template <>
		class FFT<dt_float64, edev_cpu>
		{
		public:
			using value_type = dt_float64;
			static const eDev device = edev_cpu;

			using TVctr_c = Vctr<complex<value_type>, edev_cpu>;

			FFT(): plan_forward(nullptr), plan_backward(nullptr)
			{ 
				fftw_init_threads();
			}

			FFT(const dt_int32& s0, Stream_cpu* pstream = nullptr): FFT()
			{ 
				const dt_int32 n_thread = (fcn_is_single_thread(pstream))?1:pstream->size();

				create_plan_1d(s0, n_thread);
			}

			FFT(const dt_int32& s0, const dt_int32& s1, Stream_cpu* pstream = nullptr): FFT()
			{ 
				const dt_int32 n_thread = (fcn_is_single_thread(pstream))?1:pstream->size();

				create_plan_2d(s0, s1, n_thread);
			}

			~FFT()
			{
				destroy_plan();
			}

			void cleanup()
			{
				destroy_plan();

				fftw_cleanup_threads();
			}

			void destroy_plan()
			{
				if (plan_backward == nullptr)
				{
					return;
				}

				fftw_destroy_plan(plan_forward);
				fftw_destroy_plan(plan_backward);

				plan_backward = plan_forward = nullptr;
			}

			void create_plan_1d(const dt_int32& s0, dt_int32 n_thread=1)
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

			void create_plan_1d_batch(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread=1)
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

			void create_plan_2d(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread=1)
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

			void create_plan_2d_batch(const dt_int32& s0, const dt_int32& s1, const dt_int32& s2, dt_int32 n_thread=1)
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

			void forward(TVctr_c& mx)
			{
				auto pdata = mx.data_cast<fftw_complex*>();

				fftw_execute_dft(plan_forward, pdata, pdata);
			}

			void inverse(TVctr_c& mx)
			{
				auto pdata = mx.data_cast<fftw_complex*>();

				fftw_execute_dft(plan_backward, pdata, pdata);
			}

			void forward(TVctr_c& mx_i, TVctr_c& mx_o)
			{
				if (mx_i.data() != mx_o.data())
				{
					mx_o.assign(mx_i.begin(), mx_i.end());
				}

				forward(mx_o);
			}

			void inverse(TVctr_c& mx_i, TVctr_c& mx_o)
			{
				if (mx_i.data() != mx_o.data())
				{
					mx_o.assign(mx_i.begin(), mx_i.end());
				}

				inverse(mx_o);
			}
		private:
			fftw_plan plan_forward;
			fftw_plan plan_backward;
		};
	}
	
	/* gpu fourier transform */
	namespace mt
	{
	#ifdef __CUDACC__
		template <>
		class FFT<dt_float32, edev_gpu>
		{
		public:
			using value_type = dt_float32;
			static const eDev device = edev_gpu;

			using TVctr_c = Vctr<complex<value_type>, edev_gpu>;

			FFT(): plan_forward(0), plan_backward(0) {}

			FFT(const dt_int32& s0, Stream_gpu* pstream = nullptr): FFT()
			{ 
				create_plan_1d(s0);
			}

			FFT(const dt_int32& s0, const dt_int32& s1, Stream_gpu* pstream = nullptr): FFT()
			{ 
				create_plan_2d(s0, s1);
			}

			~FFT()
			{
				destroy_plan();
			}

			void cleanup()
			{
				destroy_plan();
			}

			void destroy_plan()
			{	
				if (plan_backward == 0)
				{
					return;
				}

				cudaDeviceSynchronize();

				cufftDestroy(plan_forward);
				plan_backward = plan_forward = 0;
			}

			void create_plan_1d(const dt_int32& s0, dt_int32 n_thread=1)
			{
				destroy_plan();

				cufftPlan1d(&plan_forward, s0, CUFFT_C2C, 1);

				plan_backward = plan_forward;
			}

			void create_plan_1d_batch(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread=1)
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

			void create_plan_2d(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread=1)
			{
				destroy_plan();

				// https:// docs.nvidia.com/cuda/cufft/index.html
				cufftPlan2d(&plan_forward, s1, s0, CUFFT_C2C);

				plan_backward = plan_forward;
			}

			void create_plan_2d_batch(const dt_int32& s0, const dt_int32& s1, const dt_int32& nz, dt_int32 n_thread=1)
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

			void set_stream(cudaStream_t &stream)
			{
				cufftSetStream(plan_forward, stream);
			}

			void forward(TVctr_c& mx)
			{
				forward(mx, mx);
			}

			void inverse(TVctr_c& mx)
			{
				inverse(mx, mx);
			}

			void forward(TVctr_c& mx_i, TVctr_c& mx_o)
			{
				auto pdata_i = mx_i.data_cast<cufftComplex*>();
				auto pdata_o = mx_o.data_cast<cufftComplex*>();

				cufftExecC2C(plan_forward, pdata_i, pdata_o, CUFFT_FORWARD);
			}

			void inverse(TVctr_c& mx_i, TVctr_c& mx_o)
			{
				auto pdata_i = mx_i.data_cast<cufftComplex*>();
				auto pdata_o = mx_o.data_cast<cufftComplex*>();

				cufftExecC2C(plan_backward, pdata_i, pdata_o, CUFFT_INVERSE);
			}
		private:
			cufftHandle plan_forward;
			cufftHandle plan_backward;
		};

		template <>
		class FFT<dt_float64, edev_gpu>
		{
		public:
			using value_type = dt_float64;
			static const eDev device = edev_gpu;

			using TVctr_c = Vctr<complex<value_type>, edev_gpu>;

			FFT(): plan_forward(0), plan_backward(0) {}

			FFT(const dt_int32& s0, Stream_gpu* pstream = nullptr): FFT()
			{ 
				create_plan_1d(s0);
			}

			FFT(const dt_int32& s0, const dt_int32& s1, Stream_gpu* pstream = nullptr): FFT()
			{ 
				create_plan_2d(s0, s1);
			}

			~FFT()
			{
				destroy_plan();
			}

			void cleanup()
			{
				destroy_plan();
			}

			void destroy_plan()
			{	
				if (plan_backward == 0)
				{
					return;
				}

				cudaDeviceSynchronize();

				cufftDestroy(plan_forward);
				plan_backward = plan_forward = 0;
			}

			void create_plan_1d(const dt_int32& s0, dt_int32 n_thread=1)
			{
				destroy_plan();

				cufftPlan1d(&plan_forward, s0, CUFFT_Z2Z, 1);

				plan_backward = plan_forward;
			}

			void create_plan_1d_batch(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread=1)
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

			void create_plan_2d(const dt_int32& s0, const dt_int32& s1, dt_int32 n_thread=1)
			{
				destroy_plan();

				// https:// docs.nvidia.com/cuda/cufft/index.html
				cufftPlan2d(&plan_forward, s1, s0, CUFFT_Z2Z);

				plan_backward = plan_forward;
			}

			void create_plan_2d_batch(const dt_int32& s0, const dt_int32& s1, const dt_int32& s2, dt_int32 n_thread=1)
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

			void forward(TVctr_c& mx)
			{
				forward(mx, mx);
			}

			void inverse(TVctr_c& mx)
			{
				inverse(mx, mx);
			}

			void forward(TVctr_c& mx_i, TVctr_c& mx_o)
			{
				auto pdata_i = mx_i.data_cast<cufftDoubleComplex*>();
				auto pdata_o = mx_o.data_cast<cufftDoubleComplex*>();

				cufftExecZ2Z(plan_forward, pdata_i, pdata_o, CUFFT_FORWARD);
			}

			void inverse(TVctr_c& mx_i, TVctr_c& mx_o)
			{
				auto pdata_i = mx_i.data_cast<cufftDoubleComplex*>();
				auto pdata_o = mx_o.data_cast<cufftDoubleComplex*>();

				cufftExecZ2Z(plan_backward, pdata_i, pdata_o, CUFFT_INVERSE);
			}
		private:
			cufftHandle plan_forward;
			cufftHandle plan_backward;
		};
	#endif
	}

	/* derive classes */
	namespace mt
	{
		/***************************************************************************************/
		/************************************ alias ********************************************/
		/***************************************************************************************/
		template <class T>
		using FFT_cpu = FFT<T, edev_cpu>;

		template <class T>
		using FFT_gpu = FFT<T, edev_gpu>;

	}

#endif