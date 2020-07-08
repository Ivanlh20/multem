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

#ifndef FFT_H
#define FFT_H

#include "math.cuh"
#include "types.cuh"
#include <fftw3.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>

//#ifdef __CUDACC__
//	#pragma message("Cuda FFT_H")
//#else
//	#pragma message("nonCuda FFT_H")
//#endif

namespace mt
{
	template <class T, eDevice dev>
	struct FFT;

	template <>
	struct FFT<float, e_host>
	{
		public:
			using value_type = float;
			using TVector_c = Vector<complex<float>, e_host>;

			static const eDevice device = e_host;

			FFT(): plan_forward(nullptr), plan_backward(nullptr){ fftwf_init_threads(); }

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
				if(plan_backward == nullptr)
				{
					return;
				}

				fftwf_destroy_plan(plan_forward);
				fftwf_destroy_plan(plan_backward);

				plan_backward = plan_forward = nullptr;
			}

			void create_plan_1d(const int &nx, int nThread=1)
			{
				destroy_plan();

				fftwf_import_wisdom_from_filename("fftwf_1d.wisdom");

				fftwf_plan_with_nthreads(nThread);

				TVector_c M(nx);

				fftwf_complex *V = reinterpret_cast<fftwf_complex*>(M.data());

				plan_forward = fftwf_plan_dft_1d(nx, V, V, FFTW_FORWARD, FFTW_MEASURE);
				plan_backward = fftwf_plan_dft_1d(nx, V, V, FFTW_BACKWARD, FFTW_MEASURE);

				fftwf_export_wisdom_to_filename("fftwf_1d.wisdom");
			}

			void create_plan_1d_batch(const int &ny, const int &nx, int nThread=1)
			{
				destroy_plan();

				fftwf_import_wisdom_from_filename("fftwf_1d_batch.wisdom");

				fftwf_plan_with_nthreads(nThread);

				TVector_c M(nx*ny);

				auto V = reinterpret_cast<fftwf_complex*>(M.data());

				int rank = 1;						// Dimensionality of the transform (1, 2, or 3). 
				int n[] = {ny};						// 1d transforms of length nx*ny
				int how_many = nx;
				int idist = ny, odist = ny;			// distance between two successive input elements in the least significant
				int istride = 1, ostride = 1;		// distance between two elements in the same column
				int *inembed = n, *onembed = n;		// Pointer of size rank that indicates the storage dimensions

				plan_forward = fftwf_plan_many_dft(rank, n, how_many, V, inembed, istride, idist, V, onembed, ostride, odist, FFTW_FORWARD, FFTW_MEASURE);
				plan_backward = fftwf_plan_many_dft(rank, n, how_many, V, inembed, istride, idist, V, onembed, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);

				fftwf_export_wisdom_to_filename("fftwf_1d_batch.wisdom");
			}

			void create_plan_2d(const int &ny, const int &nx, int nThread)
			{
				destroy_plan();

				fftwf_import_wisdom_from_filename("fftwf_2d.wisdom");

				fftwf_plan_with_nthreads(nThread);

				TVector_c M(nx*ny);

				auto V = reinterpret_cast<fftwf_complex*>(M.data());

				plan_forward = fftwf_plan_dft_2d(nx, ny, V, V, FFTW_FORWARD, FFTW_MEASURE);
				plan_backward = fftwf_plan_dft_2d(nx, ny, V, V, FFTW_BACKWARD, FFTW_MEASURE);

				fftwf_export_wisdom_to_filename("fftwf_2d.wisdom");
			}

			void create_plan_2d_batch(const int &ny, const int &nx, const int &nz, int nThread)
			{
				destroy_plan();

				fftwf_import_wisdom_from_filename("fftwf_2d_batch.wisdom");

				fftwf_plan_with_nthreads(nThread);

				TVector_c M(nx*ny*nz);

				auto V = reinterpret_cast<fftwf_complex*>(M.data());

				int rank = 2;						// Dimensionality of the transform (1, 2, or 3). 
				int n[] = {nx, ny};					// 2d transforms of length nx*ny
				int how_many = nz;
				int idist = nx*ny, odist = nx*ny;	// distance between two successive input elements in the least significant
				int istride = 1, ostride = 1;		// distance between two elements in the same column
				int *inembed = n, *onembed = n;		// Pointer of size rank that indicates the storage dimensions

				plan_forward = fftwf_plan_many_dft(rank, n, how_many, V, inembed, istride, idist, V, onembed, ostride, odist, FFTW_FORWARD, FFTW_MEASURE);
				plan_backward = fftwf_plan_many_dft(rank, n, how_many, V, inembed, istride, idist, V, onembed, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);

				fftwf_export_wisdom_to_filename("fftwf_2d_batch.wisdom");
			}

			template <class TVector>
			void forward(TVector &M_io)
			{
				fftwf_complex *V_io = reinterpret_cast<fftwf_complex*>(M_io.data());
				fftwf_execute_dft(plan_forward, V_io, V_io);
			}

			template <class TVector>
			void inverse(TVector &M_io)
			{
				fftwf_complex *V_io = reinterpret_cast<fftwf_complex*>(M_io.data());
				fftwf_execute_dft(plan_backward, V_io, V_io);
			}

			template <class TVector>
			void forward(TVector &M_i, TVector &M_o)
			{
				if (M_i.data() != M_o.data())
				{
					M_o.assign(M_i.begin(), M_i.end());
				}

				forward(M_o);
			}

			template <class TVector>
			void inverse(TVector &M_i, TVector &M_o)
			{
				if (M_i.data() != M_o.data())
				{
					M_o.assign(M_i.begin(), M_i.end());
				}

				inverse(M_o);
			}
		private:
			fftwf_plan plan_forward;
			fftwf_plan plan_backward;
	};

	template <>
	struct FFT<double, e_host>
	{
		public:
			using value_type = double;
			using TVector_c = Vector<complex<double>, e_host>;

			static const eDevice device = e_host;

			FFT(): plan_forward(nullptr), plan_backward(nullptr){ fftw_init_threads(); }

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
				if(plan_backward == nullptr)
				{
					return;
				}

				fftw_destroy_plan(plan_forward);
				fftw_destroy_plan(plan_backward);

				plan_backward = plan_forward = nullptr;
			}

			void create_plan_1d(const int &nx, int nThread=1)
			{
				destroy_plan();

				fftw_import_wisdom_from_filename("fftw_1d.wisdom");

				fftw_plan_with_nthreads(nThread);

				TVector_c M(nx);

				auto V = reinterpret_cast<fftw_complex*>(M.data());

				plan_forward = fftw_plan_dft_1d(nx, V, V, FFTW_FORWARD, FFTW_MEASURE);
				plan_backward = fftw_plan_dft_1d(nx, V, V, FFTW_BACKWARD, FFTW_MEASURE);

				fftw_export_wisdom_to_filename("fftw_1d.wisdom");
			}

			void create_plan_1d_batch(const int &ny, const int &nx, int nThread=1)
			{
				destroy_plan();

				fftw_import_wisdom_from_filename("fftw_1d_batch.wisdom");

				fftw_plan_with_nthreads(nThread);

				TVector_c M(nx*ny);

				auto V = reinterpret_cast<fftw_complex*>(M.data());

				int rank = 1;						// Dimensionality of the transform (1, 2, or 3). 
				int n[] = {ny};						// 1d transforms of length nx*ny
				int how_many = nx;
				int idist = ny, odist = ny;			// distance between two successive input elements in the least significant
				int istride = 1, ostride = 1;		// distance between two elements in the same column
				int *inembed = n, *onembed = n;		// Pointer of size rank that indicates the storage dimensions

				plan_forward = fftw_plan_many_dft(rank, n, how_many, V, inembed, istride, idist, V, onembed, ostride, odist, FFTW_FORWARD, FFTW_MEASURE);
				plan_backward = fftw_plan_many_dft(rank, n, how_many, V, inembed, istride, idist, V, onembed, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);

				fftw_export_wisdom_to_filename("fftw_1d_batch.wisdom");
			}

			void create_plan_2d(const int &ny, const int &nx, int nThread)
			{
				destroy_plan();

				fftw_import_wisdom_from_filename("fftw_2d.wisdom");

				fftw_plan_with_nthreads(nThread);

				TVector_c M(nx*ny);

				auto V = reinterpret_cast<fftw_complex*>(M.data());

				plan_forward = fftw_plan_dft_2d(nx, ny, V, V, FFTW_FORWARD, FFTW_MEASURE);
				plan_backward = fftw_plan_dft_2d(nx, ny, V, V, FFTW_BACKWARD, FFTW_MEASURE);

				fftw_export_wisdom_to_filename("fftw_2d.wisdom");
			}

			void create_plan_2d_batch(const int &ny, const int &nx, const int &nz, int nThread)
			{
				destroy_plan();

				fftwf_import_wisdom_from_filename("fftw_2d_batch.wisdom");

				fftwf_plan_with_nthreads(nThread);

				TVector_c M(nx*ny*nz);

				auto V = reinterpret_cast<fftw_complex*>(M.data());

				int rank = 2;						// Dimensionality of the transform (1, 2, or 3). 
				int n[] = {nx, ny};					// 2d transforms of length nx*ny
				int how_many = nz;
				int idist = nx*ny, odist = nx*ny;	// distance between two successive input elements in the least significant
				int istride = 1, ostride = 1;		// distance between two elements in the same column
				int *inembed = n, *onembed = n;		// Pointer of size rank that indicates the storage dimensions

				plan_forward = fftw_plan_many_dft(rank, n, how_many, V, inembed, istride, idist, V, onembed, ostride, odist, FFTW_FORWARD, FFTW_MEASURE);
				plan_backward = fftw_plan_many_dft(rank, n, how_many, V, inembed, istride, idist, V, onembed, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);

				fftwf_export_wisdom_to_filename("fftwf_2d_batch.wisdom");
			}

			template <class TVector>
			void forward(TVector &M_io)
			{
				auto V_io = reinterpret_cast<fftw_complex*>(M_io.data());
				fftw_execute_dft(plan_forward, V_io, V_io);
			}

			template <class TVector>
			void inverse(TVector &M_io)
			{
				auto V_io = reinterpret_cast<fftw_complex*>(M_io.data());
				fftw_execute_dft(plan_backward, V_io, V_io);
			}

			template <class TVector>
			void forward(TVector &M_i, TVector &M_o)
			{
				if (M_i.data() != M_o.data())
				{
					M_o.assign(M_i.begin(), M_i.end());
				}

				forward(M_o);
			}

			template <class TVector>
			void inverse(TVector &M_i, TVector &M_o)
			{
				if (M_i.data() != M_o.data())
				{
					M_o.assign(M_i.begin(), M_i.end());
				}

				inverse(M_o);
			}
		private:
			fftw_plan plan_forward;
			fftw_plan plan_backward;
	};

	template <>
	struct FFT<float, e_device>
	{
		public:
			using value_type = float;
			using TVector_c = Vector<complex<float>, e_device>;

			static const eDevice device = e_device;

			FFT(): plan_forward(0), plan_backward(0){}

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
				if(plan_backward == 0)
				{
					return;
				}

				cudaDeviceSynchronize();

				cufftDestroy(plan_forward);
				plan_backward = plan_forward = 0;
			}

			void create_plan_1d(const int &nx, int nThread=1)
			{
				destroy_plan();

				cufftPlan1d(&plan_forward, nx, CUFFT_C2C, 1);

				plan_backward = plan_forward;
			}

			void create_plan_1d_batch(const int &ny, const int &nx, int nThread=1)
			{
				destroy_plan();

				int rank = 1;						// Dimensionality of the transform (1, 2, or 3). 
				int n[] = {ny};						// 1d transforms of length nx*ny
				int how_many = nx;
				int idist = ny, odist = ny;			// distance between two successive input elements in the least significant
				int istride = 1, ostride = 1;		// distance between two elements in the same column
				int *inembed = n, *onembed = n;		// Pointer of size rank that indicates the storage dimensions

				cufftPlanMany(&plan_forward, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, how_many);

				plan_backward = plan_forward;
			}

			void create_plan_2d(const int &ny, const int &nx, int nThread)
			{
				destroy_plan();

				cufftPlan2d(&plan_forward, nx, ny, CUFFT_C2C);

				plan_backward = plan_forward;
			}

			void create_plan_2d_batch(const int &ny, const int &nx, const int &nz, int nThread)
			{
				destroy_plan();

				int rank = 2;						// Dimensionality of the transform (1, 2, or 3). 
				int n[] = {nx, ny};					// 2d transforms of length nx*ny
				int how_many = nz;
				int idist = nx*ny, odist = nx*ny;	// distance between two successive input elements in the least significant
				int istride = 1, ostride = 1;		// distance between two elements in the same column
				int *inembed = n, *onembed = n;		// Pointer of size rank that indicates the storage dimensions

				cufftPlanMany(&plan_forward, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, how_many);

				plan_backward = plan_forward;
			}

			void setstream(cudaStream_t &stream)
			{
				cufftSetStream(plan_forward, stream);
			}

			void forward(TVector_c &M_io)
			{
				forward(M_io, M_io);
			}

			void inverse(TVector_c &M_io)
			{
				inverse(M_io, M_io);
			}

			void forward(TVector_c &M_i, TVector_c &M_o)
			{
				auto V_i = reinterpret_cast<cufftComplex*>(raw_pointer_cast(M_i.data()));
				auto V_o = reinterpret_cast<cufftComplex*>(raw_pointer_cast(M_o.data()));
				cufftExecC2C(plan_forward, V_i, V_o, CUFFT_FORWARD);
			}

			void inverse(TVector_c &M_i, TVector_c &M_o)
			{
				auto V_i = reinterpret_cast<cufftComplex*>(raw_pointer_cast(M_i.data()));
				auto V_o = reinterpret_cast<cufftComplex*>(raw_pointer_cast(M_o.data()));
				cufftExecC2C(plan_backward, V_i, V_o, CUFFT_INVERSE);
			}
		private:
			cufftHandle plan_forward;
			cufftHandle plan_backward;
	};

	template <>
	struct FFT<double, e_device>
	{
		public:
			using value_type = double;
			using TVector_c = Vector<complex<double>, e_device>;

			static const eDevice device = e_device;

			FFT(): plan_forward(0), plan_backward(0){}

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
				if(plan_backward == 0)
				{
					return;
				}

				cudaDeviceSynchronize();

				cufftDestroy(plan_forward);
				plan_backward = plan_forward = 0;
			}

			void create_plan_1d(const int &nx, int nThread=1)
			{
				destroy_plan();

				cufftPlan1d(&plan_forward, nx, CUFFT_Z2Z, 1);

				plan_backward = plan_forward;
			}

			void create_plan_1d_batch(const int &ny, const int &nx, int nThread=1)
			{
				destroy_plan();

				int rank = 1;						// Dimensionality of the transform (1, 2, or 3). 
				int n[] = {ny};						// 1d transforms of length nx*ny
				int how_many = nx;
				int idist = ny, odist = ny;			// distance between two successive input elements in the least significant
				int istride = 1, ostride = 1;		// distance between two elements in the same column
				int *inembed = n, *onembed = n;		// Pointer of size rank that indicates the storage dimensions

				cufftPlanMany(&plan_forward, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, how_many);

				plan_backward = plan_forward;
			}

			void create_plan_2d(const int &ny, const int &nx, int nThread)
			{
				destroy_plan();

				cufftPlan2d(&plan_forward, nx, ny, CUFFT_Z2Z);

				plan_backward = plan_forward;
			}

			void create_plan_2d_batch(const int &ny, const int &nx, const int &nz, int nThread)
			{
				destroy_plan();

				int rank = 2;						// Dimensionality of the transform (1, 2, or 3). 
				int n[] = {nx, ny};					// 2d transforms of length nx*ny
				int how_many = nz;
				int idist = nx*ny, odist = nx*ny;	// distance between two successive input elements in the least significant
				int istride = 1, ostride = 1;		// distance between two elements in the same column
				int *inembed = n, *onembed = n;		// Pointer of size rank that indicates the storage dimensions

				cufftPlanMany(&plan_forward, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, how_many);

				plan_backward = plan_forward;
			}

			void forward(TVector_c &M_io)
			{
				forward(M_io, M_io);
			}

			void inverse(TVector_c &M_io)
			{
				inverse(M_io, M_io);
			}

			void forward(TVector_c &M_i, TVector_c &M_o)
			{
				auto V_i = reinterpret_cast<cufftDoubleComplex*>(raw_pointer_cast(M_i.data()));
				auto V_o = reinterpret_cast<cufftDoubleComplex*>(raw_pointer_cast(M_o.data()));
				cufftExecZ2Z(plan_forward, V_i, V_o, CUFFT_FORWARD);
			}

			void inverse(TVector_c &M_i, TVector_c &M_o)
			{
				auto V_i = reinterpret_cast<cufftDoubleComplex*>(raw_pointer_cast(M_i.data()));
				auto V_o = reinterpret_cast<cufftDoubleComplex*>(raw_pointer_cast(M_o.data()));
				cufftExecZ2Z(plan_backward, V_i, V_o, CUFFT_INVERSE);
			}
		private:
			cufftHandle plan_forward;
			cufftHandle plan_backward;
	};
} // namespace mt

#endif