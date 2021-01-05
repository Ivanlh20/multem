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

#ifndef BLAS_H
#define BLAS_H

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"

#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>

namespace bs
{
	template <class T, mt::eDevice dev>
	struct GEAM;

	template <class T>
	struct GEAM<T, mt::e_host>
	{
		public:
			using value_type = T;
			using TVector_c = mt::Vector<T, mt::e_host>;

			const mt::eDevice device;

			GEAM(): device(mt::e_host) {}
	};

	template <class T>
	struct GEAM<T, mt::e_device>
	{
		public:
			using value_type = T;
			using TVector = mt::Vector<T, mt::e_device>;

			const mt::eDevice device;

			GEAM(): device(mt::e_device)
			{ 
				cublasCreate(&handle);
			}

			~GEAM()
			{ 
				cublasDestroy(handle); 
			}

			void operator()(mt::eOP trsa, mt::eOP trsb, int m, int n, T alpha, TVector &A, T beta, TVector &B, TVector &C)
			{
				cublasOperation_t transa = op_2_cbop(trsa);
				cublasOperation_t transb = op_2_cbop(trsb);

				geam(handle, transa, transb, m, n, alpha, A, m, beta, B, m, C, m);
			}

		private:

			cublasOperation_t op_2_cbop(mt::eOP op)
			{
				return static_cast<cublasOperation_t>(static_cast<int>(op)-1);
			}

			void geam(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb, 
			int m, int n, float alpha, device_vector<float> &A, int lda, float beta, 
			device_vector<float> &B, int ldb, device_vector<float> &C, int ldc)
			{
				float *ralpha = reinterpret_cast<float*>(&alpha);
				float *rbeta = reinterpret_cast<float*>(&beta);

				float *rA = reinterpret_cast<float*>(raw_pointer_cast(A.data()));
				float *rB = reinterpret_cast<float*>(raw_pointer_cast(B.data()));
				float *rC = reinterpret_cast<float*>(raw_pointer_cast(C.data()));

				auto result = cublasSgeam(handle, transa, transb, m, n, ralpha, rA, lda, rbeta, rB, ldb, rC, ldc);
			}

			void geam(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb, 
			int m, int n, double alpha, device_vector<double> &A, int lda, double beta, 
			device_vector<double> &B, int ldb, device_vector<double> &C, int ldc)
			{
				double *ralpha = reinterpret_cast<double*>(&alpha);
				double *rbeta = reinterpret_cast<double*>(&beta);

				double *rA = reinterpret_cast<double*>(raw_pointer_cast(A.data()));
				double *rB = reinterpret_cast<double*>(raw_pointer_cast(B.data()));
				double *rC = reinterpret_cast<double*>(raw_pointer_cast(C.data()));

				cublasDgeam(handle, transa, transb, m, n, ralpha, rA, lda, rbeta, rB, ldb, rC, ldc);
			}

			void geam(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb, 
			int m, int n, complex<float> alpha, device_vector<complex<float>> &A, int lda, complex<float> beta, 
			device_vector<complex<float>> &B, int ldb, device_vector<complex<float>> &C, int ldc)
			{
				cuComplex *ralpha = reinterpret_cast<cuComplex*>(&alpha);
				cuComplex *rbeta = reinterpret_cast<cuComplex*>(&beta);

				cuComplex *rA = reinterpret_cast<cuComplex*>(raw_pointer_cast(A.data()));
				cuComplex *rB = reinterpret_cast<cuComplex*>(raw_pointer_cast(B.data()));
				cuComplex *rC = reinterpret_cast<cuComplex*>(raw_pointer_cast(C.data()));

				cublasCgeam(handle, transa, transb, m, n, ralpha, rA, lda, rbeta, rB, ldb, rC, ldc);
			}

			void geam(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb, 
			int m, int n, complex<double> alpha, device_vector<complex<double>> &A, int lda, complex<double> beta, 
			device_vector<complex<double>> &B, int ldb, device_vector<complex<double>> &C, int ldc)
			{
				cuDoubleComplex *ralpha = reinterpret_cast<cuDoubleComplex*>(&alpha);
				cuDoubleComplex *rbeta = reinterpret_cast<cuDoubleComplex*>(&beta);

				cuDoubleComplex *rA = reinterpret_cast<cuDoubleComplex*>(raw_pointer_cast(A.data()));
				cuDoubleComplex *rB = reinterpret_cast<cuDoubleComplex*>(raw_pointer_cast(B.data()));
				cuDoubleComplex *rC = reinterpret_cast<cuDoubleComplex*>(raw_pointer_cast(C.data()));

				cublasZgeam(handle, transa, transb, m, n, ralpha, rA, lda, rbeta, rB, ldb, rC, ldc);
			}

			cublasHandle_t handle;
	};
} // namespace blas

#endif