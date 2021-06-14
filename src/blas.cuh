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

#ifndef BLAS_H
	#define BLAS_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include "math.cuh"
	#include "types.cuh"
	#include "type_traits_gen.cuh"

	#include <cuda.h>
	#include <cuda_runtime.h>
	#include <cublas_v2.h>

	namespace blass
	{
		template <class T, mt::eDev Dev>
		struct GEAM;

		template <class T>
		struct GEAM<T, mt::edev_cpu>
		{
			public:
				using value_type = T;
				using TVctr_c = mt::Vctr<T, mt::edev_cpu>;

				const mt::eDev device;

				GEAM(): device(mt::edev_cpu) {}
		};

		template <class T>
		struct GEAM<T, mt::edev_gpu>
		{
			public:
				using value_type = T;
				using TVctr = mt::Vctr<T, mt::edev_gpu>;

				const mt::eDev device;

				GEAM(): device(mt::edev_gpu)
				{ 
					cublasCreate(&handle);
				}

				~GEAM()
				{ 
					cublasDestroy(handle);
				}

				void operator()(mt::eOP trsa, mt::eOP trsb, dt_int32 m, dt_int32 n, T alpha, TVctr& A, T beta, TVctr& B, TVctr& C)
				{
					cublasOperation_t transa = op_2_cbop(trsa);
					cublasOperation_t transb = op_2_cbop(trsb);

					geam(handle, transa, transb, m, n, alpha, A, m, beta, B, m, C, m);
				}

			private:

				cublasOperation_t op_2_cbop(mt::eOP op)
				{
					return static_cast<cublasOperation_t>(static_cast<dt_int32>(op)-1);
				}

				void geam(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb, 
				dt_int32 m, dt_int32 n, dt_float32 alpha, mt::Vctr<dt_float32, mt::edev_gpu>& A, dt_int32 lda, dt_float32 beta, 
				mt::Vctr<T, mt::edev_cpu>& B, dt_int32 ldb, mt::Vctr<dt_float32, mt::edev_gpu>& C, dt_int32 ldc)
				{
					dt_float32 *ralpha = reinterpret_cast<dt_float32*>(&alpha);
					dt_float32 *rbeta = reinterpret_cast<dt_float32*>(&beta);

					dt_float32 *rA = reinterpret_cast<dt_float32*>(A.data());
					dt_float32 *rB = reinterpret_cast<dt_float32*>(B.data());
					dt_float32 *rC = reinterpret_cast<dt_float32*>(C.data());

					auto result = cublasSgeam(handle, transa, transb, m, n, ralpha, rA, lda, rbeta, rB, ldb, rC, ldc);
				}

				void geam(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb, 
				dt_int32 m, dt_int32 n, dt_float64 alpha, mt::Vctr<dt_float64, mt::edev_gpu>& A, dt_int32 lda, dt_float64 beta, 
				mt::Vctr<dt_float64, mt::edev_gpu>& B, dt_int32 ldb, mt::Vctr<dt_float64, mt::edev_gpu>& C, dt_int32 ldc)
				{
					dt_float64 *ralpha = reinterpret_cast<dt_float64*>(&alpha);
					dt_float64 *rbeta = reinterpret_cast<dt_float64*>(&beta);

					dt_float64 *rA = reinterpret_cast<dt_float64*>(A.data());
					dt_float64 *rB = reinterpret_cast<dt_float64*>(B.data());
					dt_float64 *rC = reinterpret_cast<dt_float64*>(C.data());

					cublasDgeam(handle, transa, transb, m, n, ralpha, rA, lda, rbeta, rB, ldb, rC, ldc);
				}

				void geam(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb, 
				dt_int32 m, dt_int32 n, dt_cfloat32 alpha, mt::Vctr<dt_cfloat32, mt::edev_gpu>& A, dt_int32 lda, dt_cfloat32 beta, 
				mt::Vctr<dt_cfloat32, mt::edev_gpu>& B, dt_int32 ldb, mt::Vctr<dt_cfloat32, mt::edev_gpu>& C, dt_int32 ldc)
				{
					cuComplex *ralpha = reinterpret_cast<cuComplex*>(&alpha);
					cuComplex *rbeta = reinterpret_cast<cuComplex*>(&beta);

					cuComplex *rA = reinterpret_cast<cuComplex*>(A.data());
					cuComplex *rB = reinterpret_cast<cuComplex*>(B.data());
					cuComplex *rC = reinterpret_cast<cuComplex*>(C.data());

					cublasCgeam(handle, transa, transb, m, n, ralpha, rA, lda, rbeta, rB, ldb, rC, ldc);
				}

				void geam(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb, 
				dt_int32 m, dt_int32 n, dt_cfloat64 alpha, mt::Vctr<dt_cfloat64, mt::edev_gpu>& A, dt_int32 lda, dt_cfloat64 beta, 
				mt::Vctr<dt_cfloat64, mt::edev_gpu>& B, dt_int32 ldb, mt::Vctr<dt_cfloat64, mt::edev_gpu>& C, dt_int32 ldc)
				{
					cuDoubleComplex *ralpha = reinterpret_cast<cuDoubleComplex*>(&alpha);
					cuDoubleComplex *rbeta = reinterpret_cast<cuDoubleComplex*>(&beta);

					cuDoubleComplex *rA = reinterpret_cast<cuDoubleComplex*>(A.data());
					cuDoubleComplex *rB = reinterpret_cast<cuDoubleComplex*>(B.data());
					cuDoubleComplex *rC = reinterpret_cast<cuDoubleComplex*>(C.data());

					cublasZgeam(handle, transa, transb, m, n, ralpha, rA, lda, rbeta, rB, ldb, rC, ldc);
				}

				cublasHandle_t handle;
		};
	} // namespace blas

#endif