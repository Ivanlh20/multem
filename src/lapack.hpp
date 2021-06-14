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

#ifndef LAPACK_H
	#define LAPACK_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#if defined(_WIN32) && defined(MATLAB_BLAS_LAPACK)
		#define ssyevr_ ssyevr
		#define dsyevr_ dsyevr
		#define sgesv_ sgesv
		#define dgesv_ dgesv
		#define sgels_ sgels
		#define dgels_ dgels
		#define sgelsd_ sgelsd
		#define dgelsd_ dgelsd
	#endif

	#include "macros.cuh"
	#include "math.cuh"
	#include "type_traits_gen.cuh"
	#include "types.cuh"
	#include "cgpu_vctr.cuh"

	namespace lapack
	{
		extern "C" void ssyevr_(
			const char* jobz, 
			const char* range, 
			const char* uplo, 
			const MTL_BL_INT* n, 
			dt_float32* a, 
			const MTL_BL_INT* lda, 
			const dt_float32* vl, 
			const dt_float32* vu, 
			const MTL_BL_INT* il, 
			const MTL_BL_INT* iu, 
			const dt_float32* abstol, 
			MTL_BL_INT* m, 
			dt_float32* w, 
			dt_float32* z, 
			const MTL_BL_INT* ldz, 
			MTL_BL_INT* isuppz, 
			dt_float32* work, 
			const MTL_BL_INT* lwork, 
			MTL_BL_INT* iwork, 
			const MTL_BL_INT* liwork, 
			MTL_BL_INT* info
		);

		extern "C" void dsyevr_(
			const char* jobz, 
			const char* range, 
			const char* uplo, 
			const MTL_BL_INT* n, 
			dt_float64* a, 
			const MTL_BL_INT* lda, 
			const dt_float64* vl, 
			const dt_float64* vu, 
			const MTL_BL_INT* il, 
			const MTL_BL_INT* iu, 
			const dt_float64* abstol, 
			MTL_BL_INT* m, 
			dt_float64* w, 
			dt_float64* z, 
			const MTL_BL_INT* ldz, 
			MTL_BL_INT* isuppz, 
			dt_float64* work, 
			const MTL_BL_INT* lwork, 
			MTL_BL_INT* iwork, 
			const MTL_BL_INT* liwork, 
			MTL_BL_INT* info
		);

		extern "C" void sgesv_(
			const MTL_BL_INT* n, 
			const MTL_BL_INT* nrhs, 
			dt_float32* a, 
			const MTL_BL_INT* lda, 
			MTL_BL_INT* ipiv, 
			dt_float32* b, 
			const MTL_BL_INT* ldb, 
			MTL_BL_INT* info
		);

		extern "C" void dgesv_(
			const MTL_BL_INT* n, 
			const MTL_BL_INT* nrhs, 
			dt_float64* a, 
			const MTL_BL_INT* lda, 
			MTL_BL_INT* ipiv, 
			dt_float64* b, 
			const MTL_BL_INT* ldb, 
			MTL_BL_INT* info
		);

		extern "C" void sgels_(
			const char	*trans, 
			const MTL_BL_INT* m, 
			const MTL_BL_INT* n, 
			const MTL_BL_INT* nrhs, 
			dt_float32 *a, 
			const MTL_BL_INT* lda, 
			dt_float32 *b, 
			const MTL_BL_INT* ldb, 
			dt_float32 *work, 
			const MTL_BL_INT* lwork, 
			MTL_BL_INT* info
		);

		extern "C" void dgels_(
			const char	*trans, 
			const MTL_BL_INT* m, 
			const MTL_BL_INT* n, 
			const MTL_BL_INT* nrhs, 
			dt_float64 *a, 
			const MTL_BL_INT* lda, 
			dt_float64 *b, 
			const MTL_BL_INT* ldb, 
			dt_float64 *work, 
			const MTL_BL_INT* lwork, 
			MTL_BL_INT* info
		);

		extern "C" void sgelsd_(
			const MTL_BL_INT* m, 
			const MTL_BL_INT* n, 
			const MTL_BL_INT* nrhs, 
			const dt_float32 *a, 
			const MTL_BL_INT* lda, 
			dt_float32 *b, 
			const MTL_BL_INT* ldb, 
			dt_float32 *s, 
			const dt_float32 *rcond, 
			MTL_BL_INT* rank, 
			dt_float32 *work, 
			const MTL_BL_INT* lwork, 
			MTL_BL_INT* iwork, 
			MTL_BL_INT* info
		);

		extern "C" void dgelsd_(
			const MTL_BL_INT* m, 
			const MTL_BL_INT* n, 
			const MTL_BL_INT* nrhs, 
			const dt_float64 *a, 
			const MTL_BL_INT* lda, 
			dt_float64 *b, 
			const MTL_BL_INT* ldb, 
			dt_float64 *s, 
			const dt_float64 *rcond, 
			MTL_BL_INT* rank, 
			dt_float64 *work, 
			const MTL_BL_INT* lwork, 
			MTL_BL_INT* iwork, 
			MTL_BL_INT* info
		);

		/*	computes eigenvalue of a real symmetric matrix A */
		template <class T>
		struct EIGR
		{
			public:
				using value_type = T;

				const mt::eDev device;

				EIGR(): jobz('N'), range('I'), uplo('U'), n(0), lda(0), vl(0), vu(0), 
				il(1), iu(1), abstol(0), m(0), ldz(0), lwork(0), liwork(0), info(0), device(mt::edev_cpu) {}

				EIGR(MTL_BL_INT n_i): jobz('N'), range('I'), uplo('U'), n(n_i), lda(0), vl(0), vu(0), 
				il(1), iu(1), abstol(0), m(0), ldz(0), lwork(0), liwork(0), info(0), device(mt::edev_cpu)
				{
					init(n_i);
				}

				// get eigenvalues
				T operator()(T* AS, MTL_BL_INT idx=1)
				{
					il = max<MTL_BL_INT>(1, idx);
					iu = min<MTL_BL_INT>(n, idx);

					std::copy(AS, AS+n*n, MS.begin());
					syevr_c(jobz, range, uplo, n, MS.data(), lda, vl, vu, il, iu, abstol, m, w.data(), z.data(), ldz, isuppz.data(), work.data(), lwork, iwork.data(), liwork, info);

					return w[0];
				}

				// get min eigenvalue
				T eig_min(T* AS)
				{
					il = 1;
					iu = 1;

					std::copy(AS, AS+n*n, MS.begin());
					syevr_c(jobz, range, uplo, n, MS.data(), lda, vl, vu, il, iu, abstol, m, w.data(), z.data(), ldz, isuppz.data(), work.data(), lwork, iwork.data(), liwork, info);

					return w[0];
				}

				// get max eigenvalue
				T eig_max(T* AS)
				{
					il = n;
					iu = n;

					std::copy(AS, AS+n*n, MS.begin());
					syevr_c(jobz, range, uplo, n, MS.data(), lda, vl, vu, il, iu, abstol, m, w.data(), z.data(), ldz, isuppz.data(), work.data(), lwork, iwork.data(), liwork, info);

					return w[0];
				}

				// get max eigenvalue
				void eig_min_max(T* AS, T& eig_min, T& eig_max)
				{
					il = 1;
					iu = n;

					std::copy(AS, AS+n*n, MS.begin());
					syevr_c(jobz, range, uplo, n, MS.data(), lda, vl, vu, il, iu, abstol, m, w.data(), z.data(), ldz, isuppz.data(), work.data(), lwork, iwork.data(), liwork, info);

					eig_min = w[0];
					eig_max = w[n-1];
				}

				void init(MTL_BL_INT n_i)
				{
					n = n_i;
					lda = n;
					vl = 0.0;
					vu = 0.0;
					il = 1;
					iu = 1;
					abstol = -1;
					ldz = n;

					// (n*n-n)/2+n
					MS.resize(n*n);
					w.resize(n);
					z.resize(n*n);
					isuppz.resize(2*n);

					// query optimal size of work array
					lwork = -1;
					liwork = -1;
					T work_query = 0;
					MTL_BL_INT iwork_query = 0;
					syevr_c(jobz, range, uplo, n, MS.data(), lda, vl, vu, il, iu, abstol, m, w.data(), z.data(), ldz, isuppz.data(), &work_query, lwork, &iwork_query, liwork, info);

					// set arrays
					lwork = max<MTL_BL_INT>(1, MTL_BL_INT(work_query));
					work.resize(lwork);

					liwork = max<MTL_BL_INT>(1, iwork_query);
					iwork.resize(liwork);
				}

			private:
				void syevr_c(const char &jobz, const char &range, const char &uplo, const MTL_BL_INT& n, dt_float32 *a, 
				const MTL_BL_INT& lda, const dt_float32 &vl, const dt_float32 &vu, const MTL_BL_INT& il, const MTL_BL_INT& iu, 
				const dt_float32 &abstol, MTL_BL_INT& m, dt_float32 *w, dt_float32 *z, const MTL_BL_INT& ldz, MTL_BL_INT* isuppz, dt_float32 *work, 
				const MTL_BL_INT& lwork, MTL_BL_INT* iwork, const MTL_BL_INT& liwork, MTL_BL_INT& info)
				{
					ssyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
				}

				void syevr_c(const char &jobz, const char &range, const char &uplo, const MTL_BL_INT& n, dt_float64 *a, 
				const MTL_BL_INT& lda, const dt_float64 &vl, const dt_float64 &vu, const MTL_BL_INT& il, const MTL_BL_INT& iu, 
				const dt_float64 &abstol, MTL_BL_INT& m, dt_float64 *w, dt_float64 *z, const MTL_BL_INT& ldz, MTL_BL_INT* isuppz, dt_float64 *work, 
				const MTL_BL_INT& lwork, MTL_BL_INT* iwork, const MTL_BL_INT& liwork, MTL_BL_INT& info)
				{
					dsyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
				}

				char jobz;
				char range;
				char uplo;
				mt::Vctr_cpu<T> MS;
				MTL_BL_INT n;
				MTL_BL_INT lda;
				T vl;
				T vu;
				MTL_BL_INT il;
				MTL_BL_INT iu;
				T abstol;
				MTL_BL_INT m;
				mt::Vctr_cpu<T> w;
				mt::Vctr_cpu<T> z;
				MTL_BL_INT ldz;
				mt::Vctr_cpu<MTL_BL_INT> isuppz;
				mt::Vctr_cpu<T> work;
				MTL_BL_INT lwork;
				mt::Vctr_cpu<MTL_BL_INT> iwork;
				MTL_BL_INT liwork;
				MTL_BL_INT info;
		};

		/* Fast least square fitting */
		template <class T>
		struct LSF_1
		{
			public:
				using value_type = T;

				const mt::eDev device;

				LSF_1(): device(mt::edev_cpu) {}

				void operator()(MTL_BL_INT A_rows, MTL_BL_INT A_cols, T* A, T* b, T* x)
				{
					using TVctr = mt::Vctr_cpu<T>;

					// get ATA
					TVctr M;
					M.reserve(A_cols*A_cols);

					for(auto ir=0; ir<A_cols; ir++)
					{
						for(auto ic=0; ic<A_cols; ic++)
						{
							T v = 0;
							for(auto k=0;k<A_rows;k++)
							{
								v += A[ir*A_rows+k]*A[ic*A_rows+k];
							}
							M.push_back(v);
						}
					}

					// get y
					TVctr y;
					y.reserve(A_cols);

					for(auto ic=0; ic<A_cols; ic++)
					{
						T v = 0;
						for(auto ir=0; ir<A_rows; ir++)
						{
							v += A[ic*A_rows+ir]*b[ir];
						}
						y.push_back(v);
					}

					MTL_BL_INT n = A_cols;
					MTL_BL_INT nrhs = 1;
					MTL_BL_INT lda = n;
					MTL_BL_INT ldb = n;
					MTL_BL_INT info = 0;

					mt::Vctr_cpu<MTL_BL_INT> ipiv(n);

					gesv_c(n, nrhs, M.data(), lda, ipiv.data(), y.data(), ldb, info);
					std::copy(y.begin(), y.end(), x);
				}

				void operator()(MTL_BL_INT A_rows, MTL_BL_INT A_cols, T* A, T* b, T* x, T lambda, T* D, T* G)
				{
					using TVctr = mt::Vctr_cpu<T>;

					// get ATA
					TVctr M;
					M.reserve(A_cols*A_cols);

					for(auto ir=0; ir<A_cols; ir++)
					{
						for(auto ic=0; ic<A_cols; ic++)
						{
							T v = 0;
							for(auto k=0;k<A_rows;k++)
							{
								v += A[ir*A_rows+k]*A[ic*A_rows+k];
							}
							T d = 0;
							if (ir == ic)
							{
								d = (v<1e-7)?lambda:lambda*v;
								D[ir] = d;
							}
							M.push_back(v+d);
						}
					}

					// get y
					TVctr y;
					y.reserve(A_cols);

					for(auto ic=0; ic<A_cols; ic++)
					{
						T v = 0;
						for(auto ir=0; ir<A_rows; ir++)
						{
							v += A[ic*A_rows+ir]*b[ir];
						}
						y.push_back(v);
						G[ic] = v;
					}

					MTL_BL_INT n = A_cols;
					MTL_BL_INT nrhs = 1;
					MTL_BL_INT lda = n;
					MTL_BL_INT ldb = n;
					MTL_BL_INT info = 0;

					mt::Vctr_cpu<MTL_BL_INT> ipiv(n);

					gesv_c(n, nrhs, M.data(), lda, ipiv.data(), y.data(), ldb, info);
					std::copy(y.begin(), y.end(), x);
				}

			private:
				void gesv_c(const MTL_BL_INT& n, const MTL_BL_INT& nrhs, dt_float32 *a, const MTL_BL_INT& lda, 
				MTL_BL_INT* ipiv, dt_float32 *b, const MTL_BL_INT& ldb, MTL_BL_INT& info)
				{
					sgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
				}

				void gesv_c(const MTL_BL_INT& n, const MTL_BL_INT& nrhs, dt_float64 *a, const MTL_BL_INT& lda, 
				MTL_BL_INT* ipiv, dt_float64 *b, const MTL_BL_INT& ldb, MTL_BL_INT& info)
				{
					dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
				}

				mt::Vctr_cpu<MTL_BL_INT> ipiv;
				mt::Vctr_cpu<T> M;
				mt::Vctr_cpu<T> y;
		};

		/* RQ minimum-norm solution */
		template <class T>
		struct MNS_QR
		{
			public:
				using value_type = T;

				const mt::eDev device;

				MNS_QR(): trans('N'), info(0), lwork(-1), m(0), n(0), device(mt::edev_cpu) {}

				MNS_QR(MTL_BL_INT A_rows, MTL_BL_INT A_cols): trans('N'), info(0), lwork(-1), m(A_rows), n(A_cols), device(mt::edev_cpu) 
				{
					init(A_rows, A_cols);
				}

				void operator()(T* A, T* b, T* x)
				{
					MTL_BL_INT nrhs = 1;
					MTL_BL_INT lda = n;
					MTL_BL_INT ldb = n;

					std::copy(A, A+m*n, M.begin());
					std::copy(b, b+n, y.begin());

					// perform minimum-norm solution
					gels_c(trans, m, n, nrhs, M.data(), lda, y.data(), ldb, work.data(), lwork, info);
					std::copy(y.begin(), y.end(), x);
				}

				T operator()(T* A, T* b, T* x, T lambda)
				{
					MTL_BL_INT nrhs = 1;
					MTL_BL_INT lda = n;
					MTL_BL_INT ldb = n;

					std::copy(A, A+m*n, M.begin());
					std::copy(b, b+n, y.begin());

					// get lambda and get maximum gradient
					T G_max = 0;
					for(auto ic=0; ic<n; ic++)
					{
						M[ic*m+ic] *= 1+lambda;
						G_max = ::fmax(G_max, fabs(y[ic]));
					}

					// perform minimum-norm solution
					gels_c(trans, m, n, nrhs, M.data(), lda, y.data(), ldb, work.data(), lwork, info);
					std::copy(y.begin(), y.end(), x);

					return G_max;
				}

				void init(MTL_BL_INT A_rows, MTL_BL_INT A_cols)
				{
					m = A_rows;
					n = A_cols;

					MTL_BL_INT nrhs = 1;
					MTL_BL_INT lda = n;
					MTL_BL_INT ldb = n;
					MTL_BL_INT info = 0;

					M.resize(m*n);
					y.resize(n);

					// query optimal size of work array
					lwork = -1;
					T work_query = 0;
					gels_c(trans, m, n, nrhs, M.data(), lda, y.data(), ldb, &work_query, lwork, info);

					// set arrays
					lwork = max<MTL_BL_INT>(1, MTL_BL_INT(work_query));
					work.resize(lwork);
				}

			private:
				void gels_c(const char &trans, const MTL_BL_INT& m, const MTL_BL_INT& n, const MTL_BL_INT& nrhs, dt_float32 *a, 
				const MTL_BL_INT& lda, dt_float32 *b, const MTL_BL_INT& ldb, dt_float32 *work, const MTL_BL_INT& lwork, MTL_BL_INT& info)
				{
					sgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
				}

				void gels_c(const char &trans, const MTL_BL_INT& m, const MTL_BL_INT& n, const MTL_BL_INT& nrhs, dt_float64 *a, 
				const MTL_BL_INT& lda, dt_float64 *b, const MTL_BL_INT& ldb, dt_float64 *work, const MTL_BL_INT& lwork, MTL_BL_INT& info)
				{
					dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
				}

				const char trans;
				MTL_BL_INT info;
				MTL_BL_INT lwork;
				MTL_BL_INT m;
				MTL_BL_INT n;
				mt::Vctr_cpu<T> work;

				mt::Vctr_cpu<T> M;
				mt::Vctr_cpu<T> y;
		};

		/* SVD minimum-norm solution */
		template <class T>
		struct MNS_SVD
		{
			public:
				using value_type = T;

				const mt::eDev device;

				MNS_SVD(): rcond(-1), 	info(0), rank(0), lwork(-1), m(0), n(0), device(mt::edev_cpu) {}
				MNS_SVD(MTL_BL_INT A_rows, MTL_BL_INT A_cols): rcond(-1), info(0), rank(0), lwork(-1), m(A_rows), n(A_cols), device(mt::edev_cpu) 
				{
					init(A_rows, A_cols);
				}

				void operator()(T* A, T* b, T* x)
				{
					std::copy(A, A+m*n, M.begin());
					std::copy(b, b+n, y.begin());

					// perform minimum-norm solution
					gelsd_c(m, n, nrhs, M.data(), lda, y.data(), ldb, S.data(), rcond, rank, work.data(), lwork, iwork.data(), info);
					std::copy(y.begin(), y.end(), x);
				}

				void init(MTL_BL_INT A_rows, MTL_BL_INT A_cols)
				{
					m = A_rows;
					n = A_cols;

					nrhs = 1;
					lda = n;
					ldb = n;

					S.resize(n);
					M.resize(m*n);
					y.resize(n);

					// query optimal size of work array
					lwork = -1;
					T work_query = 0;
					MTL_BL_INT iwork_query = 0;
					gelsd_c(m, n, nrhs, M.data(), lda, y.data(), ldb, S.data(), rcond, rank, &work_query, lwork, &iwork_query, info);

					// set arrays
					lwork = max<MTL_BL_INT>(1000, MTL_BL_INT(work_query));
					work.resize(lwork);

					liwork = max<MTL_BL_INT>(1, iwork_query);
					iwork.resize(liwork);
				}

			private:
				void gelsd_c(const MTL_BL_INT& m, const MTL_BL_INT& n, const MTL_BL_INT& nrhs, dt_float32 *a, const MTL_BL_INT& lda, 
				dt_float32 *b, const MTL_BL_INT& ldb, dt_float32 *s, const dt_float32 &rcond, MTL_BL_INT& rank, dt_float32 *work, 
				const MTL_BL_INT& lwork, MTL_BL_INT* iwork, MTL_BL_INT& info)
				{
					sgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);
				}

				void gelsd_c(const MTL_BL_INT& m, const MTL_BL_INT& n, const MTL_BL_INT& nrhs, dt_float64 *a, const MTL_BL_INT& lda, 
				dt_float64 *b, const MTL_BL_INT& ldb, dt_float64 *s, const dt_float64 &rcond, MTL_BL_INT& rank, dt_float64 *work, 
				const MTL_BL_INT& lwork, MTL_BL_INT* iwork, MTL_BL_INT& info)
				{
					dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);
				}

				MTL_BL_INT m;
				MTL_BL_INT n;
				MTL_BL_INT nrhs;
				MTL_BL_INT lda;
				MTL_BL_INT ldb;
				T rcond;
				MTL_BL_INT info;
				MTL_BL_INT rank;
				MTL_BL_INT liwork;
				MTL_BL_INT lwork;
				mt::Vctr_cpu<MTL_BL_INT> iwork;
				mt::Vctr_cpu<T> work;
				mt::Vctr_cpu<T> S;

				mt::Vctr_cpu<T> M;
				mt::Vctr_cpu<T> y;
		};

		/*	The program computes the solution to the system of linear
			equations with a square matrix A and multiple right-hand 
			sides B, where A is the coefficient matrix and b is the 
			right-hand side matrix:
		*/

		template <class T>
		struct GESV
		{
			public:
				using value_type = T;

				const mt::eDev device;

				GESV(): device(mt::edev_cpu) {}

				void operator()(MTL_BL_INT A_n, T* A, MTL_BL_INT b_cols, T* b, T* x)
				{
					MTL_BL_INT n = A_n;
					MTL_BL_INT nrhs = b_cols;
					MTL_BL_INT lda = n;
					MTL_BL_INT ldb = n;
					MTL_BL_INT info = 0;

					std::copy(b, b+n, x);
					ipiv.resize(n);

					gesv_c(n, nrhs, A, lda, ipiv.data(), x, ldb, info);
				}

			private:
				void gesv_c(const MTL_BL_INT& n, const MTL_BL_INT& nrhs, dt_float32 *a, const MTL_BL_INT& lda, 
				MTL_BL_INT* ipiv, dt_float32 *b, const MTL_BL_INT& ldb, MTL_BL_INT& info)
				{
					sgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
				}

				void gesv_c(const MTL_BL_INT& n, const MTL_BL_INT& nrhs, dt_float64 *a, const MTL_BL_INT& lda, 
				MTL_BL_INT* ipiv, dt_float64 *b, const MTL_BL_INT& ldb, MTL_BL_INT& info)
				{
					dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
				}

				mt::Vctr_cpu<MTL_BL_INT> ipiv;
		};

		/*	The program solves overdetermined or underdetermined real 
			linear systems involving an M-by-N matrix A, or its trs, 
			using a QR or LQ factorization of A. It is assumed that A has full rank.
		*/
		template <class T>
		struct GELS
		{
			public:
				using value_type = T;

				const mt::eDev device;

				GELS(): device(mt::edev_cpu) {}

				void operator()(MTL_BL_INT A_rows, MTL_BL_INT A_cols, T* A, MTL_BL_INT b_cols, T* b, T* x)
				{
					const char trans = 'N';
					MTL_BL_INT m = A_rows;
					MTL_BL_INT n = A_cols;
					MTL_BL_INT nrhs = b_cols;
					MTL_BL_INT lda = A_rows;
					MTL_BL_INT ldb = max(b_cols, max(m, n));
					MTL_BL_INT min_mn = min(m, n);
					T rcond = -1;
					MTL_BL_INT info = 0;

					mt::Vctr<T, mt::edev_cpu> bv(b, b+m*b_cols);

					// query optimal size of work array
					MTL_BL_INT lwork = -1;
					T work_query= 0;
					gels_c(trans, m, n, nrhs, A, lda, bv.data(), ldb, &work_query, lwork, info);

					// set arrays
					lwork = max<MTL_BL_INT>(1, MTL_BL_INT(work_query));
					work.resize(lwork);

					// perform minimum-norm solution
					gels_c(trans, m, n, nrhs, A, lda, bv.data(), ldb, work.data(), lwork, info);

					for(auto ix = 0; ix<nrhs; ix++)
					{
						std::copy(bv.begin()+ix*m, bv.begin()+ix*m+n, x+ix*n);
					}
				}

			private:
				void gels_c(const char &trans, const MTL_BL_INT& m, const MTL_BL_INT& n, const MTL_BL_INT& nrhs, dt_float32 *a, 
				const MTL_BL_INT& lda, dt_float32 *b, const MTL_BL_INT& ldb, dt_float32 *work, const MTL_BL_INT& lwork, MTL_BL_INT& info)
				{
					sgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
				}

				void gels_c(const char &trans, const MTL_BL_INT& m, const MTL_BL_INT& n, const MTL_BL_INT& nrhs, dt_float64 *a, 
				const MTL_BL_INT& lda, dt_float64 *b, const MTL_BL_INT& ldb, dt_float64 *work, const MTL_BL_INT& lwork, MTL_BL_INT& info)
				{
					dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
				}

				mt::Vctr<T, mt::edev_cpu> work;
		};

		/*	The program solves overdetermined or underdetermined real 
			linear systems involving an M-by-N matrix A, or its trs, 
			using the singular value decomposition (SVD) of A. A is an m-by-n matrix which may be rank-deficient.
		*/
		template <class T>
		struct GELSD
		{
			public:
				using value_type = T;

				const mt::eDev device;

				GELSD(): device(mt::edev_cpu) {}

				void operator()(MTL_BL_INT A_rows, MTL_BL_INT A_cols, T* A, MTL_BL_INT b_cols, T* b, T* x)
				{
					MTL_BL_INT m = A_rows;
					MTL_BL_INT n = A_cols;
					MTL_BL_INT nrhs = b_cols;
					MTL_BL_INT lda = A_rows;
					MTL_BL_INT ldb = max(b_cols, max(m, n));
					MTL_BL_INT min_mn = min(m, n);
					MTL_BL_INT rank = 0;
					T rcond = -1;
					MTL_BL_INT info = 0;

					mt::Vctr<T, mt::edev_cpu> bv(b, b+m*b_cols);
					S.resize(min_mn);

					// query optimal size of work array
					MTL_BL_INT lwork = -1;
					T work_query= 0;
					MTL_BL_INT iwork_query= 0;
					gelsd_c(m, n, nrhs, A, lda, bv.data(), ldb, S.data(), rcond, rank, &work_query, lwork, &iwork_query, info);

					// set arrays
					lwork = max<MTL_BL_INT>(1, MTL_BL_INT(work_query));
					work.resize(lwork);

					MTL_BL_INT liwork = max<MTL_BL_INT>(1, iwork_query);
					iwork.resize(liwork);

					// perform minimum-norm solution
					gelsd_c(m, n, nrhs, A, lda, bv.data(), ldb, S.data(), rcond, rank, work.data(), lwork, iwork.data(), info);

					for(auto ix= 0; ix<nrhs; ix++)
					{
						std::copy(bv.begin()+ix*m, bv.begin()+ix*m+n, x+ix*n);
					}
				}

			private:
				void gelsd_c(const MTL_BL_INT& m, const MTL_BL_INT& n, const MTL_BL_INT& nrhs, dt_float32 *a, const MTL_BL_INT& lda, 
				dt_float32 *b, const MTL_BL_INT& ldb, dt_float32 *s, const dt_float32 &rcond, MTL_BL_INT& rank, dt_float32 *work, 
				const MTL_BL_INT& lwork, MTL_BL_INT* iwork, MTL_BL_INT& info)
				{
					sgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);
				}

				void gelsd_c(const MTL_BL_INT& m, const MTL_BL_INT& n, const MTL_BL_INT& nrhs, dt_float64 *a, const MTL_BL_INT& lda, 
				dt_float64 *b, const MTL_BL_INT& ldb, dt_float64 *s, const dt_float64 &rcond, MTL_BL_INT& rank, dt_float64 *work, 
				const MTL_BL_INT& lwork, MTL_BL_INT* iwork, MTL_BL_INT& info)
				{
					dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);
				}

				mt::Vctr<MTL_BL_INT, mt::edev_cpu> iwork;
				mt::Vctr<T, mt::edev_cpu> work;
				mt::Vctr<T, mt::edev_cpu> S;
		};
} // namespace lapack

#endif