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

#ifndef LAPACK_H
#define LAPACK_H

#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"

namespace lapack
{
	extern "C" void sgesv_(
		const int* n, 
		const int* nrhs, 
		float* a, 
		const int* lda, 
		int* ipiv, 
 float* b, 
		const int* ldb, 
		int* info
	);

	extern "C" void dgesv_(
		const int* n, 
		const int* nrhs, 
		double* a, 
		const int* lda, 
		int* ipiv, 
 double* b, 
		const int* ldb, 
		int* info
	);

	extern "C" void sgels_(
		const char	*trans, 
		const int *m, 
		const int *n, 
		const int *nrhs, 
		float *a, 
		const int *lda, 
		float *b, 
		const int *ldb, 
		float *work, 
		const int *lwork, 
		int *info
	);

	extern "C" void dgels_(
		const char	*trans, 
		const int *m, 
		const int *n, 
		const int *nrhs, 
		double *a, 
		const int *lda, 
		double *b, 
		const int *ldb, 
		double *work, 
		const int *lwork, 
		int *info
	);

	extern "C" void sgelsd_(
		const int *m, 
		const int *n, 
		const int *nrhs, 
		const float *a, 
		const int *lda, 
		float *b, 
		const int *ldb, 
		float *s, 
		const float *rcond, 
		int *rank, 
		float *work, 
		const int *lwork, 
		int *iwork, 
		int *info
	);

	extern "C" void dgelsd_(
		const int *m, 
		const int *n, 
		const int *nrhs, 
		const double *a, 
		const int *lda, 
		double *b, 
		const int *ldb, 
		double *s, 
		const double *rcond, 
		int *rank, 
		double *work, 
		const int *lwork, 
		int *iwork, 
		int *info
	);
	/*	Fast least square fitting	*/
	template <class T>
	struct FLSF
	{
		public:
			using value_type = T;

			const mt::eDevice device;

			FLSF(): device(mt::e_host) {}

			void operator()(int A_rows, int A_cols, T *A, T *b, T *x)
			{
				using TVector = vector<T>;

				// get ATA
				TVector M;
				M.reserve(A_cols*A_cols);

				for(auto ir=0; ir<A_cols; ir++)
				{
					for(auto ic=0; ic<A_cols; ic++)
					{
						T v = 0;
						for(auto k=0; k<A_rows; k++)
						{
							v += A[ir*A_rows+k]*A[ic*A_rows+k];
						}
						M.push_back(v);
					}
				}

				// get y
				TVector y;
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

				int n = A_cols;
				int nrhs = 1;
				int lda = n;
				int ldb = n;
				int info = 0;

				vector<int> ipiv(n);

				gesv(n, nrhs, M.data(), lda, ipiv.data(), y.data(), ldb, info);
				std::copy(y.begin(), y.end(), x);
			}

			void operator()(int A_rows, int A_cols, T *A, T *b, T *x, T lambda, T *D, T *G)
			{
				using TVector = vector<T>;

				// get ATA
				TVector M;
				M.reserve(A_cols*A_cols);

				for(auto ir=0; ir<A_cols; ir++)
				{
					for(auto ic=0; ic<A_cols; ic++)
					{
						T v = 0;
						for(auto k=0; k<A_rows; k++)
						{
							v += A[ir*A_rows+k]*A[ic*A_rows+k];
						}
						T d = 0;
						if(ir == ic)
						{
							d = (v<1e-7)?lambda:lambda*v;
							D[ir] = d;
						}
						M.push_back(v+d);
					}
				}

				// get y
				TVector y;
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

				int n = A_cols;
				int nrhs = 1;
				int lda = n;
				int ldb = n;
				int info = 0;

				vector<int> ipiv(n);

				gesv(n, nrhs, M.data(), lda, ipiv.data(), y.data(), ldb, info);
				std::copy(y.begin(), y.end(), x);
			}

		private:
			void gesv(const int &n, const int &nrhs, float *a, const int &lda, 
			int *ipiv, float *b, const int &ldb, int &info)
			{
				sgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
			}

			void gesv(const int &n, const int &nrhs, double *a, const int &lda, 
			int *ipiv, double *b, const int &ldb, int &info)
			{
				dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
			}

			vector<int> ipiv;
			vector<T> M;
			vector<T> y;
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

			const mt::eDevice device;

			GESV(): device(mt::e_host) {}

			void operator()(int A_n, T *A, int b_cols, T *b, T *x)
			{
				int n = A_n;
				int nrhs = b_cols;
				int lda = n;
				int ldb = n;
				int info = 0;

				std::copy(b, b+n, x);
				ipiv.resize(n);

				gesv(n, nrhs, A, lda, ipiv.data(), x, ldb, info);
			}

		private:
			void gesv(const int &n, const int &nrhs, float *a, const int &lda, 
			int *ipiv, float *b, const int &ldb, int &info)
			{
				sgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
			}

			void gesv(const int &n, const int &nrhs, double *a, const int &lda, 
			int *ipiv, double *b, const int &ldb, int &info)
			{
				dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
			}

			vector<int> ipiv;
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

			const mt::eDevice device;

			GELS(): device(mt::e_host) {}

			void operator()(int A_rows, int A_cols, T *A, int b_cols, T *b, T *x)
			{
				const char trans = 'N';
				int m = A_rows;
				int n = A_cols;
				int nrhs = b_cols;
				int lda = A_rows;
				int ldb = max(b_cols, max(m, n));
				int min_mn = min(m, n);
				T rcond = -1;
				int info = 0;

				mt::Vector<T, mt::e_host> bv(b, b+m*b_cols);

				// query optimal size of work array
				int lwork = -1;
				T work_query= 0;
				gels(trans, m, n, nrhs, A, lda, bv.data(), ldb, &work_query, lwork, info);

				// set arrays
				lwork = max(1, int(work_query));
				work.resize(lwork);

				// perform minimum-norm solution
				gels(trans, m, n, nrhs, A, lda, bv.data(), ldb, work.data(), lwork, info);

				for(auto ix=0; ix<nrhs; ix++)
				{
					std::copy(bv.begin()+ix*m, bv.begin()+ix*m+n, x+ix*n);
				}
			}

		private:
			void gels(const char &trans, const int &m, const int &n, const int &nrhs, float *a, 
			const int &lda, float *b, const int &ldb, float *work, const int &lwork, int &info)
			{
				sgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
			}

			void gels(const char &trans, const int &m, const int &n, const int &nrhs, double *a, 
			const int &lda, double *b, const int &ldb, double *work, const int &lwork, int &info)
			{
				dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
			}

			mt::Vector<T, mt::e_host> work;
	};

	template <class T>
	struct GELSD
	{
		public:
			using value_type = T;

			const mt::eDevice device;

			GELSD(): device(mt::e_host) {}

			void operator()(int A_rows, int A_cols, T *A, int b_cols, T *b, T *x)
			{
				int m = A_rows;
				int n = A_cols;
				int nrhs = b_cols;
				int lda = A_rows;
				int ldb = max(b_cols, max(m, n));
				int min_mn = min(m, n);
				int rank = 0;
				T rcond = -1;
				int info = 0;

				mt::Vector<T, mt::e_host> bv(b, b+m*b_cols);
				S.resize(min_mn);

				// query optimal size of work array
				int lwork = -1;
				T work_query= 0;
				int iwork_query= 0;
				gelsd(m, n, nrhs, A, lda, bv.data(), ldb, S.data(), rcond, rank, &work_query, lwork, &iwork_query, info);

				// set arrays
				lwork = max(1, int(work_query));
				work.resize(lwork);

				int liwork = max(1, iwork_query);
				iwork.resize(liwork);

				// perform minimum-norm solution
				gelsd(m, n, nrhs, A, lda, bv.data(), ldb, S.data(), rcond, rank, work.data(), lwork, iwork.data(), info);

				for(auto ix= 0; ix<nrhs; ix++)
				{
					std::copy(bv.begin()+ix*m, bv.begin()+ix*m+n, x+ix*n);
				}
			}

		private:
			void gelsd(const int &m, const int &n, const int &nrhs, float *a, const int &lda, 
			float *b, const int &ldb, float *s, const float &rcond, int &rank, float *work, 
			const int &lwork, int *iwork, int &info)
			{
				sgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);
			}

			void gelsd(const int &m, const int &n, const int &nrhs, double *a, const int &lda, 
			double *b, const int &ldb, double *s, const double &rcond, int &rank, double *work, 
			const int &lwork, int *iwork, int &info)
			{
				dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);
			}

			mt::Vector<int, mt::e_host> iwork;
			mt::Vector<T, mt::e_host> work;
			mt::Vector<T, mt::e_host> S;
	};
} // namespace lapack

#endif