/*
 * This file is part of MULTEM.
 * Copyright 2016 Ivan Lobato <Ivanlh20@gmail.com>
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

namespace lapack
{
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

	/*	The program computes the solution to the system of linear
		equations with a square matrix A and multiple right-hand 
		sides B, where A is the coefficient matrix and b is the 
		right-hand side matrix:
	*/
	template<class T>
	struct GELS
	{
		public:
			using value_type = T;

			static const mt::eDevice device = mt::e_host;

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

				//query optimal size of work array
				int lwork = -1;
				T work_query= 0;
				gels(trans, m, n, nrhs, A, lda, bv.data(), ldb, &work_query, lwork, info);

				//set arrays
				lwork = max(1, int(work_query));
				work.resize(lwork);

				//perform minimum-norm solution
				gels(trans, m, n, nrhs, A, lda, bv.data(), ldb, work.data(), lwork, info);

				for(auto ix= 0; ix<nrhs; ix++)
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

	template<class T>
	struct GELSD
	{
		public:
			using value_type = T;

			static const mt::eDevice device = mt::e_host;

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

				//query optimal size of work array
				int lwork = -1;
				T work_query= 0;
				int iwork_query= 0;
				gelsd(m, n, nrhs, A, lda, bv.data(), ldb, S.data(), rcond, rank, &work_query, lwork, &iwork_query, info);

				//set arrays
				lwork = max(1, int(work_query));
				work.resize(lwork);

				int liwork = max(1, iwork_query);
				iwork.resize(liwork);

				//perform minimum-norm solution
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