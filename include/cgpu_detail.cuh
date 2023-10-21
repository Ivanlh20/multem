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

#include <type_traits>
#include <algorithm>

#include "const_enum.h"
#include "math_mt.h"
#include "type_traits_gen.h"
#include "kahan_sum.h"
#include "fcns_cgpu_gen.h"
#include "grid_1d.h"
#include "grid_2d.h"
#include "grid_3d.h"

#include <thrust/iterator/zip_iterator.h>

namespace mt
{
	/* synchronize every */
	template <eDev Dev>
	void synchronize_every(const dt_int32& i_sync, const dt_int32& n_sync)
	{

	}

	template <>
	void synchronize_every<edev_gpu>(const dt_int32& i_sync, const dt_int32& n_sync)
	{
	#ifdef __CUDACC__
		if (i_sync % n_sync == 0)
		{
			cudaDeviceSynchronize();
		}
	#endif
	}

	template <class... Args>
	auto fcn_mkzipiter_begin(Args&... args)
	{
		return thrust::make_zip_iterator(thrust::make_tuple(std::begin(args)...));
	}

	template <class... Args>
	auto fcn_mkzipiter_end(Args&... args)
	{
		return thrust::make_zip_iterator(thrust::make_tuple(std::end(args)...));
	}

	/* unrolled binary search */
	namespace cgpu_detail
	{
		template <class T> 
		CGPU_EXEC_INL 
		dt_int32 fcn_unrolled_binary_search_256(const T& x, Ctpr<T> xv)
		{
			dt_int32 i_0 = 0;
			dt_int32 i_e = 255;

			dt_int32 im = (i_0 + i_e)>>1;				// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 128
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 64
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 32
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 16
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 8
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 4
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 2
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 1
	
			return i_0;
		}

		template <class T> 
		CGPU_EXEC_INL 
		dt_int32 fcn_unrolled_binary_search_128(const T& x, Ctpr<T> xv)
		{
			dt_int32 i_0 = 0;
			dt_int32 i_e = 127;

			dt_int32 im = (i_0 + i_e)>>1;				// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 64
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 32
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 16
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 8
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 4
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 2
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 1
	
			return i_0;
		}

		template <class T> 
		CGPU_EXEC_INL 
		dt_int32 fcn_unrolled_binary_search_64(const T& x, Ctpr<T> xv)
		{
			dt_int32 i_0 = 0;
			dt_int32 i_e = 63;

			dt_int32 im = (i_0 + i_e)>>1;				// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 32
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 16
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 8
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 4
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 2
			im = (i_0 + i_e)>>1;						// divide by 2
			if (x < xv[im]) i_e = im; else i_0 = im;		// 1
	
			return i_0;
		}
	}

	/* add - assign - crop - norm_2 - fftsft */
	namespace cgpu_detail
	{
		/***************************************************************************************/
		/**************** Shift zero-frequency component to center of spectrum *****************/
		/***************************************************************************************/
		/* shift matrix respect to nx_h */
		template <class T>
		CGPU_EXEC_INL 
		void fcn_fftsft_1d(const dt_int32& ix, const iGrid_1d& igrid, T* mx_io)
		{
			const auto ix_sft = igrid.sub_2_ind(igrid.nx_h+ix);
			thrust::swap(mx_io[ix], mx_io[ix_sft]);
		}

		/* shift matrix respect to ny_h */
		template <class T>
		CGPU_EXEC_INL 
		void fcn_fftsft_bc_2d(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, T* mx_io)
		{
			const auto ixy = igrid.sub_2_ind(ix, iy);
			const auto ixy_sft = igrid.sub_2_ind(ix, igrid.ny_h+iy);
			thrust::swap(mx_io[ixy], mx_io[ixy_sft]);
		}

		/* shift matrix respect to (nx_h, ny_h) */
		template <class T>
		CGPU_EXEC_INL 
		void fcn_fftsft_2d(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, T* mx_io)
		{
			auto ixy = igrid.sub_2_ind(ix, iy);
			auto ixy_sft = igrid.sub_2_ind(igrid.nx_h+ix, igrid.ny_h+iy);
			thrust::swap(mx_io[ixy], mx_io[ixy_sft]);

			ixy = igrid.sub_2_ind(ix, igrid.ny_h+iy);
			ixy_sft = igrid.sub_2_ind(igrid.nx_h+ix, iy);
			thrust::swap(mx_io[ixy], mx_io[ixy_sft]);
		}

		/* shift matrix respect to (nx_h, ny_h) */
 		template <class T>
		CGPU_EXEC_INL 
		void fcn_fftsft_2d(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, T* mx_i, T* mx_o)
		{
			auto ixy = igrid.sub_2_ind(ix, iy);
			auto ixy_sft = igrid.sub_2_ind(igrid.nx_h+ix, igrid.ny_h+iy);
			mx_o[ixy] = mx_i[ixy_sft];
			mx_o[ixy_sft] = mx_i[ixy];

			ixy = igrid.sub_2_ind(ix, igrid.ny_h+iy);
			ixy_sft = igrid.sub_2_ind(igrid.nx_h+ix, iy);
			mx_o[ixy] = mx_i[ixy_sft];
			mx_o[ixy_sft] = mx_i[ixy];
		}

		/* shift 3d matrix respect to (nx_h, ny_h, nz_h) */
		template <class T>
		CGPU_EXEC_INL 
		void fcn_fftsft_3d(const dt_int32& ix, const dt_int32& iy, const dt_int32& iz, const iGrid_3d& igrid, T* mx_io)
		{
			auto ixy = igrid.sub_2_ind(ix, iy, iz);
			auto ixy_sft = igrid.sub_2_ind(igrid.nx_h+ix, igrid.ny_h+iy, igrid.nz_h+iz);
			thrust::swap(mx_io[ixy], mx_io[ixy_sft]);

			ixy = igrid.sub_2_ind(igrid.nx_h+ix, igrid.ny_h+iy, iz);
			ixy_sft = igrid.sub_2_ind(ix, iy, igrid.nz_h+iz);
			thrust::swap(mx_io[ixy], mx_io[ixy_sft]);

			ixy = igrid.sub_2_ind(ix, igrid.ny_h+iy, iz);
			ixy_sft = igrid.sub_2_ind(igrid.nx_h+ix, iy, igrid.nz_h+iz);
			thrust::swap(mx_io[ixy], mx_io[ixy_sft]);

			ixy = igrid.sub_2_ind(igrid.nx_h+ix, iy, iz);
			ixy_sft = igrid.sub_2_ind(ix, igrid.ny_h+iy, igrid.nz_h+iz);
			thrust::swap(mx_io[ixy], mx_io[ixy_sft]);
		}

		/* shift 3d matrix respect to (nx_h, ny_h, nz_h) */
 		template <class T>
		CGPU_EXEC_INL 
		void fcn_fftsft_3d(const dt_int32& ix, const dt_int32& iy, const dt_int32& iz, const iGrid_3d& igrid, T* mx_i, T* mx_o)
		{
			auto ixy = igrid.sub_2_ind(ix, iy, iz);
			auto ixy_sft = igrid.sub_2_ind(igrid.nx_h+ix, igrid.ny_h+iy, igrid.nz_h+iz);
			mx_o[ixy] = mx_i[ixy_sft];
			mx_o[ixy_sft] = mx_i[ixy];

			ixy = igrid.sub_2_ind(igrid.nx_h+ix, igrid.ny_h+iy, iz);
			ixy_sft = igrid.sub_2_ind(ix, iy, igrid.nz_h+iz);
			mx_o[ixy] = mx_i[ixy_sft];
			mx_o[ixy_sft] = mx_i[ixy];

			ixy = igrid.sub_2_ind(ix, igrid.ny_h+iy, iz);
			ixy_sft = igrid.sub_2_ind(igrid.nx_h+ix, iy, igrid.nz_h+iz);
			mx_o[ixy] = mx_i[ixy_sft];
			mx_o[ixy_sft] = mx_i[ixy];

			ixy = igrid.sub_2_ind(igrid.nx_h+ix, iy, iz);
			ixy_sft = igrid.sub_2_ind(ix, igrid.ny_h+iy, igrid.nz_h+iz);
			mx_o[ixy] = mx_i[ixy_sft];
			mx_o[ixy_sft] = mx_i[ixy];
		}

		/***************************************************************************************/
		/* add, scale and shift */
 		template <class T>
		CGPU_EXEC_INL 
		void fcn_add_sc_fftsft_2d(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, T* mx_i, const T& w, T* mx_o)
		{
			auto ixy = igrid.sub_2_ind(ix, iy);
			auto ixy_sft = igrid.sub_2_ind(igrid.nx_h+ix, igrid.ny_h+iy);
			mx_o[ixy] += w*mx_i[ixy_sft];
 			mx_o[ixy_sft] += w*mx_i[ixy];

			ixy = igrid.sub_2_ind(ix, igrid.ny_h+iy);
			ixy_sft = igrid.sub_2_ind(igrid.nx_h+ix, iy);
			mx_o[ixy] += w*mx_i[ixy_sft];
 			mx_o[ixy_sft] += w*mx_i[ixy];
		}

		/* add, scale, square and shift */
		template <class T, class U>
		CGPU_EXEC_INL 
		void fcn_add_sc_norm_2_fftsft_2d(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, T* mx_i, const U& w, U* mx_o)
		{
			auto ixy = igrid.sub_2_ind(ix, iy);
			auto ixy_sft = igrid.sub_2_ind(igrid.nx_h+ix, igrid.ny_h+iy);
			mx_o[ixy] += w*::norm_2(mx_i[ixy_sft]);
 			mx_o[ixy_sft] += w*::norm_2(mx_i[ixy]);

			ixy = igrid.sub_2_ind(ix, igrid.ny_h+iy);
			ixy_sft = igrid.sub_2_ind(igrid.nx_h+ix, iy);
			mx_o[ixy] += w*::norm_2(mx_i[ixy_sft]);
 			mx_o[ixy_sft] += w*::norm_2(mx_i[ixy]);
		}

		/***************************************************************************************/
 		/* Assign and crop */
 		template <class T>
		CGPU_EXEC_INL 
		void fcn_assign_crop_2d(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, T* mx_i, const iRegion_Rect_2d& iregion, T* mx_o)
		{
			if (iregion.chk_bound(ix, iy))
			{
				mx_o[iregion.sub_2_ind(ix, iy)] = mx_i[igrid.sub_2_ind(ix, iy)];
			}
		}

 		/* assign, crop and shift */
 		template <class T>
		CGPU_EXEC_INL 
		void fcn_assign_crop_fftsft_2d(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, T* mx_i, const iRegion_Rect_2d& iregion, T* mx_o)
		{
 			auto ix_i = ix;
			auto iy_i = iy;

			auto ix_s = igrid.nx_h+ix;
			auto iy_s = igrid.ny_h+iy;

			if (iregion.chk_bound(ix_i, iy_i))
			{
				mx_o[iregion.sub_2_ind(ix_i, iy_i)] = mx_i[igrid.sub_2_ind(ix_s, iy_s)];
			}

 			if (iregion.chk_bound(ix_s, iy_s))
			{
				mx_o[iregion.sub_2_ind(ix_s, iy_s)] = mx_i[igrid.sub_2_ind(ix_i, iy_i)];
			}

			/***************************************************************************************/
 			ix_i = ix;
			iy_i = igrid.ny_h+iy;

			ix_s = igrid.nx_h+ix;
			iy_s = iy;

			if (iregion.chk_bound(ix_i, iy_i))
			{
				mx_o[iregion.sub_2_ind(ix_i, iy_i)] = mx_i[igrid.sub_2_ind(ix_s, iy_s)];
			}

 			if (iregion.chk_bound(ix_s, iy_s))
			{
				mx_o[iregion.sub_2_ind(ix_s, iy_s)] = mx_i[igrid.sub_2_ind(ix_i, iy_i)];
			}
		}
			
		/* add, scale, crop and shift */
 		template <class T>
		CGPU_EXEC_INL 
		void fcn_add_sc_crop_fftsft_2d(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, T* mx_i, const iRegion_Rect_2d& iregion, const T& w, T* mx_o)
		{
 			auto ix_i = ix;
			auto iy_i = iy;

			auto ix_s = igrid.nx_h+ix;
			auto iy_s = igrid.ny_h+iy;

			if (iregion.chk_bound(ix_i, iy_i))
			{
				mx_o[iregion.sub_2_ind(ix_i, iy_i)] += w*mx_i[igrid.sub_2_ind(ix_s, iy_s)];
			}

 			if (iregion.chk_bound(ix_s, iy_s))
			{
				mx_o[iregion.sub_2_ind(ix_s, iy_s)] += w*mx_i[igrid.sub_2_ind(ix_i, iy_i)];
			}

			/***************************************************************************************/
 			ix_i = ix;
			iy_i = igrid.ny_h+iy;

			ix_s = igrid.nx_h+ix;
			iy_s = iy;

			if (iregion.chk_bound(ix_i, iy_i))
			{
				mx_o[iregion.sub_2_ind(ix_i, iy_i)] += w*mx_i[igrid.sub_2_ind(ix_s, iy_s)];
			}

 			if (iregion.chk_bound(ix_s, iy_s))
			{
				mx_o[iregion.sub_2_ind(ix_s, iy_s)] += w*mx_i[igrid.sub_2_ind(ix_i, iy_i)];
			}
		}
			
		/* add, scale, square, crop and shift */
		template <class T, class U>
		CGPU_EXEC_INL 
		void fcn_add_sc_norm_2_crop_fftsft_2d(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, T* mx_i, const iRegion_Rect_2d& iregion, const U& w, U* mx_o)
		{
 			auto ix_i = ix;
			auto iy_i = iy;

			auto ix_s = igrid.nx_h+ix;
			auto iy_s = igrid.ny_h+iy;

			if (iregion.chk_bound(ix_i, iy_i))
			{
				mx_o[iregion.sub_2_ind(ix_i, iy_i)] += w*::norm_2(mx_i[igrid.sub_2_ind(ix_s, iy_s)]);
			}

 			if (iregion.chk_bound(ix_s, iy_s))
			{
				mx_o[iregion.sub_2_ind(ix_s, iy_s)] += w*::norm_2(mx_i[igrid.sub_2_ind(ix_i, iy_i)]);
			}

			/***************************************************************************************/
 			ix_i = ix;
			iy_i = igrid.ny_h+iy;

			ix_s = igrid.nx_h+ix;
			iy_s = iy;

			if (iregion.chk_bound(ix_i, iy_i))
			{
				mx_o[iregion.sub_2_ind(ix_i, iy_i)] += w*::norm_2(mx_i[igrid.sub_2_ind(ix_s, iy_s)]);
			}

 			if (iregion.chk_bound(ix_s, iy_s))
			{
				mx_o[iregion.sub_2_ind(ix_s, iy_s)] += w*::norm_2(mx_i[igrid.sub_2_ind(ix_i, iy_i)]);
			}
		}
	}

	/* transpose - element wise matrix op vector */
	namespace cgpu_detail
	{
		/* transpose */
		template <class T>
		CGPU_EXEC_INL 
		void fcn_trs_2d(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, Ctpr<T> mx_i, Tpr<T> mx_o)
		{
			mx_o[igrid.sub_2_ind_by_d2(ix, iy)] = mx_i[igrid.sub_2_ind(ix, iy)];
		}

		/***************************************************************************************/
		/* element wise addition: matrix + vector row */
		template <class T>
		CGPU_EXEC_INL 
		void fcn_ew_add_mx_vctr_row(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& grid, Ctpr<T> vctr, Tpr<T> mx_io)
		{
			mx_io[grid.sub_2_ind(ix, iy)] += vctr[ix];
		}

		/* element wise addition: matrix + vector col */
		template <class T>
		CGPU_EXEC_INL 
		void fcn_ew_add_mx_vctr_col(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& grid, Ctpr<T> vctr, Tpr<T> mx_io)
		{
			mx_io[grid.sub_2_ind(ix, iy)] += vctr[iy];
		}

		/***************************************************************************************/
		/* element wise subtraction: matrix - vector row */
		template <class T>
		CGPU_EXEC_INL 
		void fcn_ew_sub_mx_vctr_row(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& grid, Ctpr<T> vctr, Tpr<T> mx_io)
		{
			mx_io[grid.sub_2_ind(ix, iy)] -= vctr[ix];
		}

		/* element wise subtraction: matrix - vector col */
		template <class T>
		CGPU_EXEC_INL 
		void fcn_ew_sub_mx_vctr_col(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& grid, Ctpr<T> vctr, Tpr<T> mx_io)
		{
			mx_io[grid.sub_2_ind(ix, iy)] -= vctr[iy];
		}

		/***************************************************************************************/
		/* element wise multiplication matrix X vector row */
		template <class T>
		CGPU_EXEC_INL 
		void fcn_ew_mult_mx_vctr_row(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& grid, Ctpr<T> vctr, Tpr<T> mx_io)
		{
			mx_io[grid.sub_2_ind(ix, iy)] *= vctr[ix];
		}

		/* element wise multiplication matrix X vector col */
		template <class T>
		CGPU_EXEC_INL 
		void fcn_ew_mult_mx_vctr_col(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& grid, Ctpr<T> vctr, Tpr<T> mx_io)
		{
			mx_io[grid.sub_2_ind(ix, iy)] *= vctr[iy];
		}
	}

	/* aperture functions */
	namespace cgpu_detail
	{
		template <class T, class U>
		CGPU_EXEC_INL 
		void fcn_fermi_aperture(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
			const T& g2_cut, const T& alpha, const T& w, U* mx_io)
		{
			const auto ixy = grid.sub_2_ind(ix, iy);
			const auto g2 = grid.g2_sft(ix, iy);

			mx_io[ixy] *= w*fcn_fermi_lpf(alpha, g2_cut, g2);
		}

		template <class T, class U>
		CGPU_EXEC_INL 
		void fcn_hard_aperture(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
			const T& g2_cut, const T& w, U* mx_io)
		{
			const auto ixy = grid.sub_2_ind(ix, iy);
			const auto g2 = grid.g2_sft(ix, iy);

			mx_io[ixy] = ((g2 <= g2_cut))?(w*mx_io[ixy]):T(0);
		}
	}

	/* phase shifts real space*/
	namespace cgpu_detail
	{
		// phase factor 1d
		template <class T>
		CGPU_EXEC_INL 
		void fcn_rs_exp_factor_1d(const dt_int32& ix, const Grid_1d<T>& grid, 
		complex<T>* psi_i, const T& gx, const T& w, complex<T>* psi_o)
		{
			const auto rx = grid.rx_sft(ix)-grid.rx_c();
			psi_o[ix] = w*psi_i[ix]*euler(gx*rx);
		}

		// phase factor 2d by col
		template <class T>
		CGPU_EXEC_INL 
		void fcn_rs_exp_factor_2d_bc(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
		complex<T>* psi_i, const T& alpha, Ctpr<T> gy, const T& w, complex<T>* psi_o)
		{
			const auto ixy = grid.sub_2_ind(ix, iy);
			const auto ry = grid.ry_sft(iy)-grid.ry_c();
			psi_o[ixy] = w*psi_i[ixy]*euler(alpha*gy[ix]*ry);
		}

		// phase factor 2d
		template <class T>
		CGPU_EXEC_INL 
		void fcn_rs_exp_factor_2d(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
		complex<T>* psi_i, const R_2d<T>& g, const T& w, complex<T>* psi_o)
		{
			const auto ixy = grid.sub_2_ind(ix, iy);
			const auto rv = grid.rv_sft(ix, iy)-grid.rv_c();
			psi_o[ixy] = w*psi_i[ixy]*euler(g*rv);
		}

		// phase factor 2d multipositions
		template <class T>
		CGPU_EXEC_INL 
		void fcn_rs_mul_exp_factor_2d(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
		complex<T>* psi_i, R_2d<T>* g, const dt_int32& n_g, const T& w, complex<T>* psi_o)
		{
			const auto ixy = grid.sub_2_ind(ix, iy);
			const auto rv = grid.rv_sft(ix, iy)-grid.rv_c();

			complex<T> exp_sum = 0;
			for(auto ig=0; ig<n_g; ig++)
			{
				exp_sum += euler(g[ig]*rv);
			}
			psi_o[ixy] = w*psi_i[ixy]*exp_sum;
		}
	}

	/* phase shifts fourier space*/
	namespace cgpu_detail
	{
		template <class T>
		CGPU_EXEC_INL 
		void fcn_fs_exp_factor_1d(const dt_int32& ix, const Grid_1d<T>& grid, 
		complex<T>* psi_i, const T& rx, const T& w, complex<T>* psi_o)
		{
			const auto gx = grid.gx_sft(ix);
			psi_o[ix] = w*psi_i[ix]*euler(rx*gx);
		}

		template <class T>
		CGPU_EXEC_INL 
		void fcn_fs_exp_factor_2d_bc(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
		complex<T>* psi_i, const T& alpha, T* ry, const T& w, complex<T>* psi_o)
		{
			const auto ixy = grid.sub_2_ind(ix, iy);
			const auto gy = grid.gy_sft(iy);
			psi_o[ixy] = w*psi_i[ixy]*euler(alpha*ry[ix]*gy);
		}

		template <class T>
		CGPU_EXEC_INL 
		void fcn_fs_exp_factor_2d(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
		complex<T>* psi_i, const R_2d<T>& r, const T& w, complex<T>* psi_o)
		{
			const auto ixy = grid.sub_2_ind(ix, iy);
			const auto gv = grid.gv_sft(ix, iy);
			psi_o[ixy] = w*psi_i[ixy]*euler(r*gv);
		}

		template <class T>
		CGPU_EXEC_INL 
		void fcn_fs_mul_exp_factor_2d(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, 
		complex<T>* psi_i, R_2d<T>* r, const dt_int32& n_r, const T& w, complex<T>* psi_o)
		{
			const auto ixy = grid.sub_2_ind(ix, iy);
			const auto gv = grid.gv_sft(ix, iy);

			complex<T> exp_sum = 0;
			for(auto ir=0; ir<n_r; ir++)
			{
				exp_sum += euler(r[ir]*gv);
			}
			psi_o[ixy] = w*psi_i[ixy]*exp_sum;
		}
	}

	/* gradient */
	namespace cgpu_detail
	{
		template <class T>
		CGPU_EXEC_INL 
		T fcn_grad_x(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, Ctpr<T> mx_i)
		{
			const auto ixy = igrid.sub_2_ind(ix, iy);

			if (ix == 0)
			{
				const auto ixy_f = igrid.sub_2_ind(ix+1, iy);

				return mx_i[ixy_f]-mx_i[ixy];
			}
			else if (ix == igrid.nx-1)
			{
				const auto ixy_b = igrid.sub_2_ind(ix-1, iy);

				return mx_i[ixy]-mx_i[ixy_b];
			}
			else
			{
				const auto ixy_b = igrid.sub_2_ind(ix-1, iy);
				const auto ixy_f = igrid.sub_2_ind(ix+1, iy);

				return (mx_i[ixy_f]-mx_i[ixy_b])/T(2);
			}
		}

		template <class T>
		CGPU_EXEC_INL 
		void fcn_grad_x(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, Ctpr<T> mx_i, Tpr<T> mx_o)
		{
			mx_o[igrid.sub_2_ind(ix, iy)] = fcn_grad_x(ix, iy, igrid, mx_i);
		}

		/***************************************************************************************/
		template <class T>
		CGPU_EXEC_INL 
		T fcn_grad_y(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, Ctpr<T> mx_i)
		{
			const auto ixy = igrid.sub_2_ind(ix, iy);

			if (iy == 0)
			{
				const auto ixy_f = igrid.sub_2_ind(ix, iy+1);

				return mx_i[ixy_f]-mx_i[ixy];
			}
			else if (iy == igrid.ny-1)
			{
				const auto ixy_b = igrid.sub_2_ind(ix, iy-1);

				return mx_i[ixy]-mx_i[ixy_b];
			}
			else
			{
				const auto ixy_b = igrid.sub_2_ind(ix, iy-1);
				const auto ixy_f = igrid.sub_2_ind(ix, iy+1);

				return (mx_i[ixy_f]-mx_i[ixy_b])/T(2);
			}
		}

		template <class T>
		CGPU_EXEC_INL 
		void fcn_grad_y(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, Ctpr<T> mx_i, Tpr<T> mx_o)
		{
			mx_o[igrid.sub_2_ind(ix, iy)] = fcn_grad_y(ix, iy, igrid, mx_i);
		}

		/***************************************************************************************/
		//template <class T>
		//CGPU_EXEC_INL 
		//R_2d<T> fcn_grad(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, Ctpr<T> mx_i)
		//{
		//	return {fcn_grad_x(ix, iy, igrid, mx_i), fcn_grad_y(ix, iy, igrid, mx_i)};
		//}

		template <class T>
		CGPU_EXEC_INL 
		void fcn_grad(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, Ctpr<T> mx_i, Tpr<T> dm_x, Tpr<T> dm_y)
		{
			fcn_grad_x(ix, iy, igrid, mx_i, dm_x);
			fcn_grad_y(ix, iy, igrid, mx_i, dm_y);
		}
	}

	/* function multiplication fourier space */
	namespace cgpu_detail
	{
		#define FCN_MULT_FS_FCN_CGPU(POW, DIM)														\
		template <class T, class U, class TFcn>														\
		CGPU_EXEC_INL 																				\
		void fcn_mult_fs_fcn_g##POW##_##DIM##d(IDX_CDEF_ND(DIM), const Grid_##DIM##d<T>& grid, 		\
		const TFcn& fcn, const T& w, Tpr<U> mx_io)													\
		{																							\
			const auto fg = fcn(grid.g##POW##_sft(IDX_ND(DIM)));									\
			mx_io[grid.sub_2_ind(IDX_ND(DIM))] *= w*fg;												\
		}

		FCN_MULT_FS_FCN_CGPU(, 1);		// fcn_mult_fs_fcn_g_1d
		FCN_MULT_FS_FCN_CGPU(, 2);		// fcn_mult_fs_fcn_g_2d
		FCN_MULT_FS_FCN_CGPU(, 3);		// fcn_mult_fs_fcn_g_3d

		FCN_MULT_FS_FCN_CGPU(2, 1);		// fcn_mult_fs_fcn_g2_1d
		FCN_MULT_FS_FCN_CGPU(2, 2);		// fcn_mult_fs_fcn_g2_2d
		FCN_MULT_FS_FCN_CGPU(2, 3);		// fcn_mult_fs_fcn_g2_3d
	}

	/* deconvolution */
	namespace cgpu_detail
	{
		#define FCN_DCV_FS_FCN_CGPU(POW, DIM)													\
		template <class T, class U, class TFcn>													\
		CGPU_EXEC_INL 																			\
		void fcn_dcv_fs_fcn_g##POW##_##DIM##d(IDX_CDEF_ND(DIM), const Grid_##DIM##d<T>& grid, 	\
		const TFcn& fcn, const T& psnr, const T& w, Tpr<U> mx_io)								\
		{																						\
			const auto fg = fcn(grid.g##POW##_sft(IDX_ND(DIM)));								\
			mx_io[grid.sub_2_ind(IDX_ND(DIM))] *= w*fg/(fg*fg+psnr);							\
		}

		FCN_DCV_FS_FCN_CGPU(, 1);		// fcn_dcv_fs_fcn_g_1d
		FCN_DCV_FS_FCN_CGPU(, 2);		// fcn_dcv_fs_fcn_g_2d
		FCN_DCV_FS_FCN_CGPU(, 3);		// fcn_dcv_fs_fcn_g_3d

		FCN_DCV_FS_FCN_CGPU(2, 1);		// fcn_dcv_fs_fcn_g2_1d
		FCN_DCV_FS_FCN_CGPU(2, 2);		// fcn_dcv_fs_fcn_g2_2d
		FCN_DCV_FS_FCN_CGPU(2, 3);		// fcn_dcv_fs_fcn_g2_3d
	}

	/* window functions */
	namespace cgpu_detail
	{
		//template <class T, class U, eFcn_typ Fcn_typ>
		//CGPU_EXEC_INL
		//void fcn_wd_fcn_r2_xd(const dt_int32& ix, const Grid_xd<T, edim_1>& grid,
		//const Wd_fcn_xd<T, edim_1, Fcn_typ>& fcn, Tpr<U> mx_o)
		//{
		//	mx_o[grid.sub_2_ind(ix)] = fcn.eval_r2(grid.r2(ix, fcn.r));
		//}			
		//
		//template <class T, class U, eFcn_typ Fcn_typ>
		//CGPU_EXEC_INL
		//void fcn_wd_fcn_r2_xd(const dt_int32& ix, const dt_int32& iy, const Grid_xd<T, edim_2>& grid,
		//const Wd_fcn_xd<T, edim_2, Fcn_typ>& fcn, Tpr<U> mx_o)
		//{
		//	mx_o[grid.sub_2_ind(ix, iy)] = fcn.eval_r2(grid.r2(ix, iy, fcn.r));
		//}

		//template <class T, class U, eFcn_typ Fcn_typ>
		//CGPU_EXEC_INL
		//void fcn_wd_fcn_r2_xd(const dt_int32& ix, const dt_int32& iy, const dt_int32& iz, const Grid_xd<T, edim_3>& grid,
		//const Wd_fcn_xd<T, edim_3, Fcn_typ>& fcn, Tpr<U> mx_o)
		//{
		//	mx_o[grid.sub_2_ind(ix, iy, iz)] = fcn.eval_r2(grid.r2(ix, iy, iz, fcn.r));
		//}

		///***************************************************************************************/
		//template <class T, class U, eFcn_typ Fcn_typ>
		//CGPU_EXEC_INL
		//void fcn_wd_fcn_r2_sft_xd(const dt_int32& ix, const Grid_xd<T, edim_1>& grid,
		//const Wd_fcn_xd<T, edim_1, Fcn_typ>& fcn, Tpr<U> mx_o)
		//{
		//	mx_o[grid.sub_2_ind(ix)] = fcn.eval_r2(grid.r2_sft(ix, fcn.r));
		//}			
		//
		//template <class T, class U, eFcn_typ Fcn_typ>
		//CGPU_EXEC_INL
		//void fcn_wd_fcn_r2_sft_xd(const dt_int32& ix, const dt_int32& iy, const Grid_xd<T, edim_2>& grid,
		//const Wd_fcn_xd<T, edim_2, Fcn_typ>& fcn, Tpr<U> mx_o)
		//{
		//	mx_o[grid.sub_2_ind(ix, iy)] = fcn.eval_r2(grid.r2_sft(ix, iy, fcn.r));
		//}

		//template <class T, class U, eFcn_typ Fcn_typ>
		//CGPU_EXEC_INL
		//void fcn_wd_fcn_r2_sft_xd(const dt_int32& ix, const dt_int32& iy, const dt_int32& iz, const Grid_xd<T, edim_3>& grid,
		//const Wd_fcn_xd<T, edim_3, Fcn_typ>& fcn, Tpr<U> mx_o)
		//{
		//	mx_o[grid.sub_2_ind(ix, iy, iz)] = fcn.eval_r2(grid.r2_sft(ix, iy, iz, fcn.r));
		//}

		#define FCN_WD_FCN_CGPU(DIM)															\
		template <class T, class U, class TFcn>													\
		CGPU_EXEC_INL 																			\
		void fcn_wd_fcn_r2_##DIM##d(IDX_CDEF_ND(DIM), const Grid_##DIM##d<T>& grid, 			\
		const TFcn& fcn, Tpr<U> mx_o)															\
		{																						\
			mx_o[grid.sub_2_ind(IDX_ND(DIM))] = fcn.eval_r2(grid.r2(IDX_ND(DIM), fcn.r));		\
		}

		FCN_WD_FCN_CGPU(1);			// fcn_wd_fcn_r2_1d
		FCN_WD_FCN_CGPU(2);			// fcn_wd_fcn_r2_2d
		FCN_WD_FCN_CGPU(3);			// fcn_wd_fcn_r2_3d

		#define FCN_WD_FCN_SFT_CGPU(DIM)														\
		template <class T, class U, class TFcn>													\
		CGPU_EXEC_INL 																			\
		void fcn_wd_fcn_r2_sft_##DIM##d(IDX_CDEF_ND(DIM), const Grid_##DIM##d<T>& grid, 		\
		const TFcn& fcn, Tpr<U> mx_o)															\
		{																						\
			mx_o[grid.sub_2_ind(IDX_ND(DIM))] = fcn.eval_r2(grid.r2_sft(IDX_ND(DIM), fcn.r));	\
		}

		FCN_WD_FCN_SFT_CGPU(1);		// fcn_wd_fcn_r2_sft_1d
		FCN_WD_FCN_SFT_CGPU(2);		// fcn_wd_fcn_r2_sft_2d
		FCN_WD_FCN_SFT_CGPU(3);		// fcn_wd_fcn_r2_sft_3d
	}

	/* phase correlation */
	namespace cgpu_detail
	{
		/****************** pcf data processing real space *******************/
		#define FCN_RS_PCF_XD_DP(DIM)															\
		template <class T, class U, class TFcn>													\
		CGPU_EXEC_INL 																			\
		void fcn_rs_pcf_##DIM##d_dp(IDX_CDEF_ND(DIM), const Grid_##DIM##d<T>& grid, 			\
		Ctpr<T> mx_i, const TFcn& wd, const T& w, Tpr<U> mx_o)									\
		{																						\
			const auto ind_0 = grid.sub_2_ind(IDX_ND(DIM));										\
			const auto ind_n = grid.sub_2_ind_bd(IDX_ND(DIM)+1);								\
			mx_o[ind_0] = w*(mx_i[ind_n]-mx_i[ind_0])*wd(grid.r2(IDX_ND(DIM), wd.r_c));			\
		}

		FCN_RS_PCF_XD_DP(1);		// fcn_rs_pcf_1d_dp
		FCN_RS_PCF_XD_DP(2);		// fcn_rs_pcf_2d_dp
		FCN_RS_PCF_XD_DP(3);		// fcn_rs_pcf_3d_dp

		/**************** mx_o data processing fourier space ******************/
		#define FCN_FS_PCF_XD_DP(DIM)															\
		template <class T, class U, class TFcn>													\
		CGPU_EXEC_INL 																			\
		void fcn_fs_pcf_##DIM##d_dp(IDX_CDEF_ND(DIM), const Grid_##DIM##d<T>& grid, 			\
		Ctpr<T> mx_1i, Ctpr<T> mx_2i, const TFcn& wd, const T& w, Tpr<U> mx_o)					\
		{																						\
			const auto ind_0 = grid.sub_2_ind(IDX_ND(DIM));										\
			auto z = conj(mx_1i[ind_0])*mx_2i[ind_0];											\
			mx_o[ind_0] = polar(wd(grid.g2_sft(IDX_ND(DIM))), arg(z));							\
		}

		FCN_FS_PCF_XD_DP(1);		// fcn_fs_pcf_1d_dp
		FCN_FS_PCF_XD_DP(2);		// fcn_fs_pcf_2d_dp
		FCN_FS_PCF_XD_DP(3);		// fcn_fs_pcf_3d_dp
	}

	/* optical flow */
	namespace cgpu_detail
	{
		// https:// en.wikipedia.org/wiki/Optical_flow
		template <class T>
		CGPU_EXEC_INL 
		void fcn_opt_flow(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& igrid, Ctpr<T> mx_s, Ctpr<T> mx_m, 
		const T& alpha, Tpr<T> dphi_x, Tpr<T> dphi_y)
		{
			const auto ixy = igrid.sub_2_ind(ix, iy);

			const auto sm = mx_m[ixy]-mx_s[ixy];
			const auto alpha_sm_2 = ::square(alpha*sm);

			auto dmx_s = fcn_grad(ix, iy, igrid, mx_s);
			auto dmx_m = fcn_grad(ix, iy, igrid, mx_m);

			const auto dphi = -sm*(dmx_s/(dmx_s.norm_2() + alpha_sm_2) + dmx_s/(dmx_m.norm_2() + alpha_sm_2));

			if (isfinite(dphi.x) && isfinite(dphi.y))
			{
				dphi_y[ixy] = dphi.x;
				dphi_x[ixy] = dphi.y;
			}
			else
			{
				dphi_y[ixy] = T(0);
				dphi_x[ixy] = T(0);
			}
		}
	}
	
	/* interpolation poly3 */
	namespace cgpu_detail
	{
		// calculate interpolation coefficients from its value and derivative
		template <class T>
		CGPU_EXEC_INL 
		void fcn_vd_2_coef_poly3(const T& dx_ij, const T& y_i, const T& y_j, const T& t_i, const T& t_j, T& c2, T& c3)
		{
			const auto m_i = (y_j-y_i)/dx_ij;
			const auto n_i = t_i + t_j - T(2)*m_i;

			c2 = (m_i-t_i-n_i)/dx_ij;
			c3 = n_i/(dx_ij*dx_ij);
		}

		// eval polynomial of order 3
		template <class T>
		CGPU_EXEC_INL 
		T fcn_eval_poly3(const T& x, Ctpr<T> xv, Ctpr<T> c0, Ctpr<T> c1, Ctpr<T> c2, Ctpr<T> c3)
		{
			const auto ix = fcn_unrolled_binary_search_128<T>(x, xv);
			const auto dx = x - xv[ix];

			return c0[ix] + (c1[ix] +(c2[ix] + c3[ix]*dx)*dx)*dx;
		}
	}

	/* bilinear interpolation */
	namespace cgpu_detail
	{
		/* regular grid bilinear interpolation */
		// https:// en.wikipedia.org/wiki/Bilinear_interpolation
 		template <class T>
		CGPU_EXEC_INL 
		T fcn_intrpl_bl_rg_2d(const R_2d<T>& p, const Grid_2d<T>& grid_i, Tpr<T> mx_i)
		{
			const auto ix = fcn_set_bound(grid_i.rx_2_irx_fds(p.x), 0, grid_i.nx-2);
			const auto iy = fcn_set_bound(grid_i.ry_2_iry_fds(p.y), 0, grid_i.ny-2);

			const auto f11 = mx_i[grid_i.sub_2_ind(ix, iy)];
			const auto f12 = mx_i[grid_i.sub_2_ind(ix, iy+1)];
			const auto f21 = mx_i[grid_i.sub_2_ind(ix+1, iy)];
			const auto f22 = mx_i[grid_i.sub_2_ind(ix+1, iy+1)];

			const auto x1 = grid_i.rx(ix);
			const auto x2 = grid_i.rx(ix+1);
			const auto y1 = grid_i.ry(iy);
			const auto y2 = grid_i.ry(iy+1);

			const auto dx1 = (p.x-x1)/(x2-x1);
			const auto dx2 = (x2-p.x)/(x2-x1);
			const auto dy1 = (p.y-y1)/(y2-y1);
			const auto dy2 = (y2-p.y)/(y2-y1);

			return dx2*(f11*dy2 + f12*dy1)+dx1*(f21*dy2 + f22*dy1);
		};

		template <class T>
		CGPU_EXEC_INL 
		void fcn_intrpl_bl_rg_2d(const dt_int32& ixy, const Grid_2d<T>& grid_i, Ctpr<T> mx_i, 
		Ctpr<T> vrx, Ctpr<T> vry, const T& bg, Tpr<T> mx_o)
		{
			R_2d<T> p(vrx[ixy], vry[ixy]);
			mx_o[ixy] = (grid_i.chk_bound_eps(p))?fcn_intrpl_bl_rg_2d(p, grid_i, mx_i):bg;
		};

		/***************************************************************************************/
		/* non regular grid bilinear interpolation */
		// https:// www.codesd.com/item/bilinear-interpolation-with-non-aligned-entry-points.html
 		template <class T>
		CGPU_EXEC_INL 
		T fcn_intrpl_bl_nrg_2d(const R_2d<T>& p, const R_3d<T>& p1, 
		const R_3d<T>& p2, const R_3d<T>& p3, const R_3d<T>& p4)
		{
			const auto a = -p1.x + p3.x;
			const auto b = -p1.x + p2.x;
			const auto c = p1.x - p2.x - p3.x + p4.x;
			const auto d = p.x - p1.x;
			const auto e = -p1.y + p3.y;
			const auto f = -p1.y + p2.y;
			const auto g = p1.y - p2.y - p3.y + p4.y;
			const auto h = p.y - p1.y;
	
			const auto c_x = ::sqrt(T(-4)*(c*e - a*g)*(d*f - b*h) + ::square(b*e - a*f + d*g - c*h));
			const auto c_y = T(2)*(c*e - a*g);
			const auto c_z = T(2)*(c*f - b*g);

			T alpha = -(b*e - a*f + d*g - c*h + c_x)/c_y;
			T beta = (b*e - a*f - d*g + c*h + c_x)/c_z;

			if ((alpha < T(0)) || (alpha > T(1)) || (beta < T(0)) || (beta > T(1)))
			{
				alpha = (-b*e + a*f - d*g + c*h + c_x)/c_y;
				beta = -(-b*e + a*f + d*g - c*h + c_x)/c_z;
			}

			return (T(1) - alpha)*((T(1) - beta)*p1.z + beta*p2.z) + alpha*((T(1) - beta)*p3.z + beta*p4.z);
		};

		template <class T>
		CGPU_EXEC_INL 
		void fcn_intrpl_bl_nrg_2d(const dt_int32& ixy, Ctpr<T> vrx_i, Ctpr<T> vry_i, Ctpr<T>mx_i, 
		Ctpr<T> vrx, Ctpr<T> vry, const T& bg, Tpr<T> mx_o)
		{
			R_2d<T> p(vrx[ixy], vry[ixy]);
			// neighboring search algorithm: check out tessellation or MD find neighbors
			R_3d<T> p1;
			R_3d<T> p2;
			R_3d<T> p3;
			R_3d<T> p4;

			// how to check out bound for non regular grid?
			// mx_o[ixy] = (grid_2d_i.chk_bound_eps(p))?fcn_intrpl_bl_rg_2d(p, p1, p2, p3, p4, mx_i):bg;
		};
	}

	namespace cgpu_detail
	{
		/* distort regular grid */
		template <class T>
		CGPU_EXEC_INL 
		void fcn_distort_mx_2d(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, Ctpr<T> mx, 
		Ctpr<T> dvx, Ctpr<T> dvy, const T& bg, Tpr<T> mx_o)
		{
			const R_2d<T> p(grid.rx(ix)+dvx[iy], grid.ry(iy)+dvy[iy]);
			mx_o[grid.sub_2_ind(ix, iy)] = (grid.chk_bound_eps(p))?fcn_intrpl_bl_rg_2d(p, grid, mx):bg;
		};

		/* scan distortion */
		template <class T, class TPar>
		CGPU_EXEC_INL 
		void fcn_sd_nr_2d(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid, Ctpr<T> mx, 
		TPar &parm, Tpr<T> mx_o)
		{
			const auto ixy = grid.sub_2nd(ix, iy);
			const auto ny = grid.ny-2;
			const auto nx = grid.nx-2;

			const R_2d<T> p(grid.rx(ix), grid.ry(iy));

			const auto iy_s1 = fcn_r_2r_b_by_vctr(parm.ry_s, p.y, 0, ny);
			const auto iy_s2 = iy_s1 + 1;

			if ((iy_s1 < 0)||(iy_s2 > ny))
			{
				mx_o[ixy] = (parm.bg_opt == efst_same_in)?mx[ixy]:parm.bg;
				return;
			}

			const auto iy_k1 = parm.iy[iy_s1];
			const auto ix_k1 = grid.rx_2rx_fds(p.x - parm.dx[iy_k1]);
			const auto iy_k2 = parm.iy[iy_s2];
			const auto ix_k2 = grid.rx_2rx_fds(p.x - parm.dx[iy_k2]);

			if ((ix_k1 < 0)||(ix_k1 > nx)||(ix_k2 < 0)||(ix_k2 > nx))
			{
				mx_o[ixy] = (parm.bg_opt == efst_same_in)?mx[ixy]:parm.bg;
				return;
			}

			const R_3d<T> p1(grid.rx(ix_k1) + parm.dx[iy_k1], parm.ry_s[iy_s1], mx[grid.sub_2nd(ix_k1, iy_k1)]);
			const R_3d<T> p4(grid.rx(ix_k1+1) + parm.dx[iy_k1], parm.ry_s[iy_s1], mx[grid.sub_2nd(ix_k1+1, iy_k1)]);
			const R_3d<T> p2(grid.rx(ix_k2) + parm.dx[iy_k2], parm.ry_s[iy_s2], mx[grid.sub_2nd(ix_k2, iy_k2)]);
			const R_3d<T> p3(grid.rx(ix_k2+1) + parm.dx[iy_k2], parm.ry_s[iy_s2], mx[grid.sub_2nd(ix_k2+1, iy_k2)]);

			/* this must work */
			// mx_o[ixy] = fcn_intrpl_bl_nrg_2d(p, p1, p2, p3, p4, grid);

			auto u = p2-p1;
			auto v = p4-p1;
			p = p - p1;
			T m = v.x*u.y-v.y*u.x;
			T alpha = (u.y*p.x-u.x*p.y)/m;
			T beta = (-v.y*p.x+v.x*p.y)/m;

			mx_o[ixy] = beta*(alpha*p3.z + (1-alpha)*p2.z) + (1-beta)*(alpha*p4.z + (1-alpha)*p1.z);
		};

		// /* find peak maximum 1d */
		// template <class T>
		// CGPU_EXEC_INL 
		// T fd_max_peak_pos(const Grid_1d<T>& grid, Ctpr<T> M)
		// {
		// 	const dt_int32 ix_max = thrust::distance(M.begin(), thrust::max_element(M.begin(), M.end()));

		// 	return grid.rx(ix_max);
		// };

		// /* find peak maximum 2d */
		// template <class T>
		// CGPU_EXEC FORCE_INLINE 
		// R_2d<T> fd_max_peak_pos(const Grid_2d<T>& grid, Ctpr<T>M)
		// {
		// 	const dt_int32 ixy_max = thrust::distance(M.begin(), thrust::max_element(M.begin(), M.end()));
		// 	const dt_int32 ix_max, iy_max;
		// 	grid.ind_2_sub(ixy_max, ix_max, iy_max);

		// 	return {grid.rx(ix_max), grid.ry(iy_max)};
		// };

		// template <class T>
		// CGPU_EXEC_INL 
		// R_2d<T> af_iscale(const R_2d<T>& p, const T& fxy)
		// {
		// 	return p/fxy;
		// }

		// template <class T>
		// CGPU_EXEC_INL 
		// void sc_2d(const dt_int32& ix, const dt_int32& iy, const Grid_2d<T>& grid_2d_i, TVctr& mx_i, 
		// const T& fxy, const Grid_2d<T>& grid_2d_o, TVctr& mx_o)
		// {
		// 	using T = T;

		// 	R_2d<T> p(grid_2d_o.rx(ix), grid_2d_o.ry(iy));
		// 	p = af_iscale(p, fxy);

		// 	mx_o[grid_2d_o.sub_2_ind(ix, iy)] = fcn_intrpl_bl_rg_2d(p, grid_2d_i, mx_i);
		// }

		// template <class T>
		// CGPU_EXEC_INL 
		// R_2d<T> af_irot_sca_sft(const T& theta, const R_2d<T>& p0, const T& fx, const T& fy, const R_2d<T>& ps, R_2d<T> p)
		// {
		// 	T sin_t, cos_t;
		// 	sincos(theta, &sin_t, &cos_t);
		// 	p.x = (p.x-ps.x)/fx;
		// 	p.y = (p.y-ps.y)/fy;
		// 	p -= p0;
		// 	p = R_2d<T>(cos_t*p.x+sin_t*p.y, -sin_t*p.x+cos_t*p.y);
		// 	p += p0;
		// 	return p;
		// }

		// template <class T>
		// CGPU_EXEC_INL 
		// void rot_sca_sft_2d(const dt_int32& ix, const dt_int32& iy, const TGrid& grid_2d_i, TVctr& mx_i, 
		// const T& theta, const R_2d<T>& p0, const T& fx, const T& fy, 
		// const R_2d<T>& ps, const T& bg, const TGrid& grid_2d_o, TVctr& mx_o)
		// {
		// 	using T = T;

		// 	R_2d<T> p(grid_2d_o.rx(ix), grid_2d_o.ry(iy));
		// 	p = af_irot_sca_sft(theta, p0, fx, fy, ps, p);

		// 	mx_o[grid_2d_o.sub_2_ind(ix, iy)] = (grid_2d_i.chk_bound_eps(p))?fcn_intrpl_bl_rg_2d(p, grid_2d_i, mx_i):bg;
		// };

		// template <class T>
		// CGPU_EXEC_INL 
		// void at_2d(const dt_int32& ix, const dt_int32& iy, const iGrid_2d& grid, TVctr& mx_i, 
		// const Mx_2x2<T>& A, const R_2d<T>& txy, const T& bg, 
		// TVctr& mx_o)
		// {
		// 	using T = T;

		// 	R_2d<T> p(grid.rx(ix), grid.ry(iy));
		// 	p = A*p+txy;

		// 	mx_o[grid.sub_2_ind(ix, iy)] = (grid.chk_bound_eps(p))?fcn_intrpl_bl_rg_2d(p, grid, mx_i):bg;
		// };

 	// 	template <class T>
		// CGPU_EXEC_INL 
		// R_2d<T> af_shx_scy(const R_2d<T>& p, const T& a, const T& b)
		// {
		// 	return R_2d<T>(p.x+a*p.y, b*p.y);
		// }

		// template <class T>
		// CGPU_EXEC_INL 
		// void shx_scy(const dt_int32& ix, const dt_int32& iy, const TGrid& grid_2d_i, TVctr& mx_i, 
		// const T& fx, const T& fy, const T& bg, 
		// const TGrid& grid_2d_o, TVctr& mx_o)
		// {
		// 	using T = T;

		// 	R_2d<T> p(grid_2d_o.rx(ix), grid_2d_o.ry(iy));
		// 	p = af_shx_scy(p, fx, fy);

		// 	mx_o[grid_2d_o.sub_2_ind(ix, iy)] = (grid_2d_i.chk_bound_eps(p))?fcn_intrpl_bl_rg_2d(p, grid_2d_i, mx_i):bg;
		// };


		// template <class T>
		// CGPU_EXEC_INL 
		// void dilate(const dt_int32& ix_i, const dt_int32& iy_i, TVctr& Im_i, TVctr& Im_o)
		// {
		// 	dt_int32 ix_0 = max(ix_i+nk0, 0);
		// 	dt_int32 ix_e = min(ix_i+nke, nx_i);

		// 	dt_int32 iy_0 = max(iy_i+nk0, 0);
		// 	dt_int32 iy_e = min(iy_i+nke, ny_i);

		// 	for(auto ix = ix_0; ix < ix_e; ix++)
		// 	{
		// 		for(auto iy = iy_0; iy < iy_e; iy++)
		// 		{
		// 			if (Im_i[ix*ny_i+iy]>0.5)
		// 			{	
		// 				Im_o[ix_i*ny_i+iy_i] = 1;
		// 				return;
		// 			}
		// 		}
		// 	}
		// 	Im_o[ix_i*ny_i+iy_i] = 0;
		// }
	} // cgpu_detail
}
