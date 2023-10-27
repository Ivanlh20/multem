/*
 * This file is part of Multem.
 * Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
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

#include "intrpl_rn_2d.h"

namespace mt
{
	template <class T, eDev Dev>
	Interp_rn_2d<T, Dev>::Interp_rn_2d(): pstream(nullptr) {}

	template <class T, eDev Dev>
	Interp_rn_2d<T, Dev>::Interp_rn_2d(Stream<Dev> *pstream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i, T bg_i)
	{
		set_in_data(pstream_i, grid_2d_i, bg_opt_i, bg_i);
	}

	template <class T, eDev Dev>
	inline
	void Interp_rn_2d<T, Dev>::set_in_data(Stream<Dev> *pstream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i, T bg_i)
	{
		pstream = pstream_i;
		grid_2d = grid_2d_i;
		bg_opt = bg_opt_i;
		bg = bg_i;
	}

	/***************************************** cpu *****************************************/
	template <class T, eDev Dev>
	template <eDev devn>
	enable_if_edev_cpu<devn, T>
	Interp_rn_2d<T, Dev>::operator()(TVctr& mx_i, TVctr& rx_i, TVctr& ry_i, TVctr& mx_o)
	{
		// calculate background
		const T bg = get_bg(mx_i);

		auto igrid = mx_o.igrid_2d();

		fcn_stream_exec_xd_krn<edim_1>(pstream, mx_o.size(), detail_cgpu::fcn_intrpl_bl_rg_2d<T>, grid_2d, mx_i.m_data, rx_i.m_data, ry_i.m_data, bg, mx_o.m_data);

		return bg;
	}

	/**********************Device**********************/
#ifdef __CUDACC__
	template <class T, eDev Dev>
	template <eDev devn>
	enable_if_edev_gpu<devn, T>
	Interp_rn_2d<T, Dev>::operator()(TVctr& mx_i, TVctr& rx_i, TVctr& ry_i, TVctr& mx_o)
	{
		// calculate background
		T bg = get_bg(mx_i);

		auto igrid = mx_o.igrid_2d();

		detail_gpu::fcn_intrpl_bl_rg_2d<T><<<igrid.d_grid_size(), igrid.d_blk_size()>>>(grid_2d, mx_i.ptr_32(), rx_i.ptr_32(), ry_i.ptr_32(), bg, mx_o.ptr_32());

		return bg;
	}
#endif

	template <class T, eDev Dev>
	T Interp_rn_2d<T, Dev>::get_bg(TVctr& M)
	{
		T bg_r = 0;
		switch (bg_opt) 
		{
			case efst_min:
				bg_r = fcn_min_element(M);
			break;
			case efst_max:
				bg_r = fcn_max_element(M);
			break;
			case efst_mean:
				bg_r = fcn_mean(M);
			break;
			case efst_min_mean:
				bg_r = 0.5*(fcn_mean(M) + fcn_min_element(M));
			break;
			case efst_max_mean:
				bg_r = 0.5*(fcn_mean(M) + fcn_max_element(M));
			break;
			case efst_user_def:
				bg_r = bg;
			break;
		}

		return bg_r;
	}
}