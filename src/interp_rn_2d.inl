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

#include "interp_rn_2d.h"

namespace mt
{
	template <class T, eDev Dev>
	Interp_rn_2d<T, Dev>::Interp_rn_2d(): stream(nullptr) {}

	template <class T, eDev Dev>
	Interp_rn_2d<T, Dev>::Interp_rn_2d(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, Grid_2d<T>& grid_2d_o, eFil_Sel_Typ bg_opt_i, T bg_i)
	{
		set_in_data(stream_i, grid_2d_i, grid_2d_o, bg_opt_i, bg_i);
	}

	template <class T, eDev Dev>
	inline
	void Interp_rn_2d<T, Dev>::set_in_data(Stream<Dev> *stream_i, Grid_2d<T>& grid_2d_i, Grid_2d<T>& grid_2d_o, eFil_Sel_Typ bg_opt_i, T bg_i)
	{
		stream = stream_i;
		grid_2d = grid_2d_i;
		grid_2d_mo = grid_2d_o;
		bg_opt = bg_opt_i;
		bg = bg_i;
	}

	/***************************************** cpu *****************************************/
	template <class T, eDev Dev>
	template <eDev devn = Dev>
	enable_if_edev_cpu<devn, T>
	Interp_rn_2d<T, Dev>::operator()(TVctr& mx_i, TVctr& Rx_i, TVctr& Ry_i, TVctr& mx_o)
	{
		// calculate background
		T bg = get_bg(mx_i);

		stream->set_n_stream_act(grid_2d_mo.nx);
		stream->set_grid(grid_2d_mo.nx, grid_2d_mo.ny);
		stream->exec_2d(detail_cgpu::intrpl_rg_2d<Grid_2d<T>, TVctr>, grid_2d, mx_i, Rx_i, Ry_i, grid_2d_mo, bg, mx_o);
	
		return bg;
	}

	/**********************Device**********************/
#ifdef __CUDACC__
	template <class T, eDev Dev>
	template <eDev devn = Dev>
	enable_if_edev_gpu<devn, T>
	Interp_rn_2d<T, Dev>::operator()(TVctr& mx_i, TVctr& Rx_i, TVctr& Ry_i, TVctr& mx_o)
	{
		// calculate background
		T bg = get_bg(mx_i);

		auto d_grid_blk = grid_2d_mo.d_grid_blk();
		detail_gpu::intrpl_rg_2d<Grid_2d<T>, typename TVctr::value_type><<<d_grid_blk.grid, d_grid_blk.blk>>>(grid_2d, mx_i, Rx_i, Ry_i, grid_2d_mo, bg, mx_o);
	
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