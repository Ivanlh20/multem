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

#pragma once

#include "const_enum.h"
#include "grid_2d.h"
#include "vctr_cpu.h"
#include "vctr_gpu.h"
#include "stream_cpu.h"
#include "stream_gpu.h"
#include "detail_cgpu.h"
#include "fcns_cpu.h"
#include "fcns_gpu.h"

namespace mt
{
	template <class T, eDev Dev>
	class Interp_rn_2d
	{
		public:
			using T_r = T;
			using TVctr = Vctr<T, Dev>;

			static const eDev device = Dev;

			Interp_rn_2d();

			Interp_rn_2d(Stream<Dev> *pstream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0);

			inline
			void set_in_data(Stream<Dev> *pstream_i, Grid_2d<T>& grid_2d_i, eFil_Sel_Typ bg_opt_i=efst_min, T bg_i=0);

			/***************************************** cpu *****************************************/
			template <eDev devn = Dev>
			enable_if_edev_cpu<devn, T>
			operator()(TVctr& mx_i, TVctr& rx_i, TVctr& ry_i, TVctr& mx_o);

			/**********************Device**********************/
		#ifdef __CUDACC__
			template <eDev devn = Dev>
			enable_if_edev_gpu<devn, T>
			operator()(TVctr& mx_i, TVctr& rx_i, TVctr& ry_i, TVctr& mx_o);
		#endif

		protected:
			Grid_2d<T> grid_2d;
			Stream<Dev> *pstream;
			eFil_Sel_Typ bg_opt;
			T bg;

			T get_bg(TVctr& M);
	};
}

#include "../src/intrpl_rn_2d.inl"