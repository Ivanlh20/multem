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

#ifndef CPU_DETAIL_H
	#define CPU_DETAIL_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include <thread>
	#include <type_traits>
	#include <algorithm>

	#include "macros.h"
	#include "const_enum.h"
	#include "math_mt.h"
	#include "types.cuh"
	#include "type_traits_gen.h"

	namespace mt
	{
		/* for loops */
		namespace detail_cpu
		{
			template <class TFcn, class... TArgs>
			void for_loop(const iThread_Rect_1d& range, TFcn& fcn, TArgs& ...arg)
			{
				for(auto ixy = range.ind_0; ixy < range.ind_e; ixy++)
				{
					fcn(ixy, arg...);
				}
			}			
			
			template <class TFcn, class... TArgs>
			void for_loop(const iThread_Rect_2d& range, TFcn& fcn, TArgs& ...arg)
			{
				for(auto ix = range.ix_0; ix < range.ix_e; ix++)
				{
					for(auto iy = range.iy_0; iy < range.iy_e; iy++)
					{
						fcn(ix, iy, arg...);
					}
				}
			}

			template <class TFcn, class... TArgs>
			void for_loop(const iThread_Rect_3d& range, TFcn& fcn, TArgs& ...arg)
			{
				for(auto iz = range.iz_0; iz < range.iz_e; iz++)
				{
					for(auto ix = range.ix_0; ix < range.ix_e; ix++)
					{
						for(auto iy = range.iy_0; iy < range.iy_e; iy++)
						{
							fcn(ix, iy, iz, arg...);
						}
					}
				}
			}
		}

	}

#endif
