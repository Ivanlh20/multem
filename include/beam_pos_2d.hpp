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

#ifndef BEAM_POS_2D_H
	#define BEAM_POS_2D_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include <type_traits>
	#include <vector>
	#include <algorithm>
	#include <thread>
	#include <mutex>
	#include <utility>
	#include <functional> 

	#include "const_enum_mt.cuh"
	#include "math_mt.h"
	#include "type_traits_gen.h"
	#include "fcns_cgpu_gen.h"
	#include "r_2d.h"
	#include "r_3d.h"
	#include "particles.cuh"
	#include "types.cuh"

	namespace mt
	{
		template <class T>
		struct Beam_Pos_2d
		{
			Vctr<dt_int32, edev_cpu> idx;		// index
			Vctr<R_2d<T>, edev_cpu> p;			// xy-position 

			Beam_Pos_2d() {}

			Beam_Pos_2d(dt_int32 new_size) { resize(new_size); }

			template <class TVctr>
			void set_in_data(const TVctr& x) 
			{ 
				dt_int32 n = x.cols;
				resize(n);
				for(auto ik = 0; ik<n; ik++)
				{
					idx[ik] = ik;
					p[ik] = R_2d<T>(x(0, ik), x(1, ik));
				}
			}

			template <class TVctr>
			void set_in_data(const TVctr& x, const TVctr& y) 
			{
				dt_int32 n = min(x.size(), y.size());
				resize(n);
				for(auto ik = 0; ik<n; ik++)
				{
					idx[ik] = ik;
					p[ik] = R_2d<T>(x[ik], y[ik]);
				}
			}

			typename Vctr<R_2d<T>, edev_cpu>::iterator begin() const 
			{ 
				return p.begin();
			}

			typename Vctr<R_2d<T>, edev_cpu>::iterator end() const 
			{ 
				return p.end();
			}

			dt_int32 size() const { return idx.size(); }

			void resize(dt_int32 new_size)
			{
				idx.resize(new_size);
				p.resize(new_size);
			}

			void init()
			{
				thrust::fill(idx.begin(), idx.end(), 0);
				thrust::fill(p.begin(), p.end(), R_2d<T>());
			}

			void clear()
			{
				idx.clear();
				p.clear();
			}

			void shrink_to_fit()
			{
				idx.shrink_to_fit();
				p.shrink_to_fit();
			}

			void clear_shrink_to_fit()
			{
				clear();
				shrink_to_fit();
			}

			Vctr<T, edev_cpu> x()
			{
				Vctr<T, edev_cpu> x;
				x.reserve(p.size());
				for(auto ik = 0; ik<p.size(); ik++)
				{
					x.push_back(p[ik].x);
				}

				return x;
			}

			Vctr<T, edev_cpu> y()
			{
				Vctr<T, edev_cpu> y;
				y.reserve(p.size());
				for(auto ik = 0; ik<p.size(); ik++)
				{
					y.push_back(p[ik].y);
				}

				return y;
			}

			template <class TBeam_Pos_2d>
			void assign(TBeam_Pos_2d &beam_pos)
			{
				resize(beam_pos.size());
				thrust::copy(beam_pos.idx.begin(), beam_pos.idx.end(), idx.begin());
				thrust::copy(beam_pos.p.begin(), beam_pos.p.end(), p.begin());
			}

			template <class TBeam_Pos_2d>
			Beam_Pos_2d<T>& operator=(const TBeam_Pos_2d &beam_pos)
			{
			 assign(beam_pos);
			 return *this;
			}
		};
	}

#endif