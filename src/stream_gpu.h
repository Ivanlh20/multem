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

#include <vector>

#ifdef __CUDACC__
	#include <cuda.h>
	#include <cuda_runtime.h>
#endif

#include "const_enum.h"
#include "type_traits_gen.h"
#include "kahan_sum.h"
#include "ithread_rect_1d.h"
#include "ithread_rect_2d.h"
#include "ithread_rect_3d.h"

/* forward declaration */
namespace mt
{
#ifndef STREAM_DEC
	#define STREAM_DEC
	template <eDev Dev> class Stream;
#endif
}

/* derived class */
namespace mt
{
	using Stream_gpu = Stream<edev_gpu>;
}

/* template specialization gpu stream */
namespace mt
{
#ifdef __CUDACC__
	template <>
	class Stream<edev_gpu>
	{
		public:
			static const eDev device = edev_gpu;

			Stream();

			Stream(const dt_int32& n_stream);

			Stream(const dt_int32& n_stream, const dt_shape& shape);

			~Stream();

			dt_int32 size() const;

			void resize(const dt_int32& n_stream);

			cudaStream_t& operator[](const dt_int32& iy);

			const cudaStream_t& operator[](const dt_int32& iy) const;

			void synchronize();

			void set_n_stream_act(const dt_int32& n_stream_act);

			void set_grid(const dt_shape& shape);

			void set_n_stream_act_grid(const dt_int32& n_stream_act, const dt_int32& nx);

			void set_n_stream_act_grid(const dt_int32& n_stream_act, const dt_int32& nx, const dt_int32& ny);

			void set_n_stream_act_grid(const dt_int32& n_stream_act, const dt_int32& nx, const dt_int32& ny, const dt_int32& nz);

			template <eDim Dim>
			enable_if_edim_1<Dim, iThread_Rect_xd<Dim>>
			get_ithread_rect_xd(dt_int32 istm=0);

			template <eDim Dim>
			enable_if_edim_2<Dim, iThread_Rect_xd<Dim>>
			get_ithread_rect_xd(dt_int32 istm=0);

			template <eDim Dim>
			enable_if_edim_3<Dim, iThread_Rect_xd<Dim>>
			get_ithread_rect_xd(dt_int32 istm=0);

			void cleanup();

			dt_int32 n_stream_act;

		private:
			dt_shape shape;

			dt_int32 n_stream;

			std::vector<cudaStream_t> stream;
	};
#endif
}

#include "detail/stream_gpu.inl"