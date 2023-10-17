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

#include <thread>
#include <vector>
#include <mutex>

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
	using Stream_cpu = Stream<edev_cpu>;
}

/* template specialization cpu stream */
namespace mt
{
	template <>
	class Stream<edev_cpu>
	{
		public:
			static const eDev device = edev_cpu;

			std::mutex stream_mutex;

			Stream();

			Stream(const dt_int32& n_stream);

			Stream(const dt_int32& n_stream, const dt_shape& shape);

			~Stream();

			dt_int32 size() const;

			dt_bool is_single_thread() const;

			void resize(const dt_int32& n_stream);

			std::thread& operator[](const dt_int32& iy);

			const std::thread& operator[](const dt_int32& iy) const;

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

			/***************************************************************************************/
			template <eDim Dim, class TFcn, class... TArgs>
			enable_if_edim_1<Dim, void>
			exec_xd_krn(const dt_int32& nx, TFcn& fcn, TArgs& ...arg);				

			template <eDim Dim, class TFcn, class... TArgs>
			enable_if_edim_2<Dim, void>
			exec_xd_krn(const dt_int32& nx, const dt_int32& ny, TFcn& fcn, TArgs& ...arg);

			template <eDim Dim, class TFcn, class... TArgs>
			enable_if_edim_3<Dim, void>
			exec_xd_krn(const dt_int32& nx, const dt_int32& ny, const dt_int32& nz, TFcn& fcn, TArgs& ...arg);

			/***************************************************************************************/
			template <eDim Dim, class TFcn, class... TArgs>
			enable_if_edim_1<Dim, void>
			exec_xd_fcn(const dt_int32& nx, TFcn& fcn, TArgs& ...arg);

			template <eDim Dim, class TFcn, class... TArgs>
			enable_if_edim_2<Dim, void>
			exec_xd_fcn(const dt_int32& nx, const dt_int32& ny, TFcn& fcn, TArgs& ...arg);

			template <eDim Dim, class TFcn, class... TArgs>
			enable_if_edim_3<Dim, void>
			exec_xd_fcn(const dt_int32& nx, const dt_int32& ny, const dt_int32& nz, TFcn& fcn, TArgs& ...arg);

 			void cleanup();

			dt_int32 n_stream_act;
		private:
			dt_shape shape;

			dt_int32 n_stream;
			std::thread *stream;
	};
}

#include "detail/stream_cpu.inl"