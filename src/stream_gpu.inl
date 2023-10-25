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

#include "stream_gpu.h"

/* template specialization gpu stream */
namespace mt
{
#ifdef __CUDACC__
	Stream<edev_gpu>::Stream(): shape(), n_stream(0), n_stream_act(0) {}

	Stream<edev_gpu>::Stream(const dt_int32& n_stream): Stream()
	{
		resize(n_stream);
	}

	Stream<edev_gpu>::Stream(const dt_int32& n_stream, const dt_shape& shape): Stream()
	{
		resize(n_stream);
		set_grid(shape);
	}

	Stream<edev_gpu>::~Stream()
	{ 
		cleanup();
		n_stream = n_stream_act = 0;
	}

	dt_int32 Stream<edev_gpu>::size() const
	{
		return static_cast<dt_int32>(stream.size());
	}

	void Stream<edev_gpu>::resize(const dt_int32& n_stream)
	{
		this->n_stream = max(1, n_stream);

		cleanup();

		stream.resize(this->n_stream);

		for(auto i = 0; i < this->n_stream; i++)
		{
			cudaStreamCreate(&(stream[i]));
		}

		set_n_stream_act(size());
	}

	cudaStream_t& Stream<edev_gpu>::operator[](const dt_int32& iy)
	{ 
		return stream[iy];
	}

	const cudaStream_t& Stream<edev_gpu>::operator[](const dt_int32& iy) const 
	{ 
		return stream[iy];
	}

	void Stream<edev_gpu>::synchronize()
	{
		cudaDeviceSynchronize();
	}

	void Stream<edev_gpu>::set_n_stream_act(const dt_int32& n_stream_act)
	{
		this->n_stream_act = (n_stream_act<0)?0:min(size(), n_stream_act);
	}

	void Stream<edev_gpu>::set_grid(const dt_shape& shape)
	{
		this->shape = shape;
	}

	void Stream<edev_gpu>::set_n_stream_act_grid(const dt_int32& n_stream_act, const dt_int32& nx)
	{
		set_n_stream_act(n_stream_act);
		shape = dt_shape(nx);
	}

	void Stream<edev_gpu>::set_n_stream_act_grid(const dt_int32& n_stream_act, const dt_int32& nx, const dt_int32& ny)
	{
		set_n_stream_act(n_stream_act);
		shape = dt_shape(ny, nx);
	}

	void Stream<edev_gpu>::set_n_stream_act_grid(const dt_int32& n_stream_act, const dt_int32& nx, const dt_int32& ny, const dt_int32& nz)
	{
		set_n_stream_act(n_stream_act);
		shape = dt_shape(ny, nx, nz);
	}

	template <eDim Dim>
	enable_if_edim_1<Dim, iThread_Rect_xd<Dim>>
	Stream<edev_gpu>::get_ithread_rect_xd(dt_int32 istm)
	{
		iThread_Rect_xd<Dim> ithread_rect;
				
		ithread_rect.istm = istm;
		auto sz = shape.shape_size();
		auto qsiz = sz/n_stream_act;
		ithread_rect.ind_0 = istm*qsiz;
		ithread_rect.ind_e = (istm+1)*qsiz;		

		ithread_rect.ix_0 = istm*qsiz;
		ithread_rect.ix_e = (istm+1)*qsiz;

		if (istm == n_stream_act-1)
		{
			ithread_rect.ind_e += (sz - qsiz*n_stream_act);
			ithread_rect.ix_e += (sz - qsiz*n_stream_act);
		}

		return ithread_rect;
	}

	template <eDim Dim>
	enable_if_edim_2<Dim, iThread_Rect_xd<Dim>>
	Stream<edev_gpu>::get_ithread_rect_xd(dt_int32 istm)
	{
		iThread_Rect_xd<Dim> ithread_rect;
				
		ithread_rect.istm = istm;
		auto sz = shape.shape_size();
		auto qsiz = sz/n_stream_act;
		ithread_rect.ind_0 = istm*qsiz;
		ithread_rect.ind_e = (istm+1)*qsiz;

		auto qnx = shape[1]/n_stream_act;
		ithread_rect.ix_0 = istm*qnx;
		ithread_rect.ix_e = (istm+1)*qnx;
		ithread_rect.iy_0 = 0;
		ithread_rect.iy_e = shape[0];

		if (istm == n_stream_act-1)
		{
			ithread_rect.ix_e += (shape[1] - qnx*n_stream_act);
			ithread_rect.ind_e += (sz - qsiz*n_stream_act);
		}

		return ithread_rect;
	}

	template <eDim Dim>
	enable_if_edim_3<Dim, iThread_Rect_xd<Dim>>
	Stream<edev_gpu>::get_ithread_rect_xd(dt_int32 istm)
	{
		iThread_Rect_xd<Dim> ithread_rect;
				
		ithread_rect.istm = istm;
		auto sz = shape.shape_size();
		auto qsiz = sz/n_stream_act;
		ithread_rect.ind_0 = istm*qsiz;
		ithread_rect.ind_e = (istm+1)*qsiz;

		auto qnz = shape[2]/n_stream_act;
		ithread_rect.iz_0 = istm*qnz;
		ithread_rect.iz_e = (istm+1)*qnz;
		ithread_rect.ix_0 = 0;
		ithread_rect.ix_e = shape[1];
		ithread_rect.iy_0 = 0;
		ithread_rect.iy_e = shape[0];

		if (istm == n_stream_act-1)
		{
			ithread_rect.iz_e += (shape[2] - qsiz*n_stream_act);
			ithread_rect.ind_e += (sz - qsiz*n_stream_act);
		}

		return ithread_rect;
	}

	void Stream<edev_gpu>::cleanup()
	{
		if (stream.empty())
		{
			return;
		}

		cudaDeviceSynchronize();

		for(auto istm = 0; istm < stream.size(); istm++)
		{
			cudaStreamDestroy(stream[istm]);
		}

		stream.clear();
	}
#endif
}
