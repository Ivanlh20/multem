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

#include "stream_cpu.h"

/* template specialization cpu stream */
namespace mt
{
	Stream<edev_cpu>::Stream(): shape(), n_stream(0), n_stream_act(0), stream(nullptr) {}

	Stream<edev_cpu>::Stream(const dt_int32& n_stream): Stream()
	{
		resize(n_stream);
	}

	Stream<edev_cpu>::Stream(const dt_int32& n_stream, const dt_shape& shape): Stream()
	{
		resize(n_stream);
		set_grid(shape);
	}

	Stream<edev_cpu>::~Stream(){ 
		cleanup();
		n_stream = n_stream_act = 0;
	}

	dt_int32 Stream<edev_cpu>::size() const
	{
		return n_stream;
	}

	dt_bool Stream<edev_cpu>::is_single_thread() const
	{
		return n_stream < 2;
	}

	void Stream<edev_cpu>::resize(const dt_int32& n_stream)
	{
		cleanup();

		this->n_stream = (n_stream<1)?1:n_stream;

		if (this->n_stream > 1)
		{
			stream = new std::thread[this->n_stream-1];
		}

		set_n_stream_act(size());
	}

	std::thread& Stream<edev_cpu>::operator[](const dt_int32& iy)
	{ 
		return stream[iy];
	}

	const std::thread& Stream<edev_cpu>::operator[](const dt_int32& iy) const 
	{ 
		return stream[iy];
	}

	void Stream<edev_cpu>::synchronize()
	{
		cleanup();

		if (n_stream > 1)
		{
			stream = new std::thread[n_stream-1];
		}
	}

	void Stream<edev_cpu>::set_n_stream_act(const dt_int32& n_stream_act)
	{
		this->n_stream_act = (n_stream_act<0)?0:min(size(), n_stream_act);
	}

	void Stream<edev_cpu>::set_grid(const dt_shape& shape)
	{
		this->shape = shape;
	}

	void Stream<edev_cpu>::set_n_stream_act_grid(const dt_int32& n_stream_act, const dt_int32& nx)
	{
		set_n_stream_act(n_stream_act);
		shape = dt_shape(nx);
	}

	void Stream<edev_cpu>::set_n_stream_act_grid(const dt_int32& n_stream_act, const dt_int32& nx, const dt_int32& ny)
	{
		set_n_stream_act(n_stream_act);
		shape = dt_shape(ny, nx);
	}

	void Stream<edev_cpu>::set_n_stream_act_grid(const dt_int32& n_stream_act, const dt_int32& nx, const dt_int32& ny, const dt_int32& nz)
	{
		set_n_stream_act(n_stream_act);
		shape = dt_shape(ny, nx, nz);
	}

	template <eDim Dim>
	enable_if_edim_1<Dim, iThread_Rect_xd<Dim>>
	Stream<edev_cpu>::get_ithread_rect_xd(dt_int32 istm)
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
	Stream<edev_cpu>::get_ithread_rect_xd(dt_int32 istm)
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
	Stream<edev_cpu>::get_ithread_rect_xd(dt_int32 istm)
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

	/***************************************************************************************/
	template <eDim Dim, class TFcn, class... TArgs>
	enable_if_edim_1<Dim, void>
	Stream<edev_cpu>::exec_xd_krn(const dt_int32& nx, TFcn& fcn, TArgs& ...arg)
	{
		set_n_stream_act_grid(min(size(), nx), nx);

		if (n_stream_act < 1)
		{
			return;
		}

		auto thr_fcn = [&](const iThread_Rect_xd<Dim>& ithread)
		{
			for(auto ind = ithread.ind_0; ind < ithread.ind_e; ind++)
			{
				fcn(ind, arg...);
			}
		};

		for(auto istm = 0; istm < n_stream_act-1; istm++)
		{
			stream[istm] = std::thread(thr_fcn, std::move(get_ithread_rect_xd<Dim>(istm)));
		}

		thr_fcn(std::move(get_ithread_rect_xd<Dim>(n_stream_act-1)));

		synchronize();
	}				

	template <eDim Dim, class TFcn, class... TArgs>
	enable_if_edim_2<Dim, void>
	Stream<edev_cpu>::exec_xd_krn(const dt_int32& nx, const dt_int32& ny, TFcn& fcn, TArgs& ...arg)
	{
		set_n_stream_act_grid(min(size(), nx), nx, ny);

		if (n_stream_act < 1)
		{
			return;
		}

		auto thr_fcn = [&](const iThread_Rect_xd<Dim>& ithread)
		{
			for(auto ix = ithread.ix_0; ix < ithread.ix_e; ix++)
			{
				for(auto iy = ithread.iy_0; iy < ithread.iy_e; iy++)
				{
					fcn(ix, iy, arg...);
				}
			}
		};

		for(auto istm = 0; istm < n_stream_act-1; istm++)
		{
			stream[istm] = std::thread(thr_fcn, std::move(get_ithread_rect_xd<Dim>(istm)));
		}

		thr_fcn(std::move(get_ithread_rect_xd<Dim>(n_stream_act-1)));

		synchronize();
	}

	template <eDim Dim, class TFcn, class... TArgs>
	enable_if_edim_3<Dim, void>
	Stream<edev_cpu>::exec_xd_krn(const dt_int32& nx, const dt_int32& ny, const dt_int32& nz, TFcn& fcn, TArgs& ...arg)
	{
		set_n_stream_act_grid(min(size(), nz), nx, ny, nz);

		if (n_stream_act < 1)
		{
			return;
		}

		auto thr_fcn = [&](const iThread_Rect_xd<Dim>& ithread)
		{
			for(auto iz = ithread.iz_0; iz < ithread.iz_e; iz++)
			{
				for(auto ix = ithread.ix_0; ix < ithread.ix_e; ix++)
				{
					for(auto iy = ithread.iy_0; iy < ithread.iy_e; iy++)
					{
						fcn(ix, iy, iz, arg...);
					}
				}
			}
		};

		for(auto istm = 0; istm < n_stream_act-1; istm++)
		{
			stream[istm] = std::thread(thr_fcn, std::move(get_ithread_rect_xd<Dim>(istm)));
		}

		thr_fcn(std::move(get_ithread_rect_xd<Dim>(n_stream_act-1)));

		synchronize();
	}

	/***************************************************************************************/
	template <eDim Dim, class TFcn, class... TArgs>
	enable_if_edim_1<Dim, void>
	Stream<edev_cpu>::exec_xd_fcn(const dt_int32& nx, TFcn& fcn, TArgs& ...arg)
	{
		set_n_stream_act_grid(min(this->size(), nx), nx);

		if (n_stream_act < 1)
		{
			return;
		}

		for(auto istm = 0; istm < n_stream_act-1; istm++)
		{
			stream[istm] = std::thread(std::bind(fcn, get_ithread_rect_xd<Dim>(istm), std::ref<TArgs>(arg)...));
		}

		fcn(get_ithread_rect_xd<Dim>(n_stream_act-1), std::ref<TArgs>(arg)...);

		synchronize();
	}

	template <eDim Dim, class TFcn, class... TArgs>
	enable_if_edim_2<Dim, void>
	Stream<edev_cpu>::exec_xd_fcn(const dt_int32& nx, const dt_int32& ny, TFcn& fcn, TArgs& ...arg)
	{
		set_n_stream_act_grid(min(size(), nx), nx, ny);

		if (n_stream_act < 1)
		{
			return;
		}

		for(auto istm = 0; istm < n_stream_act-1; istm++)
		{
			stream[istm] = std::thread(std::bind(fcn, get_ithread_rect_xd<Dim>(istm), std::ref<TArgs>(arg)...));
		}

		fcn(get_ithread_rect_xd<Dim>(n_stream_act-1), std::ref<TArgs>(arg)...);

		synchronize();
	}

	template <eDim Dim, class TFcn, class... TArgs>
	enable_if_edim_3<Dim, void>
	Stream<edev_cpu>::exec_xd_fcn(const dt_int32& nx, const dt_int32& ny, const dt_int32& nz, TFcn& fcn, TArgs& ...arg)
	{
		set_n_stream_act_grid(min(size(), nz), nx, ny, nz);

		if (n_stream_act < 1)
		{
			return;
		}

		for(auto istm = 0; istm < n_stream_act-1; istm++)
		{
			stream[istm] = std::thread(std::bind(fcn, get_ithread_rect_xd<Dim>(istm), std::ref<TArgs>(arg)...));
		}

		fcn(get_ithread_rect_xd<Dim>(n_stream_act-1), std::ref<TArgs>(arg)...);

		synchronize();
	}

 	void Stream<edev_cpu>::cleanup()
	{
		if ((n_stream <= 0) || (stream==nullptr))
		{
			return;
		}

		for(auto i = 0; i < n_stream-1; i++)
		{
			if (stream[i].joinable())
			{
				stream[i].join();
			}
		}

		delete [] stream;
		stream = nullptr;
	};
}
	
/* stream functions */
namespace mt
{
	dt_bool fcn_is_single_thread(Stream_cpu* pstream)
	{
		return (pstream == nullptr)?true:pstream->is_single_thread();
	}
}

/* kernel execution functions */
namespace mt
{
	template <eDim Dim, class TFcn, class... TArgs>
	enable_if_edim_1<Dim, void>
	fcn_stream_exec_xd_krn(Stream_cpu* pstream, const dt_int32& nx, TFcn& fcn, TArgs&& ...arg)
	{
		if (pstream == nullptr)
		{
			Stream_cpu stream(1);
			stream.exec_xd_krn<Dim>(nx, fcn, arg...);
		}
		else
		{
			pstream->exec_xd_krn<Dim>(nx, fcn, arg...);
		}
	}

	template <eDim Dim, class TFcn, class... TArgs>
	enable_if_edim_2<Dim, void>
	fcn_stream_exec_xd_krn(Stream_cpu* pstream, const dt_int32& nx, const dt_int32& ny, TFcn& fcn, TArgs&& ...arg)
	{
		if (pstream == nullptr)
		{
			Stream_cpu stream(1);
			stream.exec_xd_krn<Dim>(nx, ny, fcn, arg...);
		}
		else
		{
			pstream->exec_xd_krn<Dim>(nx, ny, fcn, arg...);
		}
	}

	template <eDim Dim, class TFcn, class... TArgs>
	enable_if_edim_3<Dim, void>
	fcn_stream_exec_xd_krn(Stream_cpu* pstream, const dt_int32& nx, const dt_int32& ny, const dt_int32& nz, TFcn& fcn, TArgs&& ...arg)
	{
		if (pstream == nullptr)
		{
			Stream_cpu stream(1);
			stream.exec_xd_krn<Dim>(nx, ny, nz, fcn, arg...);
		}
		else
		{
			pstream->exec_xd_krn<Dim>(nx, ny, nz, fcn, arg...);
		}
	}
}

/* kernel reduction execution functions */
namespace mt
{
	template <eDim Dim, class T, class TFcn, class TOp, class... TArgs>
	enable_if_edim_1<Dim, T>
	fcn_stream_exec_xd_krn_reduce(Stream_cpu* pstream, const dt_int32& nx, TOp&& fcn_op, T v_0, TFcn& fcn, TArgs&& ...arg)
	{
		auto thr_fcn = [](const iThread_Rect_xd<Dim>& ithread, TFcn& fcn, TArgs&& ...arg, std::vector<T>& vctr)
		{
			KS<T> v_part = vctr[ithread.istm];

			for(auto ind = ithread.ind_0; ind < ithread.ind_e; ind++)
			{
				fcn(ind, arg..., v_part);
			}

			vctr[ithread.istm] = v_part;
		};

		T v_tot = v_0;

		if (fcn_is_single_thread(pstream))
		{
			std::vector<T> vctr(1, v_0);
			thr_fcn(iThread_Rect_xd<Dim>(nx), fcn, arg..., vctr);
			v_tot = vctr[0];
		}
		else
		{
			std::vector<T> vctr(pstream->size(), v_0);
			pstream->exec_xd_fcn<Dim>(nx, thr_fcn, fcn, arg..., vctr);

			// run remainder operations
			for(auto ik=0; ik<vctr.size(); ik++)
			{
				v_tot = fcn_op(v_tot, vctr[ik]);
			}
		}

		return v_tot;
	}

	template <eDim Dim, class T, class TFcn, class TOp, class... TArgs>
	enable_if_edim_2<Dim, T>
	fcn_stream_exec_xd_krn_reduce(Stream_cpu* pstream, const dt_int32& nx, const dt_int32& ny, TOp&& fcn_op, T v_0, TFcn& fcn, TArgs&& ...arg)
	{
		auto thr_fcn = [](const iThread_Rect_xd<Dim>& ithread, TFcn& fcn, TArgs&& ...arg, std::vector<T>& vctr)
		{
			KS<T> v_part = vctr[ithread.istm];

			for(auto ix = ithread.ix_0; ix < ithread.ix_e; ix++)
			{
				for(auto iy = ithread.iy_0; iy < ithread.iy_e; iy++)
				{
					fcn(ix, iy, arg..., v_part);
				}
			}

			vctr[ithread.istm] = v_part;
		};

		T v_tot = v_0;

		if (fcn_is_single_thread(pstream))
		{
			std::vector<T> vctr(1, v_0);
			thr_fcn(iThread_Rect_xd<Dim>(nx, ny), fcn, arg..., vctr);
			v_tot = vctr[0];
		}
		else
		{
			std::vector<T> vctr(pstream->size(), v_0);
			pstream->exec_xd_fcn<Dim>(nx, ny, thr_fcn, fcn, arg..., vctr);

			// run remainder operations
			for(auto ik=0; ik<vctr.size(); ik++)
			{
				v_tot = fcn_op(v_tot, vctr[ik]);
			}
		}

		return v_tot;
	}

	template <eDim Dim, class T, class TFcn, class TOp, class... TArgs>
	enable_if_edim_3<Dim, T>
	fcn_stream_exec_xd_krn_reduce(Stream_cpu* pstream, const dt_int32& nx, const dt_int32& ny, const dt_int32& nz, TOp&& fcn_op, T v_0, TFcn& fcn, TArgs&& ...arg)
	{
		auto thr_fcn = [](const iThread_Rect_xd<Dim>& ithread, TFcn& fcn, TArgs&& ...arg, std::vector<T>& vctr)
		{
			KS<T> v_part = vctr[ithread.istm];

			for(auto iz = ithread.iz_0; iz < ithread.iz_e; iz++)
			{
				for(auto ix = ithread.ix_0; ix < ithread.ix_e; ix++)
				{
					for(auto iy = ithread.iy_0; iy < ithread.iy_e; iy++)
					{
						fcn(ix, iy, iz, arg..., v_part);
					}
				}
			}

			vctr[ithread.istm] = v_part;
		};

		T v_tot = v_0;

		if (fcn_is_single_thread(pstream))
		{
			std::vector<T> vctr(1, v_0);
			thr_fcn(iThread_Rect_xd<Dim>(nx, ny, nz), fcn, arg..., vctr);
			v_tot = vctr[0];
		}
		else
		{
			std::vector<T> vctr(pstream->size(), v_0);
			pstream->exec_xd_fcn<Dim>(nx, ny, nz, thr_fcn, fcn, arg..., vctr);

			// run remainder operations
			for(auto ik=0; ik<vctr.size(); ik++)
			{
				v_tot = fcn_op(v_tot, vctr[ik]);
			}
		}

		return v_tot;
	}
}

/* function execution functions */
namespace mt
{
	template <eDim Dim, class TFcn, class... TArgs>
	enable_if_edim_1<Dim, void>
	fcn_stream_exec_xd_fcn(Stream_cpu* pstream, const dt_int32& nx, TFcn& fcn, TArgs&& ...arg)
	{
		if (fcn_is_single_thread(pstream))
		{
			fcn(iThread_Rect_xd<Dim>(nx), std::ref<TArgs>(arg)...);
		}
		else
		{
			pstream->exec_xd_fcn<Dim>(nx, fcn, arg...);
		}
	}

	template <eDim Dim, class TFcn, class... TArgs>
	enable_if_edim_2<Dim, void>
	fcn_stream_exec_xd_fcn(Stream_cpu* pstream, const dt_int32& nx, const dt_int32& ny, TFcn& fcn, TArgs&& ...arg)
	{
		if (fcn_is_single_thread(pstream))
		{
			fcn(iThread_Rect_xd<Dim>(nx, ny), std::ref<TArgs>(arg)...);
		}
		else
		{
			pstream->exec_xd_fcn<Dim>(nx, ny, fcn, arg...);
		}
	}
		
	template <eDim Dim, class TFcn, class... TArgs>
	enable_if_edim_3<Dim, void>
	fcn_stream_exec_xd_fcn(Stream_cpu* pstream, const dt_int32& nx, const dt_int32& ny, const dt_int32& nz, TFcn& fcn, TArgs&& ...arg)
	{
		if (fcn_is_single_thread(pstream))
		{
			fcn(iThread_Rect_xd<Dim>(nx, ny, nz), std::ref<TArgs>(arg)...);
		}
		else
		{
			pstream->exec_xd_fcn<Dim>(nx, ny, nz, fcn, arg...);
		}
	}
}