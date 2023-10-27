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

#include "rnd_2d_cpu.h"

/* 2d distribution */
namespace mt
{
	template <class T, eDist Dist>
	Rnd_xd<T, edim_2, Dist, edev_cpu>::Rnd_xd(dt_uint64 seed, dt_int32 n_iter): Gen_seed(seed), 
		m_sc{1, 1}, m_b_dim{true, true}
	{
		set_seed(seed, n_iter);
	}

	template <class T, eDist Dist>
	void Rnd_xd<T, edim_2, Dist, edev_cpu>::set_in_data(const dt_uint64& seed, const R_2d<T>& sc, const R_2d<T>& sft, dt_int32 n_iter)
	{
		set_seed(seed, n_iter);
		set_sc_sft(sc, sft);
	}

	template <class T, eDist Dist>
	void Rnd_xd<T, edim_2, Dist, edev_cpu>::set_seed(dt_uint64 seed, dt_int32 n_iter)
	{
		auto gen_seed = this->seed_2d(seed, n_iter);

		rnd_x.set_seed(gen_seed.x);
		rnd_y.set_seed(gen_seed.y);
	}

	template <class T, eDist Dist>
	void Rnd_xd<T, edim_2, Dist, edev_cpu>::set_sc_sft(const R_2d<T>& sc, const R_2d<T>& sft)
	{
		m_sc = sc;
		m_sft = sft;
	}

	template <class T, eDist Dist>
	void Rnd_xd<T, edim_2, Dist, edev_cpu>::set_act_dim(const R_2d<dt_bool>& b_dim)
	{
		m_b_dim = b_dim;
	}

	template <class T, eDist Dist>
	R_2d<T> Rnd_xd<T, edim_2, Dist, edev_cpu>::operator()()
	{
		return this->operator()(m_sc, m_sft);
	}

	template <class T, eDist Dist>
	R_2d<T> Rnd_xd<T, edim_2, Dist, edev_cpu>::operator()(const R_2d<T>& sc)
	{
		return gen_rnd(sc);
	}

	template <class T, eDist Dist>
	R_2d<T> Rnd_xd<T, edim_2, Dist, edev_cpu>::operator()(const R_2d<T>& sc, const R_2d<T>& sft)
	{
		return gen_rnd(sc, sft);
	}

	template <class T, eDist Dist>
	R_2d<T> Rnd_xd<T, edim_2, Dist, edev_cpu>::gen_rnd(const R_2d<T>& sc) 
	{
		return {(m_b_dim.x)?(rnd_x(sc.x)):T(0), (m_b_dim.y)?(rnd_y(sc.y)):T(0)};
	}

	template <class T, eDist Dist>
	R_2d<T> Rnd_xd<T, edim_2, Dist, edev_cpu>::gen_rnd(const R_2d<T>& sc, const R_2d<T>& sft) 
	{
		return {(m_b_dim.x)?(rnd_x(sc.x, sft.x)):T(0), (m_b_dim.y)?(rnd_y(sc.y, sft.y)):T(0)};
	}
}