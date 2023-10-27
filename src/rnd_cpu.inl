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

#include "rnd_cpu.h"

/* class template random distribution */
namespace mt
{
	template <class T, eDist Dist>
	void Rnd_Dist<T, Dist, edev_cpu>::set_seed(const dt_uint64& seed)
	{
		gen.seed(seed);
		rnd.reset();
	}

	template <class T, eDist Dist>
	void Rnd_Dist<T, Dist, edev_cpu>::reset_rnd()
	{
		rnd.reset();
	}

	template <class T, eDist Dist>
	T Rnd_Dist<T, Dist, edev_cpu>::operator()(const T& x)
	{
		return x*rnd(gen);
	}

	template <class T, eDist Dist>
	T Rnd_Dist<T, Dist, edev_cpu>::operator()(const T& x, const T& sft)
	{
		return this->operator()(x) + sft;
	}
}

/* template specialization random distribution Poisson */
namespace mt
{
	template <class T>
	void Rnd_Dist<T, edist_p, edev_cpu>::set_seed(const dt_uint64& seed)
	{
		gen.seed(seed);
		rnd.reset();
	}

	template <class T>
	void Rnd_Dist<T, edist_p, edev_cpu>::reset_rnd()
	{
		rnd.reset();
	}

	template <class T>
	T Rnd_Dist<T, edist_p, edev_cpu>::operator()(const T& x)
	{
		rnd.param(typename Dist_type<edist_p, T>::param_type(x));

		return static_cast<T>(rnd(gen));
	}

	template <class T>
	T Rnd_Dist<T, edist_p, edev_cpu>::operator()(const T& x, const T& sft)
	{
		return this->operator()(x) + sft;
	}
}

/* class generation seed */
namespace mt
{
	Gen_seed::Gen_seed(dt_uint64 seed): m_seed(seed), rnd_u(1) 
	{
		set_seed(seed);
	}

	dt_uint64 Gen_seed::seed_1d(dt_uint64 seed, dt_int32 n_iter)
	{
		set_seed(seed);

		return iseed(n_iter);
	}

	R_2d<dt_uint64> Gen_seed::seed_2d(dt_uint64 seed, dt_int32 n_iter)
	{
		set_seed(seed);

		return {iseed(n_iter), iseed(n_iter)};
	}

	R_3d<dt_uint64> Gen_seed::seed_3d(dt_uint64 seed, dt_int32 n_iter)
	{
		set_seed(seed);

		return {iseed(n_iter), iseed(n_iter), iseed(n_iter)};
	}

	dt_uint64 Gen_seed::operator()()
	{
		return rnd_u(gen_u);
	}

	Vctr_uint64_cpu Gen_seed::operator()(const dt_int32& n_spl)
	{
		Vctr_uint64_cpu vctr;
		vctr.reserve(n_spl);
		for(auto ik=0; ik<n_spl; ik++)
		{
			vctr.push_back(rnd_u(gen_u));
		}

		return vctr;
	}


	void Gen_seed::set_seed(const dt_uint64& seed)
	{
		m_seed = dt_uint64(seed);
		gen_u.seed(m_seed);
		rnd_u.reset();
	}

	/* generate seed by iteration */
	dt_uint64 Gen_seed::iseed(const dt_int32& n_iter)
	{
		dt_uint64 seed = m_seed;

		for(auto ic = 0; ic<n_iter; ic++)
		{
			seed = rnd_u(gen_u);
		}

		return seed;
	}
}

/* traits */
namespace mt
{
	template <eDist Dist, class T>
	using enable_if_dist_p = typename std::enable_if<Dist==edist_p, T>::type;

	template <eDist Dist, class T>
	using enable_ifn_dist_p = typename std::enable_if<Dist!=edist_p, T>::type;
}