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

#include <random>

#include "const_enum.h"
#include "math_mt.h"
#include "r_2d.h"
#include "r_3d.h"
#include "vctr_cpu.h"

/* template definition */
namespace mt
{
	template <class T, eDist Dist, eDev Dev> class Rnd_Dist;
}

/* derived class */
namespace mt
{
	template <eDist Dist, class T>
	using Dist_type = typename std::conditional<Dist==edist_u, std::uniform_real_distribution<T>, 
	typename std::conditional<Dist==edist_n, std::normal_distribution<T>, \
	typename std::conditional<Dist==edist_p, std::poisson_distribution<dt_int32>, T>::type>::type>::type;
}	

/* class template random distribution */
namespace mt
{
	template <class T, eDist Dist> 
	class Rnd_Dist<T, Dist, edev_cpu>
	{
	 	public:
	 		using value_type = T;

			void set_seed(const dt_uint64& seed);

			void reset_rnd();

			T operator()(const T& x);

			T operator()(const T& x, const T& sft);
		private:
			std::mt19937_64 gen;
			Dist_type<Dist, T> rnd;
	};
}

/* template specialization random distribution Poisson */
namespace mt
{
	template <class T> 
	class Rnd_Dist<T, edist_p, edev_cpu>
	{
	 	public:
	 		using value_type = T;

			void set_seed(const dt_uint64& seed);

			void reset_rnd();

			T operator()(const T& x);

			T operator()(const T& x, const T& sft);

		private:
			std::mt19937_64 gen;
			Dist_type<edist_p, T> rnd;
	};
}

/* class generation seed */
namespace mt
{
	class Gen_seed
	{
		public:
			dt_uint64 m_seed;	// seed

			Gen_seed(dt_uint64 seed=300183);

			dt_uint64 seed_1d(dt_uint64 seed, dt_int32 n_iter=1);

			R_2d<dt_uint64> seed_2d(dt_uint64 seed, dt_int32 n_iter=1);

			R_3d<dt_uint64> seed_3d(dt_uint64 seed, dt_int32 n_iter=1);

			dt_uint64 operator()();

			Vctr_uint64_cpu operator()(const dt_int32& n_spl);

		private:
			void set_seed(const dt_uint64& seed);

			/* generate seed by iteration */
			dt_uint64 iseed(const dt_int32& n_iter);

			std::mt19937_64 gen_u;
			std::uniform_int_distribution<dt_uint64> rnd_u;
	};
}

#include "detail/rnd_cpu.inl"