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

#include "rnd_cpu.h"
#include "stream_cpu.h"

/* template definition */
namespace mt
{
#ifndef RND_XD_DEC
	#define RND_XD_DEC
	template <class T, eDim Dim, eDist Dist, eDev Dev> class Rnd_xd;
#endif
}

/* derived class */
namespace mt
{
	/* uniform distribution */
	template <class T>
	using Rndu_2d_cpu = Rnd_xd<T, edim_2, edist_u, edev_cpu>;	
		
	/* normal distribution */
	template <class T>
	using Rndn_2d_cpu = Rnd_xd<T, edim_2, edist_n, edev_cpu>;		
		
	/* poisson distribution */
	template <class T>
	using Rndp_2d_cpu = Rnd_xd<T, edim_2, edist_p, edev_cpu>;
}

/* 2d distribution */
namespace mt
{
	template <class T, eDist Dist>
	class Rnd_xd<T, edim_2, Dist, edev_cpu>: public Gen_seed
	{
		public:
			using value_type = T;

			static const eDev device = edev_cpu;

			R_2d<T> m_sc;			// scaling
			R_2d<T> m_sft;			// shift
			R_2d<dt_bool> m_b_dim;	// active dimensions

			Rnd_xd(dt_uint64 seed=300183, dt_int32 n_iter=1);

			void set_in_data(const dt_uint64& seed, const R_2d<T>& sc, const R_2d<T>& sft, dt_int32 n_iter=1);

			void set_seed(dt_uint64 seed, dt_int32 n_iter=1);

			void set_sc_sft(const R_2d<T>& sc, const R_2d<T>& sft);

			void set_act_dim(const R_2d<dt_bool>& b_dim);

			R_2d<T> operator()();

			R_2d<T> operator()(const R_2d<T>& sc);

			R_2d<T> operator()(const R_2d<T>& sc, const R_2d<T>& sft);

		private:
			R_2d<T> gen_rnd(const R_2d<T>& sc);

			R_2d<T> gen_rnd(const R_2d<T>& sc, const R_2d<T>& sft);			

			Rnd_Dist<T, Dist, edev_cpu> rnd_x;
			Rnd_Dist<T, Dist, edev_cpu> rnd_y;
	};
}

#include "../src/rnd_2d_cpu.inl"