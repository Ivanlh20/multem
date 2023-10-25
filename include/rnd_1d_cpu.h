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
	using Rndu_1d_cpu = Rnd_xd<T, edim_1, edist_u, edev_cpu>;	
		
	/* normal distribution */
	template <class T>
	using Rndn_1d_cpu = Rnd_xd<T, edim_1, edist_n, edev_cpu>;	
		
	/* poisson distribution */
	template <class T>
	using Rndp_1d_cpu = Rnd_xd<T, edim_1, edist_p, edev_cpu>;
}

/* template specialization 1d distribution */
namespace mt
{
	template <class T, eDist Dist>
	class Rnd_xd<T, edim_1, Dist, edev_cpu>: public Gen_seed
	{
		public:
			using value_type = T;

			static const eDev device = edev_cpu;

			T m_sc;			// scaling
			T m_sft;		// shift

			Rnd_xd(dt_uint64 seed=300183, dt_int32 n_iter=1);

			void set_in_data(const dt_uint64& seed, const T& sc, const T& sft, dt_int32 n_iter=1);

			void set_seed(dt_uint64 seed, dt_int32 n_iter=1);

			void set_sc_sft(const T& sc, const T& sft);

			T operator()();

			T operator()(const T& sc);

			T operator()(const T& sc, const T& sft);

		private:
			Rnd_Dist<T, Dist, edev_cpu> rnd;
	};
}

/* fcns */
namespace mt
{
	// add uniform noise
	template <class TVctr_i, class TVctr_o>
	enable_if_vctr_cpu_and_vctr_cpu<TVctr_i, TVctr_o, void>
	fcn_add_unif_nois(TVctr_i& mx_i, Value_type<TVctr_o> sc, dt_uint64 seed_i, TVctr_o& mx_o, Stream_cpu* pstream = nullptr);

	// add gaussian noise
	template <class TVctr_i, class TVctr_o>
	enable_if_vctr_cpu_and_vctr_cpu<TVctr_i, TVctr_o, void>
	fcn_add_gauss_nois(TVctr_i& mx_i, Value_type<TVctr_o> sigma, dt_uint64 seed_i, TVctr_o& mx_o, Stream_cpu* pstream = nullptr);


	// add poisson noise
	template <class TVctr_i, class TVctr_o>
	enable_if_vctr_cpu_and_vctr_cpu<TVctr_i, TVctr_o, void>
	fcn_add_poiss_nois(TVctr_i& mx_i, double sc_i, dt_uint64 seed_i, TVctr_o& mx_o, Stream_cpu* pstream = nullptr);

	//// add Poisson noise
	//template <class TVctr>
	//TVctr add_poiss_nois_by_SNR(Stream<edev_cpu>& stream, TVctr& Im_i, Value_type<TVctr> SNR_i, Value_type<TVctr>& scl_o);
}

#include "../src/rnd_1d_cpu.inl"