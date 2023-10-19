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

#include "rnd_1d_cpu.h"

/* 1d distribution */
namespace mt
{
	template <class T, eDist Dist>
	Rnd_xd<T, edim_1, Dist, edev_cpu>::Rnd_xd(dt_uint64 seed, dt_int32 n_iter): Gen_seed(seed), m_sc{1}, m_sft{0}
	{
		set_seed(seed, n_iter);
	}

	template <class T, eDist Dist>
	void Rnd_xd<T, edim_1, Dist, edev_cpu>::set_in_data(const dt_uint64& seed, const T& sc, const T& sft, dt_int32 n_iter)
	{
		set_seed(seed, n_iter);
		set_sc_sft(sc, sft);
	}

	template <class T, eDist Dist>
	void Rnd_xd<T, edim_1, Dist, edev_cpu>::set_seed(dt_uint64 seed, dt_int32 n_iter)
	{
		auto gen_seed = this->seed_1d(seed, n_iter);

		rnd.set_seed(gen_seed);
	}

	template <class T, eDist Dist>
	void Rnd_xd<T, edim_1, Dist, edev_cpu>::set_sc_sft(const T& sc, const T& sft)
	{
		m_sc = sc;
		m_sft = sft;
	}

	template <class T, eDist Dist>
	T Rnd_xd<T, edim_1, Dist, edev_cpu>::operator()()
	{
		return rnd(m_sc, m_sft);
	}

	template <class T, eDist Dist>
	T Rnd_xd<T, edim_1, Dist, edev_cpu>::operator()(const T& sc)
	{
		return rnd(sc);
	}

	template <class T, eDist Dist>
	T Rnd_xd<T, edim_1, Dist, edev_cpu>::operator()(const T& sc, const T& sft)
	{
		return rnd(sc, sft);
	}
}

/* fcns */
namespace mt
{
	// add uniform noise
	template <class TVctr_i, class TVctr_o>
	enable_if_vctr_cpu_and_vctr_cpu<TVctr_i, TVctr_o, void>
	fcn_add_unif_nois(TVctr_i& mx_i, Value_type<TVctr_o> sc, dt_uint64 seed_i, TVctr_o& mx_o, Stream_cpu* pstream)
	{	
		using T = Value_type<TVctr_o>;
		auto seed = seed_i;

		if (seed<1)
		{
			std::random_device rd;
			seed = rd();
		}

		Gen_seed gen_seed(seed);
		Vctr_uint64_cpu vctr_seed = gen_seed((pstream==nullptr)?1:pstream->size());

		const T a = -sc;
		const T b = T(2)*sc;

		auto thr_add_nois = [a, b, &mx_i, &vctr_seed, &mx_o](const iThread_Rect_1d& range)
		{
			Rndu_1d_cpu<T> rnd(vctr_seed[range.istm]);

			for(auto ind = range.ind_0; ind < range.ind_e; ind++)
			{
				mx_o[ind] = rnd(b, T(mx_i[ind]) + a);
			}
		};

		fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_add_nois);
	};

	// add gaussian noise
	template <class TVctr_i, class TVctr_o>
	enable_if_vctr_cpu_and_vctr_cpu<TVctr_i, TVctr_o, void>
	fcn_add_gauss_nois(TVctr_i& mx_i, Value_type<TVctr_o> sigma, dt_uint64 seed_i, TVctr_o& mx_o, Stream_cpu* pstream)
	{	
		using T = Value_type<TVctr_o>;
		auto seed = seed_i;

		if (seed<1)
		{
			std::random_device rd;
			seed = rd();
		}

		Gen_seed gen_seed(seed);
		Vctr_uint64_cpu vctr_seed = gen_seed((pstream==nullptr)?1:pstream->size());

		auto thr_add_nois = [sigma, &mx_i, &vctr_seed, &mx_o](const iThread_Rect_1d& range)
		{
			Rndn_1d_cpu<T> rnd(vctr_seed[range.istm]);

			for(auto ind = range.ind_0; ind < range.ind_e; ind++)
			{
				mx_o[ind] = rnd(sigma, T(mx_i[ind]));
			}
		};

		fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_add_nois);
	};


	// add poisson noise
	template <class TVctr_i, class TVctr_o>
	enable_if_vctr_cpu_and_vctr_cpu<TVctr_i, TVctr_o, void>
	fcn_add_poiss_nois(TVctr_i& mx_i, double sc_i, dt_uint64 seed_i, TVctr_o& mx_o, Stream_cpu* pstream)
	{	
		using T = Value_type<TVctr_o>;
		auto seed = seed_i;
		auto sc = T(sc_i);

		if (seed<1)
		{
			std::random_device rd;
			seed = rd();
		}

		Gen_seed gen_seed(seed);
		Vctr_uint64_cpu vctr_seed = gen_seed((pstream==nullptr)?1:pstream->size());

		auto thr_add_nois = [sc, &mx_i, &vctr_seed, &mx_o](const iThread_Rect_1d& range)
		{
			Rndp_1d_cpu<T> rnd(vctr_seed[range.istm]);

			for(auto ind = range.ind_0; ind < range.ind_e; ind++)
			{
				mx_o[ind] = rnd(T(mx_i[ind])*sc);
			}
		};

		fcn_stream_exec_xd_fcn<edim_1>(pstream, mx_o.size_32(), thr_add_nois);
	};

	//// add Poisson noise
	//template <class TVctr>
	//TVctr add_poiss_nois_by_SNR(Stream<edev_cpu>& stream, TVctr& Im_i, 
	//	Value_type<TVctr> SNR_i, Value_type<TVctr>& scl_o)
	//{	
	//	using T = Value_type<TVctr>;

	//	auto get_SNR = [](Stream<edev_cpu>& stream, TVctr& Im, T Im_std, T scf)->T
	//	{
	//		T x_mean = 0;
	//		T x_var = 0;
	//		auto thr_SNR = [&](const iThread_Rect_2d& range)
	//		{
	//			std::mt19937_64 gen;
	//			std::poisson_distribution<dt_int32> rand;

	//			T x_mean_partial = 0;
	//			T x_var_partial = 0;
	//			for(auto ixy = range.ind_0; ixy < range.ind_e; ixy++)
	//			{
	//				auto x0 = scf*Im[ixy];
	//				rand.param(std::poisson_distribution<dt_int32>::param_type(x0));
	//				auto xn = rand(gen)-x0;
	//				x_mean_partial += xn;
	//				x_var_partial += xn*xn;
	//			}

	//			stream.stream_mutex.lock();
	//			x_mean += x_mean_partial;
	//			x_var += x_var_partial;
	//			stream.stream_mutex.unlock();
	//		};

	//		stream.set_n_stream_act(Im.size());
	//		stream.set_grid(1, Im.size());
	//		stream.exec(thr_SNR);

	//		x_mean /= Im.size();
	//		x_var = x_var/Im.size()-x_mean*x_mean;

	//		return scf*Im_std/sqrt(x_var);
	//	};

	//	auto Im_i_std = sqrt(fcn_variance(stream, Im_i));

	//	T SNR_k = get_SNR(stream, Im_i, Im_i_std, 1);
	
	//	T k_0 = 1;
	//	T k_e = 1;

	//	if (SNR_k<SNR_i)
	//	{
	//		do
	//		{
	//			k_e *= 2;
	//			SNR_k = get_SNR(stream, Im_i, Im_i_std, k_e);
	//		}
	//		while (SNR_k < SNR_i);
	//		k_0 = k_e/2;
	//	}
	//	else
	//	{
	//		do
	//		{
	//			k_0 /= 2;
	//			SNR_k = get_SNR(stream, Im_i, Im_i_std, k_0);
	//		}
	//		while (SNR_k >= SNR_i);
	//		k_e = 2*k_0;
	//	}


	//	// bisection method
	//	dt_int32 ic = 0;
	//	do
	//	{
	//		scl_o = 0.5*(k_0 + k_e);
	//		auto SNR_k = get_SNR(stream, Im_i, Im_i_std, scl_o);

	//		if (SNR_k < SNR_i)
	//		{
	//			k_0 = scl_o;
	//		}
	//		else
	//		{
	//			k_e = scl_o;
	//		}
	//		ic++;
	//	} while ((fabs(SNR_i-SNR_k)>0.1) && (ic<10));

	//	// add Poisson noise
	//	return fcn_add_poiss_nois(stream, Im_i, scl_o);
	//}

}