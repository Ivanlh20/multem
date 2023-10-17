/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef CGPU_RAND_H
	#define CGPU_RAND_H

	#include <random>

	#include "macros.cuh"
	#include "const_enum.cuh"
	#include "math.cuh"
	#include "r_2d.cuh"
	#include "r_3d.cuh"
	#include "cgpu_vctr.cuh"
	#include "cgpu_stream.cuh"

	#ifdef __CUDACC__
		#include <cuda.h>
		#include <cuda_runtime.h>
	#endif

	/* traits */
	namespace mt
	{
		template <eDist Dist, class T>
		using enable_if_dist_p = typename std::enable_if<Dist==edist_p, T>::type;

		template <eDist Dist, class T>
		using enable_ifn_dist_p = typename std::enable_if<Dist!=edist_p, T>::type;

		template <eDist Dist, class T>
	 	using Dist_type = typename std::conditional<Dist==edist_u, std::uniform_real_distribution<T>, 
		typename std::conditional<Dist==edist_n, std::normal_distribution<T>, \
		typename std::conditional<Dist==edist_p, std::poisson_distribution<dt_int32>, T>::type>::type>::type;
	}

	/* distribution type */
	namespace mt
	{	
		template <class T, eDist Dist, eDev Dev> class Rnd_Dist;

		template <class T, eDist Dist> 
		class Rnd_Dist<T, Dist, edev_cpu>
		{
	 		public:
	 			using value_type = T;

				void set_seed(const dt_uint64& seed)
				{
					gen.seed(seed);
					rnd.reset();
				}

				void reset_rnd()
				{
					rnd.reset();
				}

				T operator()(const T& x)
				{
					return x*rnd(gen);
				}

				T operator()(const T& x, const T& sft)
				{
					return this->operator()(x) + sft;
				}
			private:
				std::mt19937_64 gen;
				Dist_type<Dist, T> rnd;
		};

		template <class T> 
		class Rnd_Dist<T, edist_p, edev_cpu>
		{
	 		public:
	 			using value_type = T;

				void set_seed(const dt_uint64& seed)
				{
					gen.seed(seed);
					rnd.reset();
				}

				void reset_rnd()
				{
					rnd.reset();
				}

				T operator()(const T& x)
				{
					rnd.param(Dist_type<edist_p, T>::param_type(x));

					return T(rnd(gen));
				}

				T operator()(const T& x, const T& sft)
				{
					return this->operator()(x) + sft;
				}

			private:
				std::mt19937_64 gen;
				Dist_type<edist_p, T> rnd;
		};
	}

	/* distribution */
	namespace mt
	{
		class Gen_seed
		{
			public:
				dt_uint64 m_seed;	// seed

				Gen_seed(dt_uint64 seed=300183): m_seed(seed), rand_u(1) 
				{
					set_seed(seed);
				}

				dt_uint64 seed_1d(dt_uint64 seed, dt_int32 n_iter=1)
				{
					set_seed(seed);

					return iseed(n_iter);
				}

				R_2d<dt_uint64> seed_2d(dt_uint64 seed, dt_int32 n_iter=1)
				{
					set_seed(seed);

					return {iseed(n_iter), iseed(n_iter)};
				}

				R_3d<dt_uint64> seed_3d(dt_uint64 seed, dt_int32 n_iter=1)
				{
					set_seed(seed);

					return {iseed(n_iter), iseed(n_iter), iseed(n_iter)};
				}

				dt_uint64 operator()()
				{
					return rand_u(gen_u);
				}

				Vctr_uint64_cpu operator()(const dt_int32& n_spl)
				{
					Vctr_uint64_cpu vctr;
					vctr.reserve(n_spl);
					for(auto ik=0; ik<n_spl; ik++)
					{
						vctr.push_back(rand_u(gen_u));
					}

					return vctr;
				}

			private:

				void set_seed(const dt_uint64& seed)
				{
					m_seed = dt_uint64(seed);
					gen_u.seed(m_seed);
					rand_u.reset();
				}

				/* generate seed by iteration */
				dt_uint64 iseed(const dt_int32& n_iter)
				{
					dt_uint64 seed = m_seed;

					for(auto ic = 0; ic<n_iter; ic++)
					{
						seed = rand_u(gen_u);
					}

					return seed;
				}

				std::mt19937_64 gen_u;
				std::uniform_int_distribution<dt_uint64> rand_u;
		};
	
		template <class T, eDim Dim, eDist Dist, eDev Dev> class Rand_xd;

		/* uniform distribution */
		template <class T>
		using Randu_1d_cpu = Rand_xd<T, edim_1, edist_u, edev_cpu>;		
		
		template <class T>
		using Randu_2d_cpu = Rand_xd<T, edim_2, edist_u, edev_cpu>;		

		template <class T>
		using Randu_3d_cpu = Rand_xd<T, edim_3, edist_u, edev_cpu>;		
		
		/* normal distribution */
		template <class T>
		using Randn_1d_cpu = Rand_xd<T, edim_1, edist_n, edev_cpu>;		
		
		template <class T>
		using Randn_2d_cpu = Rand_xd<T, edim_2, edist_n, edev_cpu>;		
		
		template <class T>
		using Randn_3d_cpu = Rand_xd<T, edim_3, edist_n, edev_cpu>;		
		
		/* poisson distribution */
		template <class T>
		using Randp_1d_cpu = Rand_xd<T, edim_1, edist_p, edev_cpu>;		
		
		template <class T>
		using Randp_2d_cpu = Rand_xd<T, edim_2, edist_p, edev_cpu>;		
		
		template <class T>
		using Randp_3d_cpu = Rand_xd<T, edim_3, edist_p, edev_cpu>;
	}

	/* 1d distribution */
	namespace mt
	{
		template <class T, eDist Dist>
		class Rand_xd<T, edim_1, Dist, edev_cpu>: public Gen_seed
		{
			public:
				using value_type = T;

				static const eDev device = edev_cpu;

				T m_sc;			// scaling
				T m_sft;		// shift

				Rand_xd(dt_uint64 seed=300183, dt_int32 n_iter=1): Gen_seed(seed), m_sc{1}, m_sft{0}
				{
					set_seed(seed, n_iter);
				}

				void set_in_data(const dt_uint64& seed, const T& sc, const T& sft, dt_int32 n_iter=1)
				{
					set_seed(seed, n_iter);
					set_sc_sft(sc, sft);
				}

				void set_seed(dt_uint64 seed, dt_int32 n_iter=1)
				{
					auto gen_seed = this->seed_1d(seed, n_iter);

					rnd.set_seed(gen_seed);
				}

				void set_sc_sft(const T& sc, const T& sft)
				{
					m_sc = sc;
					m_sft = sft;
				}

				T operator()()
				{
					return rnd(m_sc, m_sft);
				}

				T operator()(const T& sc)
				{
					return rnd(sc);
				}

				T operator()(const T& sc, const T& sft)
				{
					return rnd(sc, sft);
				}

			private:
				Rnd_Dist<T, Dist, edev_cpu> rnd;
		};
	}

	/* 2d distribution */
	namespace mt
	{
		template <class T, eDist Dist>
		class Rand_xd<T, edim_2, Dist, edev_cpu>: public Gen_seed
		{
			public:
				using value_type = T;

				static const eDev device = edev_cpu;

				R_2d<T> m_sc;			// scaling
				R_2d<T> m_sft;			// shift
				R_2d<dt_bool> m_b_dim;	// active dimensions

				Rand_xd(dt_uint64 seed=300183, dt_int32 n_iter=1): Gen_seed(seed), 
					m_sc{1, 1}, m_b_dim{true, true}
				{
					set_seed(seed, n_iter);
				}

				void set_in_data(const dt_uint64& seed, const R_2d<T>& sc, const R_2d<T>& sft, dt_int32 n_iter=1)
				{
					set_seed(seed, n_iter);
					set_sc_sft(sc, sft);
				}

				void set_seed(dt_uint64 seed, dt_int32 n_iter=1)
				{
					auto gen_seed = this->seed_2d(seed, n_iter);

					rnd_x.set_seed(gen_seed.x);
					rnd_y.set_seed(gen_seed.y);
				}

				void set_sc_sft(const R_2d<T>& sc, const R_2d<T>& sft)
				{
					m_sc = sc;
					m_sft = sft;
				}

				void set_act_dim(const R_2d<dt_bool>& b_dim)
				{
					m_b_dim = b_dim;
				}

				R_2d<T> operator()()
				{
					return this->operator()(m_sc, m_sft);
				}

				R_2d<T> operator()(const R_2d<T>& sc)
				{
					return gen_rnd(sc);
				}

				R_2d<T> operator()(const R_2d<T>& sc, const R_2d<T>& sft)
				{
					return gen_rnd(sc, sft);
				}

			private:
				R_2d<T> gen_rnd(const R_2d<T>& sc) 
				{
					return {(m_b_dim.x)?(rnd_x(sc.x)):T(0), (m_b_dim.y)?(rnd_y(sc.y)):T(0)};
				}

				R_2d<T> gen_rnd(const R_2d<T>& sc, const R_2d<T>& sft) 
				{
					return {(m_b_dim.x)?(rnd_x(sc.x, sft.x)):T(0), (m_b_dim.y)?(rnd_y(sc.y, sft.y)):T(0)};
				}				

				Rnd_Dist<T, Dist, edev_cpu> rnd_x;
				Rnd_Dist<T, Dist, edev_cpu> rnd_y;
		};
	}

	/* 3d distribution */
	namespace mt
	{
		template <class T, eDist Dist>
		class Rand_xd<T, edim_3, Dist, edev_cpu>: public Gen_seed
		{
			public:
				using value_type = T;

				static const eDev device = edev_cpu;

				R_3d<T> m_sc;			// scaling
				R_3d<T> m_sft;			// shift
				R_3d<dt_bool> m_b_dim;	// active dimensions

				Rand_xd(dt_uint64 seed=300183, dt_int32 n_iter=1): Gen_seed(seed), 
					m_sc{1, 1, 1}, m_sft{0, 0, 0}, m_b_dim{true, true, true}
				{
					set_seed(seed, n_iter);
				}

				void set_in_data(const dt_uint64& seed, const R_3d<T>& sc, const R_3d<T>& sft, dt_int32 n_iter=1)
				{
					set_seed(seed, n_iter);
					set_sc_sft(sc, sft);
				}

				void set_seed(dt_uint64 seed, dt_int32 n_iter=1)
				{
					auto gen_seed = this->seed_3d(seed, n_iter);

					rnd_x.set_seed(gen_seed.x);
					rnd_y.set_seed(gen_seed.y);
					rnd_z.set_seed(gen_seed.z);
				}

				void set_sc_sft(const R_3d<T>& sc, const R_3d<T>& sft)
				{
					m_sc = sc;
					m_sft = sft;
				}

				void set_act_dim(const R_3d<dt_bool>& b_dim)
				{
					m_b_dim = b_dim;
				}

				R_3d<T> operator()()
				{
					return this->operator()(m_sc, m_sft);
				}

				R_3d<T> operator()(const R_3d<T>& sc)
				{
					return gen_rnd(sc);
				}

				R_3d<T> operator()(const R_3d<T>& sc, const R_3d<T>& sft)
				{
					return gen_rnd(sc, sft);
				}

			private:
				R_3d<T> gen_rnd(const R_3d<T>& sc) 
				{
					return {(m_b_dim.x)?(rnd_x(sc.x)):T(0), (m_b_dim.y)?(rnd_y(sc.y)):T(0), (m_b_dim.z)?(rnd_z(sc.z)):T(0)};
				}

				R_3d<T> gen_rnd(const R_3d<T>& sc, const R_3d<T>& sft) 
				{
					return {(m_b_dim.x)?(rnd_x(sc.x, sft.x)):T(0), (m_b_dim.y)?(rnd_y(sc.y, sft.y)):T(0), (m_b_dim.z)?(rnd_z(sc.z, sft.z)):T(0)};
				}				

				Rnd_Dist<T, Dist, edev_cpu> rnd_x;
				Rnd_Dist<T, Dist, edev_cpu> rnd_y;
				Rnd_Dist<T, Dist, edev_cpu> rnd_z;
		};
	}

	/* fcns */
	namespace mt
	{
		// add uniform noise
		template <class TVctr_i, class TVctr_o>
		enable_if_vctr_cpu_and_vctr_cpu<TVctr_i, TVctr_o, void>
		fcn_add_unif_nois(TVctr_i& mx_i, Value_type<TVctr_o> sc, dt_uint64 seed_i, TVctr_o& mx_o, Stream_cpu* pstream = nullptr)
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
				Randu_1d_cpu<T> rnd(vctr_seed[range.istm]);

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
		fcn_add_gauss_nois(TVctr_i& mx_i, Value_type<TVctr_o> sigma, dt_uint64 seed_i, TVctr_o& mx_o, Stream_cpu* pstream = nullptr)
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
				Randn_1d_cpu<T> rnd(vctr_seed[range.istm]);

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
		fcn_add_poiss_nois(TVctr_i& mx_i, double sc_i, dt_uint64 seed_i, TVctr_o& mx_o, Stream_cpu* pstream = nullptr)
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
				Randp_1d_cpu<T> rnd(vctr_seed[range.istm]);

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

#endif