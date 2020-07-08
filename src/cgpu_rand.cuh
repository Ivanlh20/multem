/*
 * This file is part of MULTEM.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * MULTEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MULTEM is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef CGPU_RAND_H
#define CGPU_RAND_H

#include <random>

#include "math.cuh"
#include "types.cuh"
#include "lin_alg_def.cuh"

#ifdef __CUDACC__
	#include <cuda.h>
	#include <cuda_runtime.h>
#endif

namespace mt
{
	/******************** Uniform 2d distribution *******************/
	template <class T, eDevice dev>
	struct Rand_2d;

	template <class T>
	struct Rand_2d<T, e_host>
	{
		public:
			using value_type = T;

			static const eDevice device = e_host;

			Rand_2d(): m_seed(300183), m_l_x(1), m_l_y(1), 
			m_bx(true), m_by(true), rand_x(T(0), T(1)), 
			rand_y(T(0), T(1)), rand_u(1){}

			void seed(int seed_i, int n_iter=1)
			{
				m_seed = seed_i;

				gen_u.seed(m_seed);
				rand_u.reset();

				unsigned int seed_x, seed_y;
				for(auto ic = 0; ic<n_iter; ic++)
				{
					seed_x = rand_u(gen_u);
					seed_y = rand_u(gen_u);
				}

				gen_x.seed(seed_x);
				rand_x.reset();

				gen_y.seed(seed_y);
				rand_y.reset();
			}

			void set_box_size(const T &l_x, const T &l_y)
			{
				m_l_x = l_x;
				m_l_y = l_y;			
			}

			void set_activation(const bool &bx, const bool &by)
			{
				m_bx = bx;
				m_by = by;
			}

			r2d<T> operator()()
			{
				return r2d<T>(gen_rand_x(m_l_x), gen_rand_y(m_l_y));
			}

			r2d<T> operator()(const T &x, const T &y)
			{
				return r2d<T>(gen_rand_x(x), gen_rand_y(y));
			}

		private:
			T gen_rand_x(const T &x) 
			{
				return ((m_bx)?(x*rand_x(gen_x)):0); 
			}

			T gen_rand_y(const T &y) 
			{
				return ((m_by)?(y*rand_y(gen_y)):0); 
			}

			int m_seed;
			T m_l_x;
			T m_l_y;
			bool m_bx;
			bool m_by;

			std::mt19937_64 gen_u;
			std::mt19937_64 gen_x;
			std::mt19937_64 gen_y;
			std::uniform_int_distribution<int> rand_u;
			std::uniform_real_distribution<T> rand_x;
			std::uniform_real_distribution<T> rand_y;
	};

	/******************** Uniform 3d distribution *******************/
	template <class T, eDevice dev>
	struct Rand_3d;

	template <class T>
	struct Rand_3d<T, e_host>
	{
		public:
			using value_type = T;

			static const eDevice device = e_host;

			Rand_3d(): m_seed(300183), m_l_x(1), m_l_y(1), m_l_z(1), 
			m_bx(true), m_by(true), m_bz(true), rand_x(T(0), T(1)), 
			rand_y(T(0), T(1)), rand_z(T(0), T(1)), rand_u(1){}

			void seed(int seed, int n_iter=1)
			{
				n_iter = max(1, n_iter);

				m_seed = seed;

				gen_u.seed(m_seed);
				rand_u.reset();

				unsigned int seed_x, seed_y, seed_z;
				for(auto ic = 0; ic<n_iter; ic++)
				{
					seed_x = rand_u(gen_u);
					seed_y = rand_u(gen_u);
					seed_z = rand_u(gen_u);
				}

				gen_x.seed(seed_x);
				rand_x.reset();

				gen_y.seed(seed_y);
				rand_y.reset();

				gen_z.seed(seed_z);
				rand_z.reset();
			}

			void set_box_size(const T &l_x, const T &l_y, const T &l_z)
			{
				m_l_x = l_x;
				m_l_y = l_y;
				m_l_z = l_z;
			}

			void set_activation(const bool &bx, const bool &by, const bool &bz)
			{
				m_bx = bx;
				m_by = by;
				m_bz = bz;
			}

			r3d<T> operator()()
			{
				return r3d<T>(gen_rand_x(m_l_x), gen_rand_y(m_l_y), gen_rand_z(m_l_z));
			}

			r3d<T> operator()(const T &x, const T &y, const T &z)
			{
				return r3d<T>(gen_rand_x(x), gen_rand_y(y), gen_rand_z(z));
			}

		private:
			T gen_rand_x(const T &x) 
			{
				return ((m_bx)?(x*rand_x(gen_x)):0); 
			}

			T gen_rand_y(const T &y) 
			{
				return ((m_by)?(y*rand_y(gen_y)):0); 
			}

			T gen_rand_z(const T &z) 
			{
				return ((m_bz)?(z*rand_z(gen_z)):0); 
			}

			int m_seed;
			T m_l_x;
			T m_l_y;
			T m_l_z;
			bool m_bx;
			bool m_by;
			bool m_bz;
			std::mt19937_64 gen_u;
			std::mt19937_64 gen_x;
			std::mt19937_64 gen_y;
			std::mt19937_64 gen_z;
			std::uniform_int_distribution<int> rand_u;
			std::uniform_real_distribution<T> rand_x;
			std::uniform_real_distribution<T> rand_y;
			std::uniform_real_distribution<T> rand_z;
	};

	/******************** normal 2d distribution *******************/
	template <class T, eDevice dev>
	struct Randn_2d;

	template <class T>
	struct Randn_2d<T, e_host>
	{
		public:
			using value_type = T;

			static const eDevice device = e_host;

			Randn_2d(): m_seed(300183), m_l_x(1), m_l_y(1), 
			m_bx(true), m_by(true), randn_x(T(0), T(1)), 
			randn_y(T(0), T(1)), rand_u(1){}

			void seed(int seed_i, int n_iter=1)
			{
				m_seed = seed_i;

				gen_u.seed(m_seed);
				rand_u.reset();

				unsigned int seed_x, seed_y;
				for(auto ic = 0; ic<n_iter; ic++)
				{
					seed_x = rand_u(gen_u);
					seed_y = rand_u(gen_u);
				}

				gen_x.seed(seed_x);
				randn_x.reset();

				gen_y.seed(seed_y);
				randn_y.reset();
			}

			void set_box_size(const T &l_x, const T &l_y)
			{
				m_l_x = l_x;
				m_l_y = l_y;			
			}

			void set_activation(const bool &bx, const bool &by)
			{
				m_bx = bx;
				m_by = by;
			}

			r2d<T> operator()()
			{
				return r2d<T>(gen_rand_x(m_l_x), gen_rand_y(m_l_y));
			}

			r2d<T> operator()(const T &x, const T &y)
			{
				return r2d<T>(gen_rand_x(x), gen_rand_y(y));
			}

		private:
			T gen_rand_x(const T &x) 
			{
				return ((m_bx)?(x*randn_x(gen_x)):0); 
			}

			T gen_rand_y(const T &y) 
			{
				return ((m_by)?(y*randn_y(gen_y)):0); 
			}

			int m_seed;
			T m_l_x;
			T m_l_y;
			bool m_bx;
			bool m_by;

			std::mt19937_64 gen_u;
			std::mt19937_64 gen_x;
			std::mt19937_64 gen_y;
			std::uniform_int_distribution<int> rand_u;
			std::normal_distribution<T> randn_x;
			std::normal_distribution<T> randn_y;
	};

	/******************** normal 3d distribution *******************/
	template <class T, eDevice dev>
	struct Randn_3d;

	template <class T>
	struct Randn_3d<T, e_host>
	{
		public:
			using value_type = T;

			static const eDevice device = e_host;

			Randn_3d(): m_seed(300183), m_l_x(1), m_l_y(1), m_l_z(1), 
			m_bx(true), m_by(true), m_bz(true), randn_x(T(0), T(1)), 
			randn_y(T(0), T(1)), randn_z(T(0), T(1)), rand_u(1){}

			void seed(int seed, int n_iter=1)
			{
				n_iter = max(1, n_iter);

				m_seed = seed;

				gen_u.seed(m_seed);
				rand_u.reset();

				unsigned int seed_x, seed_y, seed_z;
				for(auto ic = 0; ic<n_iter; ic++)
				{
					seed_x = rand_u(gen_u);
					seed_y = rand_u(gen_u);
					seed_z = rand_u(gen_u);
				}

				gen_x.seed(seed_x);
				randn_x.reset();

				gen_y.seed(seed_y);
				randn_y.reset();

				gen_z.seed(seed_z);
				randn_z.reset();
			}

			void set_box_size(const T &l_x, const T &l_y, const T &l_z)
			{
				m_l_x = l_x;
				m_l_y = l_y;
				m_l_z = l_z;
			}

			void set_activation(const bool &bx, const bool &by, const bool &bz)
			{
				m_bx = bx;
				m_by = by;
				m_bz = bz;
			}

			r3d<T> operator()()
			{
				return r3d<T>(gen_rand_x(m_l_x), gen_rand_y(m_l_y), gen_rand_z(m_l_z));
			}

			r3d<T> operator()(const T &x, const T &y, const T &z)
			{
				return r3d<T>(gen_rand_x(x), gen_rand_y(y), gen_rand_z(z));
			}

		private:
			T gen_rand_x(const T &x) 
			{
				return ((m_bx)?(x*randn_x(gen_x)):0); 
			}

			T gen_rand_y(const T &y) 
			{
				return ((m_by)?(y*randn_y(gen_y)):0); 
			}

			T gen_rand_z(const T &z) 
			{
				return ((m_bz)?(z*randn_z(gen_z)):0); 
			}

			int m_seed;
			T m_l_x;
			T m_l_y;
			T m_l_z;
			bool m_bx;
			bool m_by;
			bool m_bz;
			std::mt19937_64 gen_u;
			std::mt19937_64 gen_x;
			std::mt19937_64 gen_y;
			std::mt19937_64 gen_z;
			std::uniform_int_distribution<int> rand_u;
			std::normal_distribution<T> randn_x;
			std::normal_distribution<T> randn_y;
			std::normal_distribution<T> randn_z;
	};

	/***************************************************************/

	// add Gaussian noise
	template <class TVector>
	void add_gauss_nois(Stream<e_host> &stream, TVector &M_i, 
	Value_type<TVector> sigma, TVector &M_o)
	{	
		using T = Value_type<TVector>;

		auto thr_gauss_nois = [&](const Range_2d &range)
		{
			std::random_device rd;
			std::mt19937_64 gen(rd());
			std::normal_distribution<T> randn;

			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				 M_o[ixy] = M_i[ixy] + sigma*randn(gen);
			}
		};

		stream.set_n_act_stream(M_i.size());
		stream.set_grid(M_i.size(), 1);
		stream.exec(thr_gauss_nois);
	};

	// add Gaussian noise
	template <class TVector>
	TVector add_gauss_nois(Stream<e_host> &stream, TVector &M_i, 
	Value_type<TVector> sigma)
	{	
		using T = Value_type<TVector>;

		TVector M(M_i.size());
		auto thr_gauss_nois = [&](const Range_2d &range)
		{
			std::random_device rd;
			std::mt19937_64 gen(rd());
			std::normal_distribution<T> randn;

			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				 M[ixy] = M_i[ixy] + sigma*randn(gen);
			}
		};

		stream.set_n_act_stream(M_i.size());
		stream.set_grid(M_i.size(), 1);
		stream.exec(thr_gauss_nois);

		return M;
	};
	/***************************************************************/
	// add Poisson noise
	template <class TVector>
	TVector add_poiss_nois(Stream<e_host> &stream, TVector &M_i, 
	Value_type<TVector> scf)
	{	
		using T = Value_type<TVector>;

		TVector M(M_i.size());
		auto thr_poiss_nois = [&](const Range_2d &range)
		{
			std::random_device rd;
			std::mt19937_64 gen(rd());
			std::poisson_distribution<int> randp;

			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				auto x0 = scf*M_i[ixy];
				randp.param(std::poisson_distribution<int>::param_type(x0));
				M[ixy] = randp(gen);
			}
		};

		stream.set_n_act_stream(M_i.size());
		stream.set_grid(M_i.size(), 1);
		stream.exec(thr_poiss_nois);

		return M;
	};

	// add Poisson noise
	template <class TVector>
	TVector add_poiss_nois_by_SNR(Stream<e_host> &stream, TVector &Im_i, 
		Value_type<TVector> SNR_i, Value_type<TVector> &scl_o)
	{	
		using T = Value_type<TVector>;

		auto get_SNR = [](Stream<e_host> &stream, TVector &Im, T Im_std, T scf)->T
		{
			T x_mean = 0;
			T x_var = 0;
			auto thr_SNR = [&](const Range_2d &range)
			{
				std::mt19937_64 gen;
				std::poisson_distribution<int> rand;

				T x_mean_partial = 0;
				T x_var_partial = 0;
				for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
				{
					auto x0 = scf*Im[ixy];
					rand.param(std::poisson_distribution<int>::param_type(x0));
					auto xn = rand(gen)-x0;
					x_mean_partial += xn;
					x_var_partial += xn*xn;
				}

				stream.stream_mutex.lock();
				x_mean += x_mean_partial;
				x_var += x_var_partial;
				stream.stream_mutex.unlock();
			};

			stream.set_n_act_stream(Im.size());
			stream.set_grid(1, Im.size());
			stream.exec(thr_SNR);

			x_mean /= Im.size();
			x_var = x_var/Im.size()-x_mean*x_mean;

			return scf*Im_std/sqrt(x_var);
		};

		auto Im_i_std = sqrt(variance(stream, Im_i));

		T SNR_k = get_SNR(stream, Im_i, Im_i_std, 1);
	
		T k_0 = 1;
		T k_e = 1;

		if(SNR_k<SNR_i)
		{
			do
			{
				k_e *= 2;
				SNR_k = get_SNR(stream, Im_i, Im_i_std, k_e);
			}
			while (SNR_k < SNR_i);
			k_0 = k_e/2;	
		}
		else
		{
			do
			{
				k_0 /= 2;
				SNR_k = get_SNR(stream, Im_i, Im_i_std, k_0);
			}
			while (SNR_k >= SNR_i);
			k_e = 2*k_0;
		}


		// bisection method
		int ic = 0;
		do
		{
			scl_o = 0.5*(k_0 + k_e);
			auto SNR_k = get_SNR(stream, Im_i, Im_i_std, scl_o);

			if(SNR_k < SNR_i)
			{
				k_0 = scl_o;
			}
			else
			{
				k_e = scl_o;
			}
			ic++;
		} while ((fabs(SNR_i-SNR_k)>0.1) && (ic<10));

		// add Poisson noise
		return add_poiss_nois(stream, Im_i, scl_o);
	}

} // namespace mt

#endif