/*
 * This file is part of MULTEM.
 * Copyright 2017 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef LIN_ALG_DEF_H
#define LIN_ALG_DEF_H

#include <random>

#include "math.cuh"

namespace mt
{
	template <class T>
	struct r3d;

	// template <class T, eDevice dev>
	// class Vector;

	template <class T>
	struct Rand_r3d
	{
		public:
			void set_seed(int seed_new)
			{	
				gen_u.seed(seed_new);
				rand_u.reset();
				rand_temp.reset();
			}

			r3d<T> operator()()
			{
				return r3d<T>(rand_u(gen_u), rand_u(gen_u), rand_u(gen_u));
			}

			T temp()
			{
				rand_temp(gen_u);
			}

		private:
			std::mt19937_64 gen_u;
			std::uniform_real_distribution<T> rand_u;
			std::uniform_real_distribution<T> rand_temp;
	};

	// /*******************************mxn matrix*********************************/
	// template <class T>
	// struct Mmxn
	// {		
	// 	using value_type = T;

	// 	int m;
	// 	int n;

	// 	Vector<T, e_host> data;


	// 	DEVICE_CALLABLE
	// 	inline Mmxn(const int &m_i = 1, const int &n_i = 1)
	// 	{
	// 		m = m_i;
	// 		n = n_i;

	// 		data.resize(m*n);
	// 	} 

	// 	DEVICE_CALLABLE
	// 	inline Mmxn<T>& operator+=(const Mmxn<T> r)
	// 	{
	// 		
	// 		for(auto i=0; i<data.size(); i++)
	// 		{
	// 			data[i] += r.data[i];		
	// 		}

	// 		return *this;
	// 	}

	// 	DEVICE_CALLABLE
	// 	inline Mmxn<T>& operator-=(const Mmxn<T> r)
	// 	{
	// 		for(auto i=0; i<data.size(); i++)
	// 		{
	// 			data[i] -= r.data[i];		
	// 		}

	// 		return *this;
	// 	}

	// 	DEVICE_CALLABLE
	// 	inline Mmxn<T>& operator*=(const T r)
	// 	{
	// 		if(this->n!=r.m)
	// 		{
	// 			return;
	// 		}

	// 		for(auto ir=0; ir<m; ir++)
	// 		{
	// 			for(auto ic=0; ic<r.n; ic++)
	// 			{
	// 				irc = r.n*ir + ic;
	// 				for(auto i=0; i<n; i++)
	// 				{
	// 					data[irc] += r.data[i];		
	// 				}	
	// 			}
	// 		}

	// 		return *this;
	// 	}

	// 	DEVICE_CALLABLE
	// 	inline Mmxn<T>& operator/=(const T r)
	// 	{
	// 		for(auto i=0; i<data.size(); i++)
	// 		{
	// 			data[i] /= r.data[i];		
	// 		}

	// 		return *this;
	// 	}

	// };

	// template <class X>
	// DEVICE_CALLABLE
	// inline Mmxn<X> operator+(const Mmxn<X> &lhs, const Mmxn<X> &rhs)
	// {
	// 	return Mmxn<X>(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
	// }

	// template <class X>
	// DEVICE_CALLABLE
	// inline Mmxn<X> operator+(const Mmxn<X> &lhs, const X &rhs)
	// {
	// 	return Mmxn<X>(lhs.x+rhs, lhs.y+rhs, lhs.z+rhs);
	// }

	// template <class X>
	// DEVICE_CALLABLE
	// inline Mmxn<X> operator+(const X &lhs, const Mmxn<X> &rhs)
	// {
	// 	return Mmxn<X>(lhs+rhs.x, lhs+rhs.y, lhs+rhs.z);
	// }


	// template <class X>
	// DEVICE_CALLABLE
	// inline Mmxn<X> operator-(const Mmxn<X> &lhs, const Mmxn<X> &rhs)
	// {
	// 	return Mmxn<X>(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
	// }

	// template <class X>
	// DEVICE_CALLABLE
	// inline Mmxn<X> operator-(const Mmxn<X> &lhs, const X &rhs)
	// {
	// 	return Mmxn<X>(lhs.x-rhs, lhs.y-rhs, lhs.z-rhs);
	// }

	// template <class X>
	// DEVICE_CALLABLE
	// inline Mmxn<X> operator-(const X &lhs, const Mmxn<X> &rhs)
	// {
	// 	return Mmxn<X>(lhs-rhs.x, lhs-rhs.y, lhs-rhs.z);
	// }

	// template <class X>
	// DEVICE_CALLABLE
	// inline Mmxn<X> operator*(const Mmxn<X> &lhs, const Mmxn<X> &rhs)
	// {
	// 	return Mmxn<X>(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z);
	// }

	// template <class X>
	// DEVICE_CALLABLE
	// inline Mmxn<X> operator*(const Mmxn<X> &lhs, const X &rhs)
	// {
	// 	return Mmxn<X>(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs);
	// }

	// template <class X>
	// DEVICE_CALLABLE
	// inline Mmxn<X> operator*(const X &lhs, const Mmxn<X> &rhs)
	// {
	// 	return Mmxn<X>(lhs*rhs.x, lhs*rhs.y, lhs*rhs.z);
	// }

	// template <class X>
	// DEVICE_CALLABLE
	// inline Mmxn<X> operator/(const Mmxn<X> &lhs, const X &rhs)
	// {
	// 	return Mmxn<X>(lhs.x/rhs, lhs.y/rhs, lhs.z/rhs);
	// }

	// template <class X>
	// DEVICE_CALLABLE
	// inline Mmxn<X> fmin(const Mmxn<X> lhs, const Mmxn<X> rhs)
	// {
	// 	return Mmxn<X>(::fmin(lhs.x, rhs.x), ::fmin(lhs.y, rhs.y), ::fmin(lhs.z, rhs.z));
	// }

	// template <class X>
	// DEVICE_CALLABLE
	// inline Mmxn<X> fmax(const Mmxn<X> lhs, const Mmxn<X> rhs)
	// {
	// 	return Mmxn<X>(::fmax(lhs.x, rhs.x), ::fmax(lhs.y, rhs.y), ::fmax(lhs.z, rhs.z));
	// }

	// template <typename X>
	// DEVICE_CALLABLE
	// inline X norm(const Mmxn<X>& r)
	// {
	// 	return r.x*r.x + r.y*r.y + r.z*r.z;
	// }

	// template <typename X>
	// DEVICE_CALLABLE
	// inline X module(const Mmxn<X>& r)
	// {
	// 	return sqrt(norm(r));
	// }

	// template <typename X>
	// DEVICE_CALLABLE
	// inline Mmxn<X> normalized(const Mmxn<X>& r)
	// {
	// 	return r/r.module();
	// }


	// // template <class TVector>
	// // DEVICE_CALLABLE
	// // inline Mmxn<T> matrix_prod(const TVector &Rm)
	// // {
	// // 	Mmxn<T> r_o;
	// // 	r_o.x = Rm[0]*x + Rm[3]*y + Rm[6]*z;
	// // 	r_o.y = Rm[1]*x + Rm[4]*y + Rm[7]*z;
	// // 	r_o.z = Rm[2]*x + Rm[5]*y + Rm[8]*z;
	// // 	return r_o;
	// // }

	// template <class X>
	// DEVICE_CALLABLE
	// inline X dot(const Mmxn<X> &lhs, const Mmxn<X> &rhs)
	// {
	// 	return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
	// }

	// template <class X>
	// DEVICE_CALLABLE
	// inline X angle(const Mmxn<X> &lhs, const Mmxn<X> &rhs)
	// {
	// 	return acos(dot(lhs, rhs)/(lhs.module()*rhs.module()));
	// }

	/******************************2d vector*********************************/
	template <class T>
	struct r2d
	{		
		using value_type = T;

		T x;
		T y;

		DEVICE_CALLABLE
		inline r2d(const T &x_i = T(), const T &y_i = T())
		{
			x = x_i;
			y = y_i;
		} 

		template <typename X> 
		DEVICE_CALLABLE
		inline r2d(const r2d<X> &r)
		{
			x = r.x;
			y = r.y;
		} 

		template <typename X> 
		friend DEVICE_CALLABLE
		inline r2d<X> operator-(const r2d<X> r);

		DEVICE_CALLABLE
		inline r2d<T>& operator+=(const r2d<T> r)
		{
			x += r.x;
			y += r.y;
			return *this;
		}

		DEVICE_CALLABLE
		inline r2d<T>& operator-=(const r2d<T> r)
		{
			x -= r.x;
			y -= r.y;
			return *this;
		}

		DEVICE_CALLABLE
		inline r2d<T>& operator*=(const T r)
		{
			x *= r;
			y *= r;
			return *this;
		}

		DEVICE_CALLABLE
		inline r2d<T>& operator/=(const T r)
		{
			x /= r;
			y /= r;
			return *this;
		}

		DEVICE_CALLABLE
		inline T norm() const
		{
			return x*x + y*y;
		}

		DEVICE_CALLABLE
		inline T module() const
		{
			return sqrt(norm());
		}

		DEVICE_CALLABLE
		inline void normalized()
		{
			*this /= module();
		}

		template <class TVector>
		DEVICE_CALLABLE
		inline r2d<T> apply_matrix(const TVector &Rm)
		{
			r2d<T> r_o;
			r_o.x = Rm[0]*x + Rm[2]*y;
			r_o.y = Rm[1]*x + Rm[3]*y;
			return r_o;
		}

		template <class TVector>
		DEVICE_CALLABLE
		inline r2d<T> rotate(const TVector &Rm, const r2d<T> &p0)
		{
			r2d<T> r(x, y); 
			r -= p0;
			return (apply_matrix(Rm) + p0);
		}

	};

	template <class X>
	DEVICE_CALLABLE
	inline r2d<X> operator-(const r2d<X> r)
	{
		return r2d<X>(-r.x, -r.y);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r2d<X> operator+(const r2d<X> &lhs, const r2d<X> &rhs)
	{
		return r2d<X>(lhs.x+rhs.x, lhs.y+rhs.y);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r2d<X> operator+(const r2d<X> &lhs, const X &rhs)
	{
		return r2d<X>(lhs.x+rhs, lhs.y+rhs);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r2d<X> operator+(const X &lhs, const r2d<X> &rhs)
	{
		return r2d<X>(lhs+rhs.x, lhs+rhs.y);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r2d<X> operator-(const r2d<X> &lhs, const r2d<X> &rhs)
	{
		return r2d<X>(lhs.x-rhs.x, lhs.y-rhs.y);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r2d<X> operator-(const r2d<X> &lhs, const X &rhs)
	{
		return r2d<X>(lhs.x-rhs, lhs.y-rhs);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r2d<X> operator-(const X &lhs, const r2d<X> &rhs)
	{
		return r2d<X>(lhs-rhs.x, lhs-rhs.y);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r2d<X> operator*(const r2d<X> &lhs, const r2d<X> &rhs)
	{
		return r2d<X>(lhs.x*rhs.x, lhs.y*rhs.y);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r2d<X> operator*(const r2d<X> &lhs, const X &rhs)
	{
		return r2d<X>(lhs.x*rhs, lhs.y*rhs);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r2d<X> operator*(const X &lhs, const r2d<X> &rhs)
	{
		return r2d<X>(lhs*rhs.x, lhs*rhs.y);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r2d<X> operator/(const r2d<X> &lhs, const X &rhs)
	{
		return r2d<X>(lhs.x/rhs, lhs.y/rhs);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r2d<X> fmin(const r2d<X> lhs, const r2d<X> rhs)
	{
		return r2d<X>(::fmin(lhs.x, rhs.x), ::fmin(lhs.y, rhs.y));
	}

	template <class X>
	DEVICE_CALLABLE
	inline r2d<X> fmax(const r2d<X> lhs, const r2d<X> rhs)
	{
		return r2d<X>(::fmax(lhs.x, rhs.x), ::fmax(lhs.y, rhs.y));
	}

	template <typename X>
	DEVICE_CALLABLE
	inline X norm(const r2d<X>& r)
	{
		return r.x*r.x + r.y*r.y;
	}

	template <typename X>
	DEVICE_CALLABLE
	inline X module(const r2d<X>& r)
	{
		return sqrt(norm(r));
	}

	template <typename X>
	DEVICE_CALLABLE
	inline r2d<X> normalized(const r2d<X>& r)
	{
		return r/module(r);
	}

	// template <class TVector>
	// DEVICE_CALLABLE
	// inline r2d<T> matrix_prod(const TVector &Rm)
	// {
	// 	r2d<T> r_o;
	// 	r_o.x = Rm[0]*x + Rm[3]*y + Rm[6]*z;
	// 	r_o.y = Rm[1]*x + Rm[4]*y + Rm[7]*z;
	// 	r_o.z = Rm[2]*x + Rm[5]*y + Rm[8]*z;
	// 	return r_o;
	// }

	template <class X>
	DEVICE_CALLABLE
	inline X dot(const r2d<X> &lhs, const r2d<X> &rhs)
	{
		return lhs.x*rhs.x + lhs.y*rhs.y;
	}

	template <class X>
	DEVICE_CALLABLE
	inline X angle(const r2d<X> &lhs, const r2d<X> &rhs)
	{
		return acos(dot(lhs, rhs)/(lhs.module()*rhs.module()));
	}

	/******************************3d vector*********************************/
	template <class T>
	struct r3d
	{		
		using value_type = T;

		T x;
		T y;
		T z;

		DEVICE_CALLABLE
		r3d(const T &x_i = T(), const T &y_i = T(), const T &z_i = T())
		{
			x = x_i;
			y = y_i;
			z = z_i;
		} 

		template <typename X> 
		DEVICE_CALLABLE
		r3d(const r3d<X> &r)
		{
			x = r.x;
			y = r.y;
			z = r.z;
		} 

		DEVICE_CALLABLE
		inline r3d<T>& operator+=(const r3d<T> r)
		{
			x += r.x;
			y += r.y;
			z += r.z;
			return *this;
		}

		DEVICE_CALLABLE
		inline r3d<T>& operator-=(const r3d<T> r)
		{
			x -= r.x;
			y -= r.y;
			z -= r.z;
			return *this;
		}

		DEVICE_CALLABLE
		inline r3d<T>& operator*=(const T r)
		{
			x *= r;
			y *= r;
			z *= r;
			return *this;
		}

		DEVICE_CALLABLE
		inline r3d<T>& operator/=(const T r)
		{
			x /= r;
			y /= r;
			z /= r;
			return *this;
		}

		DEVICE_CALLABLE
		inline T norm()
		{
			return pow(x, 2)+pow(y, 2)+pow(z, 2);
		}

		DEVICE_CALLABLE
		inline T module()
		{
			return sqrt(norm());
		}

		DEVICE_CALLABLE
		inline void normalized()
		{
			*this /= module();
		}

		template <class TVector>
		DEVICE_CALLABLE
		inline r3d<T> apply_matrix(const TVector &Rm)
		{
			r3d<T> r_o;
			r_o.x = Rm[0]*x + Rm[3]*y + Rm[6]*z;
			r_o.y = Rm[1]*x + Rm[4]*y + Rm[7]*z;
			r_o.z = Rm[2]*x + Rm[5]*y + Rm[8]*z;
			return r_o;
		}

		template <class TVector>
		DEVICE_CALLABLE
		inline r3d<T> rotate(const TVector &Rm, const r3d<T> &p0) const
		{
			r3d<T> r(x, y, z); 
			r -= p0;
			return (r.apply_matrix(Rm) + p0);
		}

	};

	template <class X>
	DEVICE_CALLABLE
	inline r3d<X> operator+(const r3d<X> &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r3d<X> operator+(const r3d<X> &lhs, const X &rhs)
	{
		return r3d<X>(lhs.x+rhs, lhs.y+rhs, lhs.z+rhs);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r3d<X> operator+(const X &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(lhs+rhs.x, lhs+rhs.y, lhs+rhs.z);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r3d<X> operator-(const r3d<X> &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r3d<X> operator-(const r3d<X> &lhs, const X &rhs)
	{
		return r3d<X>(lhs.x-rhs, lhs.y-rhs, lhs.z-rhs);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r3d<X> operator-(const X &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(lhs-rhs.x, lhs-rhs.y, lhs-rhs.z);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r3d<X> operator*(const r3d<X> &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r3d<X> operator*(const r3d<X> &lhs, const X &rhs)
	{
		return r3d<X>(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r3d<X> operator*(const X &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(lhs*rhs.x, lhs*rhs.y, lhs*rhs.z);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r3d<X> operator/(const r3d<X> &lhs, const X &rhs)
	{
		return r3d<X>(lhs.x/rhs, lhs.y/rhs, lhs.z/rhs);
	}

	template <class X>
	DEVICE_CALLABLE
	inline r3d<X> fmin(const r3d<X> lhs, const r3d<X> rhs)
	{
		return r3d<X>(::fmin(lhs.x, rhs.x), ::fmin(lhs.y, rhs.y), ::fmin(lhs.z, rhs.z));
	}

	template <class X>
	DEVICE_CALLABLE
	inline r3d<X> fmax(const r3d<X> lhs, const r3d<X> rhs)
	{
		return r3d<X>(::fmax(lhs.x, rhs.x), ::fmax(lhs.y, rhs.y), ::fmax(lhs.z, rhs.z));
	}

	template <typename X>
	DEVICE_CALLABLE
	inline X norm(const r3d<X>& r)
	{
		return r.x*r.x + r.y*r.y + r.z*r.z;
	}

	template <typename X>
	DEVICE_CALLABLE
	inline X module(const r3d<X>& r)
	{
		return sqrt(norm(r));
	}

	template <typename X>
	DEVICE_CALLABLE
	inline r3d<X> normalized(const r3d<X>& r)
	{
		return r/r.module();
	}


	// template <class TVector>
	// DEVICE_CALLABLE
	// inline r3d<T> matrix_prod(const TVector &Rm)
	// {
	// 	r3d<T> r_o;
	// 	r_o.x = Rm[0]*x + Rm[3]*y + Rm[6]*z;
	// 	r_o.y = Rm[1]*x + Rm[4]*y + Rm[7]*z;
	// 	r_o.z = Rm[2]*x + Rm[5]*y + Rm[8]*z;
	// 	return r_o;
	// }

	template <class X>
	DEVICE_CALLABLE
	inline X dot(const r3d<X> &lhs, const r3d<X> &rhs)
	{
		return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
	}

	template <class X>
	DEVICE_CALLABLE
	inline X angle(const r3d<X> &lhs, const r3d<X> &rhs)
	{
		return acos(dot(lhs, rhs)/(lhs.module()*rhs.module()));
	}
} // namespace mt
#endif