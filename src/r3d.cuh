/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef R3D_H
#define R3D_H

#include <random>

#include "math.cuh"

namespace multem
{
	template<class T>
	struct r3d;

	// template<class T, eDevice dev>
	// class Vector;

	template<class T>
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

	template<class T>
	struct r3d
	{		
		using value_type = T;

		T x;
		T y;
		T z;

		DEVICE_CALLABLE
		inline r3d(const T &x_i = T(), const T &y_i = T(), const T &z_i = T())
		{
			x = x_i;
			y = y_i;
			z = z_i;
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
		inline void normalized()
		{
			*this /= norm(*this);
		}

		template<class TVector>
		DEVICE_CALLABLE
		inline r3d<T> rotate(const TVector &Rm, const r3d<T> &p0)
		{
			r3d<T> r(x, y, z); 
			r -= p0;
			r3d<T> r_o;
			r_o.x = Rm[0]*r.x + Rm[3]*r.y + Rm[6]*r.z + p0.x;
			r_o.y = Rm[1]*r.x + Rm[4]*r.y + Rm[7]*r.z + p0.y;
			r_o.z = Rm[2]*r.x + Rm[5]*r.y + Rm[8]*r.z + p0.z;
			return r_o;
		}

	};

	template<class X>
	DEVICE_CALLABLE
	inline r3d<X> operator+(const r3d<X> &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
	}

	template<class X>
	DEVICE_CALLABLE
	inline r3d<X> operator+(const r3d<X> &lhs, const X &rhs)
	{
		return r3d<X>(lhs.x+rhs, lhs.y+rhs, lhs.z+rhs);
	}

	template<class X>
	DEVICE_CALLABLE
	inline r3d<X> operator+(const X &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(lhs+rhs.x, lhs+rhs.y, lhs+rhs.z);
	}


	template<class X>
	DEVICE_CALLABLE
	inline r3d<X> operator-(const r3d<X> &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
	}

	template<class X>
	DEVICE_CALLABLE
	inline r3d<X> operator-(const r3d<X> &lhs, const X &rhs)
	{
		return r3d<X>(lhs.x-rhs, lhs.y-rhs, lhs.z-rhs);
	}

	template<class X>
	DEVICE_CALLABLE
	inline r3d<X> operator-(const X &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(lhs-rhs.x, lhs-rhs.y, lhs-rhs.z);
	}

	template<class X>
	DEVICE_CALLABLE
	inline r3d<X> operator*(const r3d<X> &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z);
	}

	template<class X>
	DEVICE_CALLABLE
	inline r3d<X> operator*(const r3d<X> &lhs, const X &rhs)
	{
		return r3d<X>(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs);
	}

	template<class X>
	DEVICE_CALLABLE
	inline r3d<X> operator*(const X &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(lhs*rhs.x, lhs*rhs.y, lhs*rhs.z);
	}

	template<class X>
	DEVICE_CALLABLE
	inline r3d<X> operator/(const r3d<X> &lhs, const X &rhs)
	{
		return r3d<X>(lhs.x/rhs, lhs.y/rhs, lhs.z/rhs);
	}

	template<class X>
	DEVICE_CALLABLE
	inline r3d<X> fmin(const r3d<X> lhs, const r3d<X> rhs)
	{
		return r3d<X>(::fmin(lhs.x, rhs.x), ::fmin(lhs.y, rhs.y), ::fmin(lhs.z, rhs.z));
	}

	template<class X>
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
	inline X abs(const r3d<X>& r)
	{
		return sqrt(norm(r));
	}

	template <typename X>
	DEVICE_CALLABLE
	inline X normalized(const r3d<X>& r)
	{
		return r/norm(r);
	}

	template<class X>
	DEVICE_CALLABLE
	inline X dot(const r3d<X> &lhs, const r3d<X> &rhs)
	{
		return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
	}
} // namespace multem
#endif