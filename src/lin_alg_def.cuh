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

#ifndef LIN_ALG_DEF_H
#define LIN_ALG_DEF_H

#include <random>
#include <vector>

#include "math.cuh"

namespace mt
{
	/************************************************************************/
	/********************************2d vector*******************************/
	/************************************************************************/

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
			r_o.x = Rm[0] * x + Rm[2] * y;
			r_o.y = Rm[1] * x + Rm[3] * y;
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
		return r2d<X>(lhs.x + rhs.x, lhs.y + rhs.y);
	}

	template <class X>
	DEVICE_CALLABLE
		inline r2d<X> operator+(const r2d<X> &lhs, const X &rhs)
	{
		return r2d<X>(lhs.x + rhs, lhs.y + rhs);
	}

	template <class X>
	DEVICE_CALLABLE
		inline r2d<X> operator+(const X &lhs, const r2d<X> &rhs)
	{
		return r2d<X>(lhs + rhs.x, lhs + rhs.y);
	}

	template <class X>
	DEVICE_CALLABLE
		inline r2d<X> operator-(const r2d<X> &lhs, const r2d<X> &rhs)
	{
		return r2d<X>(lhs.x - rhs.x, lhs.y - rhs.y);
	}

	template <class X>
	DEVICE_CALLABLE
		inline r2d<X> operator-(const r2d<X> &lhs, const X &rhs)
	{
		return r2d<X>(lhs.x - rhs, lhs.y - rhs);
	}

	template <class X>
	DEVICE_CALLABLE
		inline r2d<X> operator-(const X &lhs, const r2d<X> &rhs)
	{
		return r2d<X>(lhs - rhs.x, lhs - rhs.y);
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
		return r2d<X>(lhs.x / rhs, lhs.y / rhs);
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
		return r / module(r);
	}

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
		return acos(dot(lhs, rhs) / (lhs.module()*rhs.module()));
	}

	/************************************************************************/
	/********************************3d vector*******************************/
	/************************************************************************/

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
			return pow(x, 2) + pow(y, 2) + pow(z, 2);
		}

		DEVICE_CALLABLE
			inline T module()
		{
			return sqrt(norm());
		}

		DEVICE_CALLABLE
			inline void normalized()
		{
			auto rm = module();
			*this /= module();
		}

		template <class TVector>
		DEVICE_CALLABLE
			inline r3d<T> apply_matrix(const TVector &Rm)
		{
			r3d<T> r_o;
			r_o.x = Rm[0] * x + Rm[3] * y + Rm[6] * z;
			r_o.y = Rm[1] * x + Rm[4] * y + Rm[7] * z;
			r_o.z = Rm[2] * x + Rm[5] * y + Rm[8] * z;
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
		return r3d<X>(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
	}

	template <class X>
	DEVICE_CALLABLE
		inline r3d<X> operator+(const r3d<X> &lhs, const X &rhs)
	{
		return r3d<X>(lhs.x + rhs, lhs.y + rhs, lhs.z + rhs);
	}

	template <class X>
	DEVICE_CALLABLE
		inline r3d<X> operator+(const X &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(lhs + rhs.x, lhs + rhs.y, lhs + rhs.z);
	}

	template <class X>
	DEVICE_CALLABLE
		inline r3d<X> operator-(const r3d<X> &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
	}

	template <class X>
	DEVICE_CALLABLE
		inline r3d<X> operator-(const r3d<X> &lhs, const X &rhs)
	{
		return r3d<X>(lhs.x - rhs, lhs.y - rhs, lhs.z - rhs);
	}

	template <class X>
	DEVICE_CALLABLE
		inline r3d<X> operator-(const X &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(lhs - rhs.x, lhs - rhs.y, lhs - rhs.z);
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
		return r3d<X>(lhs.x / rhs, lhs.y / rhs, lhs.z / rhs);
	}

	template <class X>
	DEVICE_CALLABLE
		inline r3d<X> fmin(const r3d<X> &lhs, const r3d<X> &rhs)
	{
		return r3d<X>(::fmin(lhs.x, rhs.x), ::fmin(lhs.y, rhs.y), ::fmin(lhs.z, rhs.z));
	}

	template <class X>
	DEVICE_CALLABLE
		inline r3d<X> fmax(const r3d<X> &lhs, const r3d<X> &rhs)
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
		return r / r.module();
	}

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
		return acos(dot(lhs, rhs) / (lhs.module()*rhs.module()));
	}

	/************************************************************************/
	/********************************3x3 matrix******************************/
	/************************************************************************/

	template <class T>
	struct M3x3
	{
		using value_type = T;

		const int m;
		const int n;
		const int m_size;

		std::vector<T> data;

		DEVICE_CALLABLE
			inline M3x3() : m(3), n(3), m_size(m*n)
		{
			data.resize(m_size);
		}

		template<class TVector>
		DEVICE_CALLABLE
			inline M3x3(const TVector &v) : m(3), n(3), m_size(m*n)
		{
			data.resize(m_size);

			if (v.size() >= m_size)
			{
				std::copy(v.begin(), v.begin() + m_size, data.begin());
			}
		}

		int size() const
		{
			return m_size;
		}

		DEVICE_CALLABLE FORCE_INLINE
			T& operator[](const int idx) { return data[idx]; }

		DEVICE_CALLABLE FORCE_INLINE
			const T& operator[](const int idx) const { return data[idx]; }

		DEVICE_CALLABLE FORCE_INLINE
			T& operator()(const int ir, const int ic) { return data[ir - 1 + m*(ic - 1)]; }

		DEVICE_CALLABLE FORCE_INLINE
			const T& operator()(const int ir, const int ic) const { return data[ir - 1 + m*(ic - 1)]; }

		DEVICE_CALLABLE
			inline M3x3<T>& operator+=(const T &r)
		{
			for (auto ik = 0; ik < m_size; ik++)
			{
				data[ik] += r;
			}

			return *this;
		}

		DEVICE_CALLABLE
			inline M3x3<T>& operator+=(const M3x3<T> &r)
		{
			for (auto ik = 0; ik < m_size; ik++)
			{
				data[ik] += r.data[ik];
			}

			return *this;
		}

		DEVICE_CALLABLE
			inline M3x3<T>& operator-=(const T &r)
		{
			for (auto ik = 0; ik < m_size; ik++)
			{
				data[ik] -= r;
			}

			return *this;
		}

		DEVICE_CALLABLE
			inline M3x3<T>& operator-=(const M3x3<T> &r)
		{
			for (auto ik = 0; ik < m_size; ik++)
			{
				data[ik] -= r.data[ik];
			}

			return *this;
		}

		DEVICE_CALLABLE
			inline M3x3<T>& operator*=(const T &r)
		{
			for (auto ik = 0; ik < m_size; ik++)
			{
				data[ik] *= r;
			}

			return *this;
		}

		DEVICE_CALLABLE
			inline M3x3<T>& operator*=(const M3x3<T> &r)
		{
			T a11 = data[0] * r.data[0] + data[3] * r.data[1] + data[6] * r.data[2];
			T a21 = data[0] * r.data[3] + data[3] * r.data[4] + data[6] * r.data[5];
			T a31 = data[0] * r.data[6] + data[3] * r.data[7] + data[6] * r.data[8];

			T a12 = data[1] * r.data[0] + data[4] * r.data[1] + data[7] * r.data[2];
			T a22 = data[1] * r.data[3] + data[4] * r.data[4] + data[7] * r.data[5];
			T a32 = data[1] * r.data[6] + data[4] * r.data[7] + data[7] * r.data[8];

			T a13 = data[2] * r.data[0] + data[5] * r.data[1] + data[8] * r.data[2];
			T a23 = data[2] * r.data[3] + data[5] * r.data[4] + data[8] * r.data[5];
			T a33 = data[2] * r.data[6] + data[5] * r.data[7] + data[8] * r.data[8];

			data[0] = a11;
			data[1] = a21;
			data[2] = a31;

			data[3] = a12;
			data[4] = a22;
			data[5] = a32;

			data[6] = a13;
			data[7] = a23;
			data[8] = a33;

			return *this;
		}

		DEVICE_CALLABLE
			inline M3x3<T>& operator/=(const T r)
		{
			for (auto ik = 0; ik < m_size; ik++)
			{
				data[ik] /= r;
			}

			return *this;
		}
	};

	template <class X>
	DEVICE_CALLABLE
		inline M3x3<X> operator+(const M3x3<X> &lhs, const X &rhs)
	{
		M3x3<X> M;

		for (auto ik = 0; ik < M.size(); ik++)
		{
			M[ik] = lhs[ik] + rhs;
		}

		return M;
	}

	template <class X>
	DEVICE_CALLABLE
		inline M3x3<X> operator+(const X &lhs, const M3x3<X> &rhs)
	{
		M3x3<X> M;

		for (auto ik = 0; ik < M.size(); ik++)
		{
			M[ik] = lhs + rhs[ik];
		}

		return M;
	}

	template <class X>
	DEVICE_CALLABLE
		inline M3x3<X> operator+(const M3x3<X> &lhs, const M3x3<X> &rhs)
	{
		M3x3<X> M;

		for (auto ik = 0; ik < M.size(); ik++)
		{
			M[ik] = lhs[ik] + rhs[ik];
		}

		return M;
	}

	template <class X>
	DEVICE_CALLABLE
		inline M3x3<X> operator-(const M3x3<X> &lhs, const X &rhs)
	{
		M3x3<X> M;

		for (auto ik = 0; ik < M.size(); ik++)
		{
			M[ik] = lhs[ik] - rhs;
		}

		return M;
	}

	template <class X>
	DEVICE_CALLABLE
		inline M3x3<X> operator-(const X &lhs, const M3x3<X> &rhs)
	{
		M3x3<X> M;

		for (auto ik = 0; ik < M.size(); ik++)
		{
			M[ik] = lhs - rhs[ik];
		}

		return M;
	}

	template <class X>
	DEVICE_CALLABLE
		inline M3x3<X> operator-(const M3x3<X> &lhs, const M3x3<X> &rhs)
	{
		M3x3<X> M;

		for (auto ik = 0; ik < M.size(); ik++)
		{
			M[ik] = lhs[ik] - rhs[ik];
		}

		return M;
	}

	template <class X>
	DEVICE_CALLABLE
		inline M3x3<X> operator*(const M3x3<X> &lhs, const X &rhs)
	{
		M3x3<X> M;

		for (auto ik = 0; ik < M.size(); ik++)
		{
			M[ik] = lhs[ik] * rhs;
		}

		return M;
	}

	template <class X>
	DEVICE_CALLABLE
		inline M3x3<X> operator*(const X &lhs, const M3x3<X> &rhs)
	{
		M3x3<X> M;

		for (auto ik = 0; ik < M.size(); ik++)
		{
			M[ik] = lhs*rhs[ik];
		}

		return M;
	}

	template <class X>
	DEVICE_CALLABLE
		inline r3d<X> operator*(const M3x3<X> &lhs, const r3d<X> &rhs)
	{
		r3d<X> V;

		V.x = lhs[0] * rhs.x + lhs[3] * rhs.y + lhs[6] * rhs.z;
		V.y = lhs[1] * rhs.x + lhs[4] * rhs.y + lhs[7] * rhs.z;
		V.z = lhs[2] * rhs.x + lhs[5] * rhs.y + lhs[8] * rhs.z;

		return V;
	}

	template <class X>
	DEVICE_CALLABLE
		inline r3d<X> operator*(const r3d<X> &lhs, const M3x3<X> &rhs)
	{
		r3d<X> V;

		V.x = lhs[0] * rhs.x + lhs[1] * rhs.y + lhs[2] * rhs.z;
		V.y = lhs[3] * rhs.x + lhs[4] * rhs.y + lhs[5] * rhs.z;
		V.z = lhs[6] * rhs.x + lhs[7] * rhs.y + lhs[8] * rhs.z;

		return V;
	}

	template <class X>
	DEVICE_CALLABLE
		inline M3x3<X> operator*(const M3x3<X> &lhs, const M3x3<X> &rhs)
	{
		M3x3<X> M;

		M[0] = lhs[0] * rhs[0] + lhs[3] * rhs[1] + lhs[6] * rhs[2];
		M[1] = lhs[0] * rhs[3] + lhs[3] * rhs[4] + lhs[6] * rhs[5];
		M[2] = lhs[0] * rhs[6] + lhs[3] * rhs[7] + lhs[6] * rhs[8];

		M[3] = lhs[1] * rhs[0] + lhs[4] * rhs[1] + lhs[7] * rhs[2];
		M[4] = lhs[1] * rhs[3] + lhs[4] * rhs[4] + lhs[7] * rhs[5];
		M[5] = lhs[1] * rhs[6] + lhs[4] * rhs[7] + lhs[7] * rhs[8];

		M[6] = lhs[2] * rhs[0] + lhs[5] * rhs[1] + lhs[8] * rhs[2];
		M[7] = lhs[2] * rhs[3] + lhs[5] * rhs[4] + lhs[8] * rhs[5];
		M[8] = lhs[2] * rhs[6] + lhs[5] * rhs[7] + lhs[8] * rhs[8];

		return M;
	}

	template <class X>
	DEVICE_CALLABLE
		inline M3x3<X> operator/(const M3x3<X> &lhs, const X &rhs)
	{
		M3x3<X> M;

		for (auto ik = 0; ik < M.size(); ik++)
		{
			M[ik] = lhs[ik] / rhs;
		}

		return M;
	}

	template <class X>
	DEVICE_CALLABLE
		inline M3x3<X> operator/(const X &lhs, const M3x3<X> &rhs)
	{
		M3x3<X> M;

		for (auto ik = 0; ik < M.size(); ik++)
		{
			M[ik] = lhs / lhs[ik];
		}

		return M;
	}

	template <class X>
	DEVICE_CALLABLE
	inline M3x3<X> fmin(const M3x3<X> &lhs, const M3x3<X> &rhs)
	{
		M3x3<X> M;

		for (auto ik = 0; ik < M.size(); ik++)
		{
			M[ik] = ::fmin(lhs[ik], rhs[ik]);
		}

		return M;
	}

	template <class X>
	DEVICE_CALLABLE
	inline M3x3<X> fmax(const M3x3<X> &lhs, const M3x3<X> &rhs)
	{
		M3x3<X> M;

		for (auto ik = 0; ik < M.size(); ik++)
		{
			M[ik] = ::fmax(lhs[ik], rhs[ik]);
		}

		return M;
	}

	template <typename X>
	DEVICE_CALLABLE
	inline X det(const M3x3<X>& M)
	{
		X d = M(1, 1)*M(2, 2)*M(3, 3) + M(1, 2)*M(2, 3)*M(3, 1) + M(1, 3)*M(2, 1)*M(3, 2);
		d -= M(1, 3)*M(2, 2)*M(3, 1) + M(1, 2)*M(2, 1)*M(3, 3) + M(1, 1)*M(2, 3)*M(3, 2);

		return d;
	}

	template <typename X>
	DEVICE_CALLABLE
	inline M3x3<X> inv(const M3x3<X>& M)
	{
		M3x3<X> Mo;

		Mo(1, 1) = M(2, 2)*M(3, 3) - M(2, 3)*M(3, 2);
		Mo(2, 1) = M(2, 3)*M(3, 1) - M(2, 1)*M(3, 3);
		Mo(3, 1) = M(2, 1)*M(3, 2) - M(2, 2)*M(3, 1);

		Mo(1, 2) = M(1, 3)*M(3, 2) - M(1, 2)*M(3, 3);
		Mo(2, 2) = M(1, 1)*M(3, 3) - M(1, 3)*M(3, 1);
		Mo(3, 2) = M(1, 2)*M(3, 1) - M(1, 1)*M(3, 2);

		Mo(1, 3) = M(1, 2)*M(2, 3) - M(1, 3)*M(2, 2);
		Mo(3, 3) = M(1, 1)*M(2, 2) - M(1, 2)*M(2, 1);
		Mo(2, 3) = M(1, 3)*M(2, 1) - M(1, 1)*M(2, 3);

		Mo /= det(M);

		return Mo;
	}
} // namespace mt
#endif