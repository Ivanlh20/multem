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

#include "r_2d.h"

/* forward declarations */
namespace mt
{
	template <class T>
	CGPU_EXEC
	dt_bool fcn_is_zero(const R_2d<T>& r);

	template <class T>
	CGPU_EXEC
	dt_bool fcn_is_nzero(const R_2d<T>& r);

	template <class T>
	CGPU_EXEC
	dt_bool fcn_is_one(const R_2d<T>& r);

	template <class T>
	CGPU_EXEC
	R_2d<T> floor(const R_2d<T>& r);

	template <class T>
	CGPU_EXEC
	R_2d<T> ceil(const R_2d<T>& r);

	template <class T>
	CGPU_EXEC
	R_2d<T> round(const R_2d<T>& r);

	template <class T>
	CGPU_EXEC
	T fmin(const R_2d<T>& r);

	template <class T>
	CGPU_EXEC
	T fmax(const R_2d<T>& r);

	template <class T>
	CGPU_EXEC
	T norm_2(const R_2d<T>& r);

	template <class T>
	CGPU_EXEC
	T norm(const R_2d<T>& r);
	
	template <class T>
	CGPU_EXEC
	R_2d<T> normalize(const R_2d<T>& r);
}

/* R_2d */
namespace mt
{
	/************************************* constructors ************************************/
	/* Default constructor */
	template <class T>
	CGPU_EXEC
	R_xtd<T, edim_2>::R_xtd(): x(T()), y(T()) {}

	/* constructor by initializer list */
	template <class T>
	template <class U>
	CGPU_EXEC
	R_xtd<T, edim_2>::R_xtd(const dt_init_list<U>& list)
	{
		auto ptr = list.begin();

		x = T(ptr[0]);
		y = T(ptr[1]);
	}

	/* constructor by pointer */
	template <class T>
	template <class U, class>
	CGPU_EXEC
	R_xtd<T, edim_2>::R_xtd(U* v, const dt_int64& n_r, const dt_int64& n_c, const dt_int64& idx, dt_int64 icol)
	{
		const auto ip = icol*n_r + idx;

		x = T(v[ip + 0*n_r]);
		y = (n_c>1)?T(v[ip + 1*n_r]):T(0);
	}

	/* constructor by pointer */
	template <class T>
	template <class U, class>
	CPU_EXEC
	R_xtd<T, edim_2>::R_xtd(U* v)
	{
		x = T(v[0]);
		y = T(v[1]);
	}

	/* constructor by pointer */
	template <class T>
	template <class U, class>
	CPU_EXEC
	R_xtd<T, edim_2>::R_xtd(U* v, const dt_int32& n_v)
	{
		x = (n_v>0)?T(v[0]):T(0);
		y = (n_v>1)?T(v[1]):T(0);
	}

	/* constructor a R_2d with an x */
	template <class T>
	CGPU_EXEC
	R_xtd<T, edim_2>::R_xtd(const T& x): x(x), y(T()) {}

	/* converting constructor: R_2d with an x */
	template <class T>
	template <class U>
	CGPU_EXEC
	R_xtd<T, edim_2>::R_xtd(const U& x): x(T(x)), y(T()) {}

	/* constructor a R_2d from its x and y parts */
	template <class T>
	CGPU_EXEC
	R_xtd<T, edim_2>::R_xtd(const T& x, const T& y): x(x), y(y) {}

	/* converting constructor: R_2d from its x and y parts */ 
	template <class T>
	template <class U>
	CGPU_EXEC
	R_xtd<T, edim_2>::R_xtd(const U& x, const U& y): x(T(x)), y(T(y)) {}

	/* Copy constructor */
	template <class T>
	CGPU_EXEC
	R_xtd<T, edim_2>::R_xtd(const R_xtd<T, edim_2>& r): x(r.x), y(r.y) {}

	/* Move constructor */
	template <class T>
	CGPU_EXEC
	R_xtd<T, edim_2>::R_xtd(R_xtd<T, edim_2>&& r): x(std::move(r.x)), y(std::move(r.y)) {}

	/* converting copy constructor */
	template <class T>
	template <class U>
	CGPU_EXEC
	R_xtd<T, edim_2>::R_xtd(const R_xtd<U, edim_2>& r): x(T(r.x)), y(T(r.y)) {}

	/******************************** assignment operators *********************************/
	/*Assign x and set y to 0 */
	template <class T>
	CGPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator=(const T& x)
	{
		this->x = x;
		this->y = T(0);

		return *this;
	}

	/* Convert and assign x and set y to 0 */
	template <class T>
	template <class U>
	CGPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator=(const U& x)
	{
		this->x = T(x);
		this->y = T(0);

		return *this;
	}

	/* copy assignment operator */
	template <class T>
	CGPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator=(const R_xtd<T, edim_2>& r)
	{
		x = r.x;
		y = r.y;

		return *this;
	}

	/* Move assignment operator */
	template <class T>
	CGPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator=(R_xtd<T, edim_2>&& r)
	{
		if (this != &r)
		{
			x = std::move(r.x);
			y = std::move(r.y);
		}

		return *this;
	}

	/* converting assignment operator */
	template <class T>
	template <class U>
	CGPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator=(const R_xtd<U, edim_2>& r)
	{
		x = T(r.x);
		y = T(r.y);

		return *this;
	}

#ifdef __CUDACC__
	/* Assignment operator from thrust::complex */
	template <class T>
	CGPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator=(const thrust::complex<T>& r)
	{
			x = r.real();
			y = r.imag();

			return *this;
	}

	/* converting assignment operator from thrust::complex */
	template <class T>
	template <class U>
	CGPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator=(const thrust::complex<U>& r)
	{
			x = T(r.real());
			y = T(r.imag());

			return *this;
	}
#endif

	/* Assignment operator from std::complex */
	template <class T>
	CPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator=(const std::complex<T>& r)
	{
			x = r.real();
			y = r.imag();

			return *this;
	}

	/* converting assignment operator from std::complex */
	template <class T>
	template <class U> 
	CPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator=(const std::complex<U>& r)
	{
		x = T(r.real());
		y = T(r.imag());

		return *this;
	}

	/******************* Compound assignment operators *******************/
	template <class T>
	template <class U>
	CGPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator+=(const R_xtd<U, edim_2>& r)
	{
		*this = *this + r;

		return *this;
	}

	template <class T>
	template <class U>
	CGPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator+=(const U& r)
	{
		*this = *this + r;

		return *this;
	}

	template <class T>
	template <class U>
	CGPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator-=(const R_xtd<U, edim_2>& r)
	{
		*this = *this - r;

		return *this;
	}

	template <class T>
	template <class U>
	CGPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator-=(const U& r)
	{
		*this = *this - r;

		return *this;
	}

	template <class T>
	template <class U>
	CGPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator*=(const R_xtd<U, edim_2>& r)
	{
		*this = *this*r;

		return *this;
	}

	template <class T>
	template <class U>
	CGPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator*=(const U& r)
	{
		*this = *this*r;

		return *this;
	}

	template <class T>
	template <class U>
	CGPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator/=(const R_xtd<U, edim_2>& r)
	{
		*this = *this/r;

		return *this;
	}

	template <class T>
	template <class U>
	CGPU_EXEC
	R_xtd<T, edim_2>& R_xtd<T, edim_2>::operator/=(const U& r)
	{
		*this = *this/r;

		return *this;
	}

	/************************** Other operators **************************/
	template <class T>
	CGPU_EXEC
	dt_bool R_xtd<T, edim_2>::fcn_is_zero() const
	{
		return mt::fcn_is_zero(*this);
	}

	template <class T>
	CGPU_EXEC
	dt_bool R_xtd<T, edim_2>::is_nzero() const
	{
		return mt::fcn_is_nzero(*this);
	}

	template <class T>
	CGPU_EXEC
	dt_bool R_xtd<T, edim_2>::is_one() const
	{
		return mt::fcn_is_one(*this);
	}

	template <class T>
	CGPU_EXEC
	R_xtd<T, edim_2> R_xtd<T, edim_2>::floor() const
	{
		return mt::floor(*this);
	}

	template <class T>
	CGPU_EXEC
	R_xtd<T, edim_2> R_xtd<T, edim_2>::ceil() const
	{
		return mt::ceil(*this);
	}

	template <class T>
	CGPU_EXEC
	R_xtd<T, edim_2> R_xtd<T, edim_2>::round() const
	{
		return mt::round(*this);
	}

	template <class T>
	CGPU_EXEC
	T R_xtd<T, edim_2>::min() const
	{
		return mt::fmin(*this);
	}

	template <class T>
	CGPU_EXEC
	T R_xtd<T, edim_2>::max() const
	{
		return mt::fmax(*this);
	}

	template <class T>
	CGPU_EXEC
	T R_xtd<T, edim_2>::norm_2() const
	{
		return mt::norm_2(*this);
	}

	template <class T>
	CGPU_EXEC
	T R_xtd<T, edim_2>::norm() const
	{
		return mt::norm(*this);
	}

	template <class T>
	CGPU_EXEC
	void R_xtd<T, edim_2>::normalize()
	{
		*this = mt::normalize(*this);
	}

	template <class T>
	CGPU_EXEC
	void R_xtd<T, edim_2>::fill(const T& v)
	{
		x = y = v;
	}

	template <class T>
	CGPU_EXEC
	void R_xtd<T, edim_2>::assign_nzero(const R_xtd<T, edim_2>& r)
	{
		x = r.x;
		y = (fcn_is_nzero(r.y))?r.y:r.x;
	}
}

/* unitary operators */
namespace mt
{
	template <class T>
	CGPU_EXEC
	R_2d<T> operator-(const R_2d<T>& r)
	{
		return {-r.x, -r.y};
	}

	template <class T>
	CGPU_EXEC
	dt_bool fcn_is_zero(const R_2d<T>& r)
	{
		return fcn_is_zero(r.x, r.y);
	}

	template <class T>
	CGPU_EXEC
	dt_bool fcn_is_nzero(const R_2d<T>& r)
	{
		return fcn_is_nzero(r.x, r.y);
	}

	template <class T>
	CGPU_EXEC
	dt_bool fcn_is_one(const R_2d<T>& r)
	{
		return fcn_is_one(r.x, r.y);
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> floor(const R_2d<T>& r)
	{
		return {::floor(r.x), ::floor(r.y)};
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> ceil(const R_2d<T>& r)
	{
		return {::ceil(r.x), ::ceil(r.y)};
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> round(const R_2d<T>& r)
	{
		return {::round(r.x), ::round(r.y)};
	}

	template <class T>
	CGPU_EXEC
	T fmin(const R_2d<T>& r)
	{
		return ::fmin(r.x, r.y);
	}

	template <class T>
	CGPU_EXEC
	T fmax(const R_2d<T>& r)
	{
		return ::fmax(r.x, r.y);
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> square(const R_2d<T>& r)
	{
		return {r.x*r.x, r.y*r.y};
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> sqrt(const R_2d<T>& r)
	{
		return {::sqrt(r.x), ::sqrt(r.y)};
	}

	template <class T>
	CGPU_EXEC
	T norm_2(const R_2d<T>& r)
	{
		return r.x*r.x + r.y*r.y;
	}

	template <class T>
	CGPU_EXEC
	T norm(const R_2d<T>& r)
	{
		return ::sqrt(norm_2(r));
	}
	
	template <class T>
	CGPU_EXEC
	R_2d<T> normalize(const R_2d<T>& r)
	{
		return fcn_div(r, r.norm());
	}

	template <class T>
	CGPU_EXEC
	R_2d<T> fill_nzero(R_2d<T> r)
	{
		r.x = r.x;
		r.y = (fcn_is_nzero(r.y))?r.y:r.x;

		return r;
	}
}

/* binary operators */
namespace mt
{
	template <class T, class U, class X = sel_lg_type<T, U>>
	CGPU_EXEC
	R_2d<X> operator+(const R_2d<T>& lhs, const R_2d<U>& rhs)
	{
		return {lhs.x + rhs.x, lhs.y + rhs.y};
	}

	template <class T, class U, class = enable_if_real<U>>
	CGPU_EXEC
	R_2d<T> operator+(const R_2d<T>& lhs, const U& rhs)
	{
		return {lhs.x + rhs, lhs.y + rhs};
	}

	template <class T, class U, class = enable_if_real<T>>
	CGPU_EXEC
	R_2d<U> operator+(const T& lhs, const R_2d<U>& rhs)
	{
		return {lhs+ rhs.x, lhs + rhs.y};
	}

	template <class T, class U, class X = sel_lg_type<T, U>>
	CGPU_EXEC
	R_2d<X> operator-(const R_2d<T>& lhs, const R_2d<U>& rhs)
	{
		return {lhs.x - rhs.x, lhs.y - rhs.y};
	}

	template <class T, class U, class = enable_if_real<U>>
	CGPU_EXEC
	R_2d<T> operator-(const R_2d<T>& lhs, const U& rhs)
	{
		return {lhs.x - rhs, lhs.y - rhs};
	}

	template <class T, class U, class = enable_if_real<T>>
	CGPU_EXEC
	R_2d<U> operator-(const T& lhs, const R_2d<U>& rhs) 
	{
		return {lhs - rhs.x, lhs - rhs.y};
	}

	template <class T, class U, class X = sel_lg_type<T, U>>
	CGPU_EXEC
	R_2d<X> operator*(const R_2d<T>& lhs, const R_2d<U>& rhs)
	{
		return {lhs.x*rhs.x, lhs.y*rhs.y};
	}

	template <class T, class U, class = enable_if_real<U>>
	CGPU_EXEC
	R_2d<T> operator*(const R_2d<T>& lhs, const U& rhs)
	{
		return {lhs.x*rhs, lhs.y*rhs};
	}

	template <class T, class U, class = enable_if_real<T>>
	CGPU_EXEC
	R_2d<U> operator*(const T& lhs, const R_2d<U>& rhs)
	{
		return {lhs*rhs.x, lhs*rhs.y};
	}

	template <class T, class U, class X = sel_lg_type<T, U>>
	CGPU_EXEC
	R_2d<X> operator/(const R_2d<T>& lhs, const R_2d<U>& rhs)
	{
		return {fcn_div(lhs.x, rhs.x), fcn_div(lhs.y, rhs.y)};
	}

	template <class T, class U, class = enable_if_real<U>>
	CGPU_EXEC
	R_2d<T> operator/(const R_2d<T>& lhs, const U& rhs)
	{
		return {fcn_div(lhs.x, rhs), fcn_div(lhs.y, rhs)};
	}

	template <class T, class U, class = enable_if_real<T>>
	CGPU_EXEC
		R_2d<U> operator/(const T& lhs, const R_2d<U>& rhs)
	{
		return {fcn_div(lhs, rhs.x), fcn_div(lhs, rhs.y)};
	}

	template <class T, class U, class X = sel_lg_type<T, U>>
	CGPU_EXEC
	R_2d<X> fmin(const R_2d<T>& lhs, const R_2d<U>& rhs)
	{
		return {::fmin(lhs.x, rhs.x), ::fmin(lhs.y, rhs.y)};
	}
	
	template <class T, class U, class X = sel_lg_type<T, U>>
	CGPU_EXEC
	R_2d<X> fmax(const R_2d<T>& lhs, const R_2d<U>& rhs)
	{
		return {::fmax(lhs.x, rhs.x), ::fmax(lhs.y, rhs.y)};
	}

	template <class T, class U, class X = sel_lg_type<T, U>>
	CGPU_EXEC
	X dot(const R_2d<T>& lhs, const R_2d<U>& rhs)
	{
		return lhs.x*rhs.x + lhs.y*rhs.y;
	}

	template <class T, class U, class X = sel_lg_type<T, U>>
	CGPU_EXEC
	X angle(const R_2d<T>& lhs, const R_2d<U>& rhs)
	{
		X m = ::sqrt(lhs.norm_2()*rhs.norm_2());
		return ::acos(fcn_div(dot(lhs, rhs), m));
	}
}

/* traits */
namespace mt
{
	template <class T>
	struct is_r_2d: std::integral_constant<dt_bool, std::is_same<T, R_2d<dt_bool>>::value ||
	std::is_same<T, R_2d<dt_int8>>::value || std::is_same<T, R_2d<dt_uint8>>::value || 
	std::is_same<T, R_2d<dt_int16>>::value || std::is_same<T, R_2d<dt_uint16>>::value ||
	std::is_same<T, R_2d<dt_int32>>::value || std::is_same<T, R_2d<dt_uint32>>::value ||
	std::is_same<T, R_2d<dt_int64>>::value || std::is_same<T, R_2d<dt_uint64>>::value ||
	std::is_same<T, R_2d<dt_float64>>::value || std::is_same<T, R_2d<dt_float32>>::value> {};

	template <class T, class U=void>
	using enable_if_r_2d = typename std::enable_if<is_r_2d<T>::value, U>::type;

	/***************************************************************************************/
	/******************************** initializer list *************************************/
	/***************************************************************************************/
	using dt_init_list_r_2d_f64 = std::initializer_list<R_2d<dt_float64>>;
}

/* data copy */
namespace mt
{
	// (dst, src): cpu -> cpu: R_2d
	template <class Td, class Ts>
	enable_if_real_real<Td, Ts, void>
	memcpy_pos_cpu_cpu(R_2d<Td>* pcpu_dst, const Ts* pcpu_src, dt_uint64 n_size, dt_uint64 icol=0, Ts* pcpu_jk = nullptr)
	{
		for (dt_uint64 ik = 0; ik < n_size; ik++)
		{
			const auto ip = icol*n_size + ik;
			pcpu_dst[ik]= R_2d<Td>(pcpu_src[ik + 0*n_size], pcpu_src[ik + 1*n_size]);
		}
	}

	// (dst, src): cpu: R_2d -> cpu
	template <class Td, class Ts>
	enable_if_real_real<Td, Ts, void>
	memcpy_pos_cpu_cpu(Td* pcpu_dst, const R_2d<Ts>* pcpu_src, dt_uint64 n_size, dt_uint64 icol=0, R_2d<Ts>* pcpu_jk = nullptr)
	{
		for (dt_uint64 ik = 0; ik < n_size; ik++)
		{
			const auto ip = icol*n_size + ik;
			const auto& r_2d = pcpu_src[ik];

			pcpu_dst[ip + 0*n_size] = Td(r_2d.x);
			pcpu_dst[ip + 1*n_size] = Td(r_2d.y);
		}
	}

#ifdef __CUDACC__
	// (dst, src): gpu: R_2d -> cpu
	template <class Td, class Ts>
	enable_if_real_real<Td, Ts, void>
	memcpy_pos_gpu_cpu(Td* pcpu_dst, const R_2d<Ts>* pgpu_src, dt_uint64 n_size, dt_uint64 icol=0, R_2d<Ts>* pcpu_jk = nullptr)
	{
		auto size_bytes = n_size*dt_uint64(sizeof(R_2d<Ts>));

		if (pcpu_jk == nullptr)
		{
			R_2d<Ts>* pcpu_t = new R_2d<Ts>[n_size];
			cudaMemcpy(pcpu_t, pgpu_src, size_bytes, cudaMemcpyDeviceToHost);
			memcpy_pos_cpu_cpu(pcpu_dst, pcpu_t, n_size, icol);
			delete[] pcpu_t;
		}
		else
		{
			cudaMemcpy(pcpu_jk, pgpu_src, size_bytes, cudaMemcpyDeviceToHost);
			memcpy_pos_cpu_cpu(pcpu_dst, pcpu_jk, n_size, icol);
		}
	}

	// (dst, src): cpu -> gpu: R_2d
	template <class Td, class Ts>
	enable_if_real_real<Td, Ts, void>
	memcpy_pos_cpu_gpu(R_2d<Td>* pgpu_dst, const Ts* pcpu_src, dt_uint64 n_size, dt_uint64 icol=0, R_2d<Td>* pcpu_jk = nullptr)
	{
		auto size_bytes = n_size*dt_uint64(sizeof(R_2d<Td>));

		if (pcpu_jk == nullptr)
		{
			R_2d<Td>* pcpu_t = new R_2d<Td>[n_size];
			memcpy_pos_cpu_cpu(pcpu_t, pcpu_src, n_size, icol);
			cudaMemcpy(pgpu_dst, pcpu_t, size_bytes, cudaMemcpyHostToDevice);
			delete[] pcpu_t;
		}
		else
		{
			memcpy_pos_cpu_cpu(pcpu_jk, pcpu_src, n_size, icol);
			cudaMemcpy(pgpu_dst, pcpu_jk, size_bytes, cudaMemcpyHostToDevice);
		}
	}
#endif
}

/* other operators */
namespace mt
{
	// distance from point to line
	template <class T>
	CGPU_EXEC_INL 
	T fcn_pt_2_ln_dist(const T& a, const T& b, const T& c, const R_2d<T>& r_0)
	{
		return fcn_pt_2_ln_dist(a, b, c, r_0.x, r_0.y);
	}

	// calculate intersection from point to line
	template <class T>
	CGPU_EXEC_INL 
	R_2d<T> fcn_pt_2_ln_intxn_2d(const T& a, const T& b, const T& c, const R_2d<T>& r_0)
	{
		R_2d<T> r;
		fcn_pt_2_ln_intxn_2d(a, b, c, r_0.x, r_0.y, r.x, r.y);

		return r;
	}
}