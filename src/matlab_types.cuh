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

#ifndef MATLAB_TYPES_H
#define MATLAB_TYPES_H

#include <type_traits>
#include "math.cuh"
#include "types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "cpu_fcns.hpp"

namespace mt
{
	template <class T>
	struct complex_s
	{
	public:
		complex_s(): m_real(nullptr), m_imag(nullptr){}

		template <class U> 
		inline complex_s<T>& operator = (const complex<U> & z)
		{
			*m_real = z.real();
			*m_imag = z.imag();
			return *this;
		}

		inline void operator()(T &real, T &imag)
		{
			m_real = &real;
			m_imag = &imag;
		}

		template <class U>
		inline operator complex<U>() const 
		{
			return complex<U>(*m_real, *m_imag);
		}

		inline void real(const T &re){ *m_real = re; }

		inline void imag(const T &im){ *m_imag = im; }

		inline T real() const { return *m_real; }

		inline T imag() const { return *m_imag; }

		template <class U>
		inline complex_s<T>& operator+= (const complex<U> &z)
		{
			real(real()+z.real());
			imag(imag()+z.imag());
			return *this;
		}

	private:
		T *m_real;
		T *m_imag;
	};

	template <class T>
	std::ostream& operator<<(std::ostream& out, const complex_s<T>& z){
		return out << "("<< z.real() << ", " << z.imag() <<")";
	}

	/*********************pointer to double matrix*****************/
	struct rmatrix_r
	{
	public:
		using value_type = double;
		using size_type = std::size_t;
		static const eDevice device = e_host;

		int m_size;
		int rows;
		int cols;
		double *real;
		double *imag;
		rmatrix_r(): m_size(0), rows(0), cols(0), real(nullptr), imag(nullptr){}

		double& operator[](const int i){ return real[i]; }
		const double& operator[](const int i) const { return real[i]; }
		host_vector<double>::iterator begin() const { return real; };
		host_vector<double>::iterator end() const { return (real + m_size); };

		size_type size() const
		{
			return m_size;
		}

		bool empty() const
		{
			return size()==0;
		}

		void resize(const size_type &new_size)
		{
			m_size = static_cast<int>(new_size);
			delete [] real;
			real = new double [new_size];
		}

		void clear()
		{
			m_size = 0;
			rows = 0;
			cols = 0;
			delete [] real; 
			real = nullptr;
		}

		template <class TInput_Iterator>
		void assign(TInput_Iterator first, TInput_Iterator last)
		{
			if(real!= nullptr)
			{
				thrust::copy(first, last, real);
			}
		}

		template <class T>
		T get(const int &i) const
		{
			return T(real[i]);
		}

		template <class T>
		typename std::enable_if<is_complex<T>::value, void>::type
		set(const int &i, const T &z)
		{
			real[i] = z.real();
		}

		template <class T>
		typename std::enable_if<!is_complex<T>::value, void>::type
		set(const int &i, const T &z)
		{
			real[i] = z;
		}

		void swap(const int &i, const int &j)
		{
			thrust::swap(real[i], real[j]);
		}
	};

	/*********************pointer to complex matrix****************/
	struct rmatrix_c
	{
	public:
		using value_type = complex<double>;
		using size_type = std::size_t;
		static const eDevice device = e_host;

		int m_size;
		int rows;
		int cols;
		double *real;
		double *imag;

		rmatrix_c(): m_size(0), rows(0), cols(0), real(nullptr), imag(nullptr){}

		complex<double> operator[](const int i) const
		{ 
			return complex<double>(real[i], imag[i]); 
		}

		size_type size() const
		{
			return m_size;
		}

		bool empty() const
		{
			return size()==0;
		}

		void resize(const size_type &new_size)
		{
			m_size = static_cast<int>(new_size);
			delete [] real;
			real = new double [new_size];
			delete [] imag;
			imag = new double [new_size];
		}

		void clear()
		{
			m_size = 0;
			rows = 0;
			cols = 0;
			delete [] real; 
			real = nullptr;
			delete [] imag; 
			imag = nullptr;
		}

		template <class T>
		typename std::enable_if<is_complex<T>::value, T>::type
		get(const int &i) const
		{
			return T(real[i], imag[i]);
		}

		template <class T>
		typename std::enable_if<!is_complex<T>::value, T>::type
		get(const int &i) const
		{
			return T(real[i]);
		}

		template <class T>
		typename std::enable_if<is_complex<T>::value, void>::type
		set(const int &i, const T &z)
		{
			real[i] = z.real();
			imag[i] = z.imag();
		}

		template <class T>
		typename std::enable_if<!is_complex<T>::value, void>::type
		set(const int &i, const T &z)
		{
			real[i] = z;
		}

		void swap(const int &i, const int &j)
		{
			thrust::swap(real[i], real[j]);
			thrust::swap(imag[i], imag[j]);
		}
	};

	/***********************Matlab traits**************************/
	template <class T>
	struct is_rmatrix_r: std::integral_constant<bool, std::is_same<T, rmatrix_r>::value> {};

	template <class T>
	struct is_rmatrix_c: std::integral_constant<bool, std::is_same<T, rmatrix_c>::value> {};

	template <class T>
	struct is_rmatrix: std::integral_constant<bool, is_rmatrix_r<T>::value || is_rmatrix_c<T>::value> {};

	template <class T1, class T2>
	struct is_rmatrix_and_rmatrix: std::integral_constant<bool, is_rmatrix<T1>::value && is_rmatrix<T2>::value> {};

	template <class T1, class T2>
	struct is_rmatrix_and_host_vector: std::integral_constant<bool, is_rmatrix<T1>::value && is_host_vector<T2>::value> {};

	template <class T1, class T2>
	struct is_rmatrix_and_device_vector: std::integral_constant<bool, is_rmatrix<T1>::value && is_device_vector<T2>::value> {};

	template <class T1, class T2>
	struct is_host_vector_and_rmatrix: std::integral_constant<bool, is_host_vector<T1>::value && is_rmatrix<T2>::value> {};

	template <class T1, class T2>
	struct is_device_vector_and_rmatrix: std::integral_constant<bool, is_device_vector<T1>::value && is_rmatrix<T2>::value> {};

	template <class T, class U>
	using enable_if_rmatrix_r = typename std::enable_if<is_rmatrix_r<T>::value, U>::type;

	template <class T, class U>
	using enable_if_rmatrix_c = typename std::enable_if<is_rmatrix_c<T>::value, U>::type;

	template <class T, class U>
	using enable_if_rmatrix = typename std::enable_if<is_rmatrix<T>::value, U>::type;

	template <class T1, class T2, class U>
	using enable_if_rmatrix_and_rmatrix = typename std::enable_if<is_rmatrix_and_rmatrix<T1, T2>::value, U>::type;

	template <class T1, class T2, class U>
	using enable_if_rmatrix_and_host_vector = typename std::enable_if<is_rmatrix_and_host_vector<T1, T2>::value, U>::type;

	template <class T1, class T2, class U>
	using enable_if_rmatrix_and_device_vector = typename std::enable_if<is_rmatrix_and_device_vector<T1, T2>::value, U>::type;

	template <class T1, class T2, class U>
	using enable_if_host_vector_and_rmatrix = typename std::enable_if<is_host_vector_and_rmatrix<T1, T2>::value, U>::type;

	template <class T1, class T2, class U>
	using enable_if_device_vector_and_rmatrix = typename std::enable_if<is_device_vector_and_rmatrix<T1, T2>::value, U>::type;

	/***********************Matlab functions**************************/

	template <class TVector_1, class TVector_2>
	enable_if_rmatrix_and_host_vector<TVector_1, TVector_2, void>
	assign(TVector_1 &M_i, TVector_2 &M_o, Vector<Value_type<TVector_2>, e_host> *M_i_h = nullptr)
	{
		using value_type = Value_type<TVector_2>;

		M_o.resize(M_i.size());

		for(auto ixy = 0; ixy < M_o.size(); ixy++)
		{
			M_o[ixy] = M_i.template get<value_type>(ixy);
		}
	}

	template <class TVector_1, class TVector_2>
	enable_if_rmatrix_and_device_vector<TVector_1, TVector_2, void>
	assign(TVector_1 &M_i, TVector_2 &M_o, Vector<Value_type<TVector_2>, e_host> *M_i_h = nullptr)
	{
		Vector<Value_type<TVector_2>, e_host> M_h;
		M_i_h = (M_i_h == nullptr)?&M_h:M_i_h;

		assign(M_i, *M_i_h);
		M_o.assign(M_i_h->begin(), M_i_h->end());
	}

	template <class TVector_1, class TVector_2>
	typename std::enable_if<is_host_vector_and_rmatrix<TVector_1, TVector_2>::value && is_complex<Value_type<TVector_1>>::value, void>::type
	assign_real(TVector_1 &M_i, TVector_2 &M_o, Vector<Value_type<TVector_2>, e_host> *M_i_h = nullptr)
	{
		for(auto ixy = 0; ixy < M_o.size(); ixy++)
		{
			M_o[ixy] = M_i[ixy].real();
		}
	}
	
	template <class TVector>
	enable_if_rmatrix<TVector, void>
	fill(Stream<e_host> &stream, TVector &M_io, Value_type<TVector> value_i)
	{
		auto thr_fill = [&](const Range_2d &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				M_io.set(ixy, value_i);
			}
		};

		stream.set_n_act_stream(M_io.size());
		stream.set_grid(1, M_io.size());
		stream.exec(thr_fill);
	}

	template <class TVector_1, class TVector_2>
	enable_if_rmatrix_and_rmatrix<TVector_1, TVector_2, void>
	scale(Stream<e_host> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_scale = [&](const Range_2d &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				auto z = w_i*M_i.template get<value_type>(ixy);
				M_o.set(ixy, z);
			};
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_scale);
	}

	template <class TVector>
	enable_if_rmatrix<TVector, void>
	scale(Stream<e_host> &stream, Value_type<TVector> w_i, TVector &M_io)
	{
		scale(stream, w_i, M_io, M_io);
	}

	template <class TVector_1, class TVector_2>
	enable_if_rmatrix_and_rmatrix<TVector_1, TVector_2, void>
	square(Stream<e_host> &stream, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_square = [&](const Range_2d &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				auto z = ::norm(M_i.template get<value_type>(ixy));
				M_o.set(ixy, z);
			}
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_square);
	}

	template <class TVector_1, class TVector_2>
	enable_if_rmatrix_and_rmatrix<TVector_1, TVector_2, void>
	square_scale(Stream<e_host> &stream, Value_type<TVector_2> w_i, TVector_1 &M_i, TVector_2 &M_o)
	{
		using value_type = Value_type<TVector_2>;
		auto thr_square_scale = [&](const Range_2d &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				auto z = w_i*::norm(M_i.template get<value_type>(ixy));
				M_o.set(ixy, z);
			}
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_square_scale);
	}

	template <class TGrid, class TVector>
	enable_if_rmatrix<TVector, void>
	fft1_shift(TGrid &grid_1d, TVector &M_io)
	{
		for(auto ix = 0; ix < grid_1d.nxh; ix++)
		{
			int ix_shift = grid_1d.iRx_shift(ix);
			M_io.swap(ix, ix_shift);
		}
	}

	template <class TGrid, class TVector>
	enable_if_rmatrix<TVector, void>
	fft2_shift(Stream<e_host> &stream, TGrid &grid_2d, TVector &M_io)
	{
		auto krn_fft2_shift = [](const int &ix, const int &iy, const TGrid &grid_2d, TVector &M_io)
		{
			int ixy = grid_2d.ind_col(ix, iy); 
			int ixy_shift = grid_2d.ind_col(grid_2d.nxh+ix, grid_2d.nyh+iy);
			M_io.swap(ixy, ixy_shift);

			ixy = grid_2d.ind_col(ix, grid_2d.nyh+iy); 
			ixy_shift = grid_2d.ind_col(grid_2d.nxh+ix, iy);
			M_io.swap(ixy, ixy_shift);
		};

		auto thr_fft2_shift = [&](const Range_2d &range)
		{
			host_detail::matrix_iter(range, krn_fft2_shift, grid_2d, M_io);
		};

		stream.set_n_act_stream(grid_2d.nxh);
		stream.set_grid(grid_2d.nxh, grid_2d.nyh);
		stream.exec(thr_fft2_shift);
	}

	template <class TVector>
	enable_if_rmatrix_r<TVector, Value_type<TVector>>
	sum(Stream<e_host> &stream, TVector &M_i)
	{
		using value_type = Value_type<TVector>;

		value_type sum_total = 0;
		auto thr_sum = [&](const Range_2d &range)
		{
			auto sum_partial = thrust::reduce(M_i.begin()+range.ixy_0, M_i.begin()+range.ixy_e);

			stream.stream_mutex.lock();
			sum_total += sum_partial;
			stream.stream_mutex.unlock();
		};

		stream.set_n_act_stream(M_i.size());
		stream.set_grid(1, M_i.size());
		stream.exec(thr_sum);

		return sum_total;
	}

	template <class TVector>
	enable_if_rmatrix_r<TVector, Value_type_r<TVector>>
	mean(Stream<e_host> &stream, TVector &M_i)
	{
		return sum(stream, M_i)/M_i.size();
	}

	/***************************************************************************/
	/***************************************************************************/
	template <class TVector_i, class TVector_o>
	enable_if_host_vector_and_rmatrix<TVector_i, TVector_o, void>
	copy_to_host(Stream<e_host> &stream, TVector_i &M_i, 
	TVector_o &M_o, TVector_i *M_i_h = nullptr)
	{
		auto thr_copy_to_host = [&](const Range_2d &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				M_o.set(ixy, M_i[ixy]);
			}
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_copy_to_host);
	}

	template <class TVector_i, class TVector_o>
	enable_if_host_vector_and_rmatrix<TVector_i, TVector_o, void>
	add_scale_to_host(Stream<e_host> &stream, Value_type<TVector_i> w_i, 
	TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_host> *M_i_h = nullptr)
	{
		using value_type = Value_type<TVector_o>;
		auto thr_add_scale_to_host = [&](const Range_2d &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				auto z = M_o.template get<value_type>(ixy) + value_type(w_i*M_i[ixy]);
				M_o.set(ixy, z);
			};
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_add_scale_to_host);
	}

	template <class TVector_i, class TVector_o>
	enable_if_host_vector_and_rmatrix<TVector_i, TVector_o, void>
	add_scale_square_to_host(Stream<e_host> &stream, Value_type<TVector_o> w_i, 
	TVector_i &M_i, TVector_o &M_o, Vector<Value_type<TVector_i>, e_host> *M_i_h = nullptr)
	{
		auto thr_add_scale_square_to_host = [&](const Range_2d &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				auto z = M_o.template get<Value_type<TVector_o>>(ixy) + w_i*thrust::norm(M_i[ixy]);
				M_o.set(ixy, z);
			};
		};

		stream.set_n_act_stream(M_o.size());
		stream.set_grid(1, M_o.size());
		stream.exec(thr_add_scale_square_to_host);
	}

	template <class TVector_c_i, class TVector_r_o, class TVector_c_o>
	enable_if_host_vector_and_rmatrix<TVector_c_i, TVector_c_o, void>
	add_scale_m2psi_psi_to_host(Stream<e_host> &stream, Value_type<TVector_r_o> w_i, 
	TVector_c_i &psi_i, TVector_r_o &m2psi_o, TVector_c_o &psi_o, Vector<Value_type<TVector_c_i>, e_host> *psi_i_h = nullptr)
	{
		using T_r = Value_type<TVector_r_o>;
		using T_c = Value_type<TVector_c_o>;
		auto thr_add_scale_m2psi_psi_to_host = [&](const Range_2d &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				auto z1 = m2psi_o.template get<T_r>(ixy) + T_r(w_i)*thrust::norm(psi_i[ixy]);
				auto z2 = psi_o.template get<T_c>(ixy) + T_c(w_i)*T_c(psi_i[ixy]);
				m2psi_o.set(ixy, z1);
				psi_o.set(ixy, z2);
			}
		};

		stream.set_n_act_stream(psi_o.size());
		stream.set_grid(1, psi_o.size());
		stream.exec(thr_add_scale_m2psi_psi_to_host);
	}
} // namespace mt

#endif