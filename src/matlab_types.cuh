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
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MATLAB_TYPES_H
#define MATLAB_TYPES_H

#include <type_traits>
#include "math.cuh"
#include "types.cuh"
#include "stream.cuh"

namespace multem
{
	template<class T>
	struct complex_s
	{
	public:
		complex_s():m_real(nullptr), m_imag(nullptr){}

		template <class U> 
		inline complex_s<T>& operator=(const complex<U> & z)
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

		template<class U>
		inline operator complex<U>() const 
		{
			return complex<U>(*m_real, *m_imag);
		}

		inline void real(const T &re){ *m_real = re; }

		inline void imag(const T &im){ *m_imag = im; }

		inline T real() const { return *m_real; }

		inline T imag() const { return *m_imag; }

		template <class U>
		inline complex_s<T>& operator+=(const complex<U> &z)
		{
			real(real()+z.real());
			imag(imag()+z.imag());
			return *this;
		}

	private:
		T *m_real;
		T *m_imag;
	};

	template<class T>
	std::ostream& operator<<(std::ostream& out, const complex_s<T>& z){
		return out << "("<< z.real() << ", " << z.imag() << ")";
	}

	/*********************pointer to double matrix*****************/
	struct rmatrix_r
	{
	public:
		using value_type = double;
		using size_type = std::size_t;
		static const eDevice device = e_host;

		int size;
		int rows;
		int cols;
		double *real;
		rmatrix_r(): size(0), rows(0), cols(0), real(nullptr){}

		double& operator[](const int i){ return real[i]; }
		const double& operator[](const int i) const { return real[i]; }
		host_vector<double>::iterator begin() const { return real; };
		host_vector<double>::iterator end() const { return (real + size); };

		void resize(const size_type &new_size)
		{
			size = new_size;
			delete [] real;
			real = new double [new_size];
		}

		void clear()
		{
			size = 0;
			rows = 0;
			cols = 0;
			delete [] real; 
			real = nullptr;
		}
	};

	/*********************pointer to complex matrix****************/
	struct rmatrix_c
	{
	public:
		using value_type = complex<double>;
		using size_type = std::size_t;
		static const eDevice device = e_host;

		int size;
		int rows;
		int cols;
		double *real;
		double *imag;
		rmatrix_c(): size(0), rows(0), cols(0), real(nullptr), imag(nullptr){}

		//complex_s<double>& operator[](const int i)
		//{ 
		//	m_data(real[i], imag[i]);
		//	return m_data; 
		//}

		//complex<double> operator[](const int i) const
		//{ 
		//	return m_data; 
		//}

		//thrust::zip_iterator<thrust::tuple<double*, double*>> begin() const 
		//{ 
		//	return thrust::make_zip_iterator(thrust::make_tuple(real, imag)); 
		//};

		//thrust::zip_iterator<thrust::tuple<double*, double*>> end() const 
		//{ 
		//	return thrust::make_zip_iterator(thrust::make_tuple(real+size, imag+size)); 
		//};

		void resize(const size_type &new_size)
		{
			size = new_size;
			delete [] real;
			real = new double [new_size];
			delete [] imag;
			imag = new double [new_size];
		}

		void clear()
		{
			size = 0;
			rows = 0;
			cols = 0;
			delete [] real; 
			real = nullptr;
			delete [] imag; 
			imag = nullptr;
		}
	private:
		complex_s<double> m_data;

	};

	template<class T>
	struct is_rmatrix_r: std::integral_constant<bool, std::is_same<T, rmatrix_r>::value> {};

	template<class T>
	struct is_rmatrix_c: std::integral_constant<bool, std::is_same<T, rmatrix_c>::value> {};

	template <class T, class U>
	using enable_if_rmatrix_r = typename std::enable_if<is_rmatrix_r<T>::value, U>::type;

	template <class T, class U>
	using enable_if_rmatrix_c = typename std::enable_if<is_rmatrix_c<T>::value, U>::type;

	template <class T, class U>
	using enable_if_not_rmatrix_c = typename std::enable_if<!is_rmatrix_c<T>::value, U>::type;

	template <class T, class U>
	using enable_if_rmatrix = typename std::enable_if<is_rmatrix_r<T>::value || is_rmatrix_c<T>::value, U>::type;

	template<class TGrid>
	void fft2_shift(Stream<e_host> &stream, const TGrid &grid, rmatrix_c &M_io)
	{
		auto fft2_shift = [&](const multem::Range &range)
		{
			for(auto ix = range.ix_0; ix < range.ix_e; ix++)
			{
				for(auto iy = range.iy_0; iy < range.iy_e; iy++)
				{
					int ixy = grid.ind_col(ix, iy); 
					int ixy_shift = grid.ind_col(grid.nxh+ix, grid.nyh+iy);
					thrust::swap(M_io.real[ixy], M_io.real[ixy_shift]);
					thrust::swap(M_io.imag[ixy], M_io.imag[ixy_shift]);

					ixy = grid.ind_col(ix, grid.nyh+iy); 
					ixy_shift = grid.ind_col(grid.nxh+ix, iy);
					thrust::swap(M_io.real[ixy], M_io.real[ixy_shift]);
					thrust::swap(M_io.imag[ixy], M_io.imag[ixy_shift]);
				}
			}
		};

		stream.set_n_act_stream(grid.nxh);
		stream.set_grid(grid.nxh, grid.nyh);
		stream.exec(fft2_shift);
	}

	void fill(rmatrix_c &M_io, complex<double> value_i)
	{
		for(auto ixy = 0; ixy < M_io.size; ixy++)
		{
			M_io.real[ixy] = value_i.real();
			M_io.imag[ixy] = value_i.imag();
		}
	}

	template<class T>
	void assign(rmatrix_c &M_i, host_vector<T> &M_o, host_vector<T> *M_i_h=nullptr)
	{
		M_o.resize(M_i.size);
		for(auto ixy = 0; ixy < M_o.size(); ixy++)
		{
			M_o[ixy].real(M_i.real[ixy]);
			M_o[ixy].imag(M_i.imag[ixy]);
		}
	}

	template<class T>
	void assign(rmatrix_c &M_i, device_vector<T> &M_o, host_vector<T> *M_i_h=nullptr)
	{
		host_vector<T> M_h;
		M_i_h = (M_i_h==nullptr)?&M_h:M_i_h;

		assign(M_i, *M_i_h);
		M_o.assign(M_i_h->begin(), M_i_h->end());
	}

	void assign_square(rmatrix_c &M_i, rmatrix_r &M_o)
	{
		for(auto i = 0; i < M_i.size; i++)
		{
			M_o[i] = M_i.real[i]*M_i.real[i] + M_i.imag[i]*M_i.imag[i];
		}
	}

	template<class TGrid, class T>
	void copy_to_host(Stream<e_host> &stream, const TGrid &grid, host_vector<T> &M_i, rmatrix_c &M_o, host_vector<T> *M_i_h=nullptr)
	{
		auto copy_to_host = [&](const multem::Range &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				M_o.real[ixy] = M_i[ixy].real();
				M_o.imag[ixy] = M_i[ixy].imag();
			}
		};

		stream.set_n_act_stream(grid.nxy());
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(copy_to_host);
	}

	template<class TGrid, class T>
	void copy_to_host(Stream<e_host> &stream, const TGrid &grid, device_vector<T> &M_i, rmatrix_c &M_o, host_vector<T> *M_i_h=nullptr)
	{
		host_vector<T> M_h;
		M_i_h = (M_i_h==nullptr)?&M_h:M_i_h;

		M_i_h->assign(M_i.begin(), M_i.end());
		copy_to_host(stream, grid, *M_i_h, M_o);
	}

	template<class TGrid, class T>
	void add_scale_to_host(Stream<e_host> &stream, TGrid &grid, typename TGrid::value_type w_i, host_vector<T> &M_i, rmatrix_c &M_o, host_vector<T> *M_i_h=nullptr)
	{
		auto add_scale_to_host = [&](const multem::Range &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				M_o.real[ixy] += w_i*M_i[ixy].real();
				M_o.imag[ixy] += w_i*M_i[ixy].imag();
			}
		};

		stream.set_n_act_stream(grid.nxy());
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(add_scale_to_host);
	}

	template<class TGrid, class T>
	void add_scale_to_host(Stream<e_host> &stream, TGrid &grid, typename TGrid::value_type w_i, device_vector<T> &M_i, rmatrix_c &M_o, host_vector<T> *M_i_h=nullptr)
	{
		host_vector<T> M_h;
		M_i_h = (M_i_h==nullptr)?&M_h:M_i_h;

		M_i_h->assign(M_i.begin(), M_i.end());
		add_scale_to_host(stream, grid, w_i, *M_i_h, M_o);
	}

	template<class TGrid, class T>
	void add_scale_m2psi_psi_to_host(Stream<e_host> &stream, TGrid &grid, typename TGrid::value_type w_i, host_vector<T> &psi_i, rmatrix_r &m2psi_o, rmatrix_c &psi_o, host_vector<T> *psi_i_h=nullptr)
	{
		auto add_scale_m2psi_psi_to_host = [&](const multem::Range &range)
		{
			for(auto ixy = range.ixy_0; ixy < range.ixy_e; ixy++)
			{
				m2psi_o[ixy] += w_i*thrust::norm(psi_i[ixy]);
				psi_o.real[ixy] += w_i*psi_i[ixy].real();
				psi_o.imag[ixy] += w_i*psi_i[ixy].imag();
			}
		};

		stream.set_n_act_stream(grid.nxy());
		stream.set_grid(grid.nx, grid.ny);
		stream.exec(add_scale_m2psi_psi_to_host);
	}

	template<class TGrid, class T>
	void add_scale_m2psi_psi_to_host(Stream<e_host> &stream, TGrid &grid, typename TGrid::value_type w_i, device_vector<T> &psi_i, rmatrix_r &m2psi_o, rmatrix_c &psi_o, host_vector<T> *psi_i_h=nullptr)
	{
		host_vector<T> M_h;
		psi_i_h = (psi_i_h==nullptr)?&M_h:psi_i_h;

		psi_i_h->assign(psi_i.begin(), psi_i.end());
		add_scale_m2psi_psi_to_host(stream, grid, w_i, *psi_i_h, m2psi_o, psi_o);
	}

} // namespace multem

#endif