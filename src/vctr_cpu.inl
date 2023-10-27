/*
* This file is part of Multem.
* Copyright 2022 Ivan Lobato <Ivanlh20@gmail.com>
*
* Multem is destroy software: you can redistribute it and/or modify
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

#include "vctr_cpu.h"

/* vector copy */
namespace mt
{
	template <class T, eDev Dev> class Vctr;

	// (dst, src): cpu -> cpu
	template <class Td, class Ts>
	void vctr_cpy(Vctr<Td, edev_cpu>& vctr_cpu_dst, const Vctr<Ts, edev_cpu>& vctr_cpu_src, Ts* pvctr_cpu_jk = nullptr)
	{
		memcpy_cpu_cpu<Td, Ts>(vctr_cpu_dst.data(), vctr_cpu_src.data(), vctr_cpu_src.size(), pvctr_cpu_jk);
	}

#ifdef __CUDACC__
	// (dst, src): gpu -> cpu
	template <class Td, class Ts>
	void vctr_cpy(Vctr<Td, edev_cpu>& vctr_cpu_dst, const Vctr<Ts, edev_gpu>& vctr_gpu_src, Ts* pvctr_cpu_jk = nullptr)
	{
		memcpy_gpu_cpu<Td, Ts>(vctr_cpu_dst.data(), vctr_gpu_src.data(), vctr_gpu_src.size(), pvctr_cpu_jk);
	}

	// (dst, src): cpu -> gpu
	template <class Td, class Ts>
	void vctr_cpy(Vctr<Td, edev_gpu>& vctr_gpu_dst, const Vctr<Ts, edev_cpu>& vctr_cpu_src, Td* pvctr_cpu_jk = nullptr)
	{
		memcpy_cpu_gpu<Td, Ts>(vctr_gpu_dst.data(), vctr_cpu_src.data(), vctr_cpu_src.size(), pvctr_cpu_jk);
	}

	// (dst, src): gpu -> gpu
	template <class Td, class Ts>
	void vctr_cpy(Vctr<Td, edev_gpu>& vctr_gpu_dst, const Vctr<Ts, edev_gpu>& vctr_gpu_src, Td* pvctr_cpu_jk = nullptr)
	{
		memcpy_gpu_gpu<Td, Ts>(vctr_gpu_dst.data(), vctr_gpu_src.data(), vctr_gpu_src.size(), pvctr_cpu_jk);
	}
#endif
}

/* cpu vector */
namespace mt
{
	/************************************* constructors ************************************/
	template <class T>
	Vctr<T, edev_cpu>::Vctr(): m_data(nullptr), m_s0(0), m_s1(1), m_s2(1), m_s3(1), 
	m_size(shape_size()), m_capacity(0)
	{
		set_picth();
	}

	template <class T>
	Vctr<T, edev_cpu>::Vctr(const dt_init_list_f64& data): Vctr()
	{
		assign(data.begin(), data.end());
	}

	template <class T>
	Vctr<T, edev_cpu>::Vctr(size_type s0): Vctr()
	{
		resize({ s0, size_type(1), size_type(1), size_type(1) });
	}

	template <class T>
	Vctr<T, edev_cpu>::Vctr(size_type s0, const T& value): Vctr()
	{
		resize({ s0, size_type(1), size_type(1), size_type(1) }, value);
	}

	template <class T>
	Vctr<T, edev_cpu>::Vctr(const dt_shape_st<size_type>& shape): Vctr()
	{
		resize(shape);
	}

	template <class T>
	Vctr<T, edev_cpu>::Vctr(const dt_shape_st<size_type>& shape, const T& value): Vctr()
	{
		resize(shape, value);
	}

	/* copy constructor */
	template <class T>
	Vctr<T, edev_cpu>::Vctr(const Vctr<T, edev_cpu>& vctr): Vctr()
	{
		*this = vctr;
	}

	/* Move constructor */
	template <class T>
	Vctr<T, edev_cpu>::Vctr(Vctr<T, edev_cpu>&& vctr): Vctr()
	{
		*this = std::move(vctr);
	}

	/* converting constructor */
	template <class T>
	template <class U>
	Vctr<T, edev_cpu>::Vctr(const Vctr<U, edev_cpu>& vctr): Vctr()
	{
		assign(vctr);
	}

	template <class T>
	template <class U>
	Vctr<T, edev_cpu>::Vctr(U* first, U* last): Vctr()
	{
		assign(first, last);
	}

	template <class T>
	template <class U, class V, class>
	Vctr<T, edev_cpu>::Vctr(U *p, dt_int64 n_p, dt_int64 icol): Vctr()
	{
		set_pos_from_cpu_ptr(p, n_p, icol);
	}

	template <class T>
	template <class U>
	Vctr<T, edev_cpu>::Vctr(const std::vector<U>& vctr): Vctr()
	{
		assign(vctr);
	}

	// from cpu pVctr to Vctr
	template <class T>
	template <class U, class STU>
	Vctr<T, edev_cpu>::Vctr(const pVctr<U, edev_cpu, STU>& pvctr): Vctr()
	{
		assign(pvctr);
	}

#ifdef __CUDACC__
template <class T>
	template <class U>
	Vctr<T, edev_cpu>::Vctr(const Vctr<U, edev_gpu>& vctr): Vctr()
	{
		assign(vctr);
	}

	template <class T>
	template <class U>
	Vctr<T, edev_cpu>::Vctr(const thrust::host_vector<U>& vctr): Vctr()
	{
		assign(vctr);
	}

	template <class T>
	template <class U>
	Vctr<T, edev_cpu>::Vctr(const thrust::device_vector<U>& vctr): Vctr()
	{
		assign(vctr);
	}

	// from gpu pVctr to Vctr
	template <class T>
	template <class U, class STU>
	Vctr<T, edev_cpu>::Vctr(const pVctr<U, edev_gpu, STU>& pvctr): Vctr()
	{
		assign(pvctr);
	}
#endif
	template <class T>
	Vctr<T, edev_cpu>::~Vctr()
	{
		destroy();
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T>
	Vctr<T, edev_cpu>& Vctr<T, edev_cpu>::operator=(const Vctr<T, edev_cpu>& vctr)
	{
		if (vctr.m_data == nullptr)
		{
			destroy();
		}
		else if (this != &vctr)
		{
			delete[] m_data;

			m_data = new T[vctr.m_capacity];
			memcpy_cpu_cpu(m_data, vctr.data(), vctr.size());

			m_s0 = vctr.m_s0;
			m_s1 = vctr.m_s1;
			m_s3 = vctr.m_s3;
			m_size = vctr.m_size;
			m_capacity = vctr.m_capacity;	
			m_pitch_s1 = vctr.m_pitch_s1;
			m_pitch_s2 = vctr.m_pitch_s2;
			m_pitch_s3 = vctr.m_pitch_s3;
		}

		return *this;
	}

	/* Move assignment operator */
	template <class T>
	Vctr<T, edev_cpu>& Vctr<T, edev_cpu>::operator=(Vctr<T, edev_cpu>&& vctr)
	{
		if (vctr.m_data == nullptr)
		{
			destroy();
		}
		else if (this != &vctr)
		{
			delete[] m_data;

			m_data = vctr.m_data;
			m_s0 = vctr.m_s0;
			m_s1 = vctr.m_s1;
			m_s2 = vctr.m_s2;
			m_s3 = vctr.m_s3;
			m_size = vctr.m_size;
			m_capacity = vctr.m_capacity;
			m_pitch_s1 = vctr.m_pitch_s1;
			m_pitch_s2 = vctr.m_pitch_s2;
			m_pitch_s3 = vctr.m_pitch_s3;

			vctr.m_data = nullptr;
			vctr.m_s0 = 0;
			vctr.m_s1 = 0;
			vctr.m_s2 = 0;
			vctr.m_s3 = 0;
			vctr.m_size = 0;
			vctr.m_capacity = 0;
			vctr.m_pitch_s1 = 0;
			vctr.m_pitch_s2 = 0;
			vctr.m_pitch_s3 = 0;
		}

		return *this;
	}

	/* converting assignment operator */
	template <class T>
	template <class U>
	Vctr<T, edev_cpu>& Vctr<T, edev_cpu>::operator=(const Vctr<U, edev_cpu>& vctr) 
	{
		assign(vctr);
			
		return *this;
	}

	template <class T>
	template <class U>
	Vctr<T, edev_cpu>& Vctr<T, edev_cpu>::operator=(const std::vector<U>& vctr)
	{
		resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
		memcpy_cpu_cpu(m_data, vctr.data(), vctr.size());

		return *this;
	}

#ifdef __CUDACC__
template <class T>
	template <class U>
	Vctr<T, edev_cpu>& Vctr<T, edev_cpu>::operator=(const Vctr<U, edev_gpu>& vctr)
	{
		assign(vctr);
			
		return *this;
	}

	template <class T>
	template <class U>
	Vctr<T, edev_cpu>& Vctr<T, edev_cpu>::operator=(const thrust::host_vector<U>& vctr)
	{
		assign(vctr);
			
		return *this;
	}

	template <class T>
	template <class U>
	Vctr<T, edev_cpu>& Vctr<T, edev_cpu>::operator=(const thrust::device_vector<U>& vctr)
	{
		assign(vctr);
			
		return *this;
	}
#endif

	template <class T>
	template <class U>
	void Vctr<T, edev_cpu>::assign(const Vctr<U, edev_cpu>& vctr, U* pvctr_cpu)
	{
		if ((void*)this != (void*)&vctr)
		{
			this->allocate(vctr.shape());
			vctr_cpy(*this, vctr, pvctr_cpu);
		}
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_cpu>::assign(U* first, U* last)
	{
		if ((void*)m_data != (void*)first)
		{
			auto m_size_i = std::distance(first, last);
			if (m_size_i>0)
			{
				auto m_size = size_type(m_size_i);
				resize({ m_size, size_type(1), size_type(1), size_type(1) });
				memcpy_cpu_cpu(m_data, first, m_size);
			}
		}
	}

	// from cpu pVctr to Vctr
	template <class T>
	template <class U, class STU>
	void Vctr<T, edev_cpu>::assign(const pVctr<U, edev_cpu, STU>& pvctr)
	{
		resize(pvctr.shape());
		memcpy_cpu_cpu(m_data, pvctr.data(), m_size);
	}

#ifdef __CUDACC__
	template <class T>
	template <class U>
	void Vctr<T, edev_cpu>::assign(const thrust::device_ptr<U>& first, const thrust::device_ptr<U>& last, U* pvctr_cpu)
	{
		auto m_size_i = thrust::distance(first, last);

		if (m_size_i>0)
		{
			auto m_size = size_type(m_size_i);
			resize({ m_size, size_type(1), size_type(1), size_type(1) });
			memcpy_gpu_cpu(m_data, (U*)first.get(), m_size, pvctr_cpu);
		}
	}

	// from gpu pVctr to Vctr
	template <class T>
	template <class U, class STU>
	void Vctr<T, edev_cpu>::assign(const pVctr<U, edev_gpu, STU>& pvctr)
	{
		resize(pvctr.shape());
		memcpy_gpu_cpu(m_data, pvctr.data(), m_size);
	}
#endif

	template <class T>
	template <class U>
	void Vctr<T, edev_cpu>::assign(const std::vector<U>& vctr, U* pvctr_cpu)
	{
		resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
		memcpy_cpu_cpu(m_data, vctr.data(), vctr.size(), pvctr_cpu);
	}

#ifdef __CUDACC__
	template <class T>
	template <class U>
	void Vctr<T, edev_cpu>::assign(const Vctr<U, edev_gpu>& vctr, U* pvctr_cpu)
	{
		this->allocate(vctr.shape());
		vctr_cpy(*this, vctr, pvctr_cpu);
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_cpu>::assign(const thrust::host_vector<U>& vctr, U* pvctr_cpu)
	{
		resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
		memcpy_cpu_cpu(m_data, vctr.data(), vctr.size(), pvctr_cpu);
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_cpu>::assign(const thrust::device_vector<U>& vctr, U* pvctr_cpu)
	{
		resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
		memcpy_gpu_cpu(m_data, (U*)vctr.data().get(), vctr.size(), pvctr_cpu);
	}
#endif
	/**************** user define conversion operators *******************/
	template <class T>
	pVctr_cpu_32<T> Vctr<T, edev_cpu>::ptr_32() const
	{
		return pVctr_cpu_32<T>(*this);
	}

	template <class T>
	pVctr_cpu_64<T> Vctr<T, edev_cpu>::ptr_64() const
	{
		return pVctr_cpu_64<T>(*this);
	}

	template <class T>
	Vctr<T, edev_cpu>::operator pVctr_cpu_32<T>() const
	{
		return pVctr_cpu_32<T>(*this);
	}

	template <class T>
	Vctr<T, edev_cpu>::operator pVctr_cpu_64<T>() const
	{
		return pVctr_cpu_64<T>(*this);
	}

	/* user define conversion */
	template <class T>
	Vctr<T, edev_cpu>::operator std::vector<T>() const
	{
		std::vector<T> vctr(m_size);
		memcpy_cpu_cpu(vctr.data(), m_data, m_size);

		return vctr;
	}

	/* user define conversion in which T is the complemented precision */
	template <class T>
	template <class U>
	Vctr<T, edev_cpu>::operator std::vector<U>() const
	{
		std::vector<U> vctr(m_size);
		memcpy_cpu_cpu(vctr.data(), m_data, m_size);

		return vctr;
	}

#ifdef __CUDACC__
	/* user define conversion to output type std::vector<thrust::complex<dt_float32>> */
	template <class T>
	template <class U, class>
	Vctr<T, edev_cpu>::operator std::vector<thrust::complex<dt_float32>>() const
	{
		std::vector<thrust::complex<dt_float32>> vctr(m_size);
		memcpy_cpu_cpu(vctr.data(), m_data, m_size);

		return vctr;
	}

	/* user define conversion to output type std::vector<thrust::complex<dt_float64>> */
	template <class T>
	template <class U, class>
	Vctr<T, edev_cpu>::operator std::vector<thrust::complex<dt_float64>>() const
	{
		std::vector<thrust::complex<dt_float64>>vctr(m_size);
		memcpy_cpu_cpu(vctr.data(), m_data, m_size);

		return vctr;
	}

	/***************************************************************************************/
	/* user define conversion */
	template <class T>
	Vctr<T, edev_cpu>::operator thrust::host_vector<T>() const
	{
		thrust::host_vector<T> vctr(m_size);
		memcpy_cpu_cpu(vctr.data(), m_data, m_size);

		return vctr;
	}

	/* user define conversion in which T is the complemented precision */
	template <class T>
	template <class U>
	Vctr<T, edev_cpu>::operator thrust::host_vector<U>() const
	{
		thrust::host_vector<U> vctr(m_size);
		memcpy_cpu_cpu(vctr.data(), m_data, m_size);

		return vctr;
	}

	/* user define conversion to output type thrust::host_vector<std::complex<dt_float32>> */
	template <class T>
	template <class U, class>
	Vctr<T, edev_cpu>::operator thrust::host_vector<std::complex<dt_float32>>() const
	{
		thrust::host_vector<std::complex<dt_float32>> vctr(m_size);
		memcpy_cpu_cpu(vctr.data(), m_data, m_size);

		return vctr;
	}

	/* user define conversion to output type thrust::host_vector<std::complex<dt_float64>> */
	template <class T>
	template <class U, class>
	Vctr<T, edev_cpu>::operator thrust::host_vector<std::complex<dt_float64>>() const
	{
		thrust::host_vector<std::complex<dt_float64>> vctr(m_size);
		memcpy_cpu_cpu(vctr.data(), m_data, m_size);

		return vctr;
	}

	/***************************************************************************************/
	/* user define conversion */
	template <class T>
	Vctr<T, edev_cpu>::operator thrust::device_vector<T>() const
	{
		thrust::device_vector<T> vctr(m_size);
		memcpy_cpu_gpu((T*)vctr.data().get(), m_data, m_size);

		return vctr;
	}

	/* user define conversion in which T is the complemented precision */
	template <class T>
	template <class U>
	Vctr<T, edev_cpu>::operator thrust::device_vector<U>() const
	{
		thrust::device_vector<U> vctr(m_size);
		memcpy_cpu_gpu((U*)vctr.data().get(), m_data, m_size);

		return vctr;
	}
#endif
	/***************************************************************************************/
	template <class T>
	template <class U>
	void Vctr<T, edev_cpu>::cpy_to_cpu_ptr(U* pdata, size_type n_data, T* pvctr_cpu)
	{
		n_data = std::min(m_size, n_data);
		memcpy_cpu_cpu(pdata, m_data, n_data, pvctr_cpu);
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_cpu>::cpy_to_cpu_ptr(U* first, U* last, T* pvctr_cpu)
	{
		auto n_data = std::distance(first, last);
		memcpy_cpu_cpu(first, m_data, n_data, pvctr_cpu);
	}

#ifdef __CUDACC__
	template <class T>
	template <class U>
	void Vctr<T, edev_cpu>::cpy_to_gpu_ptr(U* pdata, size_type n_data, U* pvctr_cpu)
	{
		n_data = std::min(m_size, n_data);
		memcpy_cpu_gpu(pdata, m_data, n_data, pvctr_cpu);
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_cpu>::cpy_to_gpu_ptr(U* first, U* last, U* pvctr_cpu)
	{
		auto n_data = std::distance(first, last);
		memcpy_cpu_gpu(first, m_data, n_data, pvctr_cpu);
	}
#endif

	/***************************************************************************************/
	template <class T>
	template <class U, class V, class >
	void Vctr<T, edev_cpu>::cpy_real_to_cpu_ptr(U* pdata, size_type n_data, T* pvctr_cpu)
	{
		n_data = std::min(m_size, n_data);
		memcpy_real_cpu_cpu(pdata, m_data, n_data, pvctr_cpu);
	}

	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_cpu>::cpy_real_to_cpu_ptr(U* first, U* last, T* pvctr_cpu)
	{
		auto n_data = std::distance(first, last);
		memcpy_real_cpu_cpu(first, m_data, n_data, pvctr_cpu);
	}

#ifdef __CUDACC__
	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_cpu>::cpy_real_to_gpu_ptr(U* pdata, size_type n_data, U* pvctr_cpu)
	{
		n_data = std::min(m_size, n_data);
		memcpy_real_cpu_gpu(pdata, m_data, n_data, pvctr_cpu);
	}

	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_cpu>::cpy_real_to_gpu_ptr(U* first, U* last, U* pvctr_cpu)
	{
		auto n_data = std::distance(first, last);
		memcpy_real_cpu_gpu(first, m_data, n_data, pvctr_cpu);
	}
#endif

	/***************************************************************************************/
	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_cpu>::cpy_imag_to_cpu_ptr(U* pdata, size_type n_data, T* pvctr_cpu)
	{
		n_data = std::min(m_size, n_data);
		memcpy_imag_cpu_cpu(pdata, m_data, n_data, pvctr_cpu);
	}

	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_cpu>::cpy_imag_to_cpu_ptr(U* first, U* last, T* pvctr_cpu)
	{
		auto n_data = std::distance(first, last);
		memcpy_imag_cpu_cpu(first, m_data, n_data, pvctr_cpu);
	}

#ifdef __CUDACC__
	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_cpu>::cpy_imag_to_gpu_ptr(U* pdata, size_type n_data, U* pvctr_cpu)
	{
		n_data = std::min(m_size, n_data);
		memcpy_imag_cpu_gpu(pdata, m_data, n_data, pvctr_cpu);
	}

	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_cpu>::cpy_imag_to_gpu_ptr(U* first, U* last, U* pvctr_cpu)
	{
		auto n_data = std::distance(first, last);
		memcpy_imag_cpu_gpu(first, m_data, n_data, pvctr_cpu);
	}
#endif

	/***************************************************************************************/
	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_cpu>::set_pos_from_cpu_ptr(U *p, size_type n_p, dt_int64 icol)
	{
		resize({n_p});

		memcpy_pos_cpu_cpu(m_data, p, size(), icol);
	}

	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_cpu>::cpy_pos_to_cpu_ptr(U *p, size_type icol)
	{
		memcpy_pos_cpu_cpu(p, m_data, size(), icol);
	}

	/***************************************************************************************/
	template <class T>
	void Vctr<T, edev_cpu>::resize(const dt_shape_st<size_type>& shape)
	{
		this->allocate(shape);
	}

	template <class T>
	void Vctr<T, edev_cpu>::resize(const dt_shape_st<size_type>& shape, const T& value)
	{
		auto m_size_t = m_size;
		this->allocate(shape);
		std::fill(this->begin()+m_size_t, this->end(), value);
	}

	template <class T>
	void Vctr<T, edev_cpu>::reserve(const dt_shape_st<size_type>& shape)
	{
		this->allocate(shape, true);
	}

	template <class T>
	void Vctr<T, edev_cpu>::shrink_to_fit()
	{
		if (m_size < m_capacity)
		{
			if (m_size > 0)
			{
				/* allocate memory and transfer data */
				T* data_tmp = new T[m_size];
				memcpy_cpu_cpu(data_tmp, m_data, m_size);
				delete[] m_data;

				m_data = data_tmp;
				data_tmp = nullptr;
			}
			else
			{
				delete[] m_data;
				m_data = nullptr;
			}

			m_capacity = m_size;

			if (shape().dim()==2)
			{
				if (m_s1==1)
					m_s0 = m_size;

				if (m_s0==1)
					m_s1 = m_size;
			}
		}
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_cpu>::push_back(const U& val)
	{
		if (m_size >= m_capacity)
		{
			dt_shape_st<size_type> shape{m_s0, m_s1, m_s2, m_s3};

			if (m_s3 == 1)
			{
				if (m_s2 == 1)
				{
					if (m_s1 == 1)
						shape[0] = max(size_type(1), 2*m_size);
					else
						shape[1] *= size_type(2);
				}
				else
				{
					shape[2] *= size_type(2);
				}
			}
			else
			{
				shape[3] *= size_type(2);
			}

			this->allocate(shape, true);
		}

		m_data[m_size++] = T(val);
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_cpu>::push_back(const Vctr<U, edev_cpu>& vctr)
	{
		const size_type size_new = m_size + vctr.size();

		if (size_new >= m_capacity)
		{
			dt_shape_st<size_type> shape{m_s0, m_s1, m_s2, m_s3};

			if (m_s3 == 1)
			{
				if (m_s2 == 1)
				{
					if (m_s1 == 1)
						shape[0] = max(size_new, 2*shape[0]);
					else
						shape[1] = max(size_new/shape[0], 2*shape[1]);
				}
				else
				{
					shape[2] = max(size_new/(shape[0]*shape[1]), 2*shape[2]);
				}
			}
			else
			{
				shape[3] = max(size_new/(shape[0]*shape[1]*shape[2]), 2*shape[3]);
			}

			this->allocate(shape, true);
		}

		memcpy_cpu_cpu(m_data + m_size, vctr.data(), vctr.size());

		m_size = size_new;
	}

	template <class T>
	void Vctr<T, edev_cpu>::pop_back()
	{
		if (m_size >= 1)
		{
			m_size--;
		}
	}

	template <class T>
	void Vctr<T, edev_cpu>::fill(T val)
	{
		std::fill(this->begin(), this->end(), val);
	}

	/***************************************************************************************/
	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::s0() const
	{
		return m_s0;
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::s1() const
	{
		return m_s1;
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::s2() const
	{
		return m_s2;
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::s3() const
	{
		return m_s2;
	}			
			
	template <class T>
	dt_int32 Vctr<T, edev_cpu>::s0_32() const
	{
		return static_cast<dt_int32>(m_s0);
	}

	template <class T>
	dt_int32 Vctr<T, edev_cpu>::s1_32() const
	{
		return static_cast<dt_int32>(m_s1);
	}

	template <class T>
	dt_int32 Vctr<T, edev_cpu>::s2_32() const
	{
		return static_cast<dt_int32>(m_s2);
	}

	template <class T>
	dt_int32 Vctr<T, edev_cpu>::s3_32() const
	{
		return static_cast<dt_int32>(m_s3);
	}
			
	template <class T>
	dt_int64 Vctr<T, edev_cpu>::s0_64() const
	{
		return static_cast<dt_int64>(m_s0);
	}

	template <class T>
	dt_int64 Vctr<T, edev_cpu>::s1_64() const
	{
		return static_cast<dt_int64>(m_s1);
	}

	template <class T>
	dt_int64 Vctr<T, edev_cpu>::s2_64() const
	{
		return static_cast<dt_int64>(m_s2);
	}

	template <class T>
	dt_int64 Vctr<T, edev_cpu>::s3_64() const
	{
		return static_cast<dt_int64>(m_s3);
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::s0h() const
	{
		return m_s0/size_type(2);
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::s1h() const
	{
		return m_s1/size_type(2);
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::s2h() const
	{
		return m_s2/size_type(2);
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::s3h() const
	{
		return m_s3/size_type(2);
	}

	template <class T>
	dt_shape_st<typename Vctr<T, edev_cpu>::size_type> Vctr<T, edev_cpu>::shape() const
	{
		return {m_s0, m_s1, m_s2, m_s3};
	}

	template <class T>
	dt_shape_st<typename Vctr<T, edev_cpu>::size_type> Vctr<T, edev_cpu>::shape_2d_trs() const
	{
		return {m_s1, m_s0, m_s2, m_s3};
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::shape_size() const 
	{
		return m_s0*max(m_s1, size_type(1))*max(m_s2, size_type(1))*max(m_s3, size_type(1));
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::pitch_s1() const
	{
		return m_pitch_s1;
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::pitch_s2() const
	{
		return m_pitch_s2;
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::pitch_s3() const
	{
		return m_pitch_s3;
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::size() const
	{
		return m_size;
	}

	template <class T>
	dt_int32 Vctr<T, edev_cpu>::size_32() const
	{
		return static_cast<dt_int32>(m_size);
	}

	template <class T>
	dt_int64 Vctr<T, edev_cpu>::size_64() const
	{
		return static_cast<dt_int64>(m_size);
	}

	template <class T>
	iGrid_1d Vctr<T, edev_cpu>::igrid_1d() const
	{
		return {size_32()};
	}

	template <class T>
	iGrid_2d Vctr<T, edev_cpu>::igrid_2d() const
	{
		return {s1_32(), s0_32()};
	}

	template <class T>
	iGrid_3d Vctr<T, edev_cpu>::igrid_3d() const
	{
		return {s1_32(), s0_32(), s2_32()};
	}

	template <class T>
	iGrid_1d_64 Vctr<T, edev_cpu>::igrid_1d_64() const
	{
		return {size_64()};
	}

	template <class T>
	iGrid_2d_64 Vctr<T, edev_cpu>::igrid_2d_64() const
	{
		return {s1_64(), s0_64()};
	}

	template <class T>
	iGrid_3d_64 Vctr<T, edev_cpu>::igrid_3d_64() const
	{
		return {s1_64(), s0_64(), s2_64()};
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::capacity() const
	{
		return m_capacity;
	}

	template <class T>
	dt_bool Vctr<T, edev_cpu>::empty() const
	{
		return m_size == 0;
	}

	template <class T>
	dt_bool Vctr<T, edev_cpu>::is_1d() const
	{
		return (m_s0 == 1) || (m_s1 == 1);
	}

	template <class T>
	void Vctr<T, edev_cpu>::clear()
	{
		m_size = 0;
	}

	template <class T>
	void Vctr<T, edev_cpu>::clear_shrink_to_fit()
	{
		destroy();
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::sub_2_ind(const size_type& ix_0) const 
	{ 
		return ix_0;
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::sub_2_ind(const size_type& ix_0, const size_type& ix_1) const 
	{ 
		return ix_0 + m_s0*ix_1;
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::sub_2_ind(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2) const 
	{ 
		return ix_0 + m_s0*(ix_1 + m_s1*ix_2);
	}

	template <class T>
	typename Vctr<T, edev_cpu>::size_type Vctr<T, edev_cpu>::sub_2_ind(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2, const size_type& ix_3) const 
	{ 
		return ix_0 + m_s0*(ix_1 + m_s1*(ix_2 + m_s2*ix_3));
	}

	template <class T>
	T& Vctr<T, edev_cpu>::operator[](const size_type& iy)
	{ 
		return m_data[iy];
	}

	template <class T>
	const T& Vctr<T, edev_cpu>::operator[](const size_type& iy) const 
	{ 
		return m_data[iy];
	}

	template <class T>
	T& Vctr<T, edev_cpu>::operator()(const size_type& iy)
	{ 
		return m_data[iy];
	}

	template <class T>
	const T& Vctr<T, edev_cpu>::operator()(const size_type& iy) const 
	{ 
		return m_data[iy];
	}

	template <class T>
	T& Vctr<T, edev_cpu>::operator()(const size_type& ix_0, const size_type& ix_1)
	{ 
		return m_data[sub_2_ind(ix_0, ix_1)];
	}

	template <class T>
	const T& Vctr<T, edev_cpu>::operator()(const size_type& ix_0, const size_type& ix_1) const 
	{ 
		return m_data[sub_2_ind(ix_0, ix_1)];
	}

	template <class T>
	T& Vctr<T, edev_cpu>::operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2)
	{ 
		return m_data[sub_2_ind(ix_0, ix_1, ix_2)];
	}

	template <class T>
	const T& Vctr<T, edev_cpu>::operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2) const 
	{ 
		return m_data[sub_2_ind(ix_0, ix_1, ix_2)];
	}

	template <class T>
	T& Vctr<T, edev_cpu>::operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2, const size_type& ix_3)
	{ 
		return m_data[sub_2_ind(ix_0, ix_1, ix_2, ix_3)];
	}

	template <class T>
	const T& Vctr<T, edev_cpu>::operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2, const size_type& ix_3) const 
	{ 
		return m_data[sub_2_ind(ix_0, ix_1, ix_2, ix_3)];
	}

	template <class T>
	T* Vctr<T, edev_cpu>::begin()
	{ 
		return m_data;
	}

	template <class T>
	const T* Vctr<T, edev_cpu>::begin() const 
	{ 
		return m_data;
	}

	template <class T>
	T* Vctr<T, edev_cpu>::end()
	{ 
		return m_data + m_size;
	}

	template <class T>
	const T* Vctr<T, edev_cpu>::end() const
	{ 
		return m_data + m_size;
	}

	template <class T>
	T* Vctr<T, edev_cpu>::data()
	{
		return m_data;
	}

	template <class T>
	const T* Vctr<T, edev_cpu>::data() const
	{
		return m_data;
	}

	template <class T>
	template <class U>
	U Vctr<T, edev_cpu>::data_cast()
	{
		return reinterpret_cast<U>(m_data);
	}

	template <class T>
	template <class U>
	const U Vctr<T, edev_cpu>::data_cast() const
	{
		return reinterpret_cast<U>(m_data);
	}

	template <class T>
	T& Vctr<T, edev_cpu>::front()
	{ 
		return m_data[0];
	}

	template <class T>
	const T& Vctr<T, edev_cpu>::front() const 
	{ 
		return m_data[0];
	}

	template <class T>
	T& Vctr<T, edev_cpu>::back()
	{ 
		return m_data[m_size-1];
	}

	template <class T>
	const T& Vctr<T, edev_cpu>::back() const 
	{ 
		return m_data[m_size-1];
	}

	// set shape
	template <class T>
	void Vctr<T, edev_cpu>::set_shape(const dt_shape_st<size_type>& shape)
	{
		m_s0 = max(shape[0], size_type(0));
		m_s1 = max(shape[1], size_type(1));
		m_s2 = max(shape[2], size_type(1));
		m_s3 = max(shape[3], size_type(1));

		set_picth();
	}

	template <class T>
	void Vctr<T, edev_cpu>::trs_shape_2d()
	{
		set_shape({m_s1, m_s0, m_s2, m_s3});
	}

	template <class T>
	template <class U, eDev Dev>
	void Vctr<T, edev_cpu>::set_shape(const Vctr<U, Dev>& vctr, dt_bool bb_size)
	{
		this->set_shape(vctr.shape(), bb_size);
	}

#ifdef __CUDACC__
	FCNS_IMP_GPU_GRID_BLK_VCTR(template <class T>, Vctr<T COMMA edev_cpu>);
#endif

	template <class T>
	void Vctr<T, edev_cpu>::set_picth()
	{
		m_pitch_s1 = m_s0*sizeof(T);
		m_pitch_s2 = m_pitch_s1*m_s1;
		m_pitch_s3 = m_pitch_s2*m_s2;
	}

	template <class T>
	void Vctr<T, edev_cpu>::set_capacity(size_type size_r)
	{
		m_capacity = max(m_capacity, size_r);
	}

	template <class T>
	void Vctr<T, edev_cpu>::set_shape_cstr(dt_shape_st<size_type>& shape)
	{
		shape[0] = max(shape[0], size_type(0));
		shape[1] = max(shape[1], size_type(1));
		shape[2] = max(shape[2], size_type(1));
		shape[3] = max(shape[3], size_type(1));
	}

	// reallocate and copy memory
	template <class T>
	void Vctr<T, edev_cpu>::allocate(dt_shape_st<size_type> shape, dt_bool bb_reserve)
	{
		set_shape_cstr(shape);

		auto size_r = shape[0]*shape[1]*shape[2]*shape[3];
			
		if (size_r == 0)
		{
			return;
		}

		if (m_size == 0)
		{
			if (m_data != nullptr)
			{
				delete[] m_data;
			}

			m_data = new T[size_r];
		}
		else if (size_r>m_capacity)
		{
			/*allocate	memory and transfer data */
			T* data_tmp = new T[size_r];
			memcpy_cpu_cpu(data_tmp, m_data, m_size);

			if (m_data != nullptr)
			{
				delete[] m_data;
			}

			m_data = data_tmp;
			data_tmp = nullptr;
		}

		this->set_shape(shape);
		if (!bb_reserve)
		{
			m_size = size_r;
		}
		this->set_capacity(size_r);
	}

	// destroy memory on the device
	template <class T>
	void Vctr<T, edev_cpu>::init()
	{
		m_data = nullptr;
		m_s0 = 0; m_s1 = m_s2 = m_s3 = 1;
		m_pitch_s1 = m_pitch_s2 = m_pitch_s3 = 0;
		m_size = 0;
		m_capacity = 0;
	}

	// destroy memory on the device
	template <class T>
	void Vctr<T, edev_cpu>::destroy()
	{
		if (m_data != 0)
		{
			delete[] m_data;
		}

		init();
	}

}

/* traits */
namespace mt
{
	template <class T>
	struct is_vctr_cpu_r_2d: std::integral_constant<dt_bool, is_vctr_cpu<T>::value && is_r_2d<typename T::value_type>::value> {};

	template <class T, class U>
	struct is_vctr_cpu_r_2d_and_vctr_cpu: std::integral_constant<dt_bool, is_vctr_cpu_r_2d<T>::value && is_vctr_cpu<U>::value> {};

	/***************************************************************************************/
	template <class T, class U=void>
	using enable_if_vctr_cpu_r_2d = typename std::enable_if<is_vctr_cpu_r_2d<T>::value, U>::type;

	template <class T, class U, class V=void>
	using enable_if_vctr_cpu_r_2d_and_vctr_cpu = typename std::enable_if<is_vctr_cpu_r_2d_and_vctr_cpu<T, U>::value, V>::type;

	/***************************************************************************************/
	template <class T>
	struct is_vctr_cpu_r_3d: std::integral_constant<dt_bool, is_vctr_cpu<T>::value && is_r_3d<typename T::value_type>::value> {};

	template <class T, class U>
	struct is_vctr_cpu_r_3d_and_vctr_cpu: std::integral_constant<dt_bool, is_vctr_cpu_r_3d<T>::value && is_vctr_cpu<U>::value> {};

	/***************************************************************************************/
	template <class T, class U=void>
	using enable_if_vctr_cpu_r_3d = typename std::enable_if<is_vctr_cpu_r_3d<T>::value, U>::type;

	template <class T, class U, class V=void>
	using enable_if_vctr_cpu_r_3d_and_vctr_cpu = typename std::enable_if<is_vctr_cpu_r_3d_and_vctr_cpu<T, U>::value, V>::type;
		
	/***************************************************************************************/
	template <class T>
	struct is_vctr_cpu_r_nd: std::integral_constant<dt_bool, is_vctr_cpu<T>::value && is_r_nd<typename T::value_type>::value> {};

	template <class T, class U>
	struct is_vctr_cpu_r_nd_and_vctr_cpu: std::integral_constant<dt_bool, is_vctr_cpu_r_nd<T>::value && is_vctr_cpu<U>::value> {};

	/***************************************************************************************/
	template <class T, class U=void>
	using enable_if_vctr_cpu_r_nd = typename std::enable_if<is_vctr_cpu_r_nd<T>::value, U>::type;

	template <class T, class U, class V=void>
	using enable_if_vctr_cpu_r_nd_and_vctr_cpu = typename std::enable_if<is_vctr_cpu_r_nd_and_vctr_cpu<T, U>::value, V>::type;
}