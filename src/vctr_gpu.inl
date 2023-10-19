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

#include "vctr_gpu.h"

/* gpu vector */
namespace mt
{
	/************************************* constructors ************************************/
	template <class T>
	Vctr<T, edev_gpu>::Vctr(): m_data(nullptr), m_s0(0), m_s1(1), m_s2(1), m_s3(1), 
	m_size(shape_size()), m_capacity(0)
	{
		set_picth();
	}

	template <class T>
	Vctr<T, edev_gpu>::Vctr(const dt_init_list_f64& data): Vctr()
	{
		assign(data.begin(), data.end());
	}

	template <class T>
	Vctr<T, edev_gpu>::Vctr(size_type s0): Vctr()
	{
		resize({ s0, size_type(1), size_type(1), size_type(1) });
	}

	template <class T>
	Vctr<T, edev_gpu>::Vctr(size_type s0, const T& value): Vctr()
	{
		resize({ s0, size_type(1), size_type(1), size_type(1)}, value);
	}

	template <class T>
	Vctr<T, edev_gpu>::Vctr(const dt_shape_st<size_type>& shape): Vctr()
	{
		resize(shape);
	}

	template <class T>
	Vctr<T, edev_gpu>::Vctr(const dt_shape_st<size_type>& shape, const T& value): Vctr()
	{
		resize(shape, value);
	}

	/* copy constructor */
	template <class T>
	Vctr<T, edev_gpu>::Vctr(const Vctr<T, edev_gpu>& vctr): Vctr()
	{
		*this = vctr;
	}

	/* Move constructor */
	template <class T>
	Vctr<T, edev_gpu>::Vctr(Vctr<T, edev_gpu>&& vctr): Vctr()
	{
		*this = std::move(vctr);
	}

	/* converting constructor */
	template <class T>
	template <class U>
	Vctr<T, edev_gpu>::Vctr(const Vctr<U, edev_gpu>& vctr): Vctr()
	{
		assign(vctr);
	}

	template <class T>
	template <class U>
	Vctr<T, edev_gpu>::Vctr(U* first, U* last): Vctr()
	{
		assign(first, last);
	}

	template <class T>
	template <class U, class V, class>
	Vctr<T, edev_gpu>::Vctr(U *p, dt_int64 n_p, size_type icol): Vctr()
	{
		set_pos_from_cpu_ptr(p, n_p, icol);
	}

	template <class T>
	template <class U>
	Vctr<T, edev_gpu>::Vctr(const std::vector<U>& vctr): Vctr()
	{
		assign(vctr);
	}

	// from gpu pVctr to Vctr
	template <class T>
	template <class U, class STU>
	Vctr<T, edev_gpu>::Vctr(const pVctr<U, edev_gpu, STU>& pvctr): Vctr()
	{
		assign(pvctr);
	}

	// from cpu pVctr to Vctr
	template <class T>
	template <class U, class STU>
	Vctr<T, edev_gpu>::Vctr(const pVctr<U, edev_cpu, STU>& pvctr): Vctr()
	{
		assign(pvctr);
	}

	template <class T>
	template <class U>
	Vctr<T, edev_gpu>::Vctr(const Vctr<U, edev_cpu>& vctr): Vctr()
	{
		assign(vctr);
	}

	template <class T>
	template <class U>
	Vctr<T, edev_gpu>::Vctr(const thrust::host_vector<U>& vctr): Vctr()
	{
		assign(vctr);
	}

	template <class T>
	template <class U>
	Vctr<T, edev_gpu>::Vctr(const thrust::device_vector<U>& vctr): Vctr()
	{
		assign(vctr);
	}

	template <class T>
	Vctr<T, edev_gpu>::~Vctr()
	{
		destroy();
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T>
	Vctr<T, edev_gpu>& Vctr<T, edev_gpu>::operator=(const Vctr<T, edev_gpu>& vctr)
	{
		if (vctr.m_data == nullptr)
		{
			destroy();
		}
		else if (this != &vctr)
		{
			fcn_cuda_free(m_data);

			fcn_cuda_malloc(m_data, vctr.m_capacity);

			memcpy_gpu_gpu(m_data, vctr.data(), vctr.size());

			m_s0 = vctr.m_s0;
			m_s1 = vctr.m_s1;
			m_s2 = vctr.m_s2;
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
	Vctr<T, edev_gpu>& Vctr<T, edev_gpu>::operator=(Vctr<T, edev_gpu>&& vctr)
	{
		if (vctr.m_data == nullptr)
		{
			destroy();
		}
		else if (this != &vctr)
		{
			fcn_cuda_free(m_data);

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
	Vctr<T, edev_gpu>& Vctr<T, edev_gpu>::operator=(const Vctr<U, edev_gpu>& vctr)
	{
		assign(vctr);
			
		return *this;
	}

	template <class T>
	template <class U>
	Vctr<T, edev_gpu>& Vctr<T, edev_gpu>::operator=(const std::vector<U>& vctr)
	{
		resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
		memcpy_cpu_gpu(m_data, vctr.data(), vctr.size());

		return *this;
	}

	template <class T>
	template <class U>
	Vctr<T, edev_gpu>& Vctr<T, edev_gpu>::operator=(const Vctr<U, edev_cpu>& vctr)
	{
		assign(vctr);
			
		return *this;
	}

	template <class T>
	template <class U>
	Vctr<T, edev_gpu>& Vctr<T, edev_gpu>::operator=(const thrust::host_vector<U>& vctr)
	{
		assign(vctr);
			
		return *this;
	}

	template <class T>
	template <class U>
	Vctr<T, edev_gpu>& Vctr<T, edev_gpu>::operator=(const thrust::device_vector<U>& vctr)
	{
		assign(vctr);
			
		return *this;
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_gpu>::assign(const Vctr<U, edev_gpu>& vctr, T* pvctr_cpu)
	{
		if ((void*)this != (void*)&vctr)
		{
			this->allocate(vctr.shape());
			vctr_cpy(*this, vctr, pvctr_cpu);
		}
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_gpu>::assign(U* first, U* last)
	{
		if ((void*)m_data != (void*)first) 
		{
			auto m_size_i = std::distance(first, last);
			if (m_size_i>0)
			{
				auto m_size = size_type(m_size_i);
				resize({ m_size, size_type(1), size_type(1), size_type(1) });
				memcpy_cpu_gpu(m_data, first, m_size);
			}
		}
	}
			
	// from gpu pVctr to Vctr
	template <class T>
	template <class U, class STU>
	void Vctr<T, edev_gpu>::assign(const pVctr<U, edev_gpu, STU>& pvctr)
	{
		resize(pvctr.shape());
		memcpy_gpu_gpu(m_data, pvctr.data(), m_size);
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_gpu>::assign(const thrust::device_ptr<U>& first, const thrust::device_ptr<U>& last, T* pvctr_cpu)
	{
		auto m_size_i = thrust::distance(first, last);
		if (m_size_i>0)
		{
			auto m_size = size_type(m_size_i);
			resize({ m_size, size_type(1), size_type(1), size_type(1) });
			memcpy_gpu_gpu(m_data, (U*)first.get(), m_size, pvctr_cpu);
		}
	}

	// from cpu pVctr to Vctr
	template <class T>
	template <class U, class STU>
	void Vctr<T, edev_gpu>::assign(const pVctr<U, edev_cpu, STU>& pvctr)
	{
		resize(pvctr.shape());
		memcpy_cpu_gpu(m_data, pvctr.data(), m_size);
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_gpu>::assign(const std::vector<U>& vctr, T* pvctr_cpu)
	{
		resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
		memcpy_cpu_gpu(m_data, vctr.data(), vctr.size(), pvctr_cpu);
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_gpu>::assign(const Vctr<U, edev_cpu>& vctr, T* pvctr_cpu)
	{
		this->allocate(vctr.shape());
		vctr_cpy(*this, vctr, pvctr_cpu);
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_gpu>::assign(const thrust::host_vector<U>& vctr, T* pvctr_cpu)
	{
		resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
		memcpy_cpu_gpu(m_data, vctr.data(), vctr.size(), pvctr_cpu);
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_gpu>::assign(const thrust::device_vector<U>& vctr, T* pvctr_cpu)
	{
		resize({ size_type(vctr.size()), size_type(1), size_type(1), size_type(1) });
		memcpy_gpu_gpu(m_data, (T*)vctr.data().get(), vctr.size(), pvctr_cpu);
	}

	/**************** user define conversion operators *******************/
	template <class T>
	pVctr_gpu_32<T> Vctr<T, edev_gpu>::ptr_32() const
	{
		
		return pVctr_gpu_32<T>(*this);
	}

	template <class T>
	pVctr_gpu_64<T> Vctr<T, edev_gpu>::ptr_64() const
	{
		
		return pVctr_gpu_64<T>(*this);
	}

	/* user define conversion for pointer Vctr */
	template <class T>
	Vctr<T, edev_gpu>::operator pVctr_gpu_32<T>() const
	{
		return pVctr_gpu_32<T>(*this);
	}

	template <class T>
	Vctr<T, edev_gpu>::operator pVctr_gpu_64<T>() const
	{
		return pVctr_gpu_64<T>(*this);
	}

	/* user define conversion */
	template <class T>
	Vctr<T, edev_gpu>::operator std::vector<T>() const
	{
		std::vector<T> vctr(m_size);
		memcpy_gpu_cpu(vctr.data(), m_data, m_size);

		return vctr;
	}

	/* user define conversion in which T is the complemented precision */
	template <class T>
	template <class U>
	Vctr<T, edev_gpu>::operator std::vector<U>() const
	{
		std::vector<U> vctr(m_size);
		memcpy_gpu_cpu(vctr.data(), m_data, m_size);

		return vctr;
	}

	/* user define conversion to output type std::vector<thrust::complex<dt_float32>> */
	template <class T>
	template <class U, class>
	Vctr<T, edev_gpu>::operator std::vector<thrust::complex<dt_float32>>() const
	{
		std::vector<thrust::complex<dt_float32>> vctr(m_size);
		memcpy_gpu_cpu(vctr.data(), m_data, m_size);
		return vctr;
	}

	/* user define conversion to output type std::vector<thrust::complex<dt_float64>> */
	template <class T>
	template <class U, class>
	Vctr<T, edev_gpu>::operator std::vector<thrust::complex<dt_float64>>() const
	{
		std::vector<thrust::complex<dt_float64>>vctr(m_size);
		memcpy_gpu_cpu(vctr.data(), m_data, m_size);
		return vctr;
	}

	/* user define conversion */
	template <class T>
	Vctr<T, edev_gpu>::operator thrust::host_vector<T>() const
	{
		thrust::host_vector<T> vctr(m_size);
		memcpy_gpu_cpu(vctr.data(), m_data, m_size);

		return vctr;
	}

	/* user define conversion in which T is the complemented precision */
	template <class T>
	template <class U>
	Vctr<T, edev_gpu>::operator thrust::host_vector<U>() const
	{
		thrust::host_vector<U> vctr(m_size);
		memcpy_gpu_cpu(vctr.data(), m_data, m_size);

		return vctr;
	}

	/* user define conversion to output type thrust::host_vector<std::complex<dt_float32>> */
	template <class T>
	template <class U, class>
	Vctr<T, edev_gpu>::operator thrust::host_vector<std::complex<dt_float32>>() const
	{
		thrust::host_vector<std::complex<dt_float32>> vctr(m_size);
		memcpy_gpu_cpu(vctr.data(), m_data, m_size);

		return vctr;
	}

	/* user define conversion to output type thrust::host_vector<std::complex<dt_float64>> */
	template <class T>
	template <class U, class>
	Vctr<T, edev_gpu>::operator thrust::host_vector<std::complex<dt_float64>>() const
	{
		thrust::host_vector<std::complex<dt_float64>> vctr(m_size);
		memcpy_gpu_cpu(vctr.data(), m_data, m_size);

		return vctr;
	}

	/***************************************************************************************/
	/* user define conversion */
	template <class T>
	Vctr<T, edev_gpu>::operator thrust::device_vector<T>() const
	{
		thrust::host_vector<T> vctr(m_size);
		memcpy_gpu_gpu((T*)vctr.data().get(), m_data, m_size);
		return vctr;
	}

	/* user define conversion in which T is the complemented precision */
	template <class T>
	template <class U>
	Vctr<T, edev_gpu>::operator thrust::device_vector<U>() const
	{
		thrust::host_vector<U> vctr(m_size);
		memcpy_gpu_gpu((U*)vctr.data().get(), m_data, m_size);

		return vctr;
	}

	/***************************************************************************************/
	template <class T>
	template <class U>
	void Vctr<T, edev_gpu>::cpy_to_cpu_ptr(U* pdata, size_type n_data, T* pvctr_cpu)
	{
		n_data = std::min(m_size, n_data);
		memcpy_gpu_cpu(pdata, m_data, n_data, pvctr_cpu);
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_gpu>::cpy_to_cpu_ptr(U* first, U* last, T* pvctr_cpu)
	{
		auto n_data = std::distance(first, last);
		memcpy_gpu_cpu(first, m_data, n_data, pvctr_cpu);
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_gpu>::cpy_to_gpu_ptr(U* pdata, size_type n_data, U* pvctr_cpu)
	{
		n_data = std::min(m_size, n_data);
		memcpy_gpu_gpu(pdata, m_data, n_data, pvctr_cpu);
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_gpu>::cpy_to_gpu_ptr(U* first, U* last, U* pvctr_cpu)
	{
		auto n_data = thrust::distance(first, last);
		memcpy_gpu_gpu(first, m_data, n_data, pvctr_cpu);
	}

	/***************************************************************************************/
	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_gpu>::cpy_real_to_cpu_ptr(U* pdata, size_type n_data, T* pvctr_cpu)
	{
		n_data = std::min(m_size, n_data);
		memcpy_real_gpu_cpu(pdata, m_data, n_data, pvctr_cpu);
	}

	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_gpu>::cpy_real_to_cpu_ptr(U* first, U* last, T* pvctr_cpu)
	{
		auto n_data = std::distance(first, last);
		memcpy_real_gpu_cpu(first, m_data, n_data, pvctr_cpu);
	}

	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_gpu>::cpy_real_to_gpu_ptr(U* pdata, size_type n_data, U* pvctr_cpu)
	{
		n_data = std::min(m_size, n_data);
		memcpy_real_gpu_gpu(pdata, m_data, n_data, pvctr_cpu);
	}

	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_gpu>::cpy_real_to_gpu_ptr(U* first, U* last, U* pvctr_cpu)
	{
		auto n_data = std::distance(first, last);
		memcpy_real_gpu_gpu(first, m_data, n_data, pvctr_cpu);
	}

	/***************************************************************************************/
	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_gpu>::cpy_imag_to_cpu_ptr(U* pdata, size_type n_data, T* pvctr_cpu)
	{
		n_data = std::min(m_size, n_data);
		memcpy_imag_gpu_cpu(pdata, m_data, n_data, pvctr_cpu);
	}

	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_gpu>::cpy_imag_to_cpu_ptr(U* first, U* last, T* pvctr_cpu)
	{
		auto n_data = std::distance(first, last);
		memcpy_imag_gpu_cpu(first, m_data, n_data, pvctr_cpu);
	}

	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_gpu>::cpy_imag_to_gpu_ptr(U* pdata, size_type n_data, U* pvctr_cpu)
	{
		n_data = std::min(m_size, n_data);
		memcpy_imag_gpu_gpu(pdata, m_data, n_data, pvctr_cpu);
	}

	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_gpu>::cpy_imag_to_gpu_ptr(U* first, U* last, U* pvctr_cpu)
	{
		auto n_data = std::distance(first, last);
		memcpy_imag_gpu_gpu(first, m_data, n_data, pvctr_cpu);
	}

	/***************************************************************************************/
	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_gpu>::set_pos_from_cpu_ptr(U *p, size_type n_p, dt_int64 icol)
	{
		resize({n_p});

		memcpy_pos_cpu_gpu(m_data, p, size(), icol);
	}

	template <class T>
	template <class U, class V, class>
	void Vctr<T, edev_gpu>::cpy_pos_to_cpu_ptr(U *p, size_type icol)
	{
		memcpy_pos_gpu_cpu(p, m_data, size(), icol);
	}

	/***************************************************************************************/
	template <class T>
	void Vctr<T, edev_gpu>::resize(const dt_shape_st<size_type>& shape)
	{
		this->allocate(shape);
	}

	template <class T>
	void Vctr<T, edev_gpu>::resize(const dt_shape_st<size_type>& shape, const T& value)
	{
		auto m_size_t = m_size;
		this->allocate(shape);
		thrust::fill(this->begin()+m_size_t, this->end(), value);
	}

	template <class T>
	void Vctr<T, edev_gpu>::reserve(const dt_shape_st<size_type>& shape)
	{
		this->allocate(shape, true);
	}

	template <class T>
	void Vctr<T, edev_gpu>::shrink_to_fit()
	{
		if (m_size < m_capacity)
		{
			if (m_size > 0)
			{
				/* allocate memory and transfer data */
				T* data_tmp = nullptr;
				fcn_cuda_malloc(data_tmp, m_size);
				memcpy_gpu_gpu(data_tmp, m_data, m_size);
				fcn_cuda_free(m_data);

				m_data = data_tmp;
				data_tmp = nullptr;
			}
			else
			{
				fcn_cuda_free(m_data);
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
	void Vctr<T, edev_gpu>::push_back(const U& val)
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

		memcpy_cpu_gpu(m_data + m_size, &val, 1);

		m_size++;
	}

	template <class T>
	template <class U>
	void Vctr<T, edev_gpu>::push_back(const Vctr<U, edev_cpu>& vctr)
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

		memcpy_cpu_gpu(m_data + m_size, vctr.data(), vctr.size());

		m_size = size_new;
	}

	template <class T>
	void Vctr<T, edev_gpu>::pop_back()
	{
		if (m_size >= 1)
		{
			m_size--;
		}
	}

	template <class T>
	void Vctr<T, edev_gpu>::fill(T val)
	{
		thrust::fill(this->begin(), this->end(), val);
	}

	/***************************************************************************************/
	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::s0() const
	{
		return m_s0;
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::s1() const
	{
		return m_s1;
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::s2() const
	{
		return m_s2;
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::s3() const
	{
		return m_s2;
	}

	template <class T>
	dt_int32 Vctr<T, edev_gpu>::s0_32() const
	{
		return static_cast<dt_int32>(m_s0);
	}

	template <class T>
	dt_int32 Vctr<T, edev_gpu>::s1_32() const
	{
		return static_cast<dt_int32>(m_s1);
	}

	template <class T>
	dt_int32 Vctr<T, edev_gpu>::s2_32() const
	{
		return static_cast<dt_int32>(m_s2);
	}

	template <class T>
	dt_int32 Vctr<T, edev_gpu>::s3_32() const
	{
		return static_cast<dt_int32>(m_s3);
	}
			
	template <class T>
	dt_int64 Vctr<T, edev_gpu>::s0_64() const
	{
		return static_cast<dt_int64>(m_s0);
	}

	template <class T>
	dt_int64 Vctr<T, edev_gpu>::s1_64() const
	{
		return static_cast<dt_int64>(m_s1);
	}

	template <class T>
	dt_int64 Vctr<T, edev_gpu>::s2_64() const
	{
		return static_cast<dt_int64>(m_s2);
	}

	template <class T>
	dt_int64 Vctr<T, edev_gpu>::s3_64() const
	{
		return static_cast<dt_int64>(m_s3);
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::s0h() const
	{
		return m_s0/size_type(2);
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::s1h() const
	{
		return m_s1/size_type(2);
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::s2h() const
	{
		return m_s2/size_type(2);
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::s3h() const
	{
		return m_s3/size_type(2);
	}

	template <class T>
	dt_shape_st<typename Vctr<T, edev_gpu>::size_type> Vctr<T, edev_gpu>::shape() const
	{
		return {m_s0, m_s1, m_s2, m_s3};
	}

	template <class T>
	dt_shape_st<typename Vctr<T, edev_gpu>::size_type> Vctr<T, edev_gpu>::shape_2d_trs() const
	{
		return {m_s1, m_s0, m_s2, m_s3};
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::shape_size() const 
	{
		return m_s0*max(m_s1, size_type(1))*max(m_s2, size_type(1))*max(m_s3, size_type(1));
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::pitch_s1() const
	{
		return m_pitch_s1;
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::pitch_s2() const
	{
		return m_pitch_s2;
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::pitch_s3() const
	{
		return m_pitch_s3;
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::size() const
	{
		return m_size;
	}

	template <class T>
	dt_int32 Vctr<T, edev_gpu>::size_32() const
	{
		return static_cast<dt_int32>(m_size);
	}

	template <class T>
	dt_int64 Vctr<T, edev_gpu>::size_64() const
	{
		return static_cast<dt_int64>(m_size);
	}

	template <class T>
	iGrid_1d Vctr<T, edev_gpu>::igrid_1d() const
	{
		return {size_32()};
	}

	template <class T>
	iGrid_2d Vctr<T, edev_gpu>::igrid_2d() const
	{
		return {s1_32(), s0_32()};
	}

	template <class T>
	iGrid_3d Vctr<T, edev_gpu>::igrid_3d() const
	{
		return {s1_32(), s0_32(), s2_32()};
	}

	template <class T>
	iGrid_1d_64 Vctr<T, edev_gpu>::igrid_1d_64() const
	{
		return {size_64()};
	}

	template <class T>
	iGrid_2d_64 Vctr<T, edev_gpu>::igrid_2d_64() const
	{
		return {s1_64(), s0_64()};
	}

	template <class T>
	iGrid_3d_64 Vctr<T, edev_gpu>::igrid_3d_64() const
	{
		return {s1_64(), s0_64(), s2_64()};
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::capacity() const
	{
		return m_capacity;
	}

	template <class T>
	dt_bool Vctr<T, edev_gpu>::empty() const
	{
		return m_size == 0;
	}

	template <class T>
	dt_bool Vctr<T, edev_gpu>::is_1d() const
	{
		return (m_s0 == 1) || (m_s1 == 1);
	}

	template <class T>
	void Vctr<T, edev_gpu>::clear()
	{
		m_size = 0;
	}

	template <class T>
	void Vctr<T, edev_gpu>::clear_shrink_to_fit()
	{
		destroy();
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::sub_2_ind(const size_type& ix_0) const 
	{ 
		return ix_0;
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::sub_2_ind(const size_type& ix_0, const size_type& ix_1) const 
	{ 
		return ix_0 + m_s0*ix_1;
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::sub_2_ind(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2) const 
	{ 
		return ix_0 + m_s0*(ix_1 + m_s1*ix_2);
	}

	template <class T>
	typename Vctr<T, edev_gpu>::size_type Vctr<T, edev_gpu>::sub_2_ind(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2, const size_type& ix_3) const 
	{ 
		return ix_0 + m_s0*(ix_1 + m_s1*(ix_2 + m_s2*ix_3));
	}

	//template <class T>
	// T& Vctr<T, edev_gpu>::operator[](const size_type& iy)
	// { 
	// 	return m_data[iy];
	// }

	//template <class T>
	// const T& Vctr<T, edev_gpu>::operator[](const size_type& iy) const 
	// { 
	// 	return m_data[iy];
	// }

	//template <class T>
	// T& Vctr<T, edev_gpu>::operator()(const size_type& iy)
	// { 
	// 	return m_data[iy];
	// }

	//template <class T>
	// const T& Vctr<T, edev_gpu>::operator()(const size_type& iy) const 
	// { 
	// 	return m_data[iy];
	// }

	//template <class T>
	// T& Vctr<T, edev_gpu>::operator()(const size_type& ix_0, const size_type& ix_1)
	// { 
	// 	return m_data[sub_2_ind(ix_0, ix_1)];
	// }

	//template <class T>
	// const T& Vctr<T, edev_gpu>::operator()(const size_type& ix_0, const size_type& ix_1) const 
	// { 
	// 	return m_data[sub_2_ind(ix_0, ix_1)];
	// }

	//template <class T>
	// T& Vctr<T, edev_gpu>::operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2)
	// { 
	// 	return m_data[sub_2_ind(ix_0, ix_1, ix_2)];
	// }

	//template <class T>
	// const T& Vctr<T, edev_gpu>::operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2) const 
	// { 
	// 	return m_data[sub_2_ind(ix_0, ix_1, ix_2)];
	// }

	// T& Vctr<T, edev_gpu>::operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2, const size_type& ix_3)
	// { 
	// 	return m_data[sub_2_ind(ix_0, ix_1, ix_2, ix_3)];
	// }

	//template <class T>
	// const T& Vctr<T, edev_gpu>::operator()(const size_type& ix_0, const size_type& ix_1, const size_type& ix_2, const size_type& ix_3) const 
	// { 
	// 	return m_data[sub_2_ind(ix_0, ix_1, ix_2, ix_3)];
	// }

	template <class T>
	thrust::device_ptr<T> Vctr<T, edev_gpu>::begin() noexcept 
	{ 
		return thrust::device_ptr<T>(m_data);
	}

	template <class T>
	const thrust::device_ptr<T> Vctr<T, edev_gpu>::begin() const
	{ 
		return thrust::device_ptr<T>(m_data);
	}

	template <class T>
	thrust::device_ptr<T> Vctr<T, edev_gpu>::end() noexcept 
	{ 
		return thrust::device_ptr<T>(m_data) + m_size;
	}

	template <class T>
	const thrust::device_ptr<T> Vctr<T, edev_gpu>::end() const
	{ 
		return thrust::device_ptr<T>(m_data) + m_size;
	}

	template <class T>
	T* Vctr<T, edev_gpu>::data()
	{
		return m_data;
	}

	template <class T>
	const T* Vctr<T, edev_gpu>::data() const
	{
		return m_data;
	}

	template <class T>
	template <class U>
	U Vctr<T, edev_gpu>::data_cast()
	{
		return reinterpret_cast<U>(m_data);
	}

	template <class T>
	template <class U>
	const U Vctr<T, edev_gpu>::data_cast() const
	{
		return reinterpret_cast<U>(m_data);
	}

	//template <class T>
	// T& Vctr<T, edev_gpu>::front()
	// { 
	// 	return m_data[0];
	// }

	//template <class T>
	// const T& Vctr<T, edev_gpu>::front() const 
	// { 
	// 	return m_data[0];
	// }

	//template <class T>
	// T& Vctr<T, edev_gpu>::back(){ 
	// 	return m_data[m_size-1];
	// }

	//template <class T>
	// const T& Vctr<T, edev_gpu>::back() const 
	// { 
	// 	return m_data[m_size-1];
	// }

	// set shape
	template <class T>
	void Vctr<T, edev_gpu>::set_shape(const dt_shape_st<size_type>& shape)
	{
		m_s0 = max(shape[0], size_type(0));
		m_s1 = max(shape[1], size_type(1));
		m_s2 = max(shape[2], size_type(1));
		m_s3 = max(shape[3], size_type(1));

		set_picth();
	}

	template <class T>
	void Vctr<T, edev_gpu>::trs_shape_2d()
	{
		set_shape({m_s1, m_s0, m_s2, m_s3});
	}

	template <class T>
	template <class U, eDev Dev>
	void Vctr<T, edev_gpu>::set_shape(const Vctr<U, Dev>& vctr, dt_bool bb_size)
	{
		this->set_shape(vctr.shape(), bb_size);
	}

	FCNS_IMP_GPU_GRID_BLK_VCTR(template <class T>, Vctr<T COMMA edev_gpu>);

	template <class T>
	void Vctr<T, edev_gpu>::set_picth()
	{
		m_pitch_s1 = m_s0*sizeof(T);
		m_pitch_s2 = m_pitch_s1*m_s1;
		m_pitch_s3 = m_pitch_s2*m_s2;
	}

	template <class T>
	void Vctr<T, edev_gpu>::set_capacity(size_type size_r)
	{
		m_capacity = max(m_capacity, size_r);
	}

	template <class T>
	void Vctr<T, edev_gpu>::set_shape_cstr(dt_shape_st<size_type>& shape)
	{
		shape[0] = max(shape[0], size_type(0));
		shape[1] = max(shape[1], size_type(1));
		shape[2] = max(shape[2], size_type(1));
		shape[3] = max(shape[3], size_type(1));
	}

	// reallocate and copy memory
	template <class T>
	void Vctr<T, edev_gpu>::allocate(dt_shape_st<size_type> shape, dt_bool bb_reserve)
	{
		set_shape_cstr(shape);

		auto size_r = shape[0]*shape[1]*shape[2]*shape[3];
			
		if (size_r == 0)
		{
			return;
		}

		if (m_size == 0)
		{
			fcn_cuda_free(m_data);

			fcn_cuda_malloc(m_data, size_r);
		}
		else if (size_r>m_capacity)
		{
			// allocate	memory and transfer data
			T* data_tmp = nullptr;
			fcn_cuda_malloc(data_tmp, size_r);
			memcpy_gpu_gpu(data_tmp, m_data, m_size);
			fcn_cuda_free(m_data);

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

	// initialization
	template <class T>
	void Vctr<T, edev_gpu>::init()
	{
		m_data = nullptr;
		m_s0 = 0; m_s1 = m_s2 = m_s3 = 1;
		m_pitch_s1 = m_pitch_s2 = m_pitch_s3 = 0;
		m_size = 0;
		m_capacity = 0;
	}

	// destroy memory on the device
	template <class T>
	void Vctr<T, edev_gpu>::destroy()
	{
		if (m_data != nullptr)
		{
			fcn_cuda_free(m_data);
		}

		init();
	}
#endif
}

/* derived class */
namespace mt
{
#ifdef __CUDACC__
	template <class T>
	struct is_vctr_gpu_r_2d: std::integral_constant<dt_bool, is_vctr_gpu<T>::value && is_r_2d<typename T::value_type>::value> {};

	template <class T, class U>
	struct is_vctr_gpu_r_2d_and_vctr_gpu: std::integral_constant<dt_bool, is_vctr_gpu_r_2d<T>::value && is_vctr_gpu<U>::value> {};

	/***************************************************************************************/
	template <class T, class U=void>
	using enable_if_vctr_gpu_r_2d = typename std::enable_if<is_vctr_gpu_r_2d<T>::value, U>::type;

	template <class T, class U, class V=void>
	using enable_if_vctr_gpu_r_2d_and_vctr_gpu = typename std::enable_if<is_vctr_gpu_r_2d_and_vctr_gpu<T, U>::value, V>::type;

	/***************************************************************************************/
	template <class T>
	struct is_vctr_gpu_r_3d: std::integral_constant<dt_bool, is_vctr_gpu<T>::value && is_r_3d<typename T::value_type>::value> {};

	template <class T, class U>
	struct is_vctr_gpu_r_3d_and_vctr_gpu: std::integral_constant<dt_bool, is_vctr_gpu_r_3d<T>::value && is_vctr_gpu<U>::value> {};

	/***************************************************************************************/
	template <class T, class U=void>
	using enable_if_vctr_gpu_r_3d = typename std::enable_if<is_vctr_gpu_r_3d<T>::value, U>::type;

	template <class T, class U, class V=void>
	using enable_if_vctr_gpu_r_3d_and_vctr_gpu = typename std::enable_if<is_vctr_gpu_r_3d_and_vctr_gpu<T, U>::value, V>::type;

	/***************************************************************************************/
	template <class T>
	struct is_vctr_gpu_r_nd: std::integral_constant<dt_bool, is_vctr_gpu<T>::value && is_r_nd<typename T::value_type>::value> {};

	template <class T, class U>
	struct is_vctr_gpu_r_nd_and_vctr_gpu: std::integral_constant<dt_bool, is_vctr_gpu_r_nd<T>::value && is_vctr_gpu<U>::value> {};

	/***************************************************************************************/
	template <class T, class U=void>
	using enable_if_vctr_gpu_r_nd = typename std::enable_if<is_vctr_gpu_r_nd<T>::value, U>::type;

	template <class T, class U, class V=void>
	using enable_if_vctr_gpu_r_nd_and_vctr_cpu = typename std::enable_if<is_vctr_gpu_r_nd_and_vctr_gpu<T, U>::value, V>::type;
}