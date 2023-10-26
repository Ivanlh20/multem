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

#include "poly_coef_1d.h"	

// class linear polynomial coefficients
namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDev Dev>
	Poly_Coef_1d<T, Dev>::Poly_Coef_1d() {}

	template <class T, eDev Dev>
	Poly_Coef_1d<T, Dev>::Poly_Coef_1d(const Vctr_cpu<T>& c0, const Vctr_cpu<T>& c1): 
		c0(c0), c1(c1){}

	template <class T, eDev Dev>
	Poly_Coef_1d<T, Dev>::Poly_Coef_1d(const dt_init_list_f64& c0, const dt_init_list_f64& c1): 
		c0(c0), c1(c1) {}

	template <class T, eDev Dev>
	Poly_Coef_1d<T, Dev>::Poly_Coef_1d(const size_type& new_size)
	{
		resize(new_size);
	}

	template <class T, eDev Dev>
	Poly_Coef_1d<T, Dev>::Poly_Coef_1d(const size_type& new_size, const T& value)
	{
		resize(new_size, value);
	}

	/* copy constructor */
	template <class T, eDev Dev>
	Poly_Coef_1d<T, Dev>::Poly_Coef_1d(const Poly_Coef_1d<T, Dev>& coef_poly1)
	{
		*this = coef_poly1;
	}

	/* converting constructor */
	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	Poly_Coef_1d<T, Dev>::Poly_Coef_1d(const Poly_Coef_1d<U, Dev_u>& coef_poly1)
	{
		*this = coef_poly1;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDev Dev>
	Poly_Coef_1d<T, Dev>& Poly_Coef_1d<T, Dev>::operator=(const Poly_Coef_1d<T, Dev>& coef_poly1)
	{
		this->assign(coef_poly1);
	
		return *this;
	}

	/* converting assignment operator */
	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	Poly_Coef_1d<T, Dev>& Poly_Coef_1d<T, Dev>::operator=(const Poly_Coef_1d<U, Dev_u>& coef_poly1)
	{
		this->assign(coef_poly1);
	
		return *this;
	}

	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	void Poly_Coef_1d<T, Dev>::assign(const Poly_Coef_1d<U, Dev_u>& coef_poly1)
	{ 
		c0.assign(coef_poly1.c0);
		c1.assign(coef_poly1.c1);
	}

	/**************** user define conversion operators *******************/
	template <class T, eDev Dev>
	pPoly_Coef_1d<T, Dev> Poly_Coef_1d<T, Dev>::ptr() const
	{
		return pPoly_Coef_1d<T, Dev>(*this);
	}

	/* user define conversion for pointer Vctr */
	template <class T, eDev Dev>
	Poly_Coef_1d<T, Dev>::operator pPoly_Coef_1d<T, Dev>() const
	{
		return pPoly_Coef_1d<T, Dev>(*this);
	}

	template <class T, eDev Dev>
	void Poly_Coef_1d<T, Dev>::fill(const T& val_c0, const T& val_c1)
	{
		c0.fill(val_c0);
		c1.fill(val_c1);
	}

	template <class T, eDev Dev>
	typename Poly_Coef_1d<T, Dev>::size_type Poly_Coef_1d<T, Dev>::size() const
	{
		return c0.size();
	}

	template <class T, eDev Dev>
	void Poly_Coef_1d<T, Dev>::clear()
	{
		c0.clear();
		c1.clear();
	}

	template <class T, eDev Dev>
	void Poly_Coef_1d<T, Dev>::reserve(const size_type& new_size)
	{
		c0.reserve(new_size);
		c1.reserve(new_size);
	}

	template <class T, eDev Dev>
	void Poly_Coef_1d<T, Dev>::resize(const size_type& new_size)
	{
		c0.resize(new_size);
		c1.resize(new_size);
	}

	template <class T, eDev Dev>
	void Poly_Coef_1d<T, Dev>::resize(const size_type& new_size, const T& value)
	{
		c0.resize(new_size, value);
		c1.resize(new_size, value);
	}

	template <class T, eDev Dev>
	void Poly_Coef_1d<T, Dev>::shrink_to_fit()
	{
		c0.shrink_to_fit();
		c1.shrink_to_fit();
	}
}

// pointer class linear polynomial coefficients
namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDev Dev>
	CGPU_EXEC
	pPoly_Coef_1d<T, Dev>::pPoly_Coef_1d(): c0(nullptr), c1(nullptr), m_size(0) {}

	template <class T, eDev Dev>
	CGPU_EXEC
	pPoly_Coef_1d<T, Dev>::pPoly_Coef_1d(T* c0, T* c1, const size_type& size): c0(c0), c1(c1), m_size(size) {}

	/* copy constructor */
	template <class T, eDev Dev>
	CGPU_EXEC
	pPoly_Coef_1d<T, Dev>::pPoly_Coef_1d(const pPoly_Coef_1d<T, Dev>& pcoef_poly1)
	{
		*this = pcoef_poly1;
	}

	/* constructor from Poly_Coef_1d */
	template <class T, eDev Dev>
	pPoly_Coef_1d<T, Dev>::pPoly_Coef_1d(const Poly_Coef_1d<T, Dev>& coef_poly1)
	{
		*this = coef_poly1;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDev Dev>
	CGPU_EXEC
	pPoly_Coef_1d<T, Dev>& pPoly_Coef_1d<T, Dev>::operator=(const pPoly_Coef_1d<T, Dev>& pcoef_poly1)
	{
		if (this != &pcoef_poly1)
		{
			c0 = pcoef_poly1.c0;
			c1 = pcoef_poly1.c1;
			m_size = pcoef_poly1.m_size;
		}

		return *this;
	}

	/* Assignment operator: Poly_Coef_1d -> pPoly_Coef_1d */
	template <class T, eDev Dev>
	CPU_EXEC
	pPoly_Coef_1d<T, Dev>& pPoly_Coef_1d<T, Dev>::operator=(const Poly_Coef_1d<T, Dev>& coef_poly1)
	{
		c0 = coef_poly1.c0.data();
		c1 = coef_poly1.c1.data();
		m_size = size_type(coef_poly1.size());

		return *this;
	}

	template <class T, eDev Dev>
	typename pPoly_Coef_1d<T, Dev>::size_type pPoly_Coef_1d<T, Dev>::size() const
	{
		return m_size;
	}
}