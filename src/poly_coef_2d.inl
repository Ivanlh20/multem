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

#include "intrpl_coef.h"	

/***************************************************************************************/
/************************* quadratic polynomial coefficients ***************************/
/***************************************************************************************/
namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDev Dev>
	CGPU_EXEC
	pPoly_Coef_2d<T, Dev>::pPoly_Coef_2d(): c0(nullptr), c1(nullptr), c2(nullptr), m_size(0) {}

	template <class T, eDev Dev>
	CGPU_EXEC
	pPoly_Coef_2d<T, Dev>::pPoly_Coef_2d(T* c0, T* c1, T* c2, const size_type& size): c0(c0), c1(c1), c2(c2), m_size(size) {}

	/* copy constructor */
	template <class T, eDev Dev>
	CGPU_EXEC
	pPoly_Coef_2d<T, Dev>::pPoly_Coef_2d(const pPoly_Coef_2d<T, Dev>& pcoef_poly2)
	{
		*this = pcoef_poly2;
	}

	/* constructor from Poly_Coef_2d */
	template <class T, eDev Dev>
	explicit pPoly_Coef_2d<T, Dev>::pPoly_Coef_2d(const Poly_Coef_2d<T, Dev>& coef_poly2)
	{
		*this = coef_poly2;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDev Dev>
	CGPU_EXEC
	pPoly_Coef_2d<T, Dev>& pPoly_Coef_2d<T, Dev>::operator=(const pPoly_Coef_2d<T, Dev>& pcoef_poly2)
	{
		if (this != &pcoef_poly2)
		{
			c0 = pcoef_poly2.c0;
			c1 = pcoef_poly2.c1;
			c2 = pcoef_poly2.c2;
			m_size = pcoef_poly2.m_size;
		}

		return *this;
	}

	/* Assignment operator: Poly_Coef_2d -> pPoly_Coef_2d */
	template <class T, eDev Dev>
	CPU_EXEC
	pPoly_Coef_2d<T, Dev>& pPoly_Coef_2d<T, Dev>::operator=(const Poly_Coef_2d<T, Dev>& coef_poly2)
	{
		c0 = coef_poly2.c0.data();
		c1 = coef_poly2.c1.data();
		c2 = coef_poly2.c2.data();
		m_size = size_type(coef_poly2.size());

		return *this;
	}

	template <class T, eDev Dev>
	size_type pPoly_Coef_2d<T, Dev>::size() const
	{
		return m_size;
	}
}

namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDev Dev>
	Poly_Coef_1d<T, Dev>::Poly_Coef_2d() {}

	template <class T, eDev Dev>
	Poly_Coef_1d<T, Dev>::Poly_Coef_2d(const Vctr_cpu<T>& c0, const Vctr_cpu<T>& c1, const Vctr_cpu<T>& c2): 
		c0(c0), c1(c1), c2(c2){}

	template <class T, eDev Dev>
	Poly_Coef_1d<T, Dev>::Poly_Coef_2d(const dt_init_list_f64& c0, const dt_init_list_f64& c1, const dt_init_list_f64& c2): 
		c0(c0), c1(c1), c2(c2) {}

	template <class T, eDev Dev>
	Poly_Coef_1d<T, Dev>::Poly_Coef_2d(const size_type& new_size)
	{
		resize(new_size);
	}

	template <class T, eDev Dev>
	Poly_Coef_1d<T, Dev>::Poly_Coef_2d(const size_type& new_size, const T& value)
	{
		resize(new_size, value);
	}

	/* copy constructor */
	template <class T, eDev Dev>
	Poly_Coef_1d<T, Dev>::Poly_Coef_2d(const Poly_Coef_2d<T, Dev>& coef_poly2)
	{
		*this = coef_poly2;
	}

	/* converting constructor */
	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	Poly_Coef_1d<T, Dev>::Poly_Coef_2d(const Poly_Coef_2d<U, Dev_u>& coef_poly2)
	{
		*this = coef_poly2;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDev Dev>
	Poly_Coef_2d<T, Dev>& Poly_Coef_1d<T, Dev>::operator=(const Poly_Coef_2d<T, Dev>& coef_poly2)
	{
		this->assign(coef_poly2);
	
		return *this;
	}

	/* converting assignment operator */
	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	Poly_Coef_2d<T, Dev>& Poly_Coef_1d<T, Dev>::operator=(const Poly_Coef_2d<U, Dev_u>& coef_poly2)
	{
		this->assign(coef_poly2);
	
		return *this;
	}

	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	voidPoly_Coef_1d<T, Dev>::assign(const Poly_Coef_2d<U, Dev_u>& coef_poly2)
	{ 
		c0.assign(coef_poly2.c0);
		c1.assign(coef_poly2.c1);
		c2.assign(coef_poly2.c2);
	}

	/**************** user define conversion operators *******************/
	template <class T, eDev Dev>
	pPoly_Coef_2d<T, Dev> Poly_Coef_1d<T, Dev>::ptr() const
	{
		return pPoly_Coef_2d<T, Dev>(*this);
	}

	/* user define conversion for pointer Vctr */
	template <class T, eDev Dev>
	Poly_Coef_1d<T, Dev>::operator pPoly_Coef_2d<T, Dev>() const
	{
		return pPoly_Coef_2d<T, Dev>(*this);
	}

	template <class T, eDev Dev>
	void Poly_Coef_1d<T, Dev>::fill(const T& val_c0, const T& val_c1, const T& val_c2)
	{
		c0.fill(val_c0);
		c1.fill(val_c1);
		c2.fill(val_c2);
	}

	template <class T, eDev Dev>
	size_type Poly_Coef_1d<T, Dev>::size() const
	{
		return c0.size();
	}

	template <class T, eDev Dev>
	void Poly_Coef_1d<T, Dev>::clear()
	{
		c0.clear();
		c1.clear();
		c2.clear();
	}

	template <class T, eDev Dev>
	void Poly_Coef_1d<T, Dev>::reserve(const size_type& new_size)
	{
		c0.reserve(new_size);
		c1.reserve(new_size);
		c2.reserve(new_size);
	}

	template <class T, eDev Dev>
	void Poly_Coef_1d<T, Dev>::resize(const size_type& new_size)
	{
		c0.resize(new_size);
		c1.resize(new_size);
		c2.resize(new_size);
	}

	template <class T, eDev Dev>
	void Poly_Coef_1d<T, Dev>::resize(const size_type& new_size, const T& value)
	{
		c0.resize(new_size, value);
		c1.resize(new_size, value);
		c2.resize(new_size, value);
	}

	template <class T, eDev Dev>
	void Poly_Coef_1d<T, Dev>::shrink_to_fit()
	{
		c0.shrink_to_fit();
		c1.shrink_to_fit();
		c2.shrink_to_fit();
	}
}

/***************************************************************************************/
/************************** cubic interpolation coefficients ***************************/
/***************************************************************************************/
namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDev Dev> 
	CGPU_EXEC
	pPoly_Coef_3d<T, Dev>::pPoly_Coef_3d(): c0(nullptr), c1(nullptr), c2(nullptr), c3(nullptr), m_size(0) {}

	template <class T, eDev Dev> 
	CGPU_EXEC
	pPoly_Coef_3d<T, Dev>::pPoly_Coef_3d(T* c0, T* c1, T* c2, T* c3, const size_type& size): c0(c0), c1(c1), c2(c2), c3(c3), m_size(size) {}

	/* copy constructor */
	template <class T, eDev Dev>
	CGPU_EXEC
	pPoly_Coef_3d<T, Dev>::pPoly_Coef_3d(const pPoly_Coef_3d<T, Dev>& pcoef_poly3)
	{
		*this = pcoef_poly3;
	}

	/* constructor from Poly_Coef_3d */
	template <class T, eDev Dev>
	explicit pPoly_Coef_3d<T, Dev>::pPoly_Coef_3d(const Poly_Coef_3d<T, Dev>& coef_poly3)
	{
		*this = coef_poly3;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDev Dev>
	CGPU_EXEC
	pPoly_Coef_3d<T, Dev>& pPoly_Coef_3d<T, Dev>::operator=(const pPoly_Coef_3d<T, Dev>& pcoef_poly3)
	{
		if (this != &pcoef_poly3)
		{
			c0 = pcoef_poly3.c0;
			c1 = pcoef_poly3.c1;
			c2 = pcoef_poly3.c2;
			c3 = pcoef_poly3.c3;
			m_size = pcoef_poly3.m_size;
		}

		return *this;
	}

	/* Assignment operator: Poly_Coef_3d -> pPoly_Coef_3d */
	template <class T, eDev Dev>
	CPU_EXEC
	pPoly_Coef_3d<T, Dev>& pPoly_Coef_3d<T, Dev>::operator=(const Poly_Coef_3d<T, Dev>& coef_poly3)
	{
		c0 = coef_poly3.c0.data();
		c1 = coef_poly3.c1.data();
		c2 = coef_poly3.c2.data();
		c3 = coef_poly3.c3.data();
		m_size = size_type(coef_poly3.size());

		return *this;
	}

	template <class T, eDev Dev>
	size_type pPoly_Coef_3d<T, Dev>::size() const
	{
		return m_size;
	}
}

namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDev Dev>
	Poly_Coef_3d<T, Dev>::Poly_Coef_3d() {}

	template <class T, eDev Dev>
	Poly_Coef_3d<T, Dev>::Poly_Coef_3d(const Vctr_cpu<T>& c0, const Vctr_cpu<T>& c1, const Vctr_cpu<T>& c2, const Vctr_cpu<T>& c3): 
		c0(c0), c1(c1), c2(c2), c3(c3){}

	template <class T, eDev Dev>
	Poly_Coef_3d<T, Dev>::Poly_Coef_3d(const dt_init_list_f64& c0, const dt_init_list_f64& c1, const dt_init_list_f64& c2, const dt_init_list_f64& c3): 
		c0(c0), c1(c1), c2(c2), c3(c3) {}

	template <class T, eDev Dev>
	Poly_Coef_3d<T, Dev>::Poly_Coef_3d(const size_type& new_size)
	{
		resize(new_size);
	}

	template <class T, eDev Dev>
	Poly_Coef_3d<T, Dev>::Poly_Coef_3d(const size_type& new_size, const T& value)
	{
		resize(new_size, value);
	}

	/* copy constructor */
	template <class T, eDev Dev>
	Poly_Coef_3d<T, Dev>::Poly_Coef_3d(const Poly_Coef_3d<T, Dev>& coef_poly3)
	{
		*this = coef_poly3;
	}

	/* converting constructor */
	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	Poly_Coef_3d<T, Dev>::Poly_Coef_3d(const Poly_Coef_3d<U, Dev_u>& coef_poly3)
	{
		*this = coef_poly3;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDev Dev>
	Poly_Coef_3d<T, Dev>& Poly_Coef_3d<T, Dev>::operator=(const Poly_Coef_3d<T, Dev>& coef_poly3)
	{
		this->assign(coef_poly3);
	
		return *this;
	}

	/* converting assignment operator */
	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	Poly_Coef_3d<T, Dev>& Poly_Coef_3d<T, Dev>::operator=(const Poly_Coef_3d<U, Dev_u>& coef_poly3)
	{
		this->assign(coef_poly3);
	
		return *this;
	}

	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	void Poly_Coef_3d<T, Dev>::assign(const Poly_Coef_3d<U, Dev_u>& coef_poly3)
	{ 
		c0.assign(coef_poly3.c0);
		c1.assign(coef_poly3.c1);
		c2.assign(coef_poly3.c2);
		c3.assign(coef_poly3.c3);
	}

	/**************** user define conversion operators *******************/
	template <class T, eDev Dev>
	pPoly_Coef_3d<T, Dev> Poly_Coef_3d<T, Dev>::ptr() const
	{
		return pPoly_Coef_3d<T, Dev>(*this);
	}

	/* user define conversion for pointer Vctr */
	template <class T, eDev Dev>
	Poly_Coef_3d<T, Dev>::operator pPoly_Coef_3d<T, Dev>() const
	{
		return pPoly_Coef_3d<T, Dev>(*this);
	}

	template <class T, eDev Dev>
	void Poly_Coef_3d<T, Dev>::fill(const T& val_c0, const T& val_c1, const T& val_c2, const T& val_c3)
	{
		c0.fill(val_c0);
		c1.fill(val_c1);
		c2.fill(val_c2);
		c3.fill(val_c3);
	}

	template <class T, eDev Dev>
	size_type Poly_Coef_3d<T, Dev>::size() const
	{
		return c0.size();
	}

	template <class T, eDev Dev>
	void Poly_Coef_3d<T, Dev>::clear()
	{
		c0.clear();
		c1.clear();
		c2.clear();
		c3.clear();
	}

	template <class T, eDev Dev>
	void Poly_Coef_3d<T, Dev>::reserve(const size_type& new_size)
	{
		c0.reserve(new_size);
		c1.reserve(new_size);
		c2.reserve(new_size);
		c3.reserve(new_size);
	}

	template <class T, eDev Dev>
	void Poly_Coef_3d<T, Dev>::resize(const size_type& new_size)
	{
		c0.resize(new_size);
		c1.resize(new_size);
		c2.resize(new_size);
		c3.resize(new_size);
	}

	template <class T, eDev Dev>
	void Poly_Coef_3d<T, Dev>::resize(const size_type& new_size, const T& value)
	{
		c0.resize(new_size, value);
		c1.resize(new_size, value);
		c2.resize(new_size, value);
		c3.resize(new_size, value);
	}

	template <class T, eDev Dev>
	void Poly_Coef_3d<T, Dev>::shrink_to_fit()
	{
		c0.shrink_to_fit();
		c1.shrink_to_fit();
		c2.shrink_to_fit();
		c3.shrink_to_fit();
	}
}