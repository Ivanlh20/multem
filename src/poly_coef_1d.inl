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
/************************* linear and non-linear coefficients **************************/
/***************************************************************************************/
namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDev Dev>
	CGPU_EXEC
	pLNL_Coef<T, eDev Dev>::pLNL_Coef(): cl(nullptr), cnl(nullptr), m_size(0) {}

	template <class T, eDev Dev>
	CGPU_EXEC
	pLNL_Coef<T, eDev Dev>::pLNL_Coef(T* cl, T* cnl, const size_type& size): cl(cl), cnl(cnl), m_size(size) {}

	/* copy constructor */
	template <class T, eDev Dev>
	CGPU_EXEC
	pLNL_Coef<T, eDev Dev>::pLNL_Coef(const pLNL_Coef<T, Dev>& pcoef_lnl)
	{
		*this = pcoef_lnl;
	}

	/* constructor from LNL_Coef */
	template <class T, eDev Dev>
	explicit pLNL_Coef<T, eDev Dev>::pLNL_Coef(const LNL_Coef<T, Dev>& coef_lnl)
	{
		*this = coef_lnl;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDev Dev>
	CGPU_EXEC
	pLNL_Coef<T, Dev>& pLNL_Coef<T, eDev Dev>::operator=(const pLNL_Coef<T, Dev>& pcoef_lnl)
	{
		if (this != &pcoef_lnl)
		{
			cl = pcoef_lnl.cl;
			cnl = pcoef_lnl.cnl;
			m_size = pcoef_lnl.m_size;
		}

		return *this;
	}

	/* Assignment operator: LNL_Coef -> pLNL_Coef */
	template <class T, eDev Dev>
	CPU_EXEC
	pLNL_Coef<T, Dev>& pLNL_Coef<T, eDev Dev>::operator=(const LNL_Coef<T, Dev>& coef_lnl)
	{
		cl = coef_lnl.cl.data();
		cnl = coef_lnl.cnl.data();
		m_size = size_type(coef_lnl.size());

		return *this;
	}

	template <class T, eDev Dev>
	size_type pLNL_Coef<T, eDev Dev>::size() const
	{
		return m_size;
	}
}

namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDev Dev>
	LNL_Coef<T, eDev Dev>::LNL_Coef() {}

	template <class T, eDev Dev>
	LNL_Coef<T, eDev Dev>::LNL_Coef(const Vctr_cpu<T>& cl, const Vctr_cpu<T>& cnl): 
		cl(cl), cnl(cnl){}

	template <class T, eDev Dev>
	LNL_Coef<T, eDev Dev>::LNL_Coef(const dt_init_list_f64& cl, const dt_init_list_f64& cnl): 
		cl(cl), cnl(cnl) {}

	template <class T, eDev Dev>
	LNL_Coef<T, eDev Dev>::LNL_Coef(const size_type& new_size)
	{
		resize(new_size);
	}

	template <class T, eDev Dev>
	LNL_Coef<T, eDev Dev>::LNL_Coef(const size_type& new_size, const T& value)
	{
		resize(new_size, value);
	}

	/* copy constructor */
	template <class T, eDev Dev>
	LNL_Coef<T, eDev Dev>::LNL_Coef(const LNL_Coef<T, Dev>& coef_lnl)
	{
		*this = coef_lnl;
	}

	/* converting constructor */
	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	LNL_Coef<T, eDev Dev>::LNL_Coef(const LNL_Coef<U, Dev_u>& coef_lnl)
	{
		*this = coef_lnl;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDev Dev>
	LNL_Coef<T, Dev>& LNL_Coef<T, eDev Dev>::operator=(const LNL_Coef<T, Dev>& coef_lnl)
	{
		this->assign(coef_lnl);
	
		return *this;
	}

	/* converting assignment operator */
	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	LNL_Coef<T, Dev>& LNL_Coef<T, eDev Dev>::operator=(const LNL_Coef<U, Dev_u>& coef_lnl)
	{
		this->assign(coef_lnl);
	
		return *this;
	}

	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	void LNL_Coef<T, eDev Dev>::assign(const LNL_Coef<U, Dev_u>& coef_lnl)
	{ 
		cl.assign(coef_lnl.cl);
		cnl.assign(coef_lnl.cnl);
	}

	/**************** user define conversion operators *******************/
	pLNL_Coef<T, Dev> LNL_Coef<T, eDev Dev>::ptr() const
	{
		return pLNL_Coef<T, Dev>(*this);
	}

	/* user define conversion for pointer */
	template <class T, eDev Dev>
	LNL_Coef<T, eDev Dev>::operator pLNL_Coef<T, Dev>() const
	{
		return pLNL_Coef<T, Dev>(*this);
	}

	template <class T, eDev Dev>
	void LNL_Coef<T, eDev Dev>::fill(const T& val_l, const T& val_nl)
	{
		cl.fill(val_l);
		cnl.fill(val_nl);
	}

	template <class T, eDev Dev>
	size_type LNL_Coef<T, eDev Dev>::size() const
	{
		return size_type(cl.size());
	}

	template <class T, eDev Dev>
	void LNL_Coef<T, eDev Dev>::clear()
	{
		cl.clear();
		cnl.clear();
	}

	template <class T, eDev Dev>
	void LNL_Coef<T, eDev Dev>::reserve(const size_type& new_size)
	{
		cl.reserve(new_size);
		cnl.reserve(new_size);
	}

	template <class T, eDev Dev>
	void LNL_Coef<T, eDev Dev>::resize(const size_type& new_size)
	{
		cl.resize(new_size);
		cnl.resize(new_size);
	}

	template <class T, eDev Dev>
	void LNL_Coef<T, eDev Dev>::resize(const size_type& new_size, const T& value)
	{
		cl.resize(new_size, value);
		cnl.resize(new_size, value);
	}

	template <class T, eDev Dev>
	void LNL_Coef<T, eDev Dev>::shrink_to_fit()
	{
		cl.shrink_to_fit();
		cnl.shrink_to_fit();
	}
}

/***************************************************************************************/
/*************************** linear polynomial coefficients ****************************/
/***************************************************************************************/
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
	explicit pPoly_Coef_1d<T, Dev>::pPoly_Coef_1d(const Poly_Coef_1d<T, Dev>& coef_poly1)
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
	size_type pPoly_Coef_1d<T, Dev>::size() const
	{
		return m_size;
	}
}

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
	Poly_Coef_1d(const Poly_Coef_1d<U, Dev_u>& coef_poly1)
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
	size_type Poly_Coef_1d<T, Dev>::size() const
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
	void reserve(const size_type& new_size)
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