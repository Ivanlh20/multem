/*
 * This file is part of Multem.
 * Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
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

#include "lnl_coef.h"	

// class linear and non-linear coefficients
namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDev Dev>
	LNL_Coef<T, Dev>::LNL_Coef() {}

	template <class T, eDev Dev>
	LNL_Coef<T, Dev>::LNL_Coef(const Vctr_cpu<T>& cl, const Vctr_cpu<T>& cnl): 
		cl(cl), cnl(cnl){}

	template <class T, eDev Dev>
	LNL_Coef<T, Dev>::LNL_Coef(const dt_init_list_f64& cl, const dt_init_list_f64& cnl): 
		cl(cl), cnl(cnl) {}

	template <class T, eDev Dev>
	LNL_Coef<T, Dev>::LNL_Coef(const size_type& new_size)
	{
		resize(new_size);
	}

	template <class T, eDev Dev>
	LNL_Coef<T, Dev>::LNL_Coef(const size_type& new_size, const T& value)
	{
		resize(new_size, value);
	}

	/* copy constructor */
	template <class T, eDev Dev>
	LNL_Coef<T, Dev>::LNL_Coef(const LNL_Coef<T, Dev>& coef_lnl)
	{
		*this = coef_lnl;
	}

	/* converting constructor */
	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	LNL_Coef<T, Dev>::LNL_Coef(const LNL_Coef<U, Dev_u>& coef_lnl)
	{
		*this = coef_lnl;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDev Dev>
	LNL_Coef<T, Dev>& LNL_Coef<T, Dev>::operator=(const LNL_Coef<T, Dev>& coef_lnl)
	{
		this->assign(coef_lnl);
	
		return *this;
	}

	/* converting assignment operator */
	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	LNL_Coef<T, Dev>& LNL_Coef<T, Dev>::operator=(const LNL_Coef<U, Dev_u>& coef_lnl)
	{
		this->assign(coef_lnl);
	
		return *this;
	}

	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	void LNL_Coef<T, Dev>::assign(const LNL_Coef<U, Dev_u>& coef_lnl)
	{ 
		cl.assign(coef_lnl.cl);
		cnl.assign(coef_lnl.cnl);
	}

	/**************** user define conversion operators *******************/
	template <class T, eDev Dev>
	pLNL_Coef<T, Dev> LNL_Coef<T, Dev>::ptr() const
	{
		return pLNL_Coef<T, Dev>(*this);
	}

	/* user define conversion for pointer */
	template <class T, eDev Dev>
	LNL_Coef<T, Dev>::operator pLNL_Coef<T, Dev>() const
	{
		return pLNL_Coef<T, Dev>(*this);
	}

	template <class T, eDev Dev>
	void LNL_Coef<T, Dev>::fill(const T& val_l, const T& val_nl)
	{
		cl.fill(val_l);
		cnl.fill(val_nl);
	}

	template <class T, eDev Dev>
	size_type LNL_Coef<T, Dev>::size() const
	{
		return size_type(cl.size());
	}

	template <class T, eDev Dev>
	void LNL_Coef<T, Dev>::clear()
	{
		cl.clear();
		cnl.clear();
	}

	template <class T, eDev Dev>
	void LNL_Coef<T, Dev>::reserve(const size_type& new_size)
	{
		cl.reserve(new_size);
		cnl.reserve(new_size);
	}

	template <class T, eDev Dev>
	void LNL_Coef<T, Dev>::resize(const size_type& new_size)
	{
		cl.resize(new_size);
		cnl.resize(new_size);
	}

	template <class T, eDev Dev>
	void LNL_Coef<T, Dev>::resize(const size_type& new_size, const T& value)
	{
		cl.resize(new_size, value);
		cnl.resize(new_size, value);
	}

	template <class T, eDev Dev>
	void LNL_Coef<T, Dev>::shrink_to_fit()
	{
		cl.shrink_to_fit();
		cnl.shrink_to_fit();
	}
}

// pointer class linear and non-linear coefficients
namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDev Dev>
	CGPU_EXEC
	pLNL_Coef<T, Dev>::pLNL_Coef(): cl(nullptr), cnl(nullptr), m_size(0) {}

	template <class T, eDev Dev>
	CGPU_EXEC
	pLNL_Coef<T, Dev>::pLNL_Coef(T* cl, T* cnl, const size_type& size): cl(cl), cnl(cnl), m_size(size) {}

	/* copy constructor */
	template <class T, eDev Dev>
	CGPU_EXEC
	pLNL_Coef<T, Dev>::pLNL_Coef(const pLNL_Coef<T, Dev>& pcoef_lnl)
	{
		*this = pcoef_lnl;
	}

	/* constructor from LNL_Coef */
	template <class T, eDev Dev>
	pLNL_Coef<T, Dev>::pLNL_Coef(const LNL_Coef<T, Dev>& coef_lnl)
	{
		*this = coef_lnl;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDev Dev>
	CGPU_EXEC
	pLNL_Coef<T, Dev>& pLNL_Coef<T, Dev>::operator=(const pLNL_Coef<T, Dev>& pcoef_lnl)
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
	pLNL_Coef<T, Dev>& pLNL_Coef<T, Dev>::operator=(const LNL_Coef<T, Dev>& coef_lnl)
	{
		cl = coef_lnl.cl.data();
		cnl = coef_lnl.cnl.data();
		m_size = size_type(coef_lnl.size());

		return *this;
	}

	template <class T, eDev Dev>
	size_type pLNL_Coef<T, Dev>::size() const
	{
		return m_size;
	}
}