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

#include "quad_coef_2d.h"

/* class pQuad_Coef_2d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDev Dev>
	CGPU_EXEC
	pQuad_Coef_2d<T, Dev>::pQuad_Coef_2d(): x(nullptr), y(nullptr), w(nullptr), m_size(0) {}

	template <class T, eDev Dev>
	CGPU_EXEC
	pQuad_Coef_2d<T, Dev>::pQuad_Coef_2d(T* x, T* y, T* w, const size_type& size): x(x), y(y), w(w), m_size(size) {}

	/* copy constructor */
	template <class T, eDev Dev>
	CGPU_EXEC
	pQuad_Coef_2d<T, Dev>::pQuad_Coef_2d(const pQuad_Coef_2d<T, Dev>& pcoef_quad_2d)
	{
		*this = pcoef_quad_2d;
	}

	/* constructor from Quad_Coef_2d */
	template <class T, eDev Dev>
	pQuad_Coef_2d<T, Dev>::pQuad_Coef_2d(const Quad_Coef_2d<T, Dev>& coef_quad_2d)
	{
		*this = coef_quad_2d;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDev Dev>
	CGPU_EXEC
	pQuad_Coef_2d<T, Dev>& pQuad_Coef_2d<T, Dev>::operator=(const pQuad_Coef_2d<T, Dev>& pcoef_quad_2d)
	{
		if (this != &pcoef_quad_2d)
		{
			x = pcoef_quad_2d.x;
			y = pcoef_quad_2d.y;
			w = pcoef_quad_2d.w;
			m_size = pcoef_quad_2d.m_size;
		}

		return *this;
	}

	/* Assignment operator: Quad_Coef_2d -> pQuad_Coef_2d */
	template <class T, eDev Dev>
	CPU_EXEC
	pQuad_Coef_2d<T, Dev>& pQuad_Coef_2d<T, Dev>::operator=(Quad_Coef_2d<T, Dev>& coef_quad_2d)
	{
		x = coef_quad_2d.x.data();
		y = coef_quad_2d.y.data();
		w = coef_quad_2d.w.data();
		m_size = size_type(coef_quad_2d.size());

		return *this;
	}

	template <class T, eDev Dev>
	pQuad_Coef_2d<T, Dev>::size_type pQuad_Coef_2d<T, Dev>::size() const
	{
		return m_size;
	}
}

/* class Quad_Coef_2d */
namespace mt
{
	/************************************* constructors ************************************/
	template <class T, eDev Dev>
	Quad_Coef_2d<T, Dev>::Quad_Coef_2d(const Vctr_cpu<T>& x, const Vctr_cpu<T>& y, const Vctr_cpu<T>& w): 
		x(x), y(y), w(w){}

	template <class T, eDev Dev>
	Quad_Coef_2d<T, Dev>::Quad_Coef_2d(const dt_init_list_f64& x, const dt_init_list_f64& y, const dt_init_list_f64& w): 
		x(x), y(y), w(w) {}

	template <class T, eDev Dev>
	Quad_Coef_2d<T, Dev>::Quad_Coef_2d(const size_type& new_size)
	{
		resize(new_size);
	}

	template <class T, eDev Dev>
	Quad_Coef_2d<T, Dev>::Quad_Coef_2d(const size_type& new_size, const T& value)
	{
		resize(new_size, value);
	}

	/* copy constructor */
	template <class T, eDev Dev>
	Quad_Coef_2d<T, Dev>::Quad_Coef_2d(const Quad_Coef_2d<T, Dev>& coef_quad_2d)
	{
		*this = coef_quad_2d;
	}

	/* converting constructor */
	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	Quad_Coef_2d<T, Dev>::Quad_Coef_2d(const Quad_Coef_2d<U, Dev_u>& coef_quad_2d)
	{
		*this = coef_quad_2d;
	}

	/******************************** assignment operators *********************************/
	/* copy assignment operator */
	template <class T, eDev Dev>
	Quad_Coef_2d<T, Dev>& Quad_Coef_2d<T, Dev>::operator=(const Quad_Coef_2d<T, Dev>& coef_quad_2d)
	{
		this->assign(coef_quad_2d);
			
		return *this;
	}

	/* converting assignment operator */
	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	Quad_Coef_2d<T, Dev>& Quad_Coef_2d<T, Dev>::operator=(const Quad_Coef_2d<U, Dev_u>& coef_quad_2d)
	{
		this->assign(coef_quad_2d);
			
		return *this;
	}

	template <class T, eDev Dev>
	template <class U, eDev Dev_u> 
	void Quad_Coef_2d<T, Dev>::assign(const Quad_Coef_2d<U, Dev_u>& coef_quad_2d)
	{ 
		x.assign(coef_quad_2d.x);
		y.assign(coef_quad_2d.y);
		w.assign(coef_quad_2d.w);
	}

	/**************** user define conversion operators *******************/
	template <class T, eDev Dev>
	pQuad_Coef_2d<T, Dev> Quad_Coef_2d<T, Dev>::ptr() const
	{
		return pQuad_Coef_2d<T, Dev>(*this);
	}

	/* user define conversion for pointer Vctr */
	template <class T, eDev Dev>
	Quad_Coef_2d<T, Dev>::operator pQuad_Coef_2d<T, Dev>() const
	{
		return pQuad_Coef_2d<T, Dev>(*this);
	}

	template <class T, eDev Dev>
	void Quad_Coef_2d<T, Dev>::fill(const T& val_x, const T& val_y, const T& val_w)
	{
		x.fill(val_x);
		y.fill(val_y);
		w.fill(val_w);
	}

	template <class T, eDev Dev>
	Quad_Coef_2d<T, Dev>::size_type Quad_Coef_2d<T, Dev>::size() const
	{
		return x.size();
	}

	template <class T, eDev Dev>
	void Quad_Coef_2d<T, Dev>::clear()
	{
		x.clear();
		y.clear();
		w.clear();
	}

	template <class T, eDev Dev>
	void Quad_Coef_2d<T, Dev>::reserve(const size_type& new_size)
	{
		x.reserve(new_size);
		y.reserve(new_size);
		w.reserve(new_size);
	}

	template <class T, eDev Dev>
	void Quad_Coef_2d<T, Dev>::resize(const size_type& new_size)
	{
		x.resize(new_size);
		y.resize(new_size);
		w.resize(new_size);
	}

	template <class T, eDev Dev>
	void Quad_Coef_2d<T, Dev>::resize(const size_type& new_size, const T& value)
	{
		x.resize(new_size, value);
		y.resize(new_size, value);
		w.resize(new_size, value);
	}

	template <class T, eDev Dev>
	void Quad_Coef_2d<T, Dev>::shrink_to_fit()
	{
		x.shrink_to_fit();
		y.shrink_to_fit();
		w.shrink_to_fit();
	}
}