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

#pragma once

#include "macros.h"

/***************************************************************************************/
/**************************** Kahan summation algorithm ********************************/
/*           https:// en.wikipedia.org/wiki/Kahan_summation_algorithmnnn               */
/***************************************************************************************/
namespace mt
{
	template <class T>
	class KS
	{
	public:
		using value_type = T;
		
		T m_sum;
		T m_error;

		/************************************* constructors ************************************/
		/* Default constructor */
		CGPU_EXEC
		KS();

		/* Copy constructor */
		CGPU_EXEC
		KS(const T& sum, T error=0);

		/* Copy constructor */
		CGPU_EXEC
		KS(const KS<T>& ks);

		/* Move constructor */
		CGPU_EXEC
		KS(KS<T>&& ks);

		/* converting copy constructor */
		template <class U>
		CGPU_EXEC
		KS(const KS<U>& ks);

		/******************************** assignment operators *********************************/
		/* constructor from r */
		CGPU_EXEC
		KS<T>& operator=(const T& r);

		/* copy assignment operator */
		CGPU_EXEC
		KS<T>& operator=(const KS<T>& ks);

		/* Move assignment operator */
		CGPU_EXEC
		KS<T>& operator=(KS<T>&& ks);

		/* converting assignment operator */
		template <class U>
		CGPU_EXEC
		KS<T>& operator=(const KS<U>& ks);

		// user define conversion
		CGPU_EXEC
		T operator()();
			
		// user define conversion
		CGPU_EXEC
		operator T();

		/******************* Compound assignment operators *******************/
		CGPU_EXEC
		KS<T>& operator+=(T r);
		
		CGPU_EXEC
		KS<T>& operator+=(const KS<T>& r);
		
		CGPU_EXEC
		KS<T>& operator-=(T r);	

		CGPU_EXEC
		KS<T>& operator-=(const KS<T>& r);

		CGPU_EXEC
		KS<T>& operator/=(T r);

		CGPU_EXEC
		KS<T>& operator/=(const KS<T>& r);

		CGPU_EXEC
		void add(T r);
	};
}

#include "../src/kahan_sum.inl"