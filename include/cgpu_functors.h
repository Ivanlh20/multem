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

#pragma once

#include "const_enum.h"
#include "math_mt.h"

namespace mt
{
	namespace cgpu_fctr
	{
		// Assigns the real part of a complex number
		template <class T>
		struct assign_real
		{
			template <class U>
			CGPU_EXEC
			T operator()(const U& x) const { return T(x.real()); }
		};
	}

	namespace cgpu_fctr
	{
		// Assigns the maximum real part based on a limit
		template <class T>
		struct assign_max_real
		{
			const T v_lim;  // Value limit
			const T val;    // Replacement value
			assign_max_real(const T& v_lim, const T& val): v_lim(v_lim), val(val) {}

			CGPU_EXEC
			T operator()(const complex<T>& x) const 
			{ 
				auto x_r = T(x.real()); 
				return (x_r >= v_lim) ? x_r : val; 
			}
		};
	}

	namespace cgpu_fctr
	{
		// Assigns the absolute value of the real part of a complex number
		template <class T>
		struct assign_abs_real
		{
			CGPU_EXEC  
			T operator()(const complex<T>& x) const { return ::abs(x.real()); }
		};
	}

	namespace cgpu_fctr
	{
		// Assigns the absolute value of a given number
		template <class T>
		struct assign_abs
		{
			template <class U>
			CGPU_EXEC
			T operator()(const U& x) const { return ::abs(x); }
		};
	}

	namespace cgpu_fctr
	{
		// Scales a given number by a constant factor 'w'
		template <class T>
		struct scale
		{
			const T w;  // Scaling factor
			scale(T w): w(w) {}

			template <class U>
			CGPU_EXEC
			T operator()(const U& x) const { return w * x; }
		};

	}

	namespace cgpu_fctr
	{
		// Computes the 2-norm of a given number using the ::norm_2 function
		template <class T>
		struct norm_2
		{
			template <class U>
			CGPU_EXEC
			T operator()(const U& x) const { return ::norm_2(x); }
		};
	}

	namespace cgpu_fctr
	{
		// Scales the 2-norm of a given number by a constant factor 'w'
		template <class T>
		struct scale_norm_2
		{
			const T w;  // Scaling factor
			scale_norm_2(T w): w(w) {}

			template <class U>
			CGPU_EXEC
			T operator()(const U& x) const { return w * ::norm_2(x); }
		};

	}

	namespace cgpu_fctr
	{
		// Adds two given numbers and returns the result
		template <class T>
		struct add
		{
			template <class U, class V>
			CGPU_EXEC
			T operator()(const U& lhs, const V &rhs) const { return lhs + rhs; }
		};
	}

	namespace cgpu_fctr
	{
		// Subtracts the second number from the first and returns the result
		template <class T>
		struct sub
		{
			template <class U, class V>
			CGPU_EXEC
			T operator()(const U& lhs, const V &rhs) const { return lhs - rhs; }
		};
	}

	namespace cgpu_fctr
	{
		// Scales the first number by 'w' and adds it to the second number
		template <class T>
		struct add_scale
		{
			const T w;  // Scaling factor
			add_scale(T w): w(w) {}

			template <class U, class V>
			CGPU_EXEC
			T operator()(const U& lhs, const V &rhs) const { return w * lhs + rhs; }
		};

	}

	namespace cgpu_fctr
	{
		// Scales the first number by 'w1' and the second number by 'w2', then adds them together
		template <class T>
		struct add_scale_i
		{
			const T w1;  // Scaling factor for the first number
			const T w2;  // Scaling factor for the second number
			add_scale_i(T w1, T w2): w1(w1), w2(w2) {}

			template <class U, class V>
			CGPU_EXEC
			T operator()(const U& lhs, const V &rhs) const { return w1 * lhs + w2 * rhs; }
		};
	}

	namespace cgpu_fctr
	{
		// Adds the the square of the norm of the first number to the second number
		template <class T>
		struct add_norm_2
		{
			template <class U, class V>
			CGPU_EXEC
			T operator()(const U& lhs, const V &rhs) const { return ::norm_2(lhs) + rhs; }
		};
	}

	namespace cgpu_fctr
	{
		// Adds the square of the norm of the first number to the square of the 2-norm of the second number
		template <class T>
		struct add_norm_2_i
		{
			template <class U, class V>
			CGPU_EXEC
			T operator()(const U& lhs, const V &rhs) const { return ::norm_2(lhs) + ::norm_2(rhs); }
		};

	}

	namespace cgpu_fctr
	{
		// Scales the square of the norm of the first number by 'w' and adds it to the second number
		template <class T>
		struct add_scale_norm_2
		{
			const T w;  // Scaling factor
			add_scale_norm_2(T w): w(w) {}

			template <class U, class V>
			CGPU_EXEC
			T operator()(const U& lhs, const V &rhs) const { return w * ::norm_2(lhs) + rhs; }
		};
	}

	namespace cgpu_fctr
	{
		// Scales the square of the norm of the first number by 'w1' and the square of the norm of the second number by 'w2', then adds them together
		template <class T>
		struct add_scale_norm_2_i
		{
			const T w1;  // Scaling factor for the square of the norm of the first number
			const T w2;  // Scaling factor for the square of the norm of the second number
			add_scale_norm_2_i(T w1, T w2): w1(w1), w2(w2) {}

			template <class U, class V>
			CGPU_EXEC
			T operator()(const U& lhs, const V &rhs) const { return w1 * ::norm_2(lhs) + w2 * ::norm_2(rhs); }
		};

	}

	namespace cgpu_fctr
	{
		// Multiplies the first number by the second number and returns the result
		template <class T>
		struct mult
		{
			template <class U, class V>
			CGPU_EXEC
			T operator()(const U& lhs, const V &rhs) const { return lhs * rhs; }
		};
	}

	namespace cgpu_fctr
	{
		// Computes the square of the norm after a shift value 'x_sft' is applied
		template <class T, class Tr>
		struct norm_2_sft
		{
			const T x_sft;
			norm_2_sft(T x_sft): x_sft(x_sft) {}

			template <class U>
			CGPU_EXEC
			Tr operator()(const U& x) const { return ::norm_2(x - x_sft); }
		};
	}

	namespace cgpu_fctr
	{
		// Computes the absolute value after a shift value 'x_sft' is applied
		template <class T, class Tr>
		struct abs_sft
		{
			const T x_sft;
			abs_sft(T x_sft): x_sft(x_sft) {}

			template <class U>
			CGPU_EXEC
			Tr operator()(const U& x) const { return ::fabs(x - x_sft); }
		};
	}

	namespace cgpu_fctr
	{
		// Divides the input by 'x_div' after subtracting a shift value 'x_sft'
		template <class T>
		struct div_sft
		{
			const T x_sft;
			const T x_div;
			div_sft(T x_sft, T x_div): x_sft(x_sft), x_div(x_div) {}

			template <class U>
			CGPU_EXEC
			T operator()(const U& x) const { return (x - x_sft) / x_div; }
		};
	}

	namespace cgpu_fctr
	{
		// Compares two elements and returns true if the first is less than the second
		template<typename T>
		struct less
		{
			CGPU_EXEC
			bool operator()(const T &lhs, const T &rhs) const
			{
				return lhs < rhs;
			}
		};
	}

	namespace cgpu_fctr
	{
		// Compares two tuples based on the element at index 'idx' and returns true if the first is less than the second
		template<dt_int32 idx>
		struct less_soa
		{
			template <class TTuple1, class TTuple2>
			bool operator()(const TTuple1& lhs, const TTuple2& rhs)
			{
				return thrust::get<idx>(lhs) < thrust::get<idx>(rhs);
			}
		};
	}

	namespace cgpu_fctr
	{
		// Finds the element closest to a given value 'm_val'
		template <class T>
		struct closest_element
		{
			const T m_val;
			closest_element(T val): m_val(val) {}

			dt_bool operator()(const T& a, const T& b)
			{	
				const T da = 1e-4;
				return fabs(a - m_val - da) < fabs(b - m_val);
			};
		};
	}

	namespace cgpu_fctr
	{
		// Binarizes the input based on a threshold 'thr'
		template <class T>
		struct binarize
		{
			const T thr;
			binarize(T thr): thr(thr) {}

			template <class U>
			CGPU_EXEC
			T operator()(const U& x) const { return (x < thr) ? T(0) : T(1); }
		};
	}

	namespace cgpu_fctr
	{
		// Applies a maximum threshold 'thr' to the input
		template <class T>
		struct threshold_max
		{
			const T thr;
			threshold_max(T thr): thr(thr) {}

			template <class U>
			CGPU_EXEC
			T operator()(const U& x) const { return fcn_max(thr, x); }
		};
	}

	namespace cgpu_fctr
	{
		// Applies a minimum threshold 'thr' to the input
		template <class T>
		struct threshold_min
		{
			const T thr;
			threshold_min(T thr): thr(thr) {}

			template <class U>
			CGPU_EXEC
			T operator()(const U& x) const { return fcn_min(thr, x); }
		};
	}

	namespace cgpu_fctr
	{
		// Applies the Anscombe forward transformation to the input
		template <class T>
		struct anscombe_fwd 
		{
			const T xs;
			anscombe_fwd(): xs(T(3.0 / 8.0)) {}

			template <class U>
			CGPU_EXEC
			T operator()(const U& x) const 
			{ 
				return fcn_max(T(0), T(2) * ::sqrt(x + xs)); 
			}
		};
	}

	namespace cgpu_fctr
	{
		// Applies the Anscombe inverse transformation to the input
		template <class T>
		struct anscombe_inv 
		{
			const T a;
			const T b;
			const T c;
			const T d;
			const T e;
			anscombe_inv(): a(T(1.0 / 4.0)), b(T(::sqrt(3.0 / 2.0) / 4.0)), 
				c(T(-11.0 / 8.0)), d(T(5.0 * ::sqrt(3.0 / 2.0) / 8)), e(T(-1.0 / 8.0)) {}

			template <class U>
			CGPU_EXEC
			T operator()(const U& x) const 
			{ 
				if (fcn_is_zero(x))
				{
					return fcn_max(T(0), a * x + e);
				}
				else
				{
					const T ix = T(1) / x;
					return fcn_max(T(0), a * x * x + e + ix * (b + ix * (c + ix * d)));
				}
			}
		};
	}

	namespace cgpu_fctr
	{
		// Calculates the radial distance index by division
		template <class T>
		struct rad_dist_ind_by_div
		{
			T r_0;
			T dr;
			dt_int32 ix_min;
			dt_int32 ix_max;
			rad_dist_ind_by_div(T r_0, T dr, dt_int32 ix_min, dt_int32 ix_max): r_0(r_0), dr(dr), ix_min(ix_min), ix_max(ix_max) {}

			CGPU_EXEC
			dt_int32 operator()(const T& r) const 
			{ 
				return fcn_bcfloor<dt_int32>((r - r_0) / dr, ix_min, ix_max); 
			}
		};
	}

	namespace cgpu_fctr
	{
		// Calculates the radial distance index by searching within a vector 'rv'
		template <class T>
		struct rad_dist_ind_by_srch
		{
			T* rv;
			dt_int32 ix_min;
			dt_int32 ix_max;
			rad_dist_ind_by_srch(T* rv, dt_int32 ix_min, dt_int32 ix_max): rv(rv), ix_min(ix_min), ix_max(ix_max) {}

			CGPU_EXEC
			dt_int32 operator()(const T& r) const 
			{ 
				return fcn_r_2_ir_by_vctr(rv, r, ix_min, ix_max); 
			}
		};
	}
} // namespace mt