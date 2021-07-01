/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef CONST_ENUM_H
	#define CONST_ENUM_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include <type_traits>
	#include <stdint.h>
	#include <float.h>
	#include <complex>
	#include <initializer_list>

	#ifdef __CUDACC__
		#include <cuda.h>
		#include <cuda_runtime.h>
		#include <vector_types.h>
		#include <thrust/complex.h>
		using thrust::complex;
	#else
		using std::complex;
	#endif

	#include "macros.cuh"
	#include "shape_t.hpp"

	/*data type */
	using dt_bool = bool;
	using dt_int8 = int8_t;
	using dt_uint8 = uint8_t;
	using dt_int16 = int16_t;
	using dt_uint16 = uint16_t;
	using dt_int32 = int32_t;
	using dt_uint32 = uint32_t;
	using dt_int64 = int64_t;
	using dt_uint64 = uint64_t;
	using dt_float32 = float;
	using dt_float64 = double;

	using dt_cint8 = complex<dt_int8>;
	using dt_cuint8 = complex<dt_uint8>;
	using dt_cint16 = complex<dt_int16>;
	using dt_cuint16 = complex<dt_uint16>;
	using dt_cint32 = complex<dt_int32>;
	using dt_cuint32 = complex<dt_uint32>;
	using dt_cint64 = complex<dt_int64>;
	using dt_cuint64 = complex<dt_uint64>;
	using dt_cfloat32 = complex<dt_float32>;
	using dt_cfloat64 = complex<dt_float64>;

	using dt_std_cint8 = std::complex<dt_int8>;
	using dt_std_cuint8 = std::complex<dt_uint8>;
	using dt_std_cint16 = std::complex<dt_int16>;
	using dt_std_cuint16 = std::complex<dt_uint16>;
	using dt_std_cint32 = std::complex<dt_int32>;
	using dt_std_cuint32 = std::complex<dt_uint32>;
	using dt_std_cint64 = std::complex<dt_int64>;
	using dt_std_cuint64 = std::complex<dt_uint64>;
	using dt_std_cfloat32 = std::complex<dt_float32>;
	using dt_std_cfloat64 = std::complex<dt_float64>;

	#ifdef __CUDACC__
		using dt_thr_cint8 = thrust::complex<dt_int8>;
		using dt_thr_cuint8 = thrust::complex<dt_uint8>;
		using dt_thr_cint16 = thrust::complex<dt_int16>;
		using dt_thr_cuint16 = thrust::complex<dt_uint16>;
		using dt_thr_cint32 = thrust::complex<dt_int32>;
		using dt_thr_cuint32 = thrust::complex<dt_uint32>;
		using dt_thr_cint64 = thrust::complex<dt_int64>;
		using dt_thr_cuint64 = thrust::complex<dt_int64>;
		using dt_thr_cfloat32 = thrust::complex<dt_float32>;
		using dt_thr_cfloat64 = thrust::complex<dt_float64>;
	#endif

	/* data type */
	enum eData_Typ
	{
		edt_none = 0, edt_bool = 1, 
		edt_int8 = 2, edt_uint8 = 3, edt_int16 = 4, edt_uint16 = 5, edt_int32 = 6, 
		edt_uint32 = 7, edt_int64 = 8, edt_uint64 = 9, edt_float32 = 10, edt_float64 = 11, 
		edt_cint8 = 12, edt_cuint8 = 13, edt_cint16 = 14, edt_cuint16 = 15, edt_cint32 = 16, 
		edt_cuint32 = 17, edt_cint64 = 18, edt_cuint64 = 19, edt_cfloat32 = 20, edt_cfloat64 = 21, 
		edt_std_cint8 = 22, edt_std_cuint8 = 23, edt_std_cint16 = 24, edt_std_cuint16 = 25, edt_std_cint32 = 26, 
		edt_std_cuint32 = 27, edt_std_cint64 = 28, edt_std_cuint64 = 29, edt_std_cfloat32 = 30, edt_std_cfloat64 = 31, 
		edt_thr_cint8 = 32, edt_thr_cuint8 = 33, edt_thr_cint16 = 34, edt_thr_cuint16 = 35, edt_thr_cint32 = 36, 
		edt_thr_cuint32 = 37, edt_thr_cint64 = 38, edt_thr_cuint64 = 39, edt_thr_cfloat32 = 40, edt_thr_cfloat64 = 41
	};

	template <class T> 
	using Ctpr = const T* __restrict__;		// const template pointer restrict

	template <class T> 
	using Tpr = T* __restrict__;			// template pointer restrict

	template <class T>
	using dt_init_list = std::initializer_list<T>;

	using dt_init_list_int32 = std::initializer_list<dt_int32>;
	using dt_init_list_int64 = std::initializer_list<dt_int64>;

	using dt_init_list_f32 = std::initializer_list<dt_float32>;
	using dt_init_list_f64 = std::initializer_list<dt_float64>;

	/* physical constants */
	namespace mt
	{
		template<class T>
		const T c_Ha = T(27.2113850656389L);					// Hartree to electron-Volt
		
		template<class T>
		const T c_a0 = T(0.52917721077817892L);					// Bohr radius
		
		template<class T>
		const T c_pot_factor = T(47.877645145863056L);			// 
		
		template<class T>
		const T c_2Pi2a0 = T(10.445539456905012L);				// 2*pi^2*a0
		
		template<class T>
		const T c_H = T(6.62606876e-34L);						// Planck's constant - J s

		template<class T>
		const T c_bH = T(1.054571596e-34L);						// h/(2*pi) - J s
		
		template<class T>
		const T c_C = T(2.99792458e+8L);						// Velocity of light - m s^-1 
		
		template<class T>
		const T c_Qe = T(1.602176462e-19L);						// Elementary charge
		
		template<class T>
		const T c_me = T(9.10938291e-31L);						// Electron rest mass [kg]
		
		template<class T>
		const T c_mp = T(1.672621637e-27L);						// Proton rest mass [kg]
		
		template<class T>
		const T c_KB = T(1.3806504e-23L);						// Boltzmann's constant - J K^-1
		
		template<class T>
		const T c_Na = T(6.0221415e+23L);						// Avogadro's Number - mol^-1
		
		template<class T>
		const T c_E0 = T(8.854187817e-12L);						// Vacuum permittivity
	}
	
	/* pi and power constants */
	namespace mt
	{
		template<class T>
		const T c_E = T(2.7182818284590452354L);				// e (base of natural log)

		template<class T>
		const T c_pi = T(3.141592653589793238463L);				// pi
		
		template<class T>
		const T c_ipi = T(0.3183098861837906715378L);			// 1.0/pi
		
		template<class T>
		const T c_i2pi = T(1.570796326794896619231L);			// pi/2
		
		template<class T>
		const T c_i3pi = T(1.047197551196597746154L);			// pi/3
		
		template<class T>
		const T c_i4pi = T(0.7853981633974483096157L);			// pi/4
		
		template<class T>
		const T c_2pi = T(6.283185307179586476925L);			// 2*pi
		
		template<class T>
		const T c_3pi = T(9.424777960769379715388L);			// 3*pi
		
		template<class T>
		const T c_4pi = T(12.56637061435917295385L);			// 4*pi
		
		template<class T>
		const T c_pi2 = T(9.869604401089358618834L);			// pi^2
		
		template<class T>
		const T c_pi3 = T(31.00627668029982017548L);			// pi^3
		
		template<class T>
		const T c_pi4 = T(97.4090910340024372364L);				// pi^4
		
		template<class T>
		const T c_pii2 = T(1.772453850905516027298L);			// pi^(1/2)
		
		template<class T>
		const T c_pii3 = T(1.46459188756152326302L);			// pi^(1/3)
		
		template<class T>
		const T c_pii4 = T(1.331335363800389712798L);			// pi^(1/4)

		template<class T>
		const T c_2i2 = T(1.414213562373095048802L);			// 2^(1/2)
		
		template<class T>
		const T c_3i2 = T(1.732050807568877293527L);			// 3^(1/2)
		
		template<class T>
		const T c_5i2 = T(2.236067977499789696409L);			// 5^(1/2)
		
		template<class T>
		const T c_7i2 = T(2.645751311064590590502L);			// 7^(1/2)
	}

	/* others constants */
	namespace mt
	{
		const dt_int32 c_cSynCPU = 5;

		const dt_int32 c_n_atom_typ = 103;
		const dt_int32 c_n_atom_ions = 15;
		const dt_int32 c_nqz = 128;
		const dt_int32 c_nR = 128;

		const dt_float64 c_vr_min = 0.001;							// before was 0.015V

		const dt_float64 c_dflt_pos_ee = 1e-04;						// position error

		const dt_float32 c_dflt_rms3d = 0.085f;						// default 3d root mean squared displacement
		const dt_float32 c_dflt_occ = 1.0f;							// default occupancy
		const dt_int32 c_dflt_tag = 0;								// default tag
		const dt_int32 c_dflt_charge = 0;							// default charge

		const dt_int32 cSizeofI = sizeof(dt_int32);
		const dt_int32 cSizeofRD = sizeof(dt_float64);
		const dt_int32 cSizeofRF = sizeof(dt_float32);
		const dt_int32 cSizeofCD = 2*cSizeofRD;

		const dt_uint64 c_bytes_2_kb = 1024;
		const dt_uint64 c_bytes_2_mb = 1024*1024;
		const dt_uint64 c_bytes_2_gb = 1024*1024*1024;

		template<class T>
		const T c_hwhm_2_sigma = T(0.84932180028801907L);		// hwhm to sigma 1/(sqrt(2*log(2)))
		
		template<class T>
		const T c_fwhm_2_sigma = T(0.42466090014400953L);		// fwhm to sigma 1/(2*sqrt(2*log(2)))
		
		template<class T>
		const T c_iehwgd_2_sigma = T(0.70710678118654746L);		// iehwgd to sigma 1/sqrt(2)

		template<class T>
		const T c_mrad_2_rad = T(1.0e-03L);						// mrad-->rad
		
		template<class T>
		const T c_deg_2_rad = T(0.01745329251994329576924L);	// degrees-->rad
		
		template<class T>
		const T c_mm_2_angs = T(1.0e+07L);						// mm-->Angstrom
		
		template<class T>
		const T c_eV_2_keV = T(1e-03L);							// ev-->keV

		template<class T>
		const T c_cm3_A3 = T(1e24L);							// cm^3 --> A^3
	}

	/* error class and fcns */
	namespace mt
	{
		template <class T>
		CGPU_EXEC 
		T epsilon_eps(){ return 0; }

		template <>
		CGPU_EXEC 
		dt_float64 epsilon_eps(){ return 10.0*DBL_EPSILON; }

		template <>
		CGPU_EXEC 
		dt_float32 epsilon_eps(){ return 10.0*FLT_EPSILON; }

		template <class T>
		CGPU_EXEC 
		T epsilon_abs(){ return 0; }

		template <>
		CGPU_EXEC 
		dt_float64 epsilon_abs(){ return 1e-13; }

		template <>
		CGPU_EXEC 
		dt_float32 epsilon_abs(){ return 1e-5f; }

		template <class T>
		CGPU_EXEC 
		T epsilon_rel(){ return 0; }

		template <>
		CGPU_EXEC 
		dt_float64 epsilon_rel(){ return 1e-8; }

		template <>
		CGPU_EXEC 
		dt_float32 epsilon_rel(){ return 1e-4f; }


		template <class T>
		struct Epsilon
		{
			static const T eps;
			static const T abs;
			static const T rel;
		};

		template <class T>
		const T Epsilon<T>::eps = 0;

		template <class T>
		const T Epsilon<T>::abs = 0;

		template <class T>
		const T Epsilon<T>::rel = 0;

		template <>
		const dt_float64 Epsilon<dt_float64>::eps = 10.0*DBL_EPSILON;

		template <>
		const dt_float64 Epsilon<dt_float64>::abs = 1e-13;

		template <>
		const dt_float64 Epsilon<dt_float64>::rel = 1e-8;

		template <>
		const dt_float32 Epsilon<dt_float32>::eps = 10.0*FLT_EPSILON;

		template <> 
		const dt_float32 Epsilon<dt_float32>::abs = 1e-5f;

		template <> 
		const dt_float32 Epsilon<dt_float32>::rel = 1e-4f;
	}

	/* enumerations */
	namespace mt
	{
		/* Dimension*/
		enum eDim
		{
			edim_1 = 1, edim_2 = 2, edim_3 = 3
		};

		/* Distribution */
		// 1: uniform, 2: normal, 3: poisson
		enum eDist
		{
			edist_u = 1, edist_n = 2, edist_p = 3
		};

		/* Device type */
		enum eDev
		{
			edev_cpu = 1, edev_gpu = 2, edev_cpu_gpu = 3
		};

		/* functions type */
		enum eFcn_typ
		{
			efcn_cos_tap = 1, efcn_gauss = 2, efcn_exp = 3, efcn_fermi = 5, efcn_butwth = 6, efcn_hann = 7
		};

		/* modify vector */
		enum eModify_Vector
		{
			eMV_yes = 1, eMV_no = 2
		};

		/* Multem type */
		enum ePrecision
		{
			eprc_float32 = 1, eprc_float64 = 2
		};

		/* Show Data type */
		enum eShow_CData
		{
			escd_creal = 1, escs_cimag = 2, escd_cmod = 3, escd_cphase = 4
		};

		/** operation mode */
		enum eOperation_Mode
		{
			eOM_Normal = 1, eOM_Advanced = 2
		};

		enum eRot_Ctr_Typ
		{
			erct_none = 0, erct_geometric_ctr = 1, erct_user_def = 2
		};

		/* real or fourier space */
		enum eSpace
		{
			esp_real = 1, esp_fourier = 2
		};

		/* match boder */
		enum eMatch_Bdr
		{
			emb_none = 0, emb_min = 1, emb_max = 2, emb_minmax = 3
		};

		/*** output type */
		enum eOutput_Typ
		{
			eot_matlab = 1, eot_vector = 2
		};

		/* data sel type */
		enum eDat_Sel_Typ
		{
			edst_closest = 1, edst_less_than = 2, edst_greater_than = 3, edst_eless_than = 4, edst_egreater_than = 5
		};

		/* fill sel type */
		enum eFil_Sel_Typ
		{
			efst_min = 1, efst_max = 2, efst_mean = 3, efst_min_mean = 4, efst_max_mean = 5, efst_user_def = 6, efst_same_in = 7
		};

		/* structuring element */
		enum eStr_Ele
		{
			ese_disk = 1, ese_square = 2
		};

		/** operation */
		enum eOP {
			eOP_N = 1, eOP_T = 2, eOP_C = 3 
		};

		/* Exec type */
		enum eET {
			eET_Matrix = 1, eET_Vector = 2 
		};
			
		/* quad_data coefficients */

		// 1: eqt_tanh_sinh_int_n1_p1 -> int_-1^1 f(x) dx
		// 2: eqt_exp_sinh_int_0_pinfty -> int_0^infty f(x) dx
		// 3: eqt_exp_exp_int_0_pinfty -> int_0^infty f(x)exp(-x) dx
		// 4: eqt_sinh_sinh_int_ninfty_pinfty -> int_-infty^infty f(x) dx

		// 5: eqt_fourier_sin_int_0_pinfty -> int_0^infty f(x)sin(wx) dx
		// 6: eqt_fourier_cos_int_0_pinfty -> int_0^infty f(x)Cos(wx) dx

		// 7: eqt_gauss_legendre_int_n1_p1 -> int_-1^1 f(x) dx

		// 8: eqt_gauss_hermite_x0_int_ninfty_pinfty -> int_-infty^infty f(x) x^0 Exp[-x^2] dx
		// 9: eqt_gauss_hermite_x1_int_ninfty_pinfty -> int_-infty^infty f(x) |x|^1 Exp[-x^2] dx
		// 10: eqt_gauss_hermite_x2_int_ninfty_pinfty -> int_-infty^infty f(x) |x|^2 Exp[-x^2] dx

		// 11: eqt_gauss_laguerre_x0_int_0_pinfty -> int_0^infty f(x) x^0 Exp[-x] dx
		// 12: eqt_gauss_laguerre_x1_int_0_pinfty -> int_0^infty f(x) x^1 Exp[-x] dx
		// 13: eqt_gauss_laguerre_x2_int_0_pinfty -> int_0^infty f(x) x^2 Exp[-x] dx
		// 14: eqt_gauss_laguerre_xi2_int_0_pinfty -> int_0^infty f(x) Exp[-x]/Sqrt[x] dx

		// 15: eqt_legendre -> legendre, (a, b)
		// 16: eqt_chebyshev -> (a, b) ((b-x)*(x-a))^(-0.5)
		// 17: eqt_gegenbauer -> (a, b) ((b-x)*(x-a))^alpha
		// 18: eqt_jacobi -> (a, b) (b-x)^alpha*(x-a)^beta
		// 19: eqt_laguerre -> (a, inf) (x-a)^alpha*exp(-b*(x-a))
		// 20: eqt_hermite -> (-inf, inf) |x-a|^alpha*exp(-b*(x-a)^2)
		// 21: eqt_exponential -> (a, b) |x-(a+b)/2.0|^alpha
		// 22: eqt_rational -> (a, inf) (x-a)^alpha*(x+b)^beta

		enum eQuad_Typ
		{
			eqt_none = 0, eqt_tanh_sinh_int_n1_p1 = 1, eqt_exp_sinh_int_0_pinfty = 2, eqt_exp_exp_int_0_pinfty = 3, 
			eqt_sinh_sinh_int_ninfty_pinfty = 4, eqt_fourier_sin_int_0_pinfty = 5, eqt_fourier_cos_int_0_pinfty = 6, 
			eqt_gauss_legendre_int_n1_p1 = 7, eqt_gauss_hermite_x0_int_ninfty_pinfty = 8, eqt_gauss_hermite_x1_int_ninfty_pinfty = 9, 
			eqt_gauss_hermite_x2_int_ninfty_pinfty = 10, eqt_gauss_laguerre_x0_int_0_pinfty = 11, eqt_gauss_laguerre_x1_int_0_pinfty = 12, 
			eqt_gauss_laguerre_x2_int_0_pinfty = 13, eqt_gauss_laguerre_xi2_int_0_pinfty = 14, eqt_legendre = 15, 
			eqt_chebyshev = 16, eqt_gegenbauer = 17, eqt_jacobi = 18, 
			eqt_laguerre = 19, eqt_hermite = 20, eqt_exponential = 21, 
			eqt_rational = 22
		};
	}
	
	/* constant - enumeration comparison */
	namespace mt
	{
		inline
		dt_bool is_rot_ctr_none(const eRot_Ctr_Typ &type)
		{
			return type == mt::erct_none;
		}

		inline
		dt_bool is_rot_ctr_geometric_ctr(const eRot_Ctr_Typ &type)
		{
			return type == mt::erct_geometric_ctr;
		}		
		
		inline
		dt_bool is_rot_ctr_user_def(const eRot_Ctr_Typ &type)
		{
			return type == mt::erct_user_def;
		}
	}

	/* R_xd */
	namespace mt
	{
		template <class T, eDim Dim> class R_xtd;

		template <class T>
		using R_1d = T;

		template<class T, eDim Dim>
		using R_xd = typename std::conditional<Dim == edim_1, T, R_xtd<T, Dim>>::type;
	}

	/* pointers functions */
	namespace mt
	{
		template <class T>
		dt_bool fcn_is_null_ptr(T* ptr) 
		{ 
			return ptr==nullptr; 
		};
	}

	/* cuda grid and block constants */
	namespace mt
	{
 		const dt_int32 c_blk_x = 4;
		const dt_int32 c_blk_y = 16;

		const dt_int32 c_thr_1d = 256;

		const dt_int32 c_thr_2d_x = 32;
		const dt_int32 c_thr_2d_y = 8;

		const dt_int32 c_thr_3d_x = 32;
		const dt_int32 c_thr_3d_y = 4;
		const dt_int32 c_thr_3d_z = 2;
	}

	#ifdef __CUDACC__
		class D_Grid_Blk
		{
		public:
			D_Grid_Blk(): grid(), blk() {}

			D_Grid_Blk(const dim3& grid, const dim3& blk): grid(grid), blk(blk) {}

			D_Grid_Blk(const dt_int32& grid, const dt_int32& blk): grid(grid), blk(blk) {}

			dt_int32 grid_size() 
			{
				return grid.x*grid.y*grid.z; 
			};

			dt_int32 blk_size() 
			{ 
				return blk.x*blk.y*blk.z; 
			};

			// shared memory size reduction
			dt_int32 smems_red() 
			{ 
				// when there is only one warp per block, we need to allocate two warps
				// worth of shared memory so that we don't index shared memory out of bounds
				auto threads = blk_size();
				dt_int32 smems = (threads <= 32) ? 2*threads:threads;

				return smems;
			}

			dim3 blk;		// Threads
			dim3 grid;		// Blocks
		};

		// cuda dimension grid

		template <class ST>
		dt_uint32 fcn_cdg(const ST& n, const dt_int32& dn)
		{
			return static_cast<dt_uint32>((n + ST(dn)-ST(1))/ST(dn));
		}

		/***************************************************************************************/
 		dim3 fcn_cdb_size()
		{
			return dim3(mt::c_thr_1d);
		}

		template <class ST>
		dim3 fcn_cdg_size(const ST& n)
		{
			return dim3(fcn_cdg(n, mt::c_thr_1d));
		}

		/***************************************************************************************/
 		dim3 fcn_cdb_1d()
		{
			return dim3(mt::c_thr_1d);
		}

		template <class ST>
		dim3 fcn_cdg_1d(const ST& n)
		{
			return dim3(fcn_cdg(n, mt::c_thr_1d));
		}

		/***************************************************************************************/
 		dim3 fcn_cdb_2d()
		{
			return dim3(mt::c_thr_2d_x, mt::c_thr_2d_y);
		}

		template <class ST>
		dt_uint32 fcn_cdg_2d_x(const ST& n)
		{
			return fcn_cdg(n, mt::c_thr_2d_x);
		}

		template <class ST>
		dt_uint32 fcn_cdg_2d_y(const ST& n)
		{
			return fcn_cdg(n, mt::c_thr_2d_y);
		}

		template <class ST>
		dim3 fcn_cdg_2d(const ST& nx, const ST& ny)
		{
			return dim3(fcn_cdg_2d_x(nx), fcn_cdg_2d_y(ny));
		}

		/***************************************************************************************/
 		dim3 fcn_cdb_3d()
		{
			return dim3(mt::c_thr_3d_x, mt::c_thr_3d_y, mt::c_thr_3d_z);
		}

		template <class ST>
		dt_uint32 fcn_cdg_3d_x(const ST& n)
		{
			return fcn_cdg(n, mt::c_thr_3d_x);
		}

		template <class ST>
		dt_uint32 fcn_cdg_3d_y(const ST& n)
		{	
			return fcn_cdg(n, mt::c_thr_3d_y);
		}

		template <class ST>
		dt_uint32 fcn_cdg_3d_z(const ST& n)
		{	
			return fcn_cdg(n, mt::c_thr_3d_z);
		}

		template <class ST>
		dim3 fcn_cdg_3d(const ST& nx, const ST& ny, const ST& nz)
		{	
			return dim3(fcn_cdg_3d_x(nx), fcn_cdg_3d_y(ny), fcn_cdg_3d_z(nz));
		}
		
		/***************************************************************************************/
		template <class T, class ST>
		void fcn_cuda_malloc(T*& data, const ST& n_data)
		{
			cudaMalloc((void**)&data, n_data*sizeof(T));
		}

		template <class T>
		void fcn_cuda_free(T*& data)
		{
			if (data != nullptr)
			{
				cudaFree(data);
				data = nullptr;
			}
		}

	#endif
#endif