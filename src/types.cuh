/*
 * This file is part of MULTEM.
 * Copyright 2016 Ivan Lobato <Ivanlh20@gmail.com>
 *
 * MULTEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MULTEM is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
 */

#ifndef TYPES_H
#define TYPES_H

#ifndef DEVICE_CALLABLE
	#ifdef __CUDACC__
		#define DEVICE_CALLABLE __host__ __device__
		#define FORCE_INLINE __forceinline__
	#else
		#define DEVICE_CALLABLE
		#define FORCE_INLINE inline
	#endif
#endif

#include <cfloat>
#include <type_traits>
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <thread>
#include <mutex>

#include "math.cuh"
#include "lin_alg_def.cuh"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/detail/raw_pointer_cast.h>
#include <thrust/tuple.h>

using std::vector;
using thrust::device_vector;
using thrust::host_vector;
using thrust::raw_pointer_cast;

namespace mt
{
	const double c_Ha = 27.2113850656389; 					// Hartree to electron-Volt
	const double c_a0 = 0.52917721077817892; 				// Bohr radius
	const double c_Potf = 47.877645145863056; 				// 
	const double c_2Pi2a0 = 10.445539456905012; 			// 2*pi^2*a0

	const double c_H = 6.62606876e-34; 						// Planck's constant - J s
	const double c_bH = 1.054571596e-34; 					// h/(2*pi) - J s
	const double c_C = 2.99792458e+8; 						// Velocity of light - m s^-1 
	const double c_Qe = 1.602176462e-19; 					// Elementary charge
	const double c_me = 9.10938291e-31; 					// Electron rest mass [kg]
	const double c_mp = 1.672621637e-27; 					// Proton rest mass [kg]
	const double c_KB = 1.3806504e-23; 						// Boltzmann's constant - J K^-1
	const double c_Na = 6.0221415e+23; 						// Avogadro's Number - mol^-1

	const double c_E = 2.7182818284590452354; 				// e (base of natural log)

	const double c_Pi = 3.141592653589793238463; 			// pi
	const double c_iPi = 0.3183098861837906715378; 			// 1.0/pi
	const double c_i2Pi = 1.570796326794896619231; 			// pi/2
	const double c_i3Pi = 1.047197551196597746154; 			// pi/3
	const double c_i4Pi = 0.7853981633974483096157; 		// pi/4
	const double c_2Pi = 6.283185307179586476925; 			// 2*pi
	const double c_3Pi = 9.424777960769379715388; 			// 3*pi
	const double c_4Pi = 12.56637061435917295385; 			// 4*pi
	const double c_Pi2 = 9.869604401089358618834; 			// pi^2
	const double c_Pi3 = 31.00627668029982017548; 			// pi^3
	const double c_Pi4 = 97.4090910340024372364; 			// pi^4
	const double c_Pii2 = 1.772453850905516027298; 			// pi^(1/2)
	const double c_Pii3 = 1.46459188756152326302; 			// pi^(1/3)
	const double c_Pii4 = 1.331335363800389712798; 			// pi^(1/4)

	const double c_2i2 = 1.414213562373095048802; 			// 2^(1/2)
	const double c_3i2 = 1.732050807568877293527; 			// 3^(1/2)
	const double c_5i2 = 2.236067977499789696409; 			// 5^(1/2)
	const double c_7i2 = 2.645751311064590590502; 			// 7^(1/2)	

	const double c_fwhm2sigma = 0.42466090014400953; 		// fwhm to sigma 

	const double c_mrad_2_rad = 1.0e-03; 					// mrad-->rad
	const double c_deg_2_rad = 0.01745329251994329576924;	// degrees-->rad
	const double c_mm_2_Ags = 1.0e+07; 						// mm-->Angstrom
	const double c_meV_2_keV = 1e-03; 						// ev-->keV

	const int c_cSynCPU = 5;

	const int c_nAtomsTypes = 103;
	const int c_nAtomsIons = 15;
	const int c_nqz = 128;
	const int c_nR = 128;

	const int c_thrnxny = 16;
	const int c_thrnxy = 256;
	const double c_Vrl = 0.015;

	#define IsFS(i, nh)	((i < nh)?i:i-2*nh)
	#define IsRS(i, nh)	((i < nh)?i+nh:i-nh)

	const int cSizeofI = sizeof(int);
	const int cSizeofRD = sizeof(double);
	const int cSizeofRF = sizeof(float);
	const int cSizeofCD = 2*cSizeofRD;

	/******************************modify vector******************************/
	enum eInput_Atoms
	{
		eIA_yes = 1, eIA_no = 2
	};
	/******************************modify vector******************************/
	enum eModify_Vector
	{
		eMV_yes = 1, eMV_no = 2
	};

	/******************************e_device type******************************/
	enum eDevice
	{
		e_host = 1, e_device = 2, e_host_device = 3
	};

	/******************************Slice memory type******************************/
	enum eSlice_Memory_Type
	{
		eSMT_Transmission = 1, eSMT_Potential = 2, eSMT_none = 3
	};

	/******************************Microscope effects*****************************/
	enum eIllumination_Model
	{
		eIM_Coherent = 1, eIM_Partial_Coherent = 2, eIM_Trans_Cross_Coef = 3, eIM_Full_Integration = 4, eIM_none = 5
	};

	/******************************Spatial and temporal***************************/
	enum eTemporal_Spatial_Incoh
	{
		eTSI_Temporal_Spatial = 1, eTSI_Temporal = 2, eTSI_Spatial = 3, eTSI_none = 4
	};

	/********************************MULTEM type**********************************/
	enum ePrecision
	{
		eP_float = 1, eP_double = 2
	};

	/**********************************Input File*********************************/
	enum eInput_File
	{
		eIF_txt = 1, eIF_pdb = 2
	};

	/***********************************Operation mode*********************************/
	enum eOperation_Mode
	{
		eOM_Normal = 1, eOM_Advanced = 2
	};

	/*****************************lens variable type*********************************/
	enum eLens_Var_Type
	{
		eLVT_off = 0, eLVT_m = 1, eLVT_f = 2, eLVT_Cs3 = 3, eLVT_Cs5 = 4, 
		eLVT_mfa2 = 5, eLVT_afa2 = 6, eLVT_mfa3 = 7, eLVT_afa3 = 8,
		eLVT_inner_aper_ang = 9, eLVT_outer_aper_ang = 10
	};

	/*****************************simulation type*********************************/
	enum eTEM_Sim_Type
	{
		eTEMST_STEM = 11, eTEMST_ISTEM = 12, 
		eTEMST_CBED = 21, eTEMST_CBEI = 22, 
		eTEMST_ED = 31, eTEMST_HRTEM = 32, 
		eTEMST_PED = 41, eTEMST_HCI = 42, 
		eTEMST_EWFS = 51, eTEMST_EWRS = 52, 
		eTEMST_EELS = 61, eTEMST_EFTEM = 62, 
		eTEMST_IWFS = 71, eTEMST_IWRS = 72, 
		eTEMST_PPFS = 81, eTEMST_PPRS = 82, 			// projected potential
		eTEMST_TFFS = 91, eTEMST_TFRS = 92, 			// transmission function
		eTEMST_PropFS = 101, eTEMST_PropRS = 102		// propagate
	};

	/***************************Projected_Potential model************************/
	enum ePhonon_Model
	{
		ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
	};

	/******************Electron specimen interaction model**********************/
	enum eElec_Spec_Int_Model
	{
		eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
	};

	/**************************Projected_Potential Slicing Type***************************/
	enum ePotential_Slicing
	{
		ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
	};

	/***************************Projected_Potential parameterization************************/
	enum ePotential_Type
	{
		ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, 
		ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6
	};

	enum eRot_Point_Type
	{
		eRPT_geometric_center = 1, eRPT_User = 2
	};

	/***************************Real or Fourier space***************************/
	enum eSpace
	{
		eS_Real = 1, eS_Reciprocal = 2
	};

	/****************************Defocus plane type*****************************/
	enum eZero_Defocus_Type
	{
		eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
	};

	/**********************************Incident Wave Type********************************/
	enum eIncident_Wave_Type
	{
		eIWT_Plane_Wave = 1, eIWT_Convergent_Wave = 2, eIWT_User_Define_Wave = 3, eIWT_Auto = 4
	};

	/******************************thickness Type*******************************/
	enum eThickness_Type
	{
		eTT_Whole_Specimen = 1, eTT_Through_Thickness = 2, eTT_Through_Slices = 3
	};

	/******************************Scanning Type********************************/
	enum eScanning_Type
	{
		eST_Line = 1, eST_Area = 2
	};
	/******************************grid Type********************************/
	enum eGrid_Type
	{
		eGT_Regular = 1, eGT_Quadratic = 2
	};

	/****************************Detector type*****************************/
	enum eDetector_Type
	{
		eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3
	};

	/**********************************Channelling type********************************/
	enum eChannelling_Type
	{
		eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3
	};

	/****************************Output type*****************************/
	enum eOutput_Type
	{
		eOT_Matlab = 1, eOT_Vector = 2
	};

	/******************Data selection type**********************/
	enum eDat_Sel_Type
	{
		eDST_Closest = 1, eDST_Less_Than = 2, eDST_Greater_Than = 3
	};

	/******************structuring element**********************/
	enum eStr_Ele
	{
		eSE_Disk = 1, eSE_Square = 2
	};

	template<class T>
	void get_bn(const T &R, const int &nR, const T &dR, const T &R_max, const bool &pbc, int &iR0, int &iRn);
	
	/*********************************Epsilon***********************************/
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
	const double Epsilon<double>::eps = 10.0*DBL_EPSILON;

	template <>
	const double Epsilon<double>::abs = 1e-13;

	template <>
	const double Epsilon<double>::rel = 1e-8;

	template <>
	const float Epsilon<float>::eps = 10.0*FLT_EPSILON;

	template <> 
	const float Epsilon<float>::abs = 1e-5f;

	template <> 
	const float Epsilon<float>::rel = 1e-4f;

	struct GridBT
	{
		dim3 Blk; 	// Blocks
		dim3 Thr; 	// Threads

	};

	/*******************forward declarations********************/
	template<class T>
	struct is_fundamental;

	template<class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T get_lambda(const T &E_0);

	template<class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T get_sigma(const T &E_0);

	template<class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T get_gamma(const T &E_0);

	/**************************vector**************************/
	template<class T, eDevice dev>
	using Vector = typename std::conditional<dev == e_host, typename std::conditional<std::is_fundamental<T>::value || 
	std::is_same<T, complex<float>>::value || std::is_same<T, complex<double>>::value, host_vector<T>, vector<T>>::type, device_vector<T>>::type;

	template<class T>
	struct rVector
	{
	public:
		using value_type = T;
		using size_type = std::size_t;

		int m_size;
		T *V;

		rVector(): m_size(0), V(nullptr){}

		size_type size() const
		{
			return m_size;
		}

		template<class TVector>
		rVector(TVector &vector)
		{
			m_size = vector.size();
			V = raw_pointer_cast(vector.data());
		}

		DEVICE_CALLABLE FORCE_INLINE
		T& operator[](const int i){ return V[i]; }

		DEVICE_CALLABLE FORCE_INLINE
		const T& operator[](const int i) const { return V[i]; }
	};

	template<class T>
	DEVICE_CALLABLE FORCE_INLINE
	double sizeMb(const int &n)
	{
		return static_cast<double>(n*sizeof(T)/1048576.0);
	}

	// static member function are not support for the cuda compiler
	template<class T>
	DEVICE_CALLABLE FORCE_INLINE
	bool isEqual(const T &a, const T &b);

	template<>
	DEVICE_CALLABLE FORCE_INLINE
	bool isEqual<int>(const int &a, const int &b)
	{
		return a == b;
	}

	template<>
	DEVICE_CALLABLE FORCE_INLINE
	bool isEqual<float>(const float &a, const float &b)
	{
		const float eps_abs = 1e-5f;
		const float eps_rel = 1e-4f;

		// Check if the numbers are really close -- needed when comparing numbers near zero.
		float diff = fabs(a - b);
		if (diff <= eps_abs)
			return true;
 
		// Otherwise fall back to Knuth's algorithm
		return diff <= ((fabs(a)<fabs(b)?fabs(b):fabs(a))*eps_rel);
	}

	template<>
	DEVICE_CALLABLE FORCE_INLINE
	bool isEqual<double>(const double &a, const double &b)
	{
		const double eps_abs = 1e-13;
		const double eps_rel = 1e-8;

		// Check if the numbers are really close -- needed when comparing numbers near zero.
		double diff = fabs(a - b);
		if (diff <= eps_abs)
			return true;
 
		// Otherwise fall back to Knuth's algorithm
		return diff <= ((fabs(a)<fabs(b)?fabs(b):fabs(a))*eps_rel);
	}

	template<class T>
	DEVICE_CALLABLE FORCE_INLINE
	bool isZero(const T &x)
	{
		return isEqual<T>(x, 0);
	}

	template<class T, class U>
	DEVICE_CALLABLE FORCE_INLINE
	bool isZero(const T &x, const U &y)
	{
		return isEqual<T>(x, 0) && isEqual<U>(y, 0);
	}

	template<class T>
	DEVICE_CALLABLE FORCE_INLINE
	bool nonZero(const T &x)
	{
		return !isEqual<T>(x, 0);
	}

	template<class T, class U>
	DEVICE_CALLABLE FORCE_INLINE
	bool nonZero(const T &x, const U &y)
	{
		return !(isEqual<T>(x, 0) || isEqual<U>(y, 0));
	}

	template<class T, class U>
	DEVICE_CALLABLE FORCE_INLINE
	T Div(const T &x, const U &y)
	{
		return (isEqual<U>(y, 0))?0:static_cast<T>(x)/static_cast<T>(y);
	}

	template<class T>
	DEVICE_CALLABLE FORCE_INLINE
	bool Check_Bound(const T &x, const T &x_min, const T &x_max)
	{
		return (x_min <= x)&&(x <= x_max);
	}

	template<class T, class U>
	DEVICE_CALLABLE FORCE_INLINE
	bool Check_Bound(const T &x, const T &x_min, const T &x_max, const U &y, const U &y_min, const U &y_max)
	{
		return (x_min <= x)&&(x <= x_max)&&(y_min <= y)&&(y <= y_max);
	}

	/************************thickness*************************/
	template<class T, eDevice dev>
	struct Thickness
	{
		using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = dev;

		size_type size() const
		{
			return islice.size();
		}

		void resize(const size_type &new_size, const value_type &value = value_type())
		{
			islice.resize(new_size, value);
			islice.shrink_to_fit();

			iatom_e.resize(new_size, value);
			iatom_e.shrink_to_fit();

			z.resize(new_size, value);
			z.shrink_to_fit();

			z_zero_def_plane.resize(new_size, value);
			z_zero_def_plane.shrink_to_fit();

			z_back_prop.resize(new_size, value);
			z_back_prop.shrink_to_fit();
		}

		template<class TThickness> 
		void assign(TThickness &thickness)
		{ 
			islice.assign(thickness.islice.begin(), thickness.islice.end());
			iatom_e.assign(thickness.iatom_e.begin(), thickness.iatom_e.end());
			z.assign(thickness.z.begin(), thickness.z.end());
			z_zero_def_plane.assign(thickness.z_zero_def_plane.begin(), thickness.z_zero_def_plane.end());
			z_back_prop.assign(thickness.z_back_prop.begin(), thickness.z_back_prop.end());
		}

		Vector<int, dev> islice; 			// slice position
		Vector<int, dev> iatom_e; 			// Last atom index
		Vector<T, dev> z; 					// z
		Vector<T, dev> z_zero_def_plane; 	// z: Zero defocus
		Vector<T, dev> z_back_prop; 		// z: Back propagation
	};

	/*************************slice thickness*****************/
	template<class T, eDevice dev>
	struct Slice
	{
		using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = dev;

		size_type size() const
		{
			return z_0.size();
		}

		void resize(const size_type &new_size, const value_type &value = value_type())
		{
			z_0.resize(new_size, value);
			z_e.resize(new_size, value);

			z_int_0.resize(new_size, value);
			z_int_e.resize(new_size, value);

			iatom_0.resize(new_size, value);
			iatom_e.resize(new_size, value);

			ithk.resize(new_size, -1);
		}

		template<class TSlice> 
		void assign(TSlice &slice)
		{ 
			z_0.assign(slice.z_0.begin(), slice.z_0.end());
			z_e.assign(slice.z_e.begin(), slice.z_e.end());
			z_int_0.assign(slice.z_int_0.begin(), slice.z_int_0.end());
			z_int_e.assign(slice.z_int_e.begin(), slice.z_int_e.end());
			iatom_0.assign(slice.iatom_0.begin(), slice.iatom_0.end());
			iatom_e.assign(slice.iatom_e.begin(), slice.iatom_e.end());
			ithk.assign(slice.ithk.begin(), slice.ithk.end());
		}

		T dz(const int &islice_0, const int &islice_e)
		{
			return (islice_e<size())?(z_e[islice_e]-z_0[islice_0]):0.0;
		}

		T dz(const int &islice)
		{
			return dz(islice, islice);
		}

		T z_m(const int &islice)
		{
			return (islice<size())?(0.5*fabs(z_e[islice]+z_0[islice])):0.0;
		}

		T dz_m(const int &islice_0, const int &islice_e)
		{
			return fabs(z_m(islice_e) - z_m(islice_0));
		}

		Vector<T, dev> z_0; 		// Initial z-position
		Vector<T, dev> z_e; 		// Final z-position
		Vector<T, dev> z_int_0; 	// Initial z-position
		Vector<T, dev> z_int_e; 	// Final z-position
		Vector<int, dev> iatom_0; 	// Index to initial z-position
		Vector<int, dev> iatom_e; 	// Index to final z-position
		Vector<int, dev> ithk;		// thickness index
	};

	/************************quadrature x**********************/
	template<class T, eDevice dev>
	struct Q1
	{
		using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = dev;

		size_type size() const
		{
			return x.size();
		}

		void resize(const size_type &new_size, const value_type &value = value_type())
		{
			x.resize(new_size, value);
			w.resize(new_size, value);
		}

		template<class TQ1> 
		void assign(TQ1 &q1)
		{ 
			x.assign(q1.x.begin(), q1.x.end());
			w.assign(q1.w.begin(), q1.w.end());
		}

		Vector<T, dev> x;
		Vector<T, dev> w;
	};

	template<class T>
	struct rQ1
	{
		using value_type = T;

		rQ1(): m_size(0), x(nullptr), w(nullptr){}

		template<class TQ1> 
		rQ1<T>& operator = (TQ1 &q1)
		{
			m_size = q1.size();
			x = raw_pointer_cast(q1.x.data());
			w = raw_pointer_cast(q1.w.data());
			return *this; 
		}

		template<class TQ1>
		rQ1(TQ1 &q1)
		{
			*this = q1;
		}

		int m_size;
		T *x;
		T *w;
	};

	/***********************quadrature xy**********************/
	template<class T, eDevice dev>
	struct Q2
	{
		using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = dev;

		size_type size() const
		{
			return x.size();
		}

		void resize(const size_type &new_size, const value_type &value = value_type())
		{
			x.resize(new_size, value);
			y.resize(new_size, value);
			w.resize(new_size, value);
		}

		template<class TQ2> 
		void assign(TQ2 &q2)
		{ 
			x.assign(q2.x.begin(), q2.x.end());
			y.assign(q2.y.begin(), q2.y.end());
			w.assign(q2.w.begin(), q2.w.end());
		}

		Vector<T, dev> x;
		Vector<T, dev> y;
		Vector<T, dev> w;
	};

	template<class T>
	struct rQ2
	{
		using value_type = T;

		rQ2(): m_size(0), x(nullptr), y(nullptr), w(nullptr){}

		template<class TQ2> 
		rQ2<T>& operator = (TQ2 &q2)
		{
			m_size = q2.size();
			x = raw_pointer_cast(q2.x.data());
			y = raw_pointer_cast(q2.y.data());
			w = raw_pointer_cast(q2.w.data());
			return *this; 
		}

		template<class TQ2>
		rQ2(TQ2 &q2)
		{
			*this = q2;
		}

		int m_size;
		T *x;
		T *y;
		T *w;
	};

	/***********Lineal and non-Lineal Coefficients************/
	template<class T, eDevice dev>
	struct PP_Coef
	{	
		using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = dev;

		size_type size() const
		{
			return cl.size();
		}

		void fill(const value_type &value = value_type())
		{
			thrust::fill(cl.begin(), cl.end(), value);
			thrust::fill(cnl.begin(), cnl.end(), value);
		}

		void resize(const size_type &new_size, const value_type &value = value_type())
		{
			cl.resize(new_size, value);
			cnl.resize(new_size, value);
		}

		template<class TPP_Coef> 
		void assign(TPP_Coef &pp_coef)
		{
			cl.assign(pp_coef.cl.begin(), pp_coef.cl.end());
			cnl.assign(pp_coef.cnl.begin(), pp_coef.cnl.end());
		}

		Vector<T, dev> cl; 	// Lineal coefficients fep
		Vector<T, dev> cnl; // Non-Lineal coefficients fep

	};

	template<class T>
	struct rPP_Coef
	{
		using value_type = T;

		rPP_Coef(): m_size(0), cl(nullptr), cnl(nullptr){}

		template<class TPP_Coef> 
		rPP_Coef<T>& operator = (TPP_Coef &rhs)
		{ 
			m_size = rhs.size();
			cl = raw_pointer_cast(rhs.cl.data());
			cnl = raw_pointer_cast(rhs.cnl.data());
			return *this; 
		}

		template<class TPP_Coef>
		rPP_Coef(TPP_Coef &pp_coef)
		{
			*this = pp_coef;
		}

		int m_size;
		T *cl;
		T *cnl;
	};

	/************Cubic interpolation coefficients*************/
	template<class T, eDevice dev>
	struct CI_Coef
	{
		using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = dev;

		size_type size() const
		{
			return c0.size();
		}

		void resize(const size_type &new_size, const value_type &value = value_type())
		{
			c0.resize(new_size, value);
			c1.resize(new_size, value);
			c2.resize(new_size, value);
			c3.resize(new_size, value);
		}
 
		template<class TCI_Coef> 
		void assign(TCI_Coef &ci_coef)
		{
			c0.assign(ci_coef.c0.begin(), ci_coef.c0.end());
			c1.assign(ci_coef.c1.begin(), ci_coef.c1.end());
			c2.assign(ci_coef.c2.begin(), ci_coef.c2.end());
			c3.assign(ci_coef.c3.begin(), ci_coef.c3.end());
		}

		Vector<T, dev> c0; 	// zero coefficient
		Vector<T, dev> c1; 	// first coefficient
		Vector<T, dev> c2; 	// second coefficient
		Vector<T, dev> c3; 	// third coefficient
	};

	template<class T>
	struct rCI_Coef
	{
		using value_type = T;

		rCI_Coef(): m_size(0), c0(nullptr), c1(nullptr), c2(nullptr), c3(nullptr){}

		template<class TCI_Coef> 
		rCI_Coef<T>& operator = (TCI_Coef &ci_coef)
		{
			m_size = ci_coef.size();
			c0 = raw_pointer_cast(ci_coef.c0.data());
			c1 = raw_pointer_cast(ci_coef.c1.data());
			c2 = raw_pointer_cast(ci_coef.c2.data());
			c3 = raw_pointer_cast(ci_coef.c3.data());
			return *this; 
		}

		template<class TCI_Coef>
		rCI_Coef(TCI_Coef &ci_coef)
		{
			*this = ci_coef;
		}

		int m_size;
		T *c0;
		T *c1;
		T *c2;
		T *c3;
	};

	/**********************STEM Detector**********************/
	template<class T, eDevice dev>
	struct Detector
	{
		using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = dev;

		Detector(): type(eDT_Circular){}

		size_type size() const
		{
			size_type size_out = 0;
			switch (type)
			{
				case eDT_Circular:
				{
					size_out = g_inner.size();
				}
				break;
				case eDT_Radial:
				{
					size_out = fx.size();
				}
					break;
				case eDT_Matrix:
				{
					size_out = fR.size();
				}
				break;
			}
			return size_out;
		}

		void resize(const size_type &new_size)
		{
			switch (type)
			{
				case eDT_Circular:
				{
					g_inner.resize(new_size);
					g_outer.resize(new_size);
				}
				break;
				case eDT_Radial:
				{
					fx.resize(new_size);
				}
					break;
				case eDT_Matrix:
				{
					fR.resize(new_size);
				}
				break;
			}

		}

		template<class TDetector> 
		void assign(TDetector &detector)
		{
			type = detector.type;
			g_inner.assign(detector.g_inner.begin(), detector.g_inner.end());
			g_outer.assign(detector.g_outer.begin(), detector.g_outer.end());

			fx.resize(detector.fx.size());
			for(auto i= 0; i<detector.fx.size(); i++)
			{
				fx[i].assign(detector.fx[i].begin(), detector.fx[i].end());
			}

			fR.resize(detector.fR.size());
			for(auto i= 0; i<detector.fR.size(); i++)
			{
				fR[i].assign(detector.fR[i].begin(), detector.fR[i].end());
			}
		}

		eDetector_Type type;					// eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3
		Vector<T, e_host> g_inner;				// Inner aperture Ang^-1
		Vector<T, e_host> g_outer;				// Outer aperture Ang^-1
		Vector<Vector<T, dev>, e_host> fx;		// radial sensitivity value
		Vector<Vector<T, dev>, e_host> fR;		// 2D sensitivity value
	};

	/********************STEM Intensity***********************/
	template<class TVector>
	struct Det_Int
	{
		using value_type = typename TVector::value_type;
		using size_type = std::size_t;

		static const eDevice device = e_host;

		size_type size() const
		{
			return image.size();
		}

		Vector<TVector, e_host> image;
	};

	/******************************range************************************/
	struct Range
	{
		int ix_0; 	// initial index
		int ix_e; 	// final index
		int iy_0; 	// initial index
		int iy_e; 	// final index
		int ixy_0; 	// initial index
		int ixy_e; 	// final index
		Range(): ix_0(0), ix_e(0), iy_0(0), iy_e(0), ixy_0(0), ixy_e(0){}

		template<class TGrid>
		Range(const TGrid &grid){ set_grid(grid); }

		template<class TGrid>
		void set_grid(const TGrid &grid)
		{
			ix_0 = 0;
			ix_e = grid.nx;
			iy_0 = 0;
			iy_e = grid.ny;
			ixy_0 = 0;
			ixy_e = grid.nxy();
		}
	};

	/******************************pair<int, val>************************************/
	template<class T>
	struct sPair
	{
		using value_type = T;

		int idx; 	// index
		T value; 	// Value

	};

	/******************************Dim 3************************************/
	struct FP_Dim
	{
		bool x;
		bool y;
		bool z;

		FP_Dim(): x(true), y(true), z(false){}

		void set(const int &Dim)
		{
			switch(Dim)
			{
				case 111:
				{
					x = true;
					y = true;
					z = true;
				}
				break;
				case 110:
				{
					x = true;
					y = true;
					z = false;
				}
				break;
				case 101:
				{
					x = true;
					y = false;
					z = true;
				}
				break;
				case 11:
				{
					x = false;
					y = true;
					z = true;
				}
				break;
				case 100:
				{
					x = true;
					y = false;
					z = false;
				}
				break;
				case 10:
				{
					x = false;
					y = true;
					z = false;
				}
				break;
				case 1:
				{
					x = false;
					y = false;
					z = true;
				}
				break;
			}
		}
	};

	/**********************closest prime factorization**********************/
	struct Prime_Num
	{
		public:
			Prime_Num()
			{
				int64_t b_2 = 2, b_3 = 3, b_5 = 5, b_7 = 7;
				int np2 = 16, np3 = 7, np5 = 5, np7 = 4;
				int64_t prime_0 = 64;
				int64_t prime_e = static_cast<int64_t>(std::pow(b_2, 16));


				number.reserve(np2*np3*np5*np7);
				for(auto ie=1; ie<7; ie++)
				{
					auto p = static_cast<int64_t>(std::pow(b_2, ie));
					number.push_back(p);
				}

				for(auto ip7=0; ip7<np7; ip7++)
				{
					auto p7 = static_cast<int64_t>(std::pow(b_7, ip7));

					for(auto ip5=0; ip5<np5; ip5++)
					{
						auto p5 = static_cast<int64_t>(std::pow(b_5, ip5));

						for(auto ip3=0; ip3<np3; ip3++)
						{
							auto p3 = static_cast<int64_t>(std::pow(b_3, ip3));

							for(auto ip2=0; ip2<np2; ip2++)
							{
								auto p2 = static_cast<int64_t>(std::pow(b_2, ip2));

								auto p = p7*p5*p3*p2;

								if((prime_0<p)&&(p<=prime_e))
								{
									number.push_back(p);
								}
							}
						}
					}
				}
				number.shrink_to_fit();
				std::sort(number.begin(), number.end()); 
			}

			int operator()(int64_t n, eDat_Sel_Type dst = eDST_Closest)
			{
				auto idx = std::min_element(number.begin(), number.end(), [&n](int64_t p0, int64_t pe){ return std::abs(n-p0)<std::abs(n-pe);});

				auto pn = static_cast<int>(*idx);

				switch(dst)
				{
					case eDST_Less_Than:
					{
						if(pn>n)
						{
							pn = static_cast<int>(*(idx-1));
						}
					}
					break;
					case eDST_Greater_Than:
					{
						if(pn<n)
						{
							pn = static_cast<int>(*(idx+1));
						}
					}
					break;
				}
				return pn;
			}

		private:
			std::vector<int64_t> number;
	};

	/***************************grid****************************/
	template<class T>
	struct Grid
	{
		using value_type = T;

		int nx; 				// Number of pixels in x direction
		int ny; 				// Number of pixels in y direction

		int nxh; 				// Half number of pixels in x direction
		int nyh; 				// Half number of pixels in y direction

		T inxy; 				// 1.0/nxy

		T lx; 					// Box m_size in x direction(Angstroms)
		T ly; 					// Box m_size in y direction(Angstroms)
		T dz; 					// slice thicknes

		bool bwl; 				// Band-width limit
		bool pbc_xy; 			// Peridic boundary contions

		T dRx; 					// x-sampling in real space
		T dRy; 					// y-sampling in real space

		T dgx; 					// x-sampling in reciprocal space
		T dgy; 					// y-sampling in reciprocal space

		T gl2_max; 				// Squared of the maximun limited frequency
		T alpha;				// 1/(1+exp(alpha*(x^2-x_c^2)))

		inline
		Grid(): nx(0), ny(0), nxh(0), nyh(0), inxy(0), 
			lx(0), ly(0), dz(0), pbc_xy(true), bwl(true), 
			dRx(0), dRy(0), dgx(0), dgy(0), gl2_max(0){}

		Grid(int nx, int ny)
		{
			lx = nx;
			ly = ny;
			dz = 0.5;
			bwl = false;
			pbc_xy = false;
			set_input_data(nx, ny, lx, ly, dz, bwl, pbc_xy);
		}

		Grid(int nx, int ny, T lx, T ly, T dz = 0.5, bool bwl = false, bool pbc_xy = false)
		{
			set_input_data(nx, ny, lx, ly, dz, bwl, pbc_xy);
		}

		inline
		void set_input_data(int nx_i, int ny_i, T lx_i, T ly_i, T dz_i, bool BWL_i, bool PBC_xy_i)
		{
			nx = nx_i;
			ny = ny_i;
			nxh = nx/2;
			nyh = ny/2;
			inxy = static_cast<T>(mt::Div(1.0, nxy()));
			lx = lx_i;
			ly = ly_i;
			dz = dz_i;
			bwl = BWL_i;
			pbc_xy = PBC_xy_i;
			dRx = mt::Div(lx, nx);
			dRy = mt::Div(ly, ny);
			dgx = mt::Div(1.0, lx);
			dgy = mt::Div(1.0, ly);
			gl2_max = pow(gl_max(), 2);

			// y = 1/(1+exp(alpha*(x^2-x_c^2)))
			T dg0 = 0.25, fg0 = 1e-02;
			alpha = log(1.0/fg0-1.0)/(pow(gl_max()+dg0, 2)-gl2_max);
		}

		template<class TGrid> 
		void assign(TGrid &grid)
		{
			set_input_data(grid.nx, grid.ny, grid.lx, grid.ly, grid.dz, grid.bwl, grid.pbc_xy);
		}

		template<class TGrid> 
		Grid<T>& operator=(TGrid &grid)
		{
			assign(grid);
			return *this; 
		}
		// Maximun limited frequency
		DEVICE_CALLABLE FORCE_INLINE
		T gl_max() const
		{
			return 2.0*g_max()/3.0;
		}	

		DEVICE_CALLABLE FORCE_INLINE
		int nxy() const { return nx*ny; }

		DEVICE_CALLABLE FORCE_INLINE
		int nx_ny_min() const { return min(nx, ny); }

		DEVICE_CALLABLE FORCE_INLINE
		int nx_ny_max() const { return max(nx, ny); }

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		int lx_ly_min() const { return min(lx, ly); }

		DEVICE_CALLABLE FORCE_INLINE
		int lx_ly_max() const { return max(lx, ly); }

		DEVICE_CALLABLE FORCE_INLINE
		T lxh() const { return 0.5*lx; }

		DEVICE_CALLABLE FORCE_INLINE
		T lyh() const { return 0.5*ly; }

		DEVICE_CALLABLE FORCE_INLINE
		int lxh_lyh_min() const { return min(lxh(), lyh()); }

		DEVICE_CALLABLE FORCE_INLINE
		int lxh_lyh_max() const { return max(lxh(), lyh()); }

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		int floor_dR_min(const T &x) const { return static_cast<int>(floor(x/dR_min())); }

		DEVICE_CALLABLE FORCE_INLINE
		int floor_dRx(const T &x) const { return static_cast<int>(floor(x/dRx)); }

		DEVICE_CALLABLE FORCE_INLINE
		int floor_dRy(const T &y) const { return static_cast<int>(floor(y/dRy)); }

		DEVICE_CALLABLE FORCE_INLINE
		int ceil_dR_min(const T &x) const { return static_cast<int>(ceil(x/dR_min())); }

		DEVICE_CALLABLE FORCE_INLINE
		int ceil_dRx(const T &x) const { return static_cast<int>(ceil(x/dRx)); }

		DEVICE_CALLABLE FORCE_INLINE
		int ceil_dRy(const T &y) const { return static_cast<int>(ceil(y/dRy)); }

		// lower bound
		DEVICE_CALLABLE FORCE_INLINE
		int lb_index_x(const T &x) const 
		{ 
			return max(0, static_cast<int>(floor(x/dRx)));
		}

		// upper bound
		DEVICE_CALLABLE FORCE_INLINE
		int ub_index_x(const T &x) const 
		{ 
			return min(nx, static_cast<int>(ceil(x/dRx)));
		}

		// lower bound
		DEVICE_CALLABLE FORCE_INLINE
		int lb_index_y(const T &y) const 
		{ 
			return max(0, static_cast<int>(floor(y/dRy)));
		}

		// upper bound
		DEVICE_CALLABLE FORCE_INLINE
		int ub_index_y(const T &y) const 
		{ 
			return min(ny, static_cast<int>(ceil(y/dRy)));
		}
		/*********************************************************/
		// Maximun frequency
		DEVICE_CALLABLE FORCE_INLINE
		T g_max() const { return ::fmin(static_cast<T>(nxh)*dgx, static_cast<T>(nyh)*dgy); }

		// Squared of the maximun frequency
		DEVICE_CALLABLE FORCE_INLINE
		T g2_max() const { return pow(g_max(), 2); }

		DEVICE_CALLABLE FORCE_INLINE
		T dR_min() const { return ::fmin(dRx, dRy); }

		DEVICE_CALLABLE FORCE_INLINE
		T dg_min() const { return ::fmin(dgx, dgy); }

		DEVICE_CALLABLE FORCE_INLINE
		int nx_dRx(const T &lx) const { return static_cast<int>(ceil(lx/dRx)); }

		DEVICE_CALLABLE FORCE_INLINE
		int ny_dRy(const T &ly) const { return static_cast<int>(ceil(ly/dRy)); }

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		int igx(const int &ix) const { return ix-nxh; }

		DEVICE_CALLABLE FORCE_INLINE
		int igy(const int &iy) const { return iy-nyh; }

		DEVICE_CALLABLE FORCE_INLINE
		T gx(const int &ix) const { return igx(ix)*dgx; }

		DEVICE_CALLABLE FORCE_INLINE
		T gy(const int &iy) const { return igy(iy)*dgy; }

		DEVICE_CALLABLE FORCE_INLINE
		T g2(const int &ix, const int &iy) const 
		{ 
			T gxi = gx(ix);
			T gyi = gy(iy);
			return (gxi*gxi + gyi*gyi);
		}

		DEVICE_CALLABLE FORCE_INLINE
		T g(const int &ix, const int &iy) const { return sqrt(g2(ix, iy)); }

		DEVICE_CALLABLE FORCE_INLINE
		T gx(const int &ix, const T &x) const { return gx(ix)-x; }

		DEVICE_CALLABLE FORCE_INLINE
		T gy(const int &iy, const T &y) const { return gy(iy)-y; }

		DEVICE_CALLABLE FORCE_INLINE
		T g2(const int &ix, const int &iy, const T &x, const T &y) const 
		{ 
			T gxi = gx(ix, x);
			T gyi = gy(iy, y);
			return (gxi*gxi + gyi*gyi);
		}

		DEVICE_CALLABLE FORCE_INLINE
		T g(const int &ix, const int &iy, const T &x, const T &y) const { return sqrt(g2(ix, iy, x, y)); }

		DEVICE_CALLABLE FORCE_INLINE
		int iRx(const int &ix) const { return ix; }

		DEVICE_CALLABLE FORCE_INLINE
		int iRy(const int &iy) const { return iy; }

		DEVICE_CALLABLE FORCE_INLINE
		T Rx(const int &ix) const { return iRx(ix)*dRx; }

		DEVICE_CALLABLE FORCE_INLINE
		T Ry(const int &iy) const { return iRy(iy)*dRy; }

		DEVICE_CALLABLE FORCE_INLINE
		T R2(const int &ix, const int &iy) const 
		{ 
			T Rxi = Rx(ix);
			T Ryi = Ry(iy);
			return (Rxi*Rxi + Ryi*Ryi);
		}

		DEVICE_CALLABLE FORCE_INLINE
		T R(const int &ix, const int &iy) const { return sqrt(R2(ix, iy)); }

		DEVICE_CALLABLE FORCE_INLINE
		T Rx(const int &ix, const T &x) const { return Rx(ix)-x; }

		DEVICE_CALLABLE FORCE_INLINE
		T Ry(const int &iy, const T &y) const { return Ry(iy)-y; }

		DEVICE_CALLABLE FORCE_INLINE
		T R2(const int &ix, const int &iy, const T &x, const T &y) const 
		{ 
			T Rxi = Rx(ix, x);
			T Ryi = Ry(iy, y);
			return (Rxi*Rxi + Ryi*Ryi);
		}

		DEVICE_CALLABLE FORCE_INLINE
		T R(const int &ix, const int &iy, const T &x, const T &y) const { return sqrt(R2(ix, iy, x, y)); }

		DEVICE_CALLABLE FORCE_INLINE
		T bwl_factor(const int &ix, const int &iy) const 
		{ 
			return inxy/(1.0+exp(alpha*(g2_shift(ix, iy)-gl2_max)));
		}

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		int igx_shift(const int &ix) const { return (ix<nxh)?ix:ix-nx; }

		DEVICE_CALLABLE FORCE_INLINE
		int igy_shift(const int &iy) const { return (iy<nyh)?iy:iy-ny; }

		DEVICE_CALLABLE FORCE_INLINE
		T gx_shift(const int &ix) const { return static_cast<T>(igx_shift(ix))*dgx; }

		DEVICE_CALLABLE FORCE_INLINE
		T gy_shift(const int &iy) const { return static_cast<T>(igy_shift(iy))*dgy; }

		DEVICE_CALLABLE FORCE_INLINE
		T g2_shift(const int &ix, const int &iy) const 
		{ 
			T gxi = gx_shift(ix);
			T gyi = gy_shift(iy);
			return gxi*gxi + gyi*gyi;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T g_shift(const int &ix, const int &iy) const { return sqrt(g2_shift(ix, iy)); }

		DEVICE_CALLABLE FORCE_INLINE
		int iRx_shift(const int &ix) const { return (ix<nxh)?ix+nxh:ix-nxh; }

		DEVICE_CALLABLE FORCE_INLINE
		int iRy_shift(const int &iy) const { return (iy<nyh)?iy+nyh:iy-nyh; }

		DEVICE_CALLABLE FORCE_INLINE
		T Rx_shift(const int &ix) const { return static_cast<T>(iRx_shift(ix))*dRx; }

		DEVICE_CALLABLE FORCE_INLINE
		T Ry_shift(const int &iy) const { return static_cast<T>(iRy_shift(iy))*dRy; }

		DEVICE_CALLABLE FORCE_INLINE
		T R2_shift(const int &ix, const int &iy) const 
		{ 
			T Rxi = Rx_shift(ix);
			T Ryi = Ry_shift(iy);
			return Rxi*Rxi + Ryi*Ryi;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T R_shift(const int &ix, const int &iy) const { return sqrt(R2_shift(ix, iy)); }

		DEVICE_CALLABLE FORCE_INLINE
		T Rx_shift(const int &ix, const T &x) const { return Rx_shift(ix)-x; }

		DEVICE_CALLABLE FORCE_INLINE
		T Ry_shift(const int &iy, const T &y) const { return Ry_shift(iy)-y; }

		DEVICE_CALLABLE FORCE_INLINE
		T R2_shift(const int &ix, const int &iy, const T &x, const T &y) const 
		{ 
			T Rxi = Rx_shift(ix, x);
			T Ryi = Ry_shift(iy, y);
			return (Rxi*Rxi + Ryi*Ryi);
		}

		DEVICE_CALLABLE FORCE_INLINE
		T R_shift(const int &ix, const int &iy, const T &x, const T &y) const { return sqrt(R2_shift(ix, iy, x, y)); }

		DEVICE_CALLABLE FORCE_INLINE
		T bwl_factor_shift(const int &ix, const int &iy) const 
		{ 
			if(bwl)
			{
				return inxy/(1.0+exp(alpha*(g2_shift(ix, iy)-gl2_max)));
			}
			else
			{
				return inxy;
			}
		}
		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		int ind_row(const int &ix, const int &iy) const { return iy*nx+ix; }

		DEVICE_CALLABLE FORCE_INLINE
		int ind_col(const int &ix, const int &iy) const { return ix*ny+iy; }

		DEVICE_CALLABLE FORCE_INLINE
		void row_col(const int &ixy, int &ix, int &iy) const 
		{ 
			ix = ixy/ny;
			iy = ixy - ix*ny;
		}

	};

	/****************************lens***************************/
	template<class T>
	struct Lens
	{
		using value_type = T;

		int m; 					// Momentum of the vortex
		T f; 					// Defocus (Å)
		T Cs3; 					// Third order spherical aberration (Å)
		T Cs5; 					// Fifth order spherical aberration (Å)
		T mfa2; 				// Twofold astigmatism (Å)
		T afa2; 				// Azimuthal angle of the twofold astigmatism (rad)
		T mfa3; 				// Threefold astigmatism (Å)
		T afa3; 				// Azimuthal angle of the threefold astigmatism (rad)
		T inner_aper_ang; 		// Inner aperture (rad);
		T outer_aper_ang; 		// Outer aperture (rad);

		T sf; 					// Defocus Spread (Å)
		int nsf; 				// Number of integration steps for the defocus Spread

		T beta; 				// Divergence semi-angle (rad)
		int nbeta; 				// Number of integration steps for the divergence semi-angle

		eZero_Defocus_Type zero_defocus_type; 	// Defocus type: eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
		T zero_defocus_plane; 	// plane

		T gamma; 				// Relativistic factor
		T lambda; 				// wavelength(Angstrom)
		T lambda2; 				// wavelength(Angstrom)^2

		T cf; 					// pi*f*lambda
		T cCs3; 				// -0.5*pi*Cs3*lambda^3
		T cCs5; 				// -pi*Cs5*lambda^5/3
		T cmfa2; 				// -pi*lambda
		T cmfa3; 				// -2*pi*lambda^3/3
		T g2_min; 				// inner_aper_ang/lambda
		T g2_max; 				// outer_aper_ang/lambda

		T sggs; 				// standard deviation
		int ngxs; 				// Number of source sampling points x
		int ngys; 				// Number of source sampling points y
		T dgxs; 				// source sampling m_size;
		T dgys; 				// source sampling m_size;
		T g2_maxs; 				// q maximum square;

		Lens(): gamma(0), lambda(0), m(0), f(0), Cs3(0), Cs5(0), 
				mfa2(0), afa2(0), mfa3(0), afa3(0), inner_aper_ang(0), 
				outer_aper_ang(0), sf(0), nsf(0), beta(0), nbeta(0), 
				lambda2(0), cf(0), cCs3(0), 	cCs5(0), cmfa2(0), 
				cmfa3(0), g2_min(0), g2_max(0), sggs(0), ngxs(0), 
				ngys(0), dgxs(0), dgys(0), g2_maxs(0){}

		void set_input_data(T E_0, Grid<T> &grid)
		{
			gamma = get_gamma(E_0);

			lambda = get_lambda(E_0);
			lambda2 = pow(lambda, 2);

			cf = (isZero(f))?0:c_Pi*f*lambda;
			cCs3 = (isZero(Cs3))?0:-c_Pi*Cs3*pow(lambda, 3)/2.0;
			cCs5 = (isZero(Cs5))?0:-c_Pi*Cs5*pow(lambda, 5)/3.0;
			cmfa2 = (isZero(mfa2))?0:-c_Pi*mfa2*lambda;
			cmfa3 = (isZero(mfa3))?0:-2.0*c_Pi*mfa3*pow(lambda, 2)/3.0;
			g2_min = (isZero(inner_aper_ang)||(inner_aper_ang<0))?0:pow(sin(inner_aper_ang)/lambda, 2);
			g2_max = (isZero(outer_aper_ang)||(outer_aper_ang<0))?grid.g2_max(): pow(sin(outer_aper_ang)/lambda, 2);

			T g0s = sin(beta)/lambda;
			sggs = g0s/c_2i2;
			T gmaxs = 3.5*sggs;
			g2_maxs = gmaxs*gmaxs;
			T dgs = gmaxs/static_cast<T>(nbeta);

			int n;
			n = (dgs<grid.dgx)?static_cast<int>(floor(grid.dgx/dgs)+1):1;
			ngxs = static_cast<int>(floor(n*gmaxs/grid.dgx)) + 1;
			dgxs = gmaxs/ngxs;

			n = (dgs<grid.dgy)?static_cast<int>(floor(grid.dgy/dgs)+1):1;
			ngys = static_cast<int>(floor(n*gmaxs/grid.dgy)) + 1;
			dgys = gmaxs/ngys;
		}

		void set_defocus(T f_i)
		{
			f = f_i;
			cf = (isZero(f))?0:c_Pi*f*lambda;
		}

		T get_zero_defocus_plane(const T &z_min, const T &z_max)
		{
			T z = 0;
			switch(zero_defocus_type)
			{
				case eZDT_First:
				{
					z = z_min;
				}
				break;
				case eZDT_Middle:
				{
					z = 0.5*(z_min + z_max);
				}
				break;
				case eZDT_Last:
				{
					z = z_max;
				}
				break;
				default:
				{
					z = zero_defocus_plane;
				}
			}
			return z;
		};

		bool is_zero_defocus_type_First()
		{
			return zero_defocus_type == eZDT_First;
		}

		bool is_zero_defocus_type_Middle()
		{
			return zero_defocus_type == eZDT_Middle;
		}

		bool is_zero_defocus_type_Last()
		{
			return zero_defocus_type == eZDT_Last;
		}

		bool is_zero_defocus_type_User_Define()
		{
			return zero_defocus_type == eZDT_User_Define;
		}

		template<class TLens> 
		void assign(TLens &lens)
		{
			m = lens.m;
			f = lens.f;
			Cs3 = lens.Cs3;
			Cs5 = lens.Cs5;
			mfa2 = lens.mfa2;
			afa2 = lens.afa2;
			mfa3 = lens.mfa3;
			afa3 = lens.afa3;
			inner_aper_ang = lens.inner_aper_ang;
			outer_aper_ang = lens.outer_aper_ang;

			sf = lens.sf;
			nsf = lens.nsf;

			beta = lens.beta;
			nbeta = lens.nbeta;

			gamma = lens.gamma;
			lambda = lens.lambda;
			lambda2 = lens.lambda2;

			cf = lens.cf;
			cCs3 = lens.cCs3;
			cCs5 = lens.cCs5;
			cmfa2 = lens.cmfa2;
			cmfa3 = lens.cmfa3;
			g2_min = lens.g2_min;
			g2_max = lens.g2_max;

			sggs = lens.sggs;
			ngxs = lens.ngxs;
			ngys = lens.ngys;
			dgxs = lens.dgxs;
			dgys = lens.dgys;
			g2_maxs = lens.g2_maxs;
		}

		template<class TLens> 
		Lens<T>& operator=(TLens &lens)
		{
			assign(lens);
			return *this; 
		}

		inline
		T gxs(const int &ix) const { return static_cast<T>(ix)*dgxs; }

		inline
		T gys(const int &iy) const { return static_cast<T>(iy)*dgys; }

		inline
		T g2s(const int &ix, const int &iy) const 
		{ 
			T gxi = gxs(ix);
			T gyi = gys(iy);
			return gxi*gxi + gyi*gyi;
		}
	};

	/****************************EELS***************************/
	template<class T>
	struct EELS
	{
		using value_type = T;

		inline
		EELS(): space(eS_Real), E_0(0), E_loss(0), ge(0), ge2(0), gc(0), gc2(0), 
		m_selection(0), collection_angle(0), channelling_type(eCT_Double_Channelling), 
		g_collection(0), factor(0), Z(0), x(0), y(0), occ(0){}

		inline
		void set_input_data(eSpace space_i, T E_0_i, T E_loss_i, int m_selection_i, T collection_angle_i, eChannelling_Type channelling_type_i, int Z_i)
		{
			space = space_i;
			E_0 = E_0_i;
			E_loss = E_loss_i;
			T gamma = get_gamma(E_0);
			T lambda = get_lambda(E_0);

			T theta = get_theta(E_loss, E_0);

			ge = theta/(lambda*gamma*gamma);
			ge2 = pow(gamma*ge, 2);
			gc = sqrt(2.0*theta)/lambda;
			gc2 = pow(gc, 2);

			m_selection = m_selection_i;
			collection_angle = collection_angle_i;
			g_collection = collection_angle/lambda;
			channelling_type = channelling_type_i;

			Z = Z_i;
		}

		template<class TEELS> 
		void assign(TEELS &eels)
		{
			set_input_data(eels.space, eels.E_0, eels.E_loss, eels.m_selection, eels.collection_angle, eels.channelling_type, eels.Z);
			factor = eels.factor;
			x = eels.x;
			y = eels.y;
			occ = eels.occ;
		}

		template<class TEELS> 
		EELS<T>& operator=(TEELS &eels)
		{
			assign(eels);
			return *this; 
		}

		// effective scattering angle
		inline
		T get_theta(const T &E_loss, const T &E_0)
		{
			T emass = 510.99906;
			T x = (emass + E_0)/(2*emass + E_0);
			return E_loss*x/E_0;
		}

		bool is_Single_Channelling() const
		{
			return channelling_type == eCT_Single_Channelling;
		}

		bool is_Mixed_Channelling() const
		{
			return channelling_type == eCT_Mixed_Channelling;
		}

		bool is_Double_Channelling() const
		{
			return channelling_type == eCT_Double_Channelling;
		}

		eSpace space;
		T E_0;
		T E_loss; 				// Energy loss
		T ge; 					// relativistic corrected effective scattering momentum transfer
		T ge2; 					// ge square
		T gc; 					// ge cut-off at the Bethe Ridge
		T gc2; 					// gc square
		int m_selection;		// selection rule
		T collection_angle; 	// Collection half angle(rad)
		eChannelling_Type channelling_type;

		T factor;
		int Z;					// Atomic type
		T x;
		T y;
		T occ;
		T g_collection;
	};

	/**************************atoms for Superposition**************************/
	template<class T>
	class Atom_Data_Sp{
		public:
			using value_type = T;
			using size_type = std::size_t;

			Atom_Data_Sp(): l_x(0), l_y(0), 
				x_min(0), x_max(0), 
				y_min(0), y_max(0),
				a_min(0), a_max(0), 
				sigma_min(0), sigma_max(0), 
				b_min(0), b_max(0), 
				x_mean(0), y_mean(0), 
				x_std(0), y_std(0), 
				s_x(0), s_y(0){}

			size_type size() const
			{
				return x.size();
			}

			bool empty() const
			{
				return size() == 0;
			}

			// resize number of atoms
			void resize(const size_type &new_size, const value_type &value = value_type())
			{
				x.resize(new_size, value);
				x.shrink_to_fit();

				y.resize(new_size, value);
				y.shrink_to_fit();

				a.resize(new_size, value);
				a.shrink_to_fit();

				sigma.resize(new_size, value);
				sigma.shrink_to_fit();

				b.resize(new_size, value);
				b.shrink_to_fit();
			}

			// set atoms
			void set_Atoms(const size_type &nr_atoms_i, const size_type &nc_atoms_i, double *atoms_i, T l_x_i = 0, T l_y_i = 0, bool PBC_xy_i =false)
			{
				resize(nr_atoms_i);

				l_x = l_x_i;
				l_y = l_y_i;

				T dl = 1e-04;
				T lx_b = l_x - dl;
				T ly_b = l_y - dl;

				size_type j = 0;
				for(auto i = 0; i < size(); i++)
				{
					auto atom = read_atom(nr_atoms_i, nc_atoms_i, atoms_i, i);
					if((!PBC_xy_i)||((atom.x<lx_b)&&(atom.y<ly_b)))
					{
						x[j] = atom.x; 				// x-position
						y[j] = atom.y; 				// y-position
						a[j] = atom.a; 				// height
						sigma[j] = atom.sigma;		// standard deviation
						b[j] = atom.b;				// background
						j++;
					}
				}

				resize(j);

				get_Statistic();
			}

			// set atoms
			void set_Atoms(const Atom_Data_Sp<T> &atoms, bool PBC_xy_i=false)
			{
				resize(atoms.size());

				l_x = atoms.l_x;
				l_y = atoms.l_y;

				T dl = 1e-04;
				T lx_b = l_x - dl;
				T ly_b = l_y - dl;

				size_type j = 0;
				for(auto i = 0; i < size(); i++)
				{
					if((!PBC_xy_i)||((atoms.x[i]<lx_b)&&(atoms.y[i]<ly_b)))
					{
						x[j] = atoms.x[i]; 				// x-position
						y[j] = atoms.y[i]; 				// y-position
						a[j] = atoms.a[i];				// height
						sigma[j] = atoms.sigma[i];		// standard deviation
						b[j] = atoms.b[i];				// background
						j++;
					}
				}

				resize(j);

				get_Statistic();
			}
		
			// get statistic
			void get_Statistic()
			{
				if(empty())
				{
					return;
				}

				x_min = x_max = x[0];
				y_min = y_max = y[0];
				a_min = a_max = a[0];
				sigma_min = sigma_max = sigma[0];
				b_min = b_max = b[0];

				x_mean = y_mean = 0.0;
				x_std = y_std = 0.0;

				for(auto iAtom = 0; iAtom < size(); iAtom++)
				{
					x_min = min(x[iAtom], x_min);
					x_max = max(x[iAtom], x_max);

					y_min = min(y[iAtom], y_min);
					y_max = max(y[iAtom], y_max);

					a_min = min(a[iAtom], a_min);
					a_max = max(a[iAtom], a_max);

					sigma_min = min(sigma[iAtom], sigma_min);
					sigma_max = max(sigma[iAtom], sigma_max);

					b_min = min(b[iAtom], b_min);
					b_max = max(b[iAtom], b_max);

					x_mean += x[iAtom];
					y_mean += y[iAtom];

					x_std += x[iAtom]*x[iAtom];
					y_std += y[iAtom]*y[iAtom];
				}

				T nAtoms = static_cast<T>(size());

				x_mean /= nAtoms;
				y_mean /= nAtoms;

				x_std = sqrt(x_std/nAtoms - x_mean*x_mean);
				y_std = sqrt(y_std/nAtoms - y_mean*y_mean);

				s_x = x_max - x_min;
				s_y = y_max - y_min;

				if(isZero(l_x))
				{
					l_x = s_x;
				}

				if(isZero(l_y))
				{
					l_y = s_y;
				}
			}

			T l_x; 			// Box m_size-x
			T l_y; 			// Box m_size-y

			Vector<T, e_host> x;
			Vector<T, e_host> y;
			Vector<T, e_host> a;
			Vector<T, e_host> sigma;
			Vector<T, e_host> b;

			T x_min;
			T x_max;

			T y_min;
			T y_max;

			T a_min;
			T a_max;

			T sigma_min;
			T sigma_max;

			T b_min;
			T b_max;

			T x_mean;
			T y_mean;

			T x_std;
			T y_std;

			T s_x; 			// m_size-x
			T s_y; 			// m_size-y

		private:
			struct Atom
			{
				T x;
				T y;
				T a;
				T sigma;
				T b;

				Atom():x(0), y(0), a(0), sigma(0), b(0){};
			};

			template<class TIn>
			Atom read_atom(const int &nr, const int &nc, TIn *atoms, const int &iatom)
			{
				Atom atom;
				atom.x = atoms[0*nr + iatom]; 						// x-position
				atom.y = atoms[1*nr + iatom]; 						// y-position
				atom.a = (nc>2)?atoms[2*nr + iatom]:1.0;			// height
				atom.sigma = (nc>3)?atoms[3*nr + iatom]:1.0;		// standard deviation
				atom.b = (nc>4)?atoms[4*nr + iatom]:0.0; 			// background

				return atom;
			}
	};

	template<class T, eDevice dev>
	struct Atom_Sp
	{
		public:
			using value_type = T;
			static const eDevice device = dev;

			T x;
			T y;
			T occ;
			T R2_max;
			T *R2;
			T *c3;
			T *c2;
			T *c1;
			T *c0;

			int ix0;
			int ixn;
			int iy0;
			int iyn;

			T a;
			T alpha;
			T b;
			T dtR;
			T R2_tap;
			T tap_cf;

			int *iv;
			T *v;

			Atom_Sp(): x(0), y(0), occ(1), R2_max(0), R2(nullptr), c3(nullptr), c2(nullptr), 
			c1(nullptr), c0(nullptr), ix0(1), ixn(0), iy0(0), a(0), alpha(0), b(0), 
			dtR(0), R2_tap(0), tap_cf(0), iyn(0), iv(nullptr), v(nullptr){}

			inline
			void set_ix0_ixn(const Grid<T> &grid, const T &R_max)
			{
				get_bn(x, grid.nx, grid.dRx, R_max, grid.pbc_xy, ix0, ixn);
			}

			inline
			void set_iy0_iyn(const Grid<T> &grid, const T &R_max)
			{
				get_bn(y, grid.ny, grid.dRy, R_max, grid.pbc_xy, iy0, iyn);
			}

			inline
			GridBT get_eval_cubic_poly_gridBT()
			{
				GridBT gridBT;
				gridBT.Blk = dim3((iyn+c_thrnxny-1)/c_thrnxny, (ixn+c_thrnxny-1)/c_thrnxny);
				gridBT.Thr = dim3(c_thrnxny, c_thrnxny);
				return gridBT;
			}
	};

	/***********************atoms for simulated annealing***********************/
	template<class T>
	class Atom_Data_Sa{
		public:
			using value_type = T;
			using size_type = std::size_t;

			size_type size() const
			{
				return Z.size();
			}

			bool empty() const
			{
				return size() == 0;
			}

			template<class TAtom_SA> 
			void assign(TAtom_SA &atom_sa)
			{ 
				Z.assign(atom_sa.Z.begin(), atom_sa.Z.end());

				r_min.assign(atom_sa.r_min.begin(), atom_sa.r_min.end());
				r_max.assign(atom_sa.r_max.begin(), atom_sa.r_max.end());
				r_0.assign(atom_sa.r_0.begin(), atom_sa.r_0.end());
				r_d.assign(atom_sa.r_d.begin(), atom_sa.r_d.end());

				r.assign(atom_sa.r.begin(), atom_sa.r.end());
				r_n.assign(atom_sa.r_n.begin(), atom_sa.r_n.end());
				r_opt.assign(atom_sa.r_opt.begin(), atom_sa.r_opt.end());

				chi2.assign(atom_sa.chi2.begin(), atom_sa.chi2.end());
				chi2_n.assign(atom_sa.chi2_n.begin(), atom_sa.chi2_n.end());
				chi2_opt.assign(atom_sa.chi2_opt.begin(), atom_sa.chi2_opt.end());

				df.assign(atom_sa.df.begin(), atom_sa.df.end());
			}

			// resize number of atoms
			void resize(const size_type &new_size, const value_type &value = value_type())
			{
				Z.resize(new_size, value);

				r_min.resize(new_size, value);
				r_max.resize(new_size, value);
				r_0.resize(new_size, value);
				r_d.resize(new_size, value);

				r.resize(new_size, value);
				r_n.resize(new_size, value);
				r_opt.resize(new_size, value);

				chi2.resize(new_size, value);
				chi2_n.resize(new_size, value);
				chi2_opt.resize(new_size, value);

				df.resize(new_size, value);

			}

			// set atoms
			void set_Atoms(const size_type &natoms_i, double *atoms_i, double *atoms_min_i, double *atoms_max_i)
			{
				resize(natoms_i);

				for(auto iatoms = 0; iatoms < size(); iatoms++)
				{
					Z[iatoms] = static_cast<int>(atoms_i[0*natoms_i + iatoms]); 		// Atomic number
					r[iatoms].x = atoms_i[1*natoms_i + iatoms]; 						// x-position
					r[iatoms].y = atoms_i[2*natoms_i + iatoms]; 						// y-position
					r[iatoms].z = atoms_i[3*natoms_i + iatoms]; 						// z-position

					r_min[iatoms].x = atoms_min_i[0*natoms_i + iatoms]; 				// x-position
					r_min[iatoms].y = atoms_min_i[1*natoms_i + iatoms]; 				// y-position
					r_min[iatoms].z = atoms_min_i[2*natoms_i + iatoms]; 				// z-position

					r_max[iatoms].x = atoms_max_i[0*natoms_i + iatoms]; 				// x-position
					r_max[iatoms].y = atoms_max_i[1*natoms_i + iatoms]; 				// y-position
					r_max[iatoms].z = atoms_max_i[2*natoms_i + iatoms]; 				// z-position

					r_0[iatoms] = r_min[iatoms];
					r_d[iatoms] = r_max[iatoms]-r_min[iatoms];

					df[iatoms] = 1;
				}
			}

			// set atoms
			void set_Atoms(const size_type &natoms_i, double *atoms_i, r3d<T> d_i)
			{
				resize(natoms_i);

				for(auto iatoms = 0; iatoms < size(); iatoms++)
				{
					Z[iatoms] = static_cast<int>(atoms_i[0*natoms_i + iatoms]); 		// Atomic number

					r[iatoms].x = atoms_i[0*natoms_i + iatoms]; 						// x-position
					r[iatoms].y = atoms_i[1*natoms_i + iatoms]; 						// y-position
					r[iatoms].z = atoms_i[2*natoms_i + iatoms]; 						// z-position

					r_min[iatoms] = r - d_i;
					r_max[iatoms] = r + d_i;

					r_0[iatoms] = r_min[iatoms];
					r_d[iatoms] = r_max[iatoms]-r_min[iatoms];

					df[iatoms] = 1;
				}
			}

			void set_range(int Z_i, r3d<T> r_min_i, r3d<T> r_max_i)
			{
				for(auto iatoms = 0; iatoms < size(); iatoms++)
				{
					Z[iatoms] = Z_i;

					r_min[iatoms] = r_min_i;
					r_max[iatoms] = r_max_i;
					r_0[iatoms] = r_min[iatoms];
					r_d[iatoms] = r_max[iatoms]-r_min[iatoms];
					df[iatoms] = 1;
				}
			}

			Vector<int, e_host> Z;

			Vector<r3d<T>, e_host> r_min;
			Vector<r3d<T>, e_host> r_max;
			Vector<r3d<T>, e_host> r_0;
			Vector<r3d<T>, e_host> r_d;

			Vector<r3d<T>, e_host> r;
			Vector<r3d<T>, e_host> r_n;
			Vector<r3d<T>, e_host> r_opt;

			Vector<T, e_host> chi2;
			Vector<T, e_host> chi2_n;
			Vector<T, e_host> chi2_opt;

			Vector<T, e_host> df;
	};

	template<class T>
	struct Atom_Sa
	{
		public:
			using value_type = T;

			T x;
			T y;
			T R2_max;
			T *R2;
			T *c3;
			T *c2;
			T *c1;
			T *c0;

			int ix0;
			int ixn;
			int iy0;
			int iyn;

			int *iv;
			T *v;

			Atom_Sa(): x(0), y(0), R2_max(0), R2(nullptr), c3(nullptr), c2(nullptr), c1(nullptr), 
			c0(nullptr), ix0(1), ixn(0), iy0(0), iyn(0), iv(nullptr), v(nullptr){}

			inline
			void set_ix0_ixn(const Grid<T> &grid, const T &R_max)
			{
				get_bn(x, grid.nx, grid.dRx, R_max, grid.pbc_xy, ix0, ixn);
			}

			inline
			void set_iy0_iyn(const Grid<T> &grid, const T &R_max)
			{
				get_bn(y, grid.ny, grid.dRy, R_max, grid.pbc_xy, iy0, iyn);
			}

			inline
			GridBT get_eval_cubic_poly_gridBT()
			{
				GridBT gridBT;
				gridBT.Blk = dim3((iyn+c_thrnxny-1)/c_thrnxny, (ixn+c_thrnxny-1)/c_thrnxny);
				gridBT.Thr = dim3(c_thrnxny, c_thrnxny);
				return gridBT;
			}
	};

	/*****************************Potential**************************/
	template<class T, eDevice dev>
	struct Atom_Vp
	{
	public:
		using value_type = T;
		static const eDevice device = dev;

		int charge;
		T x;
		T y;
		T z0h;
		T zeh;
		bool split;
		T occ;
		T R2_min;
		T R2_max;
		T *R2;
		T *cl;
		T *cnl;
		T *c3;
		T *c2;
		T *c1;
		T *c0;
		int ix0;
		int ixn;
		int iy0;
		int iyn;

		T R2_tap;
		T tap_cf;

		int *iv;
		T *v;

		Atom_Vp() : charge(0), x(0), y(0), z0h(0), zeh(0), split(false), occ(1), R2_min(0), R2_max(0),
			R2(nullptr), cl(nullptr), cnl(nullptr), c3(nullptr), c2(nullptr), c1(nullptr),
			c0(nullptr), ix0(1), ixn(0), iy0(0), iyn(0), R2_tap(0), tap_cf(0), iv(nullptr), v(nullptr){}

		inline
			void set_ix0_ixn(const Grid<T> &grid, const T &R_max)
		{
			get_bn(x, grid.nx, grid.dRx, R_max, grid.pbc_xy, ix0, ixn);
		}

		inline
			void set_iy0_iyn(const Grid<T> &grid, const T &R_max)
		{
			get_bn(y, grid.ny, grid.dRy, R_max, grid.pbc_xy, iy0, iyn);
		}

		inline
			GridBT get_eval_cubic_poly_gridBT()
		{
			GridBT gridBT;
			gridBT.Blk = dim3((iyn + c_thrnxny - 1) / c_thrnxny, (ixn + c_thrnxny - 1) / c_thrnxny);
			gridBT.Thr = dim3(c_thrnxny, c_thrnxny);
			return gridBT;
		}

	};

	/*****************************Atomic Coefficients**************************/
	template<class T, eDevice dev>
	struct Atom_Coef
	{
		using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = dev;

		Atom_Coef(): charge(0), tag(0), R_min(0), R_max(0), R_tap(0), tap_cf(0){}

		template<class TAtom_Coef> 
		void assign(TAtom_Coef &atom_coef)
		{ 
			charge = atom_coef.charge;
			tag = atom_coef.tag;

			R_min = atom_coef.R_min;
			R_max = atom_coef.R_max;

			R_tap = atom_coef.R_tap;
			tap_cf = atom_coef.tap_cf;

			feg.assign(atom_coef.feg);
			fxg.assign(atom_coef.fxg);
			Pr.assign(atom_coef.Pr);
			Vr.assign(atom_coef.Vr);
			VR.assign(atom_coef.VR);

			R.assign(atom_coef.R.begin(), atom_coef.R.end());
			R2.assign(atom_coef.R2.begin(), atom_coef.R2.end());
			ciVR.assign(atom_coef.ciVR);
		}

		template<class TAtom_Coef> 
		Atom_Coef<T, dev>& operator=(TAtom_Coef &atom_coef)
		{
			assign(atom_coef);
			return *this; 
		}

		// Minimum interaction radius squared
		T R2_min() const { return pow(R_min, 2); }

		// Maximum interaction radius squared
		T R2_max() const { return pow(R_max, 2); }

		// Tapering radius squared
		T R2_tap() const { return pow(R_tap, 2); }

		int charge; 				// Charge
		T tag; 						// tag

		T R_min; 					// Minimum interaction radius
		T R_max; 					// Maximum interaction radius
		T R_tap; 					// Tapering radius
		T tap_cf; 					// Tapering cosine factor

		PP_Coef<T, dev> feg; 		// Electron scattering factor coefficients
		PP_Coef<T, dev> fxg; 		// X-ray scattering factor coefficients
		PP_Coef<T, dev> Pr; 		// Projected_Potential coefficients
		PP_Coef<T, dev> Vr; 		// Projected_Potential coefficients
		PP_Coef<T, dev> VR; 		// Projected potential coefficients

		Vector<T, dev> R; 			// R
		Vector<T, dev> R2; 			// R2
		CI_Coef<T, dev> ciVR; 		// Look up table - Projected potential coefficients

	};

	/********************************Atomic type*******************************/
	template<class T, eDevice dev>
	struct Atom_Type
	{
		using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = dev;

		Atom_Type(): Z(0), m(0), A(0), rn_e(0), rn_c(0), ra_e(0), ra_c(0){}

		template<class TAtom_Type> 
		void assign(TAtom_Type &atom_type)
		{ 
			Z = atom_type.Z;
			m = atom_type.m;
			A = atom_type.A;
			rn_e = atom_type.rn_e;
			rn_c = atom_type.rn_c;
			ra_e = atom_type.ra_e;
			ra_c = atom_type.ra_c;

			coef.resize(atom_type.coef.size());
			for(auto i= 0; i<atom_type.coef.size(); i++)
			{
				coef[i].assign(atom_type.coef[i]);
			}
		}

		template<class TAtom_Type> 
		Atom_Type<T, dev>& operator=(TAtom_Type &atom_type)
		{
			assign(atom_type);
			return *this; 
		}

		int check_charge(const int &charge) const
		{ 
			for(auto i= 0; i<coef.size(); i++)
			{
				if(coef[i].charge==charge)
				{
					return charge;
				}
			}
			return 0;
		};

		int charge_to_idx(const int &charge) const
		{ 
			int icharge = 0;
			for(auto i= 0; i<coef.size(); i++)
			{
				if(coef[i].charge==charge)
				{
					icharge = i;
					break;
				}
			}
			return icharge;
		};

		PP_Coef<T, dev>* feg(const int &charge)
		{
			int icharge = charge_to_idx(charge);
			return &(coef[icharge].feg);
		};

		PP_Coef<T, dev>* fxg(const int &charge)
		{
			int icharge = charge_to_idx(charge);
			return &(coef[icharge].fxg);
		};

		PP_Coef<T, dev>* Pr(const int &charge)
		{
			int icharge = charge_to_idx(charge);
			return &(coef[icharge].Pr);
		};

		PP_Coef<T, dev>* Vr(const int &charge)
		{
			int icharge = charge_to_idx(charge);
			return &(coef[icharge].Vr);
		};

		PP_Coef<T, dev>* VR(const int &charge)
		{
			int icharge = charge_to_idx(charge);
			return &(coef[icharge].VR);
		};

		int Z; 										// Atomic number
		T m; 										// Atomic mass
		int A; 										// Mass number
		T rn_e; 									// Experimental Nuclear radius
		T rn_c; 									// Calculated Nuclear radius
		T ra_e; 									// Experimental atomic radius
		T ra_c; 									// Calculated atomic radius

		Vector<Atom_Coef<T, dev>, e_host> coef;		// atomic coefficients
	};

	/********************************Scanning**********************************/
	template<class T>
	struct Scanning
	{
		public:
			using value_type = T;
			using size_type = std::size_t;

			eScanning_Type type;			// 1: Line, 2: Area, 
			eGrid_Type grid_type;			// 1: regular, 2: quadratic
			bool pbc;						// periodic boundary conditions
			int ns; 						// Number of sampling points
			int nx;
			int ny;
			T x0; 							// Initial scanning position in x
			T y0; 							// Initial scanning in y
			T xe; 							// final scanning position in x
			T ye; 							// final scanning position in y
			T dRx;
			T dRy;

			Vector<T, e_host> x;
			Vector<T, e_host> y;
			Vector<T, e_host> r;

			size_type size() const
			{
				return x.size();
			}

			Scanning(): type(eST_Line), grid_type(eGT_Regular), pbc(false), ns(1), 
				nx(0), dRx(0), dRy(0), ny(0), x0(0), y0(0), xe(0), ye(0) {};

			template<class TScanning> 
			void assign(TScanning &scanning)
			{
				type = scanning.type;
				grid_type = scanning.grid_type;
				pbc = scanning.pbc; 
				ns = scanning.ns;
				nx = scanning.nx;
				ny = scanning.ny;
				x0 = scanning.x0;
				y0 = scanning.y0;
				xe = scanning.xe;
				ye = scanning.ye;
				dRx = scanning.dRx;
				dRy = scanning.dRy;

				x = scanning.x;
				y = scanning.y;
				r = scanning.r;
			}

			template<class TScanning> 
			Scanning<T>& operator=(TScanning &scanning)
			{
				assign(scanning);
				return *this; 
			}

			void set_default()
			{
				type = eST_Line;
				grid_type = eGT_Regular;
				pbc = false;
				ns = 1;
				x0 = y0 = 0;
				xe = ye = 0;
			}

			int nxy() const { return nx*ny; }

			T Rx(const int &ix) const 
			{
				T x = 0;
				switch (grid_type)
				{
					case eGT_Regular:
					{
						x = x0 + ix*dRx;
					}
					break;
					case eGT_Quadratic:
					{
						x = x0 + pow(ix*dRx, 2);
					}
					break;
				}
				return x; 
			}

			T Ry(const int &iy) const
			{ 
				T y = 0;
				switch (grid_type)
				{
					case eGT_Regular:
					{
						y = y0 + iy*dRy;
					}
					break;
					case eGT_Quadratic:
					{
						y = y0 + pow(iy*dRy, 2);
					}
					break;
				}
				return y; 
			}

			void set_grid()
			{
				if(ns <= 0)
				{
					ns = nx = ny = 0;
					x.clear();
					y.clear();
					r.clear();
					return;
				}

				nx = ny = ns;
				if(is_line())
				{
					T xu = xe-x0;
					T yu = ye-y0;
					T ds = sqrt(yu*yu+xu*xu);
					T theta = atan2(yu, xu);
					T cos_theta = cos(theta);
					cos_theta = (isZero(cos_theta))?0:cos_theta;
					T sin_theta = sin(theta);
					theta = (isZero(theta))?0:theta;

					switch (grid_type)
					{
						case eGT_Regular:
						{
							dRx = ds*cos_theta/((pbc)?ns:(ns-1));
							dRy = ds*sin_theta/((pbc)?ns:(ns-1));
						}
						break;
						case eGT_Quadratic:
						{
							dRx = sqrt(ds*cos_theta)/((pbc)?ns:(ns-1));
							dRy = sqrt(ds*sin_theta)/((pbc)?ns:(ns-1));
						}
						break;
					}

					x.resize(ns);
					y.resize(ns);
					r.resize(ns);

					for(auto i = 0; i < ns; i++)
					{
						x[i] = Rx(i);
						y[i] = Ry(i);
						r[i] = sqrt(pow(x[i]-x0, 2)+pow(y[i]-y0, 2));
					}
				}
				else
				{
					T xu = xe-x0;
					T yu = ye-y0;
					if(fabs(xu)>fabs(yu))
					{
						dRx = xu/((pbc)?ns:(ns-1));
						dRy = std::copysign(dRx, yu);
						ny = int(floor(yu/dRy+Epsilon<T>::rel+0.5));
						ny += (pbc)?0:1;
					}
					else
					{
						dRy = yu/((pbc)?ns:(ns-1));
						dRx = std::copysign(dRy, xu);
						nx = int(floor(xu/dRx+Epsilon<T>::rel+0.5));
						nx += (pbc)?0:1;
					}

					x.resize(nxy());
					y.resize(nxy());

					for(auto ix = 0; ix<nx; ix++)
					{
						for(auto iy = 0; iy<ny; iy++)
						{
							x[ix*ny+iy] = Rx(ix);
							y[ix*ny+iy] = Ry(iy);
						}
					}
				}

				x.shrink_to_fit();
				y.shrink_to_fit();
				r.shrink_to_fit();
			}

			bool is_line() const
			{
				return type == eST_Line;
			}

			bool is_area() const
			{
				return type == eST_Area;
			}

		private:
	};

	/*************************Radial Schrodinger equation**********************/
	template<class T>
	struct In_Rad_Schr
	{
		T E_0; 					// Acceleration Voltage
		ePotential_Type potential_type; 			// Parameterization type
		int n; 					// Principal quantum number
		int nr; 				// Number of gridBT points
		int natomsM; 			// Number of atoms
		T *atomsM; 			// atoms
	};

	/*****************************e_device properties**************************/
	struct Device_Properties
	{
		int id;
		std::string name;
		int compute_capability;
		double total_memory_size;		// Mb
		double free_memory_size;		// Mb

		Device_Properties(): id(0), name(""), compute_capability(0), 
			total_memory_size(0), free_memory_size(0){}
	};

	/*****************************e_device properties**************************/
	struct Host_Properties
	{
		int nprocessors;
		int nthreads;
		double total_memory_size;		// Mb
		double free_memory_size;		// Mb

		Host_Properties(): nprocessors(0), nthreads(0), total_memory_size(0), 
			free_memory_size(0){}
	};

} // namespace mt

#endif