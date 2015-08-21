/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
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
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TYPES_H
#define TYPES_H

#ifndef DEVICE_CALLABLE
	#ifdef __CUDACC__
		#define DEVICE_CALLABLE __host__ __device__
		#define FORCEINLINE __forceinline__
	#else
		#define DEVICE_CALLABLE
		#define FORCEINLINE inline
	#endif
#endif

#include <cfloat>
#include <type_traits>
#include <stdio.h>
#include <string>
#include <vector>
#include <thread>
#include <mutex>

#include "math.cuh"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/detail/raw_pointer_cast.h>
#include <thrust/tuple.h>

using std::vector;
using thrust::device_vector;
using thrust::host_vector;
using thrust::raw_pointer_cast;

namespace multem
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

	const double c_mrad_2_rad = 1.0e-03; 					// mrad-->rad
	const double c_deg_2_rad = 0.01745329251994329576924;	// degrees-->rad
	const double c_mm_2_Ags = 1.0e+07; 						// mm-->Angstrom
	const double c_meV_2_keV = 1e-03; 						// ev-->keV

	const int c_cSynCPU = 5;

	const int c_nAtomsTypes = 103;
	const int c_nqz = 128;
	const int c_nR = 128;

	const int c_thrnxny = 16;
	const int c_thrnxy = 256;
	const double c_Vrl = 0.015;

	#define IsFS(i,nh)	((i < nh)?i:i-2*nh)
	#define IsRS(i,nh)	((i < nh)?i+nh:i-nh)

	const int cSizeofI = sizeof(int);
	const int cSizeofRD = sizeof(double);
	const int cSizeofRF = sizeof(float);
	const int cSizeofCD = 2*cSizeofRD;

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
	enum eMicroscope_Effect
	{
		eME_Partial_Coherent = 1, eME_Transmission_Cross_Coefficient = 2, eME_Full_Calculation = 3, eME_none = 4
	};

	/******************************Spatial and temporal***************************/
	enum eSpatial_Temporal_Effect
	{
		eSTE_Spatial_Temporal = 1, eSTE_Temporal = 2, eSTE_Spatial = 3, eSTE_none = 4
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

	/*****************************simulation type*********************************/
	enum eSimulation_Type
	{
		eST_STEM=11, eST_ISTEM=12, 
		eST_CBED=21, eST_CBEI=22, 
		eST_ED=31, eST_HRTEM=32, 
		eST_PED=41, eST_HCI=42, 
		eST_EWFS=51, eST_EWRS=52,
		eST_EELS=61, eST_EFTEM=62,
		eST_ProbeFS=71, eST_ProbeRS=72,
		eST_PPFS=81, eST_PPRS=82,			//projected potential
		eST_TFFS=91, eST_TFRS=92,			//transmission function
		eST_EWSFS=101, eST_EWSRS=102,		//exit wave single configuration
		eST_PropFS=111, eST_PropRS=112		//propagate
	};

	/********************************Potential model****************************/
	enum ePhonon_Model
	{
		ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
	};

	/******************Electron specimen interaction model**********************/
	enum eElec_Spec_Int_Model
	{
		eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
	};

	/**************************Potential Slicing Type***************************/
	enum ePotential_Slicing
	{
		ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
	};

	/***************************Potential parameterization************************/
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
		eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User = 4
	};

	/*******************************Surface Type********************************/
	enum eSurface_Type
	{
		eST_Bottom = 1, eST_Top = 2
	};

	/******************************thickness Type*******************************/
	enum eThickness_Type
	{
		eTT_Whole_Specimen = 1, eTT_Through_Thickness = 2, eTT_Through_Slices = 3
	};

	/***************************Input wave function*****************************/
	enum eInput_Wave_Type
	{
		eIWT_Automatic = 1, eIWT_User_Define = 2
	};

	/**********************************Beam type********************************/
	enum eBeam_Type
	{
		eBT_Plane_Wave = 1, eBT_Convergent = 2, eBT_User_Define = 3
	};

	/******************************Scanning Type********************************/
	enum eScanning_Type
	{
		eST_Line = 1, eST_Area = 2
	};


	/**********************************Channelling type********************************/
	enum eChannelling_Type
	{
		eCT_Single_Channelling = 1, eCT_Double_Channelling = 2, eCT_Double_Channelling_FOMS = 3, eCT_Double_Channelling_SOMS = 4
	};

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

	template<class T>
	using Pos_2d = thrust::tuple<T, T>;

	template<class T>
	using Pos_3d = thrust::tuple<T, T, T>;

	/*******************forward declarations********************/
	template<class T>
	struct is_fundamental;

	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	T get_lambda(const T &E_0);

	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	T get_sigma(const T &E_0);

	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
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

		int size;
		T *V;

		rVector(): size(0), V(nullptr){}

		template<class TVector>
		rVector(TVector &vector)
		{
			size = vector.size();
			V = raw_pointer_cast(vector.data());
		}

		DEVICE_CALLABLE FORCEINLINE
		T& operator[](const int i){ return V[i]; }

		DEVICE_CALLABLE FORCEINLINE
		const T& operator[](const int i) const { return V[i]; }
	};

	template<class T>
	DEVICE_CALLABLE FORCEINLINE
	double sizeMb(const int &n)
	{
		return static_cast<double>(n*sizeof(T)/1048576.0);
	}

	//static member function are not support for the cuda compiler
	template<class T>
	DEVICE_CALLABLE FORCEINLINE
	bool isEqual(const T &a, const T &b);

	template<>
	DEVICE_CALLABLE FORCEINLINE
	bool isEqual<int>(const int &a, const int &b)
	{
		return a == b;
	}

	template<>
	DEVICE_CALLABLE FORCEINLINE
	bool isEqual<float>(const float &a, const float &b)
	{
		const float eps_abs = 1e-5f;
		const float eps_rel = 1e-4f;

		// Check if the numbers are really close -- needed when comparing numbers near zero.
		float diff = abs(a - b);
		if (diff <= eps_abs)
			return true;
 
		// Otherwise fall back to Knuth's algorithm
		return diff <= ((abs(a)<abs(b)?abs(b):abs(a))*eps_rel);
	}

	template<>
	DEVICE_CALLABLE FORCEINLINE
	bool isEqual<double>(const double &a, const double &b)
	{
		const double eps_abs = 1e-13;
		const double eps_rel = 1e-8;

		// Check if the numbers are really close -- needed when comparing numbers near zero.
		double diff = abs(a - b);
		if (diff <= eps_abs)
			return true;
 
		// Otherwise fall back to Knuth's algorithm
		return diff <= ((abs(a)<abs(b)?abs(b):abs(a))*eps_rel);
	}

	template<class T>
	DEVICE_CALLABLE FORCEINLINE
	bool isZero(const T &x)
	{
		return isEqual<T>(x, 0);
	}

	template<class T, class U>
	DEVICE_CALLABLE FORCEINLINE
	bool isZero(const T &x, const U &y)
	{
		return isEqual<T>(x, 0) && isEqual<U>(y, 0);
	}

	template<class T>
	DEVICE_CALLABLE FORCEINLINE
	bool nonZero(const T &x)
	{
		return !isEqual<T>(x, 0);
	}

	template<class T, class U>
	DEVICE_CALLABLE FORCEINLINE
	bool nonZero(const T &x, const U &y)
	{
		return !(isEqual<T>(x, 0) || isEqual<U>(y, 0));
	}

	template<class T, class U>
	DEVICE_CALLABLE FORCEINLINE
	T Div(const T &x, const U &y)
	{
		return (isEqual<U>(y, 0))?0:static_cast<T>(x)/static_cast<T>(y);
	}

	template<class T>
	DEVICE_CALLABLE FORCEINLINE
	bool Check_Bound(const T &x, const T &x_min, const T &x_max)
	{
		return (x_min <= x)&&(x <= x_max);
	}

	template<class T, class U>
	DEVICE_CALLABLE FORCEINLINE
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
		Vector<int, dev> iatom_e; 			// Last Atom index
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
			return (islice<size())?(0.5*abs(z_e[islice]+z_0[islice])):0.0;
		}

		T dz_m(const int &islice_0, const int &islice_e)
		{
			return abs(z_m(islice_e) - z_m(islice_0));
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

		rQ1(): size(0), x(nullptr), w(nullptr){}

		template<class TQ1> 
		rQ1<T>& operator=(TQ1 &q1)
		{
			size = q1.size();
			x = raw_pointer_cast(q1.x.data());
			w = raw_pointer_cast(q1.w.data());
			return *this; 
		}

		template<class TQ1>
		rQ1(TQ1 &q1)
		{
			*this = q1;
		}

		int size;
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

		rQ2(): size(0), x(nullptr), y(nullptr), w(nullptr){}

		template<class TQ2> 
		rQ2<T>& operator=(TQ2 &q2)
		{
			size = q2.size();
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

		int size;
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

		rPP_Coef(): size(0), cl(nullptr), cnl(nullptr){}

		template<class TPP_Coef> 
		rPP_Coef<T>& operator=(TPP_Coef &rhs)
		{ 
			size = rhs.size();
			cl = raw_pointer_cast(rhs.cl.data());
			cnl = raw_pointer_cast(rhs.cnl.data());
			return *this; 
		}

		template<class TPP_Coef>
		rPP_Coef(TPP_Coef &pp_coef)
		{
			*this = pp_coef;
		}

		int size;
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

		rCI_Coef(): size(0), c0(nullptr), c1(nullptr), c2(nullptr), c3(nullptr){}

		template<class TCI_Coef> 
		rCI_Coef<T>& operator=(TCI_Coef &ci_coef)
		{
			size = ci_coef.size();
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

		int size;
		T *c0;
		T *c1;
		T *c2;
		T *c3;
	};

	/**********************STEM Detector**********************/
	template<class T>
	struct Det_Cir
	{
		using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = e_host;

		inline
		void set_input_data(value_type E_0)
		{
			lambda = get_lambda(E_0);
		}		
		
		size_type size() const
		{
			return ang_inner.size();
		}

		void resize(const size_type &new_size, const value_type &value = value_type())
		{
			ang_inner.resize(new_size, value);
			ang_outer.resize(new_size, value);
		}

		template<class TDet_Cir> 
		void assign(TDet_Cir &det_cir)
		{
			ang_inner.assign(det_cir.ang_inner.begin(), det_cir.ang_inner.end());
			ang_outer.assign(det_cir.ang_outer.begin(), det_cir.ang_outer.end());
		}

		value_type g_inner(const int & idx) const
		{
			return ang_inner[idx]/lambda;
		}

		value_type g_outer(const int & idx) const
		{
			return ang_outer[idx]/lambda;
		}

		Vector<T, e_host> ang_inner; // Inner aperture (rad)
		Vector<T, e_host> ang_outer; // Outer aperture (rad)
		value_type lambda;	 // lambda
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
		Range():ix_0(0), ix_e(0), iy_0(0), iy_e(0), ixy_0(0), ixy_e(0){}

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

		T lx; 					// Box size in x direction(Angstroms)
		T ly; 					// Box size in y direction(Angstroms)
		T dz; 					// slice thicknes

		bool bwl; 				// Band-width limit
		bool pbc_xy; 			// Peridic boundary contions

		T dRx; 					// x-sampling in real space
		T dRy; 					// y-sampling in real space

		T dgx; 					// x-sampling in reciprocal space
		T dgy; 					// y-sampling in reciprocal space

		T gl_max; 				// Maximun limited frequency 
		T gl2_max; 				// Squared of the maximun limited frequency
		T alpha;

		inline
		Grid(): nx(0), ny(0), nxh(0), nyh(0), inxy(0), 
			lx(0), ly(0), dz(0), pbc_xy(true), bwl(true), 
			dRx(0), dRy(0), dgx(0), dgy(0), gl_max(0), gl2_max(0){}

		inline
		void set_input_data(int nx_i, int ny_i, T lx_i, T ly_i, T dz_i, bool BWL_i, bool PBC_xy_i)
		{
			nx = nx_i;
			ny = ny_i;
			nxh = nx/2;
			nyh = ny/2;
			inxy = multem::Div(1.0, nxy());
			lx = lx_i;
			ly = ly_i;
			dz = dz_i;
			bwl = BWL_i;
			pbc_xy = PBC_xy_i;
			dRx = multem::Div(lx, nx);
			dRy = multem::Div(ly, ny);
			dgx = multem::Div(1.0, lx);
			dgy = multem::Div(1.0, ly);
			gl_max = 2.0*g_max()/3.0;
			gl2_max = pow(gl_max, 2);

			T dg0 = 0.02, fg0 = 0.99;
			alpha = log(1.0/fg0-1.0)/(pow(gl_max-dg0, 2)-gl2_max);
		}
	
		DEVICE_CALLABLE FORCEINLINE
		int nxy() const { return nx*ny; }

		DEVICE_CALLABLE FORCEINLINE
		int nx_ny_min() const { return min(nx, ny); }

		DEVICE_CALLABLE FORCEINLINE
		int nx_ny_max() const { return max(nx, ny); }

		/*********************************************************/
		DEVICE_CALLABLE FORCEINLINE
		int lx_ly_min() const { return min(lx, ly); }

		DEVICE_CALLABLE FORCEINLINE
		int lx_ly_max() const { return max(lx, ly); }

		DEVICE_CALLABLE FORCEINLINE
		T lxh() const { return 0.5*lx; }

		DEVICE_CALLABLE FORCEINLINE
		T lyh() const { return 0.5*ly; }

		/*********************************************************/
		// Maximun frequency
		DEVICE_CALLABLE FORCEINLINE
		T g_max() const { return fmin(static_cast<T>(nxh)*dgx, static_cast<T>(nyh)*dgy); }

		// Squared of the maximun frequency
		DEVICE_CALLABLE FORCEINLINE
		T g2_max() const { return pow(g_max(), 2); }

		DEVICE_CALLABLE FORCEINLINE
		T dR_min() const { return fmin(dRx, dRy); }

		DEVICE_CALLABLE FORCEINLINE
		T dg_min() const { return fmin(dgx, dgy); }

		DEVICE_CALLABLE FORCEINLINE
		int nx_dRx(const T &lx) const { return static_cast<int>(ceil(lx/dRx)); }

		DEVICE_CALLABLE FORCEINLINE
		int ny_dRy(const T &ly) const { return static_cast<int>(ceil(ly/dRy)); }

		/*********************************************************/
		DEVICE_CALLABLE FORCEINLINE
		int igx(const int &ix) const { return ix-nxh; }

		DEVICE_CALLABLE FORCEINLINE
		int igy(const int &iy) const { return iy-nyh; }

		DEVICE_CALLABLE FORCEINLINE
		T gx(const int &ix) const { return igx(ix)*dgx; }

		DEVICE_CALLABLE FORCEINLINE
		T gy(const int &iy) const { return igy(iy)*dgy; }

		DEVICE_CALLABLE FORCEINLINE
		T g2(const int &ix, const int &iy) const 
		{ 
			T gxi = gx(ix);
			T gyi = gy(iy);
			return (gxi*gxi + gyi*gyi);
		}

		DEVICE_CALLABLE FORCEINLINE
		T g(const int &ix, const int &iy) const { return sqrt(g2(ix, iy)); }


		DEVICE_CALLABLE FORCEINLINE
		T gx(const int &ix, const T &x) const { return gx(ix)-x; }

		DEVICE_CALLABLE FORCEINLINE
		T gy(const int &iy, const T &y) const { return gy(iy)-y; }

		DEVICE_CALLABLE FORCEINLINE
		T g2(const int &ix, const int &iy, const T &x, const T &y) const 
		{ 
			T gxi = gx(ix, x);
			T gyi = gy(iy, y);
			return (gxi*gxi + gyi*gyi);
		}

		DEVICE_CALLABLE FORCEINLINE
		T g(const int &ix, const int &iy, const T &x, const T &y) const { return sqrt(g2(ix, iy, x, y)); }

		DEVICE_CALLABLE FORCEINLINE
		int iRx(const int &ix) const { return ix; }

		DEVICE_CALLABLE FORCEINLINE
		int iRy(const int &iy) const { return iy; }

		DEVICE_CALLABLE FORCEINLINE
		T Rx(const int &ix) const { return iRx(ix)*dRx; }

		DEVICE_CALLABLE FORCEINLINE
		T Ry(const int &iy) const { return iRy(iy)*dRy; }

		DEVICE_CALLABLE FORCEINLINE
		T R2(const int &ix, const int &iy) const 
		{ 
			T Rxi = Rx(ix);
			T Ryi = Ry(iy);
			return (Rxi*Rxi + Ryi*Ryi);
		}

		DEVICE_CALLABLE FORCEINLINE
		T R(const int &ix, const int &iy) const { return sqrt(R2(ix, iy)); }

		DEVICE_CALLABLE FORCEINLINE
		T Rx(const int &ix, const T &x) const { return Rx(ix)-x; }

		DEVICE_CALLABLE FORCEINLINE
		T Ry(const int &iy, const T &y) const { return Ry(iy)-y; }

		DEVICE_CALLABLE FORCEINLINE
		T R2(const int &ix, const int &iy, const T &x, const T &y) const 
		{ 
			T Rxi = Rx(ix, x);
			T Ryi = Ry(iy, y);
			return (Rxi*Rxi + Ryi*Ryi);
		}

		DEVICE_CALLABLE FORCEINLINE
		T R(const int &ix, const int &iy, const T &x, const T &y) const { return sqrt(R2(ix, iy, x, y)); }

		DEVICE_CALLABLE FORCEINLINE
		T bwl_factor(const int &ix, const int &iy) const 
		{ 
			if(bwl)
			{
				return inxy/(1.0+exp(alpha*(g2(ix, iy)-gl2_max)));
			}
			else
			{
				return inxy;
			}
		}

		/*********************************************************/
		DEVICE_CALLABLE FORCEINLINE
		int igx_shift(const int &ix) const { return (ix<nxh)?ix:ix-nx; }

		DEVICE_CALLABLE FORCEINLINE
		int igy_shift(const int &iy) const { return (iy<nyh)?iy:iy-ny; }

		DEVICE_CALLABLE FORCEINLINE
		T gx_shift(const int &ix) const { return static_cast<T>(igx_shift(ix))*dgx; }

		DEVICE_CALLABLE FORCEINLINE
		T gy_shift(const int &iy) const { return static_cast<T>(igy_shift(iy))*dgy; }

		DEVICE_CALLABLE FORCEINLINE
		T g2_shift(const int &ix, const int &iy) const 
		{ 
			T gxi = gx_shift(ix);
			T gyi = gy_shift(iy);
			return gxi*gxi + gyi*gyi;
		}

		DEVICE_CALLABLE FORCEINLINE
		T g_shift(const int &ix, const int &iy) const { return sqrt(g2_shift(ix, iy)); }

		DEVICE_CALLABLE FORCEINLINE
		int iRx_shift(const int &ix) const { return (ix<nxh)?ix+nxh:ix-nxh; }

		DEVICE_CALLABLE FORCEINLINE
		int iRy_shift(const int &iy) const { return (iy<nyh)?iy+nyh:iy-nyh; }

		DEVICE_CALLABLE FORCEINLINE
		T Rx_shift(const int &ix) const { return static_cast<T>(iRx_shift(ix))*dRx; }

		DEVICE_CALLABLE FORCEINLINE
		T Ry_shift(const int &iy) const { return static_cast<T>(iRy_shift(iy))*dRy; }

		DEVICE_CALLABLE FORCEINLINE
		T R2_shift(const int &ix, const int &iy) const 
		{ 
			T Rxi = Rx_shift(ix);
			T Ryi = Ry_shift(iy);
			return Rxi*Rxi + Ryi*Ryi;
		}

		DEVICE_CALLABLE FORCEINLINE
		T R_shift(const int &ix, const int &iy) const { return sqrt(R2_shift(ix, iy)); }

		DEVICE_CALLABLE FORCEINLINE
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
		DEVICE_CALLABLE FORCEINLINE
		int ind_row(const int &ix, const int &iy) const { return iy*nx+ix; }

		DEVICE_CALLABLE FORCEINLINE
		int ind_col(const int &ix, const int &iy) const { return ix*ny+iy; }

		DEVICE_CALLABLE FORCEINLINE
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

		int m; 		// vortex momentum
		T f; 		// defocus
		T Cs3; 	// spherical aberration(Angstrom)
		T Cs5; 	// spherical aberration(Angstrom)
		T mfa2; 	// magnitude 2-fold astigmatism
		T afa2; 	// angle 2-fold astigmatism(rad)
		T mfa3; 	// magnitude 3-fold astigmatism
		T afa3; 	// angle 3-fold astigmatism(rad)
		T aobjl; 	// lower objective aperture(rad);
		T aobju; 	// upper objective aperture(rad);

		T sf; 		// defocus spread
		int nsf; 	// Number of defocus sampling points

		T beta; 	// semi-convergence angle
		int nbeta; 	// Number of semi-convergence angle sampling points

		T gamma; 	// Relativistic factor
		T lambda; 	// wavelength(Angstrom)
		T lambda2; 	// wavelength(Angstrom)^2

		T cf; 		// pi*f*lambda
		T cCs3; 	// -0.5*pi*Cs3*lambda^3
		T cCs5; 	// -pi*Cs5*lambda^5/3
		T cmfa2; 	// -pi*lambda
		T cmfa3; 	// -2*pi*lambda^3/3
		T g2_min; 	// aobjl/lambda
		T g2_max; 	// aobju/lambda

		T sggs; 	// Standard deviation
		int ngxs; 	// Number of source sampling points x
		int ngys; 	// Number of source sampling points y
		T dgxs; 	// source sampling size;
		T dgys; 	// source sampling size;
		T g2_maxs; 	// q maximum square;

		Lens(): gamma(0), lambda(0), m(0), f(0), Cs3(0), Cs5(0), 
				mfa2(0),afa2(0), mfa3(0), afa3(0), aobjl(0),
				aobju(0), sf(0), nsf(0), beta(0), nbeta(0),
				lambda2(0), cf(0), cCs3(0),	cCs5(0), cmfa2(0),
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
			g2_min = (isZero(aobjl)||(aobjl<0))?0:pow(aobjl/lambda, 2);
			g2_max = (isZero(aobju)||(aobju<0))?grid.g2_max():pow(aobju/lambda, 2);

			T g0s = beta/lambda;
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
			cf = (isZero(f))?0:c_Pi*f_i*lambda;
		}

		inline
		T prop_factor(const T &z) const { return -c_Pi*lambda*z; }

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
		m_selection(0), collection_angle(0), channelling_type(eCT_Double_Channelling), factor(0), Z(0), x(0), y(0), occ(0){}

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

		bool is_Double_Channelling() const
		{
			return channelling_type == eCT_Double_Channelling;
		}

		bool is_Double_Channelling_FOMS() const
		{
			return channelling_type == eCT_Double_Channelling_FOMS;
		}

		bool is_Double_Channelling_SOMS() const
		{
			return channelling_type == eCT_Double_Channelling_SOMS;
		}

		bool is_Double_Channelling_FOMS_SOMS() const
		{
			return is_Double_Channelling_FOMS() || is_Double_Channelling_SOMS();
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

	/*****************************Atomic type****************************/
	template<class T, eDevice dev>
	struct Atom_Type
	{
		using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = dev;

		Atom_Type(): Z(0), A(0), rn_e(0), rn_c(0), ra_e(0), ra_c(0), 
						R_min(0), R_max(0), R_min2(0), R_max2(0){}

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
			R_min = atom_type.R_min;
			R_max = atom_type.R_max;
			R_min2 = atom_type.R_min2;
			R_max2 = atom_type.R_max2;

			feg.assign(atom_type.feg);
			fxg.assign(atom_type.fxg);
			Pr.assign(atom_type.Pr);
			Vr.assign(atom_type.Vr);
			VR.assign(atom_type.VR);

			R.assign(atom_type.R.begin(), atom_type.R.end());
			R2.assign(atom_type.R2.begin(), atom_type.R2.end());
			ciVR.assign(atom_type.ciVR);
		}

		int Z; 						// Atomic number
		T m; 						// Atomic mass
		int A; 						// Mass number
		T rn_e; 					// Experimental Nuclear radius
		T rn_c; 					// Calculated Nuclear radius
		T ra_e; 					// Experimental atomic radius
		T ra_c; 					// Calculated atomic radius
		T R_min; 					// Minimum interaction radius
		T R_max; 					// Maximum interaction radius
		T R_min2; 					// Minimum interaction radius squared
		T R_max2; 					// Maximum interaction radius squared

		PP_Coef<T, dev> feg; 		// Electron scattering factor coefficients
		PP_Coef<T, dev> fxg; 		// X-ray scattering factor coefficients
		PP_Coef<T, dev> Pr; 		// Potential coefficients
		PP_Coef<T, dev> Vr; 		// Potential coefficients
		PP_Coef<T, dev> VR; 		// Projected potential coefficients

		Vector<T, dev> R; 			// R gridBT
		Vector<T, dev> R2; 			// R2 gridBT
		CI_Coef<T, dev> ciVR; 		// Look up table - Projected potential coefficients

	};

	/******************************Atom_Vp*******************************/
	template<class T>
	struct Atom_Vp
	{
		public:
			using value_type = T;

			T x;
			T y;
			T z0h;
			T zeh;
			bool split;
			T occ;
			T R_min2;
			T R_max2;
			T *R2;
			T *cl;
			T *cnl;
			T *c0;
			T *c1;
			T *c2;
			T *c3;
			int ix0;
			int ixn;
			int iy0;
			int iyn;

			int *iv;
			T *v;

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

		private:
			inline
			void get_bn(const T &R, const int &nR, const T &dR, const T &R_max, const bool &pbc, int &iR0, int &iRn)
			{
				int iR_0 = static_cast<int>(floor((R-R_max)/dR));
				int iR_e = static_cast<int>(ceil((R+R_max)/dR));

				if(!pbc)
				{
					auto set_Bound = [](const int &i, const int &n)->int{ return (i<0)?0:((i>=n)?n-1:i); };
					iR_0 = set_Bound(iR_0, nR);
					iR_e = set_Bound(iR_e, nR);
				}

				iR0 = iR_0;
				iRn = (iR_0 == iR_e)?0:iR_e-iR_0+1;
			};
	};

	/******************************Scanning******************************/
	template<class T>
	struct Scanning
	{
		using value_type = T;
		using size_type = std::size_t;

		eScanning_Type type;	// 1: Line, 2: Area, 
		int ns; 						// Sampling points
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

		size_type size() const
		{
			return x.size();
		}

		Scanning():type(eST_Line), ns(1), nx(0), dRx(0), dRy(0),
			ny(0), x0(0), y0(0), xe(0), ye(0){};

		void set_default()
		{
			type = eST_Line;
			ns = 1;
			x0 = y0 = 0;
			xe = ye = 0;
		}

		void get_grid_dim()
		{
			if(ns <= 0)
			{
				nx = ny = 0;
				return;
			}

			nx = ny = ns;
			if(type == eST_Area)
			{
				T lx = fabs(xe-x0);
				T ly = fabs(ye-y0);
				if(lx>ly)
				{
					ny = static_cast<int>(ceil(ns*ly/lx));
				}
				else
				{
					nx = static_cast<int>(ceil(ns*lx/ly));
				}
			}
		}

		int nxy() const { return nx*ny; }

		T Rx(const int &ix) const { return x0 + ix*dRx; }

		T Ry(const int &iy) const { return y0 + iy*dRy; }

		void set_grid()
		{
			get_grid_dim();

			if(ns <= 0)
			{				
				return;
			}

			x.resize(nx);
			x.shrink_to_fit();

			y.resize(ny);
			y.shrink_to_fit();

			if(type == eST_Line)
			{
				T xu = xe-x0;
				T yu = ye-y0;
				T ds = sqrt(yu*yu+xu*xu)/static_cast<T>(ns);
				T theta = atan2(yu, xu);
				T cos_theta = cos(theta);
				T sin_theta = sin(theta);

				dRx = ds*cos_theta;
				dRy = ds*sin_theta;

				x.resize(ns);
				y.resize(ns);

				for(auto i=0; i < ns; i++)
				{
					x[i] = Rx(i);
					y[i] = Ry(i);
				}
			}
			else
			{
				T xu = xe-x0;
				T yu = ye-y0;

				dRx = xu/static_cast<T>(nx);
				dRy = yu/static_cast<T>(ny);

				x.resize(nxy());
				y.resize(nxy());

				for(auto ix=0; ix<nx; ix++)
				{
					for(auto iy=0; iy<ny; iy++)
					{
						x[ix*ny+iy] = Rx(ix);
						y[ix*ny+iy] = Ry(iy);
					}
				}
			}

			x.shrink_to_fit();
			y.shrink_to_fit();
		}

		bool is_line() const
		{
			return type == eST_Line;
		}

		bool is_area() const
		{
			return type == eST_Area;
		}
	};

	/********************Radial Schrodinger equation*********************/
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

	/********************************hrtem******************************/
	template<class T>
	struct HRTEM
	{
		int xx;
		HRTEM():xx(0){};
	};

	/*****************************CBED/CBEI*****************************/
	template<class T>
	struct CBE_FR
	{
		CBE_FR(): x0(0), y0(0){};

		T x0;
		T y0;
	};

	/***************************Exit Wave FS/RS*************************/
	template<class T>
	struct EW_FR
	{
		EW_FR():convergent_beam(false), x0(0), y0(0){};

		bool convergent_beam;
		T x0;
		T y0;
	};

	/***************************e_device properties*************************/
	struct Device_Properties
	{
		int id;
		std::string name;
		int compute_capability;
		double total_memory_size;		// Mb
		double free_memory_size;		// Mb

		Device_Properties():id(0), name(""), compute_capability(0), 
			total_memory_size(0), free_memory_size(0){}
	};

	/***************************e_device properties*************************/
	struct Host_Properties
	{
		int nprocessors;
		int nthreads;
		double total_memory_size;		// Mb
		double free_memory_size;		// Mb

		Host_Properties(): nprocessors(0), nthreads(0), total_memory_size(0), 
			free_memory_size(0){}
	};

} // namespace multem

#endif