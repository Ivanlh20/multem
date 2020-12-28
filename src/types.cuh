/*
 * This file is part of MULTEM.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
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

//#ifdef __CUDACC__
//	#pragma message("Cuda TYPES_H")
//#else
//	#pragma message("nonCuda TYPES_H")
//#endif

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
#include <thrust/reduce.h>

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

	const double c_hwhm_2_sigma = 0.84932180028801907; 		// hwhm to sigma 1/(sqrt(2*log(2)))
	const double c_fwhm_2_sigma = 0.42466090014400953; 		// fwhm to sigma 1/(2*sqrt(2*log(2)))
	const double c_iehwgd_2_sigma = 0.70710678118654746; 	// iehwgd to sigma 1/sqrt(2)

	const double c_mrad_2_rad = 1.0e-03; 					// mrad-->rad
	const double c_deg_2_rad = 0.01745329251994329576924;	// degrees-->rad
	const double c_mm_2_Angs = 1.0e+07; 					// mm-->Angstrom
	const double c_eV_2_keV = 1e-03; 						// ev-->keV

	const int c_cSynCPU = 5;

	const int c_nAtomsTypes = 103;
	const int c_nAtomsIons = 15;
	const int c_nqz = 128;
	const int c_nR = 128;

	const int c_thrnxny = 16;
	const int c_thrnxy = 256;
	const double c_Vrl = 0.015;

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

	/*************************************data type******************************/
	enum eData_Type
	{
		eDT_float = 1, eDT_double = 2, eDT_cfloat = 3, eDT_cdouble = 4
	};

	/*****************************Show Data Type**********************************/
	enum eShow_CData
	{
		eSCD_CReal = 1, eSCD_CImag = 2, eSCD_CMod = 3, eSCD_CPhase = 4
	};

	/*********************************Operation mode******************************/
	enum eOperation_Mode
	{
		eOM_Normal = 1, eOM_Advanced = 2
	};

	/****************************lens variable type******************************/
	enum eLens_Var_Type
	{
		eLVT_off = 0, eLVT_m = 1, eLVT_f = 2, eLVT_Cs3 = 3, eLVT_Cs5 = 4,
		eLVT_mfa2 = 5, eLVT_afa2 = 6, eLVT_mfa3 = 7, eLVT_afa3 = 8,
		eLVT_inner_aper_ang = 9, eLVT_outer_aper_ang = 10
	};

	/*****************************simulation type********************************/
	enum eTEM_Sim_Type
	{
		eTEMST_STEM = 11, eTEMST_ISTEM = 12,
		eTEMST_CBED = 21, eTEMST_CBEI = 22,
		eTEMST_ED = 31, eTEMST_HRTEM = 32,
		eTEMST_PED = 41, eTEMST_HCTEM = 42,
		eTEMST_EWFS = 51, eTEMST_EWRS = 52,
		eTEMST_EELS = 61, eTEMST_EFTEM = 62,
		eTEMST_IWFS = 71, eTEMST_IWRS = 72,
		eTEMST_PPFS = 81, eTEMST_PPRS = 82, 			// projected potential
		eTEMST_TFFS = 91, eTEMST_TFRS = 92, 			// transmission function
		eTEMST_PropFS = 101, eTEMST_PropRS = 102		// propagate
	};

	/*************************simulation data output*****************************/
	enum eTEM_Output_Type
	{
		eTEMOT_image_tot_coh = 1, eTEMOT_image_tot = 2,
		eTEMOT_m2psi_tot_coh = 3, eTEMOT_m2psi_tot = 4,
		eTEMOT_m2psi_tot_psi_coh = 5, eTEMOT_psi_coh = 6,
		eTEMOT_psi_0 = 7, eTEMOT_V = 8, eTEMOT_trans = 9
	};

	/******************Electron specimen interaction model**********************/
	enum eElec_Spec_Int_Model
	{
		eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
	};

	/*****************************Frozen lattice model**************************/
	enum ePhonon_Model
	{
		ePM_Still_Atom = 1, ePM_Absorptive_Model = 2, ePM_Frozen_Phonon = 3
	};

	/*******************************Extract data********************************/
	enum ePhonon_Model_Output
	{
		eFMO_Total = 1, eFMO_Coherent = 2
	};

	/*********************Projected_Potential Slicing Type**********************/
	enum ePotential_Slicing
	{
		ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
	};

	/********************Projected_Potential parameterization******************/
	enum ePotential_Type
	{
		ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3,
		ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6, ePT_none = 0
	};

	/***************************Incident Wave Type******************************/
	enum eIncident_Wave_Type
	{
		eIWT_Plane_Wave = 1, eIWT_Convergent_Wave = 2, eIWT_User_Define_Wave = 3, eIWT_Auto = 4
	};

	enum eRot_Point_Type
	{
		eRPT_geometric_center = 1, eRPT_User_Define = 2
	};

	/*****************************Real or Fourier space**************************/
	enum eSpace
	{
		eS_Real = 1, eS_Reciprocal = 2
	};

	/****************************Defocus plane type*****************************/
	enum eMatch_Border
	{
		eMB_Min = 1, eMB_Max = 2, eMB_MinMax = 3
	};

	/****************************Amorphous layer Type***************************/
	enum eAmorp_Lay_Type
	{
		eALT_Top = 1, eALT_Bottom = 2, eALT_Middle= 3, eALT_none = 4
	};

	/****************************Defocus plane type*****************************/
	enum eZero_Defocus_Type
	{
		eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
	};

	/*******************************thick Type*********************************/
	enum eThick_Type
	{
		eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
	};

	/******************************Scanning Type*******************************/
	enum eScanning_Type
	{
		eST_Line = 1, eST_Area = 2
	};
	/******************************grid_2d Type*******************************/
	enum eGrid_Type
	{
		eGT_Regular = 1, eGT_Quadratic = 2
	};

	/******************************Detector type*******************************/
	enum eDetector_Type
	{
		eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3
	};

	/********************************Channelling type*************************/
	enum eChannelling_Type
	{
		eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3
	};

	/*******************************Output type*******************************/
	enum eOutput_Type
	{
		eOT_Matlab = 1, eOT_Vector = 2
	};

	/****************************Data selection type**************************/
	enum eDat_Sel_Type
	{
		eDST_Closest = 1, eDST_Less_Than = 2, eDST_Greater_Than = 3
	};

	/**************************structuring element****************************/
	enum eStr_Ele
	{
		eSE_Disk = 1, eSE_Square = 2
	};

	/******************************operation**********************************/
	enum eOP {
		eOP_N=1, eOP_T=2, eOP_C=3
	};

	/********************************Exec type********************************/
	enum eET {
		eET_Matrix=1, eET_Vector=2
	};

	template <class T>
	void get_bn(const T &R, const int &nR, const T &dR, const T &R_max, const bool &pbc, int &iR0, int &iRn);

	template <class T>
	struct r2d;

	template <typename X>
	X norm(const r2d<X>& r);

	template <class T>
	struct r3d;

	template <typename X>
	X norm(const r3d<X>& r);

	namespace host_device_detail
	{
		template <class T>
		DEVICE_CALLABLE FORCE_INLINE
		T tapering(const T &x_tap, const T &alpha, const T &x);

		template <class T>
		DEVICE_CALLABLE FORCE_INLINE
		void kh_sum(T &sum_v, T v, T &error);
	}

	template <class T>
	class Atom_Data;

	/*********************************Epsilon***********************************/
	template<class T>
	DEVICE_CALLABLE FORCE_INLINE
	T epsilon_eps(){ return 0; }

	template<>
	DEVICE_CALLABLE FORCE_INLINE
	double epsilon_eps(){ return 10.0*DBL_EPSILON; }

	template<>
	DEVICE_CALLABLE FORCE_INLINE
	float epsilon_eps(){ return 10.0*FLT_EPSILON; }

	template<class T>
	DEVICE_CALLABLE FORCE_INLINE
	T epsilon_abs(){ return 0; }

	template<>
	DEVICE_CALLABLE FORCE_INLINE
	double epsilon_abs(){ return 1e-13; }

	template<>
	DEVICE_CALLABLE FORCE_INLINE
	float epsilon_abs(){ return 1e-5f; }

	template<class T>
	DEVICE_CALLABLE FORCE_INLINE
	T epsilon_rel(){ return 0; }

	template<>
	DEVICE_CALLABLE FORCE_INLINE
	double epsilon_rel(){ return 1e-8; }

	template<>
	DEVICE_CALLABLE FORCE_INLINE
	float epsilon_rel(){ return 1e-4f; }

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

	struct Grid_BT
	{
		dim3 Blk; 	// Blocks
		dim3 Thr; 	// Threads

	};

	/*******************forward declarations********************/
	template <class T>
	struct is_fundamental;

	template <class T>
	DEVICE_CALLABLE FORCE_INLINE
	T get_lambda(const T &E_0);

	template <class T>
	DEVICE_CALLABLE FORCE_INLINE
	T get_sigma(const T &E_0);

	template <class T>
	DEVICE_CALLABLE FORCE_INLINE
	T get_gamma(const T &E_0);

	/************************vector type***********************/
	//template<class TVector>
	//eData_Type Vector_Type_to_Data_Type()
	//{
	//	if(is_float<TVector::value_type>)
	//	{
	//		return eDT_float;
	//	}
	//	else if(is_double<TVector::value_type>)
	//	{
	//		return eDT_double;
	//	}
	//	else if(is_cfloat<TVector::value_type>)
	//	{
	//		return eDT_cfloat;
	//	}
	//	else(is_cdouble<TVector::value_type>)
	//	{
	//		return eDT_cdouble;
	//	}
	//}

	/**************************vector**************************/
	template <class T, eDevice dev>
	using Vector = typename std::conditional<dev == e_host, typename std::conditional<std::is_fundamental<T>::value ||
	std::is_same<T, complex<float>>::value || std::is_same<T, complex<double>>::value, host_vector<T>, vector<T>>::type, device_vector<T>>::type;

	template <class T>
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

		rVector(const rVector<T> &vector)
		{
			m_size = vector.m_size;
			V = vector.V;
		}

		rVector(Vector<T, e_host> &vector)
		{
			m_size = vector.size();
			V = raw_pointer_cast(vector.data());
		}

		rVector(Vector<T, e_device> &vector)
		{
			m_size = vector.size();
			V = raw_pointer_cast(vector.data());
		}

		//rVector(host_vector<T> &vector)
		//{
		//	m_size = vector.size();
		//	V = raw_pointer_cast(vector.data());
		//}

		//rVector(device_vector<T> &vector)
		//{
		//	m_size = vector.size();
		//	V = raw_pointer_cast(vector.data());
		//}
		DEVICE_CALLABLE FORCE_INLINE
		T& operator[](const int i){ return V[i]; }

		DEVICE_CALLABLE FORCE_INLINE
		const T& operator[](const int i) const { return V[i]; }
	};

	template <class T>
	DEVICE_CALLABLE FORCE_INLINE
	double sizeMb(const int &n)
	{
		return static_cast<double>(n*sizeof(T)/1048576.0);
	}

	// static member function are not supported for the cuda compiler
	template <class T>
	DEVICE_CALLABLE FORCE_INLINE
	bool isEqual(const T &a, const T &b);

	template <>
	DEVICE_CALLABLE FORCE_INLINE
	bool isEqual<int>(const int &a, const int &b)
	{
		return a == b;
	}

	template <>
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

	template <>
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

	template <class T>
	DEVICE_CALLABLE FORCE_INLINE
	bool isZero(const T &x)
	{
		return isEqual<T>(x, 0);
	}

	template <class T, class U>
	DEVICE_CALLABLE FORCE_INLINE
	bool isZero(const T &x, const U &y)
	{
		return isEqual<T>(x, 0) && isEqual<U>(y, 0);
	}

	template <class T>
	DEVICE_CALLABLE FORCE_INLINE
	bool isZero(const r2d<T> &r)
	{
		return isZero<T, T>(r.x, r.y);
	}

	template <class T>
	DEVICE_CALLABLE FORCE_INLINE
	bool nonZero(const T &x)
	{
		return !isEqual<T>(x, 0);
	}

	template <class T, class U>
	DEVICE_CALLABLE FORCE_INLINE
	bool nonZero(const T &x, const U &y)
	{
		return !(isEqual<T>(x, 0) || isEqual<U>(y, 0));
	}

	template <class T>
	DEVICE_CALLABLE FORCE_INLINE
	bool nonZero(const r2d<T> &r)
	{
		return nonZero<T, T>(r.x, r.y);
	}

	template <class T, class U>
	DEVICE_CALLABLE FORCE_INLINE
	T Div(const T &x, const U &y)
	{
		return (isEqual<U>(y, 0))?0:static_cast<T>(x)/static_cast<T>(y);
	}

	template <class T>
	DEVICE_CALLABLE FORCE_INLINE
	bool Check_Bound(const T &x, const T &x_min, const T &x_max)
	{
		return (x_min <= x) && (x <= x_max);
	}

	template <class T, class U>
	DEVICE_CALLABLE FORCE_INLINE
	bool Check_Bound(const T &x, const T &x_min, const T &x_max, const U &y, const U &y_min, const U &y_max)
	{
		return (x_min <= x) && (x <= x_max) && (y_min <= y) && (y <= y_max);
	}


	/******************************Ranges************************************/
	struct Range_1d
	{
		int ix_0; 	// initial index
		int ix_e; 	// final index
		int ixy_0; 	// initial index
		int ixy_e; 	// final index
		Range_1d(): ix_0(0), ix_e(0), ixy_0(0), ixy_e(0){}

		Range_1d(int ix_0i, int ix_ei): ix_0(ix_0i), ix_e(ix_ei)
		{
			ixy_0 = 0;
			ixy_e = nx();
		}

		template <class TGrid>
		Range_1d(const TGrid &grid_1d){ set_grid(grid_1d); }

		void clear()
		{
			ix_0 = 0;
			ix_e = 0;
			ixy_0 = 0;
			ixy_e = 0;
		}

		int nx() const { return (ix_e-ix_0);}

		void clip_ix(int ix_0i, int ix_ei)
		{
			ix_0 = min(ix_ei, max(ix_0i, ix_0));
			ix_e = min(ix_ei, max(ix_0i, ix_e));
		}

		void set_ascending_index()
		{
			if(ix_0>ix_e)
			{
				std::swap(ix_0, ix_e);
			}
		}

		template <class TGrid>
		void set_grid(const TGrid &grid_1d)
		{
			ix_0 = 0;
			ix_e = grid_1d.nx;
			ixy_0 = 0;
			ixy_e = grid_1d.nx;
		}
	};

	struct Range_2d
	{
		int ix_0; 	// initial index
		int ix_e; 	// final index
		int iy_0; 	// initial index
		int iy_e; 	// final index
		int ixy_0; 	// initial index
		int ixy_e; 	// final index
		Range_2d(): ix_0(0), ix_e(0),
		iy_0(0), iy_e(0), ixy_0(0), ixy_e(0){}

		Range_2d(int ix_0i, int ix_ei, int iy_0i, int iy_ei): ix_0(ix_0i), ix_e(ix_ei),
		iy_0(iy_0i), iy_e(iy_ei)
		{
			ixy_0 = 0;
			ixy_e = nxy();
		}

		template <class TGrid>
		Range_2d(const TGrid &grid_2d){ set_grid(grid_2d); }

		void clear()
		{
			ix_0 = 0;
			ix_e = 0;
			iy_0 = 0;
			iy_e = 0;
			ixy_0 = 0;
			ixy_e = 0;
		}

		template <class TGrid>
		void set_grid(const TGrid &grid_2d)
		{
			ix_0 = 0;
			ix_e = grid_2d.nx;
			iy_0 = 0;
			iy_e = grid_2d.ny;
			ixy_0 = 0;
			ixy_e = grid_2d.nxy();
		}

		DEVICE_CALLABLE FORCE_INLINE
		int nx() const { return (ix_e-ix_0);}

		DEVICE_CALLABLE FORCE_INLINE
		int ny() const { return (iy_e-iy_0);}

		DEVICE_CALLABLE FORCE_INLINE
		int nxy() const { return (nx()*ny());}

		DEVICE_CALLABLE FORCE_INLINE
		void clip_ix(int ix_0i, int ix_ei)
		{
			ix_0 = min(ix_ei, max(ix_0i, ix_0));
			ix_e = min(ix_ei, max(ix_0i, ix_e));
		}

		DEVICE_CALLABLE FORCE_INLINE
		void clip_iy(int iy_0i, int iy_ei)
		{
			iy_0 = min(iy_ei, max(iy_0i, iy_0));
			iy_e = min(iy_ei, max(iy_0i, iy_e));
		}

		DEVICE_CALLABLE FORCE_INLINE
		void set_ascending_index()
		{
			if(ix_0>ix_e)
			{
				std::swap(ix_0, ix_e);
			}

			if(iy_0>iy_e)
			{
				std::swap(iy_0, iy_e);
			}
		}

		DEVICE_CALLABLE FORCE_INLINE
 		bool chk_ix_bound(const int &ix) const
		{
			return (ix_0<=ix) && (ix<ix_e);
		}

		DEVICE_CALLABLE FORCE_INLINE
 		bool chk_iy_bound(const int &iy) const
		{
			return (iy_0<=iy) && (iy<iy_e);
		}

		DEVICE_CALLABLE FORCE_INLINE
		bool chk_bound(const int &ix, const int &iy) const
		{
			return chk_ix_bound(ix)&&chk_iy_bound(iy);
		}

 		DEVICE_CALLABLE FORCE_INLINE
		int ind_col_o(const int &ix, const int &iy) const
		{
			return ((ix-ix_0)*ny()+(iy-iy_0));
		}
	};

	/**************************borders**************************/
	template <class T>
	struct Border_1d
	{
		T lx;

		T xb_0;
		T xb_e;

		Border_1d(T lxi=T(), T xb_0i=T(), T xb_ei=T()): lx(lxi), xb_0(xb_0i), xb_e(xb_ei){}

		template <class X>
		Border_1d(T lxi, int nptr, X *ptr): lx(lxi)
		{
			if(nptr>0)
			{
				xb_0 = ptr[0];
				xb_e = (nptr>1)?ptr[1]:0;
			}
			else
			{
				xb_0 = 0;
				xb_e = 0;
			}
		}

		Border_1d(T lxi, T dx): lx(lxi)
		{
			set_bd(dx);
		}

		template <class X>
		Border_1d<T>& operator=(Border_1d<X> &border_1d)
		{
			lx = border_1d.lx;
			xb_0 = border_1d.xb_0;
			xb_e = border_1d.xb_e;
			return *this;
		}

		inline
		void set_bd(const T &dx)
		{
			xb_0 = max(xb_0, max(dx, T(0)));
			xb_e = max(xb_e, max(-dx, T(0)));
		}

		T lx_wb() const
		{
			return ::fmax(x_e()-x_0(), T(0));
		}

		T x_c() const
		{
			return 0.5*(x_e()+x_0());
		}

		T x_0() const
		{
			return xb_0;
		}

		T x_e() const
		{
			return lx-xb_e;
		}

		T xb_max() const
		{
			return max(xb_0, xb_e);
		}

		T xb_min() const
		{
			return min(xb_0, xb_e);
		}

		// percentaje of the total length
		T radius_ptl(T p) const
		{
			return (1-0.5*p)*lx_wb()/2;
		}

		bool chk_bound(const T &x) const
		{
			return (x_0()<=x) && (x<=x_e());
		}

		void shift(T dx)
		{
			xb_0 = ::fmax(xb_0+dx, T(0));
			xb_e = ::fmax(xb_e-dx, T(0));
		}

		void clear()
		{
			xb_0 = xb_e = 0;
		}
	};

	template <class T>
	struct Border_2d
	{
		T lx;
		T ly;

		T xb_0;
		T xb_e;
		T yb_0;
		T yb_e;

		Border_2d(T lxi=T(), T lyi=T(), T xb_0i=T(), T xb_ei=T(),
		T yb_0i=T(), T yb_ei=T()): lx(lxi), ly(lyi),
		xb_0(xb_0i), xb_e(xb_ei), yb_0(yb_0i), yb_e(yb_ei){}

		template <class X>
		Border_2d(T lxi, T lyi, int nptr, X *ptr, bool b_swap = false): lx(lxi), ly(lyi)
		{
			if(nptr>0)
			{
				xb_0 = ptr[0];
				xb_e = (nptr>1)?ptr[1]:0;
				yb_0 = (nptr>2)?ptr[2]:0;
				yb_e = (nptr>3)?ptr[3]:0;
			}
			else
			{
				xb_0 = 0;
				xb_e = 0;
				yb_0 = 0;
				yb_e = 0;
			}

			if(b_swap)
			{
				thrust::swap(xb_0, yb_0);
				thrust::swap(xb_e, yb_e);
			}
		}

		Border_2d(T lxi, T lyi, r2d<T> rd): lx(lxi), ly(lyi)
		{
			set_bd(rd);
		}

		template <class X>
		Border_2d<T>& operator=(Border_2d<X> &border_2d)
		{
			lx = border_2d.lx;
			ly = border_2d.ly;
			xb_0 = border_2d.xb_0;
			xb_e = border_2d.xb_e;
			yb_0 = border_2d.yb_0;
			yb_e = border_2d.yb_e;
			return *this;
		}

		inline
		void set_bd(const r2d<T> &rd)
		{
			xb_0 = max(xb_0, max(rd.x, T(0)));
			xb_e = max(xb_e, max(-rd.x, T(0)));

			yb_0 = max(yb_0, max(rd.y, T(0)));
			yb_e = max(yb_e, max(-rd.y, T(0)));
		}

		T lx_wb() const
		{
			return ::fmax(x_e()-x_0(), T(0));
		}

		T ly_wb() const
		{
			return ::fmax(y_e()-y_0(), T(0));
		}

		T x_c() const
		{
			return 0.5*(x_e()+x_0());
		}

		T y_c() const
		{
			return 0.5*(y_e()+y_0());
		}

		T x_0() const
		{
			return xb_0;
		}

		T y_0() const
		{
			return yb_0;
		}

		T x_e() const
		{
			return lx-xb_e;
		}

		T y_e() const
		{
			return ly-yb_e;
		}

		T xb_max() const
		{
			return max(xb_0, xb_e);
		}

		T xb_min() const
		{
			return min(xb_0, xb_e);
		}

		T yb_max() const
		{
			return max(yb_0, yb_e);
		}

		T yb_min() const
		{
			return min(yb_0, yb_e);
		}

		// percentaje of the total length
		T radius_ptl(T p) const
		{
			return (1-0.5*p)*::fmin(lx_wb(), ly_wb())/2;
		}

		bool chk_x_bound(const T &x) const
		{
			return (x_0()<=x) && (x<=x_e());
		}

		bool chk_y_bound(const T &y) const
		{
			return (y_0()<=y) && (y<=y_e());
		}

		bool chk_bound(const T &x, const T &y) const
		{
			return chk_x_bound(x)&&chk_y_bound(y);
		}

		void shift(r2d<T> dr)
		{
			xb_0 = ::fmax(xb_0+dr.x, T(0));
			xb_e = ::fmax(xb_e-dr.x, T(0));

			yb_0 = ::fmax(yb_0+dr.y, T(0));
			yb_e = ::fmax(yb_e-dr.y, T(0));
		}

		void clear()
		{
			xb_0 = xb_e = yb_0 = yb_e = 0;
		}
	};

	/************************Amorphous layer information*************************/
	template <class T>
	struct Amorp_Lay_Info
	{
		Amorp_Lay_Info(): z_0(0), z_e(0), dz(2), region(0), type(eALT_none) {};

		T z_0; 							// Initial z-position
		T z_e; 							// Final z-position
		T dz;							// Slice thickness
		int region; 					// Region
		eAmorp_Lay_Type type; 			// amorphous layer type

		template <class TAmorp_Lay_Info>
		Amorp_Lay_Info<T>& operator=(TAmorp_Lay_Info &amorp_lay_info)
		{
			z_0 = amorp_lay_info.z_0;
			z_e = amorp_lay_info.z_e;
			dz = amorp_lay_info.dz;
			region = amorp_lay_info.region;
			type = amorp_lay_info.type;
			return *this;
		}

		T lz() const { return fabs(z_e-z_0); }
		bool is_at_top() const { return type==eALT_Top; };
		bool is_at_bottom() const { return type==eALT_Bottom; };
		bool is_at_middle() const { return type==eALT_Middle; };

		void set_region(Atom_Data<T> &atoms)
		{
			if(lz()<1e-4)
			{
				region = 0;
				return;
			}

			int f_region = 0;
			int c_region = 0;
			for(auto iatoms = 0; iatoms < atoms.size(); iatoms++)
			{
				auto z = atoms.z[iatoms];
				if((z_0<z) && (z<z_e))
				{
					f_region += atoms.region[iatoms];
					c_region++;
				}
			}
			c_region = max(1, c_region);
			region = static_cast<int>(round(T(f_region)/T(c_region)));
		}
	};

	/*************************slice thickness*****************/
	template <class T>
	struct Slice
	{
		Slice(): z_0(0), z_e(0), z_int_0(0),
		z_int_e(0), iatom_0(1), iatom_e(0), ithk(-1){}

		T z_0; 			// Initial z-position
		T z_e; 			// Final z-position
		T z_int_0; 		// Initial z-position
		T z_int_e; 		// Final z-position
		int iatom_0; 	// Index to initial z-position
		int iatom_e; 	// Index to final z-position
		int ithk;		// thick index

		T dz() const { return fabs(z_e-z_0); }

		T z_m() const { return 0.5*(z_e+z_0); }
	};

	/************************Thickness*************************/
	template <class T>
	struct Thick
	{
		Thick(): z(0), z_zero_def_plane(0), z_back_prop(0),
		islice(0), iatom_e(0) {}

		T z; 					// z
		T z_zero_def_plane; 	// z: Zero defocus
		T z_back_prop; 			// z: Back propagation

		int islice; 			// slice position
		int iatom_e; 			// Last atom index
	};

	/**********************Identify planes*********************/
	template <class T>
	struct Identify_Planes
	{
		using TVector = Vector<T, e_host>;
		using TVector_I = Vector<int, e_host>;

		public:
			Identify_Planes(): dv(0.1){}

			// Identify planes: Require v to be sorted
			TVector operator()(TVector &v)
			{
				TVector v_plane;

				if(v.size()==0)
				{
					return v_plane;
				}

				// min and max element
				T v_min = v.front();
				T v_max = v.back();

				// calculate hist and correct it
				auto v_hist = hist(v, dv, v_min, v_max);

				if(v_hist.size()==1)
				{
					v_plane.push_back(thrust::reduce(v.begin(), v.end())/T(v.size()));
					return v_plane;
				}

				// calculate layer limits
				TVector v_lim;
				v_lim.reserve(v_hist.size());

				for(auto iz = 0; iz < v_hist.size()-1; iz++)
				{
					if((v_hist[iz]>0) && (v_hist[iz+1]==0))
					{
						v_lim.push_back(v.front()+(iz+1)*dv);
					}
				}
				v_lim.push_back(v.back()+dv);

				// calculate planes
				v_plane.reserve(v_lim.size());

				T v_m = v.front();
				T v_m_ee = 0;
				int v_c = 1;
				int izl = 0;
				for(auto iz = 0; iz < v.size(); iz++)
				{
					auto v_v = v[iz];

					if(v_v<v_lim[izl])
					{
						host_device_detail::kh_sum(v_m, v_v, v_m_ee);
						v_c++;
					}
					else
					{
						v_plane.push_back(v_m/v_c);
						v_m = v_v;
						v_m_ee = 0;
						v_c = 1;
						izl++;
					}
				}

				v_plane.push_back(v_m/v_c);
				v_plane.shrink_to_fit();

				return v_plane;
			}

			// calculate planes
			TVector operator()(T v_min, T v_max, T dv, eMatch_Border mb=eMB_MinMax)
			{
				const T v_eps = 1e-4;

				TVector v_plane;

				if(fabs(v_max-v_min)<v_eps)
				{
					v_plane.resize(1);
					v_plane[0] = 0.5*(v_min+v_max);
					return v_plane;
				}

				if(v_max<v_min)
				{
					return v_plane;
				}

				auto quot = [v_eps](const T &a, const T &b)->int
				{
					return static_cast<int>(floor(a/b+v_eps));
				};

				T s_v = v_max-v_min;
				const int nv = max(1, quot(s_v, dv))+1;

				switch (mb)
				{
					case eMB_Min:
					{
						v_plane.resize(nv);
						for(auto iv=0; iv<nv; iv++)
						{
							v_plane[iv] = v_min + iv*dv;
						}
					}
					break;
					case eMB_Max:
					{
						v_plane.resize(nv);
						for(auto iv=0; iv<nv; iv++)
						{
							v_plane[nv-1-iv] = v_max - iv*dv;
						}
					}
					break;
					case eMB_MinMax:
					{
						const auto dv_b = dv + 0.5*(s_v-(nv-1)*dv);

						v_plane.resize(nv);
						v_plane[0] = v_min;
						for(auto iv=1; iv<nv; iv++)
						{
							const auto dv_t = ((iv==1)||(iv==nv-1))?dv_b:dv;
							v_plane[iv] = v_plane[iv-1] + dv_t;
						}
					}
					break;
				}


				return v_plane;
			}

		private:
			// calculate corrected histogram
			TVector_I hist(TVector &v, double dv, double v_min, double v_max)
			{
				const auto v_l = ::fmax(v_max-v_min, dv);
				const int nbins = static_cast<int>(ceil(v_l/dv));

				TVector_I v_hist(nbins, 0);
				for(auto iv = 0; iv< v.size(); iv++)
				{
	 auto v_id = double(v[iv]);
					auto ih = static_cast<int>(floor((v_id-v_min)/dv));
	 auto v_imin = v_min + (ih-1)*dv;
	 auto v_imax = v_imin + dv;

	 if(v_id<v_imin)
	 {
	 for (auto ik = ih; ik >= 0; ik--)
	 {
	  v_imin = v_min + (ik-1)*dv;
	  v_imax = v_imin + dv;
	  if((v_imin<=v_id) && (v_id<v_imax))
	  {
		ih = ik-1;
		break;
	  }
	 }
	 }
	 else if(v_id>v_imax)
	 {
	 for (auto ik = ih; ik < nbins; ik++)
	 {
	  v_imin = v_min + ik*dv;
	  v_imax = v_imin + dv;
	  if((v_imin<=v_id) && (v_id<v_imax))
	  {
		ih = ik;
		break;
	  }
	 }
	 }
	 ih = max(0, min(nbins-1, ih));
					v_hist[ih]++;
				}

				while(v_hist.back()==0)
				{
					v_hist.pop_back();
				}

				for(auto ih = 1; ih < v_hist.size()-1; ih++)
				{
					bool bn = (ih<v_hist.size()-2)?(0<v_hist[ih+2]):false;
					bn = (0<v_hist[ih-1]) && ((0<v_hist[ih+1])||bn);
					if((v_hist[ih]==0) && bn)
					{
						v_hist[ih] = 1;
					}
				}

				return v_hist;
			}

			double dv;
	};

	/************************quadrature x**********************/
	template <class T, eDevice dev>
	struct Q1
	{
		using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = dev;

		size_type size() const
		{
			return x.size();
		}

		void clear()
		{
			x.clear();
			w.clear();
		}

		void reserve(const size_type &new_size)
		{
			x.reserve(new_size);
			w.reserve(new_size);
		}

		void resize(const size_type &new_size, const value_type &value = value_type())
		{
			x.resize(new_size, value);
			w.resize(new_size, value);
		}

		void shrink_to_fit()
		{
			x.shrink_to_fit();
			w.shrink_to_fit();
		}

		template <class TQ1>
		void assign(TQ1 &q1)
		{
			x.assign(q1.x.begin(), q1.x.end());
			w.assign(q1.w.begin(), q1.w.end());
		}

		Vector<T, dev> x;
		Vector<T, dev> w;
	};

	template <class T>
	struct rQ1
	{
		using value_type = T;

		rQ1(): m_size(0), x(nullptr), w(nullptr){}

		template <class TQ1>
		rQ1<T>& operator = (TQ1 &q1)
		{
			m_size = q1.size();
			x = raw_pointer_cast(q1.x.data());
			w = raw_pointer_cast(q1.w.data());
			return *this;
		}

		template <class TQ1>
		rQ1(TQ1 &q1)
		{
			*this = q1;
		}

		int m_size;
		T *x;
		T *w;
	};

	/***********************quadrature xy**********************/
	template <class T, eDevice dev>
	struct Q2
	{
		using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = dev;

		size_type size() const
		{
			return x.size();
		}

		void clear()
		{
			x.clear();
			y.clear();
			w.clear();
		}

		void reserve(const size_type &new_size)
		{
			x.reserve(new_size);
			y.reserve(new_size);
			w.reserve(new_size);
		}

		void resize(const size_type &new_size, const value_type &value = value_type())
		{
			x.resize(new_size, value);
			y.resize(new_size, value);
			w.resize(new_size, value);
		}

		void shrink_to_fit()
		{
			x.shrink_to_fit();
			y.shrink_to_fit();
			w.shrink_to_fit();
		}

		template <class TQ2>
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

	template <class T>
	struct rQ2
	{
		using value_type = T;

		rQ2(): m_size(0), x(nullptr), y(nullptr), w(nullptr){}

		template <class TQ2>
		rQ2<T>& operator = (TQ2 &q2)
		{
			m_size = q2.size();
			x = raw_pointer_cast(q2.x.data());
			y = raw_pointer_cast(q2.y.data());
			w = raw_pointer_cast(q2.w.data());
			return *this;
		}

		template <class TQ2>
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
	template <class T, eDevice dev>
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

		template <class TPP_Coef>
		void assign(TPP_Coef &pp_coef)
		{
			cl.assign(pp_coef.cl.begin(), pp_coef.cl.end());
			cnl.assign(pp_coef.cnl.begin(), pp_coef.cnl.end());
		}

		Vector<T, dev> cl; 	// Lineal coefficients fep
		Vector<T, dev> cnl; // Non-Lineal coefficients fep

	};

	template <class T>
	struct rPP_Coef
	{
		using value_type = T;

		rPP_Coef(): m_size(0), cl(nullptr), cnl(nullptr){}

		template <class TPP_Coef>
		rPP_Coef<T>& operator = (TPP_Coef &rhs)
		{
			m_size = rhs.size();
			cl = raw_pointer_cast(rhs.cl.data());
			cnl = raw_pointer_cast(rhs.cnl.data());
			return *this;
		}

		template <class TPP_Coef>
		rPP_Coef(TPP_Coef &pp_coef)
		{
			*this = pp_coef;
		}

		int m_size;
		T *cl;
		T *cnl;
	};

	/************Cubic interpolation coefficients*************/
	template <class T, eDevice dev>
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

		template <class TCI_Coef>
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

	template <class T>
	struct rCI_Coef
	{
		using value_type = T;

		rCI_Coef(): m_size(0), c0(nullptr), c1(nullptr), c2(nullptr), c3(nullptr){}

		template <class TCI_Coef>
		rCI_Coef<T>& operator = (TCI_Coef &ci_coef)
		{
			m_size = ci_coef.size();
			c0 = raw_pointer_cast(ci_coef.c0.data());
			c1 = raw_pointer_cast(ci_coef.c1.data());
			c2 = raw_pointer_cast(ci_coef.c2.data());
			c3 = raw_pointer_cast(ci_coef.c3.data());
			return *this;
		}

		template <class TCI_Coef>
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

	/******************************pair<int, val>************************************/
	template <class T>
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
				load_number();
			}

			int operator()(int64_t n, eDat_Sel_Type dst = eDST_Closest)
			{
				auto p_idx = std::min_element(number.begin(), number.end(), [&n](int64_t p0, int64_t pe){ return std::abs(n-p0)<std::abs(n-pe);});

				auto pn = static_cast<int>(*p_idx);

				switch(dst)
				{
					case eDST_Less_Than:
					{
						if(pn>=n)
						{
							pn = static_cast<int>(*(p_idx-1));
						}
					}
					break;
					case eDST_Greater_Than:
					{
						if(pn<=n)
						{
							pn = static_cast<int>(*(p_idx+1));
						}
					}
					break;
				}
				return pn;
			}

		private:
			std::vector<int64_t> number;

			void load_number()
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

							for(auto ip2=1; ip2<np2; ip2++)
							{
								auto p2 = static_cast<int64_t>(std::pow(b_2, ip2));

								auto p = p7*p5*p3*p2;

								if((prime_0<p) && (p<=prime_e))
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
	};

	/***************************grids****************************/

	/***************************grid_1d****************************/
	template <class T>
	struct Grid_1d
	{
		using value_type = T;

		int nx; 				// Number of pixels in x direction

		int nxh; 				// Half number of pixels in x direction

		T lx; 					// Box m_size in x direction(Angstroms)
		T dz; 					// slice thicknes

		bool bwl; 				// Band-width limit
		bool pbc_x; 			// Peridic boundary contions

		T Rx_0; 				// starting Rx

		T dRx; 					// x-sampling in real space

		T dgx; 					// x-sampling in reciprocal space

		T gl2_max; 				// Squared of the maximun limited frequency
		T alpha;				// 1/(1+exp(alpha*(x^2-x_c^2)))

		inline
		Grid_1d(): nx(0), nxh(0), lx(0), dz(0),
		pbc_x(true), bwl(true), Rx_0(0), dRx(0), dgx(0), gl2_max(0){}

		Grid_1d(int nx_i)
		{
			lx = nx_i;
			set_input_data(nx_i, lx);
		}

		Grid_1d(int nx_i, T lx_i, T dz_i = 0.5, bool bwl_i = false, bool pbc_x_i = false, T Rx_0_i=0)
		{
			set_input_data(nx_i, lx_i, dz_i, bwl_i, pbc_x_i, Rx_0_i);
		}

		inline
		void set_input_data(int nx_i, T lx_i, T dz_i = 0.5, bool bwl_i = false, bool pbc_x_i = false, T Rx_0_i=0)
		{
			nx = nx_i;
			nxh = nx/2;
			lx = lx_i;
			dz = dz_i;
			bwl = bwl_i;
			pbc_x = pbc_x_i;
			Rx_0 = Rx_0_i;
			dRx = mt::Div(lx, nx);
			dgx = mt::Div(1.0, lx);
			gl2_max = ::square(gl_max());

			// y = 1/(1+exp(alpha*(x^2-x_c^2)))
			T dg0 = 0.25, fg0 = 1e-02;
			alpha = log(1.0/fg0-1.0)/(::square(gl_max()+dg0)-gl2_max);
		}

		inline
		void set_R_0(T Rx_0_i)
		{
			Rx_0 = Rx_0_i;
		}

		template <class TGrid>
		void assign(TGrid &grid_1d)
		{
			set_input_data(grid_1d.nx, grid_1d.lx, grid_1d.dz, grid_1d.bwl, grid_1d.pbc_x, grid_1d.Rx_0);
		}

		template <class TGrid>
		Grid_1d<T>& operator=(TGrid &grid_1d)
		{
			assign(grid_1d);
			return *this;
		}

		// Maximun limited frequency
		DEVICE_CALLABLE FORCE_INLINE
		T gl_max() const
		{
			return 2.0*g_max()/3.0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T nx_r() const { return T(nx); }

		DEVICE_CALLABLE FORCE_INLINE
		T lxh() const { return 0.5*lx; }

		DEVICE_CALLABLE FORCE_INLINE
		int Rx_c() const { return (Rx_0 + lxh()); }

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		int floor_dRx(const T &x) const { return static_cast<int>(floor((x-Rx_0)/dRx)); }

		DEVICE_CALLABLE FORCE_INLINE
		int ceil_dRx(const T &x) const { return static_cast<int>(ceil((x-Rx_0)/dRx)); }

		// lower bound
		DEVICE_CALLABLE FORCE_INLINE
		int lb_index_x(const T &x) const
		{
			return min(nx, max(0, floor_dRx(x)));
		}

		// upper bound
		DEVICE_CALLABLE FORCE_INLINE
		int ub_index_x(const T &x) const
		{
			return min(nx, max(0, ceil_dRx(x)));
		}

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		bool ckb_bound(const T &x) const
		{
			const T eps = epsilon_abs<T>();
			const T x_0 = Rx_first()-eps;
			const T x_e = Rx_last()+eps;
			return ((x<x_0)||(x>x_e))?false:true;
		}
		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		void set_x_bound(T &x) const
		{
			x = ::fmin(::fmax(x, Rx_first()), Rx_last());
		}
		/*********************************************************/
		// index
		DEVICE_CALLABLE FORCE_INLINE
		int ix(const T &x) const
		{
			return lb_index_x(x);
		}

		// range
		DEVICE_CALLABLE FORCE_INLINE
		Range_1d index_range(T x, const T &radius) const
		{
			Range_1d r;

			r.ix_0 = lb_index_x(x - radius);
			r.ix_e = ub_index_x(x + radius);

			r.ixy_0 = 0;
			r.ixy_e = r.ix_e-r.ix_0;

			return r;
		}

		/*********************************************************/
		// Maximun frequency
		DEVICE_CALLABLE FORCE_INLINE
		T g_max() const { return static_cast<T>(nxh)*dgx; }

		// Squared of the maximun frequency
		DEVICE_CALLABLE FORCE_INLINE
		T g2_max() const { return pow(g_max(), 2); }

		DEVICE_CALLABLE FORCE_INLINE
		int nx_dRx(const T &lx) const { return ceil_dRx(lx); }

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		int igx(const int &ix) const { return ix-nxh; }

		DEVICE_CALLABLE FORCE_INLINE
		T gx(const int &ix, T gx_0=T()) const { return (igx(ix)*dgx-gx_0); }

		DEVICE_CALLABLE FORCE_INLINE
		T g2(const int &ix, T gx_0=T()) const { return ::square(gx(ix, gx_0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T g(const int &ix, T gx_0=T()) const { return fabs(gx(ix, gx_0)); }


		DEVICE_CALLABLE FORCE_INLINE
		T gx_first() const { return gx(0); }

		DEVICE_CALLABLE FORCE_INLINE
		T gx_last() const { return gx(nx-1); }

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		int iRx(const int &ix) const { return ix; }

		DEVICE_CALLABLE FORCE_INLINE
		T Rx(const int &ix, T x0=T()) const { return (iRx(ix)*dRx+Rx_0-x0); }

		DEVICE_CALLABLE FORCE_INLINE
		T R2(const int &ix, T x0=T()) const { return ::square(Rx(ix, x0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T R(const int &ix, T x0=T()) const { return fabs(Rx(ix, x0)); }


		DEVICE_CALLABLE FORCE_INLINE
		T Rx_first() const { return Rx(0); }

		DEVICE_CALLABLE FORCE_INLINE
		T Rx_last() const { return Rx(nx-1); }

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		int igx_shift(const int &ix) const { return (ix<nxh)?ix:ix-nx; }

		DEVICE_CALLABLE FORCE_INLINE
		T gx_shift(const int &ix, T gx_0=T()) const { return (igx_shift(ix)*dgx-gx_0); }

		DEVICE_CALLABLE FORCE_INLINE
		T g2_shift(const int &ix, T gx_0=T()) const { return ::square(gx_shift(ix, gx_0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T g_shift(const int &ix, T gx_0=T()) const { return fabs(gx_shift(ix, gx_0)); }

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		int iRx_shift(const int &ix) const { return (ix<nxh)?ix+nxh:ix-nxh; }

		DEVICE_CALLABLE FORCE_INLINE
		T Rx_shift(const int &ix, T x0=T()) const { return (iRx_shift(ix)*dRx+Rx_0-x0); }

		DEVICE_CALLABLE FORCE_INLINE
		T R2_shift(const int &ix, T x0=T()) const { return ::square(Rx_shift(ix, x0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T R_shift(const int &ix, T x0=T()) const { return fabs(Rx_shift(ix, x0)); }

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		T bwl_factor(const int &ix) const
		{
			return (bwl)?1.0/(1.0+exp(alpha*(g2(ix)-gl2_max))):1;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T bwl_factor_shift(const int &ix) const
		{
			return (bwl)?1.0/(1.0+exp(alpha*(g2_shift(ix)-gl2_max))):1;
		}

		#ifdef __CUDACC__
		Grid_BT cuda_grid(const dim3 Blk_max = dim3(0, 0, 0))
		{
			Grid_BT grid_bt;
			grid_bt.Blk = dim3((nx+c_thrnxy-1)/c_thrnxy);
			grid_bt.Thr = dim3(c_thrnxy);

			if(Blk_max.x != 0)
			{
				grid_bt.Blk.x = min(Blk_max.x, grid_bt.Blk.x);
			}
			if(Blk_max.y != 0)
			{
				grid_bt.Blk.y = min(Blk_max.y, grid_bt.Blk.y);
			}
			if(Blk_max.z != 0)
			{
				grid_bt.Blk.z = min(Blk_max.z, grid_bt.Blk.z);
			}

			return grid_bt;
		}

		Grid_BT cuda_grid_h(const dim3 Blk_max = dim3(0, 0, 0))
		{
			return cuda_grid(dim3((nxh+c_thrnxy-1)/c_thrnxy));
		}
		#endif
	};

	template <class T>
	struct Grid_2d
	{
		using value_type = T;

		int nx; 				// Number of pixels in x direction
		int ny; 				// Number of pixels in y direction

		int nxh; 				// Half number of pixels in x direction
		int nyh; 				// Half number of pixels in y direction

		T lx; 					// Box m_size in x direction(Angstroms)
		T ly; 					// Box m_size in y direction(Angstroms)
		T dz; 					// slice thicknes

		bool bwl; 				// Band-width limit
		bool pbc_xy; 			// Peridic boundary contions

		T Rx_0; 				// starting Rx
		T Ry_0; 				// starting Ry

		T dRx; 					// x-sampling in real space
		T dRy; 					// y-sampling in real space

		T dgx; 					// x-sampling in reciprocal space
		T dgy; 					// y-sampling in reciprocal space

		T gl2_max; 				// Squared of the maximun limited frequency
		T alpha;				// 1/(1+exp(alpha*(x^2-x_c^2)))

		inline
		Grid_2d(): nx(0), ny(0), nxh(0), nyh(0),
			lx(0), ly(0), dz(0), pbc_xy(true), bwl(false),
			Rx_0(0), Ry_0(0), dRx(0), dRy(0), dgx(0), dgy(0), gl2_max(0){}

		Grid_2d(int nx_i, int ny_i)
		{
			lx = nx_i;
			ly = ny_i;
			set_input_data(nx_i, ny_i, lx, ly);
		}

		Grid_2d(int nx_i, int ny_i, T lx_i, T ly_i, T dz_i = 0.5, bool bwl_i = false,
		bool pbc_xy_i = false, T Rx_0_i=0, T Ry_0_i=0)
		{
			set_input_data(nx_i, ny_i, lx_i, ly_i, dz_i, bwl_i, pbc_xy_i, Rx_0_i, Ry_0_i);
		}

		inline
		void set_input_data(int nx_i, int ny_i, T lx_i, T ly_i, T dz_i = 0.5,
		bool bwl_i = false, bool pbc_xy_i = false, T Rx_0_i=0, T Ry_0_i=0)
		{
			nx = nx_i;
			ny = ny_i;
			nxh = nx/2;
			nyh = ny/2;
			lx = lx_i;
			ly = ly_i;
			dz = dz_i;
			bwl = bwl_i;
			pbc_xy = pbc_xy_i;
			Rx_0 = Rx_0_i;
			Ry_0 = Ry_0_i;
			dRx = mt::Div(lx, nx);
			dRy = mt::Div(ly, ny);
			dgx = mt::Div(T(1.0), lx);
			dgy = mt::Div(T(1.0), ly);
			gl2_max = ::square(gl_max());

			// y = 1/(1+exp(alpha*(x^2-x_c^2)))
			T dg0 = 0.25, fg0 = 1e-02;
			alpha = log(1.0/fg0-1.0)/(::square(gl_max()+dg0)-gl2_max);
		}

		inline
		void set_R_0(T Rx_0_i, T Ry_0_i)
		{
			Rx_0 = Rx_0_i;
			Ry_0 = Ry_0_i;
		}

		template <class TGrid>
		void assign(TGrid &grid_2d)
		{
			set_input_data(grid_2d.nx, grid_2d.ny, grid_2d.lx, grid_2d.ly, grid_2d.dz, grid_2d.bwl, grid_2d.pbc_xy, grid_2d.Rx_0, grid_2d.Ry_0);
		}

		template <class TGrid>
		Grid_2d<T>& operator=(TGrid &grid_2d)
		{
			assign(grid_2d);
			return *this;
		}

		/************index periodic boundary conditions************/
		// https:// en.wikipedia.org/wiki/Periodic_boundary_conditions

		DEVICE_CALLABLE FORCE_INLINE
		int iRx_pbc(const int &ix) const
		{
			return (ix - static_cast<int>(floor(Rx(ix)/lx))*nx);
		}

		DEVICE_CALLABLE FORCE_INLINE
		int iRy_pbc(const int &iy) const
		{
			return (iy - static_cast<int>(floor(Ry(iy)/ly))*ny);
		}

		DEVICE_CALLABLE FORCE_INLINE
		void iRx_iRy_pbc(int &ix, int &iy)
		{
			ix = iRx_pbc(ix);
			iy = iRy_pbc(iy);
		}

		DEVICE_CALLABLE FORCE_INLINE
		int ind_row_pbc(const int &ix, const int &iy) const { return iRy_pbc(iy)*nx+iRx_pbc(ix); }

		DEVICE_CALLABLE FORCE_INLINE
		int ind_col_pbc(const int &ix, const int &iy) const { return iRx_pbc(ix)*ny+iRy_pbc(iy); }

		DEVICE_CALLABLE FORCE_INLINE
		int ind_col_pbc_shift(int ix, int iy)
		{
			iRx_iRy_pbc(ix, iy);
			iRx_iRy_shift(ix, iy);
			return ix*ny+iy;
		}

		/*********************************************************/
		// Maximun limited frequency
		DEVICE_CALLABLE FORCE_INLINE
		T gl_max() const
		{
			return 2.0*g_max()/3.0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		int nxy() const { return nx*ny; }

		T inxy() const { return 1.0/nxy_r(); }

		DEVICE_CALLABLE FORCE_INLINE
		T nx_r() const { return T(nx); }

		DEVICE_CALLABLE FORCE_INLINE
		T ny_r() const { return T(ny); }

		DEVICE_CALLABLE FORCE_INLINE
		T nxy_r() const { return T(nx*ny); }

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

		DEVICE_CALLABLE FORCE_INLINE
		int Rx_c() const { return (Rx_0 + lxh()); }

		DEVICE_CALLABLE FORCE_INLINE
		int Ry_c() const { return (Ry_0 + lyh()); }

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		int floor_dR_min(const T &x) const { return static_cast<int>(floor((x-R_0_min())/dR_min())); }

		DEVICE_CALLABLE FORCE_INLINE
		int floor_dRx(const T &x) const { return static_cast<int>(floor((x-Rx_0)/dRx)); }

		DEVICE_CALLABLE FORCE_INLINE
		int floor_dRy(const T &y) const { return static_cast<int>(floor((y-Ry_0)/dRy)); }

		DEVICE_CALLABLE FORCE_INLINE
		int ceil_dR_min(const T &x) const { return static_cast<int>(ceil((x-R_0_min())/dR_min())); }

		DEVICE_CALLABLE FORCE_INLINE
		int ceil_dRx(const T &x) const { return static_cast<int>(ceil((x-Rx_0)/dRx)); }

		DEVICE_CALLABLE FORCE_INLINE
		int ceil_dRy(const T &y) const { return static_cast<int>(ceil((y-Ry_0)/dRy)); }

		// lower bound
		DEVICE_CALLABLE FORCE_INLINE
		int lb_index_x(const T &x) const
		{
			return min(nx, max(0, floor_dRx(x)));
		}

		// upper bound
		DEVICE_CALLABLE FORCE_INLINE
		int ub_index_x(const T &x) const
		{
			return min(nx, max(0, ceil_dRx(x)));
		}

		// lower bound
		DEVICE_CALLABLE FORCE_INLINE
		int lb_index_y(const T &y) const
		{
			return min(ny, max(0, floor_dRy(y)));
		}

		// upper bound
		DEVICE_CALLABLE FORCE_INLINE
		int ub_index_y(const T &y) const
		{
			return min(ny, max(0, floor_dRy(y)));
		}
		/*********************************************************/

		DEVICE_CALLABLE FORCE_INLINE
		bool ckb_x_bound(const T &x) const
		{
			const T eps = epsilon_abs<T>();
			const T x_0 = Rx_first()-eps;
			const T x_e = Rx_last()+eps;
			return ((x<x_0)||(x>x_e))?false:true;
		}

		DEVICE_CALLABLE FORCE_INLINE
		bool ckb_y_bound(const T &y) const
		{
			const T eps = epsilon_abs<T>();
			const T y_0 = Ry_first()-eps;
			const T y_e = Ry_last()+eps;
			return ((y<y_0)||(y>y_e))?false:true;
		}

		DEVICE_CALLABLE FORCE_INLINE
		bool ckb_bound(const r2d<T> &p) const
		{
			return (ckb_x_bound(p.x) && ckb_y_bound(p.y));
		}

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		void set_x_bound(T &x) const
		{
			x = ::fmin(::fmax(x, Rx_first()), Rx_last());
		}

		DEVICE_CALLABLE FORCE_INLINE
		void set_y_bound(T &y) const
		{
			y = ::fmin(::fmax(y, Ry_first()), Ry_last());
		}

		DEVICE_CALLABLE FORCE_INLINE
		void set_bound(r2d<T> &p) const
		{
			set_x_bound(p.x);
			set_y_bound(p.y);
		}
		/*********************************************************/
		// index
		DEVICE_CALLABLE FORCE_INLINE
		int ixy(const T &x, const T &y) const
		{
			auto ix = lb_index_x(x);
			auto iy = lb_index_y(y);
			return ind_col(ix, iy);
		}

		// range
		DEVICE_CALLABLE FORCE_INLINE
		Range_2d index_range(r2d<T> p, T radius) const
		{
			Range_2d r;

			r.ix_0 = lb_index_x(p.x - radius);
			r.ix_e = ub_index_x(p.x + radius);

			r.iy_0 = lb_index_y(p.y - radius);
			r.iy_e = ub_index_y(p.y + radius);

			r.ixy_0 = 0;
			r.ixy_e = (r.ix_e-r.ix_0)*(r.iy_e-r.iy_0);

			return r;
		}

		Range_2d index_range(r2d<T> p, T f0, T a, T b, T c)
		{
			Range_2d r;

			T d = log(f0);
			T dd = c*c-4*a*b;

			T radius_x = sqrt(4*b*d/dd);
			T radius_y = sqrt(4*a*d/dd);

			r.ix_0 = ub_index_x(p.x - radius_x);
			r.ix_e = lb_index_x(p.x + radius_x);

			r.iy_0 = ub_index_y(p.y - radius_y);
			r.iy_e = lb_index_y(p.y + radius_y);

			r.ixy_0 = 0;
			r.ixy_e = (r.ix_e-r.ix_0)*(r.iy_e-r.iy_0);

			return r;
		};

		/*********************************************************/
		// Maximun frequency
		DEVICE_CALLABLE FORCE_INLINE
		T g_max() const { return ::fmin(static_cast<T>(nxh)*dgx, static_cast<T>(nyh)*dgy); }

		// Squared of the maximun frequency
		DEVICE_CALLABLE FORCE_INLINE
		T g2_max() const { return ::square(g_max()); }

		DEVICE_CALLABLE FORCE_INLINE
		T R_0_min() const { return ::fmin(Rx_0, Ry_0); }

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
		int igx(const int &ix) const { return (ix-nxh); }

		DEVICE_CALLABLE FORCE_INLINE
		int igy(const int &iy) const { return (iy-nyh); }

		DEVICE_CALLABLE FORCE_INLINE
		T gx(const int &ix, T gx_0=T()) const { return (igx(ix)*dgx-gx_0); }

		DEVICE_CALLABLE FORCE_INLINE
		T gy(const int &iy, T gy_0=T()) const { return (igy(iy)*dgy-gy_0); }

		DEVICE_CALLABLE FORCE_INLINE
		T gx2(const int &ix, T gx_0=T()) const { return ::square(gx(ix, gx_0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T gy2(const int &iy, T gy_0=T()) const { return ::square(gy(iy, gy_0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T g2(const int &ix, const int &iy, T gx_0=T(), T gy_0=T())const { return (gx2(ix, gx_0)+gy2(iy, gy_0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T g(const int &ix, const int &iy, T gx_0=T(), T gy_0=T()) const { return sqrt(g2(ix, iy, gx_0, gy_0)); }


		DEVICE_CALLABLE FORCE_INLINE
		T gx_first() const { return gx(0); }

		DEVICE_CALLABLE FORCE_INLINE
		T gx_last() const { return gx(nx-1); }

		DEVICE_CALLABLE FORCE_INLINE
		T gy_first() const { return gy(0); }

		DEVICE_CALLABLE FORCE_INLINE
		T gy_last() const { return gy(ny-1); }

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		int iRx(const int &ix) const { return ix; }

		DEVICE_CALLABLE FORCE_INLINE
		int iRy(const int &iy) const { return iy; }

		DEVICE_CALLABLE FORCE_INLINE
		T Rx(const int &ix, T x0=T()) const { return (iRx(ix)*dRx+Rx_0-x0); }

		DEVICE_CALLABLE FORCE_INLINE
		T Ry(const int &iy, T y0=T()) const { return (iRy(iy)*dRy+Ry_0-y0); }

		DEVICE_CALLABLE FORCE_INLINE
		T Rx2(const int &ix, T x0=T()) const { return ::square(Rx(ix, x0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T Ry2(const int &iy, T y0=T()) const { return ::square(Ry(iy, y0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T R2(const int &ix, const int &iy, T x0=T(), T y0=T()) const { return (Rx2(ix, x0)+Ry2(iy, y0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T R(const int &ix, const int &iy, T x0=T(), T y0=T()) const { return sqrt(R2(ix, iy, x0, y0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T Rx_first() const { return Rx(0); }

		DEVICE_CALLABLE FORCE_INLINE
		T Rx_last() const { return Rx(nx-1); }

		DEVICE_CALLABLE FORCE_INLINE
		T Ry_first() const { return Ry(0); }

		DEVICE_CALLABLE FORCE_INLINE
		T Ry_last() const { return Ry(ny-1); }

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		int igx_shift(const int &ix) const { return (ix<nxh)?ix:(ix-nx); }

		DEVICE_CALLABLE FORCE_INLINE
		int igy_shift(const int &iy) const { return (iy<nyh)?iy:(iy-ny); }

		DEVICE_CALLABLE FORCE_INLINE
		T gx_shift(const int &ix, T gx_0=T()) const { return (igx_shift(ix)*dgx-gx_0); }

		DEVICE_CALLABLE FORCE_INLINE
		T gy_shift(const int &iy, T gy_0=T()) const { return (igy_shift(iy)*dgy-gy_0); }

		DEVICE_CALLABLE FORCE_INLINE
		T gx2_shift(const int &ix, T gx_0=T()) const { return ::square(gx_shift(ix, gx_0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T gy2_shift(const int &iy, T gy_0=T()) const { return ::square(gy_shift(iy, gy_0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T g2_shift(const int &ix, const int &iy, T gx_0=T(), T gy_0=T()) const { return (gx2_shift(ix, gx_0)+gy2_shift(iy, gy_0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T g_shift(const int &ix, const int &iy, T gx_0=T(), T gy_0=T()) const { return sqrt(g2_shift(ix, iy, gx_0, gy_0)); }

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		int iRx_shift(const int &ix) const { return ((ix<nxh)?ix+nxh:ix-nxh); }

		DEVICE_CALLABLE FORCE_INLINE
		int iRy_shift(const int &iy) const { return ((iy<nyh)?iy+nyh:iy-nyh); }

		DEVICE_CALLABLE FORCE_INLINE
		void iRx_iRy_shift(int &ix, int &iy)
		{
			ix = iRx_shift(ix);
			iy = iRy_shift(iy);
		}

		DEVICE_CALLABLE FORCE_INLINE
		T Rx_shift(const int &ix, T x0=T()) const { return (iRx_shift(ix)*dRx+Rx_0-x0); }

		DEVICE_CALLABLE FORCE_INLINE
		T Ry_shift(const int &iy, T y0=T()) const { return (iRy_shift(iy)*dRy+Ry_0-y0); }

		DEVICE_CALLABLE FORCE_INLINE
		T Rx2_shift(const int &ix, T x0=T()) const { return ::square(Rx_shift(ix, x0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T Ry2_shift(const int &iy, T y0=T()) const { return ::square(Ry_shift(iy, y0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T R2_shift(const int &ix, const int &iy, T x0=T(), T y0=T()) const { return (Rx2_shift(ix, x0)+Ry2_shift(iy, y0)); }

		DEVICE_CALLABLE FORCE_INLINE
		T R_shift(const int &ix, const int &iy, T x0=T(), T y0=T()) const { return sqrt(R2_shift(ix, iy, x0, y0)); }

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		T bwl_factor(const int &ix, const int &iy) const
		{
			return (bwl)?1.0/(1.0+exp(alpha*(g2(ix, iy)-gl2_max))):1;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T bwl_factor_shift(const int &ix, const int &iy) const
		{
			return (bwl)?1.0/(1.0+exp(alpha*(g2_shift(ix, iy)-gl2_max))):1.0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T exp_factor_Rx(const T &x) const
		{
			const T c_2Pi = 6.283185307179586476925;
			return -c_2Pi*(x-lxh());
		}

		DEVICE_CALLABLE FORCE_INLINE
		T exp_factor_Ry(const T &y) const
		{
			const T c_2Pi = 6.283185307179586476925;
			return -c_2Pi*(y-lyh());
		}

		/*********************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		int ind_row(const int &ix, const int &iy) const { return (iy*nx+ix); }

		DEVICE_CALLABLE FORCE_INLINE
		int ind_col(const int &ix, const int &iy) const { return (ix*ny+iy); }

		DEVICE_CALLABLE FORCE_INLINE
		void col_row(const int &ixy, int &ix, int &iy) const
		{
			ix = ixy/ny;
			iy = ixy - ix*ny;
		}

		#ifdef __CUDACC__
		Grid_BT cuda_grid(const dim3 Blk_max = dim3(0, 0, 0))
		{
			Grid_BT grid_bt;
			grid_bt.Blk = dim3((ny+c_thrnxny-1)/c_thrnxny, (nx+c_thrnxny-1)/c_thrnxny);
			grid_bt.Thr = dim3(c_thrnxny, c_thrnxny);

			if(Blk_max.x != 0)
			{
				grid_bt.Blk.x = min(Blk_max.x, grid_bt.Blk.x);
			}
			if(Blk_max.y != 0)
			{
				grid_bt.Blk.y = min(Blk_max.y, grid_bt.Blk.y);
			}
			if(Blk_max.z != 0)
			{
				grid_bt.Blk.z = min(Blk_max.z, grid_bt.Blk.z);
			}

			return grid_bt;
		}

		Grid_BT cuda_grid_h(const dim3 Blk_max = dim3(0, 0, 0))
		{
			return cuda_grid(dim3((nyh+c_thrnxny-1)/c_thrnxny, (nxh+c_thrnxny-1)/c_thrnxny));
		}
		#endif
	};


	/**********************STEM Detector**********************/
	inline
	bool is_detector_circular(const eDetector_Type &det_type)
	{
		return det_type == mt::eDT_Circular;
	}

	inline
	bool is_detector_radial(const eDetector_Type &det_type)
	{
		return det_type == mt::eDT_Radial;
	}

	inline
	bool is_detector_matrix(const eDetector_Type &det_type)
	{
		return det_type == mt::eDT_Matrix;
	}

	template <class T, eDevice dev>
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

		void clear()
		{
			g_inner.clear();
			g_outer.clear();
			fx.clear();
			fR.clear();
			fn.clear();
			grid_1d.clear();
			grid_2d.clear();
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
					fn.resize(new_size);
					grid_1d.resize(new_size);
				}
					break;
				case eDT_Matrix:
				{
					fR.resize(new_size);
					fn.resize(new_size);
					grid_2d.resize(new_size);
				}
				break;
			}
		}

		template <class TDetector>
		void assign(TDetector &detector)
		{
			type = detector.type;
			g_inner.assign(detector.g_inner.begin(), detector.g_inner.end());
			g_outer.assign(detector.g_outer.begin(), detector.g_outer.end());

			fx.resize(detector.fx.size());
			for(auto i= 0; i<detector.fx.size(); i++)
			{
				fx[i].assign(detector.fx[i].begin(), detector.fx[i].end());
				//fn[i] = detector.fn[i];
				//grid_1d[i] = detector.grid_1d[i];
			}

			fR.resize(detector.fR.size());
			for(auto i= 0; i<detector.fR.size(); i++)
			{
				fR[i].assign(detector.fR[i].begin(), detector.fR[i].end());
				//fn[i] = detector.fn[i];
				//grid_2d[i] = detector.grid_2d[i];
			}
		}

		template <class TDetector>
		Detector<T, dev>& operator=(TDetector &detector)
		{
			assign(detector);
			return *this;
		}

		bool is_detector_circular() const
		{
			return mt::is_detector_circular(type);
		}

		bool is_detector_radial() const
		{
			return mt::is_detector_radial(type);
		}

		bool is_detector_matrix() const
		{
			return mt::is_detector_matrix(type);
		}

		eDetector_Type type;					// eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3
		Vector<T, e_host> g_inner;				// Inner aperture Ang^-1
		Vector<T, e_host> g_outer;				// Outer aperture Ang^-1
		Vector<Vector<T, dev>, e_host> fx;		// radial sensitivity value
		Vector<Vector<T, dev>, e_host> fR;		// 2D sensitivity value
		std::vector<Grid_2d<T>> grid_1d;		// grid_1d
		std::vector<Grid_2d<T>> grid_2d;		// grid_2d
		std::vector<std::string> fn;			// file names
	};

	/********************STEM Intensity***********************/
	template <class TVector>
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

	/****************************lens***************************/
	template <class T>
	struct Lens
	{
		using value_type = T;

		int m; 					// Momentum of the vortex

		T c_10; 				// Defocus (�)
		T c_12; 				// 2-fold astigmatism (�)
		T phi_12; 				// Azimuthal angle of 2-fold astigmatism (rad)

		T c_21; 				// Axial coma (�)
		T phi_21; 				// Azimuthal angle of axial coma (rad)
		T c_23; 				// 3-fold astigmatism (�)
		T phi_23; 				// Azimuthal angle of 3-fold astigmatism (rad)

		T c_30; 				// 3rd order spherical aberration (�)
		T c_32; 				// Axial star aberration (�)
		T phi_32; 				// Azimuthal angle of axial star aberration (rad)
		T c_34; 				// 4-fold astigmatism (�)
		T phi_34; 				// Azimuthal angle of 4-fold astigmatism (rad)

		T c_41; 				// 4th order axial coma (�)
		T phi_41; 				// Azimuthal angle of 4th order axial coma (rad)
		T c_43; 				// 3-lobe aberration (�)
		T phi_43; 				// Azimuthal angle of 3-lobe aberration (rad)
		T c_45; 				// 5-fold astigmatism (�)
		T phi_45; 				// Azimuthal angle of 5-fold astigmatism (rad)

		T c_50; 				// 5th order spherical aberration (�)
		T c_52; 				// 5th order axial star aberration (�)
		T phi_52; 				// Azimuthal angle of 5th order axial star aberration (rad)
		T c_54; 				// 5th order rosette aberration (�)
		T phi_54; 				// Azimuthal angle of 5th order rosette aberration (rad)
		T c_56; 				// 6-fold astigmatism (�)
		T phi_56; 				// Azimuthal angle of 6-fold astigmatism (rad)

		T inner_aper_ang; 		// Inner aperture (rad);
		T outer_aper_ang; 		// Outer aperture (rad);

		T ti_a; 				// Height proportion of a normalized Gaussian [0, 1]
		T ti_sigma; 			// Standard deviation of the defocus spread function for the Gaussian component:
		T ti_beta; 				// Standard deviation of the defocus spread function for the Exponential component:
		int ti_npts; 			// Number of integration points of the defocus spread function

		T ti_iehwgd; 			// e^-1 half-width value of the Gaussian distribution

		T si_a; 				// Height proportion of a normalized Gaussian [0, 1]
		T si_sigma; 			// Standard deviation of the source spread function for the Gaussian component: For parallel ilumination(�^-1); otherwise (�)
		T si_beta; 				// Standard deviation of the source spread function for the Exponential component: For parallel ilumination(�^-1); otherwise (�)
		int si_rad_npts; 		// Number of radial integration points
		int si_azm_npts; 		// Number of azimuth integration points

		T si_iehwgd; 			// e^-1 half-width value of the Gaussian distribution
		T si_theta_c; 			// divergence semi-angle (rad)

		eZero_Defocus_Type zero_defocus_type; 	// Defocus type: eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
		T zero_defocus_plane; 	// plane

		T gamma; 				// Relativistic factor
		T lambda; 				// wavelength(Angstrom)
		T lambda2; 				// wavelength(Angstrom)^2

		T c_c_10; 				// -pi*c_10*lambda
		T c_c_12; 				// -pi*c_12*lambda

		T c_c_21; 				// -2*pi*c_21*lambda^2/3
		T c_c_23; 				// -2*pi*c_23*lambda^2/3

		T c_c_30; 				// -pi*c_30*lambda^3/2
		T c_c_32; 				// -pi*c_32*lambda^3/2
		T c_c_34; 				// -pi*c_34*lambda^3/2

		T c_c_41; 				// -2*pi*c_41*lambda^4/5
		T c_c_43; 				// -2*pi*c_43*lambda^4/5
		T c_c_45; 				// -2*pi*c_45*lambda^4/5

		T c_c_50; 				// -pi*c_50*lambda^5/3
		T c_c_52; 				// -pi*c_52*lambda^5/3
		T c_c_54; 				// -pi*c_54*lambda^5/3
		T c_c_56; 				// -pi*c_56*lambda^5/3

		T g2_min; 				// inner_aper_ang/lambda
		T g2_max; 				// outer_aper_ang/lambda
		int ngxs; 				// Number of source sampling points x
		int ngys; 				// Number of source sampling points y
		T dgxs; 				// source sampling m_size;
		T dgys; 				// source sampling m_size;
		T g2_maxs; 				// q maximum square;

		Lens(): m(0), c_10(0), c_12(0), phi_12(0), c_21(0), phi_21(0), c_23(0), phi_23(0),
			c_30(0), c_32(0), phi_32(0), c_34(0), phi_34(0),
			c_41(0), phi_41(0), c_43(0), phi_43(0), c_45(0), phi_45(0),
			c_50(0), c_52(0), phi_52(0), c_54(0), phi_54(0), c_56(0), phi_56(0),
			inner_aper_ang(0), outer_aper_ang(0), ti_a(1.0), ti_sigma(0), ti_beta(0), ti_npts(0), ti_iehwgd(0),
			si_a(1.0), si_sigma(0), si_beta(0), si_rad_npts(0), si_azm_npts(0), si_iehwgd(0), si_theta_c(0),
			zero_defocus_plane(0), gamma(0), lambda(0), lambda2(0),
			g2_min(0), g2_max(0), ngxs(0), ngys(0), dgxs(0), dgys(0), g2_maxs(0),
			c_c_10(0), c_c_12(0), c_c_21(0), c_c_23(0), c_c_30(0), c_c_32(0), c_c_34(0), c_c_41(0),
			c_c_43(0), c_c_45(0), c_c_50(0), c_c_52(0), c_c_54(0), c_c_56(0) {}

		void set_input_data(T E_0, Grid_2d<T> &grid_2d)
		{
			gamma = get_gamma(E_0);

			lambda = get_lambda(E_0);
			lambda2 = pow(lambda, 2);

			c_c_10 = (isZero(c_10))?0:-c_Pi*c_10*lambda;
			c_c_12 = (isZero(c_12))?0:-c_Pi*c_12*lambda;

			c_c_21 = (isZero(c_21))?0:-2.0*c_Pi*c_21*pow(lambda, 2)/3.0;
			c_c_23 = (isZero(c_23))?0:-2.0*c_Pi*c_23*pow(lambda, 2)/3.0;

			c_c_30 = (isZero(c_30))?0:-c_Pi*c_30*pow(lambda, 3)/2.0;
			c_c_32 = (isZero(c_32))?0:-c_Pi*c_32*pow(lambda, 3)/2.0;
			c_c_34 = (isZero(c_34))?0:-c_Pi*c_34*pow(lambda, 3)/2.0;

			c_c_41 = (isZero(c_41))?0:-2.0*c_Pi*c_41*pow(lambda, 4)/5.0;
			c_c_43 = (isZero(c_43))?0:-2.0*c_Pi*c_43*pow(lambda, 4)/5.0;
			c_c_45 = (isZero(c_45))?0:-2.0*c_Pi*c_45*pow(lambda, 4)/5.0;

			c_c_50 = (isZero(c_50))?0:-c_Pi*c_50*pow(lambda, 5)/3.0;
			c_c_52 = (isZero(c_52))?0:-c_Pi*c_52*pow(lambda, 5)/3.0;
			c_c_54 = (isZero(c_54))?0:-c_Pi*c_54*pow(lambda, 5)/3.0;
			c_c_56 = (isZero(c_56))?0:-c_Pi*c_56*pow(lambda, 5)/3.0;

			g2_min = (isZero(inner_aper_ang)||(inner_aper_ang<0))?0:pow(sin(inner_aper_ang)/lambda, 2);
			g2_max = (isZero(outer_aper_ang)||(outer_aper_ang<0))?grid_2d.g2_max(): pow(sin(outer_aper_ang)/lambda, 2);

			ti_a = max(T(0), min(T(1), ti_a));
			set_ti_sigma(ti_sigma);
			ti_beta = max(T(0), ti_beta);
			ti_npts = max(1, ti_npts);

			si_a = max(T(0), min(T(1), si_a));
			set_si_sigma(si_sigma);
			si_beta = max(T(0), si_beta);
			si_rad_npts = max(1, si_rad_npts);
			si_azm_npts = max(1, si_azm_npts);

			T gmaxs = 3.5*si_sigma;
			g2_maxs = gmaxs*gmaxs;
			T dgs = gmaxs/static_cast<T>(si_rad_npts);

			int n;
			n = (dgs<grid_2d.dgx)?static_cast<int>(floor(grid_2d.dgx/dgs)+1):1;
			ngxs = static_cast<int>(floor(n*gmaxs/grid_2d.dgx)) + 1;
			dgxs = gmaxs/ngxs;

			n = (dgs<grid_2d.dgy)?static_cast<int>(floor(grid_2d.dgy/dgs)+1):1;
			ngys = static_cast<int>(floor(n*gmaxs/grid_2d.dgy)) + 1;
			dgys = gmaxs/ngys;
		}

		void set_ti_sigma(T ti_sigma_i)
		{
			ti_sigma = max(T(0), ti_sigma_i);
			ti_iehwgd = c_2i2*ti_sigma;
		}

		void set_si_sigma(T si_sigma_i)
		{
			si_sigma = max(T(0), si_sigma_i);
			si_iehwgd = c_2i2*si_sigma;
			si_theta_c = asin(si_iehwgd*lambda);
		}

		void set_defocus(T f_i)
		{
			c_10 = f_i;
			c_c_10 = (isZero(c_10))?0:-c_Pi*c_10*lambda;
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

		template <class TLens>
		void assign(TLens &lens)
		{
			m = lens.m;
			c_10	= lens.c_10;
			c_12	= lens.c_12;
			phi_12	= lens.phi_12;

			c_21	= lens.c_21;
			phi_21	= lens.phi_21;
			c_23	= lens.c_23;
			phi_23	= lens.phi_23;

			c_30	= lens.c_30;
			c_32	= lens.c_32;
			phi_32	= lens.phi_32;
			c_34	= lens.c_34;
			phi_34	= lens.phi_34;

			c_41	= lens.c_41;
			phi_41	= lens.phi_41;
			c_43	= lens.c_43;
			phi_43	= lens.phi_43;
			c_45	= lens.c_45;
			phi_45	= lens.phi_45;

			c_50	= lens.c_50;
			c_52	= lens.c_52;
			phi_52	= lens.phi_52;
			c_54	= lens.c_54;
			phi_54	= lens.phi_54;
			c_56	= lens.c_56;
			phi_56	= lens.phi_56;

			inner_aper_ang = lens.inner_aper_ang;
			outer_aper_ang = lens.outer_aper_ang;

			ti_sigma = lens.ti_sigma;
			ti_npts = lens.ti_npts;
			ti_iehwgd = lens.ti_iehwgd;

			si_sigma = lens.si_sigma;
			si_rad_npts = lens.si_rad_npts;
			si_iehwgd = lens.si_iehwgd;
			si_theta_c = lens.si_theta_c;

			gamma = lens.gamma;
			lambda = lens.lambda;
			lambda2 = lens.lambda2;

			c_c_10	= lens.c_c_10;
			c_c_12	= lens.c_c_12;

			c_c_21	= lens.c_c_21;
			c_c_23	= lens.c_c_23;

			c_c_30	= lens.c_c_30;
			c_c_32	= lens.c_c_32;
			c_c_34	= lens.c_c_34;

			c_c_41	= lens.c_c_41;
			c_c_43	= lens.c_c_43;
			c_c_45	= lens.c_c_45;

			c_c_50	= lens.c_c_50;
			c_c_52	= lens.c_c_52;
			c_c_54	= lens.c_c_54;
			c_c_56	= lens.c_c_56;

			g2_min = lens.g2_min;
			g2_max = lens.g2_max;

			ngxs = lens.ngxs;
			ngys = lens.ngys;
			dgxs = lens.dgxs;
			dgys = lens.dgys;
			g2_maxs = lens.g2_maxs;
		}

		template <class TLens>
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

		/************************************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		T ti_sigma2() const
		{
			return ti_sigma*ti_sigma;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T ti_beta2() const
		{
			return ti_beta*ti_beta;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T si_sigma2() const
		{
			return si_sigma*si_sigma;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T si_beta2() const
		{
			return si_beta*si_beta;
		}

		/************************************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		bool is_phi_required() const
		{
			auto bb = nonZero(m)||nonZero(c_12)||nonZero(c_21)||nonZero(c_23)||nonZero(c_32)||nonZero(c_34);
			bb = bb||nonZero(c_41)||nonZero(c_43)||nonZero(c_45)||nonZero(c_52)||nonZero(c_54)||nonZero(c_56);
			return bb;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T eval_m(const T &phi) const
		{
			return (nonZero(phi))?(m*phi):0;
		}

		/************************************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_10(const T &g2) const
		{
			return (nonZero(c_10))?(c_c_10*g2):0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_12(const T &g2, const T &phi) const
		{
			return (nonZero(c_12))?(c_c_12*g2*sin(2*(phi-phi_12))):0;
		}

		/************************************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_21(const T &g3, const T &phi) const
		{
			return (nonZero(c_21))?(c_c_21*g3*sin(phi-phi_21)):0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_23(const T &g3, const T &phi) const
		{
			return (nonZero(c_23))?(c_c_23*g3*sin(3*(phi-phi_23))):0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_21_c_23(const T &g3, const T &phi) const
		{
			return (eval_c_21(g3, phi) + eval_c_23(g3, phi));
		}
		/************************************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_30(const T &g4) const
		{
			return (nonZero(c_30))?(c_c_30*g4):0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_32(const T &g4, const T &phi) const
		{
			return (nonZero(c_32))?(c_c_32*g4*sin(2*(phi-phi_32))):0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_34(const T &g4, const T &phi) const
		{
			return (nonZero(c_34))?(c_c_34*g4*sin(4*(phi-phi_34))):0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_32_c_34(const T &g4, const T &phi) const
		{
			return (eval_c_32(g4, phi) + eval_c_34(g4, phi));
		}
		/************************************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_41(const T &g5, const T &phi) const
		{
			return (nonZero(c_41))?(c_c_41*g5*sin(phi-phi_41)):0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_43(const T &g5, const T &phi) const
		{
			return (nonZero(c_43))?(c_c_43*g5*sin(3*(phi-phi_43))):0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_45(const T &g5, const T &phi) const
		{
			return (nonZero(c_45))?(c_c_45*g5*sin(5*(phi-phi_45))):0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_41_c_43_c_45(const T &g5, const T &phi) const
		{
			return (eval_c_41(g5, phi) + eval_c_43(g5, phi) + eval_c_45(g5, phi));
		}
		/************************************************************************/
		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_50(const T &g6) const
		{
			return (nonZero(c_50))?(c_c_50*g6):0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_52(const T &g6, const T &phi) const
		{
			return (nonZero(c_52))?(c_c_52*g6*sin(2*(phi-phi_52))):0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_54(const T &g6, const T &phi) const
		{
			return (nonZero(c_54))?(c_c_54*g6*sin(4*(phi-phi_54))):0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_56(const T &g6, const T &phi) const
		{
			return (nonZero(c_56))?(c_c_56*g6*sin(6*(phi-phi_56))):0;
		}

		DEVICE_CALLABLE FORCE_INLINE
		T eval_c_52_c_54_c_56(const T &g6, const T &phi) const
		{
			return (eval_c_52(g6, phi) + eval_c_54(g6, phi) + eval_c_56(g6, phi));
		}
	};

	/****************************EELS***************************/
	template <class T>
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

		template <class TEELS>
		void assign(TEELS &eels)
		{
			set_input_data(eels.space, eels.E_0, eels.E_loss, eels.m_selection, eels.collection_angle, eels.channelling_type, eels.Z);
			factor = eels.factor;
			x = eels.x;
			y = eels.y;
			occ = eels.occ;
		}

		template <class TEELS>
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

	/**************************atoms for Suppos**************************/
	template <class T>
	class Atom_Data_Sp{
		public:
			using value_type = T;
			using size_type = std::size_t;

			Atom_Data_Sp(): l_x(0), l_y(0),
				x_min(0), x_max(0),
				y_min(0), y_max(0),
				a_min(0), a_max(0),
				sigma_min(0), sigma_max(0),
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
			}

			// reserve
			void reserve(const size_type &new_size)
			{
				x.reserve(new_size);
				y.reserve(new_size);
				a.reserve(new_size);
				sigma.reserve(new_size);
			}

			// reserve
			void push_back(T x_i, T y_i, T a_i, T sigma_i)
			{
				x.push_back(x_i);
				y.push_back(y_i);
				a.push_back(a_i);
				sigma.push_back(sigma_i);
			}

			// set atoms
			void set_atoms(const size_type &nr_atoms_i, const size_type &nc_atoms_i, double *atoms_i, T l_x_i = 0, T l_y_i = 0, bool pbc_xy_i = false)
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
					if((!pbc_xy_i)||((atom.x<lx_b) && (atom.y<ly_b)))
					{
						x[j] = atom.x; 				// x-position
						y[j] = atom.y; 				// y-position
						a[j] = atom.a; 				// height
						sigma[j] = atom.sigma;		// standard deviation
						j++;
					}
				}

				resize(j);

				get_statistic();
			}

			// set atoms
			void set_atoms(const Atom_Data_Sp<T> &atoms, bool pbc_xy_i = false)
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
					if((!pbc_xy_i)||((atoms.x[i]<lx_b) && (atoms.y[i]<ly_b)))
					{
						x[j] = atoms.x[i]; 				// x-position
						y[j] = atoms.y[i]; 				// y-position
						a[j] = atoms.a[i];				// height
						sigma[j] = atoms.sigma[i];		// standard deviation
						j++;
					}
				}

				resize(j);

				get_statistic();
			}

			// get statistic
			void get_statistic()
			{
				if(empty())
				{
					return;
				}

				x_min = x_max = x[0];
				y_min = y_max = y[0];
				a_min = a_max = a[0];
				sigma_min = sigma_max = sigma[0];

				x_mean = y_mean = 0.0;
				x_std = y_std = 0.0;

				for(auto iatoms = 0; iatoms < size(); iatoms++)
				{
					x_min = min(x[iatoms], x_min);
					x_max = max(x[iatoms], x_max);

					y_min = min(y[iatoms], y_min);
					y_max = max(y[iatoms], y_max);

					a_min = min(a[iatoms], a_min);
					a_max = max(a[iatoms], a_max);

					sigma_min = min(sigma[iatoms], sigma_min);
					sigma_max = max(sigma[iatoms], sigma_max);

					x_mean += x[iatoms];
					y_mean += y[iatoms];

					x_std += x[iatoms]*x[iatoms];
					y_std += y[iatoms]*y[iatoms];
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

			T x_min;
			T x_max;

			T y_min;
			T y_max;

			T a_min;
			T a_max;

			T sigma_min;
			T sigma_max;

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

				Atom():x(0), y(0), a(0), sigma(0){};
			};

			template <class TIn>
			Atom read_atom(const int &nr, const int &nc, TIn *atoms, const int &iatoms)
			{
				Atom atom;
				atom.x = atoms[0*nr + iatoms]; 						// x-position
				atom.y = atoms[1*nr + iatoms]; 						// y-position
				atom.a = (nc>2)?atoms[2*nr + iatoms]:1.0;			// height
				atom.sigma = (nc>3)?atoms[3*nr + iatoms]:1.0;		// standard deviation

				return atom;
			}
	};

	template <class T, eDevice dev>
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

			int ix_0;
			int nx;
			int iy_0;
			int ny;

			T a;
			T alpha;
			T dtR;
			T R2_tap;
			T tap_cf;

			int *iv;
			T *v;

			Atom_Sp(): x(0), y(0), occ(1), R2_max(0), R2(nullptr), c3(nullptr), c2(nullptr),
			c1(nullptr), c0(nullptr), ix_0(0), nx(1), iy_0(0), ny(1), a(0), alpha(0),
			dtR(0), R2_tap(0), tap_cf(0), iv(nullptr), v(nullptr){}

			inline
			void set_ix0_ixn(const Grid_2d<T> &grid_2d, const T &R_max)
			{
				get_bn(x, grid_2d.nx, grid_2d.dRx, R_max, grid_2d.pbc_xy, ix_0, nx);
			}

			inline
			void set_iy0_iyn(const Grid_2d<T> &grid_2d, const T &R_max)
			{
				get_bn(y, grid_2d.ny, grid_2d.dRy, R_max, grid_2d.pbc_xy, iy_0, ny);
			}

			inline
			Grid_BT get_eval_cubic_poly_gridBT()
			{
				Grid_BT grid_bt;
				grid_bt.Blk = dim3((ny+c_thrnxny-1)/c_thrnxny, (nx+c_thrnxny-1)/c_thrnxny);
				grid_bt.Thr = dim3(c_thrnxny, c_thrnxny);
				return grid_bt;
			}
	};

	/**************************atoms for Superpositon************************/
	template <class T>
	struct Gauss_Sp
	{
		public:
			using value_type = T;

			T x;
			T y;
			T R2_max;

			int ix_0;
			int iy_0;
			int nx;
			int ny;

			T a;
			T alpha;
			T R2_tap;
			T tap_cf;

			int *iv;
			T *v;

			Gauss_Sp(): x(0), y(0), R2_max(0), ix_0(0), iy_0(0), nx(1), ny(1),
			a(0), alpha(0), R2_tap(0), tap_cf(0), iv(nullptr), v(nullptr){}


			DEVICE_CALLABLE FORCE_INLINE
			T operator()(const T &R2) const
			{
				const T tap = host_device_detail::tapering(R2_tap, tap_cf, R2);
				const T y = a*exp(-alpha*R2)*tap;
				return y;
			}

			inline
			void set_ix0_ixn(const Grid_2d<T> &grid_2d, const T &R_max)
			{
				get_bn(x, grid_2d.nx, grid_2d.dRx, R_max, grid_2d.pbc_xy, ix_0, nx);
			}

			inline
			void set_iy0_iyn(const Grid_2d<T> &grid_2d, const T &R_max)
			{
				get_bn(y, grid_2d.ny, grid_2d.dRy, R_max, grid_2d.pbc_xy, iy_0, ny);
			}
	};

	template <class T>
	class Input_Gauss_Spt
	{
		public:
			using value_type = T;

			T ff_sigma;
			Grid_2d<T> grid_2d;
			Atom_Data_Sp<T> atoms;

			Input_Gauss_Spt(): ff_sigma(3){};

			T alpha(const int &iatoms) const
			{
				return 0.5/pow(atoms.sigma[iatoms], 2);
			}

			T R_max(const int &iatoms) const
			{
				return ff_sigma*atoms.sigma[iatoms];
			}

			// maximum number of pixels
			int get_nv()
			{
				auto l_x = atoms.l_x + 2*ff_sigma*atoms.sigma_max;
				auto l_y = atoms.l_y + 2*ff_sigma*atoms.sigma_max;
				return max(grid_2d.nx_dRx(l_x), grid_2d.ny_dRy(l_y));
			}
	};

	/***********************atoms for simulated annealing***********************/
	template <class T>
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

			template <class TAtom_SA>
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
			void set_atoms(const size_type &natoms_i, double *atoms_i, double *atoms_min_i, double *atoms_max_i)
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
			void set_atoms(const size_type &natoms_i, double *atoms_i, r3d<T> d_i)
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

			inline
			T norm(const int &iatoms, const r3d<T> &r)
			{
				auto rd = r_n[iatoms]-r;
				return mt::norm(rd);
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

	template <class T>
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

			int ix_0;
			int ixn;
			int iy_0;
			int iyn;

			int *iv;
			T *v;

			Atom_Sa(): x(0), y(0), R2_max(0), R2(nullptr), c3(nullptr), c2(nullptr), c1(nullptr),
			c0(nullptr), ix_0(1), ixn(0), iy_0(0), iyn(0), iv(nullptr), v(nullptr){}

			inline
			void set_ix0_ixn(const Grid_2d<T> &grid_2d, const T &R_max)
			{
				get_bn(x, grid_2d.nx, grid_2d.dRx, R_max, grid_2d.pbc_xy, ix_0, ixn);
			}

			inline
			void set_iy0_iyn(const Grid_2d<T> &grid_2d, const T &R_max)
			{
				get_bn(y, grid_2d.ny, grid_2d.dRy, R_max, grid_2d.pbc_xy, iy_0, iyn);
			}

			inline
			Grid_BT get_eval_cubic_poly_gridBT()
			{
				Grid_BT grid_bt;
				grid_bt.Blk = dim3((iyn+c_thrnxny-1)/c_thrnxny, (ixn+c_thrnxny-1)/c_thrnxny);
				grid_bt.Thr = dim3(c_thrnxny, c_thrnxny);
				return grid_bt;
			}
	};

	/*****************************Potential**************************/
	template <class T>
	struct Atom_Vp
	{
	public:
		using value_type = T;

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
		int ix_0;
		int nx;
		int iy_0;
		int ny;

		T R2_tap;
		T tap_cf;

		int *iv;
		T *v;

		Atom_Vp(): charge(0), x(0), y(0), z0h(0), zeh(0), split(false), occ(1), R2_min(0), R2_max(0),
			R2(nullptr), cl(nullptr), cnl(nullptr), c3(nullptr), c2(nullptr), c1(nullptr),
			c0(nullptr), ix_0(1), nx(0), iy_0(0), ny(0), R2_tap(0), tap_cf(0), iv(nullptr), v(nullptr){}

		void set_ix0_ixn(const Grid_2d<T> &grid_2d, const T &R_max)
		{
			get_bn(x, grid_2d.nx, grid_2d.dRx, R_max, grid_2d.pbc_xy, ix_0, nx);
		}

		void set_iy0_iyn(const Grid_2d<T> &grid_2d, const T &R_max)
		{
			get_bn(y, grid_2d.ny, grid_2d.dRy, R_max, grid_2d.pbc_xy, iy_0, ny);
		}
	};

	/*****************************Atomic Coefficients**************************/
	template <class T, eDevice dev>
	struct Atom_Coef
	{
		using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = dev;

		Atom_Coef(): charge(0), tag(0), R_min(0), R_max(0), R_tap(0), tap_cf(0){}

		template <class TAtom_Coef>
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

		template <class TAtom_Coef>
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
	template <class T, eDevice dev>
	struct Atom_Type
	{
		using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = dev;

		Atom_Type(): Z(0), m(0), A(0), rn_e(0), rn_c(0), ra_e(0), ra_c(0){}

		template <class TAtom_Type>
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

		template <class TAtom_Type>
		Atom_Type<T, dev>& operator=(TAtom_Type &atom_type)
		{
			assign(atom_type);
			return *this;
		}

		int check_charge(const int &charge) const
		{
			for(auto i= 0; i<coef.size(); i++)
			{
				if(coef[i].charge == charge)
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
				if(coef[i].charge == charge)
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
	template <class T>
	struct Scanning
	{
		public:
			using value_type = T;
			using size_type = std::size_t;

			eScanning_Type type;			// 1: Line, 2: Area,
			eGrid_Type grid_type;			// 1: regular, 2: quadratic
			bool pbc;						// periodic boundary conditions
			bool spxs;						// square pixel size
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

			Scanning(): type(eST_Line), grid_type(eGT_Regular), pbc(false), spxs(true), ns(1),
				nx(0), dRx(0), dRy(0), ny(0), x0(0), y0(0), xe(0), ye(0) {};

			template <class TScanning>
			void assign(TScanning &scanning)
			{
				type = scanning.type;
				grid_type = scanning.grid_type;
				pbc = scanning.pbc;
				spxs = scanning.spxs;
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

			template <class TScanning>
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
				spxs = true;
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

						if (!spxs)
						{
							dRy = yu/((pbc)?ny:(ny-1));
						}
					}
					else
					{
						dRy = yu/((pbc)?ns:(ns-1));
						dRx = std::copysign(dRy, xu);
						nx = int(floor(xu/dRx+Epsilon<T>::rel+0.5));
						nx += (pbc)?0:1;

						if (!spxs)
						{
							dRx = xu/((pbc)?nx:(nx-1));
						}
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
	template <class T>
	struct In_Rad_Schr
	{
		T E_0; 					// Acceleration Voltage
		ePotential_Type potential_type; 			// Parameterization type
		int n; 					// Principal quantum number
		int nr; 				// Number of grid points
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

	/*******************forward declarations********************/
	bool is_gpu_available();

	int number_of_gpu_available();

	/************************Host device configuration************************/
	class System_Configuration
	{
		public:
			ePrecision precision;
			eDevice device; 									// eP_float = 1, eP_double = 2
			int cpu_ncores; 									// Number of Cores CPU
			int cpu_nthread; 									// Number of threads
			int gpu_device; 									// GPU device
			int gpu_nstream; 									// Number of streams

			int nstream;
			bool active;

			System_Configuration(): precision(eP_double), device(e_host), cpu_ncores(1),
				cpu_nthread(4), gpu_device(0), gpu_nstream(8), nstream(1), active(true){};

			void validate_parameters()
			{
				// check precision
				if(!(is_float() || is_double()))
				{
					precision = eP_float;
				}

				// check cpu or gpu
				if(!(is_host() || is_device()))
				{
					device = e_host;
				}
				if(is_device())
				{
					#ifdef __CUDACC__
						if(!is_gpu_available())
						{
							device = mt::e_host;
							n_gpu = 0;
						}
						else
						{
							n_gpu = number_of_gpu_available();
							gpu_device = min(max(0, gpu_device), n_gpu-1);
						}
					#endif
				}

				cpu_nthread = max(1, cpu_nthread);
				gpu_nstream = max(1, gpu_nstream);
				nstream = (is_host())?cpu_nthread:gpu_nstream;
			}

			void set_device()
			{
				if(is_device())
				{
					#ifdef __CUDACC__
						cudaSetDevice(gpu_device);
					#endif
				}
				else
				{

					device = mt::e_host;
				}
			}

 			int get_device()
			{
				int idx_dev = -1;
				if(is_device())
				{
					#ifdef __CUDACC__
						cudaGetDevice(&idx_dev);
					#endif
				}

				return idx_dev;
			}

			bool is_host() const
			{
				return device == mt::e_host;
			}

			bool is_device() const
			{
				return device == mt::e_device;
			}

			bool is_float() const
			{
				return precision == mt::eP_float;
			}

			bool is_double() const
			{
				return precision == mt::eP_double;
			}

			bool is_float_host() const
			{
				return is_float() && is_host();
			}

			bool is_double_host() const
			{
				return is_double() && is_host();
			}

			bool is_float_device() const
			{
				return is_float() && is_device();
			}

			bool is_double_device() const
			{
				return is_double() && is_device();
			}

		private:
			int n_gpu;
	};

	/************************regions*************************/
	template <class TVector>
	struct Region
	{
		public:
			using T = typename TVector::value_type;
			using size_type = std::size_t;

			TVector Rx;
			TVector Ry;
			TVector R2;
			TVector Ixy;

			T R_max;
			T Rx_sf;
			T Ry_sf;
			T Rxy_sc;

			T Ixy_sf;
			T Ixy_sc;

			Region(): Rx_sf(0), Ry_sf(0), Rxy_sc(1), Ixy_sf(0), Ixy_sc(1)
			{
			}

			void clear()
			{
				Rx.clear();
				Ry.clear();
				R2.clear();
				Ixy.clear();
			}

			void reserve(const size_type &new_size)
			{
				Rx.reserve(new_size);
				Ry.reserve(new_size);
				R2.reserve(new_size);
				Ixy.reserve(new_size);
			}

			void shrink_to_fit()
			{
				Rx.shrink_to_fit();
				Ry.shrink_to_fit();
				R2.shrink_to_fit();
				Ixy.shrink_to_fit();
			}

			TVector sft_Ixy(T bg)
			{
				TVector Ixy_s;
				Ixy_s.reserve(Ixy.size());

				for(auto ixy=0; ixy<Ixy.size(); ixy++)
				{
					Ixy_s.push_back(Ixy[ixy]-bg);
				}
				return Ixy_s;
			}

			TVector sub_region_to_Ixy(Grid_2d<T> &grid_2d, TVector &Im_s, T x, T y)
			{
				TVector v = Ixy;

				T R2_max = pow(R_max, 2);

				r2d<T> p(x, y);
				auto range = grid_2d.index_range(p, R_max);
				int iv = 0;
				for (auto ix = range.ix_0; ix < range.ix_e; ix++)
				{
					for (auto iy = range.iy_0; iy < range.iy_e; iy++)
					{
						T r2 = grid_2d.R2(ix, iy, p.x, p.y);
						if (r2 < R2_max)
						{
							v[iv++] -= Im_s[grid_2d.ind_col(ix, iy)]/Ixy_sc;
						}
					}
				}
				return v;
			}

			TVector sft_Ixy(Grid_2d<T> &grid_2d, TVector &Im_s, T x, T y, T a, T s)
			{
				TVector v = sub_region_to_Ixy(grid_2d, Im_s, x*Rxy_sc+Rx_sf, y*Rxy_sc+Ry_sf);

				T alpha = 0.5/pow(s, 2);
				T r2_l = pow(4.0*s, 2);
				for(auto im = 0; im < v.size(); im++)
				{
					T rx = Rx[im]-x;
					T ry = Ry[im]-y;
					T r2 = rx*rx+ry*ry;
					if(r2<r2_l)
					{
						v[im] += a*exp(-alpha*r2);
					}
				}

				return v;
			}

			TVector sft_Ixy(Grid_2d<T> &grid_2d, TVector &Im_s, TVector &x, TVector &y, TVector &A, TVector &S)
			{
				TVector v = sub_region_to_Ixy(grid_2d, Im_s, x[0]*Rxy_sc+Rx_sf, y[0]*Rxy_sc+Ry_sf);

				for(auto ip = 0; ip < x.size(); ip++)
				{
					T a = A[ip];
					T b = S[ip];
					T alpha = 0.5/pow(b, 2);
					T r2_l = pow(4.0*b, 2);
					for(auto im = 0; im < v.size(); im++)
					{
						T rx = Rx[im]-x[ip];
						T ry = Ry[im]-y[ip];
						T r2 = rx*rx+ry*ry;
						if(r2<r2_l)
						{
							v[im] += a*exp(-alpha*r2);
						}
					}
				}

				return v;
			}

			T sft_Rx(T x) const
			{
				return (x-Rx_sf)/Rxy_sc;
			}

			T sft_Ry(T y) const
			{
				return (y-Ry_sf)/Rxy_sc;
			}

			r2d<T> sft_x_y(T x, T y) const
			{
				x = sft_Rx(x);
				y = sft_Ry(y);
				return r2d<T>(x, y);
			}
	};

	/***********************************************/
	template <class T>
	class Gauss_1d{
		public:
			using value_type = T;

			T sigma;
			T alpha;
			T Rx2_l;
			T x_c;
			T k;
			bool b_norm;

			Gauss_1d(): sigma(1), Rx2_l(0),
			x_c(0), k(0), b_norm(false){}

			Gauss_1d(T x_c_i, T sigma_i)
			{
				set_input_data(x_c_i, sigma_i);
			}

			Gauss_1d(Border_1d<T> &bd_i, T sigma_i, bool b_norm_i = true)
			{
				set_input_data(bd_i, sigma_i, b_norm_i);
			}

			inline
			void set_input_data(T x_c_i, T sigma_i)
			{
				sigma = sigma_i;
				alpha = 0.5/(sigma*sigma);
				Rx2_l = 1e+6*sigma;
				x_c = x_c_i;
				k = 0;
				b_norm = false;
			}

			inline
			void set_input_data(Border_1d<T> &bd_i, T sigma_i, bool b_norm_i = true)
			{
				sigma = sigma_i;
				alpha = 0.5/(sigma*sigma);
				Rx2_l = pow(bd_i.lx_wb()/2, 2);
				x_c = bd_i.x_c();
				k = exp(-alpha*Rx2_l);
				b_norm = b_norm_i;
			}

			DEVICE_CALLABLE FORCE_INLINE
			T eval_norm(const T &Rx2) const
			{
				return (Rx2<=Rx2_l)?((b_norm)?::fmax(T(0), (this->operator()(Rx2)-k)/(1-k)):this->operator()(Rx2)):0;
			}

			DEVICE_CALLABLE FORCE_INLINE
			T operator()(const T &Rx2) const
			{
				return exp(-alpha*Rx2);
			}
	};

	template <class T>
	class Gauss_2d{
		public:
			using value_type = T;

			T sigma;
			T alpha;
			T R2_l;
			T x_c;
			T y_c;
			T k;
			bool b_norm;

			Gauss_2d(): sigma(1), R2_l(0),
			x_c(0), y_c(0), k(0), b_norm(false){}

			Gauss_2d(T x_c_i, T y_c_i, T sigma_i)
			{
				set_input_data(x_c_i, y_c_i, sigma_i);
			}

			Gauss_2d(Border_2d<T> &bd_i, T sigma_i, bool b_norm_i = true)
			{
				set_input_data(bd_i, sigma_i, b_norm_i);
			}

			inline
			void set_input_data(T x_c_i, T y_c_i, T sigma_i)
			{
				sigma = sigma_i;
				alpha = 0.5/(sigma*sigma);
				R2_l = 1e+6*sigma;
				x_c = x_c_i;
				y_c = y_c_i;
				k = 0;
				b_norm = false;
			}

			inline
			void set_input_data(Border_2d<T> &bd_i, T sigma_i, bool b_norm_i = true)
			{
				sigma = sigma_i;
				alpha = 0.5/(sigma*sigma);
				R2_l = pow(min(bd_i.lx_wb(), bd_i.ly_wb())/2, 2);
				x_c = bd_i.x_c();
				y_c = bd_i.y_c();
				k = exp(-alpha*R2_l);
				b_norm = b_norm_i;
			}

			DEVICE_CALLABLE FORCE_INLINE
			T eval_norm(const T &R2) const
			{
				return (R2<=R2_l)?((b_norm)?::fmax(T(0), (this->operator()(R2)-k)/(1-k)):this->operator()(R2)):0;
			}

			DEVICE_CALLABLE FORCE_INLINE
			T operator()(const T &R2) const
			{
				return exp(-alpha*R2);
			}
	};

	/***********************************************/
	template <class T>
	class Hanning_1d{
		public:
			using value_type = T;

			T lx;
			T cx;
			T Rx_l;
			T x_c;
			T k;
			bool b_norm;

			Hanning_1d(): lx(1), cx(0), Rx_l(0),
			x_c(0), k(0), b_norm(false){}

			Hanning_1d(T x_c_i, T lx_i, T k_i)
			{
				set_input_data(x_c_i, lx_i, k_i);
			}

			Hanning_1d(Border_1d<T> &bd_i, T k_i, bool b_norm_i = true)
			{
				set_input_data(bd_i, k_i, b_norm_i);
			}

			inline
			void set_input_data(T x_c_i, T lx_i, T k_i)
			{
				lx = lx_i;
				cx = c_Pi/lx;
				Rx_l = lx/2;
				x_c = x_c_i;
				k = (k_i>1)?1.0/k_i:pow(sin(0.5*c_Pi*k_i), 2);
				b_norm = false;
			}

			inline
			void set_input_data(Border_1d<T> &bd_i, T k_i, bool b_norm_i = true)
			{
				lx = bd_i.lx_wb();
				cx = c_Pi/lx;
				Rx_l = lx/2;
				x_c = bd_i.x_c();
				k = (k_i>1)?1.0/k_i:pow(sin(0.5*c_Pi*k_i), 2);
				b_norm = b_norm_i;
			}

			DEVICE_CALLABLE FORCE_INLINE
			T eval_norm(const T &Rx) const
			{
				T v = (fabs(Rx)<=Rx_l)?((b_norm)?::fmax(T(0), this->operator()(Rx)):this->operator()(Rx)):0;
				return (v>k)?1.0:v/k;
			}

			DEVICE_CALLABLE FORCE_INLINE
			T operator()(const T &Rx) const
			{
				return ::square(cos(cx*Rx));
			}
	};

	template <class T>
	class Hanning_2d{
		public:
			using value_type = T;

			T lx;
			T ly;
			T cxy;
			T R_l;
			T x_c;
			T y_c;
			T k;
			bool b_norm;

			Hanning_2d(): lx(1), ly(1), cxy(0), R_l(0),
			x_c(0), y_c(0), k(0), b_norm(false){}

			Hanning_2d(T x_c_i, T y_c_i, T ly_i, T lx_i, T k_i)
			{
				set_input_data(x_c_i, y_c_i, lx_i, ly_i, k_i);
			}

			Hanning_2d(Border_2d<T> &bd_i, T k_i, bool b_norm_i = true)
			{
				set_input_data(bd_i, k_i, b_norm_i);
			}

			inline
			void set_input_data(T x_c_i, T y_c_i, T ly_i, T lx_i, T k_i)
			{
				lx = lx_i;
				ly = lx_i;
				cxy = c_Pi/min(lx, ly);
				R_l = min(lx, ly)/2;
				x_c = x_c_i;
				y_c = y_c_i;
				k = (k_i>1)?1.0/k_i:pow(sin(0.5*c_Pi*k_i), 2);
				b_norm = false;
			}

			inline
			void set_input_data(Border_2d<T> &bd_i, T k_i, bool b_norm_i = true)
			{
				lx = bd_i.lx_wb();
				ly = bd_i.ly_wb();
				cxy = c_Pi/min(lx, ly);
				R_l = min(lx, ly)/2;
				x_c = bd_i.x_c();
				y_c = bd_i.y_c();
				k = (k_i>1)?1.0/k_i:pow(sin(0.5*c_Pi*k_i), 2);
				b_norm = b_norm_i;
			}

			DEVICE_CALLABLE FORCE_INLINE
			T eval_norm(const T &R) const
			{
				T v = (fabs(R)<=R_l)?((b_norm)?::fmax(T(0), this->operator()(R)):this->operator()(R)):0;
				return (v>k)?1.0:v/k;
			}

			DEVICE_CALLABLE FORCE_INLINE
			T operator()(const T &R) const
			{
				return ::square(cos(cxy*R));
			}
	};

	/***********************************************/
	template <class T>
	class Butterworth_1d{
		public:
			using value_type = T;

			T radius;
			T R02;
			T Rx2_l;
			T x_c;
			T k;
			int n;
			bool b_norm;

			Butterworth_1d(): radius(0), R02(0), Rx2_l(0),
			x_c(0), k(0), n(16), b_norm(false){}

			Butterworth_1d(T x_c_i, T radius_i, int n_i)
			{
				set_input_data(x_c_i, radius_i, n_i);
			}

			Butterworth_1d(Border_1d<T> &bd_i, T radius_i, int n_i, bool b_norm_i = true)
			{
				set_input_data(bd_i, radius_i, n_i, b_norm_i);
			}

			inline
			void set_input_data(T x_c_i, T radius_i, int n_i)
			{
				n = n_i;
				radius = radius_i;
				R02 = ::square(radius);
				Rx2_l = 1e+6*radius;
				x_c = x_c_i;
				k = 0;
				b_norm = false;
			}

			inline
			void set_input_data(Border_1d<T> &bd_i, T radius_i, int n_i, bool b_norm_i = true)
			{
				n = n_i;
				radius = radius_i;
				R02 = pow(radius, 2);
				Rx2_l = pow(bd_i.lx_wb()/2, 2);
				x_c = bd_i.x_c();
				k = 1.0/(1.0+pow(Rx2_l/R02, n));
				b_norm = b_norm_i;
			}

			/*********************************************************/
			DEVICE_CALLABLE FORCE_INLINE
			T eval_norm(const T &Rx2) const
			{
				return (Rx2<=Rx2_l)?((b_norm)?::fmax(T(0), (this->operator()(Rx2)-k)/(1-k)):this->operator()(Rx2)):0;
			}

			DEVICE_CALLABLE FORCE_INLINE
			T operator()(const T &Rx2) const
			{
				return 1.0/(1.0+pow(Rx2/R02, n));
			}
	};

	template <class T>
	class Butterworth_2d{
		public:
			using value_type = T;

			T radius;
			T R02;
			T R2_l;
			T x_c;
			T y_c;
			T k;
			int n;
			bool b_norm;

			Butterworth_2d(): radius(0), R02(0), R2_l(0),
			x_c(0), y_c(0), k(0), n(16), b_norm(false){}

			Butterworth_2d(T x_c_i, T y_c_i, T radius_i, int n_i)
			{
				set_input_data(x_c_i, y_c_i, radius_i, n_i);
			}

			Butterworth_2d(Border_2d<T> &bd_i, T radius_i, int n_i, bool b_norm_i = true)
			{
				set_input_data(bd_i, radius_i, n_i, b_norm_i);
			}

			inline
			void set_input_data(T x_c_i, T y_c_i, T radius_i, int n_i)
			{
				n = n_i;
				radius = radius_i;
				R02 = pow(radius, 2);
				R2_l = 1e+6*radius;
				x_c = x_c_i;
				y_c = y_c_i;
				k = 0;
				b_norm = false;
			}

			inline
			void set_input_data(Border_2d<T> &bd_i, T radius_i, int n_i, bool b_norm_i = true)
			{
				n = n_i;
				radius = radius_i;
				R02 = pow(radius, 2);
				R2_l = pow(min(bd_i.lx_wb(), bd_i.ly_wb())/2, 2);
				x_c = bd_i.x_c();
				y_c = bd_i.y_c();
				k = 1.0/(1.0+pow(R2_l/R02, n));
				b_norm = b_norm_i;
			}

			/*********************************************************/
			DEVICE_CALLABLE FORCE_INLINE
			T eval_norm(const T &R2) const
			{
				return (R2<=R2_l)?((b_norm)?::fmax(T(0), (this->operator()(R2)-k)/(1-k)):this->operator()(R2)):0;
			}

			DEVICE_CALLABLE FORCE_INLINE
			T operator()(const T &R2) const
			{
				return (1.0/(1.0+pow(R2/R02, n)));
			}
	};

	/********Affine parameters shx and scy*********/
	template <class T>
	struct Afp_2
	{
		r2d<T> f;
		r2d<T> ds;
		T chi2;

		Afp_2():f(T(0), T(1)), ds(T(0), T(0)), chi2(0){}

		Afp_2(const r2d<T> &f_i, const r2d<T> &ds_i, const T &chi2_i):f(f_i), ds(ds_i), chi2(chi2_i){}

		template <class TAfp_2>
		Afp_2<T>& operator=(TAfp_2 &afp_2)
		{
			f = afp_2.f;
			ds = afp_2.ds;
			chi2 = afp_2.chi2;

			return *this;
		}
	};

} // namespace mt

#endif