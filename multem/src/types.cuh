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

#include <multem/constants.h>
#include "math.cuh"
#include "lin_alg_def.cuh"
#include "safe_types.cuh"

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

	/************************Host device configuration************************/


  System_Configuration::System_Configuration()
    : precision(eP_double), 
      device(e_host), 
      cpu_ncores(1),
			cpu_nthread(1), 
      gpu_device(0), 
      gpu_nstream(1), 
      nstream(1), 
      active(true) {};

  void System_Configuration::validate_parameters()
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

  void System_Configuration::set_device()
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

  int System_Configuration::get_device()
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

  bool System_Configuration::is_host() const
  {
    return device == mt::e_host;
  }

  bool System_Configuration::is_device() const
  {
    return device == mt::e_device;
  }

  bool System_Configuration::is_float() const
  {
    return precision == mt::eP_float;
  }

  bool System_Configuration::is_double() const
  {
    return precision == mt::eP_double;
  }

  bool System_Configuration::is_float_host() const
  {
    return is_float() && is_host();
  }

  bool System_Configuration::is_double_host() const
  {
    return is_double() && is_host();
  }

  bool System_Configuration::is_float_device() const
  {
    return is_float() && is_device();
  }

  bool System_Configuration::is_double_device() const
  {
    return is_double() && is_device();
  }

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


  template <typename T>
  void Lens<T>::set_input_data(T E_0, Grid_2d<T> &grid_2d)
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

  template <typename T>
	void EELS<T>::set_input_data(eSpace space_i, T E_0_i, T E_loss_i, int m_selection_i, T collection_angle_i, eChannelling_Type channelling_type_i, int Z_i)
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

  template <typename T>
  struct Detector<T, e_device> {
		
    using value_type = T;
		using size_type = std::size_t;

		static const eDevice device = e_device;

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
		Detector<T, e_device>& operator=(TDetector &detector)
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
    std::vector<T> g_inner;				// Inner aperture Ang^-1
    std::vector<T> g_outer;				// Outer aperture Ang^-1
    std::vector<Vector<T, e_device>> fx;		// radial sensitivity value
    std::vector<Vector<T, e_device>> fR;		// 2D sensitivity value
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

    // Need this to stop compiler warnings about calling __host__ function from
    // __host__ __device__ function
    __host__
    ~Det_Int() {}

		static const eDevice device = e_host;

		size_type size() const
		{
			return image.size();
		}

		Vector<TVector, e_host> image;
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
	

} // namespace mt

#endif
