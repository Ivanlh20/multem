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

#ifndef TYPES_MT_H
	#define TYPES_MT_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include <type_traits>
	#include <vector>
	#include <algorithm>
	#include <thread>
	#include <mutex>
	#include <utility>
	#include <functional> 

	#include "const_enum_mt.cuh"
	#include "math.cuh"
	#include "type_traits_gen.cuh"
	#include "cgpu_fcns_gen.cuh"
	#include "r_2d.cuh"
	#include "r_3d.cuh"
	#include "particles.cuh"
	#include "types.cuh"

	namespace mt
	{
		/******************* spec layer information *********************/
		template <class T>
		class Spec_Lay_Info
		{
			public:
				R_3d<T> bs;				// box size
				R_3d<T> r_0;			// initial position
				dt_int32 tag;			// tag
				eSpec_Lay_Pos type;		// specimen layer type
				T sli_thick;			// slice thickness

				Spec_Lay_Info(): bs(), r_0(), tag(c_dflt_tag), type(eslp_none), sli_thick(2) {};

				Spec_Lay_Info(const R_3d<T>& bs, R_3d<T> r_0 = R_3d<T>(), dt_int32 tag = c_dflt_tag, 
				eSpec_Lay_Pos type = eslp_none, T sli_thick = 2.0): bs(bs), r_0(r_0), tag(tag), type(type), sli_thick(sli_thick) {};

				Spec_Lay_Info<T>& operator=(Spec_Lay_Info<T> &spec_lay_info)
				{
					if (this != &spec_lay_info)
					{
						bs = spec_lay_info.bs;
						r_0 = spec_lay_info.r_0;
						sli_thick = spec_lay_info.sli_thick;
						tag = spec_lay_info.tag;
						type = spec_lay_info.type;
					}

					return *this;
				}

				template <class U> 
				Spec_Lay_Info<T>& operator=(Spec_Lay_Info<U> &spec_lay_info)
				{
					assign(spec_lay_info);
			
					return *this;
				}

				template <class U> 
				void assign(const Spec_Lay_Info<U>& spec_lay_info)
				{ 
					if ((void*)this != (void*)&spec_lay_info)
					{
						bs = spec_lay_info.bs;
						r_0 = spec_lay_info.r_0;
						r_e = spec_lay_info.r_e;
						sli_thick = T(spec_lay_info.sli_thick);
						tag = spec_lay_info.tag;
						type = spec_lay_info.type;
					}
				}
				
				R_3d<T> r_e() const 
				{ 
					return r_0 + bs; 
				};	
				
				T z_e() const 
				{ 
					return r_0.z + bs.z; 
				};

				dt_bool is_spec_lay_top() const 
				{ 
					return is_spec_lay_top(type); 
				};

				dt_bool is_spec_lay_bottom() const 
				{ 
					return is_spec_lay_bottom(type); 
				};

				dt_bool is_spec_lay_middle() const 
				{ 
					return is_spec_lay_middle(type); 
				};

				dt_bool is_spec_lay_user_def() const 
				{ 
					return is_spec_lay_user_def(type); 
				};

				void set_region(Ptc_Atom<T>& atoms)
				{
					if (bs.z<1e-4)
					{
						tag = 0;
						return;
					}

					dt_int32 f_region = 0;
					dt_int32 c_region = 0;
					const T z_0 = r_0.z;
					const T z_e = z_e();
					for(auto iatoms = 0; iatoms < atoms.size(); iatoms++)
					{
						auto z = atoms.z[iatoms];
						if ((z_0<z) && (z<z_e))
						{
							f_region += atoms.tag[iatoms];
							c_region++;
						}
					}
					c_region = max(1, c_region);
					tag = static_cast<dt_int32>(::round(T(f_region)/T(c_region)));
				}															
		};
		
		template <class T>
		using Vctr_Spec_Lay_Info = Vctr_std<Spec_Lay_Info<T>>;


		/**************************** thickness ******************************/
		template <class T>
		struct Thick
		{
			Thick(): z(0), z_zero_def_plane(0), z_back_prop(0), 
			islice(0), iatom_e(0) {}

			T z;					// z
			T z_zero_def_plane;		// z: Zero defocus
			T z_back_prop;			// z: Back propagation

			dt_int32 islice;		// slice position
			dt_int32 iatom_e;		// Last atom index
		};

		/************************** stem fetector ****************************/
		template <class T, eDev Dev>
		struct Detector
		{
			using value_type = T;
			using size_type = dt_int32;

			static const eDev device = Dev;

			Detector(): type(mt::edt_circular) {}

			size_type size() const
			{
				size_type size_out = 0;
				switch (type)
				{
					case mt::edt_circular:
					{
						size_out = g_inner.size();
					}
					break;
					case mt::edt_radial:
					{
						size_out = fx.size();
					}
						break;
					case mt::edt_matrix:
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
				fcn.clear();
				grid_1d.clear();
				grid_2d.clear();
			}

			void resize(const size_type& new_size)
			{
				switch (type)
				{
					case mt::edt_circular:
					{
						g_inner.resize(new_size);
						g_outer.resize(new_size);
					}
					break;
					case mt::edt_radial:
					{
						fx.resize(new_size);
						fcn.resize(new_size);
						grid_1d.resize(new_size);
					}
						break;
					case mt::edt_matrix:
					{
						fR.resize(new_size);
						fcn.resize(new_size);
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
				for(auto i = 0; i<detector.fx.size(); i++)
				{
					fx[i].assign(detector.fx[i].begin(), detector.fx[i].end());
					fcn[i] = detector.fcn[i];
					grid_1d[i] = detector.grid_1d[i];
				}

				fR.resize(detector.fR.size());
				for(auto i = 0; i<detector.fR.size(); i++)
				{
					fR[i].assign(detector.fR[i].begin(), detector.fR[i].end());
					fcn[i] = detector.fcn[i];
					grid_2d[i] = detector.grid_2d[i];
				}
			}

			template <class TDetector> 
			Detector<T, Dev>& operator=(TDetector &detector)
			{
				assign(detector);
				return *this;
			}

			dt_bool is_detector_circular() const
			{
				return mt::is_detector_circular(type);
			}

			dt_bool is_detector_radial() const
			{
				return mt::is_detector_radial(type);
			}

			dt_bool is_detector_matrix() const
			{
				return mt::is_detector_matrix(type);
			}

			eDetector_Typ type;		// mt::edt_circular = 1, mt::edt_radial = 2, mt::edt_matrix = 3
			Vctr<T, edev_cpu> g_inner;		// Inner aperture Ang^-1
			Vctr<T, edev_cpu> g_outer;		// Outer aperture Ang^-1
			Vctr<Vctr<T, Dev>, edev_cpu> fx;		// radial sensitivity value
			Vctr<Vctr<T, Dev>, edev_cpu> fR;		// 2D sensitivity value
			std::vector<Grid_2d<T>> grid_1d;		// grid_1d
			std::vector<Grid_2d<T>> grid_2d;		// grid_2d
			std::vector<std::string> fcn;		// file names
		};

		/************************* stem intensity ****************************/
		template <class TVctr>
		struct Det_Int
		{
			using value_type = typename TVctr::value_type;
			using size_type = dt_int32;

			static const eDev device = edev_cpu;

			size_type size() const
			{
				return image.size();
			}

			Vctr<TVctr, edev_cpu> image;
		};

		/******************** radial schrodinger equation ********************/
		template <class T>
		class In_Rad_Schr
		{
		public:
			T E_0;		// Acceleration Voltage
			eAtomic_Pot_Parm_Typ atomic_pot_parm_typ;		// Parameterization type
			dt_int32 n;		// Principal quantum number
			dt_int32 nr;		// number of grid points
			dt_int32 natomsM;		// number of atoms
			T* atomsM;		// atoms
		};
	}

#endif