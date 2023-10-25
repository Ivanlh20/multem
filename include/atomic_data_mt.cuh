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

#ifndef ATOMIC_DATA_MT_H
	#define ATOMIC_DATA_MT_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include <string>

	#include "const_enum_mt.cuh"
	#include "fcns_cgpu_gen.h"
	#include "intrpl_coef.cuh"

	namespace mt
	{
		/************************* atomic coefficients ***********************/
		template <class T, eDev Dev>
		class Atomic_Coef
		{
		public:
			using value_type = T;
			using size_type = dt_int32;
			static const eDev device = Dev;

			eAtomic_Pot_Parm_Typ atomic_pot_parm_typ;			// feg parameterization
			dt_int32 Z;							// atomic number
			dt_int32 charge;					// charge

			mutable LNL_Coef<T, Dev> feg;		// electron scattering factor coefficients
			mutable LNL_Coef<T, Dev> fxg;		// x-ray scattering factor coefficients
			mutable LNL_Coef<T, Dev> pr;		// electron density coefficients
			mutable LNL_Coef<T, Dev> vr;		// potential coefficients
			mutable LNL_Coef<T, Dev> vzp;		// projected potential coefficients

			/************************************* constructors ************************************/
			Atomic_Coef():atomic_pot_parm_typ(eappt_lobato_0_12), Z(0), charge(0) {}

			Atomic_Coef(const size_type& new_size): Atomic_Coef()
			{
				resize(new_size);
			}

			Atomic_Coef(const size_type& new_size, const T& value): Atomic_Coef()
			{
				resize(new_size, value);
			}

			/* copy constructor */
			Atomic_Coef(const Atomic_Coef<T, Dev>& atomic_coef): Atomic_Coef()
			{
				*this = atomic_coef;
			}

			/* converting constructor */
			template <class U, eDev Dev_u> 
			Atomic_Coef(const Atomic_Coef<U, Dev_u>& atomic_coef)
			{
				*this = atomic_coef;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			Atomic_Coef<T, Dev>& operator=(const Atomic_Coef<T, Dev>& atomic_coef)
			{
				this->assign(atomic_coef);
			
				return *this;
			}

			/* converting assignment operator */
			template <class U, eDev Dev_u> 
			Atomic_Coef<T, Dev>& operator=(const Atomic_Coef<U, Dev_u>& atomic_coef)
			{
				this->assign(atomic_coef);
			
				return *this;
			}

			template <class U, eDev Dev_u> 
			void assign(const Atomic_Coef<U, Dev_u>& atomic_coef)
			{ 
				atomic_pot_parm_typ = atomic_coef.atomic_pot_parm_typ;
				Z = atomic_coef.Z;
				charge = atomic_coef.charge;

				feg.assign(atomic_coef.feg);
				fxg.assign(atomic_coef.fxg);
				pr.assign(atomic_coef.pr);
				vr.assign(atomic_coef.vr);
				vzp.assign(atomic_coef.vzp);
			}

			void fill(const T& val)
			{
				feg.fill(val);
				fxg.fill(val);
				pr.fill(val);
				vr.fill(val);
				vzp.fill(val);
			}

			size_type size() const
			{
				return size_type(feg.size());
			}

			template <class ST>
			ST size() const
			{
				return static_cast<ST>(feg.size());
			}

			void clear()
			{
				atomic_pot_parm_typ = eappt_lobato_0_12;
				Z = 0;
				charge = 0;

				feg.clear();
				fxg.clear();
				pr.clear();
				vr.clear();
				vzp.clear();
			}

			void resize(const size_type& new_size)
			{
				feg.resize(new_size);
				fxg.resize(new_size);
				pr.resize(new_size);
				vr.resize(new_size);
				vzp.resize(new_size);
			}

			void resize(const size_type& new_size, const T& value)
			{
				feg.resize(new_size, value);
				fxg.resize(new_size, value);
				pr.resize(new_size, value);
				vr.resize(new_size, value);
				vzp.resize(new_size, value);
			}

			void shrink_to_fit()
			{
				feg.shrink_to_fit();
				fxg.shrink_to_fit();
				pr.shrink_to_fit();
				vr.shrink_to_fit();
				vzp.shrink_to_fit();
			}

			T Z_diff()
			{
				return T(Z - charge);
			}
		};

		template <class T>
		using Atomic_Coef_cpu = Atomic_Coef<T, edev_cpu>;

		template <class T>
		using Atomic_Coef_gpu = Atomic_Coef<T, edev_gpu>;

		/*************************** atomic info *****************************/
		template <class T, eDev Dev>
		class Atomic_Info
		{
		public:
			using value_type = T;
			using size_type = dt_int32;
 
			std::string name;									// atom name
			dt_int32 Z;											// atomic number
			dt_int32 A;											// mass number
			T m;												// atomic mass
			T rn;												// experimental nuclear radius (Angs.)
			T ra;												// experimental atomic radius (Angs.)

			mutable Vctr_cpu<dt_float32> eels_maj_edg;			// major eels edges
			mutable Vctr_cpu<dt_float32> eels_min_edg;			// minor eels edges
			mutable Vctr_cpu<Atomic_Coef<T, Dev>> coef;			// atomic coefficients

			/************************************* constructors ************************************/
			Atomic_Info(): name(""), Z(0), A(0), m(0), rn(0), ra(0){}			
			
			Atomic_Info(const std::string& name, const dt_int32& Z, const dt_int32& A, const T& m, 
			const T& rn, const T& ra): name(name), Z(Z), A(A), m(m), rn(rn), ra(ra) {}

			Atomic_Info(const std::string& name, const dt_int32& Z, const dt_int32& A, const T& m, 
			const T& rn, const T& ra, const dt_init_list_f64& eels_maj_edg, const dt_init_list_f64& eels_min_edg): Atomic_Info(name, Z, A, m, rn, ra)
			{
				this->eels_maj_edg.assign(eels_maj_edg.begin(), eels_maj_edg.end());
				this->eels_min_edg.assign(eels_min_edg.begin(), eels_min_edg.end());
			}

			/* copy constructor */
			Atomic_Info(const Atomic_Info<T, Dev>& atomic_info)
			{
				*this = atomic_info;
			}

			/* converting constructor */
			template <class U, eDev Dev_u> 
			Atomic_Info(const Atomic_Info<U, Dev_u>& atomic_info)
			{
				*this = atomic_info;
			}

			/******************************** assignment operators *********************************/
			/* copy assignment operator */
			Atomic_Info<T, Dev>& operator=(const Atomic_Info<T, Dev>& atomic_info)
			{
				if (this != &atomic_info)
				{
					this->assign(atomic_info);
				}
			
				return *this;
			}

			/* converting assignment operator */
			template <class U, eDev Dev_u> 
			Atomic_Info<T, Dev>& operator=(const Atomic_Info<U, Dev_u>& atomic_info)
			{
				this->assign(atomic_info);
			
				return *this;
			}

			template <class U, eDev Dev_u> 
			void assign(const Atomic_Info<U, Dev_u>& atomic_info)
			{ 
				name = atomic_info.name;
				Z = atomic_info.Z;
				A = atomic_info.A;
				m = T(atomic_info.m);
				rn = T(atomic_info.rn);
				ra = T(atomic_info.ra);

				eels_maj_edg.assign(atomic_info.eels_maj_edg);
				eels_min_edg.assign(atomic_info.eels_min_edg);

				coef.resize(atomic_info.coef.size());
				for(auto ik = 0; ik < atomic_info.coef.size(); ik++)
				{
					coef[ik].assign(atomic_info.coef[ik]);
				}

			}

			template <class U>
			void add_atomic_coef(const Atomic_Coef_cpu<U>& atomic_coef)
			{
				if (atomic_coef.size()==0)
				{
					return;
				}
				coef.resize(coef.size()+1);
				coef.back() = atomic_coef;
			}

			T Z_r()
			{
				return T(Z);
			}

			int mass_number() const 
			{ 
				return A;
			}

			T atomic_mass() const
			{
				return m;
			}

			T nuclear_radius() const
			{
				return rn;
			}

			T atomic_radius() const
			{
				return ra;
			}

			T nuclear_radius_cal() const 
			{ 
				return 1.2e-05*pow(T(A), T(1.0/3.0));
			}
		};

		template <class T>
		using Atomic_Info_cpu = Atomic_Info<T, edev_cpu>;

		template <class T>
		using Atomic_Info_gpu = Atomic_Info<T, edev_gpu>;

		/**************************** atomic data ****************************/
		class Atomic_Data{
		public:
			using T = double;
			using value_type = T;
			using size_type = dt_int32;
			static const eDev device = edev_cpu;

			Atomic_Data() {}

			Atomic_Data(const dt_int32& Z) 
			{
				atomic_info = operator()(Z);
			}

			Atomic_Data(const dt_int32& Z, const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, dt_int32 charge=0) 
			{
				atomic_info = operator()(Z, atomic_pot_parm_typ, charge);
			}

			Atomic_Data(const dt_int32& Z, const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const Vctr_int32_cpu& charge)
			{
				atomic_info = operator()(Z, atomic_pot_parm_typ, charge);
			}

			// conversion operator
			operator Atomic_Info_cpu<dt_float32>() const
			{
				return atomic_info;
			}

			// conversion operator
			operator Atomic_Info_cpu<dt_float64>() const
			{
				return atomic_info;
			}

			// get atomic info
			Atomic_Info_cpu<T> operator()(const dt_int32& Z)		
			{
				return operator()(Z, eappt_lobato_0_12);
			}

			// get atomic info
			Atomic_Info_cpu<T> operator()(const dt_int32& Z, const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, dt_int32 charge=0)		
			{
				auto atomic_info = load_atomic_info_1(Z);

				switch(atomic_pot_parm_typ)
				{
					case eappt_doyle_0_4:
					{
						atomic_info.add_atomic_coef(atomic_coef_doyle_neutral_0_4(Z));
					}
					break;
					case eappt_peng_0_4:
					{
						atomic_info.add_atomic_coef(atomic_coef_peng_neutral_0_4(Z));
					}
					break;
					case eappt_peng_0_12:
					{
						atomic_info.add_atomic_coef(atomic_coef_peng_neutral_0_12(Z));
					}
					break;
					case eappt_kirkland_0_12:
					{
						atomic_info.add_atomic_coef(atomic_coef_kirkland_neutral_0_12(Z));
					}
					break;
					case eappt_weickenmeier_0_12:
					{
						atomic_info.add_atomic_coef(atomic_coef_weickenmeier_neutral_0_12(Z));
					}
					break;
					case eappt_lobato_0_12:
					{
						atomic_info.add_atomic_coef(atomic_coef_lobato_neutral_0_12(Z));
					}
					break;
					case eappt_peng_ion_0_4:
					{
						auto atomic_coef = atomic_coef_peng_ion_0_4(Z, charge);

						if (atomic_coef.size()==0)
						{
							atomic_coef = atomic_coef_peng_neutral_0_4(Z);
						}

						atomic_info.add_atomic_coef(atomic_coef);
					}
					break;
				}

				return atomic_info;
			}

			// get atomic info
			Atomic_Info_cpu<T> operator()(const dt_int32& Z, const eAtomic_Pot_Parm_Typ& atomic_pot_parm_typ, const Vctr_int32_cpu& charge)		
			{
				auto atomic_info = load_atomic_info_1(Z);

				for(auto ik=0; ik < charge.size(); ik++)
				{
					switch(atomic_pot_parm_typ)
					{
						case eappt_doyle_0_4:
						{
							atomic_info.add_atomic_coef(atomic_coef_doyle_neutral_0_4(Z));
						}
						break;
						case eappt_peng_0_4:
						{
							atomic_info.add_atomic_coef(atomic_coef_peng_neutral_0_4(Z));
						}
						break;
						case eappt_peng_0_12:
						{
							atomic_info.add_atomic_coef(atomic_coef_peng_neutral_0_12(Z));
						}
						break;
						case eappt_kirkland_0_12:
						{
							atomic_info.add_atomic_coef(atomic_coef_kirkland_neutral_0_12(Z));
						}
						break;
						case eappt_weickenmeier_0_12:
						{
							atomic_info.add_atomic_coef(atomic_coef_weickenmeier_neutral_0_12(Z));
						}
						break;
						case eappt_lobato_0_12:
						{
							atomic_info.add_atomic_coef(atomic_coef_lobato_neutral_0_12(Z));
						}
						break;
						case eappt_peng_ion_0_4:
						{
							auto atomic_coef = atomic_coef_peng_ion_0_4(Z, charge[ik]);

							if (atomic_coef.size()==0)
							{
								atomic_coef = atomic_coef_peng_neutral_0_4(Z);
							}

							atomic_info.add_atomic_coef(atomic_coef);
						}
						break;
					}
				}
				
				return atomic_info;
			}

		private:
			Atomic_Info_cpu<double> atomic_info;

			/***************************************************************************************/
			// load atomic info
			Atomic_Info_cpu<double> load_atomic_info_1(const dt_int32& Z)
			{
				switch(Z)
				{
					case 1:
						return {"H", 1, 1, 1.0080e+00, 8.7770e-06, 7.9000e-01, {13.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 2:
						return {"He", 2, 4, 4.0026e+00, 1.6753e-05, 4.9000e-01, {22.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 3:
						return {"Li", 3, 7, 6.9410e+00, 2.4173e-05, 2.0500e+00, {55.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 4:
						return {"Be", 4, 9, 9.0122e+00, 2.5180e-05, 1.4000e+00, {111.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 5:
						return {"B", 5, 11, 1.0811e+01, 2.4060e-05, 1.1700e+00, {188.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 6:
						return {"C", 6, 12, 1.2011e+01, 2.4702e-05, 9.1000e-01, {284.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 7:
						return {"N", 7, 14, 1.4007e+01, 2.5582e-05, 7.5000e-01, {401.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 8:
						return {"O", 8, 16, 1.5999e+01, 2.6991e-05, 6.5000e-01, {532.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 9:
						return {"F", 9, 19, 1.8998e+01, 2.8976e-05, 5.7000e-01, {685.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 10:
						return {"Ne", 10, 20, 2.0180e+01, 3.0058e-05, 5.1000e-01, {867.00, 18.00, 0.00, 0.00, 0.00}, {45.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 11:
						return {"Na", 11, 23, 2.2990e+01, 2.9935e-05, 2.2300e+00, {1072.00, 31.00, 0.00, 0.00, 0.00}, {63.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 12:
						return {"Mg", 12, 24, 2.4305e+01, 3.0570e-05, 1.7200e+00, {1305.00, 51.00, 0.00, 0.00, 0.00}, {89.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 13:
						return {"Al", 13, 27, 2.6982e+01, 3.0610e-05, 1.6200e+00, {1560.00, 73.00, 0.00, 0.00, 0.00}, {118.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 14:
						return {"Si", 14, 28, 2.8085e+01, 3.1224e-05, 1.4400e+00, {1839.00, 99.00, 0.00, 0.00, 0.00}, {149.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 15:
						return {"P", 15, 31, 3.0974e+01, 3.1889e-05, 1.2300e+00, {2146.00, 132.00, 0.00, 0.00, 0.00}, {189.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 16:
						return {"S", 16, 32, 3.2066e+01, 3.2847e-05, 1.0900e+00, {2472.00, 165.00, 0.00, 0.00, 0.00}, {229.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 17:
						return {"Cl", 17, 35, 3.5453e+01, 3.3654e-05, 9.7000e-01, {2822.00, 200.00, 0.00, 0.00, 0.00}, {270.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 18:
						return {"Ar", 18, 40, 3.9948e+01, 3.4276e-05, 8.8000e-01, {3203.00, 245.00, 12.00, 0.00, 0.00}, {320.00, 25.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 19:
						return {"K", 19, 39, 3.9098e+01, 3.4350e-05, 2.7700e+00, {3607.00, 296.00, 294.00, 18.00, 0.00}, {377.00, 34.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 20:
						return {"Ca", 20, 40, 4.0078e+01, 3.4777e-05, 2.2300e+00, {4038.00, 350.00, 346.00, 25.00, 0.00}, {438.00, 44.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 21:
						return {"Sc", 21, 45, 4.4956e+01, 3.5460e-05, 2.0900e+00, {4493.00, 407.00, 402.00, 32.00, 0.00}, {500.00, 54.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 22:
						return {"Ti", 22, 48, 4.7880e+01, 3.5922e-05, 2.0000e+00, {4966.00, 462.00, 456.00, 35.00, 0.00}, {564.00, 60.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 23:
						return {"V", 23, 51, 5.0941e+01, 3.6002e-05, 1.9200e+00, {521.00, 513.00, 38.00, 0.00, 0.00}, {628.00, 66.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 24:
						return {"Cr", 24, 52, 5.1996e+01, 3.6452e-05, 1.8500e+00, {584.00, 575.00, 42.00, 0.00, 0.00}, {695.00, 74.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 25:
						return {"Mn", 25, 55, 5.4938e+01, 3.7057e-05, 1.7900e+00, {651.00, 640.00, 49.00, 0.00, 0.00}, {769.00, 84.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 26:
						return {"Fe", 26, 56, 5.5847e+01, 3.7384e-05, 1.7200e+00, {721.00, 708.00, 54.00, 0.00, 0.00}, {846.00, 93.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 27:
						return {"Co", 27, 59, 5.8933e+01, 3.7875e-05, 1.6700e+00, {794.00, 779.00, 60.00, 0.00, 0.00}, {926.00, 101.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 28:
						return {"Ni", 28, 58, 5.8693e+01, 3.7757e-05, 1.6200e+00, {872.00, 855.00, 68.00, 0.00, 0.00}, {1008.00, 112.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 29:
						return {"Cu", 29, 63, 6.3456e+01, 3.8823e-05, 1.5700e+00, {951.00, 931.00, 74.00, 0.00, 0.00}, {1096.00, 120.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 30:
						return {"Zn", 30, 64, 6.5390e+01, 3.9279e-05, 1.5300e+00, {1043.00, 1020.00, 0.00, 0.00, 0.00}, {1194.00, 136.00, 87.00, 0.00, 0.00, 0.00, 0.00}};
					case 31:
						return {"Ga", 31, 69, 6.9723e+01, 3.9973e-05, 1.8100e+00, {1142.00, 1115.00, 0.00, 0.00, 0.00}, {1298.00, 158.00, 103.00, 0.00, 0.00, 0.00, 0.00}};
					case 32:
						return {"Ge", 32, 74, 7.2610e+01, 4.0742e-05, 1.5200e+00, {1248.00, 1217.00, 29.00, 0.00, 0.00}, {1414.00, 180.00, 121.00, 0.00, 0.00, 0.00, 0.00}};
					case 33:
						return {"As", 33, 75, 7.4922e+01, 4.0968e-05, 1.3300e+00, {1359.00, 1323.00, 41.00, 0.00, 0.00}, {1526.00, 203.00, 140.00, 0.00, 0.00, 0.00, 0.00}};
					case 34:
						return {"Se", 34, 80, 7.8960e+01, 4.1400e-05, 1.2200e+00, {1476.00, 1436.00, 57.00, 0.00, 0.00}, {1654.00, 231.00, 162.00, 0.00, 0.00, 0.00, 0.00}};
					case 35:
						return {"Br", 35, 79, 7.9904e+01, 4.1629e-05, 1.1200e+00, {1596.00, 1150.00, 69.00, 0.00, 0.00}, {1782.00, 256.00, 181.00, 0.00, 0.00, 0.00, 0.00}};
					case 36:
						return {"Kr", 36, 84, 8.3800e+01, 4.1883e-05, 1.0300e+00, {1727.00, 1675.00, 89.00, 11.00, 0.00}, {1921.00, 289.00, 214.00, 24.00, 0.00, 0.00, 0.00}};
					case 37:
						return {"Rb", 37, 85, 8.5468e+01, 4.2037e-05, 2.9800e+00, {1864.00, 1804.00, 110.00, 14.00, 0.00}, {2065.00, 322.00, 247.00, 238.00, 29.00, 0.00, 0.00}};
					case 38:
						return {"Sr", 38, 88, 8.7620e+01, 4.2246e-05, 2.4500e+00, {2007.00, 1940.00, 133.00, 20.00, 0.00}, {2216.00, 357.00, 280.00, 269.00, 38.00, 0.00, 0.00}};
					case 39:
						return {"Y", 39, 89, 8.8906e+01, 4.2438e-05, 2.2700e+00, {2155.00, 2080.00, 157.00, 26.00, 0.00}, {2372.00, 394.00, 312.00, 300.00, 45.00, 0.00, 0.00}};
					case 40:
						return {"Zr", 40, 90, 9.1224e+01, 4.2690e-05, 2.1600e+00, {2307.00, 2222.00, 180.00, 29.00, 0.00}, {2532.00, 430.00, 344.00, 330.00, 51.00, 0.00, 0.00}};
					case 41:
						return {"Nb", 41, 93, 9.2906e+01, 4.3240e-05, 2.0800e+00, {2465.00, 2371.00, 205.00, 34.00, 0.00}, {2698.00, 468.00, 378.00, 363.00, 58.00, 0.00, 0.00}};
					case 42:
						return {"Mo", 42, 98, 9.5940e+01, 4.4090e-05, 2.0100e+00, {2625.00, 2520.00, 227.00, 35.00, 0.00}, {2865.00, 505.00, 410.00, 392.00, 62.00, 0.00, 0.00}};
					case 43:
						return {"Tc", 43, 98, 9.8000e+01, 4.4130e-05, 1.9500e+00, {2793.00, 2677.00, 253.00, 39.00, 0.00}, {3042.00, 544.00, 445.00, 425.00, 68.00, 0.00, 0.00}};
					case 44:
						return {"Ru", 44, 102, 1.0107e+02, 4.4811e-05, 1.8900e+00, {2967.00, 2838.00, 279.00, 43.00, 0.00}, {3224.00, 585.00, 483.00, 461.00, 75.00, 0.00, 0.00}};
					case 45:
						return {"Rh", 45, 103, 1.0291e+02, 4.4945e-05, 1.8300e+00, {3146.00, 3004.00, 307.00, 48.00, 0.00}, {3412.00, 627.00, 521.00, 496.00, 81.00, 0.00, 0.00}};
					case 46:
						return {"Pd", 46, 106, 1.0642e+02, 4.5324e-05, 1.7900e+00, {3330.00, 3173.00, 335.00, 51.00, 0.00}, {3604.00, 670.00, 559.00, 531.00, 86.00, 0.00, 0.00}};
					case 47:
						return {"Ag", 47, 107, 1.0787e+02, 4.5464e-05, 1.7500e+00, {3524.00, 3351.00, 367.00, 0.00, 0.00}, {3806.00, 717.00, 602.00, 571.00, 95.00, 56.00, 0.00}};
					case 48:
						return {"Cd", 48, 114, 1.1241e+02, 4.6068e-05, 1.7100e+00, {3727.00, 3538.00, 404.00, 0.00, 0.00}, {4018.00, 770.00, 651.00, 616.00, 108.00, 67.00, 0.00}};
					case 49:
						return {"In", 49, 115, 1.1482e+02, 4.6155e-05, 2.0000e+00, {3938.00, 3730.00, 443.00, 0.00, 0.00}, {4237.00, 826.00, 702.00, 664.00, 122.00, 77.00, 0.00}};
					case 50:
						return {"Sn", 50, 120, 1.1871e+02, 4.6525e-05, 1.7200e+00, {4156.00, 3929.00, 485.00, 24.00, 0.00}, {4465.00, 884.00, 756.00, 714.00, 136.00, 89.00, 0.00}};
					case 51:
						return {"Sb", 51, 121, 1.2176e+02, 4.6802e-05, 1.5300e+00, {4380.00, 4132.00, 528.00, 31.00, 0.00}, {944.00, 812.00, 766.00, 152.00, 98.00, 0.00, 0.00}};
					case 52:
						return {"Te", 52, 130, 1.2760e+02, 4.7420e-05, 1.4200e+00, {4612.00, 4341.00, 40.00, 572.00, 0.00}, {1006.00, 870.00, 819.00, 168.00, 110.00, 0.00, 0.00}};
					case 53:
						return {"I", 53, 127, 1.2690e+02, 4.7500e-05, 1.3200e+00, {4852.00, 4557.00, 619.00, 50.00, 0.00}, {1072.00, 930.00, 875.00, 186.00, 123.00, 0.00, 0.00}};
					case 54:
						return {"Xe", 54, 132, 1.3129e+02, 4.7850e-05, 1.2400e+00, {4782.00, 672.00, 63.00, 0.00, 0.00}, {1145.00, 999.00, 937.00, 208.00, 147.00, 0.00, 0.00}};
					case 55:
						return {"Cs", 55, 133, 1.3291e+02, 4.8040e-05, 3.3400e+00, {740.00, 726.00, 76.00, 0.00, 0.00}, {1217.00, 1065.00, 998.00, 231.00, 162.00, 0.00, 0.00}};
					case 56:
						return {"Ba", 56, 138, 1.3733e+02, 4.8370e-05, 2.7600e+00, {796.00, 781.00, 90.00, 15.00, 0.00}, {1293.00, 1137.00, 1062.00, 253.00, 180.00, 0.00, 0.00}};
					case 57:
						return {"La", 57, 139, 1.3891e+02, 4.8549e-05, 2.7400e+00, {849.00, 832.00, 99.00, 14.00, 0.00}, {1361.00, 1204.00, 1123.00, 270.00, 191.00, 0.00, 0.00}};
					case 58:
						return {"Ce", 58, 140, 1.4012e+02, 4.8773e-05, 2.7000e+00, {901.00, 883.00, 110.00, 20.00, 0.00}, {1435.00, 1273.00, 1185.00, 290.00, 207.00, 0.00, 0.00}};
					case 59:
						return {"Pr", 59, 141, 1.4091e+02, 4.8919e-05, 2.6700e+00, {951.00, 931.00, 113.00, 22.00, 0.00}, {1511.00, 1337.00, 1242.00, 305.00, 218.00, 0.00, 0.00}};
					case 60:
						return {"Nd", 60, 142, 1.4424e+02, 4.9115e-05, 2.6400e+00, {1000.00, 978.00, 118.00, 21.00, 0.00}, {1575.00, 1403.00, 1297.00, 315.00, 225.00, 0.00, 0.00}};
					case 61:
						return {"Pm", 61, 145, 1.4500e+02, 4.9530e-05, 2.6200e+00, {1052.00, 1027.00, 120.00, 24.00, 0.00}, {1649.00, 1471.00, 1357.00, 331.00, 236.00, 0.00, 0.00}};
					case 62:
						return {"Sm", 62, 152, 1.5036e+02, 5.0823e-05, 2.5900e+00, {1106.00, 1080.00, 129.00, 21.00, 0.00}, {1723.00, 1541.00, 1420.00, 346.00, 247.00, 0.00, 0.00}};
					case 63:
						return {"Eu", 63, 153, 1.5197e+02, 5.1093e-05, 2.5600e+00, {1161.00, 1131.00, 133.00, 22.00, 0.00}, {1800.00, 1614.00, 1481.00, 360.00, 257.00, 0.00, 0.00}};
					case 64:
						return {"Gd", 64, 158, 1.5725e+02, 5.1578e-05, 2.5400e+00, {1217.00, 1185.00, 140.00, 20.00, 0.00}, {1881.00, 1688.00, 1544.00, 376.00, 271.00, 0.00, 0.00}};
					case 65:
						return {"Tb", 65, 159, 1.5893e+02, 5.0600e-05, 2.5100e+00, {1275.00, 1241.00, 147.00, 25.00, 0.00}, {1967.00, 1768.00, 1611.00, 398.00, 285.00, 0.00, 0.00}};
					case 66:
						return {"Dy", 66, 164, 1.6250e+02, 5.2240e-05, 2.4900e+00, {1332.00, 1295.00, 154.00, 26.00, 0.00}, {2047.00, 1842.00, 1676.00, 416.00, 293.00, 0.00, 0.00}};
					case 67:
						return {"Ho", 67, 165, 1.6493e+02, 5.2022e-05, 2.4700e+00, {1391.00, 1351.00, 161.00, 20.00, 0.00}, {2128.00, 1923.00, 1741.00, 436.00, 307.00, 51.00, 0.00}};
					case 68:
						return {"Er", 68, 166, 1.6726e+02, 5.2527e-05, 2.4500e+00, {1453.00, 1409.00, 168.00, 29.00, 0.00}, {2206.00, 2006.00, 1812.00, 449.00, 320.00, 60.00, 0.00}};
					case 69:
						return {"Tm", 69, 169, 1.6893e+02, 5.2256e-05, 2.4200e+00, {1515.00, 1468.00, 180.00, 32.00, 0.00}, {2307.00, 2090.00, 1884.00, 472.00, 337.00, 53.00, 0.00}};
					case 70:
						return {"Yb", 70, 174, 1.7304e+02, 5.3105e-05, 2.4000e+00, {1576.00, 1528.00, 185.00, 24.00, 0.00}, {2398.00, 2173.00, 1950.00, 487.00, 343.00, 54.00, 0.00}};
					case 71:
						return {"Lu", 71, 175, 1.7497e+02, 5.3700e-05, 2.2500e+00, {1639.00, 1588.00, 195.00, 28.00, 0.00}, {2491.00, 2263.00, 2024.00, 506.00, 359.00, 57.00, 0.00}};
					case 72:
						return {"Hf", 72, 180, 1.7849e+02, 5.3482e-05, 2.1600e+00, {1716.00, 1662.00, 31.00, 0.00, 0.00}, {2601.00, 2365.00, 2108.00, 538.00, 380.00, 214.00, 65.00}};
					case 73:
						return {"Ta", 73, 181, 1.8095e+02, 5.3507e-05, 2.0900e+00, {1793.00, 1735.00, 36.00, 0.00, 0.00}, {2708.00, 2469.00, 2194.00, 565.00, 404.00, 229.00, 71.00}};
					case 74:
						return {"W", 74, 184, 1.8385e+02, 5.3646e-05, 2.0200e+00, {1872.00, 1809.00, 36.00, 0.00, 0.00}, {2820.00, 2575.00, 2281.00, 595.00, 425.00, 245.00, 77.00}};
					case 75:
						return {"Re", 75, 187, 1.8621e+02, 5.3697e-05, 1.9700e+00, {1949.00, 1883.00, 35.00, 0.00, 0.00}, {2932.00, 2682.00, 2367.00, 625.00, 444.00, 260.00, 83.00}};
					case 76:
						return {"Os", 76, 192, 1.9020e+02, 5.4122e-05, 1.9200e+00, {2031.00, 1960.00, 45.00, 0.00, 0.00}, {2792.00, 2457.00, 654.00, 468.00, 273.00, 46.00, 84.00}};
					case 77:
						return {"Ir", 77, 193, 1.9222e+02, 5.4020e-05, 1.8700e+00, {2116.00, 2040.00, 50.00, 0.00, 0.00}, {2909.00, 2551.00, 690.00, 494.00, 295.00, 60.00, 95.00}};
					case 78:
						return {"Pt", 78, 195, 1.9508e+02, 5.4250e-05, 1.8300e+00, {2202.00, 2122.00, 52.00, 0.00, 0.00}, {3026.00, 2645.00, 722.00, 519.00, 313.00, 102.00, 71.00}};
					case 79:
						return {"Au", 79, 197, 1.9697e+02, 5.4379e-05, 1.7900e+00, {2291.00, 2206.00, 54.00, 0.00, 0.00}, {3148.00, 2743.00, 759.00, 545.00, 334.00, 108.00, 83.00}};
					case 80:
						return {"Hg", 80, 202, 2.0059e+02, 5.4637e-05, 1.7600e+00, {2385.00, 2295.00, 58.00, 0.00, 0.00}, {3278.00, 2847.00, 800.00, 571.00, 360.00, 120.00, 98.00}};
					case 81:
						return {"Tl", 81, 205, 2.0438e+02, 5.4757e-05, 2.0800e+00, {2485.00, 2389.00, 75.00, 0.00, 0.00}, {3416.00, 2957.00, 845.00, 609.00, 386.00, 136.00, 118.00}};
					case 82:
						return {"Pb", 82, 208, 2.0720e+02, 5.5010e-05, 1.8100e+00, {2586.00, 2484.00, 86.00, 19.00, 0.00}, {3554.00, 3066.00, 894.00, 644.00, 413.00, 147.00, 138.00}};
					case 83:
						return {"Bi", 83, 209, 2.0898e+02, 5.5210e-05, 1.6300e+00, {2688.00, 2580.00, 93.00, 24.00, 0.00}, {3696.00, 3177.00, 938.00, 679.00, 440.00, 159.00, 157.00}};
					case 84:
						return {"Po", 84, 209, 2.0900e+02, 5.5262e-05, 1.5300e+00, {3491.00, 3332.00, 0.00, 0.00, 0.00}, {4046.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 85:
						return {"At", 85, 210, 2.1000e+02, 5.5310e-05, 1.4300e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 86:
						return {"Rn", 86, 222, 2.2200e+02, 5.6547e-05, 1.3400e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 87:
						return {"Fr", 87, 223, 2.2300e+02, 5.6584e-05, 2.7000e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 88:
						return {"Ra", 88, 226, 2.2603e+02, 5.6840e-05, 2.2300e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 89:
						return {"Ac", 89, 227, 2.2700e+02, 5.7120e-05, 1.8800e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 90:
						return {"Th", 90, 232, 2.3204e+02, 5.7800e-05, 1.8000e+00, {3491.00, 3332.00, 344.00, 335.00, 88.00}, {4046.00, 1329.00, 967.00, 714.00, 676.00, 290.00, 182.00}};
					case 91:
						return {"Pa", 91, 231, 2.3104e+02, 5.7660e-05, 1.6100e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 92:
						return {"U", 92, 238, 2.3803e+02, 5.8571e-05, 1.3800e+00, {3728.00, 3552.00, 391.00, 381.00, 96.00}, {4303.00, 1441.00, 1045.00, 780.00, 738.00, 324.00, 195.00}};
					case 93:
						return {"Np", 93, 237, 2.3705e+02, 5.8410e-05, 1.3000e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 94:
						return {"Pu", 94, 244, 2.4400e+02, 5.8948e-05, 1.5100e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 95:
						return {"Am", 95, 243, 2.4300e+02, 5.9042e-05, 1.8400e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 96:
						return {"Cm", 96, 247, 2.4700e+02, 5.9630e-05, 1.8400e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 97:
						return {"Bk", 97, 247, 2.4700e+02, 5.9300e-05, 1.8400e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 98:
						return {"Cf", 98, 251, 2.5100e+02, 6.0200e-05, 1.8400e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 99:
						return {"Es", 99, 252, 2.5200e+02, 6.0340e-05, 1.8400e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 100:
						return {"Fm", 100, 257, 2.5700e+02, 6.1070e-05, 1.8400e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 101:
						return {"Md", 101, 258, 2.5800e+02, 6.1220e-05, 1.8400e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 102:
						return {"No", 102, 259, 2.5900e+02, 6.1370e-05, 1.8400e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
					case 103:
						return {"Lr", 103, 262, 2.6000e+02, 6.1820e-05, 1.8400e+00, {0.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}};
				}
				return {};
			}

			// 1: doyle and Turner parameterization - 4 Gaussians - [0, 4]
			LNL_Coef_cpu<T> load_feg_doyle_neutral_0_4(const dt_int32& Z)
			{
				switch(Z)
				{
					case 2:
						return {{9.0600e-02, 1.8140e-01, 1.0950e-01, 3.6200e-02}, {1.8183e+01, 6.2109e+00, 1.8026e+00, 2.8440e-01}};
					case 3:
						return {{1.6108e+00, 1.2460e+00, 3.2570e-01, 9.8600e-02}, {1.0764e+02, 3.0480e+01, 4.5331e+00, 4.9510e-01}};
					case 4:
						return {{1.2498e+00, 1.3335e+00, 3.6030e-01, 1.0550e-01}, {6.0804e+01, 1.8591e+01, 3.6534e+00, 4.1570e-01}};
					case 5:
						return {{9.4460e-01, 1.3120e+00, 4.1880e-01, 1.1590e-01}, {4.6444e+01, 1.4178e+01, 3.2228e+00, 3.7670e-01}};
					case 6:
						return {{7.3070e-01, 1.1951e+00, 4.5630e-01, 1.2470e-01}, {3.6995e+01, 1.1297e+01, 2.8139e+00, 3.4560e-01}};
					case 7:
						return {{5.7170e-01, 1.0425e+00, 4.6470e-01, 1.3110e-01}, {2.8846e+01, 9.0542e+00, 2.4213e+00, 3.1670e-01}};
					case 8:
						return {{4.5400e-01, 9.1730e-01, 4.7190e-01, 1.3840e-01}, {2.3780e+01, 7.6220e+00, 2.1440e+00, 2.9590e-01}};
					case 9:
						return {{3.6860e-01, 8.1090e-01, 4.7510e-01, 1.4590e-01}, {2.0239e+01, 6.6093e+00, 1.9310e+00, 2.7930e-01}};
					case 10:
						return {{3.0250e-01, 7.2020e-01, 4.7510e-01, 1.5340e-01}, {1.7639e+01, 5.8604e+00, 1.7623e+00, 2.6560e-01}};
					case 11:
						return {{2.2406e+00, 1.3326e+00, 9.0700e-01, 2.8630e-01}, {1.0800e+02, 2.4505e+01, 3.3914e+00, 4.3460e-01}};
					case 12:
						return {{2.2692e+00, 1.8025e+00, 8.3940e-01, 2.8920e-01}, {7.3670e+01, 2.0175e+01, 3.0131e+00, 4.0460e-01}};
					case 13:
						return {{2.2756e+00, 2.4280e+00, 8.5780e-01, 3.2660e-01}, {7.2322e+01, 1.9773e+01, 3.0799e+00, 4.0760e-01}};
					case 14:
						return {{2.1293e+00, 2.5333e+00, 8.3490e-01, 3.2160e-01}, {5.7775e+01, 1.6476e+01, 2.6796e+00, 3.8600e-01}};
					case 15:
						return {{1.8882e+00, 2.4685e+00, 8.0460e-01, 3.2040e-01}, {4.4876e+01, 1.3538e+01, 2.6424e+00, 3.6030e-01}};
					case 16:
						return {{1.6591e+00, 2.3863e+00, 7.8990e-01, 3.2080e-01}, {3.6650e+01, 1.1488e+01, 2.4686e+00, 3.4030e-01}};
					case 17:
						return {{1.4524e+00, 2.2926e+00, 7.8740e-01, 3.2170e-01}, {3.0935e+01, 9.9798e+00, 2.3336e+00, 3.2280e-01}};
					case 18:
						return {{1.2736e+00, 2.1894e+00, 7.9270e-01, 3.2250e-01}, {2.6682e+01, 8.8130e+00, 2.2186e+00, 3.0710e-01}};
					case 19:
						return {{3.9507e+00, 2.5452e+00, 1.9795e+00, 4.0170e-01}, {1.3707e+02, 2.2402e+01, 4.5319e+00, 4.3400e-01}};
					case 20:
						return {{4.4696e+00, 2.9703e+00, 1.9696e+00, 4.8180e-01}, {9.9523e+01, 2.2696e+01, 4.1954e+00, 4.1650e-01}};
					case 21:
						return {{3.9659e+00, 2.9169e+00, 1.9254e+00, 4.8020e-01}, {3.8960e+01, 2.0606e+01, 3.8557e+00, 3.9880e-01}};
					case 22:
						return {{3.5653e+00, 2.8181e+00, 1.8930e+00, 4.8250e-01}, {8.1982e+01, 1.9049e+01, 3.5904e+00, 3.8550e-01}};
					case 23:
						return {{3.2449e+00, 2.6978e+00, 1.8597e+00, 4.8640e-01}, {7.6379e+01, 1.7726e+01, 3.3632e+00, 3.7430e-01}};
					case 24:
						return {{2.3066e+00, 2.3339e+00, 1.8226e+00, 4.9010e-01}, {7.8405e+01, 1.5785e+01, 3.1566e+00, 3.6360e-01}};
					case 25:
						return {{2.7467e+00, 2.4556e+00, 1.7923e+00, 4.9840e-01}, {6.7786e+01, 1.5674e+01, 2.9998e+00, 3.5690e-01}};
					case 26:
						return {{2.5440e+00, 2.3434e+00, 1.7588e+00, 5.0620e-01}, {6.4424e+01, 1.4531e+01, 2.8539e+00, 3.5020e-01}};
					case 27:
						return {{2.3668e+00, 2.2361e+00, 1.7243e+00, 5.1480e-01}, {6.1431e+01, 1.4179e+01, 2.7247e+00, 3.4420e-01}};
					case 28:
						return {{2.2104e+00, 2.1342e+00, 1.6891e+00, 5.2380e-01}, {5.8727e+01, 1.3553e+01, 2.6094e+00, 3.3880e-01}};
					case 29:
						return {{1.5792e+00, 1.8197e+00, 1.6576e+00, 5.3230e-01}, {6.2940e+01, 1.2453e+01, 2.5042e+00, 3.3310e-01}};
					case 30:
						return {{1.9418e+00, 1.9501e+00, 1.6192e+00, 5.4340e-01}, {5.4162e+01, 1.2518e+01, 2.4164e+00, 3.2950e-01}};
					case 31:
						return {{2.3205e+00, 2.4955e+00, 1.6879e+00, 5.9920e-01}, {6.5602e+01, 1.5458e+01, 2.5806e+00, 3.5100e-01}};
					case 32:
						return {{2.4467e+00, 2.7015e+00, 1.6157e+00, 6.0090e-01}, {5.5893e+01, 1.4393e+01, 2.4461e+00, 3.4150e-01}};
					case 33:
						return {{2.3989e+00, 2.7898e+00, 1.5288e+00, 5.9360e-01}, {4.5718e+01, 3.2817e+01, 2.2799e+00, 3.2770e-01}};
					case 34:
						return {{2.2980e+00, 2.8541e+00, 1.4555e+00, 5.8950e-01}, {3.8830e+01, 1.1536e+01, 2.1463e+00, 3.1630e-01}};
					case 35:
						return {{2.1659e+00, 2.9037e+00, 1.3951e+00, 5.8860e-01}, {3.3899e+01, 1.0500e+01, 2.0413e+00, 3.0700e-01}};
					case 36:
						return {{2.0338e+00, 2.9271e+00, 1.3425e+00, 5.8880e-01}, {2.9999e+01, 9.5977e+00, 1.9520e+00, 2.9860e-01}};
					case 37:
						return {{4.7760e+00, 3.8588e+00, 2.2339e+00, 8.6830e-01}, {1.4078e+02, 1.8991e+01, 3.7010e+00, 4.1940e-01}};
					case 38:
						return {{5.8478e+00, 4.0026e+00, 2.3420e+00, 8.7950e-01}, {1.0497e+02, 1.9367e+01, 3.7368e+00, 4.1420e-01}};
					case 42:
						return {{3.1199e+00, 3.9061e+00, 2.3615e+00, 8.5040e-01}, {7.2464e+01, 1.4642e+01, 3.2370e+00, 3.6620e-01}};
					case 47:
						return {{2.0355e+00, 3.2716e+00, 2.5105e+00, 8.3720e-01}, {6.1497e+01, 1.1324e+01, 2.8456e+00, 3.2710e-01}};
					case 48:
						return {{2.5737e+00, 3.2536e+00, 2.5468e+00, 8.3790e-01}, {5.5675e+01, 1.1838e+01, 2.7842e+00, 3.2170e-01}};
					case 49:
						return {{3.1528e+00, 3.5565e+00, 2.8180e+00, 8.8420e-01}, {6.6649e+01, 1.4449e+01, 2.9758e+00, 3.3450e-01}};
					case 50:
						return {{3.4495e+00, 3.7349e+00, 2.7779e+00, 8.7860e-01}, {5.9104e+01, 1.4179e+01, 2.8548e+00, 3.2700e-01}};
					case 51:
						return {{3.5644e+00, 3.8437e+00, 2.6366e+00, 8.6380e-01}, {5.0487e+01, 1.3316e+01, 2.6909e+00, 3.1610e-01}};
					case 53:
						return {{3.4728e+00, 4.0602e+00, 2.5215e+00, 8.3980e-01}, {3.9441e+01, 1.1816e+01, 2.4148e+00, 2.9760e-01}};
					case 54:
						return {{3.3656e+00, 4.1468e+00, 2.4430e+00, 8.2930e-01}, {3.5509e+01, 1.1117e+01, 2.2940e+00, 2.8920e-01}};
					case 55:
						return {{6.0620e+00, 5.9861e+00, 3.3033e+00, 1.0958e+00}, {1.5583e+02, 1.9695e+01, 3.3354e+00, 3.7930e-01}};
					case 56:
						return {{7.8212e+00, 6.0040e+00, 3.2803e+00, 1.1030e+00}, {1.1766e+02, 1.8778e+01, 3.2634e+00, 3.7600e-01}};
					case 57:
						return {{6.2661e+00, 4.8440e+00, 3.2023e+00, 1.2009e+00}, {1.0030e+02, 1.6066e+01, 2.9803e+00, 3.6740e-01}};
					case 79:
						return {{2.3880e+00, 4.2259e+00, 2.6886e+00, 1.2551e+00}, {4.2866e+01, 9.7430e+00, 2.2641e+00, 3.0670e-01}};
					case 80:
						return {{2.6817e+00, 4.2414e+00, 2.7549e+00, 1.2706e+00}, {4.2822e+01, 9.8557e+00, 2.2951e+00, 3.0670e-01}};
					case 82:
						return {{3.5099e+00, 4.5523e+00, 3.1539e+00, 1.3591e+00}, {5.2914e+01, 1.1884e+01, 2.5713e+00, 3.2050e-01}};
					case 83:
						return {{3.8412e+00, 4.6784e+00, 3.1924e+00, 1.3625e+00}, {5.0261e+01, 1.1999e+01, 2.5598e+00, 3.1770e-01}};
					case 86:
						return {{4.0779e+00, 4.9778e+00, 3.0955e+00, 1.3259e+00}, {2.8406e+01, 1.1020e+01, 2.3549e+00, 2.9910e-01}};
					case 92:
						return {{6.7668e+00, 6.7287e+00, 4.0135e+00, 1.5607e+00}, {8.5951e+01, 1.5642e+01, 2.9364e+00, 3.3480e-01}};
				}

				return {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
			}

			// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
			LNL_Coef_cpu<T> load_feg_peng_neutral_0_4(const dt_int32& Z)
			{
				switch(Z)
				{
					case 1:
						return {{3.490000e-02, 1.201000e-01, 1.970000e-01, 5.730000e-02, 1.195000e-01}, {5.347000e-01, 3.586700e+00, 1.234710e+01, 1.895250e+01, 3.862690e+01}};
					case 2:
						return {{3.170000e-02, 8.380000e-02, 1.526000e-01, 1.334000e-01, 1.640000e-02}, {2.507000e-01, 1.475100e+00, 4.493800e+00, 1.266460e+01, 3.116530e+01}};
					case 3:
						return {{7.500000e-02, 2.249000e-01, 5.548000e-01, 1.495400e+00, 9.354000e-01}, {3.864000e-01, 2.938300e+00, 1.538290e+01, 5.355450e+01, 1.387337e+02}};
					case 4:
						return {{7.800000e-02, 2.210000e-01, 6.740000e-01, 1.386700e+00, 6.925000e-01}, {3.131000e-01, 2.238100e+00, 1.015170e+01, 3.090610e+01, 7.832730e+01}};
					case 5:
						return {{9.090000e-02, 2.551000e-01, 7.738000e-01, 1.213600e+00, 4.606000e-01}, {2.995000e-01, 2.115500e+00, 8.381600e+00, 2.412920e+01, 6.313140e+01}};
					case 6:
						return {{8.930000e-02, 2.563000e-01, 7.570000e-01, 1.048700e+00, 3.575000e-01}, {2.465000e-01, 1.710000e+00, 6.409400e+00, 1.861130e+01, 5.025230e+01}};
					case 7:
						return {{1.022000e-01, 3.219000e-01, 7.982000e-01, 8.197000e-01, 1.715000e-01}, {2.451000e-01, 1.748100e+00, 6.192500e+00, 1.738940e+01, 4.814310e+01}};
					case 8:
						return {{9.740000e-02, 2.921000e-01, 6.910000e-01, 6.990000e-01, 2.039000e-01}, {2.067000e-01, 1.381500e+00, 4.694300e+00, 1.271050e+01, 3.247260e+01}};
					case 9:
						return {{1.083000e-01, 3.175000e-01, 6.487000e-01, 5.846000e-01, 1.421000e-01}, {2.057000e-01, 1.343900e+00, 4.278800e+00, 1.139320e+01, 2.878810e+01}};
					case 10:
						return {{1.269000e-01, 3.535000e-01, 5.582000e-01, 4.674000e-01, 1.460000e-01}, {2.200000e-01, 1.377900e+00, 4.020300e+00, 9.493400e+00, 2.312780e+01}};
					case 11:
						return {{2.142000e-01, 6.853000e-01, 7.692000e-01, 1.658900e+00, 1.448200e+00}, {3.334000e-01, 2.344600e+00, 1.008300e+01, 4.830370e+01, 1.382700e+02}};
					case 12:
						return {{2.314000e-01, 6.866000e-01, 9.677000e-01, 2.188200e+00, 1.133900e+00}, {3.278000e-01, 2.272000e+00, 1.092410e+01, 3.928980e+01, 1.019748e+02}};
					case 13:
						return {{2.390000e-01, 6.573000e-01, 1.201100e+00, 2.558600e+00, 1.231200e+00}, {3.138000e-01, 2.106300e+00, 1.041630e+01, 3.445520e+01, 9.853440e+01}};
					case 14:
						return {{2.519000e-01, 6.372000e-01, 1.379500e+00, 2.508200e+00, 1.050000e+00}, {3.075000e-01, 2.017400e+00, 9.674600e+00, 2.937440e+01, 8.047320e+01}};
					case 15:
						return {{2.548000e-01, 6.106000e-01, 1.454100e+00, 2.320400e+00, 8.477000e-01}, {2.908000e-01, 1.874000e+00, 8.517600e+00, 2.434340e+01, 6.329960e+01}};
					case 16:
						return {{2.497000e-01, 5.628000e-01, 1.389900e+00, 2.186500e+00, 7.715000e-01}, {2.681000e-01, 1.671100e+00, 7.026700e+00, 1.953770e+01, 5.038880e+01}};
					case 17:
						return {{2.443000e-01, 5.397000e-01, 1.391900e+00, 2.019700e+00, 6.621000e-01}, {2.468000e-01, 1.524200e+00, 6.153700e+00, 1.666870e+01, 4.230860e+01}};
					case 18:
						return {{2.385000e-01, 5.017000e-01, 1.342800e+00, 1.889900e+00, 6.079000e-01}, {2.289000e-01, 1.369400e+00, 5.256100e+00, 1.409280e+01, 3.553610e+01}};
					case 19:
						return {{4.115000e-01, 1.403100e+00, 2.278400e+00, 2.674200e+00, 2.216200e+00}, {3.703000e-01, 3.387400e+00, 1.310290e+01, 6.895920e+01, 1.944329e+02}};
					case 20:
						return {{4.054000e-01, 1.388000e+00, 2.160200e+00, 3.753200e+00, 2.206300e+00}, {3.499000e-01, 3.099100e+00, 1.196080e+01, 5.393530e+01, 1.423892e+02}};
					case 21:
						return {{3.787000e-01, 1.218100e+00, 2.059400e+00, 3.261800e+00, 2.387000e+00}, {3.133000e-01, 2.585600e+00, 9.581300e+00, 4.176880e+01, 1.167282e+02}};
					case 22:
						return {{3.825000e-01, 1.259800e+00, 2.000800e+00, 3.061700e+00, 2.069400e+00}, {3.040000e-01, 2.486300e+00, 9.278300e+00, 3.907510e+01, 1.094583e+02}};
					case 23:
						return {{3.876000e-01, 1.275000e+00, 1.910900e+00, 2.831400e+00, 1.897900e+00}, {2.967000e-01, 2.378000e+00, 8.798100e+00, 3.595280e+01, 1.017201e+02}};
					case 24:
						return {{4.046000e-01, 1.369600e+00, 1.894100e+00, 2.080000e+00, 1.219600e+00}, {2.986000e-01, 2.395800e+00, 9.140600e+00, 3.747010e+01, 1.137121e+02}};
					case 25:
						return {{3.796000e-01, 1.209400e+00, 1.781500e+00, 2.542000e+00, 1.593700e+00}, {2.699000e-01, 2.045500e+00, 7.472600e+00, 3.106040e+01, 9.156220e+01}};
					case 26:
						return {{3.946000e-01, 1.272500e+00, 1.703100e+00, 2.314000e+00, 1.479500e+00}, {2.717000e-01, 2.044300e+00, 7.600700e+00, 2.997140e+01, 8.622650e+01}};
					case 27:
						return {{4.118000e-01, 1.316100e+00, 1.649300e+00, 2.193000e+00, 1.283000e+00}, {2.742000e-01, 2.037200e+00, 7.720500e+00, 2.996800e+01, 8.493830e+01}};
					case 28:
						return {{3.860000e-01, 1.176500e+00, 1.545100e+00, 2.073000e+00, 1.381400e+00}, {2.478000e-01, 1.766000e+00, 6.310700e+00, 2.522040e+01, 7.431460e+01}};
					case 29:
						return {{4.314000e-01, 1.320800e+00, 1.523600e+00, 1.467100e+00, 8.562000e-01}, {2.694000e-01, 1.922300e+00, 7.347400e+00, 2.898920e+01, 9.062460e+01}};
					case 30:
						return {{4.288000e-01, 1.264600e+00, 1.447200e+00, 1.829400e+00, 1.093400e+00}, {2.593000e-01, 1.799800e+00, 6.750000e+00, 2.558600e+01, 7.352840e+01}};
					case 31:
						return {{4.818000e-01, 1.403200e+00, 1.656100e+00, 2.460500e+00, 1.105400e+00}, {2.825000e-01, 1.978500e+00, 8.754600e+00, 3.252380e+01, 9.855230e+01}};
					case 32:
						return {{4.655000e-01, 1.301400e+00, 1.608800e+00, 2.699800e+00, 1.300300e+00}, {2.647000e-01, 1.792600e+00, 7.607100e+00, 2.655410e+01, 7.752380e+01}};
					case 33:
						return {{4.517000e-01, 1.222900e+00, 1.585200e+00, 2.795800e+00, 1.263800e+00}, {2.493000e-01, 1.643600e+00, 6.815400e+00, 2.236810e+01, 6.203900e+01}};
					case 34:
						return {{4.477000e-01, 1.167800e+00, 1.584300e+00, 2.808700e+00, 1.195600e+00}, {2.405000e-01, 1.544200e+00, 6.323100e+00, 1.946100e+01, 5.202330e+01}};
					case 35:
						return {{4.798000e-01, 1.194800e+00, 1.869500e+00, 2.695300e+00, 8.203000e-01}, {2.504000e-01, 1.596300e+00, 6.965300e+00, 1.984920e+01, 5.032330e+01}};
					case 36:
						return {{4.546000e-01, 1.099300e+00, 1.769600e+00, 2.706800e+00, 8.672000e-01}, {2.309000e-01, 1.427900e+00, 5.944900e+00, 1.667520e+01, 4.222430e+01}};
					case 37:
						return {{1.016000e+00, 2.852800e+00, 3.546600e+00, -7.780400e+00, 1.211480e+01}, {4.853000e-01, 5.092500e+00, 2.578510e+01, 1.304515e+02, 1.386775e+02}};
					case 38:
						return {{6.703000e-01, 1.492600e+00, 3.336800e+00, 4.460000e+00, 3.150100e+00}, {3.190000e-01, 2.228700e+00, 1.035040e+01, 5.232910e+01, 1.512216e+02}};
					case 39:
						return {{6.894000e-01, 1.547400e+00, 3.245000e+00, 4.212600e+00, 2.976400e+00}, {3.189000e-01, 2.290400e+00, 1.000620e+01, 4.407710e+01, 1.250120e+02}};
					case 40:
						return {{6.719000e-01, 1.468400e+00, 3.166800e+00, 3.955700e+00, 2.892000e+00}, {3.036000e-01, 2.124900e+00, 8.923600e+00, 3.684580e+01, 1.082049e+02}};
					case 41:
						return {{6.123000e-01, 1.267700e+00, 3.034800e+00, 3.384100e+00, 2.368300e+00}, {2.709000e-01, 1.768300e+00, 7.248900e+00, 2.794650e+01, 9.856240e+01}};
					case 42:
						return {{6.773000e-01, 1.479800e+00, 3.178800e+00, 3.082400e+00, 1.838400e+00}, {2.920000e-01, 2.060600e+00, 8.112900e+00, 3.053360e+01, 1.000658e+02}};
					case 43:
						return {{7.082000e-01, 1.639200e+00, 3.199300e+00, 3.432700e+00, 1.871100e+00}, {2.976000e-01, 2.210600e+00, 8.524600e+00, 3.314560e+01, 9.663770e+01}};
					case 44:
						return {{6.735000e-01, 1.493400e+00, 3.096600e+00, 2.725400e+00, 1.559700e+00}, {2.773000e-01, 1.971600e+00, 7.324900e+00, 2.668910e+01, 9.055810e+01}};
					case 45:
						return {{6.413000e-01, 1.369000e+00, 2.985400e+00, 2.695200e+00, 1.543300e+00}, {2.580000e-01, 1.772100e+00, 6.385400e+00, 2.325490e+01, 8.515170e+01}};
					case 46:
						return {{5.904000e-01, 1.177500e+00, 2.651900e+00, 2.287500e+00, 8.689000e-01}, {2.324000e-01, 1.501900e+00, 5.159100e+00, 1.554280e+01, 4.682130e+01}};
					case 47:
						return {{6.377000e-01, 1.379000e+00, 2.829400e+00, 2.363100e+00, 1.455300e+00}, {2.466000e-01, 1.697400e+00, 5.765600e+00, 2.009430e+01, 7.673720e+01}};
					case 48:
						return {{6.364000e-01, 1.424700e+00, 2.780200e+00, 2.597300e+00, 1.788600e+00}, {2.407000e-01, 1.682300e+00, 5.658800e+00, 2.072190e+01, 6.911090e+01}};
					case 49:
						return {{6.768000e-01, 1.658900e+00, 2.774000e+00, 3.183500e+00, 2.132600e+00}, {2.522000e-01, 1.854500e+00, 6.293600e+00, 2.514570e+01, 8.454480e+01}};
					case 50:
						return {{7.224000e-01, 1.961000e+00, 2.716100e+00, 3.560300e+00, 1.897200e+00}, {2.651000e-01, 2.060400e+00, 7.301100e+00, 2.754930e+01, 8.133490e+01}};
					case 51:
						return {{7.106000e-01, 1.924700e+00, 2.614900e+00, 3.832200e+00, 1.889900e+00}, {2.562000e-01, 1.964600e+00, 6.885200e+00, 2.476480e+01, 6.891680e+01}};
					case 52:
						return {{6.947000e-01, 1.869000e+00, 2.535600e+00, 4.001300e+00, 1.895500e+00}, {2.459000e-01, 1.854200e+00, 6.441100e+00, 2.217300e+01, 5.922060e+01}};
					case 53:
						return {{7.047000e-01, 1.948400e+00, 2.594000e+00, 4.152600e+00, 1.505700e+00}, {2.455000e-01, 1.863800e+00, 6.763900e+00, 2.180070e+01, 5.643950e+01}};
					case 54:
						return {{6.737000e-01, 1.790800e+00, 2.412900e+00, 4.210000e+00, 1.705800e+00}, {2.305000e-01, 1.689000e+00, 5.821800e+00, 1.839280e+01, 4.724960e+01}};
					case 55:
						return {{1.270400e+00, 3.801800e+00, 5.661800e+00, 9.205000e-01, 4.810500e+00}, {4.356000e-01, 4.205800e+00, 2.343420e+01, 1.367783e+02, 1.717561e+02}};
					case 56:
						return {{9.049000e-01, 2.607600e+00, 4.849800e+00, 5.160300e+00, 4.738800e+00}, {3.066000e-01, 2.436300e+00, 1.218210e+01, 5.461350e+01, 1.619978e+02}};
					case 57:
						return {{8.405000e-01, 2.386300e+00, 4.613900e+00, 5.151400e+00, 4.794900e+00}, {2.791000e-01, 2.141000e+00, 1.034000e+01, 4.191480e+01, 1.320204e+02}};
					case 58:
						return {{8.551000e-01, 2.391500e+00, 4.577200e+00, 5.027800e+00, 4.511800e+00}, {2.805000e-01, 2.120000e+00, 1.018080e+01, 4.206330e+01, 1.309893e+02}};
					case 59:
						return {{9.096000e-01, 2.531300e+00, 4.526600e+00, 4.637600e+00, 4.369000e+00}, {2.939000e-01, 2.247100e+00, 1.082660e+01, 4.888420e+01, 1.476020e+02}};
					case 60:
						return {{8.807000e-01, 2.418300e+00, 4.444800e+00, 4.685800e+00, 4.172500e+00}, {2.802000e-01, 2.083600e+00, 1.003570e+01, 4.745060e+01, 1.469976e+02}};
					case 61:
						return {{9.471000e-01, 2.546300e+00, 4.352300e+00, 4.478900e+00, 3.908000e+00}, {2.977000e-01, 2.227600e+00, 1.057620e+01, 4.936190e+01, 1.453580e+02}};
					case 62:
						return {{9.699000e-01, 2.583700e+00, 4.277800e+00, 4.457500e+00, 3.598500e+00}, {3.003000e-01, 2.244700e+00, 1.064870e+01, 5.079940e+01, 1.464179e+02}};
					case 63:
						return {{8.694000e-01, 2.241300e+00, 3.919600e+00, 3.969400e+00, 4.549800e+00}, {2.653000e-01, 1.859000e+00, 8.399800e+00, 3.673970e+01, 1.257089e+02}};
					case 64:
						return {{9.673000e-01, 2.470200e+00, 4.114800e+00, 4.497200e+00, 3.209900e+00}, {2.909000e-01, 2.101400e+00, 9.706700e+00, 4.342700e+01, 1.259474e+02}};
					case 65:
						return {{9.325000e-01, 2.367300e+00, 3.879100e+00, 3.967400e+00, 3.799600e+00}, {2.761000e-01, 1.951100e+00, 8.929600e+00, 4.159370e+01, 1.310122e+02}};
					case 66:
						return {{9.505000e-01, 2.370500e+00, 3.821800e+00, 4.047100e+00, 3.445100e+00}, {2.773000e-01, 1.946900e+00, 8.886200e+00, 4.309380e+01, 1.331396e+02}};
					case 67:
						return {{9.248000e-01, 2.242800e+00, 3.618200e+00, 3.791000e+00, 3.791200e+00}, {2.660000e-01, 1.818300e+00, 7.965500e+00, 3.311290e+01, 1.018139e+02}};
					case 68:
						return {{1.037300e+00, 2.482400e+00, 3.655800e+00, 3.892500e+00, 3.005600e+00}, {2.944000e-01, 2.079700e+00, 9.415600e+00, 4.580560e+01, 1.327720e+02}};
					case 69:
						return {{1.007500e+00, 2.378700e+00, 3.544000e+00, 3.693200e+00, 3.175900e+00}, {2.816000e-01, 1.948600e+00, 8.716200e+00, 4.184200e+01, 1.250320e+02}};
					case 70:
						return {{1.034700e+00, 2.391100e+00, 3.461900e+00, 3.655600e+00, 3.005200e+00}, {2.855000e-01, 1.967900e+00, 8.761900e+00, 4.233040e+01, 1.256499e+02}};
					case 71:
						return {{9.927000e-01, 2.243600e+00, 3.355400e+00, 3.781300e+00, 3.099400e+00}, {2.701000e-01, 1.807300e+00, 7.811200e+00, 3.448490e+01, 1.033526e+02}};
					case 72:
						return {{1.029500e+00, 2.291100e+00, 3.411000e+00, 3.949700e+00, 2.492500e+00}, {2.761000e-01, 1.862500e+00, 8.096100e+00, 3.427120e+01, 9.852950e+01}};
					case 73:
						return {{1.019000e+00, 2.229100e+00, 3.409700e+00, 3.925200e+00, 2.267900e+00}, {2.694000e-01, 1.796200e+00, 7.694400e+00, 3.109420e+01, 9.110890e+01}};
					case 74:
						return {{9.853000e-01, 2.116700e+00, 3.357000e+00, 3.798100e+00, 2.279800e+00}, {2.569000e-01, 1.674500e+00, 7.009800e+00, 2.692340e+01, 8.139100e+01}};
					case 75:
						return {{9.914000e-01, 2.085800e+00, 3.453100e+00, 3.881200e+00, 1.852600e+00}, {2.548000e-01, 1.651800e+00, 6.884500e+00, 2.672340e+01, 8.172150e+01}};
					case 76:
						return {{9.813000e-01, 2.032200e+00, 3.366500e+00, 3.623500e+00, 1.974100e+00}, {2.487000e-01, 1.597300e+00, 6.473700e+00, 2.328170e+01, 7.092540e+01}};
					case 77:
						return {{1.019400e+00, 2.064500e+00, 3.442500e+00, 3.491400e+00, 1.697600e+00}, {2.554000e-01, 1.647500e+00, 6.596600e+00, 2.322690e+01, 7.002720e+01}};
					case 78:
						return {{9.148000e-01, 1.809600e+00, 3.213400e+00, 3.295300e+00, 1.575400e+00}, {2.263000e-01, 1.381300e+00, 5.324300e+00, 1.759870e+01, 6.001710e+01}};
					case 79:
						return {{9.674000e-01, 1.891600e+00, 3.399300e+00, 3.052400e+00, 1.260700e+00}, {2.358000e-01, 1.471200e+00, 5.675800e+00, 1.871190e+01, 6.152860e+01}};
					case 80:
						return {{1.003300e+00, 1.946900e+00, 3.439600e+00, 3.154800e+00, 1.418000e+00}, {2.413000e-01, 1.529800e+00, 5.800900e+00, 1.945200e+01, 6.057530e+01}};
					case 81:
						return {{1.068900e+00, 2.103800e+00, 3.603900e+00, 3.492700e+00, 1.828300e+00}, {2.540000e-01, 1.671500e+00, 6.350900e+00, 2.315310e+01, 7.870990e+01}};
					case 82:
						return {{1.089100e+00, 2.186700e+00, 3.616000e+00, 3.803100e+00, 1.899400e+00}, {2.552000e-01, 1.717400e+00, 6.513100e+00, 2.391700e+01, 7.470390e+01}};
					case 83:
						return {{1.100700e+00, 2.230600e+00, 3.568900e+00, 4.154900e+00, 2.038200e+00}, {2.546000e-01, 1.735100e+00, 6.494800e+00, 2.364640e+01, 7.037800e+01}};
					case 84:
						return {{1.156800e+00, 2.435300e+00, 3.645900e+00, 4.406400e+00, 1.717900e+00}, {2.648000e-01, 1.878600e+00, 7.174900e+00, 2.517660e+01, 6.928210e+01}};
					case 85:
						return {{1.090900e+00, 2.197600e+00, 3.383100e+00, 4.670000e+00, 2.127700e+00}, {2.466000e-01, 1.670700e+00, 6.019700e+00, 2.076570e+01, 5.726630e+01}};
					case 86:
						return {{1.075600e+00, 2.163000e+00, 3.317800e+00, 4.885200e+00, 2.048900e+00}, {2.402000e-01, 1.616900e+00, 5.764400e+00, 1.945680e+01, 5.250090e+01}};
					case 87:
						return {{1.428200e+00, 3.508100e+00, 5.676700e+00, 4.196400e+00, 3.894600e+00}, {3.183000e-01, 2.688900e+00, 1.348160e+01, 5.438660e+01, 2.008321e+02}};
					case 88:
						return {{1.312700e+00, 3.124300e+00, 5.298800e+00, 5.389100e+00, 5.413300e+00}, {2.887000e-01, 2.289700e+00, 1.082760e+01, 4.353890e+01, 1.456109e+02}};
					case 89:
						return {{1.312800e+00, 3.102100e+00, 5.338500e+00, 5.961100e+00, 4.756200e+00}, {2.861000e-01, 2.250900e+00, 1.052870e+01, 4.177960e+01, 1.282973e+02}};
					case 90:
						return {{1.255300e+00, 2.917800e+00, 5.086200e+00, 6.120600e+00, 4.712200e+00}, {2.701000e-01, 2.063600e+00, 9.305100e+00, 3.459770e+01, 1.079200e+02}};
					case 91:
						return {{1.321800e+00, 3.144400e+00, 5.437100e+00, 5.644400e+00, 4.010700e+00}, {2.827000e-01, 2.225000e+00, 1.024540e+01, 4.111620e+01, 1.244449e+02}};
					case 92:
						return {{1.338200e+00, 3.204300e+00, 5.455800e+00, 5.483900e+00, 3.634200e+00}, {2.838000e-01, 2.245200e+00, 1.025190e+01, 4.172510e+01, 1.249023e+02}};
					case 93:
						return {{1.519300e+00, 4.005300e+00, 6.532700e+00, -1.402000e-01, 6.748900e+00}, {3.213000e-01, 2.820600e+00, 1.488780e+01, 6.891030e+01, 8.172570e+01}};
					case 94:
						return {{1.351700e+00, 3.293700e+00, 5.321300e+00, 4.646600e+00, 3.571400e+00}, {2.813000e-01, 2.241800e+00, 9.995200e+00, 4.279390e+01, 1.321739e+02}};
					case 95:
						return {{1.213500e+00, 2.796200e+00, 4.754500e+00, 4.573100e+00, 4.478600e+00}, {2.483000e-01, 1.843700e+00, 7.542100e+00, 2.938410e+01, 1.124579e+02}};
					case 96:
						return {{1.293700e+00, 3.110000e+00, 5.039300e+00, 4.754600e+00, 3.503100e+00}, {2.638000e-01, 2.034100e+00, 8.710100e+00, 3.529920e+01, 1.094972e+02}};
					case 97:
						return {{1.291500e+00, 3.102300e+00, 4.930900e+00, 4.600900e+00, 3.466100e+00}, {2.611000e-01, 2.002300e+00, 8.437700e+00, 3.415590e+01, 1.058911e+02}};
					case 98:
						return {{1.208900e+00, 2.739100e+00, 4.348200e+00, 4.004700e+00, 4.649700e+00}, {2.421000e-01, 1.748700e+00, 6.726200e+00, 2.321530e+01, 8.031080e+01}};
				}

				return {{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}};
			}

			// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
			LNL_Coef_cpu<T> load_feg_peng_neutral_0_12(const dt_int32& Z)
			{
				switch(Z)
				{
					case 1:
						return {{8.800000e-03, 4.490000e-02, 1.481000e-01, 2.356000e-01, 9.140000e-02}, {1.152000e-01, 1.086700e+00, 4.975500e+00, 1.655910e+01, 4.327430e+01}};
					case 2:
						return {{8.400000e-03, 4.430000e-02, 1.314000e-01, 1.671000e-01, 6.660000e-02}, {5.960000e-02, 5.360000e-01, 2.427400e+00, 7.785200e+00, 2.031260e+01}};
					case 3:
						return {{4.780000e-02, 2.048000e-01, 5.253000e-01, 1.522500e+00, 9.853000e-01}, {2.258000e-01, 2.103200e+00, 1.293490e+01, 5.075010e+01, 1.366280e+02}};
					case 4:
						return {{4.230000e-02, 1.874000e-01, 6.019000e-01, 1.431100e+00, 7.891000e-01}, {1.445000e-01, 1.418000e+00, 8.116500e+00, 2.797050e+01, 7.486840e+01}};
					case 5:
						return {{4.360000e-02, 1.898000e-01, 6.788000e-01, 1.327300e+00, 5.544000e-01}, {1.207000e-01, 1.159500e+00, 6.247400e+00, 2.104600e+01, 5.936190e+01}};
					case 6:
						return {{4.890000e-02, 2.091000e-01, 7.537000e-01, 1.142000e+00, 3.555000e-01}, {1.140000e-01, 1.082500e+00, 5.428100e+00, 1.788110e+01, 5.113410e+01}};
					case 7:
						return {{2.670000e-02, 1.328000e-01, 5.301000e-01, 1.102000e+00, 4.215000e-01}, {5.410000e-02, 5.165000e-01, 2.820700e+00, 1.062970e+01, 3.437640e+01}};
					case 8:
						return {{3.650000e-02, 1.729000e-01, 5.805000e-01, 8.814000e-01, 3.121000e-01}, {6.520000e-02, 6.184000e-01, 2.944900e+00, 9.629800e+00, 2.821940e+01}};
					case 9:
						return {{3.820000e-02, 1.822000e-01, 5.972000e-01, 7.707000e-01, 2.130000e-01}, {6.130000e-02, 5.753000e-01, 2.685800e+00, 8.821400e+00, 2.566680e+01}};
					case 10:
						return {{3.800000e-02, 1.785000e-01, 5.494000e-01, 6.942000e-01, 1.918000e-01}, {5.540000e-02, 5.087000e-01, 2.263900e+00, 7.331600e+00, 2.169120e+01}};
					case 11:
						return {{1.260000e-01, 6.442000e-01, 8.893000e-01, 1.819700e+00, 1.298800e+00}, {1.684000e-01, 1.715000e+00, 8.838600e+00, 5.082650e+01, 1.472073e+02}};
					case 12:
						return {{1.130000e-01, 5.575000e-01, 9.046000e-01, 2.158000e+00, 1.473500e+00}, {1.356000e-01, 1.357900e+00, 6.925500e+00, 3.231650e+01, 9.211380e+01}};
					case 13:
						return {{1.165000e-01, 5.504000e-01, 1.017900e+00, 2.629500e+00, 1.571100e+00}, {1.295000e-01, 1.261900e+00, 6.824200e+00, 2.845770e+01, 8.847500e+01}};
					case 14:
						return {{5.670000e-02, 3.365000e-01, 8.104000e-01, 2.496000e+00, 2.118600e+00}, {5.820000e-02, 6.155000e-01, 3.252200e+00, 1.679290e+01, 5.767670e+01}};
					case 15:
						return {{1.005000e-01, 4.615000e-01, 1.066300e+00, 2.585400e+00, 1.272500e+00}, {9.770000e-02, 9.084000e-01, 4.965400e+00, 1.854710e+01, 5.436480e+01}};
					case 16:
						return {{9.150000e-02, 4.312000e-01, 1.084700e+00, 2.467100e+00, 1.085200e+00}, {8.380000e-02, 7.788000e-01, 4.346200e+00, 1.558460e+01, 4.463650e+01}};
					case 17:
						return {{7.990000e-02, 3.891000e-01, 1.003700e+00, 2.333200e+00, 1.050700e+00}, {6.940000e-02, 6.443000e-01, 3.535100e+00, 1.250580e+01, 3.586330e+01}};
					case 18:
						return {{1.044000e-01, 4.551000e-01, 1.423200e+00, 2.153300e+00, 4.459000e-01}, {8.530000e-02, 7.701000e-01, 4.468400e+00, 1.458640e+01, 4.124740e+01}};
					case 19:
						return {{2.149000e-01, 8.703000e-01, 2.499900e+00, 2.359100e+00, 3.031800e+00}, {1.660000e-01, 1.690600e+00, 8.744700e+00, 4.678250e+01, 1.656923e+02}};
					case 20:
						return {{2.355000e-01, 9.916000e-01, 2.395900e+00, 3.725200e+00, 2.564700e+00}, {1.742000e-01, 1.832900e+00, 8.840700e+00, 4.745830e+01, 1.349613e+02}};
					case 21:
						return {{4.636000e-01, 2.080200e+00, 2.900300e+00, 1.419300e+00, 2.432300e+00}, {3.682000e-01, 4.031200e+00, 2.264930e+01, 7.182000e+01, 1.033691e+02}};
					case 22:
						return {{2.123000e-01, 8.960000e-01, 2.176500e+00, 3.043600e+00, 2.443900e+00}, {1.399000e-01, 1.456800e+00, 6.753400e+00, 3.311680e+01, 1.018238e+02}};
					case 23:
						return {{2.369000e-01, 1.077400e+00, 2.189400e+00, 3.082500e+00, 1.719000e+00}, {1.505000e-01, 1.639200e+00, 7.569100e+00, 3.687410e+01, 1.078517e+02}};
					case 24:
						return {{1.970000e-01, 8.228000e-01, 2.020000e+00, 2.171700e+00, 1.751600e+00}, {1.197000e-01, 1.198500e+00, 5.409700e+00, 2.523610e+01, 9.442900e+01}};
					case 25:
						return {{1.943000e-01, 8.190000e-01, 1.929600e+00, 2.496800e+00, 2.062500e+00}, {1.135000e-01, 1.131300e+00, 5.034100e+00, 2.417980e+01, 8.055980e+01}};
					case 26:
						return {{1.929000e-01, 8.239000e-01, 1.868900e+00, 2.369400e+00, 1.906000e+00}, {1.087000e-01, 1.080600e+00, 4.763700e+00, 2.285000e+01, 7.673090e+01}};
					case 27:
						return {{2.186000e-01, 9.861000e-01, 1.854000e+00, 2.325800e+00, 1.468500e+00}, {1.182000e-01, 1.230000e+00, 5.417700e+00, 2.576020e+01, 8.085420e+01}};
					case 28:
						return {{2.313000e-01, 1.065700e+00, 1.822900e+00, 2.260900e+00, 1.188300e+00}, {1.210000e-01, 1.269100e+00, 5.687000e+00, 2.709170e+01, 8.302850e+01}};
					case 29:
						return {{3.501000e-01, 1.655800e+00, 1.958200e+00, 2.134000e-01, 1.410900e+00}, {1.867000e-01, 1.991700e+00, 1.133960e+01, 5.326190e+01, 6.325200e+01}};
					case 30:
						return {{1.780000e-01, 8.096000e-01, 1.674400e+00, 1.949900e+00, 1.449500e+00}, {8.760000e-02, 8.650000e-01, 3.861200e+00, 1.887260e+01, 6.470160e+01}};
					case 31:
						return {{2.135000e-01, 9.768000e-01, 1.666900e+00, 2.566200e+00, 1.679000e+00}, {1.020000e-01, 1.021900e+00, 4.627500e+00, 2.287420e+01, 8.015350e+01}};
					case 32:
						return {{2.135000e-01, 9.761000e-01, 1.655500e+00, 2.893800e+00, 1.635600e+00}, {9.890000e-02, 9.845000e-01, 4.552700e+00, 2.155630e+01, 7.039030e+01}};
					case 33:
						return {{2.059000e-01, 9.518000e-01, 1.637200e+00, 3.049000e+00, 1.475600e+00}, {9.260000e-02, 9.182000e-01, 4.329100e+00, 1.929960e+01, 5.893290e+01}};
					case 34:
						return {{1.574000e-01, 7.614000e-01, 1.483400e+00, 3.001600e+00, 1.797800e+00}, {6.860000e-02, 6.808000e-01, 3.116300e+00, 1.434580e+01, 4.404550e+01}};
					case 35:
						return {{1.899000e-01, 8.983000e-01, 1.635800e+00, 3.184500e+00, 1.151800e+00}, {8.100000e-02, 7.957000e-01, 3.905400e+00, 1.577010e+01, 4.561240e+01}};
					case 36:
						return {{1.742000e-01, 8.447000e-01, 1.594400e+00, 3.150700e+00, 1.133800e+00}, {7.230000e-02, 7.123000e-01, 3.519200e+00, 1.377240e+01, 3.911480e+01}};
					case 37:
						return {{3.781000e-01, 1.490400e+00, 3.575300e+00, 3.003100e+00, 3.327200e+00}, {1.557000e-01, 1.534700e+00, 9.994700e+00, 5.142510e+01, 1.859828e+02}};
					case 38:
						return {{3.723000e-01, 1.459800e+00, 3.512400e+00, 4.461200e+00, 3.303100e+00}, {1.480000e-01, 1.464300e+00, 9.232000e+00, 4.988070e+01, 1.480937e+02}};
					case 39:
						return {{3.234000e-01, 1.273700e+00, 3.211500e+00, 4.056300e+00, 3.796200e+00}, {1.244000e-01, 1.194800e+00, 7.275600e+00, 3.414300e+01, 1.112079e+02}};
					case 40:
						return {{2.997000e-01, 1.187900e+00, 3.107500e+00, 3.974000e+00, 3.576900e+00}, {1.121000e-01, 1.063800e+00, 6.389100e+00, 2.870810e+01, 9.742890e+01}};
					case 41:
						return {{1.680000e-01, 9.370000e-01, 2.730000e+00, 3.815000e+00, 3.005300e+00}, {5.970000e-02, 6.524000e-01, 4.431700e+00, 1.955400e+01, 8.550110e+01}};
					case 42:
						return {{3.069000e-01, 1.171400e+00, 3.229300e+00, 3.425400e+00, 2.122400e+00}, {1.101000e-01, 1.022200e+00, 5.961300e+00, 2.519650e+01, 9.358310e+01}};
					case 43:
						return {{2.928000e-01, 1.126700e+00, 3.167500e+00, 3.661900e+00, 2.594200e+00}, {1.020000e-01, 9.481000e-01, 5.471300e+00, 2.381530e+01, 8.289910e+01}};
					case 44:
						return {{2.604000e-01, 1.044200e+00, 3.076100e+00, 3.217500e+00, 1.944800e+00}, {8.870000e-02, 8.240000e-01, 4.827800e+00, 1.989770e+01, 8.045660e+01}};
					case 45:
						return {{2.713000e-01, 1.055600e+00, 3.141600e+00, 3.045100e+00, 1.717900e+00}, {9.070000e-02, 8.324000e-01, 4.770200e+00, 1.978620e+01, 8.025400e+01}};
					case 46:
						return {{2.003000e-01, 8.779000e-01, 2.613500e+00, 2.859400e+00, 1.025800e+00}, {6.590000e-02, 6.111000e-01, 3.556300e+00, 1.276380e+01, 4.442830e+01}};
					case 47:
						return {{2.739000e-01, 1.050300e+00, 3.156400e+00, 2.754300e+00, 1.432800e+00}, {8.810000e-02, 8.028000e-01, 4.445100e+00, 1.870110e+01, 7.926330e+01}};
					case 48:
						return {{3.072000e-01, 1.130300e+00, 3.204600e+00, 2.932900e+00, 1.656000e+00}, {9.660000e-02, 8.856000e-01, 4.627300e+00, 2.067890e+01, 7.347230e+01}};
					case 49:
						return {{3.564000e-01, 1.301100e+00, 3.242400e+00, 3.483900e+00, 2.045900e+00}, {1.091000e-01, 1.045200e+00, 5.090000e+00, 2.465780e+01, 8.805130e+01}};
					case 50:
						return {{2.966000e-01, 1.115700e+00, 3.097300e+00, 3.815600e+00, 2.528100e+00}, {8.960000e-02, 8.268000e-01, 4.224200e+00, 2.069000e+01, 7.133990e+01}};
					case 51:
						return {{2.725000e-01, 1.065100e+00, 2.994000e+00, 4.069700e+00, 2.568200e+00}, {8.090000e-02, 7.488000e-01, 3.871000e+00, 1.888000e+01, 6.064990e+01}};
					case 52:
						return {{2.422000e-01, 9.692000e-01, 2.811400e+00, 4.150900e+00, 2.816100e+00}, {7.080000e-02, 6.472000e-01, 3.360900e+00, 1.607520e+01, 5.017240e+01}};
					case 53:
						return {{2.617000e-01, 1.032500e+00, 2.809700e+00, 4.480900e+00, 2.319000e+00}, {7.490000e-02, 6.914000e-01, 3.463400e+00, 1.636030e+01, 4.825220e+01}};
					case 54:
						return {{2.334000e-01, 9.496000e-01, 2.638100e+00, 4.468000e+00, 2.502000e+00}, {6.550000e-02, 6.050000e-01, 3.038900e+00, 1.408090e+01, 4.100050e+01}};
					case 55:
						return {{5.713000e-01, 2.486600e+00, 4.979500e+00, 4.019800e+00, 4.440300e+00}, {1.626000e-01, 1.821300e+00, 1.110490e+01, 4.905680e+01, 2.029987e+02}};
					case 56:
						return {{5.229000e-01, 2.287400e+00, 4.724300e+00, 5.080700e+00, 5.638900e+00}, {1.434000e-01, 1.601900e+00, 9.451100e+00, 4.276850e+01, 1.484969e+02}};
					case 57:
						return {{5.461000e-01, 2.385600e+00, 5.065300e+00, 5.760100e+00, 4.046300e+00}, {1.479000e-01, 1.655200e+00, 1.000590e+01, 4.732450e+01, 1.458464e+02}};
					case 58:
						return {{2.227000e-01, 1.076000e+00, 2.948200e+00, 5.849600e+00, 7.183400e+00}, {5.710000e-02, 5.946000e-01, 3.202200e+00, 1.642530e+01, 9.570300e+01}};
					case 59:
						return {{5.237000e-01, 2.291300e+00, 4.616100e+00, 4.723300e+00, 4.817300e+00}, {1.360000e-01, 1.506800e+00, 8.821300e+00, 4.195360e+01, 1.412424e+02}};
					case 60:
						return {{5.368000e-01, 2.330100e+00, 4.605800e+00, 4.662100e+00, 4.462200e+00}, {1.378000e-01, 1.514000e+00, 8.871900e+00, 4.359670e+01, 1.418065e+02}};
					case 61:
						return {{5.232000e-01, 2.262700e+00, 4.455200e+00, 4.478700e+00, 4.507300e+00}, {1.317000e-01, 1.433600e+00, 8.308700e+00, 4.060100e+01, 1.359196e+02}};
					case 62:
						return {{5.162000e-01, 2.230200e+00, 4.344900e+00, 4.359800e+00, 4.429200e+00}, {1.279000e-01, 1.381100e+00, 7.962900e+00, 3.912130e+01, 1.327846e+02}};
					case 63:
						return {{5.272000e-01, 2.284400e+00, 4.336100e+00, 4.317800e+00, 4.090800e+00}, {1.285000e-01, 1.394300e+00, 8.108100e+00, 4.096310e+01, 1.341233e+02}};
					case 64:
						return {{9.664000e-01, 3.405200e+00, 5.080300e+00, 1.499100e+00, 4.252800e+00}, {2.641000e-01, 2.658600e+00, 1.622130e+01, 8.020600e+01, 9.253590e+01}};
					case 65:
						return {{5.110000e-01, 2.157000e+00, 4.030800e+00, 3.993600e+00, 4.246600e+00}, {1.210000e-01, 1.270400e+00, 7.136800e+00, 3.503540e+01, 1.235062e+02}};
					case 66:
						return {{4.974000e-01, 2.109700e+00, 3.890600e+00, 3.810000e+00, 4.308400e+00}, {1.157000e-01, 1.210800e+00, 6.737700e+00, 3.241500e+01, 1.169225e+02}};
					case 67:
						return {{4.679000e-01, 1.969300e+00, 3.719100e+00, 3.963200e+00, 4.243200e+00}, {1.069000e-01, 1.099400e+00, 5.976900e+00, 2.714910e+01, 9.631190e+01}};
					case 68:
						return {{5.034000e-01, 2.108800e+00, 3.823200e+00, 3.729900e+00, 3.896300e+00}, {1.141000e-01, 1.176900e+00, 6.608700e+00, 3.343320e+01, 1.164913e+02}};
					case 69:
						return {{4.839000e-01, 2.026200e+00, 3.685100e+00, 3.587400e+00, 4.003700e+00}, {1.081000e-01, 1.101200e+00, 6.111400e+00, 3.037280e+01, 1.105988e+02}};
					case 70:
						return {{5.221000e-01, 2.169500e+00, 3.756700e+00, 3.668500e+00, 3.427400e+00}, {1.148000e-01, 1.186000e+00, 6.752000e+00, 3.568070e+01, 1.180692e+02}};
					case 71:
						return {{4.680000e-01, 1.946600e+00, 3.542800e+00, 3.849000e+00, 3.659400e+00}, {1.015000e-01, 1.019500e+00, 5.605800e+00, 2.748990e+01, 9.528460e+01}};
					case 72:
						return {{4.048000e-01, 1.737000e+00, 3.339900e+00, 3.944800e+00, 3.729300e+00}, {8.680000e-02, 8.585000e-01, 4.637800e+00, 2.169000e+01, 8.024080e+01}};
					case 73:
						return {{3.835000e-01, 1.674700e+00, 3.298600e+00, 4.046200e+00, 3.430300e+00}, {8.100000e-02, 8.020000e-01, 4.354500e+00, 1.996440e+01, 7.363370e+01}};
					case 74:
						return {{3.661000e-01, 1.619100e+00, 3.245500e+00, 4.085600e+00, 3.206400e+00}, {7.610000e-02, 7.543000e-01, 4.095200e+00, 1.828860e+01, 6.809670e+01}};
					case 75:
						return {{3.933000e-01, 1.697300e+00, 3.420200e+00, 4.127400e+00, 2.615800e+00}, {8.060000e-02, 7.972000e-01, 4.423700e+00, 1.956920e+01, 6.874770e+01}};
					case 76:
						return {{3.854000e-01, 1.655500e+00, 3.412900e+00, 4.111100e+00, 2.410600e+00}, {7.870000e-02, 7.638000e-01, 4.244100e+00, 1.837000e+01, 6.510710e+01}};
					case 77:
						return {{3.510000e-01, 1.562000e+00, 3.294600e+00, 4.061500e+00, 2.438200e+00}, {7.060000e-02, 6.904000e-01, 3.826600e+00, 1.608120e+01, 5.876380e+01}};
					case 78:
						return {{3.083000e-01, 1.415800e+00, 2.966200e+00, 3.934900e+00, 2.170900e+00}, {6.090000e-02, 5.993000e-01, 3.192100e+00, 1.252850e+01, 4.976750e+01}};
					case 79:
						return {{3.055000e-01, 1.394500e+00, 2.961700e+00, 3.899000e+00, 2.002600e+00}, {5.960000e-02, 5.827000e-01, 3.103500e+00, 1.196930e+01, 4.791060e+01}};
					case 80:
						return {{3.593000e-01, 1.573600e+00, 3.523700e+00, 3.810900e+00, 1.695300e+00}, {6.940000e-02, 6.758000e-01, 3.845700e+00, 1.562030e+01, 5.666140e+01}};
					case 81:
						return {{3.511000e-01, 1.548900e+00, 3.567600e+00, 4.090000e+00, 2.525100e+00}, {6.720000e-02, 6.522000e-01, 3.742000e+00, 1.597910e+01, 6.513540e+01}};
					case 82:
						return {{3.540000e-01, 1.545300e+00, 3.597500e+00, 4.315200e+00, 2.774300e+00}, {6.680000e-02, 6.465000e-01, 3.696800e+00, 1.620560e+01, 6.149090e+01}};
					case 83:
						return {{3.530000e-01, 1.525800e+00, 3.581500e+00, 4.553200e+00, 3.071400e+00}, {6.610000e-02, 6.324000e-01, 3.590600e+00, 1.599620e+01, 5.757600e+01}};
					case 84:
						return {{3.673000e-01, 1.577200e+00, 3.707900e+00, 4.858200e+00, 2.844000e+00}, {6.780000e-02, 6.527000e-01, 3.739600e+00, 1.706680e+01, 5.597890e+01}};
					case 85:
						return {{3.547000e-01, 1.520600e+00, 3.562100e+00, 5.018400e+00, 3.007500e+00}, {6.490000e-02, 6.188000e-01, 3.469600e+00, 1.560900e+01, 4.948180e+01}};
					case 86:
						return {{4.586000e-01, 1.778100e+00, 3.987700e+00, 5.727300e+00, 1.546000e+00}, {8.310000e-02, 7.840000e-01, 4.359900e+00, 2.001280e+01, 6.215350e+01}};
					case 87:
						return {{8.282000e-01, 2.994100e+00, 5.659700e+00, 4.929200e+00, 4.288900e+00}, {1.515000e-01, 1.616300e+00, 9.775200e+00, 4.284800e+01, 1.907366e+02}};
					case 88:
						return {{1.412900e+00, 4.426900e+00, 7.046000e+00, -1.057300e+00, 8.643000e+00}, {2.921000e-01, 3.138100e+00, 1.967670e+01, 1.020436e+02, 1.139798e+02}};
					case 89:
						return {{7.169000e-01, 2.571000e+00, 5.179100e+00, 6.348400e+00, 5.647400e+00}, {1.263000e-01, 1.290000e+00, 7.368600e+00, 3.244900e+01, 1.180558e+02}};
					case 90:
						return {{6.958000e-01, 2.493600e+00, 5.126900e+00, 6.698800e+00, 5.079900e+00}, {1.211000e-01, 1.224700e+00, 6.939800e+00, 3.009910e+01, 1.051960e+02}};
					case 91:
						return {{1.250200e+00, 4.228400e+00, 7.048900e+00, 1.139000e+00, 5.822200e+00}, {2.415000e-01, 2.644200e+00, 1.633130e+01, 7.357570e+01, 9.194010e+01}};
					case 92:
						return {{6.410000e-01, 2.264300e+00, 4.871300e+00, 5.928700e+00, 5.393500e+00}, {1.097000e-01, 1.064400e+00, 5.790700e+00, 2.502610e+01, 1.013899e+02}};
					case 93:
						return {{6.938000e-01, 2.465200e+00, 5.122700e+00, 5.596500e+00, 4.854300e+00}, {1.171000e-01, 1.175700e+00, 6.405300e+00, 2.752170e+01, 1.030482e+02}};
					case 94:
						return {{6.902000e-01, 2.450900e+00, 5.128400e+00, 5.033900e+00, 4.857500e+00}, {1.153000e-01, 1.154500e+00, 6.229100e+00, 2.707410e+01, 1.113150e+02}};
					case 95:
						return {{7.577000e-01, 2.726400e+00, 5.418400e+00, 4.819800e+00, 4.101300e+00}, {1.257000e-01, 1.304400e+00, 7.103500e+00, 3.246490e+01, 1.188647e+02}};
					case 96:
						return {{7.567000e-01, 2.756500e+00, 5.436400e+00, 5.191800e+00, 3.564300e+00}, {1.239000e-01, 1.297900e+00, 7.079800e+00, 3.278710e+01, 1.101512e+02}};
					case 97:
						return {{7.492000e-01, 2.726700e+00, 5.352100e+00, 5.036900e+00, 3.532100e+00}, {1.217000e-01, 1.265100e+00, 6.810100e+00, 3.160880e+01, 1.064853e+02}};
					case 98:
						return {{8.100000e-01, 3.000100e+00, 5.463500e+00, 4.175600e+00, 3.506600e+00}, {1.310000e-01, 1.403800e+00, 7.605700e+00, 3.401860e+01, 9.052260e+01}};
				}

				return {{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}};
			}

			// 4: Kirkland parameterization - 3 yukawa + 3 Gaussians - [0, 12]
			LNL_Coef_cpu<T> load_feg_kirkland_neutral_0_12(const dt_int32& Z)
			{
				switch(Z)
				{
					case 1:
						return {{4.20298320e-03, 6.27762505e-02, 3.00907347e-02, 6.77756695e-02, 3.56609240e-03, 2.76135815e-02}, {2.25350888e-01, 2.25366950e-01, 2.25331756e-01, 4.38854001e+00, 4.03884823e-01, 1.44490166e+00}};
					case 2:
						return {{1.87544000e-05, 4.10595800e-04, 1.96300059e-01, 8.36015740e-03, 2.95102022e-02, 4.65900000e-07}, {2.12427997e-01, 3.32212279e-01, 5.17325152e-01, 3.66668239e-01, 1.37171827e+00, 3.75768025e+04}};
					case 3:
						return {{7.45843816e-02, 7.15382250e-02, 1.45315229e-01, 1.12125769e+00, 2.51736530e-03, 3.58434971e-01}, {8.81151424e-01, 4.59142904e-02, 8.81301714e-01, 1.88483665e+01, 1.59189995e-01, 6.12371000e+00}};
					case 4:
						return {{6.11642897e-02, 1.25755034e-01, 2.00831548e-01, 7.87242876e-01, 1.58847850e-03, 2.73962031e-01}, {9.90182132e-02, 9.90272412e-02, 1.87392509e+00, 9.32794929e+00, 8.91900236e-02, 3.20687658e+00}};
					case 5:
						return {{1.25716066e-01, 1.73314452e-01, 1.84774811e-01, 1.95250221e-01, 5.29642075e-01, 1.08230500e-03}, {1.48258830e-01, 1.48257216e-01, 3.34227311e+00, 1.97339463e+00, 5.70035553e+00, 5.64857237e-02}};
					case 6:
						return {{2.12080767e-01, 1.99811865e-01, 1.68254385e-01, 1.42048360e-01, 3.63830672e-01, 8.35012000e-04}, {2.08605417e-01, 2.08610186e-01, 5.57870773e+00, 1.33311887e+00, 3.80800263e+00, 4.03982620e-02}};
					case 7:
						return {{5.33015554e-01, 5.29008883e-02, 9.24159648e-02, 2.61799101e-01, 8.80262100e-04, 1.10166555e-01}, {2.90952515e-01, 1.03547896e+01, 1.03540028e+01, 2.76252723e+00, 3.47681236e-02, 9.93421736e-01}};
					case 8:
						return {{3.39969204e-01, 3.07570172e-01, 1.30369072e-01, 8.83326058e-02, 1.96586700e-01, 9.96220000e-04}, {3.81570280e-01, 3.81571436e-01, 1.91919745e+01, 7.60635525e-01, 2.07401094e+00, 3.03266869e-02}};
					case 9:
						return {{2.30560593e-01, 5.26889648e-01, 1.24346755e-01, 1.24616890e-03, 7.20452555e-02, 1.53075777e-01}, {4.80754213e-01, 4.80763895e-01, 3.95306720e+01, 2.62181803e-02, 5.92495593e-01, 1.59127671e+00}};
					case 10:
						return {{4.08371771e-01, 4.54418858e-01, 1.44564923e-01, 5.91531395e-02, 1.24003718e-01, 1.64986040e-03}, {5.88228627e-01, 5.88288655e-01, 1.21246013e+02, 4.63963540e-01, 1.23413025e+00, 2.05869217e-02}};
					case 11:
						return {{1.36471662e-01, 7.70677865e-01, 1.56862014e-01, 9.96821513e-01, 3.80304670e-02, 1.27685089e-01}, {4.99965301e-02, 8.81899664e-01, 1.61768579e+01, 2.00132610e+01, 2.60516254e-01, 6.99559329e-01}};
					case 12:
						return {{3.04384121e-01, 7.56270563e-01, 1.01164809e-01, 3.45203403e-02, 9.71751327e-01, 1.20593012e-01}, {8.42014377e-02, 1.64065598e+00, 2.97142975e+01, 2.16596094e-01, 1.21236852e+01, 5.60865838e-01}};
					case 13:
						return {{7.77419424e-01, 5.78312036e-02, 4.26386499e-01, 1.13407220e-01, 7.90114035e-01, 3.23293496e-02}, {2.71058227e+00, 7.17532098e+01, 9.13331555e-02, 4.48867451e-01, 8.66366718e+00, 1.78503463e-01}};
					case 14:
						return {{1.06543892e+00, 1.20143691e-01, 1.80915263e-01, 1.12065620e+00, 3.05452816e-02, 1.59963502e+00}, {1.04118455e+00, 6.87113368e+01, 8.87533926e-02, 3.70062619e+00, 2.14097897e-01, 9.99096638e+00}};
					case 15:
						return {{1.05284447e+00, 2.99440284e-01, 1.17460748e-01, 9.60643452e-01, 2.63555748e-02, 1.38059330e+00}, {1.31962590e+00, 1.28460520e-01, 1.02190163e+02, 2.87477555e+00, 1.82076844e-01, 7.49165526e+00}};
					case 16:
						return {{1.01646916e+00, 4.41766748e-01, 1.21503863e-01, 8.27966670e-01, 2.33022533e-02, 1.18302846e+00}, {1.69181965e+00, 1.74180288e-01, 1.67011091e+02, 2.30342810e+00, 1.56954150e-01, 5.85782891e+00}};
					case 17:
						return {{9.44221116e-01, 4.37322049e-01, 2.54547926e-01, 5.47763323e-02, 8.00087488e-01, 1.07488641e-02}, {2.40052374e-01, 9.30510439e+00, 9.30486346e+00, 1.68655688e-01, 2.97849774e+00, 6.84240646e-02}};
					case 18:
						return {{1.06983288e+00, 4.24631786e-01, 2.43897949e-01, 4.79446296e-02, 7.64958952e-01, 8.23128430e-03}, {2.87791022e-01, 1.24156957e+01, 1.24158868e+01, 1.36979796e-01, 2.43940729e+00, 5.27258749e-02}};
					case 19:
						return {{6.92717865e-01, 9.65161085e-01, 1.48466588e-01, 2.64645027e-02, 1.80883768e+00, 5.43900018e-01}, {7.10849990e+00, 3.57532901e-01, 3.93763275e-02, 1.03591321e-01, 3.22845199e+01, 1.67791374e+00}};
					case 20:
						return {{3.66902871e-01, 8.66378999e-01, 6.67203300e-01, 4.87743636e-01, 1.82406314e+00, 2.20248453e-02}, {6.14274129e-02, 5.70881727e-01, 7.82965639e+00, 1.32531318e+00, 2.10056032e+01, 9.11853450e-02}};
					case 21:
						return {{3.78871777e-01, 9.00022505e-01, 7.15288914e-01, 1.88640973e-02, 4.07945949e-01, 1.61786540e+00}, {6.98910162e-02, 5.21061541e-01, 7.87707920e+00, 8.17512708e-02, 1.11141388e+00, 1.80840759e+01}};
					case 22:
						return {{3.62383267e-01, 9.84232966e-01, 7.41715642e-01, 3.62555269e-01, 1.49159390e+00, 1.61659509e-02}, {7.54707114e-02, 4.97757309e-01, 8.17659391e+00, 9.55524906e-01, 1.62221677e+01, 7.33140839e-02}};
					case 23:
						return {{3.52961378e-01, 7.46791014e-01, 1.08364068e+00, 1.39013610e+00, 3.31273356e-01, 1.40422612e-02}, {8.19204103e-02, 8.81189511e+00, 5.10646075e-01, 1.48901841e+01, 8.38543079e-01, 6.57432678e-02}};
					case 24:
						return {{1.34348379e+00, 5.07040328e-01, 4.26358955e-01, 1.17241826e-02, 5.11966516e-01, 3.38285828e-01}, {1.25814353e+00, 1.15042811e+01, 8.53660389e-02, 6.00177061e-02, 1.53772451e+00, 6.62418319e-01}};
					case 25:
						return {{3.26697613e-01, 7.17297000e-01, 1.33212464e+00, 2.80801702e-01, 1.15499241e+00, 1.11984488e-02}, {8.88813083e-02, 1.11300198e+01, 5.82141104e-01, 6.71583145e-01, 1.26825395e+01, 5.32334467e-02}};
					case 26:
						return {{3.13454847e-01, 6.89290016e-01, 1.47141531e+00, 1.03298688e+00, 2.58280285e-01, 1.03460690e-02}, {8.99325756e-02, 1.30366038e+01, 6.33345291e-01, 1.16783425e+01, 6.09116446e-01, 4.81610627e-02}};
					case 27:
						return {{3.15878278e-01, 1.60139005e+00, 6.56394338e-01, 9.36746624e-01, 9.77562650e-03, 2.38378578e-01}, {9.46683246e-02, 6.99436449e-01, 1.56954403e+01, 1.09392410e+01, 4.37446816e-02, 5.56286483e-01}};
					case 28:
						return {{1.72254630e+00, 3.29543044e-01, 6.23007200e-01, 9.43496510e-03, 8.54063515e-01, 2.21073515e-01}, {7.76606908e-01, 1.02262360e-01, 1.94156207e+01, 3.98684596e-02, 1.04078166e+01, 5.10869330e-01}};
					case 29:
						return {{3.58774531e-01, 1.76181348e+00, 6.36905053e-01, 7.44930670e-03, 1.89002347e-01, 2.29619589e-01}, {1.06153463e-01, 1.01640995e+00, 1.53659093e+01, 3.85345989e-02, 3.98427790e-01, 9.01419843e-01}};
					case 30:
						return {{5.70893973e-01, 1.98908856e+00, 3.06060585e-01, 2.35600223e-01, 3.97061102e-01, 6.85657230e-03}, {1.26534614e-01, 2.17781965e+00, 3.78619003e+01, 3.67019041e-01, 8.66419596e-01, 3.35778823e-02}};
					case 31:
						return {{6.25528464e-01, 2.05302901e+00, 2.89608120e-01, 2.07910594e-01, 3.45079617e-01, 6.55634300e-03}, {1.10005650e-01, 2.41095786e+00, 4.78685736e+01, 3.27807224e-01, 7.43139061e-01, 3.09411369e-02}};
					case 32:
						return {{5.90952690e-01, 5.39980660e-01, 2.00626188e+00, 7.49705041e-01, 1.83581347e-01, 9.52190740e-03}, {1.18375976e-01, 7.18937433e+01, 1.39304889e+00, 6.89943350e+00, 3.64667232e-01, 2.69888650e-02}};
					case 33:
						return {{7.77875218e-01, 5.93848150e-01, 1.95918751e+00, 1.79880226e-01, 8.63267222e-01, 9.59053430e-03}, {1.50733157e-01, 1.42882209e+02, 1.74750339e+00, 3.31800852e-01, 5.85490274e+00, 2.33777569e-02}};
					case 34:
						return {{9.58390681e-01, 6.03851342e-01, 1.90828931e+00, 1.73885956e-01, 9.35265145e-01, 8.62254660e-03}, {1.83775557e-01, 1.96819224e+02, 2.15082053e+00, 3.00006024e-01, 4.92471215e+00, 2.12308108e-02}};
					case 35:
						return {{1.14136170e+00, 5.18118737e-01, 1.85731975e+00, 1.68217399e-01, 9.75705606e-01, 7.24187870e-03}, {2.18708710e-01, 1.93916682e+02, 2.65755396e+00, 2.71719918e-01, 4.19482500e+00, 1.99325718e-02}};
					case 36:
						return {{3.24386970e-01, 1.31732163e+00, 1.79912614e+00, 4.29961430e-03, 1.00429433e+00, 1.62188197e-01}, {6.31317973e+01, 2.54706036e-01, 3.23668394e+00, 1.98965610e-02, 3.61094513e+00, 2.45583672e-01}};
					case 37:
						return {{2.90445351e-01, 2.44201329e+00, 7.69435449e-01, 1.58687000e+00, 2.81617590e-03, 1.28663830e-01}, {3.68420227e-02, 1.16013332e+00, 1.69591472e+01, 2.53082574e+00, 1.88577417e-02, 2.10753969e-01}};
					case 38:
						return {{1.37373086e-02, 1.97548672e+00, 1.59261029e+00, 1.73263882e-01, 4.66280378e+00, 1.61265060e-03}, {1.87469061e-02, 6.36079230e+00, 2.21992482e-01, 2.01624958e-01, 2.53027803e+01, 1.53610568e-02}};
					case 39:
						return {{6.75302747e-01, 4.70286720e-01, 2.63497677e+00, 1.09621746e-01, 9.60348773e-01, 5.28921560e-03}, {6.54331847e-02, 1.06108709e+02, 2.06643540e+00, 1.93131925e-01, 1.63310938e+00, 1.66083821e-02}};
					case 40:
						return {{2.64365505e+00, 5.54225147e-01, 7.61376625e-01, 6.02946890e-03, 9.91630530e-02, 9.56782020e-01}, {2.20202699e+00, 1.78260107e+02, 7.67218745e-02, 1.55143296e-02, 1.76175995e-01, 1.54330682e+00}};
					case 41:
						return {{6.59532875e-01, 1.84545854e+00, 1.25584405e+00, 1.22253422e-01, 7.06638328e-01, 2.62381590e-03}, {8.66145490e-02, 5.94774398e+00, 6.40851475e-01, 1.66646050e-01, 1.62853268e+00, 8.26257860e-03}};
					case 42:
						return {{6.10160120e-01, 1.26544000e+00, 1.97428762e+00, 6.48028962e-01, 2.60380820e-03, 1.13887493e-01}, {9.11628054e-02, 5.06776025e-01, 5.89590381e+00, 1.46634108e+00, 7.84336310e-03, 1.55114340e-01}};
					case 43:
						return {{8.55189183e-01, 1.66219641e+00, 1.45575475e+00, 1.05445664e-01, 7.71657112e-01, 2.20992640e-03}, {1.02962151e-01, 7.64907000e+00, 1.01639987e+00, 1.42303338e-01, 1.34659349e+00, 7.90358980e-03}};
					case 44:
						return {{4.70847093e-01, 1.58180781e+00, 2.02419818e+00, 1.97036260e-03, 6.26912639e-01, 1.02641320e-01}, {9.33029874e-02, 4.52831347e-01, 7.11489023e+00, 7.56181600e-03, 1.25399858e+00, 1.33786087e-01}};
					case 45:
						return {{4.20051553e-01, 1.76266507e+00, 2.02735641e+00, 1.45487180e-03, 6.22809600e-01, 9.91529915e-02}, {9.38882628e-02, 4.64441687e-01, 8.19346046e+00, 7.82704520e-03, 1.17194153e+00, 1.24532839e-01}};
					case 46:
						return {{2.10475155e+00, 2.03884487e+00, 1.82067264e-01, 9.52040948e-02, 5.91445248e-01, 1.13328680e-03}, {8.68606470e+00, 3.78924449e-01, 1.42921634e-01, 1.17125900e-01, 1.07843808e+00, 7.80252090e-03}};
					case 47:
						return {{2.07981390e+00, 4.43170726e-01, 1.96515215e+00, 5.96130591e-01, 4.78016333e-01, 9.46458470e-02}, {9.92540297e+00, 1.04920104e-01, 6.40103839e-01, 8.89594790e-01, 1.98509407e+00, 1.12744464e-01}};
					case 48:
						return {{1.63657549e+00, 2.17927989e+00, 7.71300690e-01, 6.64193880e-01, 7.64563285e-01, 8.61126689e-02}, {1.24540381e+01, 1.45134660e+00, 1.26695757e-01, 7.77659202e-01, 1.66075210e+00, 1.05728357e-01}};
					case 49:
						return {{2.24820632e+00, 1.64706864e+00, 7.88679265e-01, 8.12579069e-02, 6.68280346e-01, 6.38467475e-01}, {1.51913507e+00, 1.30113424e+01, 1.06128184e-01, 9.94045620e-02, 1.49742063e+00, 7.18422635e-01}};
					case 50:
						return {{2.16644620e+00, 6.88691021e-01, 1.92431751e+00, 5.65359888e-01, 9.18683861e-01, 7.80542213e-02}, {1.13174909e+01, 1.10131285e-01, 6.74464853e-01, 7.33564610e-01, 1.02310312e+01, 9.31104308e-02}};
					case 51:
						return {{1.73662114e+00, 9.99871380e-01, 2.13972409e+00, 5.60566526e-01, 9.93772747e-01, 7.37374982e-02}, {8.84334719e-01, 1.38462121e-01, 1.19666432e+01, 6.72672880e-01, 8.72330411e+00, 8.78577715e-02}};
					case 52:
						return {{2.09383882e+00, 1.56940519e+00, 1.30941993e+00, 6.98067804e-02, 1.04969537e+00, 5.55594354e-01}, {1.26856869e+01, 1.21236537e+00, 1.66633292e-01, 8.30817576e-02, 7.43147857e+00, 6.17487676e-01}};
					case 53:
						return {{1.60186925e+00, 1.98510264e+00, 1.48226200e+00, 5.53807199e-01, 1.11728722e+00, 6.60720847e-02}, {1.95031538e-01, 1.36976183e+01, 1.80304795e+00, 5.67912340e-01, 6.40879878e+00, 7.86615429e-02}};
					case 54:
						return {{1.60015487e+00, 1.71644581e+00, 1.84968351e+00, 6.23813648e-02, 1.21387555e+00, 5.54051946e-01}, {2.92913354e+00, 1.55882990e+01, 2.22525983e-01, 7.45581223e-02, 5.56013271e+00, 5.21994521e-01}};
					case 55:
						return {{2.95236854e+00, 4.28105721e-01, 1.89599233e+00, 5.48012938e-02, 4.70838600e+00, 5.90356719e-01}, {6.01461952e+00, 4.64151246e+01, 1.80109756e-01, 7.12799633e-02, 4.56702799e+01, 4.70236310e-01}};
					case 56:
						return {{3.19434243e+00, 1.98289586e+00, 1.55121052e-01, 6.73222354e-02, 4.48474211e+00, 5.42674414e-01}, {9.27352241e+00, 2.28741632e-01, 3.82000231e-02, 7.30961745e-02, 2.95703565e+01, 4.08647015e-01}};
					case 57:
						return {{2.05036425e+00, 1.42114311e-01, 3.23538151e+00, 6.34683429e-02, 3.97960586e+00, 5.20116711e-01}, {2.20348417e-01, 3.96438056e-02, 9.56979169e+00, 6.92443091e-02, 2.53178406e+01, 3.83614098e-01}};
					case 58:
						return {{3.22990759e+00, 1.57618307e-01, 2.13477838e+00, 5.01907609e-01, 3.80889010e+00, 5.96625028e-02}, {9.94660135e+00, 4.15378676e-02, 2.40480572e-01, 3.66252019e-01, 2.43275968e+01, 6.59653503e-02}};
					case 59:
						return {{1.58189324e-01, 3.18141995e+00, 2.27622140e+00, 3.97705472e+00, 5.58448277e-02, 4.85207954e-01}, {3.91309056e-02, 1.04139545e+01, 2.81671757e-01, 2.61872978e+01, 6.30921695e-02, 3.54234369e-01}};
					case 60:
						return {{1.81379417e-01, 3.17616396e+00, 2.35221519e+00, 3.83125763e+00, 5.25889976e-02, 4.70090742e-01}, {4.37324793e-02, 1.07842572e+01, 3.05571833e-01, 2.54745408e+01, 6.02676073e-02, 3.39017003e-01}};
					case 61:
						return {{1.92986811e-01, 2.43756023e+00, 3.17248504e+00, 3.58105414e+00, 4.56529394e-01, 4.94812177e-02}, {4.37785970e-02, 3.29336996e-01, 1.11259996e+01, 2.46709586e+01, 3.24990282e-01, 5.76553100e-02}};
					case 62:
						return {{2.12002595e-01, 3.16891754e+00, 2.51503494e+00, 4.44080845e-01, 3.36742101e+00, 4.65652543e-02}, {4.57703608e-02, 1.14536599e+01, 3.55561054e-01, 3.11953363e-01, 2.40291435e+01, 5.52266819e-02}};
					case 63:
						return {{2.59355002e+00, 3.16557522e+00, 2.29402652e-01, 4.32257780e-01, 3.17261920e+00, 4.37958317e-02}, {3.82452612e-01, 1.17675155e+01, 4.76642249e-02, 2.99719833e-01, 2.34462738e+01, 5.29440680e-02}};
					case 64:
						return {{3.19144939e+00, 2.55766431e+00, 3.32681934e-01, 4.14243130e-02, 2.61036728e+00, 4.20526863e-01}, {1.20224655e+01, 4.08338876e-01, 5.85819814e-02, 5.06771477e-02, 1.99344244e+01, 2.85686240e-01}};
					case 65:
						return {{2.59407462e-01, 3.16177855e+00, 2.75095751e+00, 2.79247686e+00, 3.85931001e-02, 4.10881708e-01}, {5.04689354e-02, 1.23140183e+01, 4.38337626e-01, 2.23797309e+01, 4.87920992e-02, 2.77622892e-01}};
					case 66:
						return {{3.16055396e+00, 2.82751709e+00, 2.75140255e-01, 4.00967160e-01, 2.63110834e+00, 3.61333817e-02}, {1.25470414e+01, 4.67899094e-01, 5.23226982e-02, 2.67614884e-01, 2.19498166e+01, 4.68871497e-02}};
					case 67:
						return {{2.88642467e-01, 2.90567296e+00, 3.15960159e+00, 3.91280259e-01, 2.48596038e+00, 3.37664478e-02}, {5.40507687e-02, 4.97581077e-01, 1.27599505e+01, 2.58151831e-01, 2.15400972e+01, 4.50664323e-02}};
					case 68:
						return {{3.15573213e+00, 3.11519560e-01, 2.97722406e+00, 3.81563854e-01, 2.40247532e+00, 3.15224214e-02}, {1.29729009e+01, 5.81399387e-02, 5.31213394e-01, 2.49195776e-01, 2.13627616e+01, 4.33253257e-02}};
					case 69:
						return {{3.15591970e+00, 3.22544710e-01, 3.05569053e+00, 2.92845100e-02, 3.72487205e-01, 2.27833695e+00}, {1.31232407e+01, 5.97223323e-02, 5.61876773e-01, 4.16534255e-02, 2.40821967e-01, 2.10034185e+01}};
					case 70:
						return {{3.10794704e+00, 3.14091221e+00, 3.75660454e-01, 3.61901097e-01, 2.45409082e+00, 2.72383990e-02}, {6.06347847e-01, 1.33705269e+01, 7.29814740e-02, 2.32652051e-01, 2.12695209e+01, 3.99969597e-02}};
					case 71:
						return {{3.11446863e+00, 5.39634353e-01, 3.06460915e+00, 2.58563745e-02, 2.13983556e+00, 3.47788231e-01}, {1.38968881e+01, 8.91708508e-02, 6.79919563e-01, 3.82808522e-02, 1.80078788e+01, 2.22706591e-01}};
					case 72:
						return {{3.01166899e+00, 3.16284788e+00, 6.33421771e-01, 3.41417198e-01, 1.53566013e+00, 2.40723773e-02}, {7.10401889e-01, 1.38262192e+01, 9.48486572e-02, 2.14129678e-01, 1.55298698e+01, 3.67833690e-02}};
					case 73:
						return {{3.20236821e+00, 8.30098413e-01, 2.86552297e+00, 2.24813887e-02, 1.40165263e+00, 3.33740596e-01}, {1.38446369e+01, 1.18381581e-01, 7.66369118e-01, 3.52934622e-02, 1.46148877e+01, 2.05704486e-01}};
					case 74:
						return {{9.24906855e-01, 2.75554557e+00, 3.30440060e+00, 3.29973862e-01, 1.09916444e+00, 2.06498883e-02}, {1.28663377e-01, 7.65826479e-01, 1.34471170e+01, 1.98218895e-01, 1.35087534e+01, 3.38918459e-02}};
					case 75:
						return {{1.96952105e+00, 1.21726619e+00, 4.10391685e+00, 2.90791978e-02, 2.30696669e-01, 6.08840299e-01}, {4.98830620e+01, 1.33243809e-01, 1.84396916e+00, 2.84192813e-02, 1.90968784e-01, 1.37090356e+00}};
					case 76:
						return {{2.06385867e+00, 1.29603406e+00, 3.96920673e+00, 2.69835487e-02, 2.31083999e-01, 6.30466774e-01}, {4.05671697e+01, 1.46559047e-01, 1.82561596e+00, 2.84172045e-02, 1.79765184e-01, 1.38911543e+00}};
					case 77:
						return {{2.21522726e+00, 1.37573155e+00, 3.78244405e+00, 2.44643240e-02, 2.36932016e-01, 6.48471412e-01}, {3.24464090e+01, 1.60920048e-01, 1.78756553e+00, 2.82909938e-02, 1.70692368e-01, 1.37928390e+00}};
					case 78:
						return {{9.84697940e-01, 2.73987079e+00, 3.61696715e+00, 3.02885602e-01, 2.78370726e-01, 1.52124129e-02}, {1.60910839e-01, 7.18971667e-01, 1.29281016e+01, 1.70134854e-01, 1.49862703e+00, 2.83510822e-02}};
					case 79:
						return {{9.61263398e-01, 3.69581030e+00, 2.77567491e+00, 2.95414176e-01, 3.11475743e-01, 1.43237267e-02}, {1.70932277e-01, 1.29335319e+01, 6.89997070e-01, 1.63525510e-01, 1.39200901e+00, 2.71265337e-02}};
					case 80:
						return {{1.29200491e+00, 2.75161478e+00, 3.49387949e+00, 2.77304636e-01, 4.30232810e-01, 1.48294351e-02}, {1.83432865e-01, 9.42368371e-01, 1.46235654e+01, 1.55110144e-01, 1.28871670e+00, 2.61903834e-02}};
					case 81:
						return {{3.75964730e+00, 3.21195904e+00, 6.47767825e-01, 2.76123274e-01, 3.18838810e-01, 1.31668419e-02}, {1.35041513e+01, 6.66330993e-01, 9.22518234e-02, 1.50312897e-01, 1.12565588e+00, 2.48879842e-02}};
					case 82:
						return {{1.00795975e+00, 3.09796153e+00, 3.61296864e+00, 2.62401476e-01, 4.05621995e-01, 1.31812509e-02}, {1.17268427e-01, 8.80453235e-01, 1.47325812e+01, 1.43491014e-01, 1.04103506e+00, 2.39575415e-02}};
					case 83:
						return {{1.59826875e+00, 4.38233925e+00, 2.06074719e+00, 1.94426023e-01, 8.22704978e-01, 2.33226953e-02}, {1.56897471e-01, 2.47094692e+00, 5.72438972e+01, 1.32979109e-01, 9.56532528e-01, 2.23038435e-02}};
					case 84:
						return {{1.71463223e+00, 2.14115960e+00, 4.37512413e+00, 2.16216680e-02, 1.97843837e-01, 6.52047920e-01}, {9.79262841e+01, 2.10193717e-01, 3.66948812e+00, 1.98456144e-02, 1.33758807e-01, 7.80432104e-01}};
					case 85:
						return {{1.48047794e+00, 2.09174630e+00, 4.75246033e+00, 1.85643958e-02, 2.05859375e-01, 7.13540948e-01}, {1.25943919e+02, 1.83803008e-01, 4.19890596e+00, 1.81383503e-02, 1.33035404e-01, 7.03031938e-01}};
					case 86:
						return {{6.30022295e-01, 3.80962881e+00, 3.89756067e+00, 2.40755100e-01, 2.62868577e+00, 3.14285931e-02}, {1.40909762e-01, 3.08515540e+01, 6.51559763e-01, 1.08899672e-01, 6.42383261e+00, 2.42346699e-02}};
					case 87:
						return {{5.23288135e+00, 2.48604205e+00, 3.23431354e-01, 2.55403596e-01, 5.53607228e-01, 5.75278890e-03}, {8.60599536e+00, 3.04543982e-01, 3.87759096e-02, 1.28717724e-01, 5.36977452e-01, 1.29417790e-02}};
					case 88:
						return {{1.44192685e+00, 3.55291725e+00, 3.91259586e+00, 2.16173519e-01, 3.94191605e+00, 4.60422605e-02}, {1.18740873e-01, 1.01739750e+00, 6.31814783e+01, 9.55806441e-02, 3.50602732e+01, 2.20850385e-02}};
					case 89:
						return {{1.45864127e+00, 4.18945405e+00, 3.65866182e+00, 2.08479229e-01, 3.16528117e+00, 5.23892556e-02}, {1.07760494e-01, 8.89090649e+01, 1.05088931e+00, 9.09335557e-02, 3.13297788e+01, 2.08807697e-02}};
					case 90:
						return {{1.19014064e+00, 2.55380607e+00, 4.68110181e+00, 2.26121303e-01, 3.58250545e-01, 7.82263950e-03}, {7.73468729e-02, 6.59693681e-01, 1.28013896e+01, 1.08632194e-01, 4.56765664e-01, 1.62623474e-02}};
					case 91:
						return {{4.68537504e+00, 2.98413708e+00, 8.91988061e-01, 2.24825384e-01, 3.04444846e-01, 9.48162710e-03}, {1.44503632e+01, 5.56438592e-01, 6.69512914e-02, 1.03235396e-01, 4.27255647e-01, 1.77730611e-02}};
					case 92:
						return {{4.63343606e+00, 3.18157056e+00, 8.76455075e-01, 2.21685477e-01, 2.72917100e-01, 1.11737298e-02}, {1.63377267e+01, 5.69517868e-01, 6.88860012e-02, 9.84254550e-02, 4.09470917e-01, 1.86215410e-02}};
					case 93:
						return {{4.56773888e+00, 3.40325179e+00, 8.61841923e-01, 2.19728870e-01, 2.38176903e-01, 1.38306499e-02}, {1.90992795e+01, 5.90099634e-01, 7.03204851e-02, 9.36334280e-02, 3.93554882e-01, 1.94437286e-02}};
					case 94:
						return {{5.45671123e+00, 1.11687906e-01, 3.30260343e+00, 1.84568319e-01, 4.93644263e-01, 3.57484743e+00}, {1.01892720e+01, 3.98131313e-02, 3.14622212e-01, 1.04220860e-01, 4.63080540e-01, 2.19369542e+01}};
					case 95:
						return {{5.38321999e+00, 1.23343236e-01, 3.46469090e+00, 1.75437132e-01, 3.39800073e+00, 4.69459519e-01}, {1.07289857e+01, 4.15137806e-02, 3.39326208e-01, 9.98932346e-02, 2.11601535e+01, 4.51996970e-01}};
					case 96:
						return {{5.38402377e+00, 3.49861264e+00, 1.88039547e-01, 1.69143137e-01, 3.19595016e+00, 4.64393059e-01}, {1.11211419e+01, 3.56750210e-01, 5.39853583e-02, 9.60082633e-02, 1.80694389e+01, 4.36318197e-01}};
					case 97:
						return {{3.66090688e+00, 2.03054678e-01, 5.30697515e+00, 1.60934046e-01, 3.04808401e+00, 4.43610295e-01}, {3.84420906e-01, 5.48547131e-02, 1.17150262e+01, 9.21020329e-02, 1.73525367e+01, 4.27132359e-01}};
					case 98:
						return {{3.94150390e+00, 5.16915345e+00, 1.61941074e-01, 4.15299561e-01, 2.91761325e+00, 1.51474927e-01}, {4.18246722e-01, 1.25201788e+01, 4.81540117e-02, 4.24913856e-01, 1.90899693e+01, 8.81568925e-02}};
					case 99:
						return {{4.09780623e+00, 5.10079393e+00, 1.74617289e-01, 2.76774658e+00, 1.44496639e-01, 4.02772109e-01}, {4.46021145e-01, 1.31768613e+01, 5.02742829e-02, 1.84815393e+01, 8.46232592e-02, 4.17640100e-01}};
					case 100:
						return {{4.24934820e+00, 5.03556594e+00, 1.88920613e-01, 3.94356058e-01, 2.61213100e+00, 1.38001927e-01}, {4.75263933e-01, 1.38570834e+01, 5.26975158e-02, 4.11193751e-01, 1.78537905e+01, 8.12774434e-02}};
					case 101:
						return {{2.00942931e-01, 4.40119869e+00, 4.97250102e+00, 2.47530599e+00, 3.86883197e-01, 1.31936095e-01}, {5.48366518e-02, 5.04248434e-01, 1.45721366e+01, 1.72978308e+01, 4.05043898e-01, 7.80821071e-02}};
					case 102:
						return {{2.16052899e-01, 4.91106799e+00, 4.54862870e+00, 2.36114249e+00, 1.26277292e-01, 3.81364501e-01}, {5.83584058e-02, 1.53264212e+01, 5.34434760e-01, 1.68164803e+01, 7.50304633e-02, 3.99305852e-01}};
					case 103:
						return {{4.86738014e+00, 3.19974401e-01, 4.58872425e+00, 1.21482448e-01, 2.31639872e+00, 3.79258137e-01}, {1.60320520e+01, 6.70871138e-02, 5.77039373e-01, 7.22275899e-02, 1.41279737e+01, 3.89973484e-01}};
				}

				return {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
			}

			// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
			LNL_Coef_cpu<T> load_feg_weickenmeier_neutral_0_12(const dt_int32& Z)
			{
				switch(Z)
				{
					case 2:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {2.542e+00, 8.743e+00, 1.269e+01, 4.371e-01, 5.294e+00, 2.825e+01}};
					case 3:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {6.845e-01, 3.065e+00, 6.240e+00, 1.262e+02, 1.312e+02, 1.318e+02}};
					case 4:
						return {{3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01}, {5.400e-01, 3.388e+00, 5.562e+01, 5.078e+01, 6.701e+01, 9.637e+01}};
					case 5:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {3.314e-01, 2.975e+00, 3.401e+01, 3.598e+01, 3.668e+01, 6.081e+01}};
					case 6:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {2.946e-01, 3.934e+00, 2.498e+01, 2.528e+01, 2.547e+01, 4.670e+01}};
					case 7:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {2.393e-01, 4.935e+00, 1.812e+01, 1.570e+01, 1.582e+01, 4.024e+01}};
					case 8:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {6.376e+00, 8.037e+00, 2.721e+01, 1.116e-01, 3.869e-01, 1.090e+01}};
					case 9:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {2.180e-01, 6.770e+00, 7.051e+00, 6.675e+00, 1.238e+01, 2.808e+01}};
					case 10:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {2.006e-01, 5.498e+00, 6.281e+00, 7.192e+00, 7.548e+00, 2.326e+01}};
					case 11:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {2.190e-01, 5.300e+00, 5.319e+00, 5.283e+00, 5.285e+00, 1.282e+02}};
					case 12:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {1.976e+00, 2.809e+00, 1.639e+01, 5.490e-02, 2.061e+00, 1.217e+02}};
					case 13:
						return {{4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01}, {2.297e+00, 2.358e+00, 2.499e+01, 7.460e-02, 5.595e-01, 1.285e+02}};
					case 14:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {1.737e+00, 3.043e+00, 3.057e+01, 5.070e-02, 9.918e-01, 8.618e+01}};
					case 15:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {1.795e-01, 2.632e+00, 2.676e+00, 3.457e+01, 3.678e+01, 5.406e+01}};
					case 16:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {1.006e+00, 4.904e+00, 3.135e+01, 3.700e-02, 9.870e-01, 4.494e+01}};
					case 17:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {1.846e-01, 1.480e+00, 5.210e+00, 2.479e+01, 3.206e+01, 3.910e+01}};
					case 18:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {2.006e-01, 6.533e+00, 2.272e+01, 1.200e+00, 1.274e+00, 3.626e+01}};
					case 19:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {4.442e-01, 3.367e+00, 1.963e+01, 1.820e-02, 2.351e+01, 2.129e+02}};
					case 20:
						return {{3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01}, {1.827e-01, 2.066e+00, 1.699e+01, 1.158e+01, 1.398e+01, 1.861e+02}};
					case 21:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {1.425e-01, 1.466e+00, 1.547e+01, 4.243e+00, 9.804e+00, 1.215e+02}};
					case 22:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {1.278e-01, 1.456e+00, 1.210e+01, 4.617e+00, 1.197e+01, 1.050e+02}};
					case 23:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {1.313e-01, 1.399e+00, 8.008e+00, 7.981e+00, 1.341e+01, 9.531e+01}};
					case 24:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {1.231e-01, 2.384e+00, 9.921e+00, 1.648e+00, 1.100e+01, 6.846e+01}};
					case 25:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {4.817e-01, 3.783e+00, 8.473e+00, 4.690e-02, 8.745e+00, 7.744e+01}};
					case 26:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {4.470e-01, 6.894e+00, 6.903e+00, 5.690e-02, 3.026e+00, 7.087e+01}};
					case 27:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {1.071e-01, 3.636e+00, 7.558e+00, 1.280e+00, 5.140e+00, 6.716e+01}};
					case 28:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {1.107e-01, 1.619e+00, 6.003e+00, 5.975e+00, 6.060e+00, 5.941e+01}};
					case 29:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {1.129e-01, 1.891e+00, 5.085e+00, 5.073e+00, 5.099e+00, 4.639e+01}};
					case 30:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {1.021e-01, 1.734e+00, 4.783e+00, 4.807e+00, 5.645e+00, 5.122e+01}};
					case 31:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {1.064e-01, 1.537e+00, 5.138e+00, 4.743e+00, 5.000e+00, 6.143e+01}};
					case 32:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {9.580e-02, 1.677e+00, 4.703e+00, 2.912e+00, 7.870e+00, 6.494e+01}};
					case 33:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {9.430e-02, 2.214e+00, 3.951e+00, 1.521e+00, 1.581e+01, 5.241e+01}};
					case 34:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {9.250e-02, 1.602e+00, 3.049e+00, 3.185e+00, 1.894e+01, 4.763e+01}};
					case 35:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {9.250e-02, 1.773e+00, 3.481e+00, 1.884e+00, 2.269e+01, 4.069e+01}};
					case 36:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {4.932e-01, 2.083e+00, 1.141e+01, 3.330e-02, 2.097e+00, 4.238e+01}};
					case 37:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {1.580e-01, 1.715e+00, 9.392e+00, 1.675e+00, 2.359e+01, 1.525e+02}};
					case 38:
						return {{3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01}, {3.605e-01, 2.128e+00, 1.246e+01, 1.530e-02, 2.108e+00, 1.332e+02}};
					case 39:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {9.000e-02, 1.414e+00, 2.053e+00, 1.026e+01, 1.075e+01, 9.064e+01}};
					case 40:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {1.009e-01, 1.154e+00, 2.347e+00, 1.058e+01, 1.095e+01, 8.282e+01}};
					case 41:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {9.240e-02, 1.170e+00, 5.940e+00, 1.306e+00, 1.343e+01, 6.637e+01}};
					case 42:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {4.354e-01, 1.248e+00, 7.454e+00, 3.540e-02, 9.914e+00, 6.172e+01}};
					case 43:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {4.594e-01, 1.182e+00, 8.317e+00, 3.230e-02, 8.323e+00, 6.498e+01}};
					case 44:
						return {{4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01}, {8.600e-02, 1.396e+00, 1.170e+01, 1.396e+00, 3.452e+00, 5.556e+01}};
					case 45:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {9.210e-02, 1.113e+00, 7.658e+00, 1.126e+00, 8.325e+00, 4.838e+01}};
					case 46:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {9.010e-02, 1.125e+00, 9.698e+00, 1.085e+00, 5.709e+00, 3.349e+01}};
					case 47:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {8.940e-02, 3.191e+00, 9.100e+00, 8.090e-01, 8.144e-01, 4.134e+01}};
					case 48:
						return {{3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01}, {2.885e-01, 1.613e+00, 8.997e+00, 1.710e-02, 9.467e+00, 5.813e+01}};
					case 49:
						return {{4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01}, {8.950e-02, 1.233e+00, 8.231e+00, 1.224e+00, 7.062e+00, 5.970e+01}};
					case 50:
						return {{6.000e-01, 6.000e-01, 6.000e-01, 6.000e-01, 6.000e-01, 6.000e-01}, {7.120e-02, 8.553e-01, 6.401e+00, 1.336e+00, 6.382e+00, 5.092e+01}};
					case 51:
						return {{6.000e-01, 6.000e-01, 6.000e-01, 6.000e-01, 6.000e-01, 6.000e-01}, {3.575e-01, 1.325e+00, 6.517e+00, 3.550e-02, 6.519e+00, 5.081e+01}};
					case 52:
						return {{6.000e-01, 6.000e-01, 6.000e-01, 6.000e-01, 6.000e-01, 6.000e-01}, {5.009e-01, 3.953e+00, 7.628e+00, 3.010e-02, 5.074e-01, 4.963e+01}};
					case 53:
						return {{4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01}, {8.430e-02, 1.130e+00, 8.862e+00, 1.130e+00, 9.132e+00, 5.602e+01}};
					case 54:
						return {{4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01}, {2.780e-01, 1.621e+00, 1.145e+01, 2.030e-02, 3.275e+00, 5.144e+01}};
					case 55:
						return {{1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01}, {1.204e-01, 1.537e+00, 9.816e+00, 4.122e+01, 4.262e+01, 2.243e+02}};
					case 56:
						return {{1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01}, {1.223e-01, 1.449e+00, 9.502e+00, 4.941e+01, 7.495e+01, 2.170e+02}};
					case 57:
						return {{3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01}, {8.930e-02, 1.262e+00, 8.097e+00, 1.203e+00, 1.766e+01, 1.166e+02}};
					case 58:
						return {{3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01}, {8.500e-02, 1.283e+00, 1.122e+01, 1.327e+00, 4.610e+00, 1.122e+02}};
					case 59:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {9.810e-02, 1.526e+00, 8.590e+00, 1.239e+00, 2.249e+01, 1.400e+02}};
					case 60:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {9.410e-02, 1.266e+00, 5.988e+00, 1.779e+01, 1.814e+01, 1.326e+02}};
					case 61:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {9.450e-02, 1.251e+00, 5.912e+00, 1.629e+01, 1.673e+01, 1.279e+02}};
					case 62:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {9.060e-02, 1.593e+00, 1.064e+01, 1.789e+00, 2.221e+00, 1.246e+02}};
					case 63:
						return {{1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01}, {1.049e-01, 1.544e+00, 8.652e+00, 7.093e+00, 5.337e+01, 1.837e+02}};
					case 64:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {9.340e-02, 1.387e+00, 7.359e+00, 1.551e+00, 2.082e+01, 1.110e+02}};
					case 65:
						return {{1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01}, {1.019e-01, 1.524e+00, 7.169e+00, 2.086e+01, 4.929e+01, 1.661e+02}};
					case 66:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {8.400e-02, 1.409e+00, 7.140e+00, 1.348e+00, 1.142e+01, 1.080e+02}};
					case 67:
						return {{1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01}, {9.440e-02, 1.618e+00, 6.271e+00, 4.035e+01, 4.283e+01, 1.306e+02}};
					case 68:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {8.210e-02, 1.251e+00, 4.812e+00, 1.084e+01, 1.090e+01, 1.001e+02}};
					case 69:
						return {{1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01}, {9.660e-02, 1.602e+00, 5.675e+00, 3.059e+01, 3.113e+01, 1.387e+02}};
					case 70:
						return {{1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01}, {9.490e-02, 1.602e+00, 5.439e+00, 2.831e+01, 2.928e+01, 1.381e+02}};
					case 71:
						return {{1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01}, {9.660e-02, 1.568e+00, 5.322e+00, 3.418e+01, 3.525e+01, 1.214e+02}};
					case 72:
						return {{1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01}, {9.290e-02, 1.555e+00, 5.251e+00, 3.752e+01, 3.888e+01, 1.052e+02}};
					case 73:
						return {{4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01}, {6.300e-02, 8.195e-01, 2.891e+00, 5.543e+00, 5.981e+00, 5.442e+01}};
					case 74:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {7.900e-02, 1.371e+00, 8.234e+00, 1.383e+00, 1.392e+00, 7.712e+01}};
					case 75:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {5.270e-02, 9.072e-01, 4.438e+00, 9.459e-01, 4.375e+00, 4.398e+01}};
					case 76:
						return {{4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01}, {2.270e-01, 1.570e+00, 6.345e+00, 1.560e-02, 1.618e+00, 4.616e+01}};
					case 77:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {5.060e-02, 8.677e-01, 5.093e+00, 8.812e-01, 3.569e+00, 3.977e+01}};
					case 78:
						return {{5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01, 5.000e-01}, {5.250e-02, 8.377e-01, 3.959e+00, 8.152e-01, 6.442e+00, 3.421e+01}};
					case 79:
						return {{4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01}, {5.493e-01, 1.728e+00, 6.720e+00, 2.640e-02, 7.250e-02, 3.546e+01}};
					case 80:
						return {{4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01}, {2.194e-01, 1.416e+00, 6.682e+00, 1.470e-02, 1.576e+00, 3.716e+01}};
					case 81:
						return {{4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01}, {2.246e-01, 1.128e+00, 4.303e+00, 1.490e-02, 7.156e+00, 4.309e+01}};
					case 82:
						return {{3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01}, {6.430e-02, 1.194e+00, 7.393e+00, 1.142e+00, 1.289e+00, 5.113e+01}};
					case 83:
						return {{4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01}, {5.380e-02, 8.672e-01, 1.875e+00, 7.648e+00, 7.868e+00, 4.564e+01}};
					case 84:
						return {{4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01}, {5.011e-01, 1.638e+00, 6.786e+00, 2.190e-02, 8.600e-02, 4.673e+01}};
					case 85:
						return {{4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01}, {2.232e-01, 1.108e+00, 3.591e+00, 1.010e-02, 1.164e+01, 4.507e+01}};
					case 86:
						return {{4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01, 4.000e-01}, {2.115e-01, 1.140e+00, 3.415e+00, 1.190e-02, 1.341e+01, 4.311e+01}};
					case 87:
						return {{1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01, 1.000e-01}, {9.440e-02, 1.026e+00, 6.255e+00, 3.251e+01, 3.629e+01, 1.491e+02}};
					case 88:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {7.300e-02, 1.018e+00, 5.896e+00, 1.031e+00, 2.037e+01, 1.153e+02}};
					case 89:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {7.520e-02, 9.494e-01, 3.725e+00, 1.758e+01, 1.975e+01, 1.091e+02}};
					case 90:
						return {{3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01}, {6.390e-02, 9.019e-01, 4.657e+00, 9.025e-01, 1.571e+01, 8.370e+01}};
					case 91:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {7.560e-02, 8.492e-01, 4.010e+00, 1.695e+01, 1.779e+01, 1.002e+02}};
					case 92:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {7.140e-02, 1.149e+00, 9.212e+00, 9.592e-01, 1.203e+00, 1.043e+02}};
					case 93:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {6.920e-02, 9.810e-01, 5.954e+00, 9.909e-01, 2.206e+01, 9.098e+01}};
					case 94:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {7.140e-02, 9.577e-01, 6.132e+00, 9.744e-01, 1.567e+01, 8.987e+01}};
					case 95:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {7.300e-02, 9.327e-01, 6.348e+00, 9.103e-01, 1.326e+01, 8.686e+01}};
					case 96:
						return {{3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01}, {5.780e-02, 7.227e-01, 3.011e+00, 9.219e+00, 9.534e+00, 6.587e+01}};
					case 97:
						return {{2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01, 2.000e-01}, {7.090e-02, 7.759e-01, 6.143e+00, 1.790e+00, 1.512e+01, 8.357e+01}};
					case 98:
						return {{3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01, 3.000e-01}, {6.160e-02, 8.136e-01, 6.562e+00, 8.381e-01, 4.189e+00, 6.141e+01}};
				}

				return {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
			}

			// 6: Lobato parameterization - Hydrogen functions - [0, 12]
			LNL_Coef_cpu<double> load_feg_lobato_neutral_0_12(const dt_int32& Z)
			{
				switch(Z)
				{
					case 1:
						return {{5.29176831e-02, 5.29176831e-02, 5.29176831e-02, 5.29176831e-02, 5.29176831e-02}, {2.76376891e+00, 2.76376891e+00, 2.76376891e+00, 2.76376891e+00, 2.76376891e+00}};
					case 2:
						return {{1.01746701e-01, 9.49040428e-02, -4.11936715e-02, 5.35294712e-02, 1.34533102e-05}, {1.49190903e+00, 1.01103199e+00, 4.90177989e-01, 4.72138286e-01, 2.12409496e-01}};
					case 3:
						return {{6.19087696e-01, 1.72226405e+00, -7.92180777e-01, 8.37330297e-02, 1.01449499e-02}, {1.05523100e+01, 9.23192024e+00, 4.76025486e+00, 5.19986629e-01, 2.14180499e-01}};
					case 4:
						return {{1.07960296e+00, 7.13181376e-01, -3.91357690e-01, 1.06852397e-01, 1.86710991e-02}, {6.00019121e+00, 4.10317707e+00, 1.37123895e+00, 5.52847624e-01, 1.53891206e-01}};
					case 5:
						return {{1.31485796e+00, 2.52792597e-01, -2.02927098e-01, 2.92299408e-02, 3.54602700e-03}, {3.94681811e+00, 1.91365397e+00, 9.58502173e-01, 1.62817106e-01, 7.77361393e-02}};
					case 6:
						return {{7.39281714e-01, 7.50964522e-01, -2.54863292e-01, 1.76229905e-02, 1.34442805e-03}, {3.47855496e+00, 1.63563204e+00, 8.30155313e-01, 9.63626578e-02, 5.00164293e-02}};
					case 7:
						return {{3.26442391e-01, 8.36000979e-01, -7.21398070e-02, 1.31275300e-02, 5.07714599e-03}, {3.17208099e+00, 1.59607995e+00, 3.68962586e-01, 9.75946933e-02, 4.86027598e-02}};
					case 8:
						return {{7.85341978e-01, 2.52642393e-01, -5.85432090e-02, 1.27868997e-02, 1.97321596e-03}, {1.78015697e+00, 8.17416310e-01, 2.44515404e-01, 6.56605288e-02, 3.27121913e-02}};
					case 9:
						return {{6.07230484e-01, 4.81327087e-01, -1.94196299e-01, 7.56342197e-03, 7.61223928e-05}, {1.58934402e+00, 5.77747822e-01, 3.46372694e-01, 3.73525992e-02, 1.65033191e-02}};
					case 10:
						return {{6.13327866e-05, 4.95682299e-01, 8.31292927e-01, -5.06411910e-01, 5.32641588e-03}, {1.25770597e+01, 1.38210797e+00, 4.62477595e-01, 3.64122301e-01, 2.77422108e-02}};
					case 11:
						return {{2.38567805e+00, -5.41312218e-01, 5.80007017e-01, -4.31193784e-02, 4.69793798e-03}, {1.17757797e+01, 5.18739605e+00, 5.75649977e-01, 1.69096500e-01, 2.32259594e-02}};
					case 12:
						return {{2.48242688e+00, -4.73447889e-01, 8.96034002e-01, -3.08043092e-01, 3.78086301e-03}, {7.48703098e+00, 1.11762202e+00, 3.95489693e-01, 2.52319604e-01, 1.92819294e-02}};
					case 13:
						return {{2.84625793e+00, -5.46510279e-01, 7.44476914e-01, -1.04492702e-01, 3.32048093e-03}, {6.67600107e+00, 7.45036721e-01, 3.77064496e-01, 1.67518005e-01, 1.65107101e-02}};
					case 14:
						return {{2.85535312e+00, -4.35813814e-01, 6.98928118e-01, -2.16126293e-01, 3.03548202e-03}, {5.05003214e+00, 5.29675722e-01, 2.47670904e-01, 1.50627494e-01, 1.44170700e-02}};
					case 15:
						return {{2.83505607e+00, -2.85362303e-01, 2.51688510e-01, -5.73274195e-02, 3.27281700e-03}, {3.81684399e+00, 1.01348305e+00, 1.71625093e-01, 7.74765834e-02, 1.31826401e-02}};
					case 16:
						return {{2.71400809e+00, -3.41581792e-01, 2.40847498e-01, -3.09837908e-02, 2.47210590e-03}, {2.95212793e+00, 6.79380417e-01, 1.75853401e-01, 6.52635694e-02, 1.12200603e-02}};
					case 17:
						return {{2.53558207e+00, -1.92675203e-01, 1.20515399e-01, -3.50671113e-02, 3.20018711e-03}, {2.74467707e+00, 2.44560409e+00, 9.37200710e-02, 4.38743904e-02, 1.08148400e-02}};
					case 18:
						return {{2.39861608e+00, -6.07183397e-01, 5.24528027e-01, -2.69268993e-02, 1.87984400e-03}, {2.05899811e+00, 3.09629112e-01, 1.89788997e-01, 5.85227199e-02, 8.71843286e-03}};
					case 19:
						return {{2.46684504e+00, 1.92999494e+00, 1.54350102e-01, -8.12437236e-02, 2.15101591e-03}, {2.41689091e+01, 2.15578604e+00, 6.77528828e-02, 4.73528206e-02, 8.29782058e-03}};
					case 20:
						return {{3.55370402e+00, 1.34050202e+00, 6.76007420e-02, -1.39076998e-02, 5.67091210e-03}, {1.34759798e+01, 1.55728805e+00, 6.80236369e-02, 1.71041396e-02, 9.30336956e-03}};
					case 21:
						return {{3.44377804e+00, 1.13595498e+00, 6.00283407e-02, 1.71150803e-03, -1.41223101e-03}, {1.11026001e+01, 1.23438895e+00, 8.91299099e-02, 2.11031688e-03, 2.00439896e-03}};
					case 22:
						return {{2.24129701e+00, 2.07236791e+00, 9.29673165e-02, -4.03135009e-02, 2.26346403e-03}, {1.70091000e+01, 1.85258400e+00, 4.71124090e-02, 2.76718698e-02, 6.66805590e-03}};
					case 23:
						return {{2.92555189e+00, 1.16348600e+00, 4.43377793e-02, -1.35539500e-02, 4.69834497e-03}, {9.08362007e+00, 1.00832999e+00, 4.14920002e-02, 1.35524096e-02, 7.13791978e-03}};
					case 24:
						return {{1.64590800e+00, 1.80008602e+00, 5.35616912e-02, -2.09835991e-02, 4.11470514e-03}, {1.21226702e+01, 1.31736004e+00, 3.45182195e-02, 1.51794404e-02, 6.56798482e-03}};
					case 25:
						return {{2.39974809e+00, 1.31345499e+00, 5.21507002e-02, -4.50833105e-02, 1.73789803e-02}, {8.79271221e+00, 9.28560615e-01, 2.16320492e-02, 1.14173004e-02, 7.73268519e-03}};
					case 26:
						return {{1.84187496e+00, 1.69585705e+00, 4.24552783e-02, 2.01590403e-04, -1.46712104e-04}, {1.19059095e+01, 1.14025199e+00, 5.56822792e-02, 4.73777589e-04, 4.29703505e-04}};
					case 27:
						return {{2.13464499e+00, 1.70659602e+00, -5.68710923e-01, 1.51593298e-01, 5.05861302e-04}, {8.62808323e+00, 5.93899727e-01, 2.49489203e-01, 9.49204788e-02, 3.46111390e-03}};
					case 28:
						return {{2.31999111e+00, 1.31185806e+00, -6.97869778e-01, 3.42176288e-01, 4.83954413e-04}, {5.89068222e+00, 4.14797008e-01, 1.55506194e-01, 9.88660678e-02, 3.21499701e-03}};
					case 29:
						return {{1.23546398e+00, 1.54618704e+00, 2.77369693e-02, -1.17329899e-02, 4.77816490e-03}, {1.07732496e+01, 8.05693388e-01, 2.23376006e-02, 7.92775396e-03, 4.88925213e-03}};
					case 30:
						return {{1.99225104e+00, 1.36421597e+00, -5.34584999e-01, 2.10798606e-01, 3.86487402e-04}, {5.84002686e+00, 4.15746897e-01, 1.55540794e-01, 8.27629790e-02, 2.77534500e-03}};
					case 31:
						return {{2.66632009e+00, 1.94698906e+00, -1.18554604e+00, 1.21890597e-01, 3.09190800e-04}, {6.34141302e+00, 2.78902292e-01, 1.83026701e-01, 6.34498373e-02, 2.52509094e-03}};
					case 32:
						return {{3.44351792e+00, -1.13572204e+00, 1.36060297e+00, 2.05168296e-02, 2.72685196e-04}, {4.83926678e+00, 1.60483098e+00, 5.87697089e-01, 3.27054188e-02, 2.32771295e-03}};
					case 33:
						return {{3.23248005e+00, -1.32028103e+00, 1.74037302e+00, 1.42622404e-02, 1.96056702e-04}, {4.26419592e+00, 9.00477171e-01, 5.46430409e-01, 2.42819097e-02, 2.06261105e-03}};
					case 34:
						return {{3.63524103e+00, -2.07056499e+00, 2.03258395e+00, 1.16819199e-02, 1.22890197e-04}, {3.28155088e+00, 8.77488911e-01, 5.28532326e-01, 1.97308101e-02, 1.78012101e-03}};
					case 35:
						return {{3.69151306e+00, -1.79441595e+00, 1.62838995e+00, 7.51289679e-03, 5.47963391e-05}, {2.62788510e+00, 7.07726896e-01, 4.17159408e-01, 1.39476703e-02, 1.39450806e-03}};
					case 36:
						return {{3.31220102e+00, -5.70945084e-01, 6.97833002e-01, 9.94181167e-03, 1.29992593e-04}, {2.61677909e+00, 9.20762420e-01, 3.22337389e-01, 1.78531203e-02, 1.64635398e-03}};
					case 37:
						return {{2.98498106e+00, 2.36873889e+00, 5.00381887e-01, 7.34424405e-03, 4.27639898e-05}, {2.40395908e+01, 2.55953598e+00, 2.63484299e-01, 1.31538603e-02, 1.21721299e-03}};
					case 38:
						return {{4.12212515e+00, 1.92479300e+00, 4.70577389e-01, 7.04240007e-03, 3.69832087e-05}, {1.62667599e+01, 2.20284200e+00, 2.45777294e-01, 1.25011001e-02, 1.12981803e-03}};
					case 39:
						return {{3.74598193e+00, 2.12386703e+00, 4.42405105e-01, 5.74155105e-03, 1.75479108e-05}, {1.34013996e+01, 2.35688090e+00, 2.22515702e-01, 1.05318800e-02, 8.94146215e-04}};
					case 40:
						return {{3.38130689e+00, 2.25415707e+00, 4.40226912e-01, 6.77093910e-03, 6.23684973e-05}, {1.21639204e+01, 2.35530710e+00, 2.21562207e-01, 1.22107305e-02, 1.17974204e-03}};
					case 41:
						return {{2.17703199e+00, 2.78299403e+00, 3.90569508e-01, 4.67265584e-03, 6.22945709e-06}, {1.11494303e+01, 2.41581297e+00, 1.92431599e-01, 8.68258812e-03, 6.24530192e-04}};
					case 42:
						return {{1.67860794e+00, 3.11303711e+00, 3.46704096e-01, 4.59088618e-03, 9.16880344e-06}, {1.36886101e+01, 2.27832890e+00, 1.75135896e-01, 8.52385256e-03, 6.70720416e-04}};
					case 43:
						return {{2.88291001e+00, 2.22659802e+00, 2.98569113e-01, 4.29125316e-03, 9.66334756e-06}, {8.43309307e+00, 1.69386697e+00, 1.55848294e-01, 8.09361599e-03, 6.57025608e-04}};
					case 44:
						return {{1.31261802e+00, 3.17361403e+00, 2.85882711e-01, 3.99975199e-03, 8.73038971e-06}, {1.25337400e+01, 1.96410704e+00, 1.46394894e-01, 7.61754299e-03, 6.17400685e-04}};
					case 45:
						return {{1.44103205e+00, 2.94216609e+00, 2.32040107e-01, 3.09590693e-03, 2.05679089e-06}, {1.03466597e+01, 1.64708400e+00, 1.23005100e-01, 6.30735699e-03, 3.96098098e-04}};
					case 46:
						return {{2.02130705e-01, 3.35615706e+00, 2.30525196e-01, 3.56838107e-03, 4.97298106e-06}, {1.69279804e+01, 1.69309103e+00, 1.22913197e-01, 6.80784881e-03, 4.96210472e-04}};
					case 47:
						return {{1.21013403e+00, 2.92162895e+00, 1.98055595e-01, 3.18102399e-03, 8.53895926e-06}, {1.01067600e+01, 1.43554401e+00, 1.08165003e-01, 6.40084315e-03, 5.55367384e-04}};
					case 48:
						return {{1.60705805e+00, 2.81508899e+00, 1.91956207e-01, 2.63145100e-03, 1.84189503e-06}, {9.29942894e+00, 1.37371194e+00, 1.02014497e-01, 5.39610581e-03, 3.47371009e-04}};
					case 49:
						return {{2.70126796e+00, 2.29345989e+00, 2.16367900e-01, 3.69811291e-03, 1.95252096e-06}, {8.33631229e+00, 1.19331098e+00, 1.16310202e-01, 6.37960806e-03, 3.49309004e-04}};
					case 50:
						return {{2.87292099e+00, 2.40417600e+00, 1.57168493e-01, 1.89157098e-03, 7.62896207e-06}, {8.46732521e+00, 1.12822700e+00, 8.41042623e-02, 4.36692312e-03, 4.75006207e-04}};
					case 51:
						return {{3.50419807e+00, 1.83627403e+00, 1.50732994e-01, 3.59381200e-03, 1.05334497e-04}, {5.89424086e+00, 9.22941923e-01, 8.61572027e-02, 8.35711230e-03, 8.86740920e-04}};
					case 52:
						return {{3.21798110e+00, 2.05408001e+00, 2.15190396e-01, 2.72924802e-03, 2.24804194e-06}, {5.26334810e+00, 1.19520295e+00, 1.02286197e-01, 5.08707995e-03, 3.28804512e-04}};
					case 53:
						return {{4.13265181e+00, 1.18767595e+00, 1.35713801e-01, 1.50258699e-03, 1.33380800e-05}, {3.71807408e+00, 7.28209019e-01, 7.14013129e-02, 3.71678802e-03, 5.00523078e-04}};
					case 54:
						return {{3.51788712e+00, 1.66314995e+00, 2.16246501e-01, 2.76901806e-03, 8.36265008e-06}, {3.94756007e+00, 1.07487500e+00, 1.00129701e-01, 5.00469515e-03, 4.41977201e-04}};
					case 55:
						return {{3.81201100e+00, 3.07259011e+00, 1.22549701e+00, 8.76324698e-02, 2.59851106e-04}, {2.28876896e+01, 3.35200000e+00, 5.80402970e-01, 4.92655188e-02, 8.89977906e-04}};
					case 56:
						return {{6.35506678e+00, 1.69075298e+00, 9.57290590e-01, 7.86023363e-02, 2.44261813e-04}, {1.39797802e+01, 1.68795502e+00, 5.12502193e-01, 4.48691696e-02, 8.56465078e-04}};
					case 57:
						return {{4.28002501e+00, 3.07966208e+00, 1.45083499e+00, 8.27694386e-02, 2.41575995e-04}, {1.66275997e+01, 4.31141520e+00, 6.13337815e-01, 4.52826284e-02, 8.27733311e-04}};
					case 58:
						return {{5.13602400e+00, 2.43740010e+00, 1.03740597e+00, 6.99751526e-02, 2.29697005e-04}, {1.42508097e+01, 2.43046999e+00, 4.76604193e-01, 4.05691713e-02, 7.97950372e-04}};
					case 59:
						return {{5.35353422e+00, 2.12939692e+00, 9.15245771e-01, 7.25957975e-02, 2.31387705e-04}, {1.41875095e+01, 1.81712794e+00, 4.47230697e-01, 4.14027907e-02, 7.73125212e-04}};
					case 60:
						return {{5.32854700e+00, 1.72848594e+00, 1.11297500e+00, 7.37465993e-02, 2.25470707e-04}, {1.31087599e+01, 1.74724901e+00, 4.96911108e-01, 4.08121310e-02, 7.47739687e-04}};
					case 61:
						return {{4.36711407e+00, 1.98772097e+00, 1.65877998e+00, 7.46220574e-02, 2.23465599e-04}, {1.53489799e+01, 3.45163202e+00, 5.88005126e-01, 4.03418690e-02, 7.24034617e-04}};
					case 62:
						return {{5.13235712e+00, 1.91035604e+00, 8.09289694e-01, 6.53573796e-02, 2.11307895e-04}, {1.30506601e+01, 1.27198803e+00, 4.12205100e-01, 3.68197896e-02, 6.99623721e-04}};
					case 63:
						return {{4.25057316e+00, 2.08642793e+00, 1.35467899e+00, 6.16801083e-02, 2.06875600e-04}, {1.51125002e+01, 2.45922089e+00, 4.77150202e-01, 3.51009294e-02, 6.77573611e-04}};
					case 64:
						return {{4.78696299e+00, 1.37942696e+00, 1.41249299e+00, 5.88963702e-02, 1.96196401e-04}, {1.08546200e+01, 2.08167791e+00, 4.78227913e-01, 3.32649015e-02, 6.55558892e-04}};
					case 65:
						return {{4.72020006e+00, 1.16532505e+00, 1.50215197e+00, 6.71210736e-02, 2.09402002e-04}, {1.22190304e+01, 1.67899096e+00, 5.07295787e-01, 3.62123288e-02, 6.38881815e-04}};
					case 66:
						return {{4.22459507e+00, 1.57796502e+00, 1.44352305e+00, 6.22941405e-02, 1.98618902e-04}, {1.33922796e+01, 2.03565693e+00, 4.68585610e-01, 3.41081098e-02, 6.18939695e-04}};
					case 67:
						return {{3.68242311e+00, 1.69007802e+00, 1.73174500e+00, 6.26287907e-02, 1.94016204e-04}, {1.36538200e+01, 3.38136411e+00, 5.00155091e-01, 3.36822011e-02, 6.00808591e-04}};
					case 68:
						return {{3.65319490e+00, 1.80632806e+00, 1.50050402e+00, 4.83703911e-02, 1.80868403e-04}, {1.33160105e+01, 2.66234398e+00, 4.23525512e-01, 2.84032505e-02, 5.81026427e-04}};
					case 69:
						return {{2.81324100e+00, 2.29724407e+00, 1.71350300e+00, 5.38761392e-02, 1.82238698e-04}, {1.63576298e+01, 3.92391706e+00, 4.60097700e-01, 2.98996102e-02, 5.65685797e-04}};
					case 70:
						return {{3.86458898e+00, 1.25103104e+00, 1.50890696e+00, 4.87409793e-02, 1.76625705e-04}, {1.09293900e+01, 2.05829406e+00, 4.10406888e-01, 2.79810801e-02, 5.49290911e-04}};
					case 71:
						return {{2.78422499e+00, 2.18659902e+00, 1.62505805e+00, 5.04949987e-02, 1.76768895e-04}, {1.05224104e+01, 4.15039682e+00, 4.16781306e-01, 2.84288507e-02, 5.34839579e-04}};
					case 72:
						return {{1.82888401e+00, 3.08626103e+00, 1.56239796e+00, 4.50708307e-02, 1.67344595e-04}, {1.09088402e+01, 4.48407412e+00, 3.92598897e-01, 2.59789601e-02, 5.18937595e-04}};
					case 73:
						return {{1.34293997e+00, 3.54547095e+00, 1.43781698e+00, 4.28609811e-02, 1.62209704e-04}, {1.22230902e+01, 3.98015690e+00, 3.64430696e-01, 2.49023698e-02, 5.04579220e-04}};
					case 74:
						return {{1.11627996e+00, 3.81305695e+00, 1.27436197e+00, 3.69497910e-02, 1.52609398e-04}, {1.21982803e+01, 3.51693296e+00, 3.22148204e-01, 2.25107409e-02, 4.89797501e-04}};
					case 75:
						return {{1.15703595e+00, 3.75507998e+00, 1.18611205e+00, 3.35603990e-02, 1.41277094e-04}, {1.14396601e+01, 3.10535192e+00, 3.00034910e-01, 2.07285509e-02, 4.75228589e-04}};
					case 76:
						return {{2.95070696e+00, 1.99730504e+00, 9.92647886e-01, 3.54575291e-02, 1.46259103e-04}, {5.83780003e+00, 1.72655404e+00, 2.71308392e-01, 2.16202103e-02, 4.64655401e-04}};
					case 77:
						return {{3.32513905e+00, 1.52911794e+00, 9.36988473e-01, 2.82717198e-02, 1.28554704e-04}, {4.62521505e+00, 1.43798900e+00, 2.50209898e-01, 1.81400497e-02, 4.49112791e-04}};
					case 78:
						return {{4.98159796e-01, 4.03055811e+00, 8.52239490e-01, 2.27640904e-02, 1.15693001e-04}, {1.65708504e+01, 2.17498207e+00, 2.20648602e-01, 1.56411994e-02, 4.34708607e-04}};
					case 79:
						return {{8.70827496e-01, 3.66290498e+00, 7.21173406e-01, 2.25834902e-02, 1.15950199e-04}, {8.26474285e+00, 1.76541495e+00, 1.96660295e-01, 1.56456195e-02, 4.24816913e-04}};
					case 80:
						return {{2.23930812e+00, 2.49683499e+00, 7.04915404e-01, 2.09938809e-02, 9.59475219e-05}, {4.66368914e+00, 1.36597097e+00, 1.94375604e-01, 1.40903797e-02, 4.09117813e-04}};
					case 81:
						return {{2.71713591e+00, 2.91902208e+00, 7.08489776e-01, 1.79738607e-02, 8.45634204e-05}, {8.71158314e+00, 1.45481896e+00, 1.87859505e-01, 1.23882899e-02, 3.94706090e-04}};
					case 82:
						return {{3.29438901e+00, 2.66650391e+00, 5.41207910e-01, 1.59845799e-02, 9.09762530e-05}, {7.26526308e+00, 1.15847301e+00, 1.53856993e-01, 1.19057801e-02, 3.87974986e-04}};
					case 83:
						return {{2.90801311e+00, 2.94222403e+00, 6.40052676e-01, 1.98365692e-02, 1.03154103e-04}, {7.08211708e+00, 1.40495801e+00, 1.71940193e-01, 1.36806397e-02, 3.83776205e-04}};
					case 84:
						return {{4.17358112e+00, 1.43807697e+00, 5.06418884e-01, 1.17866602e-02, 3.70279995e-05}, {3.63485503e+00, 7.84012079e-01, 1.42200202e-01, 8.50070175e-03, 3.35267512e-04}};
					case 85:
						return {{4.40665388e+00, 1.68819594e+00, 6.11086190e-01, 1.74588896e-02, 8.35977771e-05}, {4.43530703e+00, 9.95206296e-01, 1.63246498e-01, 1.18547501e-02, 3.61005805e-04}};
					case 86:
						return {{4.44254780e+00, 1.57980704e+00, 7.00253785e-01, 1.36486702e-02, 5.30898506e-05}, {3.88471198e+00, 1.23543596e+00, 1.67034298e-01, 9.34090279e-03, 3.36916302e-04}};
					case 87:
						return {{4.44422293e+00, 3.84529591e+00, 6.47594690e-01, 1.57298408e-02, 8.11556311e-05}, {1.55769997e+01, 1.69413495e+00, 1.57903805e-01, 1.09456200e-02, 3.44481290e-04}};
					case 88:
						return {{5.84577084e+00, 3.43635201e+00, 6.62875414e-01, 1.30639505e-02, 4.20832803e-05}, {1.44227695e+01, 1.58868802e+00, 1.56904906e-01, 8.72127432e-03, 3.15123296e-04}};
					case 89:
						return {{6.56823683e+00, 3.08116007e+00, 5.87002397e-01, 1.18503897e-02, 3.66903696e-05}, {1.23099003e+01, 1.34582400e+00, 1.42942101e-01, 8.05882178e-03, 3.02906294e-04}};
					case 90:
						return {{6.75877714e+00, 2.74244094e+00, 5.82863986e-01, 1.08115301e-02, 1.49991101e-05}, {9.66876411e+00, 1.26876402e+00, 1.39298707e-01, 7.15679117e-03, 2.59494089e-04}};
					case 91:
						return {{6.37304688e+00, 2.82305002e+00, 5.40732384e-01, 1.19372196e-02, 5.09942292e-05}, {9.53488636e+00, 1.19764698e+00, 1.32480398e-01, 8.31615925e-03, 3.02768894e-04}};
					case 92:
						return {{6.01751184e+00, 3.00574708e+00, 5.04366279e-01, 9.93358810e-03, 1.37084398e-05}, {9.77997780e+00, 1.16927600e+00, 1.23806901e-01, 6.66506588e-03, 2.45529693e-04}};
					case 93:
						return {{5.20271587e+00, 3.61388612e+00, 5.22116423e-01, 8.51191860e-03, 1.63058994e-05}, {1.04132996e+01, 1.40175903e+00, 1.19989701e-01, 6.06458494e-03, 2.46154610e-04}};
					case 94:
						return {{4.97225809e+00, 3.62590599e+00, 4.27345693e-01, 8.12398084e-03, 1.48722102e-05}, {1.17354097e+01, 1.18188906e+00, 1.05401397e-01, 5.83548006e-03, 2.37660599e-04}};
					case 95:
						return {{4.50712395e+00, 3.90453696e+00, 4.45498198e-01, 7.97858369e-03, 8.28391057e-06}, {1.28308401e+01, 1.25748503e+00, 1.06835604e-01, 5.59567707e-03, 2.09513906e-04}};
					case 96:
						return {{4.95550919e+00, 3.46335196e+00, 3.79075289e-01, 5.84759377e-03, 2.35737807e-05}, {9.56314087e+00, 1.07623196e+00, 9.16181728e-02, 4.78272792e-03, 2.42163704e-04}};
					case 97:
						return {{4.87071705e+00, 3.45333791e+00, 3.44667703e-01, 6.47798879e-03, 1.13584201e-05}, {9.75995636e+00, 9.98554826e-01, 8.70402008e-02, 4.92938887e-03, 2.12656101e-04}};
					case 98:
						return {{4.32416010e+00, 3.69497395e+00, 3.30316305e-01, 5.58224786e-03, 7.66277208e-06}, {1.12858105e+01, 1.00380003e+00, 8.26271027e-02, 4.33693221e-03, 1.90465493e-04}};
					case 99:
						return {{3.96653795e+00, 3.86160207e+00, 3.56552094e-01, 7.00159417e-03, 2.81896791e-05}, {1.22629700e+01, 1.04874003e+00, 8.83088931e-02, 5.31629194e-03, 2.38887806e-04}};
					case 100:
						return {{4.27992916e+00, 3.45720100e+00, 2.69533694e-01, 4.48193215e-03, 1.74949400e-05}, {9.27273178e+00, 8.65246892e-01, 6.92421496e-02, 3.92412720e-03, 2.13503794e-04}};
					case 101:
						return {{3.94026089e+00, 3.62439609e+00, 2.89608896e-01, 6.40732795e-03, 3.19458813e-05}, {1.02414103e+01, 8.98189783e-01, 7.49766305e-02, 5.15825208e-03, 2.35140804e-04}};
					case 102:
						return {{3.96523499e+00, 3.38761902e+00, 3.31892401e-01, 5.65199787e-03, 7.42224984e-06}, {8.30294800e+00, 8.89315307e-01, 8.13126490e-02, 4.15963400e-03, 1.77448994e-04}};
					case 103:
						return {{4.46439791e+00, 3.20867705e+00, 2.26712495e-01, 4.08070907e-03, -1.60853798e-03}, {7.96693707e+00, 7.42235303e-01, 5.82330413e-02, 1.39190501e-03, 8.74621095e-04}};
				}
				return {{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}};
			}

			// 10: Peng et al. parameterization for Ions- 5 Gaussians - [0, 4]
			LNL_Coef_cpu<T> load_feg_peng_ion_0_4(const dt_int32& Z, const dt_int32& charge)
			{
				switch(Z)
				{
					case 1:
					{
						switch(charge)
						{
							case -1:
								return {{1.400e-01, 6.490e-01, 1.370e+00, 3.370e-01, 7.870e-01}, {9.840e-01, 8.670e+00, 3.890e+01, 1.110e+02, 1.660e+02}};
						}
					}
					break;
					case 3:
					{
						switch(charge)
						{
							case 1:
								return {{4.600e-03, 1.650e-02, 4.350e-02, 6.490e-02, 2.700e-02}, {3.580e-02, 2.390e-01, 8.790e-01, 2.640e+00, 7.090e+00}};
						}
					}
					break;
					case 4:
					{
						switch(charge)
						{
							case 2:
								return {{3.400e-03, 1.030e-02, 2.330e-02, 3.250e-02, 1.200e-02}, {2.670e-02, 1.620e-01, 5.310e-01, 1.480e+00, 3.880e+00}};
						}
					}
					break;
					case 8:
					{
						switch(charge)
						{
							case -1:
								return {{2.050e-01, 6.280e-01, 1.170e+00, 1.030e+00, 2.900e-01}, {3.970e-01, 2.640e+00, 8.800e+00, 2.710e+01, 9.180e+01}};
							case -2:
								return {{4.210e-02, 2.100e-01, 8.520e-01, 1.820e+00, 1.170e+00}, {6.090e-02, 5.590e-01, 2.960e+00, 1.150e+01, 3.770e+01}};
						}
					}
					break;
					case 9:
					{
						switch(charge)
						{
							case -1:
								return {{1.340e-01, 3.910e-01, 8.140e-01, 9.280e-01, 3.470e-01}, {2.280e-01, 1.470e+00, 4.680e+00, 1.320e+01, 3.600e+01}};
						}
					}
					break;
					case 11:
					{
						switch(charge)
						{
							case 1:
								return {{2.560e-02, 9.190e-02, 2.970e-01, 5.140e-01, 1.990e-01}, {3.970e-02, 2.870e-01, 1.180e+00, 3.750e+00, 1.080e+01}};
						}
					}
					break;
					case 12:
					{
						switch(charge)
						{
							case 2:
								return {{2.100e-02, 6.720e-02, 1.980e-01, 3.680e-01, 1.740e-01}, {3.310e-02, 2.220e-01, 8.380e-01, 2.480e+00, 6.750e+00}};
						}
					}
					break;
					case 13:
					{
						switch(charge)
						{
							case 3:
								return {{1.920e-02, 5.790e-02, 1.630e-01, 2.840e-01, 1.140e-01}, {3.060e-02, 1.980e-01, 7.130e-01, 2.040e+00, 5.250e+00}};
						}
					}
					break;
					case 14:
					{
						switch(charge)
						{
							case 4:
								return {{1.920e-01, 2.890e-01, 1.000e-01, -7.280e-02, 1.200e-03}, {3.590e-01, 1.960e+00, 9.340e+00, 1.110e+01, 1.340e+01}};
						}
					}
					break;
					case 17:
					{
						switch(charge)
						{
							case -1:
								return {{2.650e-01, 5.960e-01, 1.600e+00, 2.690e+00, 1.230e+00}, {2.520e-01, 1.560e+00, 6.210e+00, 1.780e+01, 4.780e+01}};
						}
					}
					break;
					case 19:
					{
						switch(charge)
						{
							case 1:
								return {{1.990e-01, 3.960e-01, 9.280e-01, 1.450e+00, 4.500e-01}, {1.920e-01, 1.100e+00, 3.910e+00, 9.750e+00, 2.340e+01}};
						}
					}
					break;
					case 20:
					{
						switch(charge)
						{
							case 2:
								return {{1.640e-01, 3.270e-01, 7.430e-01, 1.160e+00, 3.070e-01}, {1.570e-01, 8.940e-01, 3.150e+00, 7.670e+00, 1.770e+01}};
						}
					}
					break;
					case 21:
					{
						switch(charge)
						{
							case 3:
								return {{1.630e-01, 3.070e-01, 7.160e-01, 8.800e-01, 1.390e-01}, {1.570e-01, 8.990e-01, 3.060e+00, 7.050e+00, 1.610e+01}};
						}
					}
					break;
					case 22:
					{
						switch(charge)
						{
							case 2:
								return {{3.990e-01, 1.040e+00, 1.210e+00, -7.970e-02, 3.520e-01}, {3.760e-01, 2.740e+00, 8.100e+00, 1.420e+01, 2.320e+01}};
							case 3:
								return {{3.640e-01, 9.190e-01, 1.350e+00, -9.330e-01, 5.890e-01}, {3.640e-01, 2.670e+00, 8.180e+00, 1.180e+01, 1.490e+01}};
							case 4:
								return {{1.160e-01, 2.560e-01, 5.650e-01, 7.720e-01, 1.320e-01}, {1.080e-01, 6.550e-01, 2.380e+00, 5.510e+00, 1.230e+01}};
						}
					}
					break;
					case 23:
					{
						switch(charge)
						{
							case 2:
								return {{3.170e-01, 9.390e-01, 1.490e+00, -1.310e+00, 1.470e+00}, {2.690e-01, 2.090e+00, 7.220e+00, 1.520e+01, 1.760e+01}};
							case 3:
								return {{3.410e-01, 8.050e-01, 9.420e-01, 7.830e-02, 1.560e-01}, {3.210e-01, 2.230e+00, 5.990e+00, 1.340e+01, 1.690e+01}};
							case 5:
								return {{3.670e-02, 1.240e-01, 2.440e-01, 7.230e-01, 4.350e-01}, {3.300e-02, 2.220e-01, 8.240e-01, 2.800e+00, 6.700e+00}};
						}
					}
					break;
					case 24:
					{
						switch(charge)
						{
							case 2:
								return {{2.370e-01, 6.340e-01, 1.230e+00, 7.130e-01, 8.590e-02}, {1.770e-01, 1.350e+00, 4.300e+00, 1.220e+01, 3.900e+01}};
							case 3:
								return {{3.930e-01, 1.050e+00, 1.620e+00, -1.150e+00, 4.070e-01}, {3.590e-01, 2.570e+00, 8.680e+00, 1.100e+01, 1.580e+01}};
							case 4:
								return {{1.320e-01, 2.920e-01, 7.030e-01, 6.920e-01, 9.590e-02}, {1.090e-01, 6.950e-01, 2.390e+00, 5.650e+00, 1.470e+01}};
						}
					}
					break;
					case 25:
					{
						switch(charge)
						{
							case 2:
								return {{5.760e-02, 2.100e-01, 6.040e-01, 1.320e+00, 6.590e-01}, {3.980e-02, 2.840e-01, 1.290e+00, 4.230e+00, 1.450e+01}};
							case 3:
								return {{1.160e-01, 5.230e-01, 8.810e-01, 5.890e-01, 2.140e-01}, {1.170e-02, 8.760e-01, 3.060e+00, 6.440e+00, 1.430e+01}};
							case 4:
								return {{3.810e-01, 1.830e+00, -1.330e+00, 9.950e-01, 6.180e-02}, {3.540e-01, 2.720e+00, 3.470e+00, 5.470e+00, 1.610e+01}};
						}
					}
					break;
					case 26:
					{
						switch(charge)
						{
							case 2:
								return {{3.070e-01, 8.380e-01, 1.110e+00, 2.800e-01, 2.770e-01}, {2.300e-01, 1.620e+00, 4.870e+00, 1.070e+01, 1.920e+01}};
							case 3:
								return {{1.980e-01, 3.870e-01, 8.890e-01, 7.090e-01, 1.170e-01}, {1.540e-01, 8.930e-01, 2.620e+00, 6.650e+00, 1.800e+01}};
						}
					}
					break;
					case 27:
					{
						switch(charge)
						{
							case 2:
								return {{2.130e-01, 4.880e-01, 9.980e-01, 8.280e-01, 2.300e-01}, {1.480e-01, 9.390e-01, 2.780e+00, 7.310e+00, 2.070e+01}};
							case 3:
								return {{3.310e-01, 4.870e-01, 7.290e-01, 6.080e-01, 1.310e-01}, {2.670e-01, 1.410e+00, 2.890e+00, 6.450e+00, 1.580e+01}};
						}
					}
					break;
					case 28:
					{
						switch(charge)
						{
							case 2:
								return {{3.380e-01, 9.820e-01, 1.320e+00, -3.560e+00, 3.620e+00}, {2.370e-01, 1.670e+00, 5.730e+00, 1.140e+01, 1.210e+01}};
							case 3:
								return {{3.470e-01, 8.770e-01, 7.900e-01, 5.380e-02, 1.920e-01}, {2.600e-01, 1.710e+00, 4.750e+00, 7.510e+00, 1.300e+01}};
						}
					}
					break;
					case 29:
					{
						switch(charge)
						{
							case 1:
								return {{3.120e-01, 8.120e-01, 1.110e+00, 7.940e-01, 2.570e-01}, {2.010e-01, 1.310e+00, 3.800e+00, 1.050e+01, 2.820e+01}};
							case 2:
								return {{2.240e-01, 5.440e-01, 9.700e-01, 7.270e-01, 1.820e-01}, {1.450e-01, 9.330e-01, 2.690e+00, 7.110e+00, 1.940e+01}};
						}
					}
					break;
					case 30:
					{
						switch(charge)
						{
							case 2:
								return {{2.520e-01, 6.000e-01, 9.170e-01, 6.630e-01, 1.610e-01}, {1.610e-01, 1.010e+00, 2.760e+00, 7.080e+00, 1.900e+01}};
						}
					}
					break;
					case 31:
					{
						switch(charge)
						{
							case 3:
								return {{3.910e-01, 9.470e-01, 6.900e-01, 7.090e-02, 6.530e-02}, {2.640e-01, 1.650e+00, 4.820e+00, 1.070e+01, 1.520e+01}};
						}
					}
					break;
					case 32:
					{
						switch(charge)
						{
							case 4:
								return {{3.460e-01, 8.300e-01, 5.990e-01, 9.490e-02, -2.170e-02}, {2.320e-01, 1.450e+00, 4.080e+00, 1.320e+01, 2.950e+01}};
						}
					}
					break;
					case 35:
					{
						switch(charge)
						{
							case -1:
								return {{1.250e-01, 5.630e-01, 1.430e+00, 3.520e+00, 3.220e+00}, {5.300e-02, 4.690e-01, 2.150e+00, 1.110e+01, 3.890e+01}};
						}
					}
					break;
					case 37:
					{
						switch(charge)
						{
							case 1:
								return {{3.680e-01, 8.840e-01, 1.140e+00, 2.260e+00, 8.810e-01}, {1.870e-01, 1.120e+00, 3.980e+00, 1.090e+01, 2.660e+01}};
						}
					}
					break;
					case 38:
					{
						switch(charge)
						{
							case 2:
								return {{3.460e-01, 8.040e-01, 9.880e-01, 1.890e+00, 6.090e-01}, {1.760e-01, 1.040e+00, 3.590e+00, 9.320e+00, 2.140e+01}};
						}
					}
					break;
					case 39:
					{
						switch(charge)
						{
							case 3:
								return {{4.650e-01, 9.230e-01, 2.410e+00, -2.310e+00, 2.480e+00}, {2.400e-01, 1.430e+00, 6.450e+00, 9.970e+00, 1.220e+01}};
						}
					}
					break;
					case 40:
					{
						switch(charge)
						{
							case 4:
								return {{2.340e-01, 6.420e-01, 7.470e-01, 1.470e+00, 3.770e-01}, {1.130e-01, 7.360e-01, 2.540e+00, 6.720e+00, 1.470e+01}};
						}
					}
					break;
					case 41:
					{
						switch(charge)
						{
							case 3:
								return {{3.770e-01, 7.490e-01, 1.290e+00, 1.610e+00, 4.810e-01}, {1.840e-01, 1.020e+00, 3.800e+00, 9.440e+00, 2.570e+01}};
							case 5:
								return {{8.280e-02, 2.710e-01, 6.540e-01, 1.240e+00, 8.290e-01}, {3.690e-02, 2.610e-01, 9.570e-01, 3.940e+00, 9.440e+00}};
						}
					}
					break;
					case 42:
					{
						switch(charge)
						{
							case 3:
								return {{4.010e-01, 7.560e-01, 1.380e+00, 1.580e+00, 4.970e-01}, {1.910e-01, 1.060e+00, 3.840e+00, 9.380e+00, 2.460e+01}};
							case 5:
								return {{4.790e-01, 8.460e-01, 1.560e+01, -1.520e+01, 1.600e+00}, {2.410e-01, 1.460e+00, 6.790e+00, 7.130e+00, 1.040e+01}};
							case 6:
								return {{2.030e-01, 5.670e-01, 6.460e-01, 1.160e+00, 1.710e-01}, {9.710e-02, 6.470e-01, 2.280e+00, 5.610e+00, 1.240e+01}};
						}
					}
					break;
					case 44:
					{
						switch(charge)
						{
							case 3:
								return {{4.280e-01, 7.730e-01, 1.550e+00, 1.460e+00, 4.860e-01}, {1.910e-01, 1.090e+00, 3.820e+00, 9.080e+00, 2.170e+01}};
							case 4:
								return {{2.820e-01, 6.530e-01, 1.140e+00, 1.530e+00, 4.180e-01}, {1.250e-01, 7.530e-01, 2.850e+00, 7.010e+00, 1.750e+01}};
						}
					}
					break;
					case 45:
					{
						switch(charge)
						{
							case 3:
								return {{3.520e-01, 7.230e-01, 1.500e+00, 1.630e+00, 4.990e-01}, {1.510e-01, 8.780e-01, 3.280e+00, 8.160e+00, 2.070e+01}};
							case 4:
								return {{3.970e-01, 7.250e-01, 1.510e+00, 1.190e+00, 2.510e-01}, {1.770e-01, 1.010e+00, 3.620e+00, 8.560e+00, 1.890e+01}};
						}
					}
					break;
					case 46:
					{
						switch(charge)
						{
							case 2:
								return {{9.350e-01, 3.110e+00, 2.460e+01, -4.360e+01, 2.110e+01}, {3.930e-01, 4.060e+00, 4.310e+01, 5.400e+01, 6.980e+01}};
							case 4:
								return {{3.480e-01, 6.400e-01, 1.220e+00, 1.450e+00, 4.270e-01}, {1.510e-01, 8.320e-01, 2.850e+00, 6.590e+00, 1.560e+01}};
						}
					}
					break;
					case 47:
					{
						switch(charge)
						{
							case 1:
								return {{5.030e-01, 9.400e-01, 2.170e+00, 1.990e+00, 7.260e-01}, {1.990e-01, 1.190e+00, 4.050e+00, 1.130e+01, 3.240e+01}};
							case 2:
								return {{4.310e-01, 7.560e-01, 1.720e+00, 1.780e+00, 5.260e-01}, {1.750e-01, 9.790e-01, 3.300e+00, 8.240e+00, 2.140e+01}};
						}
					}
					break;
					case 48:
					{
						switch(charge)
						{
							case 2:
								return {{4.250e-01, 7.450e-01, 1.730e+00, 1.740e+00, 4.870e-01}, {1.680e-01, 9.440e-01, 3.140e+00, 7.840e+00, 2.040e+01}};
						}
					}
					break;
					case 49:
					{
						switch(charge)
						{
							case 3:
								return {{4.170e-01, 7.550e-01, 1.590e+00, 1.360e+00, 4.510e-01}, {1.640e-01, 9.600e-01, 3.080e+00, 7.030e+00, 1.610e+01}};
						}
					}
					break;
					case 50:
					{
						switch(charge)
						{
							case 2:
								return {{7.970e-01, 2.130e+00, 2.150e+00, -1.640e+00, 2.720e+00}, {3.170e-01, 2.510e+00, 9.040e+00, 2.420e+01, 2.640e+01}};
							case 4:
								return {{2.610e-01, 6.420e-01, 1.530e+00, 1.360e+00, 1.770e-01}, {9.570e-02, 6.250e-01, 2.510e+00, 6.310e+00, 1.590e+01}};
						}
					}
					break;
					case 51:
					{
						switch(charge)
						{
							case 3:
								return {{5.520e-01, 1.140e+00, 1.870e+00, 1.360e+00, 4.140e-01}, {2.120e-01, 1.420e+00, 4.210e+00, 1.250e+01, 2.900e+01}};
							case 5:
								return {{3.770e-01, 5.880e-01, 1.220e+00, 1.180e+00, 2.440e-01}, {1.510e-01, 8.120e-01, 2.400e+00, 5.270e+00, 1.190e+01}};
						}
					}
					break;
					case 53:
					{
						switch(charge)
						{
							case -1:
								return {{9.010e-01, 2.800e+00, 5.610e+00, -8.690e+00, 1.260e+01}, {3.120e-01, 2.590e+00, 1.410e+01, 3.440e+01, 3.950e+01}};
						}
					}
					break;
					case 55:
					{
						switch(charge)
						{
							case 1:
								return {{5.870e-01, 1.400e+00, 1.870e+00, 3.480e+00, 1.670e+00}, {2.000e-01, 1.380e+00, 4.120e+00, 1.300e+01, 3.180e+01}};
						}
					}
					break;
					case 56:
					{
						switch(charge)
						{
							case 2:
								return {{7.330e-01, 2.050e+00, 2.300e+01, -1.520e+02, 1.340e+02}, {2.580e-01, 1.960e+00, 1.180e+01, 1.440e+01, 1.490e+01}};
						}
					}
					break;
					case 57:
					{
						switch(charge)
						{
							case 3:
								return {{4.930e-01, 1.100e+00, 1.500e+00, 2.700e+00, 1.080e+00}, {1.670e-01, 1.110e+00, 3.110e+00, 9.610e+00, 2.120e+01}};
						}
					}
					break;
					case 58:
					{
						switch(charge)
						{
							case 3:
								return {{5.600e-01, 1.350e+00, 1.590e+00, 2.630e+00, 7.060e-01}, {1.900e-01, 1.300e+00, 3.930e+00, 1.070e+01, 2.380e+01}};
							case 4:
								return {{4.830e-01, 1.090e+00, 1.340e+00, 2.450e+00, 7.970e-01}, {1.650e-01, 1.100e+00, 3.020e+00, 8.850e+00, 1.880e+01}};
						}
					}
					break;
					case 59:
					{
						switch(charge)
						{
							case 3:
								return {{6.630e-01, 1.730e+00, 2.350e+00, 3.510e-01, 1.590e+00}, {2.260e-01, 1.610e+00, 6.330e+00, 1.100e+01, 1.690e+01}};
							case 4:
								return {{5.210e-01, 1.190e+00, 1.330e+00, 2.360e+00, 6.900e-01}, {1.770e-01, 1.170e+00, 3.280e+00, 8.940e+00, 1.930e+01}};
						}
					}
					break;
					case 60:
					{
						switch(charge)
						{
							case 3:
								return {{5.010e-01, 1.180e+00, 1.450e+00, 2.530e+00, 9.200e-01}, {1.620e-01, 1.080e+00, 3.060e+00, 8.800e+00, 1.960e+01}};
						}
					}
					break;
					case 61:
					{
						switch(charge)
						{
							case 3:
								return {{4.960e-01, 1.200e+00, 1.470e+00, 2.430e+00, 9.430e-01}, {1.560e-01, 1.050e+00, 3.070e+00, 8.560e+00, 1.920e+01}};
						}
					}
					break;
					case 62:
					{
						switch(charge)
						{
							case 3:
								return {{5.180e-01, 1.240e+00, 1.430e+00, 2.400e+00, 7.810e-01}, {1.630e-01, 1.080e+00, 3.110e+00, 8.520e+00, 1.910e+01}};
						}
					}
					break;
					case 63:
					{
						switch(charge)
						{
							case 2:
								return {{6.130e-01, 1.530e+00, 1.840e+00, 2.460e+00, 7.140e-01}, {1.900e-01, 1.270e+00, 4.180e+00, 1.070e+01, 2.620e+01}};
							case 3:
								return {{4.960e-01, 1.210e+00, 1.450e+00, 2.360e+00, 7.740e-01}, {1.520e-01, 1.010e+00, 2.950e+00, 8.180e+00, 1.850e+01}};
						}
					}
					break;
					case 64:
					{
						switch(charge)
						{
							case 3:
								return {{4.900e-01, 1.190e+00, 1.420e+00, 2.300e+00, 7.950e-01}, {1.480e-01, 9.740e-01, 2.810e+00, 7.780e+00, 1.770e+01}};
						}
					}
					break;
					case 65:
					{
						switch(charge)
						{
							case 3:
								return {{5.030e-01, 1.220e+00, 1.420e+00, 2.240e+00, 7.100e-01}, {1.500e-01, 9.820e-01, 2.860e+00, 7.770e+00, 1.770e+01}};
						}
					}
					break;
					case 66:
					{
						switch(charge)
						{
							case 3:
								return {{5.030e-01, 1.240e+00, 1.440e+00, 2.170e+00, 6.430e-01}, {1.480e-01, 9.700e-01, 2.880e+00, 7.730e+00, 1.760e+01}};
						}
					}
					break;
					case 67:
					{
						switch(charge)
						{
							case 3:
								return {{4.560e-01, 1.170e+00, 1.430e+00, 2.150e+00, 6.920e-01}, {1.290e-01, 8.690e-01, 2.610e+00, 7.240e+00, 1.670e+01}};
						}
					}
					break;
					case 68:
					{
						switch(charge)
						{
							case 3:
								return {{5.220e-01, 1.280e+00, 1.460e+00, 2.050e+00, 5.080e-01}, {1.500e-01, 9.640e-01, 2.930e+00, 7.720e+00, 1.780e+01}};
						}
					}
					break;
					case 69:
					{
						switch(charge)
						{
							case 3:
								return {{4.750e-01, 1.200e+00, 1.420e+00, 2.050e+00, 5.840e-01}, {1.320e-01, 8.640e-01, 2.600e+00, 7.090e+00, 1.660e+01}};
						}
					}
					break;
					case 70:
					{
						switch(charge)
						{
							case 2:
								return {{5.080e-01, 1.370e+00, 1.760e+00, 2.230e+00, 5.840e-01}, {1.360e-01, 9.220e-01, 3.120e+00, 8.720e+00, 2.370e+01}};
							case 3:
								return {{4.980e-01, 1.220e+00, 1.390e+00, 1.970e+00, 5.590e-01}, {1.380e-01, 8.810e-01, 2.630e+00, 6.990e+00, 1.630e+01}};
						}
					}
					break;
					case 71:
					{
						switch(charge)
						{
							case 3:
								return {{4.830e-01, 1.210e+00, 1.410e+00, 1.940e+00, 5.220e-01}, {1.310e-01, 8.450e-01, 2.570e+00, 6.880e+00, 1.620e+01}};
						}
					}
					break;
					case 72:
					{
						switch(charge)
						{
							case 4:
								return {{5.220e-01, 1.220e+00, 1.370e+00, 1.680e+00, 3.120e-01}, {1.450e-01, 8.960e-01, 2.740e+00, 6.910e+00, 1.610e+01}};
						}
					}
					break;
					case 73:
					{
						switch(charge)
						{
							case 5:
								return {{5.690e-01, 1.260e+00, 9.790e-01, 1.290e+00, 5.510e-01}, {1.610e-01, 9.720e-01, 2.760e+00, 5.400e+00, 1.090e+01}};
						}
					}
					break;
					case 74:
					{
						switch(charge)
						{
							case 6:
								return {{1.810e-01, 8.730e-01, 1.180e+00, 1.480e+00, 5.620e-01}, {1.180e-02, 4.420e-01, 1.520e+00, 4.350e+00, 9.420e+00}};
						}
					}
					break;
					case 76:
					{
						switch(charge)
						{
							case 4:
								return {{5.860e-01, 1.310e+00, 1.630e+00, 1.710e+00, 5.400e-01}, {1.550e-01, 9.380e-01, 3.190e+00, 7.840e+00, 1.930e+01}};
						}
					}
					break;
					case 77:
					{
						switch(charge)
						{
							case 3:
								return {{6.920e-01, 1.370e+00, 1.800e+00, 1.970e+00, 8.040e-01}, {1.820e-01, 1.040e+00, 3.470e+00, 8.510e+00, 2.120e+01}};
							case 4:
								return {{6.530e-01, 1.290e+00, 1.500e+00, 1.740e+00, 6.830e-01}, {1.740e-01, 9.920e-01, 3.140e+00, 7.220e+00, 1.720e+01}};
						}
					}
					break;
					case 78:
					{
						switch(charge)
						{
							case 2:
								return {{8.720e-01, 1.680e+00, 2.630e+00, 1.930e+00, 4.750e-01}, {2.230e-01, 1.350e+00, 4.990e+00, 1.360e+01, 3.300e+01}};
							case 4:
								return {{5.500e-01, 1.210e+00, 1.620e+00, 1.950e+00, 6.100e-01}, {1.420e-01, 8.330e-01, 2.810e+00, 7.210e+00, 1.770e+01}};
						}
					}
					break;
					case 79:
					{
						switch(charge)
						{
							case 1:
								return {{8.110e-01, 1.570e+00, 2.630e+00, 2.680e+00, 9.980e-01}, {2.010e-01, 1.180e+00, 4.250e+00, 1.210e+01, 3.440e+01}};
							case 3:
								return {{7.220e-01, 1.390e+00, 1.940e+00, 1.940e+00, 6.990e-01}, {1.840e-01, 1.060e+00, 3.580e+00, 8.560e+00, 2.040e+01}};
						}
					}
					break;
					case 80:
					{
						switch(charge)
						{
							case 1:
								return {{7.960e-01, 1.560e+00, 2.720e+00, 2.760e+00, 1.180e+00}, {1.940e-01, 1.140e+00, 4.210e+00, 1.240e+01, 3.620e+01}};
							case 2:
								return {{7.730e-01, 1.490e+00, 2.450e+00, 2.230e+00, 5.700e-01}, {1.910e-01, 1.120e+00, 4.000e+00, 1.080e+01, 2.760e+01}};
						}
					}
					break;
					case 81:
					{
						switch(charge)
						{
							case 1:
								return {{8.200e-01, 1.570e+00, 2.780e+00, 2.820e+00, 1.310e+00}, {1.970e-01, 1.160e+00, 4.230e+00, 1.270e+01, 3.570e+01}};
							case 3:
								return {{8.360e-01, 1.430e+00, 3.940e-01, 2.510e+00, 1.500e+00}, {2.080e-01, 1.200e+00, 2.570e+00, 4.860e+00, 1.350e+01}};
						}
					}
					break;
					case 82:
					{
						switch(charge)
						{
							case 2:
								return {{7.550e-01, 1.440e+00, 2.480e+00, 2.450e+00, 1.030e+00}, {1.810e-01, 1.050e+00, 3.750e+00, 1.060e+01, 2.790e+01}};
							case 4:
								return {{5.830e-01, 1.140e+00, 1.600e+00, 2.060e+00, 6.620e-01}, {1.440e-01, 7.960e-01, 2.580e+00, 6.220e+00, 1.480e+01}};
						}
					}
					break;
					case 83:
					{
						switch(charge)
						{
							case 3:
								return {{7.080e-01, 1.350e+00, 2.280e+00, 2.180e+00, 7.970e-01}, {1.700e-01, 9.810e-01, 3.440e+00, 9.410e+00, 2.370e+01}};
							case 5:
								return {{6.540e-01, 1.180e+00, 1.250e+00, 1.660e+00, 7.780e-01}, {1.620e-01, 9.050e-01, 2.680e+00, 5.140e+00, 1.120e+01}};
						}
					}
					break;
					case 88:
					{
						switch(charge)
						{
							case 2:
								return {{9.110e-01, 1.650e+00, 2.530e+00, 3.620e+00, 1.580e+00}, {2.040e-01, 1.260e+00, 4.030e+00, 1.260e+01, 3.000e+01}};
						}
					}
					break;
					case 89:
					{
						switch(charge)
						{
							case 3:
								return {{9.150e-01, 1.640e+00, 2.260e+00, 3.180e+00, 1.250e+00}, {2.050e-01, 1.280e+00, 3.920e+00, 1.130e+01, 2.510e+01}};
						}
					}
					break;
					case 92:
					{
						switch(charge)
						{
							case 3:
								return {{1.140e+00, 2.480e+00, 3.610e+00, 1.130e+00, 9.000e-01}, {2.500e-01, 1.840e+00, 7.390e+00, 1.800e+01, 2.270e+01}};
							case 4:
								return {{1.090e+00, 2.320e+00, 1.200e+01, -9.110e+00, 2.150e+00}, {2.430e-01, 1.750e+00, 7.790e+00, 8.310e+00, 1.650e+01}};
							case 6:
								return {{6.870e-01, 1.140e+00, 1.830e+00, 2.530e+00, 9.570e-01}, {1.540e-01, 8.610e-01, 2.580e+00, 7.700e+00, 1.590e+01}};
						}
					}
					break;
				}

				return {{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}};
			}

			/***************************************************************************************/
			Atomic_Coef_cpu<T> atomic_coef_gaussian(const LNL_Coef_cpu<T>& feg)		
			{
				Atomic_Coef_cpu<T> atomic_coef(feg.size(), T(0));

				for(auto ik = 0; ik < feg.size(); ik++)
				{
					auto cl = feg.cl[ik];
					auto cnl = (fcn_is_zero(cl))?T(0):feg.cnl[ik]/T(4);

					if (fcn_is_nzero(cl, cnl))
					{
						atomic_coef.feg.cl[ik] = cl;
						atomic_coef.feg.cnl[ik] = cnl;

						atomic_coef.fxg.cl[ik] = c_2Pi2a0<T>*cl;
						atomic_coef.fxg.cnl[ik] = cnl;

						atomic_coef.pr.cl[ik] = T(0.5)*c_2Pi2a0<T>*cl*pow(c_pi<T>/cnl, T(3.5));
						atomic_coef.pr.cnl[ik] = c_pi2<T>/cnl;

						atomic_coef.vr.cl[ik] = c_pot_factor<T>*cl*pow(c_pi<T>/cnl, T(1.5));
						atomic_coef.vr.cnl[ik] = c_pi2<T>/cnl;

						atomic_coef.vzp.cl[ik] = c_pot_factor<T>*c_pi<T>*cl/cnl;
						atomic_coef.vzp.cnl[ik] = c_pi2<T>/cnl;
					}
				}

				return atomic_coef;
			}

			// 1: doyle and Turner parameterization - 4 Gaussians - [0, 4]
			Atomic_Coef_cpu<T> atomic_coef_doyle_neutral_0_4(const dt_int32& Z)		
			{
				auto feg_coef = load_feg_doyle_neutral_0_4(Z);

				auto atomic_coef = atomic_coef_gaussian(feg_coef);
				atomic_coef.atomic_pot_parm_typ = eappt_doyle_0_4;
				atomic_coef.Z = Z;
				atomic_coef.charge = 0;

				return atomic_coef;
			}

			// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
			Atomic_Coef_cpu<T> atomic_coef_peng_neutral_0_4(const dt_int32& Z)		
			{
				auto feg_coef = load_feg_peng_neutral_0_4(Z);

				auto atomic_coef = atomic_coef_gaussian(feg_coef);

				atomic_coef.atomic_pot_parm_typ = eappt_peng_0_4;
				atomic_coef.Z = Z;
				atomic_coef.charge = 0;

				return atomic_coef;
			}

			// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
			Atomic_Coef_cpu<T> atomic_coef_peng_neutral_0_12(const dt_int32& Z)		
			{
				auto feg_coef = load_feg_peng_neutral_0_12(Z);

				if (feg_coef.size()==0)
				{
					return {};
				}

				auto atomic_coef = atomic_coef_gaussian(feg_coef);
				atomic_coef.atomic_pot_parm_typ = eappt_peng_0_12;
				atomic_coef.Z = Z;
				atomic_coef.charge = 0;

				return atomic_coef;
			}

			// 4: Kirkland parameterization - 3 yukawa + 3 Gaussians - [0, 12]			
			Atomic_Coef_cpu<T> atomic_coef_kirkland_neutral_0_12(const dt_int32& Z)		
			{
				auto feg_coef = load_feg_kirkland_neutral_0_12(Z);

				Atomic_Coef_cpu<T> atomic_coef(feg_coef.size(), T(0));
				atomic_coef.atomic_pot_parm_typ = eappt_kirkland_0_12;
				atomic_coef.Z = Z;
				atomic_coef.charge = 0;

				for(auto ik = 0; ik < 3; ik++)
				{
					auto cl = feg_coef.cl[ik];
					auto cnl = (fcn_is_zero(cl))?0.0:feg_coef.cnl[ik];

					if (fcn_is_nzero(cl, cnl))
					{
						atomic_coef.feg.cl[ik] = cl;
						atomic_coef.feg.cnl[ik] = cnl;

						atomic_coef.fxg.cl[ik] = c_2Pi2a0<T>*cl;
						atomic_coef.fxg.cnl[ik] = cnl;

						atomic_coef.pr.cl[ik] = c_2Pi2a0<T>*c_pi<T>*cl*cnl;
						atomic_coef.pr.cnl[ik] = c_2pi<T>*::sqrt(cnl);

						atomic_coef.vr.cl[ik] = c_pot_factor<T>*c_pi<T>*cl;
						atomic_coef.vr.cnl[ik] = c_2pi<T>*::sqrt(cnl);

						atomic_coef.vzp.cl[ik] = c_pot_factor<T>*c_2pi<T>*cl;
						atomic_coef.vzp.cnl[ik] = c_2pi<T>*::sqrt(cnl);
					}
				}

				for(auto ik = 3; ik < 6; ik++)
				{
					auto cl = feg_coef.cl[ik];
					auto cnl = (fcn_is_zero(cl))?T(0):feg_coef.cnl[ik];

					if (fcn_is_nzero(cl, cnl))
					{
						atomic_coef.feg.cl[ik] = cl;
						atomic_coef.feg.cnl[ik] = cnl;

						atomic_coef.fxg.cl[ik] = c_2Pi2a0<T>*cl;
						atomic_coef.fxg.cnl[ik] = cnl;

						atomic_coef.pr.cl[ik] = T(0.5)*c_2Pi2a0<T>*cl*pow(c_pi<T>/cnl, T(3.5));
						atomic_coef.pr.cnl[ik] = c_pi2<T>/cnl;

						atomic_coef.vr.cl[ik] = c_pot_factor<T>*cl*pow(c_pi<T>/cnl, T(1.5));
						atomic_coef.vr.cnl[ik] = c_pi2<T>/cnl;

						atomic_coef.vzp.cl[ik] = c_pot_factor<T>*c_pi<T>*cl/cnl;
						atomic_coef.vzp.cnl[ik] = c_pi2<T>/cnl;
					}
				}

				return atomic_coef;
			}
			
			// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
			Atomic_Coef_cpu<T> atomic_coef_weickenmeier_neutral_0_12(const dt_int32& Z)		
			{
				auto feg_coef = load_feg_weickenmeier_neutral_0_12(Z);

				Atomic_Coef_cpu<T> atomic_coef(feg_coef.size(), T(0));
				atomic_coef.atomic_pot_parm_typ = eappt_weickenmeier_0_12;
				atomic_coef.Z = Z;
				atomic_coef.charge = 0;

				for(auto ik = 0; ik < 6; ik++)
				{
					// 1.0/c_2Pi2a0<T> = 4.0*0.0239336609991378
					auto cl = Z/(c_2Pi2a0<T>*T(3)*(T(1)+feg_coef.cl[ik]));
					cl *= (ik>= 3)?feg_coef.cl[ik]:T(1);

					auto cnl = (fcn_is_zero(cl))?0.0:feg_coef.cnl[ik]/T(4);
					cl = (fcn_is_zero(cnl))?T(0):cl;

					if (fcn_is_nzero(cl, cnl))
					{
						atomic_coef.feg.cl[ik] = cl;
						atomic_coef.feg.cnl[ik] = cnl;

						atomic_coef.fxg.cl[ik] = c_2Pi2a0<T>*cl;
						atomic_coef.fxg.cnl[ik] = cnl;

						atomic_coef.pr.cl[ik] = c_2Pi2a0<T>*cl*pow(c_pi<T>/cnl, T(1.5));
						atomic_coef.pr.cnl[ik] = c_pi2<T>/cnl;

						atomic_coef.vr.cl[ik] = c_pot_factor<T>*c_pi<T>*cl;
						atomic_coef.vr.cnl[ik] = c_pi<T>/::sqrt(cnl);

						// there is not analytic expression for the projected potential
						atomic_coef.vzp.cl[ik] = atomic_coef.vr.cl[ik];
						atomic_coef.vzp.cnl[ik] = atomic_coef.vr.cnl[ik];
					}
				}

				return atomic_coef;
			}
			
			// 6: Lobato parameterization - Hydrogen functions - [0, 12]
			Atomic_Coef_cpu<T> atomic_coef_lobato_neutral_0_12(const dt_int32& Z)		
			{
				auto feg_coef = load_feg_lobato_neutral_0_12(Z);

				Atomic_Coef_cpu<T> atomic_coef(feg_coef.size(), T(0));
				atomic_coef.atomic_pot_parm_typ = eappt_lobato_0_12;
				atomic_coef.Z = Z;
				atomic_coef.charge = 0;

				for(auto ik = 0; ik < 5; ik++)
				{
					auto cl = (fcn_is_zero(feg_coef.cl[ik]))?T(0):feg_coef.cl[ik];
					auto cnl = (fcn_is_zero(cl))?T(1):feg_coef.cnl[ik];

					atomic_coef.feg.cl[ik] = cl;
					atomic_coef.feg.cnl[ik] = cnl;

					atomic_coef.fxg.cl[ik] = c_2Pi2a0<T>*cl/cnl;
					atomic_coef.fxg.cnl[ik] = cnl;

					atomic_coef.pr.cl[ik] = c_2Pi2a0<T>*c_pi2<T>*cl/pow(cnl, T(2.5));
					atomic_coef.pr.cnl[ik] = c_2pi<T>/::sqrt(cnl);

					atomic_coef.vr.cl[ik] = c_pot_factor<T>*c_pi2<T>*cl/pow(cnl, T(1.5));
					atomic_coef.vr.cnl[ik] = c_2pi<T>/::sqrt(cnl);

					atomic_coef.vzp.cl[ik] = 2.0*c_pot_factor<T>*c_pi2<T>*cl/pow(cnl, T(1.5));
					atomic_coef.vzp.cnl[ik] = c_2pi<T>/::sqrt(cnl);
				}

				return atomic_coef;
			}

			// 10: Peng et al. parameterization for ions - 5 Gaussians - [0, 4]
			Atomic_Coef_cpu<T> atomic_coef_peng_ion_0_4(const dt_int32& Z, const dt_int32& charge)		
			{
				auto feg_coef = load_feg_peng_ion_0_4(Z, charge);

				if (fcn_is_zero(feg_coef.cl[0]))
				{
					return atomic_coef_peng_neutral_0_4(Z);
				}

				Atomic_Coef_cpu<T> atomic_coef(feg_coef.size()+1, T(0));
				atomic_coef.atomic_pot_parm_typ = eappt_peng_ion_0_4;
				atomic_coef.Z = Z;
				atomic_coef.charge = charge;

				for(auto ik = 0; ik < 5; ik++)
				{
					auto cl = feg_coef.cl[ik];
					auto cnl = (fcn_is_zero(cl))?T(0):feg_coef.cnl[ik]/T(4);

					if (fcn_is_nzero(cl, cnl))
					{
						atomic_coef.feg.cl[ik] = cl;
						atomic_coef.feg.cnl[ik] = cnl;

						atomic_coef.fxg.cl[ik] = c_2Pi2a0<T>*cl;
						atomic_coef.fxg.cnl[ik] = cnl;

						atomic_coef.pr.cl[ik] = T(0.5)*c_2Pi2a0<T>*cl*pow(c_pi<T>/cnl, T(3.5));
						atomic_coef.pr.cnl[ik] = c_pi2<T>/cnl;

						atomic_coef.vr.cl[ik] = c_pot_factor<T>*cl*pow(c_pi<T>/cnl, T(1.5));
						atomic_coef.vr.cnl[ik] = c_pi2<T>/cnl;

						atomic_coef.vzp.cl[ik] = c_pot_factor<T>*c_pi<T>*cl/cnl;
						atomic_coef.vzp.cnl[ik] = c_pi2<T>/cnl;
					}
				}

				// 1/g2
				if (fcn_is_nzero(feg_coef.cl[0]))
				{
					dt_int32 ik = 5;
					T r0 = 2.0;
					auto cl = charge/c_2Pi2a0<T>;
					auto cnl = pow(c_2pi<T>*r0, -2);

					atomic_coef.feg.cl[ik] = cl;
					atomic_coef.feg.cnl[ik] = cnl;

					atomic_coef.fxg.cl[ik] = c_2Pi2a0<T>*cl;
					atomic_coef.fxg.cnl[ik] = cnl;

					atomic_coef.pr.cl[ik] = c_2Pi2a0<T>*c_pi<T>*cl*cnl;
					atomic_coef.pr.cnl[ik] = c_2pi<T>*::sqrt(cnl);

					atomic_coef.vr.cl[ik] = c_pot_factor<T>*c_pi<T>*cl;
					atomic_coef.vr.cnl[ik] = c_2pi<T>*::sqrt(cnl);

					atomic_coef.vzp.cl[ik] = c_pot_factor<T>*c_2pi<T>*cl;
					atomic_coef.vzp.cnl[ik] = c_2pi<T>*::sqrt(cnl);
				}

				return atomic_coef;
			}
		};

	}
#endif