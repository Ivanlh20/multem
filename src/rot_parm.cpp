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

#ifndef ROT_PARM_H
	#define ROT_PARM_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include "const_enum_mt.cuh"
	#include "math.cuh"
	#include "r_3d.cuh"

	namespace mt
	{
		/*********************** rotation parameters *************************/
		template <class T>
		class Rot_Parm
		{
			public:
				T theta;						// angle
				R_3d<T> u_0;					// unitary vector			
				eRot_Ctr_Typ ctr_type;			// erct_geometric_ctr, erct_user_def		
				R_3d<T> ctr_p;					// rotation point

				Rot_Parm(): theta(0), u_0(0, 0, 1), ctr_type(erct_geometric_ctr), ctr_p(1, 0, 0) {}
				
				template<class U>
				void assign(Rot_Parm<U>& rot_parm)
				{
					if (this !=  &rot_parm)
					{
						theta = rot_parm.theta;
						u_0 = rot_parm.u_0;
						ctr_type = rot_parm.ctr_type;
						ctr_p = rot_parm.ctr_p;
					}
				}

				void set_in_data(const T& theta, const R_3d<T>& u_0, const eRot_Ctr_Typ& ctr_type, const R_3d<T>& ctr_p)
				{
					this->theta = theta;
					this->u_0 = u_0;
					this->ctr_type = ctr_type;
					this->ctr_p = ctr_p;

					set_dep_var();
				}

				void set_dep_var()
				{
					u_0.normalize();
				}

				template <class U> 
				Rot_Parm<T>& operator=(Rot_Parm<U>& rot_parm)
				{
					assign(rot_parm);
					return *this;
				}

				dt_bool is_rot_ctr_none() const
				{
					return mt::is_rot_ctr_none(center_type);
				}

				dt_bool is_rot_ctr_geometric_ctr() const
				{
					return mt::is_rot_ctr_geometric_ctr(center_type);
				}

				dt_bool is_rot_ctr_user_def() const
				{
					return mt::is_rot_ctr_user_def(center_type);
				}

				dt_bool is_rot_active() const
				{
					return fcn_is_nzero(theta);
				}
		};
	}

#endif