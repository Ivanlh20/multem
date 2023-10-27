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

#ifndef SPEC_SLIC_IN_PARM_H
	#define SPEC_SLIC_IN_PARM_H

	#include "const_enum_mt.cuh"
	#include "r_2d.h"
	#include "cgpu_vctr.cuh"

	namespace mt
	{
		template <class T>
		class Spec_Slic_In_Parm
		{
			public:
				using value_type = T;

				eSpec_Slic_Typ typ;				// esst_plns_proj = 1, esst_dz_proj = 2,  esst_plns_sub = 3, esst_dz_sub = 4, esst_user_def = 5, esst_auto = 6
				T sli_thick;					// slice thickness
				eSpec_Slic_Sel_Typ sel_typ;		// essso_tag = 1, essso_z = 2
				dt_int32 sel_tag;				// tag
				dt_int32 sel_Z;					// atomic number
				R_2d<T> sel_z_lim;				// [z_0, z_e]
				Vctr_cpu<T> z_plns;				// z planes positions

				Spec_Slic_In_Parm(): typ( esst_plns_proj), sli_thick(0), sel_typ(0), 
				sel_tag(0), sel_Z(0), sel_z_lim{0, 0} {}

				template <class U> 
				void assign(Spec_Slic_In_Parm<U>& spec_slic_in_parm)
				{
					if (this !=  &spec_slic_in_parm)
					{
						typ = spec_slic_in_parm.typ;
						sli_thick = T(spec_slic_in_parm.sli_thick);
						sel_typ = spec_slic_in_parm.sel_typ;
						sel_tag = spec_slic_in_parm.sel_tag;
						sel_Z = spec_slic_in_parm.sel_Z;
						sel_z_lim = spec_slic_in_parm.sel_z_lim;
						z_plns = spec_slic_in_parm.z_plns;
					}
				}

				void set_in_data(const eSpec_Slic_Typ& typ, const T& sli_thick, const eSpec_Slic_Sel_Typ& sel_typ, 
					const dt_int32& sel_tag, const dt_int32& sel_Z, const R_2d<T>& sel_z_lim, const Vctr_cpu<T>& z_plns)
				{
					this->typ = typ;
					this->sli_thick = sli_thick;
					this->sel_typ = sel_typ;
					this->sel_tag = sel_tag;
					this->sel_Z = sel_Z;
					this->sel_z_lim = sel_z_lim;
					this->z_plns = z_plns;

					set_dep_var();
				}

				void set_dep_var()
				{
					sel_Z = fcn_max(sel_Z, 0);

					if (!is_spec_slic_by_user_def())
					{
						z_plns.clear_shrink_to_fit();
					}
					else
					{
						if (is_spec_slic_sel_typ_by_z() && (sel_z_lim.y<sel_z_lim.x))
						{
							sel_typ = essso_tag;
						}
					}
				}

				template <class U> 
				Spec_Slic_In_Parm<T>& operator=(Spec_Slic_In_Parm<U>& spec_slic_in_parm)
				{
					assign(spec_slic_in_parm);
					return *this;
				}

				void clear()
				{
					typ = esst_plns_proj;
					sli_thick = 0;
					el_typ = 0;
					sel_tag = 0;
					sel_Z = 0;
					sel_z_lim = 0;
				}

				dt_bool is_spec_slic_by_plns_proj()
				{
					return mt::is_spec_slic_by_plns_proj(typ);
				}

				dt_bool is_spec_slic_by_dz_proj()
				{
					return  mt::is_spec_slic_by_dz_proj(typ);
				}		

				dt_bool is_spec_slic_by_plns_sub()
				{
					return  mt::is_spec_slic_by_plns_sub(typ);
				}

				dt_bool is_spec_slic_by_dz_sub()
				{
					return  mt::is_spec_slic_by_dz_sub(typ);
				}

				dt_bool is_spec_slic_by_user_def()
				{
					return  mt::is_spec_slic_by_user_def(typ);
				}

				dt_bool is_spec_slic_by_auto()
				{
					return  mt::is_spec_slic_by_auto(typ);
				}

				dt_bool is_spec_slic_by_planes()
				{
					return  mt::is_spec_slic_by_planes(typ);
				}

				dt_bool is_spec_slic_sel_typ_by_tag()
				{
					return mt::is_spec_slic_sel_typ_by_tag(sel_typ);
				}

				dt_bool is_spec_slic_sel_typ_by_z()
				{
					return mt::is_spec_slic_sel_typ_by_z(sel_typ);
				}
		};
	
		template <class T>
		using Vctr_Spec_Slic_In_Parm = Vctr_cpu<Spec_Slic_In_Parm<T>>;
	}

#endif
