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

#ifndef TYPE_TRAITS_MT_H
	#define TYPE_TRAITS_MT_H

	#ifdef _MSC_VER
		#pragma once
	#endif

	#include <type_traits>

	#include "const_enum_mt.cuh"
	#include "type_traits_gen.cuh"

	namespace mt
	{
		/***************************************************************************************/
		template <ePot_Parm_Typ pot_parm_typ, class U=void>
		using enable_if_ePPT_doyle_0_4 = typename std::enable_if<pot_parm_typ==ePPT_doyle_0_4, U>::type;

		template <ePot_Parm_Typ pot_parm_typ, class U=void>
		using enable_if_ePPT_peng_0_4 = typename std::enable_if<pot_parm_typ==ePPT_peng_0_4, U>::type;

		template <ePot_Parm_Typ pot_parm_typ, class U=void>
		using enable_if_ePPT_peng_0_12 = typename std::enable_if<pot_parm_typ==ePPT_peng_0_12, U>::type;

		template <ePot_Parm_Typ pot_parm_typ, class U=void>
		using enable_if_ePPT_kirkland_0_12 = typename std::enable_if<pot_parm_typ==ePPT_kirkland_0_12, U>::type;

		template <ePot_Parm_Typ pot_parm_typ, class U=void>
		using enable_if_ePPT_weickenmeier_0_12 = typename std::enable_if<pot_parm_typ==ePPT_weickenmeier_0_12, U>::type;

		template <ePot_Parm_Typ pot_parm_typ, class U=void>
		using enable_if_ePPT_lobato_0_12 = typename std::enable_if<pot_parm_typ==ePPT_lobato_0_12, U>::type;

		template <ePot_Parm_Typ pot_parm_typ, class U=void>
		using enable_if_ePPT_peng_ion_0_4 = typename std::enable_if<pot_parm_typ==ePPT_peng_ion_0_4, U>::type;

		/****************************** types ********************************/
		template <dt_int32 simulation_type, class U=void>
		using enable_if_STEM = typename std::enable_if<simulation_type == eTEMST_STEM, U>::type;

		template <dt_int32 simulation_type, class U=void>
		using enable_if_ISTEM = typename std::enable_if<simulation_type == eTEMST_ISTEM, U>::type;

		template <dt_int32 simulation_type, class U=void>
		using enable_if_CBED = typename std::enable_if<simulation_type == eTEMST_CBED, U>::type;

		template <dt_int32 simulation_type, class U=void>
		using enable_if_CBEI = typename std::enable_if<simulation_type == eTEMST_CBEI, U>::type;

		template <dt_int32 simulation_type, class U=void>
		using enable_if_ED = typename std::enable_if<simulation_type == eTEMST_ED, U>::type;

		template <dt_int32 simulation_type, class U=void>
		using enable_if_HRTEM = typename std::enable_if<simulation_type == eTEMST_HRTEM, U>::type;

		template <dt_int32 simulation_type, class U=void>
		using enable_if_PED = typename std::enable_if<simulation_type == eTEMST_PED, U>::type;

		template <dt_int32 simulation_type, class U=void>
		using enable_if_HCTEM = typename std::enable_if<simulation_type == eTEMST_HCTEM, U>::type;

		template <dt_int32 simulation_type, class U=void>
		using enable_if_EWFS = typename std::enable_if<simulation_type == eTEMST_EWFS, U>::type;

		template <dt_int32 simulation_type, class U=void>
		using enable_if_EWRS = typename std::enable_if<simulation_type == eTEMST_EWRS, U>::type;

		template <dt_int32 simulation_type, class U=void>
		using enable_if_EELS = typename std::enable_if<simulation_type == eTEMST_STEM_EELS, U>::type;

		template <dt_int32 simulation_type, class U=void>
		using enable_if_EFTEM = typename std::enable_if<simulation_type == eTEMST_EFTEMRS, U>::type;

		template <dt_int32 simulation_type, class U=void>
		using enable_if_ProbeFS = typename std::enable_if<simulation_type == eTEMST_IWFS, U>::type;

		template <dt_int32 simulation_type, class U=void>
		using enable_if_ProbeRS = typename std::enable_if<simulation_type == eTEMST_IWRS, U>::type;

	}

#endif
