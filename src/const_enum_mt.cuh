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

#ifndef CONST_ENUM_MT_H
	#define CONST_ENUM_MT_H

	#ifdef _MSC_VER
		#pragma once
	#endif 

	#include "macros.cuh"
	#include "const_enum.cuh"
	
	/* constant -  enumeration definitions */
	namespace mt
	{
		/**************************** in atoms ****************************/
		enum eIn_Atoms
		{
			eIA_yes = 1, eIA_no = 2
		};

		/************************** slice memory type ************************/
		enum eSlice_Memory_Typ
		{
			eSMT_Transmission = 1, eSMT_Potential = 2, eSMT_none = 3
		};

		/************************ microscope effects *************************/
		enum eIllumination_Model
		{
			eIM_Coherent = 1, eIM_Partial_Coherent = 2, eIM_Trans_Cross_Coef = 3, eIM_Full_Integration = 4, eIM_none = 5
		};

		/************************ spatial and temporal ***********************/
		enum eTemporal_Spatial_Incoh
		{
			eTSI_Temporal_Spatial = 1, eTSI_Temporal = 2, eTSI_Spatial = 3, eTSI_none = 4
		};

		/************************** lens variable type ***********************/
		enum eLens_Var_Typ
		{
			eLVT_off = 0, eLVT_m = 1, eLVT_f = 2, eLVT_Cs3 = 3, eLVT_Cs5 = 4, 
			eLVT_mfa2 = 5, eLVT_afa2 = 6, eLVT_mfa3 = 7, eLVT_afa3 = 8, 
			eLVT_inner_aper_ang = 9, eLVT_outer_aper_ang = 10
		};

		/************************* simulation type ***************************/
		enum eTEM_Sim_Typ
		{
			eTEMST_STEM = 11, eTEMST_ISTEM = 12, 
			eTEMST_CBED = 21, eTEMST_CBEI = 22, 
			eTEMST_ED = 31, eTEMST_HRTEM = 32, 
			eTEMST_PED = 41, eTEMST_HCTEM = 42, 
			eTEMST_EWFS = 51, eTEMST_EWRS = 52, 
			eTEMST_STEM_EELS = 61, eTEMST_ISTEM_EELS = 62, 
			eTEMST_EFTEMFS = 71, eTEMST_EFTEMRS = 72, 
			eTEMST_IWFS = 81, eTEMST_IWRS = 82, 
			eTEMST_PPFS = 91, eTEMST_PPRS = 92, 			// projected potential
			eTEMST_TFFS = 101, eTEMST_TFRS = 102, 			// transmission function
			eTEMST_PropFS = 111, eTEMST_PropRS = 112		// fcn_propagate
		};

		/********************* simulation data output ************************/
		enum eTEM_Output_Typ
		{
			eTEMOT_image_tot_coh = 1, eTEMOT_image_tot = 2, 
			eTEMOT_m2psi_tot_coh = 3, eTEMOT_m2psi_tot = 4, 
			eTEMOT_m2psi_tot_psi_coh = 5, eTEMOT_psi_coh = 6, 
			eTEMOT_psi_0 = 7, eTEMOT_V = 8, eTEMOT_trans = 9
		};

		/************** electron specimen interaction model ******************/
		enum eElec_Spec_Int_Model
		{
			eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
		};

		/*************************** phonon model ****************************/
		enum ePhonon_Model
		{
			ePM_Still_Atom = 1, ePM_Absorptive_Model = 2, ePM_Frozen_Phonon = 3, ePM_user_def = 4
		};

		/************************ phonon model output ************************/
		enum ePhonon_Model_Output
		{
			ePMO_Total = 1, ePMO_Coherent = 2
		};

		/******************** potential slicing Type *************************/
		enum ePot_Sli_Typ
		{
			ePST_Planes = 1, ePST_dz_Proj = 2, ePST_dz_Sub = 3, ePST_Auto = 4
		};

		/********************** feg parameterization *************************/
		// 1: doyle and Turner parameterization - 4 Gaussians - [0, 4
		// 2: Peng et al. parameterization - 5 Gaussians - [0, 4]
		// 3: Peng et al. parameterization - 5 Gaussians - [0, 12]
		// 4: Kirkland parameterization - 3 lorentzian + 3 Gaussians - [0, 12]
		// 5: Weickenmeier and H.Kohl - a*(1-exp(-bg^2)/g^2 - [0, 12]
		// 6: Lobato parameterization - 5 Hydrogen feg - [0, 12]
		// 10: Peng et al. parameterization ion - 5 Gaussians - [0, 4]

		enum ePot_Parm_Typ
		{
			ePPT_none = 0, ePPT_doyle_0_4 = 1, ePPT_peng_0_4 = 2, ePPT_peng_0_12 = 3, 
			ePPT_kirkland_0_12 = 4, ePPT_weickenmeier_0_12 = 5, ePPT_lobato_0_12 = 6, 
			ePPT_peng_ion_0_4 = 10
		};

		/*********************** incident wave type **************************/
		enum eIncident_Wave_Typ
		{
			eIWT_Plane_Wave = 1, eIWT_Convergent_Wave = 2, eIWT_user_def_Wave = 3, eIWT_Auto = 4
		};

		/************************ amorphous layer type ***********************/
		enum eSpec_Lay_Typ
		{
			eSLT_none = 0, eSLT_Top = 1, eSLT_Bottom = 2, eSLT_Middle = 3, eSLT_user_def = 4
		};

		/************************* defocus plane type ************************/
		enum eZero_Defocus_Typ
		{
			eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_user_def = 4
		};

		/************************** thickness Type ***************************/
		enum eThick_Typ
		{
			eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
		};

		/**************************** scan type ******************************/
		enum eScan_Typ
		{
			eST_Line = 1, eST_Area = 2, eST_user_def = 3
		};

		/************************* detector type *****************************/
		enum eDetector_Typ
		{
			edt_Circular = 1, edt_Radial = 2, edt_Matrix = 3
		};

		/************************ channelling type **************************/
		// 1: single channelling
		// 2: mixed channelling
		// 3: double channelling
		enum eChan_Typ
		{
			eCT_Single_Chan = 1, eCT_Mixed_Chan = 2, eCT_Double_Chan = 3
		};

	}
	
	/* constant -  enumeration comparison */
	namespace mt
	{
		/***************************************************************************************/
		inline
		dt_bool is_spec_lay_top(const eSpec_Lay_Typ &type)
		{
			return type == mt::eSLT_Top;
		}		
		
		inline
		dt_bool is_spec_lay_bottom(const eSpec_Lay_Typ &type)
		{
			return type == mt::eSLT_Bottom;
		}
		
		inline
		dt_bool is_spec_lay_middle(const eSpec_Lay_Typ &type)
		{
			return type == mt::eSLT_Middle;
		}
		
		inline
		dt_bool is_spec_lay_user_def(const eSpec_Lay_Typ &type)
		{
			return type == mt::eSLT_user_def;
		}

		/***************************************************************************************/
		inline
		dt_bool is_multislice(const eElec_Spec_Int_Model &int_model)
		{
			return int_model == mt::eESIM_Multislice;
		}

		inline
		dt_bool is_phase_object(const eElec_Spec_Int_Model &int_model)
		{
			return int_model == mt::eESIM_Phase_Object;
		}

		inline
		dt_bool is_weak_phase_object(const eElec_Spec_Int_Model &int_model)
		{
			return int_model == mt::eESIM_Weak_Phase_Object;
		}

		/***************************************************************************************/
		inline
		dt_bool is_still_atom(const ePhonon_Model &pn_model)
		{
			return pn_model == ePM_Still_Atom;
		}

		inline
		dt_bool is_absorptive_model(const ePhonon_Model &pn_model)
		{
			return pn_model == ePM_Absorptive_Model;
		}

		inline
		dt_bool is_frozen_phonon(const ePhonon_Model &pn_model)
		{
			return pn_model == ePM_Frozen_Phonon;
		}

		inline
		dt_bool is_frozen_phonon_float32_conf(const ePhonon_Model &pn_model, const dt_bool &pn_float32_conf)
		{
			return is_frozen_phonon(pn_model) && pn_float32_conf;
		}

		/***************************************************************************************/
		inline
		dt_bool is_whole_spec(const eThick_Typ &thick_type)
		{
			return thick_type == eTT_Whole_Spec;
		}

		inline
		dt_bool is_through_slices(const eThick_Typ &thick_type)
		{
			return thick_type == eTT_Through_Slices;
		}

		inline
		dt_bool is_through_thick(const eThick_Typ &thick_type)
		{
			return thick_type == eTT_Through_Thick;
		}

		/***************************************************************************************/
		inline
		dt_bool is_slicing_by_planes(const eElec_Spec_Int_Model &int_model, const ePot_Sli_Typ &pot_slic)
		{
			return mt::is_multislice(int_model) && (pot_slic == mt::ePST_Planes);
		}

		inline
		dt_bool is_slicing_by_dz(const eElec_Spec_Int_Model &int_model, const ePot_Sli_Typ &pot_slic)
		{
			return mt::is_multislice(int_model) && (pot_slic == mt::ePST_dz_Proj);
		}

		inline
		dt_bool is_subslicing(const eElec_Spec_Int_Model &int_model, const ePot_Sli_Typ &pot_slic)
		{
			return mt::is_multislice(int_model) && (pot_slic == mt::ePST_dz_Sub);
		}

		inline
		dt_bool is_subslicing_whole_spec(const eElec_Spec_Int_Model &int_model, const ePot_Sli_Typ &pot_slic, const eThick_Typ &thick_type)
		{
			return mt::is_subslicing(int_model, pot_slic) && is_whole_spec(thick_type);
		}

		/***************************************************************************************/
		inline
		dt_bool is_plane_wave(eIncident_Wave_Typ iw_type)
		{
			return iw_type == eIWT_Plane_Wave;
		}

		inline
		dt_bool is_convergent_wave(eIncident_Wave_Typ iw_type)
		{
			return iw_type == eIWT_Convergent_Wave;
		}

		inline
		dt_bool is_user_define_wave(eIncident_Wave_Typ iw_type)
		{
			return iw_type == eIWT_user_def_Wave;
		}

		/***************************************************************************************/
		inline
		dt_bool is_STEM(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_STEM;
		}

		inline
		dt_bool is_ISTEM(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_ISTEM;
		}

		inline
		dt_bool is_CBED(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_CBED;
		}

		inline
		dt_bool is_CBEI(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_CBEI;
		}

		inline
		dt_bool is_ED(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_ED;
		}

		inline
		dt_bool is_HRTEM(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_HRTEM;
		}

		inline
		dt_bool is_PED(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_PED;
		}

		inline
		dt_bool is_HCTEM(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_HCTEM;
		}

		inline
		dt_bool is_EWFS(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_EWFS;
		}

		inline
		dt_bool is_EWRS(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_EWRS;
		}

		inline
		dt_bool is_EWFS_SC(const eTEM_Sim_Typ &sim_type, const ePhonon_Model &pn_model, const dt_bool &pn_float32_conf)
		{
			return is_EWFS(sim_type) && (!is_frozen_phonon(pn_model) || is_frozen_phonon_float32_conf(pn_model, pn_float32_conf));
		}

		inline
		dt_bool is_EWRS_SC(const eTEM_Sim_Typ &sim_type, const ePhonon_Model &pn_model, const dt_bool &pn_float32_conf)
		{
			return is_EWRS(sim_type) && (!is_frozen_phonon(pn_model) || is_frozen_phonon_float32_conf(pn_model, pn_float32_conf));
		}

		inline
		dt_bool is_STEM_EELS(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_STEM_EELS;
		}

		inline
		dt_bool is_ISTEM_EELS(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_ISTEM_EELS;
		}

		inline
		dt_bool is_STEM_ISTEM_EELS(const eTEM_Sim_Typ &sim_type)
		{
			return is_STEM_EELS(sim_type) || is_ISTEM_EELS(sim_type);
		}

		inline
		dt_bool is_EFTEMFS(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_EFTEMFS;
		}

		inline
		dt_bool is_EFTEMRS(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_EFTEMRS;
		}

		inline
		dt_bool is_EFTEM(const eTEM_Sim_Typ &sim_type)
		{
			return is_EFTEMFS(sim_type) || is_EFTEMRS(sim_type);
		}

		inline
		dt_bool is_IWFS(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_IWFS;
		}

		inline
		dt_bool is_IWRS(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_IWRS;
		}

		inline
		dt_bool is_PPFS(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_PPFS;
		}

		inline
		dt_bool is_PPRS(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_PPRS;
		}

		inline
		dt_bool is_TFFS(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_TFFS;
		}

		inline
		dt_bool is_TFRS(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_TFRS;
		}

		inline
		dt_bool is_PropFS(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_PropFS;
		}

		inline
		dt_bool is_PropRS(const eTEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eTEMST_PropRS;
		}

		inline
		dt_bool is_STEM_ISTEM(const eTEM_Sim_Typ &sim_type)
		{
			return is_STEM(sim_type) || is_ISTEM(sim_type);
		}

		inline
		dt_bool is_CBED_CBEI(const eTEM_Sim_Typ &sim_type)
		{
			return is_CBED(sim_type) || is_CBEI(sim_type);
		}

		inline
		dt_bool is_ED_HRTEM(const eTEM_Sim_Typ &sim_type)
		{
			return is_ED(sim_type) || is_HRTEM(sim_type);
		}

		inline
		dt_bool is_PED_HCTEM(const eTEM_Sim_Typ &sim_type)
		{
			return is_PED(sim_type) || is_HCTEM(sim_type);
		}

		inline
		dt_bool is_EWFS_EWRS(const eTEM_Sim_Typ &sim_type)
		{
			return is_EWFS(sim_type) || is_EWRS(sim_type);
		}

		inline
		dt_bool is_EWFS_EWRS_SC(const eTEM_Sim_Typ &sim_type, const ePhonon_Model &pn_model, const dt_bool &pn_float32_conf)
		{
			return is_EWFS_SC(sim_type, pn_model, pn_float32_conf) || is_EWRS_SC(sim_type, pn_model, pn_float32_conf);
		}

		inline
		dt_bool is_EWFS_convergent_wave(const eTEM_Sim_Typ &sim_type, const eIncident_Wave_Typ &iw_type)
		{
			return is_EWFS(sim_type) && is_convergent_wave(iw_type);
		}

		inline
		dt_bool is_EWRS_convergent_wave(const eTEM_Sim_Typ &sim_type, const eIncident_Wave_Typ &iw_type)
		{
			return is_EWRS(sim_type) && is_convergent_wave(iw_type);
		}

		inline
		dt_bool is_EW_convergent_wave(const eTEM_Sim_Typ &sim_type, const eIncident_Wave_Typ &iw_type)
		{
			return is_EWFS_EWRS(sim_type) && is_convergent_wave(iw_type);
		}

		inline
		dt_bool is_EELS_EFTEM(const eTEM_Sim_Typ &sim_type)
		{
			return is_STEM_ISTEM_EELS(sim_type) || is_EFTEM(sim_type);
		}

		inline
		dt_bool is_IWFS_IWRS(const eTEM_Sim_Typ &sim_type)
		{
			return is_IWFS(sim_type) || is_IWRS(sim_type);
		}

		inline
		dt_bool is_PPFS_PPRS(const eTEM_Sim_Typ &sim_type)
		{
			return is_PPFS(sim_type) || is_PPRS(sim_type);
		}

		inline
		dt_bool is_TFFS_TFRS(const eTEM_Sim_Typ &sim_type)
		{
			return is_TFFS(sim_type) || is_TFRS(sim_type);
		}

		inline
		dt_bool is_PropFS_PropRS(const eTEM_Sim_Typ &sim_type)
		{
			return is_PropFS(sim_type) || is_PropRS(sim_type);
		}

		inline
		dt_bool is_grid_FS(const eTEM_Sim_Typ &sim_type)
		{
			auto bb = is_CBED(sim_type) || is_ED(sim_type) || is_PED(sim_type) || is_EWFS(sim_type);
			bb = bb || is_EFTEMFS(sim_type) || is_IWFS(sim_type) || is_PPFS(sim_type) || is_TFFS(sim_type);
			return bb;
		}

		inline
		dt_bool is_grid_RS(const eTEM_Sim_Typ &sim_type)
		{
			return !is_grid_FS(sim_type);
		}

		inline
		dt_bool is_simulation_type_FS(const eTEM_Sim_Typ &sim_type)
		{
			auto bb = is_STEM(sim_type) || is_CBED(sim_type) || is_ED(sim_type);
			bb = bb || is_PED(sim_type) || is_EWFS(sim_type) || is_STEM_EELS(sim_type);
			bb = bb || is_EFTEMFS(sim_type) || is_IWFS(sim_type) || is_PPFS(sim_type) || is_TFFS(sim_type);
			return bb;
		}

		inline
		dt_bool is_simulation_type_RS(const eTEM_Sim_Typ &sim_type)
		{
			return !is_simulation_type_FS(sim_type);
		}

		inline
		dt_bool is_specimen_required(const eTEM_Sim_Typ &sim_type)
		{
			return !(is_IWFS_IWRS(sim_type) || is_PropFS_PropRS(sim_type));
		}

		inline
		dt_bool is_ISTEM_CBEI_HRTEM_HCTEM_EFTEMRS(const eTEM_Sim_Typ &sim_type)
		{
			return is_ISTEM(sim_type) || is_CBEI(sim_type) || is_HRTEM(sim_type) || is_HCTEM(sim_type) || is_EFTEMRS(sim_type);
		}

		inline
		dt_bool is_CBED_ED_EWFS_PED(const eTEM_Sim_Typ &sim_type)
		{
			return is_CBED(sim_type) || is_ED(sim_type) || is_EWFS(sim_type) || is_PED(sim_type);
		}

		inline
		dt_bool is_obj_lens_temp_spat(const eTEM_Sim_Typ &sim_type)
		{
			return is_ISTEM(sim_type) || is_CBEI(sim_type) || is_HRTEM(sim_type) || is_HCTEM(sim_type) || is_EFTEMRS(sim_type);
		}

		inline
		dt_bool is_cond_lens_temp_spat(const eTEM_Sim_Typ &sim_type)
		{
			return is_STEM_ISTEM(sim_type) || is_CBED_CBEI(sim_type) || is_STEM_ISTEM_EELS(sim_type) || is_EFTEM(sim_type);
		}

		inline
		eSpace get_simulation_space(const eTEM_Sim_Typ &sim_type)
		{
			return (is_simulation_type_FS(sim_type) || is_ISTEM_CBEI_HRTEM_HCTEM_EFTEMRS(sim_type))?mt::eS_Reciprocal:mt::eS_Real;
		}

		inline
		dt_bool is_scanning(const eTEM_Sim_Typ &sim_type)
		{
			return is_STEM_ISTEM(sim_type) || is_STEM_ISTEM_EELS(sim_type);
		}

		inline
		eIncident_Wave_Typ validate_incident_wave_type(const eTEM_Sim_Typ &sim_type, eIncident_Wave_Typ iw_type)
		{
			if (iw_type == eIWT_Auto)
			{
				auto bb = is_scanning(sim_type) || is_CBED_CBEI(sim_type);
				bb = bb || ((is_EWFS_EWRS(sim_type) || is_EFTEM(sim_type) || is_IWFS_IWRS(sim_type)) && is_convergent_wave(iw_type));
				iw_type = (bb)?mt::eIWT_Convergent_Wave:mt::eIWT_Plane_Wave;
			}

			return iw_type;
		}

	}
#endif