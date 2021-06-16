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
			esmt_none = 0, esmt_transmission = 1, esmt_potential = 2
		};

		/************************ microscope effects *************************/
		enum eIllumination_Model
		{
			eim_none = 0, eim_coherent = 1, eim_partial_coherent = 2, eim_trans_cross_coef = 3, eim_full_integration = 4
		};

		/************************ spatial and temporal ***********************/
		enum eTemporal_Spatial_Incoh
		{
			etsi_none = 0, etsi_temporal_spatial = 1, etsi_temporal = 2, etst_spatial = 3
		};

		/************************** lens variable type ***********************/
		enum eLens_Var_Typ
		{
			eLVT_off = 0, eLVT_m = 1, eLVT_f = 2, eLVT_Cs3 = 3, eLVT_Cs5 = 4, 
			eLVT_mfa2 = 5, eLVT_afa2 = 6, eLVT_mfa3 = 7, eLVT_afa3 = 8, 
			eLVT_inner_aper_ang = 9, eLVT_outer_aper_ang = 10
		};

		/************************* simulation type ***************************/
		enum eEM_Sim_Typ
		{
			eemst_stem = 11, eemst_istem = 12, 
			eemst_cbed = 21, eemst_cbei = 22, 
			eemst_ed = 31, eemst_hrtem = 32, 
			eemst_ped = 41, eemst_hctem = 42, 
			eemst_ewfs = 51, eemst_ewrs = 52, 
			eemst_stem_eels = 61, eemst_istem_eels = 62, 
			eemst_eftemfs = 71, eemst_eftemrs = 72, 
			eemst_iwfs = 81, eemst_iwrs = 82, 
			eemst_ppfs = 91, eemst_pprs = 92, 			// projected potential
			eemst_tffs = 101, eemst_tfrs = 102, 			// transmission function
			eemst_propfs = 111, eemst_proprs = 112		// fcn_propagate
		};

		/********************* simulation data output ************************/
		enum eEM_Output_Typ
		{
			eemot_image_tot_coh = 1, eemot_image_tot = 2, 
			eemot_m2psi_tot_coh = 3, eemot_m2psi_tot = 4, 
			eemot_m2psi_tot_psi_coh = 5, eemot_psi_coh = 6, 
			eemot_psi_0 = 7, eemot_V = 8, eemot_trans = 9
		};

		/************** electron specimen interaction model ******************/
		enum eElec_Spec_Int_Model
		{
			eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
		};

		/*************************** phonon model ****************************/
		enum ePhonon_Model
		{
			epm_still_atom = 1, epm_absorptive_model = 2, epm_frozen_phonon = 3, epm_user_def = 4
		};

		/************************ phonon model output ************************/
		enum ePhonon_Model_Output
		{
			epmo_total = 1, epmo_coherent = 2, epmo_total_coherent = 3
		};

		/******************** potential slicing Type *************************/
		enum ePot_Sli_Typ
		{
			epst_planes = 1, epst_dz_Proj = 2, epst_dz_Sub = 3, epst_auto = 4
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
			eppt_none = 0, eppt_doyle_0_4 = 1, eppt_peng_0_4 = 2, eppt_peng_0_12 = 3, 
			eppt_kirkland_0_12 = 4, eppt_weickenmeier_0_12 = 5, eppt_lobato_0_12 = 6, 
			eppt_peng_ion_0_4 = 10
		};

		/*********************** incident wave type **************************/
		enum eIncident_Wave_Typ
		{
			eiwt_plane_wave = 1, eiwt_convergent_wave = 2, eiwt_user_def_Wave = 3, eiwt_auto = 4
		};

		/************************ amorphous layer type ***********************/
		enum eSpec_Lay_Typ
		{
			eslt_none = 0, eslt_top = 1, eslt_bottom = 2, eslt_middle = 3, eslt_user_def = 4
		};

		/************************* defocus plane type ************************/
		enum eZero_Defocus_Typ
		{
			ezdt_first = 1, ezdt_middle = 2, ezdt_last = 3, ezdt_user_def = 4
		};

		/************************** thickness Type ***************************/
		enum eThick_Typ
		{
			ett_whole_spec = 1, ett_through_thick = 2, ett_through_slices = 3
		};

		/**************************** scan type ******************************/
		enum eScan_Typ
		{
			est_line = 1, est_area = 2, eST_user_def = 3
		};

		/************************* detector type *****************************/
		enum eDetector_Typ
		{
			edt_circular = 1, edt_radial = 2, edt_matrix = 3
		};

		/************************ channelling type **************************/
		// 1: single channelling
		// 2: mixed channelling
		// 3: double channelling
		enum eChan_Typ
		{
			ect_single_chan = 1, ect_mixed_chan = 2, eCT_double_chan = 3
		};

	}
	
	/* constant -  enumeration comparison */
	namespace mt
	{
		/***************************************************************************************/
		inline
		dt_bool is_spec_lay_top(const eSpec_Lay_Typ &type)
		{
			return type == mt::eslt_top;
		}		
		
		inline
		dt_bool is_spec_lay_bottom(const eSpec_Lay_Typ &type)
		{
			return type == mt::eslt_bottom;
		}
		
		inline
		dt_bool is_spec_lay_middle(const eSpec_Lay_Typ &type)
		{
			return type == mt::eslt_middle;
		}
		
		inline
		dt_bool is_spec_lay_user_def(const eSpec_Lay_Typ &type)
		{
			return type == mt::eslt_user_def;
		}

		/***************************************************************************************/
		inline
		dt_bool is_multislice(const eElec_Spec_Int_Model &int_model)
		{
			return int_model == mt::eesim_multislice;
		}

		inline
		dt_bool is_phase_object(const eElec_Spec_Int_Model &int_model)
		{
			return int_model == mt::eesim_phase_object;
		}

		inline
		dt_bool is_weak_phase_object(const eElec_Spec_Int_Model &int_model)
		{
			return int_model == mt::eesim_weak_phase_object;
		}

		/***************************************************************************************/
		inline
		dt_bool is_still_atom(const ePhonon_Model &pn_model)
		{
			return pn_model == epm_still_atom;
		}

		inline
		dt_bool is_absorptive_model(const ePhonon_Model &pn_model)
		{
			return pn_model == epm_absorptive_model;
		}

		inline
		dt_bool is_frozen_phonon(const ePhonon_Model &pn_model)
		{
			return pn_model == epm_frozen_phonon;
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
			return thick_type == ett_whole_spec;
		}

		inline
		dt_bool is_through_slices(const eThick_Typ &thick_type)
		{
			return thick_type == ett_through_slices;
		}

		inline
		dt_bool is_through_thick(const eThick_Typ &thick_type)
		{
			return thick_type == ett_through_thick;
		}

		/***************************************************************************************/
		inline
		dt_bool is_slicing_by_planes(const eElec_Spec_Int_Model &int_model, const ePot_Sli_Typ &pot_slic)
		{
			return mt::is_multislice(int_model) && (pot_slic == mt::epst_planes);
		}

		inline
		dt_bool is_slicing_by_dz(const eElec_Spec_Int_Model &int_model, const ePot_Sli_Typ &pot_slic)
		{
			return mt::is_multislice(int_model) && (pot_slic == mt::epst_dz_Proj);
		}

		inline
		dt_bool is_subslicing(const eElec_Spec_Int_Model &int_model, const ePot_Sli_Typ &pot_slic)
		{
			return mt::is_multislice(int_model) && (pot_slic == mt::epst_dz_Sub);
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
			return iw_type == eiwt_plane_wave;
		}

		inline
		dt_bool is_convergent_wave(eIncident_Wave_Typ iw_type)
		{
			return iw_type == eiwt_convergent_wave;
		}

		inline
		dt_bool is_user_define_wave(eIncident_Wave_Typ iw_type)
		{
			return iw_type == eiwt_user_def_Wave;
		}

		/***************************************************************************************/
		inline
		dt_bool is_STEM(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_stem;
		}

		inline
		dt_bool is_ISTEM(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_istem;
		}

		inline
		dt_bool is_CBED(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_cbed;
		}

		inline
		dt_bool is_CBEI(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_cbei;
		}

		inline
		dt_bool is_ED(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_ed;
		}

		inline
		dt_bool is_HRTEM(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_hrtem;
		}

		inline
		dt_bool is_PED(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_ped;
		}

		inline
		dt_bool is_HCTEM(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_hctem;
		}

		inline
		dt_bool is_EWFS(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_ewfs;
		}

		inline
		dt_bool is_EWRS(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_ewrs;
		}

		inline
		dt_bool is_EWFS_SC(const eEM_Sim_Typ &sim_type, const ePhonon_Model &pn_model, const dt_bool &pn_float32_conf)
		{
			return is_EWFS(sim_type) && (!is_frozen_phonon(pn_model) || is_frozen_phonon_float32_conf(pn_model, pn_float32_conf));
		}

		inline
		dt_bool is_EWRS_SC(const eEM_Sim_Typ &sim_type, const ePhonon_Model &pn_model, const dt_bool &pn_float32_conf)
		{
			return is_EWRS(sim_type) && (!is_frozen_phonon(pn_model) || is_frozen_phonon_float32_conf(pn_model, pn_float32_conf));
		}

		inline
		dt_bool is_STEM_EELS(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_stem_eels;
		}

		inline
		dt_bool is_ISTEM_EELS(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_istem_eels;
		}

		inline
		dt_bool is_STEM_ISTEM_EELS(const eEM_Sim_Typ &sim_type)
		{
			return is_STEM_EELS(sim_type) || is_ISTEM_EELS(sim_type);
		}

		inline
		dt_bool is_EFTEMFS(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_eftemfs;
		}

		inline
		dt_bool is_EFTEMRS(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_eftemrs;
		}

		inline
		dt_bool is_EFTEM(const eEM_Sim_Typ &sim_type)
		{
			return is_EFTEMFS(sim_type) || is_EFTEMRS(sim_type);
		}

		inline
		dt_bool is_IWFS(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_iwfs;
		}

		inline
		dt_bool is_IWRS(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_iwrs;
		}

		inline
		dt_bool is_PPFS(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_ppfs;
		}

		inline
		dt_bool is_PPRS(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_pprs;
		}

		inline
		dt_bool is_TFFS(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_tffs;
		}

		inline
		dt_bool is_TFRS(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_tfrs;
		}

		inline
		dt_bool is_PropFS(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_propfs;
		}

		inline
		dt_bool is_PropRS(const eEM_Sim_Typ &sim_type)
		{
			return sim_type == mt::eemst_proprs;
		}

		inline
		dt_bool is_STEM_ISTEM(const eEM_Sim_Typ &sim_type)
		{
			return is_STEM(sim_type) || is_ISTEM(sim_type);
		}

		inline
		dt_bool is_CBED_CBEI(const eEM_Sim_Typ &sim_type)
		{
			return is_CBED(sim_type) || is_CBEI(sim_type);
		}

		inline
		dt_bool is_ED_HRTEM(const eEM_Sim_Typ &sim_type)
		{
			return is_ED(sim_type) || is_HRTEM(sim_type);
		}

		inline
		dt_bool is_PED_HCTEM(const eEM_Sim_Typ &sim_type)
		{
			return is_PED(sim_type) || is_HCTEM(sim_type);
		}

		inline
		dt_bool is_EWFS_EWRS(const eEM_Sim_Typ &sim_type)
		{
			return is_EWFS(sim_type) || is_EWRS(sim_type);
		}

		inline
		dt_bool is_EWFS_EWRS_SC(const eEM_Sim_Typ &sim_type, const ePhonon_Model &pn_model, const dt_bool &pn_float32_conf)
		{
			return is_EWFS_SC(sim_type, pn_model, pn_float32_conf) || is_EWRS_SC(sim_type, pn_model, pn_float32_conf);
		}

		inline
		dt_bool is_EWFS_convergent_wave(const eEM_Sim_Typ &sim_type, const eIncident_Wave_Typ &iw_type)
		{
			return is_EWFS(sim_type) && is_convergent_wave(iw_type);
		}

		inline
		dt_bool is_EWRS_convergent_wave(const eEM_Sim_Typ &sim_type, const eIncident_Wave_Typ &iw_type)
		{
			return is_EWRS(sim_type) && is_convergent_wave(iw_type);
		}

		inline
		dt_bool is_EW_convergent_wave(const eEM_Sim_Typ &sim_type, const eIncident_Wave_Typ &iw_type)
		{
			return is_EWFS_EWRS(sim_type) && is_convergent_wave(iw_type);
		}

		inline
		dt_bool is_EELS_EFTEM(const eEM_Sim_Typ &sim_type)
		{
			return is_STEM_ISTEM_EELS(sim_type) || is_EFTEM(sim_type);
		}

		inline
		dt_bool is_IWFS_IWRS(const eEM_Sim_Typ &sim_type)
		{
			return is_IWFS(sim_type) || is_IWRS(sim_type);
		}

		inline
		dt_bool is_PPFS_PPRS(const eEM_Sim_Typ &sim_type)
		{
			return is_PPFS(sim_type) || is_PPRS(sim_type);
		}

		inline
		dt_bool is_TFFS_TFRS(const eEM_Sim_Typ &sim_type)
		{
			return is_TFFS(sim_type) || is_TFRS(sim_type);
		}

		inline
		dt_bool is_PropFS_PropRS(const eEM_Sim_Typ &sim_type)
		{
			return is_PropFS(sim_type) || is_PropRS(sim_type);
		}

		inline
		dt_bool is_grid_FS(const eEM_Sim_Typ &sim_type)
		{
			auto bb = is_CBED(sim_type) || is_ED(sim_type) || is_PED(sim_type) || is_EWFS(sim_type);
			bb = bb || is_EFTEMFS(sim_type) || is_IWFS(sim_type) || is_PPFS(sim_type) || is_TFFS(sim_type);
			return bb;
		}

		inline
		dt_bool is_grid_RS(const eEM_Sim_Typ &sim_type)
		{
			return !is_grid_FS(sim_type);
		}

		inline
		dt_bool is_simulation_type_FS(const eEM_Sim_Typ &sim_type)
		{
			auto bb = is_STEM(sim_type) || is_CBED(sim_type) || is_ED(sim_type);
			bb = bb || is_PED(sim_type) || is_EWFS(sim_type) || is_STEM_EELS(sim_type);
			bb = bb || is_EFTEMFS(sim_type) || is_IWFS(sim_type) || is_PPFS(sim_type) || is_TFFS(sim_type);
			return bb;
		}

		inline
		dt_bool is_simulation_type_RS(const eEM_Sim_Typ &sim_type)
		{
			return !is_simulation_type_FS(sim_type);
		}

		inline
		dt_bool is_specimen_required(const eEM_Sim_Typ &sim_type)
		{
			return !(is_IWFS_IWRS(sim_type) || is_PropFS_PropRS(sim_type));
		}

		inline
		dt_bool is_ISTEM_CBEI_HRTEM_HCTEM_EFTEMRS(const eEM_Sim_Typ &sim_type)
		{
			return is_ISTEM(sim_type) || is_CBEI(sim_type) || is_HRTEM(sim_type) || is_HCTEM(sim_type) || is_EFTEMRS(sim_type);
		}

		inline
		dt_bool is_CBED_ED_EWFS_PED(const eEM_Sim_Typ &sim_type)
		{
			return is_CBED(sim_type) || is_ED(sim_type) || is_EWFS(sim_type) || is_PED(sim_type);
		}

		inline
		dt_bool is_obj_lens_temp_spat(const eEM_Sim_Typ &sim_type)
		{
			return is_ISTEM(sim_type) || is_CBEI(sim_type) || is_HRTEM(sim_type) || is_HCTEM(sim_type) || is_EFTEMRS(sim_type);
		}

		inline
		dt_bool is_cond_lens_temp_spat(const eEM_Sim_Typ &sim_type)
		{
			return is_STEM_ISTEM(sim_type) || is_CBED_CBEI(sim_type) || is_STEM_ISTEM_EELS(sim_type) || is_EFTEM(sim_type);
		}

		inline
		eSpace get_simulation_space(const eEM_Sim_Typ &sim_type)
		{
			return (is_simulation_type_FS(sim_type) || is_ISTEM_CBEI_HRTEM_HCTEM_EFTEMRS(sim_type))?mt::esp_fourier:mt::esp_real;
		}

		inline
		dt_bool is_scanning(const eEM_Sim_Typ &sim_type)
		{
			return is_STEM_ISTEM(sim_type) || is_STEM_ISTEM_EELS(sim_type);
		}

		inline
		eIncident_Wave_Typ validate_incident_wave_type(const eEM_Sim_Typ &sim_type, eIncident_Wave_Typ iw_type)
		{
			if (iw_type == eiwt_auto)
			{
				auto bb = is_scanning(sim_type) || is_CBED_CBEI(sim_type);
				bb = bb || ((is_EWFS_EWRS(sim_type) || is_EFTEM(sim_type) || is_IWFS_IWRS(sim_type)) && is_convergent_wave(iw_type));
				iw_type = (bb)?mt::eiwt_convergent_wave:mt::eiwt_plane_wave;
			}

			return iw_type;
		}

	}
#endif