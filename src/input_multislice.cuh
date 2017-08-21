/*
 * This file is part of MULTEM.
 * Copyright 2017 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef INPUT_MULTISLICE_H
#define INPUT_MULTISLICE_H

#include <algorithm>

#include "math.cuh"
#include "types.cuh"
#include "lin_alg_def.cuh"
#include "atom_data.hpp"
#include "slicing.hpp"
#include "memory_info.cuh"

namespace mt
{
	bool is_gpu_available();

	int number_of_gpu_available();

	template <class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T get_Vr_factor(const T &E_0, const T &theta);


	/**************************************************************************************/
	inline
	bool is_multislice(const eElec_Spec_Int_Model &int_model)
	{
		return int_model == mt::eESIM_Multislice;
	}

	inline
	bool is_phase_object(const eElec_Spec_Int_Model &int_model)
	{
		return int_model == mt::eESIM_Phase_Object;
	}

	inline
	bool is_weak_phase_object(const eElec_Spec_Int_Model &int_model)
	{
		return int_model == mt::eESIM_Weak_Phase_Object;
	}

	/**************************************************************************************/
	inline
	bool is_still_atom(const ePhonon_Model &pn_model)
	{
		return pn_model == ePM_Still_Atom;
	}

	inline
	bool is_absorptive_model(const ePhonon_Model &pn_model)
	{
		return pn_model == ePM_Absorptive_Model;
	}

	inline
	bool is_frozen_phonon(const ePhonon_Model &pn_model)
	{
		return pn_model == ePM_Frozen_Phonon;
	}

	inline
	bool is_frozen_phonon_single_conf(const ePhonon_Model &pn_model, const bool &pn_single_conf)
	{
		return is_frozen_phonon(pn_model) && pn_single_conf;
	}

	/**************************************************************************************/
	inline
	bool is_whole_spec(const eThick_Type &thick_type) 
	{
		return thick_type == eTT_Whole_Spec;
	}

	inline
	bool is_through_slices(const eThick_Type &thick_type) 
	{
		return thick_type == eTT_Through_Slices;
	}

	inline
	bool is_through_thick(const eThick_Type &thick_type) 
	{
		return thick_type == eTT_Through_Thick;
	}

	/**************************************************************************************/
	inline
	bool is_slicing_by_planes(const eElec_Spec_Int_Model &int_model, const ePotential_Slicing &pot_slic)
	{
		return mt::is_multislice(int_model) && (pot_slic == mt::ePS_Planes);
	}

	inline
	bool is_slicing_by_dz(const eElec_Spec_Int_Model &int_model, const ePotential_Slicing &pot_slic)
	{
		return mt::is_multislice(int_model) && (pot_slic == mt::ePS_dz_Proj);
	}

	inline
	bool is_subslicing(const eElec_Spec_Int_Model &int_model, const ePotential_Slicing &pot_slic)
	{
		return mt::is_multislice(int_model) && (pot_slic == mt::ePS_dz_Sub);
	}

	inline
	bool is_subslicing_whole_spec(const eElec_Spec_Int_Model &int_model, const ePotential_Slicing &pot_slic, const eThick_Type &thick_type)
	{
		return mt::is_subslicing(int_model, pot_slic) && is_whole_spec(thick_type);
	}

	/**************************************************************************************/
	inline
	bool is_plane_wave(eIncident_Wave_Type iw_type)
	{
		return iw_type == eIWT_Plane_Wave;
	}

	inline
	bool is_convergent_wave(eIncident_Wave_Type iw_type)
	{
		return iw_type == eIWT_Convergent_Wave;
	}

	inline
	bool is_user_define_wave(eIncident_Wave_Type iw_type)
	{
		return iw_type == eIWT_User_Define_Wave;
	}

	/**************************************************************************************/
	inline
	bool is_STEM(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_STEM;
	}

	inline
	bool is_ISTEM(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_ISTEM;
	}

	inline
	bool is_CBED(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_CBED;
	}

	inline
	bool is_CBEI(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_CBEI;
	}

	inline
	bool is_ED(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_ED;
	}

	inline
	bool is_HRTEM(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_HRTEM;
	}

	inline
	bool is_PED(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_PED;
	}

	inline
	bool is_HCTEM(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_HCTEM;
	}

	inline
	bool is_EWFS(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_EWFS;
	}

	inline
	bool is_EWRS(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_EWRS;
	}

	inline
	bool is_EWFS_SC(const eTEM_Sim_Type &sim_type, const ePhonon_Model &pn_model, const bool &pn_single_conf)
	{
		return is_EWFS(sim_type) && (!is_frozen_phonon(pn_model) || is_frozen_phonon_single_conf(pn_model, pn_single_conf));
	}

	inline
	bool is_EWRS_SC(const eTEM_Sim_Type &sim_type, const ePhonon_Model &pn_model, const bool &pn_single_conf)
	{
		return is_EWRS(sim_type) && (!is_frozen_phonon(pn_model) || is_frozen_phonon_single_conf(pn_model, pn_single_conf));
	}

	inline
	bool is_EELS(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_EELS;
	}

	inline
	bool is_EFTEM(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_EFTEM;
	}

	inline
	bool is_IWFS(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_IWFS;
	}

	inline
	bool is_IWRS(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_IWRS;
	}

	inline
	bool is_PPFS(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_PPFS;
	}

	inline
	bool is_PPRS(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_PPRS;
	}

	inline
	bool is_TFFS(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_TFFS;
	}

	inline
	bool is_TFRS(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_TFRS;
	}

	inline
	bool is_PropFS(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_PropFS;
	}

	inline
	bool is_PropRS(const eTEM_Sim_Type &sim_type)
	{
		return sim_type == mt::eTEMST_PropRS;
	}

	inline
	bool is_STEM_ISTEM(const eTEM_Sim_Type &sim_type)
	{
		return is_STEM(sim_type) || is_ISTEM(sim_type);
	}

	inline
	bool is_CBED_CBEI(const eTEM_Sim_Type &sim_type)
	{
		return is_CBED(sim_type) || is_CBEI(sim_type);
	}

	inline
	bool is_ED_HRTEM(const eTEM_Sim_Type &sim_type)
	{
		return is_ED(sim_type) || is_HRTEM(sim_type);
	}

	inline
	bool is_PED_HCTEM(const eTEM_Sim_Type &sim_type)
	{
		return is_PED(sim_type) || is_HCTEM(sim_type);
	}

	inline
	bool is_EWFS_EWRS(const eTEM_Sim_Type &sim_type)
	{
		return is_EWFS(sim_type) || is_EWRS(sim_type);
	}

	inline
	bool is_EWFS_EWRS_SC(const eTEM_Sim_Type &sim_type, const ePhonon_Model &pn_model, const bool &pn_single_conf)
	{
		return is_EWFS_SC(sim_type, pn_model, pn_single_conf) || is_EWRS_SC(sim_type, pn_model, pn_single_conf);
	}

	inline
	bool is_EWFS_convergent_wave(const eTEM_Sim_Type &sim_type, const eIncident_Wave_Type &iw_type)
	{
		return is_EWFS(sim_type) && is_convergent_wave(iw_type);
	}

	inline
	bool is_EWRS_convergent_wave(const eTEM_Sim_Type &sim_type, const eIncident_Wave_Type &iw_type)
	{
		return is_EWRS(sim_type) && is_convergent_wave(iw_type);
	}

	inline
	bool is_EW_convergent_wave(const eTEM_Sim_Type &sim_type, const eIncident_Wave_Type &iw_type)
	{
		return is_EWFS_EWRS(sim_type) && is_convergent_wave(iw_type);
	}

	inline
	bool is_EELS_EFTEM(const eTEM_Sim_Type &sim_type)
	{
		return is_EELS(sim_type) || is_EFTEM(sim_type);
	}

	inline
	bool is_IWFS_IWRS(const eTEM_Sim_Type &sim_type)
	{
		return is_IWFS(sim_type) || is_IWRS(sim_type);
	}

	inline
	bool is_PPFS_PPRS(const eTEM_Sim_Type &sim_type)
	{
		return is_PPFS(sim_type) || is_PPRS(sim_type);
	}

	inline
	bool is_TFFS_TFRS(const eTEM_Sim_Type &sim_type)
	{
		return is_TFFS(sim_type) || is_TFRS(sim_type);
	}

	inline
	bool is_PropFS_PropRS(const eTEM_Sim_Type &sim_type)
	{
		return is_PropFS(sim_type) || is_PropRS(sim_type);
	}

	inline
	bool is_grid_FS(const eTEM_Sim_Type &sim_type)
	{
		auto bb = is_CBED(sim_type) || is_ED(sim_type) || is_PED(sim_type) || is_EWFS(sim_type);
		bb = bb || is_IWFS(sim_type) || is_PPFS(sim_type) || is_TFFS(sim_type); 
		return bb;
	}

	inline
	bool is_grid_RS(const eTEM_Sim_Type &sim_type)
	{
		return !is_grid_FS(sim_type);
	}

	inline
	bool is_simulation_type_FS(const eTEM_Sim_Type &sim_type)
	{
		auto bb = is_STEM(sim_type) || is_CBED(sim_type) || is_ED(sim_type);
		bb = bb || is_PED(sim_type) || is_EWFS(sim_type) || is_EELS(sim_type);
		bb = bb || is_IWFS(sim_type) || is_PPFS(sim_type) || is_TFFS(sim_type);
		return bb;
	}

	inline
	bool is_simulation_type_RS(const eTEM_Sim_Type &sim_type)
	{
		return !is_simulation_type_FS(sim_type);
	}

	inline
	bool is_specimen_required(const eTEM_Sim_Type &sim_type)
	{
		return !(is_IWFS_IWRS(sim_type) || is_PropFS_PropRS(sim_type));
	}

	inline
	bool is_ISTEM_CBEI_HRTEM_HCTEM_EFTEM(const eTEM_Sim_Type &sim_type)
	{
		return is_ISTEM(sim_type) || is_CBEI(sim_type) || is_HRTEM(sim_type) || is_HCTEM(sim_type) || is_EFTEM(sim_type); 
	}

	inline
	bool is_CBED_ED_EWFS_PED(const eTEM_Sim_Type &sim_type)
	{
		return is_CBED(sim_type) || is_ED(sim_type) || is_EWFS(sim_type) || is_PED(sim_type); 
	}

	inline
	eSpace get_simulation_space(const eTEM_Sim_Type &sim_type)
	{
		return (is_simulation_type_FS(sim_type) || is_ISTEM_CBEI_HRTEM_HCTEM_EFTEM(sim_type))?mt::eS_Reciprocal:mt::eS_Real; 
	}

	inline
	bool is_scanning(const eTEM_Sim_Type &sim_type)
	{
		return is_STEM_ISTEM(sim_type) || is_EELS(sim_type);
	}

	/**************************************************************************************/
	inline
	eIncident_Wave_Type validate_incident_wave_type(const eTEM_Sim_Type &sim_type, eIncident_Wave_Type iw_type)
	{
		if(iw_type == eIWT_Auto)
		{
			auto bb = is_scanning(sim_type) || is_CBED_CBEI(sim_type);
			bb = bb || ((is_EWFS_EWRS(sim_type) || is_IWFS_IWRS(sim_type)) && is_convergent_wave(iw_type));
			iw_type = (bb)?mt::eIWT_Convergent_Wave:mt::eIWT_Plane_Wave;
		}

		return iw_type;
	}

	/**************************************************************************************/
	template <class T>
	struct Slicing;

	template <class T>
	class Input_Multislice
	{
		public:
			using value_type = T;

			eElec_Spec_Int_Model interaction_model; 			// eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
			ePotential_Type potential_type;						// potential type: 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)

			ePhonon_Model pn_model; 							// 1: Still atom model, 2: Absorptive potential model, 3: Frozen phonon
			bool pn_coh_contrib;								// true, false
			bool pn_single_conf; 								// single configuration: true, false			
			FP_Dim pn_dim; 										// Phonon dimensions
			int fp_dist; 										// 1: Gaussian (Phonon distribution)
			int pn_seed; 										// Random seed(frozen phonon)
			int pn_nconf; 										// true: single phonon configuration, false: number of frozen phonon configurations
			int fp_iconf_0;										// initial configuration

			Atom_Data<T> atoms; 								// atoms
			bool is_crystal;
			
			T spec_rot_theta; 									// angle
			r3d<T> spec_rot_u0; 								// unitary vector			
			eRot_Point_Type spec_rot_center_type; 				// 1: geometric center, 2: User define		
			r3d<T> spec_rot_center_p; 							// rotation point

			eThick_Type thick_type; 							// eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
			host_vector<T> thick; 								// Array of thickes

			ePotential_Slicing potential_slicing; 				// ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4

			Grid_2d<T> grid_2d; 								// grid information

			eTEM_Sim_Type simulation_type; 						// 11: Scanning, 12: ISTEM, 21: cbed, 22: cbei, 31: ED, 32: hrtem, 41: ped, 42: hci, ... 51: EW Fourier, 52: EW real

			eIncident_Wave_Type iw_type; 						// 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
			host_vector<complex<T>> iw_psi; 					// user define incident wave
			host_vector<T> iw_x;								// x position
			host_vector<T> iw_y; 								// y position

			T E_0; 												// Acceleration volatage in KeV
			T lambda;											// lambda
			T theta; 											// incident tilt (in spherical coordinates) (rad)
			T phi; 												// incident tilt (in spherical coordinates) (rad)

			eIllumination_Model illumination_model; 			// 1: Partial coherent mode, 2: transmission cross coefficient
			eTemporal_Spatial_Incoh temporal_spatial_incoh; 	// 1: Spatial and temporal, 2: Temporal, 3: Spatial

			Lens<T> cond_lens; 									// Condenser lens
			Lens<T> obj_lens; 									// Objective lens

			Scanning<T> scanning; 								// Scanning

			Detector<T, e_host> detector; 						// STEM Detectors

			EELS<T> eels_fr; 									// EELS

			eOperation_Mode operation_mode;						// eOM_Normal = 1, eOM_Advanced = 2
			bool slice_storage;									// true, false
			bool reverse_multislice;							// true, false
			int mul_sign;										// multislice sign
			
			T Vrl; 												// Atomic potential cut-off
			int nR; 											// Number of grid_bt points
			
			int nrot; 											// Total number of rotations
			
			eLens_Var_Type cdl_var_type; 						// eLVT_off = 0, eLVT_m = 1, eLVT_f = 2, eLVT_Cs3 = 3, eLVT_Cs5 = 4, eLVT_mfa2 = 5, eLVT_afa2 = 6, eLVT_mfa3 = 7, eLVT_afa3 = 8, eLVT_inner_aper_ang = 9, eLVT_outer_aper_ang = 10
			host_vector<T> cdl_var; 					// Array of thickes

			host_vector<int> iscan;
			host_vector<T> beam_x;						// temporal variables for
			host_vector<T> beam_y;

			int islice;
			bool dp_Shift; 										// Shift diffraction pattern

			Input_Multislice():simulation_type(eTEMST_EWRS), pn_model(ePM_Still_Atom), interaction_model(eESIM_Multislice), 
						potential_slicing(ePS_Planes), potential_type(ePT_Lobato_0_12), fp_dist(1), pn_seed(300183), 
						pn_single_conf(false), pn_nconf(1), fp_iconf_0(1), spec_rot_theta(0), spec_rot_u0(0, 0, 1), 
						spec_rot_center_type(eRPT_geometric_center), spec_rot_center_p(1, 0, 0), illumination_model(eIM_Partial_Coherent), 
						temporal_spatial_incoh(eTSI_Temporal_Spatial), thick_type(eTT_Whole_Spec), 
						operation_mode(eOM_Normal), pn_coh_contrib(false), slice_storage(false), reverse_multislice(false), 
						mul_sign(1), E_0(300), lambda(0), theta(0), phi(0), nrot(1), Vrl(c_Vrl), nR(c_nR), iw_type(eIWT_Plane_Wave), 
						is_crystal(false), islice(0), dp_Shift(false){};

			template <class TInput_Multislice>
			void assign(TInput_Multislice &input_multislice)
			{ 
				interaction_model = input_multislice.interaction_model;
				potential_type = input_multislice.potential_type;

				operation_mode = input_multislice.operation_mode;
				slice_storage = input_multislice.slice_storage;
				reverse_multislice = input_multislice.reverse_multislice;
				mul_sign = input_multislice.mul_sign;
				Vrl = input_multislice.Vrl;
				nR = input_multislice.nR;

				pn_model = input_multislice.pn_model;
				pn_coh_contrib = input_multislice.pn_coh_contrib;
				pn_dim = input_multislice.pn_dim;
				fp_dist = input_multislice.fp_dist;
				pn_seed = input_multislice.pn_seed;
				pn_single_conf = input_multislice.pn_single_conf;
				pn_nconf = input_multislice.pn_nconf;

				atoms = input_multislice.atoms;
				is_crystal = input_multislice.is_crystal;

				spec_rot_theta = input_multislice.spec_rot_theta;
				spec_rot_u0 = input_multislice.spec_rot_u0;
				spec_rot_center_type = input_multislice.spec_rot_center_type;
				spec_rot_center_p = input_multislice.spec_rot_center_p;

				thick_type = input_multislice.thick_type;
				thick = input_multislice.thick;

				potential_slicing = input_multislice.potential_slicing;

				grid_2d = input_multislice.grid_2d;

				simulation_type = input_multislice.simulation_type;

				iw_type = input_multislice.iw_type;
				iw_psi = input_multislice.iw_psi;
				iw_x = input_multislice.iw_x;
				iw_y = input_multislice.iw_y;

				E_0 = input_multislice.E_0;
				theta = input_multislice.theta;
				phi = input_multislice.phi;
				nrot = input_multislice.nrot;

				illumination_model = input_multislice.illumination_model;
				temporal_spatial_incoh = input_multislice.temporal_spatial_incoh;

				cond_lens = input_multislice.cond_lens;
				obj_lens = input_multislice.obj_lens;

				scanning = input_multislice.scanning;
				detector = input_multislice.detector;

				eels_fr = input_multislice.eels_fr;

				cdl_var_type = input_multislice.cdl_var_type;
				cdl_var = input_multislice.cdl_var;

				iscan = input_multislice.iscan;
				beam_x = input_multislice.beam_x;
				beam_y = input_multislice.beam_y;

				islice = input_multislice.islice;
				dp_Shift = input_multislice.dp_Shift;
			}
	
			template <class TInput_Multislice> 
			Input_Multislice<T>& operator=(TInput_Multislice &input_multislice)
			{
				assign(input_multislice);
				return *this; 
			}

			void validate_parameters()
			{
				pn_seed = max(1, pn_seed);

				pn_nconf = (!is_frozen_phonon())?1:max(1, pn_nconf);

				fp_iconf_0 = (!is_frozen_phonon())?1:(pn_single_conf)?pn_nconf:1;

				islice = max(0, islice);

				if(isZero(Vrl))
				{
					Vrl = c_Vrl;
				}

				if(isZero(nR))
				{
					nR = c_nR;
				}

				dp_Shift = (is_PED())?true:false;

				if(!is_scanning())
				{
					scanning.set_default();
				}
				scanning.set_grid();

				lambda = get_lambda(E_0);

				// multislice sign
				mul_sign = (reverse_multislice)?-1:1;

				theta = set_incident_angle(theta);
				nrot = max(1, nrot);
				if(!is_PED_HCTEM())
				{
					nrot = 1;
				}

				// set beam position
				set_beam_position(iw_x, iw_y);

				if(is_EELS_EFTEM())
				{
					pn_coh_contrib = false;
					interaction_model = mt::eESIM_Multislice;
					illumination_model = eIM_Coherent;
				}

				if(is_EWFS_EWRS())
				{
					pn_coh_contrib = true;
					illumination_model = eIM_Coherent;
				}

				// It will be changed when you include spatial incoherences (beta)
				if(is_CBED_CBEI() && !is_illu_mod_full_integration())
				{
					illumination_model = eIM_Coherent;
				}

				if (is_incoh_temporal_spatial() && isZero(cond_lens.beta))
				{
					temporal_spatial_incoh = eTSI_Temporal;
				}

				slice_storage = is_slice_storage();

				if(!is_multislice())
				{
					islice = 0;
					potential_slicing = ePS_Planes;
					if(is_through_slices())
					{
						thick_type = eTT_Through_Thick;
					}
					slice_storage = slice_storage || !is_whole_spec();
				}

				if(is_subslicing() && is_through_thick())
				{
					thick_type = eTT_Whole_Spec;
				}

				if(!is_subslicing())
				{
					pn_dim.z = false;
				}

				if(is_spec_rot_active())
				{
					thick_type = eTT_Whole_Spec;
                    spec_rot_u0.normalized();
					// get geometric center
					if(spec_rot_center_type == eRPT_geometric_center)
					{
						spec_rot_center_p = r3d<T>(atoms.x_mean, atoms.y_mean, atoms.z_mean);
					}
					// rotate atoms
					rotate_atoms(atoms, spec_rot_theta, spec_rot_u0, spec_rot_center_p);
					// get statistic
					atoms.get_statistic();
                    // reset theta
                    spec_rot_theta = 0;
				}

				// match slicing with the require thickness
				Slicing<T> slicing;
				slicing.match_thickness(potential_slicing, atoms, thick_type, thick);

				// remove atoms outside the thickness range
				T ee_z = 0.1;
				T z_e = thick.back()+ee_z;

				if(atoms.z_max>z_e)
				{
					T z_0 = atoms.z_min-ee_z;
					remove_atoms_outside_z_range(atoms, z_0, z_e);
					// get statistic
					atoms.get_statistic();
				}

				// set condenser lens parameters
				cond_lens.set_input_data(E_0, grid_2d);
				if(atoms.empty())
				{
					cond_lens.zero_defocus_plane = 0;
				}
				else
				{
					cond_lens.zero_defocus_plane = cond_lens.get_zero_defocus_plane(atoms.z_min, atoms.z_max);
				}

				cond_lens.zero_defocus_type = eZDT_User_Define;

				// set objetive lens parameters
 				obj_lens.set_input_data(E_0, grid_2d);
			}

			/**************************************************************************************/
			void set_reverse_multislice(bool rm)
			{
				reverse_multislice = rm;
				mul_sign = (reverse_multislice)?-1:1;
			}

			/**************************************************************************************/
			inline
			int number_conf()
			{
				return (pn_nconf-fp_iconf_0+1);
			}

			int number_of_beams()
			{
				return (is_plane_wave())?0:iw_x.size();
			}

			bool is_multi_beam()
			{
				return number_of_beams()>0;
			}

			template<class TVector>
			void set_beam_position(TVector &x, TVector &y)
			{
				beam_x.assign(x.begin(), x.end());
				beam_y.assign(y.begin(), y.end());
			}

			void set_iscan_beam_position()
			{
				beam_x.resize(iscan.size());
				beam_y.resize(iscan.size());
				for(auto is=0; is<iscan.size(); is++)
				{
					auto idx = iscan[is];
					beam_x[is] = scanning.x[idx];
					beam_y[is] = scanning.y[idx];
				}
			}

			/**************************************************************************************/
			bool is_spec_rot_active() const
			{
				return mt::nonZero(spec_rot_theta);
			}

			/**************************************************************************************/
			T Rx_exp_factor()
			{
				return grid_2d.exp_factor_Rx(beam_x);
			}

			T Ry_exp_factor()
			{
				return grid_2d.exp_factor_Ry(beam_y);
			}

			T set_incident_angle(const T &theta) const
			{
				T n = round(sin(theta)/(lambda*grid_2d.dg_min()));
				return (isZero(theta))?0:asin(n*lambda*grid_2d.dg_min());
			}

			T  get_weight() const
			{
				int nconf = (!is_frozen_phonon()||pn_single_conf)?1:pn_nconf;
				return 1.0/static_cast<T>(nconf*nrot);
			}

			void set_phi(const int &irot)
			{
				phi = (irot == 0)?0.0:(c_2Pi*static_cast<T>(irot))/static_cast<T>(nrot);
			}

			inline
			T get_propagator_factor(const T &z) const 
			{ 
				return (-mul_sign*c_Pi*lambda*z); 
			}

			T Vr_factor() const
			{
				return (mul_sign*get_Vr_factor(E_0, theta));
			}

			T gx_0() const
			{
				return sin(theta)*cos(phi)/lambda;
			}

			T gy_0() const
			{
				return sin(theta)*sin(phi)/lambda;
			}

			void set_eels_fr_atom(const int &iatoms, const Atom_Data<T> &atoms)
			{
				eels_fr.x = grid_2d.exp_factor_Rx(atoms.x[iatoms]);
				eels_fr.y = grid_2d.exp_factor_Ry(atoms.y[iatoms]);
				eels_fr.occ = atoms.occ[iatoms];
			}

			/**************************************************************************************/
			bool is_multislice() const
			{
				return mt::is_multislice(interaction_model);
			}

			bool is_phase_object() const
			{
				return mt::is_phase_object(interaction_model);
			}

			bool is_weak_phase_object() const
			{
				return mt::is_weak_phase_object(interaction_model);
			}

			/**************************************************************************************/

			bool is_still_atom() const
			{
				return mt::is_still_atom(pn_model);
			}

			bool is_absorptive_model() const
			{
				return mt::is_absorptive_model(pn_model);
			}

			bool is_frozen_phonon() const
			{
				return mt::is_frozen_phonon(pn_model);
			}

			bool is_frozen_phonon_single_conf() const
			{
				return mt::is_frozen_phonon_single_conf(pn_model, pn_single_conf);
			}

			/**************************************************************************************/
			bool is_whole_spec() const
			{
				return mt::is_whole_spec(thick_type);
			}

			bool is_through_slices() const
			{
				return mt::is_through_slices(thick_type);
			}

			bool is_through_thick() const
			{
				return mt::is_through_thick(thick_type);
			}

			/**************************************************************************************/
			bool is_slicing_by_planes() const
			{
				return mt::is_slicing_by_planes(interaction_model, potential_slicing);
			}

			bool is_slicing_by_dz() const
			{
				return mt::is_slicing_by_dz(interaction_model, potential_slicing);
			}

			bool is_subslicing() const
			{
				return mt::is_subslicing(interaction_model, potential_slicing);
			}

			bool is_subslicing_whole_spec() const
			{
				return mt::is_subslicing_whole_spec(interaction_model, potential_slicing, thick_type);
			}

			/**************************************************************************************/
			bool is_plane_wave() const
			{
				return mt::is_plane_wave(iw_type);
			}

			bool is_convergent_wave() const
			{
				return mt::is_convergent_wave(iw_type);
			}

			bool is_user_define_wave() const
			{
				return mt::is_user_define_wave(iw_type);
			}

			/**********************************************************/
			bool is_STEM() const
			{
				return mt::is_STEM(simulation_type);
			}

			bool is_ISTEM() const
			{
				return mt::is_ISTEM(simulation_type);
			}

			bool is_CBED() const
			{
				return mt::is_CBED(simulation_type);
			}

			bool is_CBEI() const
			{
				return mt::is_CBEI(simulation_type);
			}

			bool is_ED() const
			{
				return mt::is_ED(simulation_type);
			}

			bool is_HRTEM() const
			{
				return mt::is_HRTEM(simulation_type);
			}

			bool is_PED() const
			{
				return mt::is_PED(simulation_type);
			}

			bool is_HCTEM() const
			{
				return mt::is_HCTEM(simulation_type);
			}

			bool is_EWFS() const
			{
				return mt::is_EWFS(simulation_type);
			}

			bool is_EWRS() const
			{
				return mt::is_EWRS(simulation_type);
			}

			bool is_EWFS_SC() const
			{
				return mt::is_EWFS_SC(simulation_type, pn_model, pn_single_conf);
			}

			bool is_EWRS_SC() const
			{
				return mt::is_EWRS_SC(simulation_type, pn_model, pn_single_conf);
			}

			bool is_EELS() const
			{
				return mt::is_EELS(simulation_type);
			}

			bool is_EFTEM() const
			{
				return mt::is_EFTEM(simulation_type);
			}

			bool is_IWFS() const
			{
				return mt::is_IWFS(simulation_type);
			}

			bool is_IWRS() const
			{
				return mt::is_IWRS(simulation_type);
			}

			bool is_PPFS() const
			{
				return mt::is_PPFS(simulation_type);
			}

			bool is_PPRS() const
			{
				return mt::is_PPRS(simulation_type);
			}

			bool is_TFFS() const
			{
				return mt::is_TFFS(simulation_type);
			}

			bool is_TFRS() const
			{
				return mt::is_TFRS(simulation_type);
			}

			bool is_PropFS() const
			{
				return mt::is_PropFS(simulation_type);
			}

			bool is_PropRS() const
			{
				return mt::is_PropRS(simulation_type);
			}

			bool is_STEM_ISTEM() const
			{
				return mt::is_STEM_ISTEM(simulation_type);
			}

			bool is_CBED_CBEI() const
			{
				return mt::is_CBED_CBEI(simulation_type);
			}

			bool is_ED_HRTEM() const
			{
				return mt::is_ED_HRTEM(simulation_type);
			}

			bool is_PED_HCTEM() const
			{
				return mt::is_PED_HCTEM(simulation_type);
			}

			bool is_EWFS_EWRS() const
			{
				return mt::is_EWFS_EWRS(simulation_type);
			}

			bool is_EWFS_EWRS_SC() const
			{
				return mt::is_EWFS_EWRS_SC(simulation_type, pn_model, pn_single_conf);
			}

			bool is_EELS_EFTEM() const
			{
				return mt::is_EELS_EFTEM(simulation_type);
			}

			bool is_IWFS_IWRS() const
			{
				return mt::is_IWFS_IWRS(simulation_type);
			}

			bool is_PPFS_PPRS() const
			{
				return mt::is_PPFS_PPRS(simulation_type);
			}

			bool is_TFFS_TFRS() const
			{
				return mt::is_TFFS_TFRS(simulation_type);
			}

			bool is_PropFS_PropRS() const
			{
				return mt::is_PropFS_PropRS(simulation_type);
			}

			bool is_grid_FS() const
			{
				return mt::is_grid_FS(simulation_type);
			}

			bool is_grid_RS() const
			{
				return mt::is_grid_RS(simulation_type);
			}

			bool is_simulation_type_FS() const
			{
				return mt::is_simulation_type_FS(simulation_type);
			}

			bool is_simulation_type_RS() const
			{
				return mt::is_simulation_type_RS(simulation_type);
			}

			bool is_specimen_required() const
			{
				return mt::is_specimen_required(simulation_type);
			}
			
			bool is_ISTEM_CBEI_HRTEM_HCTEM_EFTEM() const
			{
				return mt::is_ISTEM_CBEI_HRTEM_HCTEM_EFTEM(simulation_type);
			}

			bool is_CBED_ED_EWFS_PED() const
			{
				return mt::is_CBED_ED_EWFS_PED(simulation_type);
			}

			eSpace get_simulation_space() const
			{
				return mt::get_simulation_space(simulation_type);
			}

			bool is_scanning() const
			{
				return mt::is_scanning(simulation_type);
			}

			/**************************************************************************************/
			void set_incident_wave_type(eIncident_Wave_Type iw_type_i)
			{
				iw_type = mt::validate_incident_wave_type(simulation_type, iw_type_i);
			}

			/**************************************************************************************/
			bool is_illu_mod_coherent() const
			{
				return illumination_model == eIM_Coherent;
			}

			bool is_illu_mod_partial_coherent() const
			{
				return illumination_model == eIM_Partial_Coherent;
			}

			bool is_illu_mod_trans_cross_coef() const
			{
				return illumination_model == eIM_Trans_Cross_Coef;
			}

			bool is_illu_mod_full_integration() const
			{
				return illumination_model == eIM_Full_Integration;
			}

			/**************************************************************************************/
			bool is_incoh_temporal_spatial() const
			{
				return temporal_spatial_incoh == eTSI_Temporal_Spatial;
			}

			bool is_incoh_temporal() const
			{
				return temporal_spatial_incoh == eTSI_Temporal;
			}

			bool is_incoh_spatial() const
			{
				return temporal_spatial_incoh == eTSI_Spatial;
			}

			/**************************************************************************************/
			bool is_detector_circular() const
			{
				return detector.is_detector_circular();
			}

			bool is_detector_radial() const
			{
				return detector.is_detector_radial();
			}

			bool is_detector_matrix() const
			{
				return detector.is_detector_matrix();
			}

			/**************************************************************************************/
			bool is_slice_storage() const
			{
				bool slice_storage = is_PED_HCTEM() || is_EELS_EFTEM();
				slice_storage = slice_storage || is_ISTEM() || (is_STEM() && !pn_coh_contrib);
				slice_storage = slice_storage || (is_CBED_CBEI() && is_illu_mod_full_integration());

				return slice_storage;
			}

			bool is_operation_mode_normal() const
			{
				return operation_mode == eOM_Normal;
			}

			bool is_operation_mode_advanced() const
			{
				return operation_mode == eOM_Advanced;
			}

			/**************************************************************************************/
			bool is_lvt_off() const
			{
				return cdl_var == eLVT_off;
			}

			bool is_lvt_m() const
			{
				return cdl_var == eLVT_m;
			}

			bool is_lvt_Cs3() const
			{
				return cdl_var == eLVT_Cs3;
			}

			bool is_lvt_Cs5() const
			{
				return cdl_var == eLVT_Cs5;
			}

			bool is_lvt_mfa2() const
			{
				return cdl_var == eLVT_mfa2;
			}

			bool is_lvt_afa2() const
			{
				return cdl_var == eLVT_afa2;
			}

			bool is_lvt_mfa3() const
			{
				return cdl_var == eLVT_mfa3;
			}

			bool is_lvt_afa3 () const
			{
				return cdl_var == eLVT_afa3;
			}

			bool is_lvt_inner_aper_ang() const
			{
				return cdl_var == eLVT_inner_aper_ang;
			}

			bool is_lvt_outer_aper_ang () const
			{
				return cdl_var == eLVT_outer_aper_ang;
			}

	};

} // namespace mt

#endif