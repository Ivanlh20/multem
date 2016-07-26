/*
 * This file is part of MULTEM.
 * Copyright 2016 Ivan Lobato <Ivanlh20@gmail.com>
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
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INPUT_MULTISLICE_H
#define INPUT_MULTISLICE_H

#include <algorithm>
#include <vector>

#include "math.cuh"
#include "types.cuh"
#include "lin_alg_def.cuh"
#include "atom_data.hpp"
#include "memory_info.cuh"

namespace mt
{
	bool is_gpu_available();

	int number_of_gpu_available();

	template<class T>
	DEVICE_CALLABLE FORCE_INLINE 
	T get_Vr_factor(const T &E_0, const T &theta);

	template<class InputIterator, class TVector>
	void match_vectors(InputIterator A_base_first, InputIterator A_base_last, TVector &A_search);

	template<class T>
	class Input_Multislice
	{
		public:
			using value_type = T;

			ePrecision precision;
			eDevice device; 									// eP_float = 1, eP_double = 2
			int cpu_ncores; 									// Number of Cores CPU
			int cpu_nthread; 									// Number of threads
			int gpu_device; 									// GPU device
			int gpu_nstream; 									// Number of streams

			eTEM_Sim_Type simulation_type; 						// 11: Scanning, 12: ISTEM, 21: cbed, 22: cbei, 31: ED, 32: hrtem, 41: ped, 42: hci, ... 51: EW Fourier, 52: EW real
			ePhonon_Model phonon_model; 						// 1: Still atom model, 2: Absorptive potential model, 3: Frozen phonon
			eElec_Spec_Int_Model interaction_model; 			// eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
			ePotential_Slicing potential_slicing; 				// ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
			ePotential_Type potential_type;						// potential type: 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)

			FP_Dim fp_dim; 										// Phonon dimensions
			int fp_dist; 										// 1: Gaussian (Phonon distribution)
			int fp_seed; 										// Random seed(frozen phonon)
			bool fp_single_conf; 								// single configuration: true , false
			int fp_nconf; 										// true: single phonon configuration, false: number of frozen phonon configurations
			int fp_iconf_0;										// initial configuration

			bool tm_active;										// true or false
			T tm_theta; 										// angle
			r3d<T> tm_u0; 										// unitary vector			
			eRot_Point_Type tm_rot_point_type; 					// 1: geometric center, 2: User define		
			r3d<T> tm_p0; 										// rotation point

			eIllumination_Model illumination_model; 				// 1: Partial coherente mode, 2: transmission cross coefficient
			eTemporal_Spatial_Incoh temporal_spatial_incoh; 	// 1: Spatial and temporal, 2: Temporal, 3: Spatial

			eThickness_Type thickness_type; 					// eTT_Whole_Specimen = 1, eTT_Through_Thickness = 2, eTT_Through_Slices = 3
			Vector<T, e_host> thickness; 						// Array of thicknesses

			eOperation_Mode operation_mode;						// eOM_Normal = 1, eOM_Advanced = 2
			bool coherent_contribution;							// true , false
			bool slice_storage;									// true , false

			T E_0; 												// Acceleration volatage in KeV
			T lambda;											// lambda
			T theta; 											// incident tilt (in spherical coordinates) (rad)
			T phi; 												// incident tilt (in spherical coordinates) (rad)
			int nrot; 											// Total number of rotations

			Grid<T> grid; 										// gridBT information

			T Vrl; 												// Atomic potential cut-off
			int nR; 											// Number of gridBT points

			eIncident_Wave_Type iw_type; 						// 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
			Vector<complex<T>, e_host> iw_psi; 					// user define incident wave
			T iw_x;												// x position
			T iw_y; 											// y position

			Lens<T> cond_lens; 									// Condenser lens
			Lens<T> obj_lens; 									// Objective lens

			eLens_Var_Type cdl_var_type; 						// eLVT_off = 0, eLVT_m = 1, eLVT_f = 2, eLVT_Cs3 = 3, eLVT_Cs5 = 4, eLVT_mfa2 = 5, eLVT_afa2 = 6, eLVT_mfa3 = 7, eLVT_afa3 = 8, eLVT_inner_aper_ang = 9, eLVT_outer_aper_ang = 10
			Vector<T, e_host> cdl_var; 							// Array of thicknesses

			bool is_crystal;
			Atom_Data<T> atoms; 								// atoms

			EELS<T> eels_fr; 									// EELS

			Scanning<T> scanning; 								// Scanning

			Detector<T, e_host> detector; 						// STEM Detectors

			T beam_x;
			T beam_y;
			int iscan;
			int islice;
			bool dp_Shift; 										// Shift diffraction pattern
			int nstream;

			Input_Multislice(): precision(eP_double), device(e_host), cpu_ncores(1), cpu_nthread(4), gpu_device(0), gpu_nstream(8), 
						simulation_type(eTEMST_EWRS), phonon_model(ePM_Still_Atom), interaction_model(eESIM_Multislice), 
						potential_slicing(ePS_Planes), potential_type(ePT_Lobato_0_12), fp_dist(1), fp_seed(300183), 
						fp_single_conf(false), fp_nconf(1), fp_iconf_0(1), tm_active(false), tm_theta(0), tm_u0(0, 0, 1), 
						tm_rot_point_type(eRPT_geometric_center), tm_p0(1, 0, 0), illumination_model(eIM_Partial_Coherent), 
						temporal_spatial_incoh(eTSI_Temporal_Spatial), thickness_type(eTT_Whole_Specimen), 
						operation_mode(eOM_Normal), coherent_contribution(false), slice_storage(false), E_0(300), lambda(0), 
						theta(0), phi(0), nrot(1), Vrl(c_Vrl), nR(c_nR), iw_type(eIWT_Plane_Wave), iw_x(0), iw_y(0), 
						is_crystal(false), beam_x(0), beam_y(0), iscan(0), islice(0), dp_Shift(false), nstream(cpu_nthread){};

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

				cpu_nthread = max(1, cpu_nthread);
				gpu_nstream = max(1, gpu_nstream);
				nstream = (is_host())?cpu_nthread:gpu_nstream;

				set_device();

				fp_seed = max(0, fp_seed);

				fp_nconf = (!is_frozen_phonon())?1:max(1, fp_nconf);

				fp_iconf_0 = (!is_frozen_phonon())?1:(fp_single_conf)?fp_nconf:1;

				tm_u0.normalized();
				if(tm_rot_point_type == eRPT_geometric_center)
				{
					tm_p0 = r3d<T>(atoms.x_mean, atoms.y_mean, atoms.z_mean);
				}

				if(is_tomography())
				{
					thickness_type = eTT_Whole_Specimen;
					if(potential_slicing == ePS_Planes)
					{
						potential_slicing = ePS_dz_Proj;
					}
				}

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

				cond_lens.set_input_data(E_0, grid);
				if(!cond_lens.is_zero_defocus_type_First())
				{
					cond_lens.zero_defocus_type = eZDT_User_Define;
				}
				obj_lens.set_input_data(E_0, grid);

				theta = set_incident_angle(theta);
				nrot = max(1, nrot);
				if(!is_PED_HCI())
				{
					nrot = 1;
				}

				// set beam position
				set_beam_position(iw_x, iw_y);

				if(is_EELS() || is_EFTEM())
				{
					coherent_contribution = false;
					interaction_model = mt::eESIM_Multislice;
					illumination_model = eIM_Coherent;
				}

				if(is_EWFS_EWRS())
				{
					coherent_contribution = true;
				}

				slice_storage = false;
				if(is_PED_HCI() || is_EELS_EFTEM()|| is_ISTEM() ||(is_STEM() && !coherent_contribution))
				{
					slice_storage = true;
				}

				if(!is_multislice())
				{
					islice = 0;
					potential_slicing = ePS_Planes;
					if(is_through_slices())
					{
						thickness_type = eTT_Through_Thickness;
					}
					slice_storage = slice_storage || !is_whole_specimen();
				}

				if(is_subslicing() && is_through_thickness())
				{
					thickness_type = eTT_Whole_Specimen;
				}

				if(!is_subslicing())
				{
					fp_dim.z = false;
				}

				if(is_whole_specimen() || thickness.empty())
				{
					thickness_type = eTT_Whole_Specimen;
					thickness.resize(1);
					thickness[0] = atoms.z_max;
				}
				else if(is_through_thickness())
				{
					// for amorphous specimen it has to be modified
					thickness_type = eTT_Through_Thickness;
					std::sort(thickness.begin(), thickness.end());
					fp_dim.z = false;
					atoms.Sort_by_z();
					atoms.get_z_layer();			
					mt::match_vectors(atoms.z_layer.begin(), atoms.z_layer.end(), thickness);
				}
				else if(is_through_slices())
				{
					std::sort(thickness.begin(), thickness.end());
					atoms.Sort_by_z();
					atoms.get_z_layer();

					Vector<T, e_host> z_slice;
					atoms.get_z_slice(potential_slicing, grid.dz, atoms, z_slice);
					mt::match_vectors(z_slice.begin()+1, z_slice.end(), thickness);
				}
			}

			/**************************************************************************************/
			bool is_whole_specimen() const
			{
				return thickness_type == eTT_Whole_Specimen;
			}

			bool is_through_slices() const
			{
				return thickness_type == eTT_Through_Slices;
			}

			bool is_through_thickness() const
			{
				return thickness_type == eTT_Through_Thickness;
			}
			/**************************************************************************************/
			void set_incident_wave_type(eIncident_Wave_Type iw_type_i)
			{
				if(iw_type_i == eIWT_Auto)
				{
					if(is_scanning() || is_CBED_CBEI() || ((is_EWFS_EWRS() || is_IWFS_IWRS()) && (iw_type == eIWT_Convergent_Wave)))
					{
						iw_type = eIWT_Convergent_Wave;
					}
					else
					{
						iw_type = eIWT_Plane_Wave;
					}
				}
				else
				{
					iw_type = iw_type_i;
				}
			}

			bool is_plane_wave() const
			{
				return iw_type == eIWT_Plane_Wave;
			}

			bool is_convergent_wave() const
			{
				return iw_type == eIWT_Convergent_Wave;
			}

			bool is_user_define_wave() const
			{
				return iw_type == eIWT_User_Define_Wave;
			}

			/**************************************************************************************/
			void set_beam_position(const T &x, const T &y)
			{
				beam_x = x;
				beam_y = y;
			}

			void set_beam_position(const int &iscan_i)
			{
				iscan = iscan_i;
				set_beam_position(scanning.x[iscan], scanning.y[iscan]);
			}

			/**************************************************************************************/
			bool is_tomography() const
			{
				return tm_active;
			}

			/**************************************************************************************/
			T get_Rx_pos_shift(const T &x)
			{
				return -c_2Pi*(x-0.5*grid.lx);
			}

			T get_Ry_pos_shift(const T &y)
			{
				return -c_2Pi*(y-0.5*grid.ly);
			}

			T get_Rx_pos_shift()
			{
				return get_Rx_pos_shift(beam_x);
			}

			T get_Ry_pos_shift()
			{
				return get_Ry_pos_shift(beam_y);
			}

			T set_incident_angle(const T &theta) const
			{
				T n = round(sin(theta)/(lambda*grid.dg_min()));
				return (isZero(theta))?0:asin(n*lambda*grid.dg_min());
			}

			T get_weight() const
			{
				int nconf = (!is_frozen_phonon()||fp_single_conf)?1:fp_nconf;
				return 1.0/static_cast<T>(nconf*nrot);
			}

			void set_phi(const int &irot)
			{
				phi = (irot == 0)?0.0:(c_2Pi*static_cast<T>(irot))/static_cast<T>(nrot);
			}

			inline
			T get_propagator_factor(const T &z) const 
			{ 
				return -c_Pi*lambda*z; 
			}

			T Vr_factor() const
			{
				return get_Vr_factor(E_0, theta);
			}

			T gx_0() const
			{
				return sin(theta)*cos(phi)/lambda;
			}

			T gy_0() const
			{
				return sin(theta)*sin(phi)/lambda;
			}

			void set_eels_fr_atom(const int &iatom, const Atom_Data<T> &atoms)
			{
				eels_fr.x = get_Rx_pos_shift(atoms.x[iatom]);
				eels_fr.y = get_Ry_pos_shift(atoms.y[iatom]);
				eels_fr.occ = atoms.occ[iatom];
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

			/**************************************************************************************/

			bool is_Still_Atom() const
			{
				return phonon_model == ePM_Still_Atom;
			}

			bool is_Absorptive() const
			{
				return phonon_model == ePM_Absorptive;
			}

			bool is_frozen_phonon() const
			{
				return phonon_model == ePM_Frozen_Phonon;
			}

			bool is_frozen_phonon_single_conf() const
			{
				return is_frozen_phonon() && fp_single_conf;
			}

			bool is_multislice() const
			{
				return interaction_model == mt::eESIM_Multislice;
			}

			bool is_phase_object() const
			{
				return interaction_model == mt::eESIM_Phase_Object;
			}

			bool is_weak_phase_object() const
			{
				return interaction_model == mt::eESIM_Weak_Phase_Object;
			}

			bool is_slicing_by_planes() const
			{
				return is_multislice() && (potential_slicing == mt::ePS_Planes);
			}

			bool is_slicing_by_dz() const
			{
				return is_multislice() && (potential_slicing == mt::ePS_dz_Proj);
			}

			bool is_subslicing() const
			{
				return is_multislice() && (potential_slicing == mt::ePS_dz_Sub);
			}

			bool is_subslicing_whole_specimen() const
			{
				return is_subslicing() && is_whole_specimen();
			}

			bool is_STEM() const
			{
				return simulation_type == mt::eTEMST_STEM;
			}

			bool is_ISTEM() const
			{
				return simulation_type == mt::eTEMST_ISTEM;
			}

			bool is_CBED() const
			{
				return simulation_type == mt::eTEMST_CBED;
			}

			bool is_CBEI() const
			{
				return simulation_type == mt::eTEMST_CBEI;
			}

			bool is_ED() const
			{
				return simulation_type == mt::eTEMST_ED;
			}

			bool is_HRTEM() const
			{
				return simulation_type == mt::eTEMST_HRTEM;
			}

			bool is_PED() const
			{
				return simulation_type == mt::eTEMST_PED;
			}

			bool is_HCI() const
			{
				return simulation_type == mt::eTEMST_HCI;
			}

			bool is_EWFS() const
			{
				return simulation_type == mt::eTEMST_EWFS;
			}

			bool is_EWRS() const
			{
				return simulation_type == mt::eTEMST_EWRS;
			}

			bool is_EWFS_SC() const
			{
				return is_EWFS() && (!is_frozen_phonon() || is_frozen_phonon_single_conf());
			}

			bool is_EWRS_SC() const
			{
				return is_EWRS() && (!is_frozen_phonon() || is_frozen_phonon_single_conf());
			}

			bool is_EELS() const
			{
				return simulation_type == mt::eTEMST_EELS;
			}

			bool is_EFTEM() const
			{
				return simulation_type == mt::eTEMST_EFTEM;
			}

			bool is_IWFS() const
			{
				return simulation_type == mt::eTEMST_IWFS;
			}

			bool is_IWRS() const
			{
				return simulation_type == mt::eTEMST_IWRS;
			}

			bool is_PPFS() const
			{
				return simulation_type == mt::eTEMST_PPFS;
			}

			bool is_PPRS() const
			{
				return simulation_type == mt::eTEMST_PPRS;
			}

			bool is_TFFS() const
			{
				return simulation_type == mt::eTEMST_TFFS;
			}

			bool is_TFRS() const
			{
				return simulation_type == mt::eTEMST_TFRS;
			}

			bool is_PropFS() const
			{
				return simulation_type == mt::eTEMST_PropFS;
			}

			bool is_PropRS() const
			{
				return simulation_type == mt::eTEMST_PropRS;
			}

			bool is_STEM_ISTEM() const
			{
				return is_STEM() || is_ISTEM();
			}

			bool is_CBED_CBEI() const
			{
				return is_CBED() || is_CBEI();
			}

			bool is_ED_HRTEM() const
			{
				return is_ED() || is_HRTEM();
			}

			bool is_PED_HCI() const
			{
				return is_PED() || is_HCI();
			}

			bool is_EWFS_EWRS() const
			{
				return is_EWFS() || is_EWRS();
			}

			bool is_EWFS_EWRS_SC() const
			{
				return is_EWFS_SC() || is_EWRS_SC();
			}

			bool is_EELS_EFTEM() const
			{
				return is_EELS() || is_EFTEM();
			}

			bool is_IWFS_IWRS() const
			{
				return is_IWFS() || is_IWRS();
			}

			bool is_PPFS_PPRS() const
			{
				return is_PPFS() || is_PPRS();
			}

			bool is_TFFS_TFRS() const
			{
				return is_TFFS() || is_TFRS();
			}

			bool is_PropFS_PropRS() const
			{
				return is_PropFS() || is_PropRS();
			}

			bool is_simulation_type_FS() const
			{
				return is_STEM() || is_CBED() || is_ED() || is_PED() || is_EWFS() ||
					is_EELS() || is_IWFS() || is_PPFS() || is_TFFS(); 
			}

			bool is_simulation_type_RS() const
			{
				return !is_simulation_type_FS();
			}

			bool is_ISTEM_CBEI_HRTEM_HCI_EFTEM() const
			{
				return is_ISTEM() || is_CBEI() || is_HRTEM() || is_HCI() || is_EFTEM(); 
			}

			eSpace get_simulation_space() const
			{
				return (is_simulation_type_FS() || is_ISTEM_CBEI_HRTEM_HCI_EFTEM())?eS_Reciprocal:eS_Real; 
			}
			bool is_scanning() const
			{
				return is_STEM() || is_ISTEM() || is_EELS();
			}

			bool is_detector_circular() const
			{
				return detector.type == eDT_Circular;
			}

			bool is_detector_radial() const
			{
				return detector.type == eDT_Radial;
			}

			bool is_detector_matrix() const
			{
				return detector.type == eDT_Matrix;
			}

			void set_device()
			{
				if(is_device())
				{
					if(!is_gpu_available())
					{
						device = mt::e_host;
					} 
					else
					{
						auto ngpu = number_of_gpu_available();
						gpu_device = min(max(0, gpu_device), ngpu-1);
						cudaSetDevice(gpu_device);
					}
				}
				else
				{
					device = mt::e_host;
				}
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

			bool is_operation_mode_normal() const
			{
				return operation_mode == eOM_Normal;
			}

			bool is_operation_mode_advanced() const
			{
				return operation_mode == eOM_Advanced;
			}
	};

} // namespace mt

#endif