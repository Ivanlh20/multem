/*
 * This file is part of MULTEM.
 * Copyright 2015 Ivan Lobato <Ivanlh20@gmail.com>
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
#include "types.hpp"
#include "atom_data.hpp"

namespace multem
{
	bool is_gpu_available();

	template<class T>
	DEVICE_CALLABLE FORCEINLINE 
	T get_Vr_factor(const T &E_0, const T &theta);

	template<class InputIterator, class TVector>
	void match_vectors(InputIterator A_base_first, InputIterator A_base_last, TVector &A_search);

	template<class T, eDevice dev>
	class Input_Multislice{
		public:
			using value_type = T;

			ePrecision precision;
			eDevice device; 									// eP_float = 1, eP_double = 2
			int cpu_ncores; 									// Number of Cores CPU
			int cpu_nthread; 									// Number of threads
			int gpu_device; 									// GPU device
			int gpu_nstream; 									// Number of streams

			eSimulation_Type simulation_type; 					// 11: Scanning, 12: ISTEM, 21: cbed, 22: cbei, 31: ED, 32: hrtem, 41: ped, 42: hci, ... 51: EW Fourier, 52: EW real
			ePhonon_Model phonon_model; 						// 1: Still atom model, 2: Absorptive potential model, 3: Frozen phonon
			eElec_Spec_Int_Model interaction_model; 			// eESIM_Multislice = 1, eESIM_Phase_Object= 2, eESIM_Weak_Phase_Object = 3
			ePotential_Slicing potential_slicing; 				// ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
			ePotential_Type potential_type;						// potential type: 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)

			int fp_nconf; 										// Number of frozen phonon configurations
			FP_Dim fp_dim; 										// Phonon dimensions
			int fp_dist; 										// 1: Gaussian (Phonon distribution)
			int fp_seed; 										// Random seed(frozen phonon)
			int fp_iconf;

			eMicroscope_Effect microscope_effect; 				// 1: Partial coherente mode, 2: transmission cross coefficient
			eSpatial_Temporal_Effect spatial_temporal_effect; 	// 1: Spatial and temporal, 2: Temporal, 3: Spatial
			
			eZero_Defocus_Type zero_defocus_type; 				// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
			T zero_defocus_plane; 								// Zero defocus plane
			
			eThickness_Type thickness_type; 					// eTT_Whole_Specimen = 1, eTT_Through_Slices = 2, eTT_Through_Planes = 3
			Vector<T, e_Host> thickness; 						// Array of thicknesses

			eInput_Wave_Type input_wave_type; 					// 1: Automatic, 2: User define
			Vector<complex<T>, dev> psi_0; 						// Input wave

			eOperation_Mode operation_mode;						// eOM_Normal = 1, eOM_Advanced = 2
			bool coherent_contribution;							// true , false
			bool slice_storage;									// true , false

			T E_0; 												// Acceleration volatage in KeV
			T theta; 											// incident tilt (in spherical coordinates) (rad)
			T phi; 												// incident tilt (in spherical coordinates) (rad)

			Grid<T> grid; 										// gridBT information

			T Vrl; 												// Atomic potential cut-off
			int nR; 											// Number of gridBT points

			Lens<T> lens; 										// Aberrations

			bool is_crystal;
			Atom_Data<T> atoms; 								// atoms

			EELS<T> eels_fr; 									// EELS

			Scanning<T> scanning; 								// Scanning

			Det_Cir<T, e_Host> det_cir; 						// Circular detectors

			HRTEM<T> hrtem;

			CBE_FR<T> cbe_fr;

			PE_FR<T> pe_fr;

			EW_FR<T> ew_fr;

			eBeam_Type beam_type;
			T conv_beam_wave_x;
			T conv_beam_wave_y;

			int islice;
			int nstream;
			bool dp_Shift; 										// Shift diffraction pattern

			Input_Multislice(): precision(eP_double), device(e_Host), cpu_ncores(1), cpu_nthread(4), gpu_device(0), gpu_nstream(8), 
						simulation_type(eST_EWRS), phonon_model(ePM_Still_Atom), interaction_model(eESIM_Multislice), 
						potential_slicing(ePS_Planes), potential_type(ePT_Lobato_0_12), fp_nconf(0), fp_dist(1), 
						fp_seed(1983), microscope_effect(eME_Partial_Coherent), spatial_temporal_effect(eSTE_Spatial_Temporal), 
						zero_defocus_type(eZDT_Last), zero_defocus_plane(0), operation_mode(eOM_Normal), coherent_contribution(false), 
						slice_storage(true), thickness_type(eTT_Whole_Specimen), dp_Shift(false), E_0(300), theta(0), phi(0), Vrl(c_Vrl), 
						nR(c_nR), is_crystal(false), input_wave_type(eIWT_Automatic), beam_type(eBT_Plane_Wave), conv_beam_wave_x(0), conv_beam_wave_y(0), 
						fp_iconf(0), islice(0), nstream(cpu_nthread){ };

			void validate_parameters()
			{
				if(is_whole_specimen() || thickness.empty())
				{
					thickness_type=eTT_Whole_Specimen;
					thickness.resize(1);
					thickness[0] = atoms.z_max;
				}
				else if(is_through_planes())
				{
					std::sort(thickness.begin(), thickness.end());
					atoms.Sort_by_z();
					atoms.get_z_layer();			
					multem::match_vectors(atoms.z_layer.begin(), atoms.z_layer.end(), thickness);
				}
				else if(is_through_slices())
				{
					std::sort(thickness.begin(), thickness.end());
					atoms.Sort_by_z();
					atoms.get_z_layer();

					Vector<T, e_Host> z_slice;
					atoms.get_z_slice(potential_slicing, grid.dz, atoms, z_slice);
					multem::match_vectors(z_slice.begin()+1, z_slice.end(), thickness);
				}

				nstream = (is_Host())?cpu_nthread:gpu_nstream;

				if(!is_float() && !is_double())
					precision = eP_float;

				if(!is_Host() && !is_Device())
					device = e_Host;

				if(!is_gpu_available())
				{
					device = e_Host;
				}

				fp_nconf = (is_frozen_phonon())?max(1, fp_nconf):1;

				fp_iconf = max(1, fp_iconf);

				fp_seed = max(0, fp_seed);

				islice = max(0, islice);

				gpu_device = max(0, gpu_device);

				if(isZero(Vrl))
				{
					Vrl = c_Vrl;
				}

				if(isZero(nR))
				{
					nR = c_nR;
				}

				dp_Shift = (is_PED())?true:false;

				if(is_scanning())
				{
					scanning.set_grid();
				}

				lens.set_input_data(E_0, grid);

				det_cir.set_input_data(E_0);

				theta = set_incident_angle(theta);
				pe_fr.theta = set_incident_angle(pe_fr.theta);

				//Set beam type
				if(is_user_define_wave())
				{
					beam_type = eBT_User_Define;
				}
				else if(is_convergent_beam_wave())
				{
					beam_type = eBT_Convergent;

					if(is_CBED_CBEI())
					{
						set_beam_position(cbe_fr.x0, cbe_fr.y0);
					}
					else if(is_EWFS_EWRS())
					{
						set_beam_position(ew_fr.x0, ew_fr.y0);
					}
				}
				else if(is_plane_wave())
				{
					beam_type = eBT_Plane_Wave;
				}

				if(is_EELS() || is_EFTEM())
				{
					coherent_contribution = false;
					interaction_model = multem::eESIM_Multislice;
				}

				if(is_EWFS_EWRS())
				{
					coherent_contribution = true;
				}

				slice_storage = true;
				if(is_CBED_CBEI() || is_ED_HRTEM()|| is_EWFS_EWRS() ||(is_STEM() && coherent_contribution))
				{
					slice_storage = false;
				}

			}
			/**************************************************************************************/
			bool is_whole_specimen() const
			{
				return thickness_type==eTT_Whole_Specimen;
			}

			bool is_through_slices() const
			{
				return thickness_type==eTT_Through_Slices;
			}

			bool is_through_planes() const
			{
				return thickness_type==eTT_Through_Planes;
			}
			/**************************************************************************************/
			bool is_user_define_wave() const
			{
				return input_wave_type == eIWT_User_Define;
			}

			bool is_convergent_beam_wave() const
			{
				return is_scanning() || is_CBED_CBEI() || (is_EWFS_EWRS() && ew_fr.convergent_beam);
			}

			bool is_plane_wave() const
			{
				return !(is_user_define_wave() || is_convergent_beam_wave());
			}

			/**************************************************************************************/
			void set_beam_position(const T &x, const T &y)
			{
				conv_beam_wave_x = x;
				conv_beam_wave_y = y;
			}

			void set_stem_beam_position(const int &idx)
			{
				set_beam_position(scanning.x[idx], scanning.y[idx]);
			}

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
				return get_Rx_pos_shift(conv_beam_wave_x);
			}

			T get_Ry_pos_shift()
			{
				return get_Ry_pos_shift(conv_beam_wave_y);
			}

			T set_incident_angle(const T &theta) const
			{
				T n = round(sin(theta)/(lens.lambda*grid.dg_min()));
				return (isZero(theta))?0:asin(n*lens.lambda*grid.dg_min());
			}

			T ifp_nconf() const
			{
				return 1.0/fp_nconf;
			}

			T Vr_factor() const
			{
				return get_Vr_factor(E_0, theta);
			}

			T gx_0() const
			{
				return sin(theta)*cos(phi)/lens.lambda;
			}

			T gy_0() const
			{
				return sin(theta)*sin(phi)/lens.lambda;
			}

			bool is_frozen_phonon() const
			{
				return phonon_model == ePM_Frozen_Phonon;
			}

			bool is_multislice() const
			{
				return interaction_model == multem::eESIM_Multislice;
			}

			bool is_subslicing() const
			{
				return is_multislice() && (potential_slicing == multem::ePS_dz_Sub);
			}

			bool is_STEM() const
			{
				return simulation_type == multem::eST_STEM;
			}

			bool is_ISTEM() const
			{
				return simulation_type == multem::eST_ISTEM;
			}

			bool is_CBED() const
			{
				return simulation_type == multem::eST_CBED;
			}

			bool is_CBEI() const
			{
				return simulation_type == multem::eST_CBEI;
			}

			bool is_ED() const
			{
				return simulation_type == multem::eST_ED;
			}

			bool is_HRTEM() const
			{
				return simulation_type == multem::eST_HRTEM;
			}

			bool is_PED() const
			{
				return simulation_type == multem::eST_PED;
			}

			bool is_HCI() const
			{
				return simulation_type == multem::eST_HCI;
			}

			bool is_EWFS() const
			{
				return simulation_type == multem::eST_EWFS;
			}

			bool is_EWRS() const
			{
				return simulation_type == multem::eST_EWRS;
			}

			bool is_EELS() const
			{
				return simulation_type == multem::eST_EELS;
			}

			bool is_EFTEM() const
			{
				return simulation_type == multem::eST_EFTEM;
			}

			bool is_ProbeFS() const
			{
				return simulation_type == multem::eST_ProbeFS;
			}

			bool is_ProbeRS() const
			{
				return simulation_type == multem::eST_ProbeRS;
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

			bool is_EELS_EFTEM() const
			{
				return is_EELS() || is_EFTEM();
			}

			bool is_ProbeFS_ProbeRS() const
			{
				return is_ProbeFS() || is_ProbeRS();
			}

			bool is_simulation_type_FS() const
			{
				return is_STEM() || is_CBED() || is_ED() || is_PED() || is_EWFS() || is_EELS(); 
			}

			bool is_simulation_type_RS() const
			{
				return is_ISTEM() || is_CBEI() || is_HRTEM() || is_HCI() || is_EWRS() || is_EFTEM(); 
			}

			bool is_HRTEM_HCI_EFTEM() const
			{
				return is_HRTEM() || is_HCI() || is_EFTEM(); 
			}

			bool is_scanning() const
			{
				return is_STEM() || is_ISTEM() || is_EELS();
			}

			bool is_Host() const
			{
				return device == multem::e_Host;
			}

			bool is_Device() const
			{
				return device == multem::e_Device;
			}

			bool is_float() const
			{
				return precision == multem::eP_float;
			}

			bool is_double() const
			{
				return precision == multem::eP_double;
			}

			bool is_float_Host() const
			{
				return is_float() && is_Host();
			}

			bool is_double_Host() const
			{
				return is_double() && is_Host();
			}

			bool is_float_Device() const
			{
				return is_float() && is_Device();
			}

			bool is_double_Device() const
			{
				return is_double() && is_Device();
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

} //namespace multem

#endif