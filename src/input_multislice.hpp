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

#include <vector>

#include "math.cuh"
#include "types.hpp"
#include "host_device_functions.cuh"
#include "atom_data.hpp"

namespace multem
{
	template<class T, eDevice dev>
	class Input_Multislice{
		public:
			using value_type = typename T;

			ePrecision precision;
			eDevice device; 										// eP_float = 1, eP_double = 2
			int cpu_ncores; 										// Number of Cores CPU
			int cpu_nthread; 										// Number of threads
			int gpu_device; 										// GPU device
			int gpu_nstream; 										// Number of streams

			eSimulation_Type simulation_type; 						// 11: Scanning, 12: ISTEM, 21: cbed, 22: cbei, 31: ED, 32: hrtem, 41: ped, 42: hci, ... 51: EW Fourier, 52: EW real
			ePhonon_Model phonon_model; 							// 1: Still atom model, 2: Absorptive potential model, 3: Frozen phonon
			eElec_Spec_Int_Model interaction_model; 				// eESIM_Mulstilice = 1, eESIM_Phase_Object= 2, eESIM_Weak_Phase_Object = 3
			ePotential_Slicing potential_slicing; 					// ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
			ePotential_Type potential_type;							// potential type: 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)

			int fp_nconf; 											// Number of frozen phonon configurations
			FP_Dim fp_dim; 											// Phonon dimensions
			int fp_dist; 											// 1: Gaussian (Phonon distribution)
			int fp_seed; 											// Random seed(frozen phonon)
			int fp_iconf;

			eMicroscope_Effect microscope_effect; 					// 1: Partial coherente mode, 2: transmission_fun cross coefficient
			eSpatial_Temporal_Effect spatial_temporal_effect; 		// 1: Spatial and temporal, 2: Temporal, 3: Spatial
			
			eZero_Defocus_Type zero_defocus_type; 					// 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
			T zero_defocus_plane; 									// Zero defocus plane
			
			eThickness_Type thickness_type; 						// eTT_Whole_Specimen = 1, eTT_Through_Thickness = 2, eTT_Through_Planes = 3
			Vector<T, Host> thickness; 							// Array of thicknesses

			eInput_Wave_Type input_wave_type; 						// 1: Automatic, 2: User define
			Vector<complex<T>, dev> psi_0; 						// Input wave

			bool fast_cal; 											// true: fast calculation(high memory consumption), false :normal mode(low memory consumption)
			
			T E_0; 												// Acceleration volatage in KeV
			T theta; 												// incident tilt (in spherical coordinates) (rad)
			T phi; 												// incident tilt (in spherical coordinates) (rad)

			Grid<T> grid; 											// gridBT information

			T Vrl; 												// Atomic potential cut-off
			int nR; 												// Number of gridBT points

			Lens<T> lens; 											// Aberrations

			Atom_Data<T> atoms; 									// atoms

			EELS<T> eels_fr; 										// EELS

			Scanning<T> scanning; 										// Scanning

			Det_Cir<T, Host> det_cir; 								// Circular detectors

			HRTEM<T> hrtem;

			CBE_FR<T> cbe_fr;

			PE_FR<T> pe_fr;

			EW_FR<T> ew_fr;

			eBeam_Type beam_type;
			T conv_beam_wave_x;
			T conv_beam_wave_y;

			int islice;
			int nstream;
			bool dp_Shift; 											// Shift diffraction pattern

			Input_Multislice(): precision(eP_double), device(Host), cpu_ncores(1), cpu_nthread(4), gpu_device(0), gpu_nstream(8), 
						simulation_type(eST_EWRS), phonon_model(ePM_Still_Atom), interaction_model(eESIM_Mulstilice), 
						potential_slicing(ePS_Planes), potential_type(ePT_Lobato_0_12), fp_nconf(0), fp_dist(1), 
						fp_seed(1983), microscope_effect(eME_Partial_Coherent), spatial_temporal_effect(eSTE_Spatial_Temporal), 
						zero_defocus_type(eZDT_Last), zero_defocus_plane(0), fast_cal(true), thickness_type(eTT_Whole_Specimen), 
						dp_Shift(false), E_0(300), theta(0), phi(0), Vrl(c_Vrl), nR(c_nR), input_wave_type(eIWT_Automatic), 
						beam_type(eBT_Plane_Wave), conv_beam_wave_x(0), conv_beam_wave_y(0), fp_iconf(0), islice(0), nstream(cpu_nthread) {};

			void validate_parameters()
			{
				nstream = (device==Host)?cpu_nthread:gpu_nstream;

				if((precision!=eP_float)&&(precision!=eP_double))
					precision=eP_float;

				if((device!=Host)&&(device!=Device))
					device=Host;

				fp_nconf = (is_frozen_phonon())?max(1, fp_nconf):1;

				fp_iconf = max(1, fp_iconf);

				fp_seed = max(0, fp_seed);

				islice = max(0, islice);

				gpu_device = max(0, gpu_device);

				switch(zero_defocus_type)
				{
					case eZDT_First:
					{
						zero_defocus_plane = atoms.z_min;
					}
					break;
					case eZDT_Middle:
					{
						zero_defocus_plane = 0.5*(atoms.z_min + atoms.z_max);
					}
					break;
					case eZDT_Last:
					{
						zero_defocus_plane = atoms.z_max;
					}
					break;
				}

				thickness_type = multem::eTT_Whole_Specimen;

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

				theta = set_incident_angle(theta);
				pe_fr.theta = set_incident_angle(pe_fr.theta);

				set_beam_type();
			}

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

			void set_beam_type()
			{
				if (input_wave_type==eIWT_User_Define)
				{
					beam_type = eBT_User_Define;
				}
				else if(is_scanning() || is_CBED_CBEI() || (is_EWFS_EWRS() && ew_fr.convergent_beam))
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
				else
				{
					beam_type = eBT_Plane_Wave;
				}
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
				return multem::get_Vr_factor(E_0, theta);
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
				return phonon_model==ePM_Frozen_Phonon;
			}

			bool is_multislice() const
			{
				return interaction_model == multem::eESIM_Mulstilice;
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

			bool is_scanning() const
			{
				return is_STEM() || is_ISTEM() || is_EELS();
			}

			bool is_ProbeFS_ProbeRS() const
			{
				return is_ProbeFS() || is_ProbeRS();
			}

			bool is_Host() const
			{
				return device==multem::Host;
			}

			bool is_Device() const
			{
				return device==multem::Device;
			}

			bool is_float() const
			{
				return precision==multem::eP_float;
			}

			bool is_double() const
			{
				return precision==multem::eP_double;
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
	};

} //namespace multem

#endif