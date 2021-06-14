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

#ifndef IN_CLASSES_MT_H
	#define IN_CLASSES_MT_H

	#include <algorithm>

	#include "math.cuh"
	#include "types.cuh"
	#include "particles.cuh"
	#include "slicing.hpp"
	#include "cgpu_info.cuh"
	#include "cgpu_vctr.cuh"

	/* multem */
	namespace mt
	{
		dt_bool is_gpu_avbl();

		dt_int32 gpu_n_avbl();

		template <class T>
		class EELS;

		/**************************** Input multem *****************************/
		template <class T>
		struct Slicing;

		template <class T>
		class In_Multem
		{
		public:
			using value_type = T;

			System_Config system_config;						// System information

			eElec_Spec_Int_Model interaction_model;			// eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
			ePot_Parm_Typ pot_parm_typ;						// potential type: 1: doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)

			Phonon_Par phonon_par;							// phonon parameters

			Ptc_Atom<T> atoms;								// atoms
			dt_bool is_xtl;

			Rot_Par<T> rot_par;								// specimen rotation parameters

			eThick_Typ thick_type;							// eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
			Vctr_cpu<T> thick;								// Array of thicknesses 

			ePot_Sli_Typ potential_slicing;					// ePST_Planes = 1, ePST_dz_Proj = 2, ePST_dz_Sub = 3, ePST_Auto = 4

			Grid_2d<T> grid_2d;								// grid information

			iThread_Rect_2d output_area;					// Output region information

			eTEM_Sim_Typ simulation_type;					// 11: Scan, 12: ISTEM, 21: cbed, 22: cbei, 31: ED, 32: hrtem, 41: ped, 42: hci, ... 51: EW Fourier, 52: EW real

			eIncident_Wave_Typ iw_type;						// 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
			Vctr_cpu<complex<T>> iw_psi;					// user define incident wave

			Beam_Pos_2d<T> beam_pos_i;						// in beam positions

			T E_0;											// Acceleration volatage in KeV
			T lambda;										// lambda
			T theta;										// incident tilt (in spherical coordinates) (rad)
			T phi;											// incident tilt (in spherical coordinates) (rad)

			eIllumination_Model illumination_model;			// 1: Partial coherente mode, 2: transmission cross coefficient
			eTemporal_Spatial_Incoh temporal_spatial_incoh;	// 1: Spatial and temporal, 2: Temporal, 3: Spatial

			Lens<T> cond_lens;								// Condenser lens
			Lens<T> obj_lens;								// Objective lens

			Scan<T> scanning;								// Scan

			Detector<T, edev_cpu> detector;					// STEM Detectors

			EELS<T> eels_fr;								// EELS

			eOperation_Mode operation_mode;					// eOM_Normal = 1, eOM_Advanced = 2
			dt_bool slice_storage;							// true, false
			dt_bool reverse_multislice;						// reverse multislice (true, false)
			dt_int32 mul_sign;								// tem_simulation sign

			T Vrl;											// atomic potential cut-off
			dt_int32 nR;									// number of d_grid_blk points

			dt_int32 nrot;									// Total number of rotations

			eLens_Var_Typ cdl_var_type;						// eLVT_off = 0, eLVT_m = 1, eLVT_f = 2, eLVT_Cs3 = 3, eLVT_Cs5 = 4, eLVT_mfa2 = 5, eLVT_afa2 = 6, eLVT_mfa3 = 7, eLVT_afa3 = 8, eLVT_inner_aper_ang = 9, eLVT_outer_aper_ang = 10
			Vctr_cpu<T> cdl_var;							// Array of thicknesses 

			Beam_Pos_2d<T> beam_pos;						// beam positions

			dt_int32 islice;
			dt_bool dp_Shift;								// Shift diffraction pattern

			In_Multem():simulation_type(eTEMST_EWRS), interaction_model(eESIM_Multislice), 
				potential_slicing(ePST_Planes), pot_parm_typ(ePPT_lobato_0_12), illumination_model(eIM_Partial_Coherent), 
				temporal_spatial_incoh(eTSI_Temporal_Spatial), thick_type(eTT_Whole_Spec), 
				operation_mode(eOM_Normal), slice_storage(false), reverse_multislice(false), 
				mul_sign(1), E_0(300), lambda(0), theta(0), phi(0), nrot(1), Vrl(c_vr_min), nR(c_nR), iw_type(eIWT_Plane_Wave), 
				is_xtl(false), islice(0), dp_Shift(false) {};

			In_Multem(const In_Multem<T>& in_multem)
			{
				assign(in_multem);
			}

			template <class TIn_Multislice>
			void assign(TIn_Multislice &in_multem)
			{
				system_config = in_multem.system_config;

				interaction_model = in_multem.interaction_model;
				pot_parm_typ = in_multem.pot_parm_typ;

				operation_mode = in_multem.operation_mode;
				slice_storage = in_multem.slice_storage;
				reverse_multislice = in_multem.reverse_multislice;
				mul_sign = in_multem.mul_sign;
				Vrl = in_multem.Vrl;
				nR = in_multem.nR;

				phonon_par = in_multem.phonon_par;

				atoms = in_multem.atoms;
				is_xtl = in_multem.is_xtl;

				rot_par = in_multem.rot_par;

				thick_type = in_multem.thick_type;
				thick = in_multem.thick;

				potential_slicing = in_multem.potential_slicing;

				grid_2d = in_multem.grid_2d;

				simulation_type = in_multem.simulation_type;

				iw_type = in_multem.iw_type;
				iw_psi = in_multem.iw_psi;

				beam_pos_i = in_multem.beam_pos_i;

				E_0 = in_multem.E_0;
				theta = in_multem.theta;
				phi = in_multem.phi;
				nrot = in_multem.nrot;

				illumination_model = in_multem.illumination_model;
				temporal_spatial_incoh = in_multem.temporal_spatial_incoh;

				cond_lens = in_multem.cond_lens;
				obj_lens = in_multem.obj_lens;

				scanning = in_multem.scanning;
				detector = in_multem.detector;

				eels_fr = in_multem.eels_fr;

				cdl_var_type = in_multem.cdl_var_type;
				cdl_var = in_multem.cdl_var;

				beam_pos = in_multem.beam_pos;

				islice = in_multem.islice;
				dp_Shift = in_multem.dp_Shift;
			}

			template <class TIn_Multislice>
			In_Multem<T>& operator=(const TIn_Multislice &in_multem)
			{
				assign(in_multem);
				return *this;
			}

			void validate_parameters()
			{
				phonon_par.validate_parameter();

				rot_par.validate_parameter();

				islice = max(0, islice);

				if (fcn_is_zero(Vrl))
				{
					Vrl = c_vr_min;
				}

				if (fcn_is_zero(nR))
				{
					nR = c_nR;
				}

				dp_Shift = (is_PED())?true:false;

				/************************ incident beam ****************************/
				if ((is_plane_wave()) || (is_scanning() && !scanning.is_user_define()))
				{
					beam_pos_i.resize(1);
				}

				beam_pos = beam_pos_i;

				/************************** Scanning *******************************/
 				if (is_scanning())
				{
					if (scanning.is_user_define())
					{
						scanning.R = beam_pos_i.p;
					}
				}
				else
				{
					scanning.set_default();
				}

 				scanning.set_grid();

				/***************************************************************************************/
				lambda = fcn_lambda(E_0);

				// multislice sign
				mul_sign = (reverse_multislice)?-1:1;

				theta = set_incident_angle(theta);
				nrot = max(1, nrot);
				if (!is_PED_HCTEM())
				{
					nrot = 1;
				}

				if (is_EELS_EFTEM())
				{
					phonon_par.coh_contrib = false;
					interaction_model = mt::eESIM_Multislice;
					illumination_model = eIM_Coherent;

					auto coll_angle_max = fcn_rangs_2_rad(E_0, 2*grid_2d.g_max());
					auto coll_angle = (fcn_is_zero(eels_fr.coll_angle)||(eels_fr.coll_angle<0))?coll_angle_max:eels_fr.coll_angle;
					eels_fr.set_coll_angle(coll_angle);

					if (is_EFTEMRS())
					{
						obj_lens.inner_aper_ang = 0;
						obj_lens.outer_aper_ang = eels_fr.coll_angle;
						obj_lens.set_in_data(E_0, grid_2d);
					}
				}

				if (is_EWFS_EWRS())
				{
					phonon_par.coh_contrib = true;
					illumination_model = eIM_Coherent;
				}

				// It will be changed when you include spatial incoherences (beta)
				if (is_CBED_CBEI() && !is_illu_mod_full_integration())
				{
					illumination_model = eIM_Coherent;
				}

				slice_storage = is_slice_storage();

				if (!is_multislice())
				{
					islice = 0;
					potential_slicing = ePST_Planes;
					if (is_through_slices())
					{
						thick_type = eTT_Through_Thick;
					}
					slice_storage = slice_storage || !is_whole_spec();
				}

				if (is_subslicing() && is_through_thick())
				{
					thick_type = eTT_Whole_Spec;
				}

				if (!is_subslicing())
				{
					phonon_par.dim_z = false;
				}

				if (is_spec_rot_active())
				{
					thick_type = eTT_Whole_Spec;
					// get geometric center
					if (rot_par.is_geometric_center())
					{
						rot_par.center_p = R_3d<T>(atoms.x_mean, atoms.y_mean, atoms.z_mean);
					}
					// rotate atoms
					rotate_atoms(atoms, rot_par.theta, rot_par.u0, rot_par.center_p);
					// get statistic
					atoms.get_statistic();
					// reset theta
					rot_par.theta = 0;
				}

				// match slicing with the require thickness
				Slicing<T> slicing;
				slicing.match_thickness(potential_slicing, atoms, thick_type, thick);

				// remove atoms outside the thickness range
				T ee_z = 0.1;
				T z_e = thick.back() + ee_z;

				if (atoms.z_max > z_e)
				{
					T z_0 = atoms.z_min - ee_z;
					remove_atoms_outside_z_range(atoms, z_0, z_e);
					// get statistic
					atoms.get_statistic();
				}

				/************* verify lenses parameters **************/

				if (fcn_is_zero(cond_lens.spt_inc_sigma))
				{
					temporal_spatial_incoh = eTSI_Temporal;
				}

				if (!is_illu_mod_full_integration())
				{
					cond_lens.spt_inc_rad_npts = 1;
					cond_lens.spt_inc_azm_npts = 1;
					cond_lens.tp_inc_npts = 1;

					obj_lens.spt_inc_rad_npts = 1;
					obj_lens.spt_inc_azm_npts = 1;
					obj_lens.tp_inc_npts = 1;
				}

				if (is_incoh_temporal())
				{
					cond_lens.spt_inc_rad_npts = 1;
					cond_lens.spt_inc_azm_npts = 1;

					obj_lens.spt_inc_rad_npts = 1;
					obj_lens.spt_inc_azm_npts = 1;
				}
				else if (is_incoh_spatial())
				{
					cond_lens.tp_inc_npts = 1;
					obj_lens.tp_inc_npts = 1;
				}

				// set condenser lens parameters
				cond_lens.set_in_data(E_0, grid_2d);
				if (atoms.empty())
				{
					cond_lens.zero_defocus_plane = 0;
				}
				else
				{
					cond_lens.zero_defocus_plane = cond_lens.get_zero_defocus_plane(atoms.z_min, atoms.z_max);
				}

				cond_lens.zero_defocus_type = eZDT_user_def;

				// set objetive lens parameters
				obj_lens.set_in_data(E_0, grid_2d);

				// validate output area
				validate_output_area();
			}

			/***************************************************************************************/
			void set_reverse_multislice(dt_bool rm)
			{
				reverse_multislice = rm;
				mul_sign = (reverse_multislice)?-1:1;
			}

			/***************************************************************************************/
			inline
			dt_int32 number_pn_conf()
			{
				return phonon_par.number_conf();
			}

			dt_int32 number_of_beams()
			{
				return (is_scanning())?scanning.size():beam_pos.size();
			}

			dt_bool is_multi_beam()
			{
				return number_of_beams() > 1;
			}

			Vctr<R_2d<T>, edev_cpu> extract_beam_pos(dt_int32 idx, dt_int32 n_idx)
			{
				dt_int32 n_beams_tot = number_of_beams();
				dt_int32 n_beams = (n_beams_tot+n_idx-1)/n_idx;

				dt_int32 ix_0 = idx*n_beams;
				dt_int32 ix_e = min(ix_0+n_beams, n_beams_tot);
			
				auto it_0 = (is_scanning())?scanning.R.begin():beam_pos.p.begin();
				return Vctr<R_2d<T>, edev_cpu>(it_0+ix_0, it_0+ix_e);
			}

			void validate_output_area()
			{
				if ((is_STEM() || is_STEM_EELS()))
				{
					if (scanning.is_user_define())
					{
						output_area.ix_0 = 0;
						output_area.ix_e = scanning.R.size();

						output_area.iy_0 = 0;
						output_area.iy_e = 1;
					}
					else
					{
						output_area.ix_0 = 0;
						output_area.ix_e = scanning.nx;

						output_area.iy_0 = 0;
						output_area.iy_e = scanning.ny;
					}
				}
				else
				{
					if (output_area.ix_0==output_area.ix_e)
					{
						output_area.ix_0 = 0;
						output_area.ix_e = grid_2d.nx;
					}

					if (output_area.iy_0==output_area.iy_e)
					{
						output_area.iy_0 = 0;
						output_area.iy_e = grid_2d.ny;
					}

					output_area.set_ind_asc_order();
					output_area.apply_bound_ix(0, grid_2d.nx);
					output_area.apply_bound_iy(0, grid_2d.ny);
				}

				output_area.ind_0 = 0;
				output_area.ind_e = output_area.size();
			}

			void set_iscan_beam_position()
			{
				// beam_pos_2d.resize(beam_pos_2d_i.size());
				// for(auto is = 0; is < beam_pos_2d_i.size(); is++)
				// {
				// 	auto idx = beam_pos_2d_i.idx[is];
				// 	beam_pos_2d[is] = R_2d<T>(scanning.x[idx], scanning.y[idx]);
				// }
			}

			/***************************************************************************************/
			dt_bool is_spec_rot_active() const
			{
				return rot_par.is_rot_active();
			}

			/***************************************************************************************/
			T Rx_exp_factor()
			{
				return 0;
				// return grid_2d.factor_2pi_rx_ctr(beam_x);
			}

			T Ry_exp_factor()
			{
				return 0;
				// return grid_2d.factor_2pi_ry_ctr(beam_y);
			}

			T set_incident_angle(const T& theta) const
			{
				T n = ::round(sin(theta)/(lambda*grid_2d.dg_min()));
				return (fcn_is_zero(theta))?0:asin(n*lambda*grid_2d.dg_min());
			}

			T get_phonon_rot_weight() const
			{
				return 1.0/static_cast<T>(number_pn_conf()*nrot);
			}

			void set_phi(const dt_int32& irot)
			{
				phi = (irot == 0)?0.0:(c_2pi<T>*static_cast<T>(irot))/static_cast<T>(nrot);
			}

			inline
			T get_propagator_factor(const T& z) const
			{
				return (-mul_sign*c_pi<T>*lambda*z);
			}

			T Vr_factor() const
			{
				return (mul_sign*fcn_transf_exp_factor(E_0, theta));
			}

			T gx_0() const
			{
				return sin(theta)*cos(phi)/lambda;
			}

			T gy_0() const
			{
				return sin(theta)*sin(phi)/lambda;
			}

			R_2d<T> gu_0() const
			{
				return R_2d<T>(gx_0(), gy_0());
			}

			void set_eels_fr_atom(const dt_int32& iatoms, const Ptc_Atom<T>& atoms)
			{
				eels_fr.x = grid_2d.factor_2pi_rx_ctr(atoms.x[iatoms]);
				eels_fr.y = grid_2d.factor_2pi_ry_ctr(atoms.y[iatoms]);
				eels_fr.occ = atoms.occ[iatoms];
			}

			/***************************************************************************************/
			dt_bool is_multislice() const
			{
				return mt::is_multislice(interaction_model);
			}

			dt_bool is_phase_object() const
			{
				return mt::is_phase_object(interaction_model);
			}

			dt_bool is_weak_phase_object() const
			{
				return mt::is_weak_phase_object(interaction_model);
			}

			/***************************************************************************************/

			dt_bool is_still_atom() const
			{
				return phonon_par.is_still_atom();
			}

			dt_bool is_absorptive_model() const
			{
				return phonon_par.is_absorptive_model();
			}

			dt_bool is_frozen_phonon() const
			{
				return phonon_par.is_frozen_phonon();
			}

			dt_bool is_frozen_phonon_float32_conf() const
			{
				return phonon_par.is_frozen_phonon_float32_conf();
			}

			/***************************************************************************************/
			dt_bool is_whole_spec() const
			{
				return mt::is_whole_spec(thick_type);
			}

			dt_bool is_through_slices() const
			{
				return mt::is_through_slices(thick_type);
			}

			dt_bool is_through_thick() const
			{
				return mt::is_through_thick(thick_type);
			}

			/***************************************************************************************/
			dt_bool is_slicing_by_planes() const
			{
				return mt::is_slicing_by_planes(interaction_model, potential_slicing);
			}

			dt_bool is_slicing_by_dz() const
			{
				return mt::is_slicing_by_dz(interaction_model, potential_slicing);
			}

			dt_bool is_subslicing() const
			{
				return mt::is_subslicing(interaction_model, potential_slicing);
			}

			dt_bool is_subslicing_whole_spec() const
			{
				return mt::is_subslicing_whole_spec(interaction_model, potential_slicing, thick_type);
			}

			/***************************************************************************************/
			dt_bool is_plane_wave() const
			{
				return mt::is_plane_wave(iw_type);
			}

			dt_bool is_convergent_wave() const
			{
				return mt::is_convergent_wave(iw_type);
			}

			dt_bool is_user_define_wave() const
			{
				return mt::is_user_define_wave(iw_type);
			}

			/***************************************************************************************/
			dt_bool is_STEM() const
			{
				return mt::is_STEM(simulation_type);
			}

			dt_bool is_ISTEM() const
			{
				return mt::is_ISTEM(simulation_type);
			}

			dt_bool is_CBED() const
			{
				return mt::is_CBED(simulation_type);
			}

			dt_bool is_CBEI() const
			{
				return mt::is_CBEI(simulation_type);
			}

			dt_bool is_ED() const
			{
				return mt::is_ED(simulation_type);
			}

			dt_bool is_HRTEM() const
			{
				return mt::is_HRTEM(simulation_type);
			}

			dt_bool is_PED() const
			{
				return mt::is_PED(simulation_type);
			}

			dt_bool is_HCTEM() const
			{
				return mt::is_HCTEM(simulation_type);
			}

			dt_bool is_EWFS() const
			{
				return mt::is_EWFS(simulation_type);
			}

			dt_bool is_EWRS() const
			{
				return mt::is_EWRS(simulation_type);
			}

			dt_bool is_EWFS_SC() const
			{
				return mt::is_EWFS_SC(simulation_type, phonon_par.model, phonon_par.single_conf);
			}

			dt_bool is_EWRS_SC() const
			{
				return mt::is_EWRS_SC(simulation_type, phonon_par.model, phonon_par.single_conf);
			}

 			dt_bool is_STEM_EELS() const
			{
				return mt::is_STEM_EELS(simulation_type);
			}

 			dt_bool is_ISTEM_EELS() const
			{
				return mt::is_ISTEM_EELS(simulation_type);
			}

			dt_bool is_STEM_ISTEM_EELS() const
			{
				return mt::is_STEM_ISTEM_EELS(simulation_type);
			}

 			dt_bool is_EFTEMFS() const
			{
				return mt::is_EFTEMFS(simulation_type);
			}

 			dt_bool is_EFTEMRS() const
			{
				return mt::is_EFTEMRS(simulation_type);
			}

 			dt_bool is_EFTEM() const
			{
				return mt::is_EFTEM(simulation_type);
			}

			dt_bool is_IWFS() const
			{
				return mt::is_IWFS(simulation_type);
			}

			dt_bool is_IWRS() const
			{
				return mt::is_IWRS(simulation_type);
			}

			dt_bool is_PPFS() const
			{
				return mt::is_PPFS(simulation_type);
			}

			dt_bool is_PPRS() const
			{
				return mt::is_PPRS(simulation_type);
			}

			dt_bool is_TFFS() const
			{
				return mt::is_TFFS(simulation_type);
			}

			dt_bool is_TFRS() const
			{
				return mt::is_TFRS(simulation_type);
			}

			dt_bool is_PropFS() const
			{
				return mt::is_PropFS(simulation_type);
			}

			dt_bool is_PropRS() const
			{
				return mt::is_PropRS(simulation_type);
			}

			dt_bool is_STEM_ISTEM() const
			{
				return mt::is_STEM_ISTEM(simulation_type);
			}

			dt_bool is_CBED_CBEI() const
			{
				return mt::is_CBED_CBEI(simulation_type);
			}

			dt_bool is_ED_HRTEM() const
			{
				return mt::is_ED_HRTEM(simulation_type);
			}

			dt_bool is_PED_HCTEM() const
			{
				return mt::is_PED_HCTEM(simulation_type);
			}

			dt_bool is_EWFS_EWRS() const
			{
				return mt::is_EWFS_EWRS(simulation_type);
			}

			dt_bool is_EWFS_EWRS_SC() const
			{
				return mt::is_EWFS_EWRS_SC(simulation_type, phonon_par.model, phonon_par.single_conf);
			}

			dt_bool is_EELS_EFTEM() const
			{
				return mt::is_EELS_EFTEM(simulation_type);
			}

			dt_bool is_IWFS_IWRS() const
			{
				return mt::is_IWFS_IWRS(simulation_type);
			}

			dt_bool is_PPFS_PPRS() const
			{
				return mt::is_PPFS_PPRS(simulation_type);
			}

			dt_bool is_TFFS_TFRS() const
			{
				return mt::is_TFFS_TFRS(simulation_type);
			}

			dt_bool is_PropFS_PropRS() const
			{
				return mt::is_PropFS_PropRS(simulation_type);
			}

			dt_bool is_grid_FS() const
			{
				return mt::is_grid_FS(simulation_type);
			}

			dt_bool is_grid_RS() const
			{
				return mt::is_grid_RS(simulation_type);
			}

			dt_bool is_simulation_type_FS() const
			{
				return mt::is_simulation_type_FS(simulation_type);
			}

			dt_bool is_simulation_type_RS() const
			{
				return mt::is_simulation_type_RS(simulation_type);
			}

			dt_bool is_specimen_required() const
			{
				return mt::is_specimen_required(simulation_type);
			}

			dt_bool is_ISTEM_CBEI_HRTEM_HCTEM_EFTEMRS() const
			{
				return mt::is_ISTEM_CBEI_HRTEM_HCTEM_EFTEMRS(simulation_type);
			}

			dt_bool is_CBED_ED_EWFS_PED() const
			{
				return mt::is_CBED_ED_EWFS_PED(simulation_type);
			}

			dt_bool is_obj_lens_temp_spat() const
			{
				return mt::is_obj_lens_temp_spat(simulation_type);
			}

			dt_bool is_cond_lens_temp_spat() const
			{
				return mt::is_cond_lens_temp_spat(simulation_type);
			}

			eSpace get_simulation_space() const
			{
				return mt::get_simulation_space(simulation_type);
			}

			dt_bool is_scanning() const
			{
				return mt::is_scanning(simulation_type);
			}

			/***************************************************************************************/
			void set_incident_wave_type(eIncident_Wave_Typ iw_type_i)
			{
				iw_type = mt::validate_incident_wave_type(simulation_type, iw_type_i);
			}

			/***************************************************************************************/
			dt_bool is_illu_mod_coherent() const
			{
				return illumination_model == eIM_Coherent;
			}

			dt_bool is_illu_mod_partial_coherent() const
			{
				return illumination_model == eIM_Partial_Coherent;
			}

			dt_bool is_illu_mod_trans_cross_coef() const
			{
				return illumination_model == eIM_Trans_Cross_Coef;
			}

			dt_bool is_illu_mod_full_integration() const
			{
				return illumination_model == eIM_Full_Integration;
			}

			/***************************************************************************************/
			dt_bool is_incoh_temporal_spatial() const
			{
				return temporal_spatial_incoh == eTSI_Temporal_Spatial;
			}

			dt_bool is_incoh_temporal() const
			{
				return temporal_spatial_incoh == eTSI_Temporal;
			}

			dt_bool is_incoh_spatial() const
			{
				return temporal_spatial_incoh == eTSI_Spatial;
			}

			/***************************************************************************************/
			dt_bool is_detector_circular() const
			{
				return detector.is_detector_circular();
			}

			dt_bool is_detector_radial() const
			{
				return detector.is_detector_radial();
			}

			dt_bool is_detector_matrix() const
			{
				return detector.is_detector_matrix();
			}

			/***************************************************************************************/
			dt_bool is_slice_storage() const
			{
				dt_bool slice_storage = is_PED_HCTEM() || is_EELS_EFTEM();
				slice_storage = slice_storage || is_ISTEM() || (is_STEM() && !phonon_par.coh_contrib);
				slice_storage = slice_storage || (is_CBED_CBEI() && is_illu_mod_full_integration());

				return slice_storage;
			}

			dt_bool is_operation_mode_normal() const
			{
				return operation_mode == eOM_Normal;
			}

			dt_bool is_operation_mode_advanced() const
			{
				return operation_mode == eOM_Advanced;
			}

			/***************************************************************************************/
			dt_bool is_lvt_off() const
			{
				return cdl_var == eLVT_off;
			}

			dt_bool is_lvt_m() const
			{
				return cdl_var == eLVT_m;
			}

			dt_bool is_lvt_Cs3() const
			{
				return cdl_var == eLVT_Cs3;
			}

			dt_bool is_lvt_Cs5() const
			{
				return cdl_var == eLVT_Cs5;
			}

			dt_bool is_lvt_mfa2() const
			{
				return cdl_var == eLVT_mfa2;
			}

			dt_bool is_lvt_afa2() const
			{
				return cdl_var == eLVT_afa2;
			}

			dt_bool is_lvt_mfa3() const
			{
				return cdl_var == eLVT_mfa3;
			}

			dt_bool is_lvt_afa3() const
			{
				return cdl_var == eLVT_afa3;
			}

			dt_bool is_lvt_inner_aper_ang() const
			{
				return cdl_var == eLVT_inner_aper_ang;
			}

			dt_bool is_lvt_outer_aper_ang() const
			{
				return cdl_var == eLVT_outer_aper_ang;
			}

		};
	}

#endif