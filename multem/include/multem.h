/*
 * This file is part of MULTEM.
 * Copyright 2020 Ivan Lobato <Ivanlh20@gmail.com>
 * Copyright 2021 Diamond Light Source
 * Copyright 2021 Rosalind Franklin Institute
 *
 * Author: James Parkhurst
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

#ifndef MULTEM_H
#define MULTEM_H

#include <array>
#include <cassert>
#include <complex>
#include <string>
#include <vector>
#include <memory>
#include <multem/constants.h>
#include <lin_alg_def.cuh>
#include <safe_types.cuh>
#include <atom_data_api.h>
#include <input_multislice_api.h>

#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_DLL
    #ifdef __GNUC__
      #define DLL_PUBLIC __attribute__ ((dllexport))
    #else
      #define DLL_PUBLIC __declspec(dllexport)
    #endif
  #else
    #ifdef __GNUC__
      #define DLL_PUBLIC __attribute__ ((dllimport))
    #else
      #define DLL_PUBLIC __declspec(dllimport)
    #endif
  #endif
#else
  #if __GNUC__ >= 4
    #define DLL_PUBLIC __attribute__ ((visibility ("default")))
  #else
    #define DLL_PUBLIC
  #endif
#endif

namespace mt {
  
  /****************************************************************************
   * The Scanning interface
   ***************************************************************************/
  
  /**
   * A class to represent an image
   */
  template <typename T>
  class Image {
  public:

    typedef T value_type;
    typedef std::array<std::size_t, 2> shape_type;

    std::vector<value_type> data;
    shape_type shape;

    Image()
      : shape({ 0, 0 }) {}

    /**
     * Construct the image from the pointer
     * @param data_ The data pointer
     * @param shape_ The image shape (Y, X)
     */
    template <typename U>
    Image(const U *data_, shape_type shape_)
      : shape(shape_) {
      std::size_t size = shape[0]*shape[1];
      data.assign(data_, data_ + size);
    }
  };

  /**
   * A class to represent a complex image
   */
  template <typename T>
  class Image< std::complex<T> > {
  public:
    
    typedef std::complex<T> value_type;
    typedef std::array<std::size_t, 2> shape_type;

    std::vector<value_type> data;
    shape_type shape;

    Image()
      : shape({ 0, 0 }) {}

    /**
     * Construct the image from the pointer
     * @param data_ The data pointer
     * @param shape_ The image shape (Y, X)
     */
    template <typename U>
    Image(const U *data_, shape_type shape_)
      : shape(shape_) {
      std::size_t size = shape[0]*shape[1];
      data.resize(size);
      for (auto i = 0; i < size; ++i) {
        data[i] = value_type(data_[i].real(), data_[i].imag());
      }
    }
    
  };
  
  /****************************************************************************
   * The Detector interface
   ***************************************************************************/

  /**
   * Interface for Det_Int. This can be removed simply by using std::vector in
   * the Det_Int class rather than thrust::vector
   */
  template <typename T>
  class DetInt {
  public:
		using value_type = typename T::value_type;
		using size_type = std::size_t;

		size_type size() const {
			return image.size();
		}

    std::vector<T> image;
  };

  /*template <typename T>*/
  /*class Input {*/
  /*public:*/
  /*  struct Data;*/

  /*  DLL_PUBLIC Input();*/
  /*  DLL_PUBLIC Input(const Input&);*/
  /*  DLL_PUBLIC Input(Input&&);*/
  /*  DLL_PUBLIC Input(const Data&);*/
  /*  DLL_PUBLIC Input& operator=(const Input&);*/
  /*  DLL_PUBLIC Input& operator=(Input&&);*/
  /*  DLL_PUBLIC ~Input();*/

  /*  DLL_PUBLIC Data& internal();*/

  /*  DLL_PUBLIC System_Configuration& get_system_conf() const;*/
  /*  DLL_PUBLIC void set_system_conf(const System_Configuration& system_conf);*/

  /*  DLL_PUBLIC eElec_Spec_Int_Model get_interaction_model() const;*/
  /*  DLL_PUBLIC void set_interaction_model(eElec_Spec_Int_Model interaction_model);*/

  /*  DLL_PUBLIC ePotential_Type get_potential_type() const;*/
  /*  DLL_PUBLIC void set_potential_type(ePotential_Type potential_type);*/

  /*  DLL_PUBLIC ePhonon_Model get_pn_model() const;*/
  /*  DLL_PUBLIC void set_pn_model(ePhonon_Model pn_model);*/

  /*  DLL_PUBLIC bool get_pn_coh_contrib() const;*/
  /*  DLL_PUBLIC void set_pn_coh_contrib(bool pn_coh_contrib);*/

  /*  DLL_PUBLIC bool get_pn_single_conf() const;*/
  /*  DLL_PUBLIC void set_pn_single_conf(bool pn_single_conf);*/

  /*  DLL_PUBLIC FP_Dim& get_pn_dim() const;*/
  /*  DLL_PUBLIC void set_pn_dim(const FP_Dim &pn_dim);*/

  /*  DLL_PUBLIC int get_fp_dist() const;*/
  /*  DLL_PUBLIC void set_fp_dist(int fp_dist);*/

  /*  DLL_PUBLIC int get_pn_seed() const;*/
  /*  DLL_PUBLIC void set_pn_seed(int pn_seed);*/

  /*  DLL_PUBLIC int get_pn_nconf() const;*/
  /*  DLL_PUBLIC void set_pn_nconf(int pn_nconf);*/

  /*  DLL_PUBLIC int get_fp_iconf_0() const;*/
  /*  DLL_PUBLIC void set_fp_iconf_0(int fp_iconf_0);*/

  /*  DLL_PUBLIC Atom_Data<T>& get_atoms() const;*/
  /*  DLL_PUBLIC void set_atoms(const Atom_Data<T>& atoms);*/

  /*  DLL_PUBLIC bool get_is_crystal() const;*/
  /*  DLL_PUBLIC void set_is_crystal(bool is_crystal);*/

  /*  DLL_PUBLIC double get_spec_rot_theta() const;*/
  /*  DLL_PUBLIC void set_spec_rot_theta(double spec_rot_theta);*/

  /*  DLL_PUBLIC r3d<T>& get_spec_rot_u0() const;*/
  /*  DLL_PUBLIC void set_spec_rot_u0(const r3d<T>& spec_rot_u0);*/

  /*  DLL_PUBLIC eRot_Point_Type get_spec_rot_center_type() const;*/
  /*  DLL_PUBLIC void set_spec_rot_center_type(eRot_Point_Type spec_rot_center_type);*/

  /*  DLL_PUBLIC r3d<T>& get_spec_rot_center_p() const;*/
  /*  DLL_PUBLIC void set_spec_rot_center_p(const r3d<T>& spec_rot_center_p);*/

  /*  DLL_PUBLIC eThick_Type get_thick_type() const;*/
  /*  DLL_PUBLIC void set_thick_type(eThick_Type thick_type);*/

  /*  DLL_PUBLIC std::vector<T> get_thick() const;*/
  /*  DLL_PUBLIC void set_thick(const std::vector<T>& thick);*/

  /*  DLL_PUBLIC ePotential_Slicing get_potential_slicing() const;*/
  /*  DLL_PUBLIC void set_potential_slicing(ePotential_Slicing potential_slicing);*/

  /*  DLL_PUBLIC Grid_2d<T>& get_grid_2d() const;*/
  /*  DLL_PUBLIC void set_grid_2d(const Grid_2d<T>& grid_2d);*/

  /*  DLL_PUBLIC Range_2d& get_output_area() const;*/
  /*  DLL_PUBLIC void set_output_area(const Range_2d& output_area);*/

  /*  DLL_PUBLIC eTEM_Sim_Type get_simulation_type() const;*/
  /*  DLL_PUBLIC void set_simulation_type(eTEM_Sim_Type simulation_type);*/

  /*  DLL_PUBLIC eIncident_Wave_Type get_iw_type() const;*/
  /*  DLL_PUBLIC void set_iw_type(eIncident_Wave_Type iw_type);*/

  /*  DLL_PUBLIC std::vector<std::complex<T>> get_iw_psi() const;*/
  /*  DLL_PUBLIC void set_iw_psi(const std::vector<std::complex<T>>& iw_psi);*/

  /*  DLL_PUBLIC std::vector<T> get_iw_x() const;*/
  /*  DLL_PUBLIC void set_iw_x(const std::vector<T>& iw_x);*/

  /*  DLL_PUBLIC std::vector<T> get_iw_y() const;*/
  /*  DLL_PUBLIC void set_iw_y(const std::vector<T>& iw_y);*/

  /*  DLL_PUBLIC double get_E_0() const;*/
  /*  DLL_PUBLIC void set_E_0(double E_0);*/

  /*  DLL_PUBLIC double get_lambda() const;*/
  /*  DLL_PUBLIC void set_lambda(double lambda);*/

  /*  DLL_PUBLIC double get_theta() const;*/
  /*  DLL_PUBLIC void set_theta(double theta);*/

  /*  DLL_PUBLIC double get_phi() const;*/
  /*  DLL_PUBLIC void set_phi(double phi);*/

  /*  DLL_PUBLIC eIllumination_Model get_illumination_model() const;*/
  /*  DLL_PUBLIC void set_illumination_model(eIllumination_Model illumination_model);*/

  /*  DLL_PUBLIC eTemporal_Spatial_Incoh get_temporal_spatial_incoh() const;*/
  /*  DLL_PUBLIC void set_temporal_spatial_incoh(eTemporal_Spatial_Incoh temporal_spatial_incoh);*/

  /*  DLL_PUBLIC Lens<T>& get_cond_lens() const;*/
  /*  DLL_PUBLIC void set_cond_lens(const Lens<T>& cond_lens);*/

  /*  DLL_PUBLIC Lens<T>& get_obj_lens() const;*/
  /*  DLL_PUBLIC void set_obj_lens(const Lens<T>& obj_lens);*/

  /*  DLL_PUBLIC Scanning<T>& get_scanning() const;*/
  /*  DLL_PUBLIC void set_scanning(const Scanning<T>& scanning);*/

  /*  DLL_PUBLIC Detector<T, e_host>& get_detector() const;*/
  /*  DLL_PUBLIC void set_detector(const Detector<T, e_host>& detector);*/

  /*  DLL_PUBLIC EELS<T>& get_eels_fr() const;*/
  /*  DLL_PUBLIC void set_eels_fr(const EELS<T>& eels_fr);*/

  /*  DLL_PUBLIC eOperation_Mode get_operation_mode() const;*/
  /*  DLL_PUBLIC void set_operation_mode(eOperation_Mode operation_mode);*/

  /*  DLL_PUBLIC bool get_slice_storage() const;*/
  /*  DLL_PUBLIC void set_slice_storage(bool slice_storage);*/

  /*  DLL_PUBLIC bool get_reverse_multislice() const;*/
  /*  DLL_PUBLIC void set_reverse_multislice(bool reverse_multislice);*/

  /*  DLL_PUBLIC int get_mul_sign() const;*/
  /*  DLL_PUBLIC void set_mul_sign(int mul_sign);*/

  /*  DLL_PUBLIC double get_Vrl() const;*/
  /*  DLL_PUBLIC void set_Vrl(double Vrl);*/

  /*  DLL_PUBLIC int get_nR() const;*/
  /*  DLL_PUBLIC void set_nR(int nR);*/

  /*  DLL_PUBLIC int get_nrot() const;*/
  /*  DLL_PUBLIC void set_nrot(int nrot);*/

  /*  DLL_PUBLIC eLens_Var_Type get_cdl_var_type() const;*/
  /*  DLL_PUBLIC void set_cdl_var_type(eLens_Var_Type cdl_var_type);*/

  /*  DLL_PUBLIC std::vector<T> get_cdl_var() const;*/
  /*  DLL_PUBLIC void set_cdl_var(const std::vector<T>& cdl_var);*/

  /*  DLL_PUBLIC std::vector<int> get_iscan() const;*/
  /*  DLL_PUBLIC void set_iscan(const std::vector<int>& iscan);*/

  /*  DLL_PUBLIC std::vector<T> get_beam_x() const;*/
  /*  DLL_PUBLIC void set_beam_x(const std::vector<T>& beam_x);*/

  /*  DLL_PUBLIC std::vector<T> get_beam_y() const;*/
  /*  DLL_PUBLIC void set_beam_y(const std::vector<T>& beam_y);*/

  /*  DLL_PUBLIC int get_islice() const;*/
  /*  DLL_PUBLIC void set_islice(int islice);*/

  /*  DLL_PUBLIC bool get_dp_Shift() const;*/
  /*  DLL_PUBLIC void set_dp_Shift(bool dp_Shift);*/

  /*  DLL_PUBLIC void set_incident_wave_type(eIncident_Wave_Type iw_type);*/

  /*  DLL_PUBLIC void validate_parameters();*/

		/*DLL_PUBLIC bool is_multislice() const;*/
		/*DLL_PUBLIC bool is_phase_object() const;*/
		/*DLL_PUBLIC bool is_weak_phase_object() const;*/
		/*DLL_PUBLIC bool is_still_atom() const;*/
		/*DLL_PUBLIC bool is_absorptive_model() const;*/
		/*DLL_PUBLIC bool is_frozen_phonon() const;*/
		/*DLL_PUBLIC bool is_frozen_phonon_single_conf() const;*/
		/*DLL_PUBLIC bool is_whole_spec() const;*/
		/*DLL_PUBLIC bool is_through_slices() const;*/
		/*DLL_PUBLIC bool is_through_thick() const;*/
		/*DLL_PUBLIC bool is_slicing_by_planes() const;*/
		/*DLL_PUBLIC bool is_slicing_by_dz() const;*/
		/*DLL_PUBLIC bool is_subslicing() const;*/
		/*DLL_PUBLIC bool is_subslicing_whole_spec() const;*/
		/*DLL_PUBLIC bool is_plane_wave() const;*/
		/*DLL_PUBLIC bool is_convergent_wave() const;*/
		/*DLL_PUBLIC bool is_user_define_wave() const;*/
		/*DLL_PUBLIC bool is_STEM() const;*/
		/*DLL_PUBLIC bool is_ISTEM() const;*/
		/*DLL_PUBLIC bool is_CBED() const;*/
		/*DLL_PUBLIC bool is_CBEI() const;*/
		/*DLL_PUBLIC bool is_ED() const;*/
		/*DLL_PUBLIC bool is_HRTEM() const;*/
		/*DLL_PUBLIC bool is_PED() const;*/
		/*DLL_PUBLIC bool is_HCTEM() const;*/
		/*DLL_PUBLIC bool is_EWFS() const;*/
		/*DLL_PUBLIC bool is_EWRS() const;*/
		/*DLL_PUBLIC bool is_EWFS_SC() const;*/
		/*DLL_PUBLIC bool is_EWRS_SC() const;*/
		/*DLL_PUBLIC bool is_EELS() const;*/
		/*DLL_PUBLIC bool is_EFTEM() const;*/
		/*DLL_PUBLIC bool is_IWFS() const;*/
		/*DLL_PUBLIC bool is_IWRS() const;*/
		/*DLL_PUBLIC bool is_PPFS() const;*/
		/*DLL_PUBLIC bool is_PPRS() const;*/
		/*DLL_PUBLIC bool is_TFFS() const;*/
		/*DLL_PUBLIC bool is_TFRS() const;*/
		/*DLL_PUBLIC bool is_PropFS() const;*/
		/*DLL_PUBLIC bool is_PropRS() const;*/
		/*DLL_PUBLIC bool is_STEM_ISTEM() const;*/
		/*DLL_PUBLIC bool is_CBED_CBEI() const;*/
		/*DLL_PUBLIC bool is_ED_HRTEM() const;*/
		/*DLL_PUBLIC bool is_PED_HCTEM() const;*/
		/*DLL_PUBLIC bool is_EWFS_EWRS() const;*/
		/*DLL_PUBLIC bool is_EWFS_EWRS_SC() const;*/
		/*DLL_PUBLIC bool is_EELS_EFTEM() const;*/
		/*DLL_PUBLIC bool is_IWFS_IWRS() const;*/
		/*DLL_PUBLIC bool is_PPFS_PPRS() const;*/
		/*DLL_PUBLIC bool is_TFFS_TFRS() const;*/
		/*DLL_PUBLIC bool is_PropFS_PropRS() const;*/
		/*DLL_PUBLIC bool is_grid_FS() const;*/
		/*DLL_PUBLIC bool is_grid_RS() const;*/
		/*DLL_PUBLIC bool is_simulation_type_FS() const;*/
		/*DLL_PUBLIC bool is_simulation_type_RS() const;*/
		/*DLL_PUBLIC bool is_specimen_required() const;*/
		/*DLL_PUBLIC bool is_ISTEM_CBEI_HRTEM_HCTEM_EFTEM() const;*/
		/*DLL_PUBLIC bool is_CBED_ED_EWFS_PED() const;*/
		/*DLL_PUBLIC bool is_obj_lens_temp_spat() const;*/
		/*DLL_PUBLIC bool is_cond_lens_temp_spat() const;*/
		/*DLL_PUBLIC bool is_scanning() const;*/
		/*DLL_PUBLIC bool is_illu_mod_coherent() const;*/
		/*DLL_PUBLIC bool is_illu_mod_partial_coherent() const;*/
		/*DLL_PUBLIC bool is_illu_mod_trans_cross_coef() const;*/
		/*DLL_PUBLIC bool is_illu_mod_full_integration() const;*/
		/*DLL_PUBLIC bool is_incoh_temporal_spatial() const;*/
		/*DLL_PUBLIC bool is_incoh_temporal() const;*/
		/*DLL_PUBLIC bool is_incoh_spatial() const;*/
		/*DLL_PUBLIC bool is_detector_circular() const;*/
		/*DLL_PUBLIC bool is_detector_radial() const;*/
		/*DLL_PUBLIC bool is_detector_matrix() const;*/
		/*DLL_PUBLIC bool is_slice_storage() const;*/
		/*DLL_PUBLIC bool is_operation_mode_normal() const;*/
		/*DLL_PUBLIC bool is_operation_mode_advanced() const;*/
		/*DLL_PUBLIC bool is_lvt_off() const;*/
		/*DLL_PUBLIC bool is_lvt_m() const;*/
		/*DLL_PUBLIC bool is_lvt_Cs3() const;*/
		/*DLL_PUBLIC bool is_lvt_Cs5() const;*/
		/*DLL_PUBLIC bool is_lvt_mfa2() const;*/
		/*DLL_PUBLIC bool is_lvt_afa2() const;*/
		/*DLL_PUBLIC bool is_lvt_mfa3() const;*/
		/*DLL_PUBLIC bool is_lvt_afa3() const;*/
		/*DLL_PUBLIC bool is_lvt_inner_aper_ang() const;*/
		/*DLL_PUBLIC bool is_lvt_outer_aper_ang() const;*/

  /*private:*/
  /*  std::unique_ptr<Data> impl_;*/
  /*};*/
  
  /****************************************************************************
   * The Output interface
   ***************************************************************************/

  /**
   * A proxy for the Output_Multislice class
   */
  template <typename T>
  class Output : public Input_Multislice<T> {
  public:
    struct Data;

    typedef std::vector<T> vector_type;
    typedef std::vector<std::complex<T>> complex_vector_type;

    DLL_PUBLIC Output();
    DLL_PUBLIC Output(Output&);
    DLL_PUBLIC Output(Output&&);
    DLL_PUBLIC Output(Data&);
    DLL_PUBLIC Output& operator=(Output&);
    DLL_PUBLIC Output& operator=(Output&&);
    DLL_PUBLIC ~Output();

    DLL_PUBLIC const Data& internal() const;
    DLL_PUBLIC Data& internal();

    DLL_PUBLIC void set_input_data(Input_Multislice<T> &input);

    DLL_PUBLIC eTEM_Output_Type get_output_type() const;
    DLL_PUBLIC void set_output_type(eTEM_Output_Type output_type);

    DLL_PUBLIC int get_ndetector() const;
		DLL_PUBLIC void set_ndetector(int ndetector);

    DLL_PUBLIC int get_nx() const;
    DLL_PUBLIC void set_nx(int nx);

    DLL_PUBLIC int get_ny() const;
    DLL_PUBLIC void set_ny(int ny);

    DLL_PUBLIC T get_dx() const;
    DLL_PUBLIC void set_dx(T dx);

    DLL_PUBLIC T get_dy() const;
    DLL_PUBLIC void set_dy(T dy);

    DLL_PUBLIC T get_dr() const;
		DLL_PUBLIC void set_dr(T dr);

    DLL_PUBLIC std::vector<T> get_x() const;
    DLL_PUBLIC void set_x(const std::vector<T>& x);

    DLL_PUBLIC std::vector<T> get_y() const;
    DLL_PUBLIC void set_y(const std::vector<T>& y);

    DLL_PUBLIC std::vector<T> get_r() const;
    DLL_PUBLIC void set_r(const std::vector<T>& r);

    DLL_PUBLIC std::vector<DetInt<vector_type>> get_image_tot() const;
    DLL_PUBLIC void set_image_tot(const std::vector<DetInt<vector_type>>& image_tot);

    DLL_PUBLIC std::vector<DetInt<vector_type>> get_image_coh() const;
    DLL_PUBLIC void set_image_coh(const std::vector<DetInt<vector_type>>& image_coh);

    DLL_PUBLIC std::vector<vector_type> get_m2psi_tot() const;
    DLL_PUBLIC void set_m2psi_tot(const std::vector<vector_type>& m2psi_tot);

    DLL_PUBLIC std::vector<vector_type> get_m2psi_coh() const;
    DLL_PUBLIC void set_m2psi_coh(const std::vector<vector_type>& m2psi_coh);

    DLL_PUBLIC std::vector<complex_vector_type> get_psi_coh() const;
    DLL_PUBLIC void set_psi_coh(const std::vector<complex_vector_type>& psi_coh);

    DLL_PUBLIC std::vector<vector_type> get_V() const;
    DLL_PUBLIC void set_V(const std::vector<vector_type>& V);

    DLL_PUBLIC std::vector<complex_vector_type> get_trans() const;
    DLL_PUBLIC void set_trans(const std::vector<complex_vector_type>& trans);

    DLL_PUBLIC std::vector<complex_vector_type> get_psi_0() const;
    DLL_PUBLIC void set_psi_0(const std::vector<complex_vector_type>& psi_0);

		DLL_PUBLIC std::vector<bool> get_thk_gpu() const;
		DLL_PUBLIC void set_thk_gpu(const std::vector<bool>& thk_gpu);

    DLL_PUBLIC void gather();
    DLL_PUBLIC void clean_temporal();
  
  private:
    std::unique_ptr<Data> impl_;

  };
 
  /****************************************************************************
   * Logging functions
   ***************************************************************************/

  DLL_PUBLIC std::string to_string(const System_Configuration &sys_conf, const std::string &prefix="");

  template <typename T>
  DLL_PUBLIC std::string to_string(const Lens<T> &lens, const std::string &prefix="");

  template <typename T>
  DLL_PUBLIC std::string to_string(const Input_Multislice<T> &input, const std::string &prefix="");

  /****************************************************************************
   * Simulation functions
   ***************************************************************************/

  template <typename T>
  DLL_PUBLIC mt::Output<T> tem_simulation(Input_Multislice<T>&);

}  // namespace mt

#endif
