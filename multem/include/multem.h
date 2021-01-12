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

namespace mt {

  /****************************************************************************
   * The Atom_Data interface
   ***************************************************************************/

  template <typename T>
  class Atom {
  public:
    int Z;
    T x;
    T y;
    T z;
    T sigma;
    T occ;
    int region;
    int charge;

    Atom()
      : Z(0),
        x(0),
        y(0),
        z(0),
        sigma(0),
        occ(0),
        region(0),
        charge(0) {}

    Atom(int Z_in,
         T x_in,
         T y_in,
         T z_in,
         T sigma_in,
         T occ_in,
         T region_in,
         T charge_in)
      : Z(Z_in),
        x(x_in),
        y(y_in),
        z(z_in),
        sigma(sigma_in),
        occ(occ_in),
        region(region_in),
        charge(charge_in) {}
  };

  /**
   * A proxy for the Atom_Data class
   */
  template <typename T>
  class AtomData {
  public:
    struct Data;

    typedef std::size_t size_type;

    /* AtomData(); */
    AtomData(const AtomData&);
    AtomData(AtomData&&);
    AtomData(const Data&);
    AtomData& operator=(const AtomData&);
    AtomData& operator=(AtomData&&);
    ~AtomData();

    const Data& internal() const;

    T get_dz() const;
    void set_dz(T dz);
    
    T get_l_x() const;
    void set_l_x(T l_x);
    
    T get_l_y() const;
    void set_l_y(T l_y);
    
    T get_l_z() const;
    void set_l_z(T l_z);
    
    int get_ct_na() const;
    void set_ct_na(int ct_na);
    
    int get_ct_nb() const;
    void set_ct_nb(int ct_nb);
    
    int get_ct_nc() const;
    void set_ct_nc(int ct_nc);
    
    T get_ct_a() const;
    void set_ct_a(T ct_a);
    
    T get_ct_b() const;
    void set_ct_b(T ct_b);
    
    T get_ct_c() const;
    void set_ct_c(T ct_c);
    
    T get_ct_x0() const;
    void set_ct_x0(T ct_x0);
    
    T get_ct_y0() const;
    void set_ct_y0(T ct_y0);

    std::vector<Amorp_Lay_Info<T>> get_amorphous_parameters() const;
    void set_amorphous_parameters(const std::vector<Amorp_Lay_Info<T>> &amorp_lay_info);
    
    std::vector<Atom<T>> get_spec_atoms() const;
    void set_spec_atoms(const std::vector<Atom<T>>& spec_atoms);
  
    void get_statistic();
    void clear();
    bool empty() const;
    size_type size() const;
  
  private:
    std::unique_ptr<Data> impl_;
  };
  
  /****************************************************************************
   * The Scanning interface
   ***************************************************************************/

  /**
   * A proxy for the Scanning class
   */
  template <typename T>
  class ScanningData {
  public:
    struct Data;

    ScanningData(const ScanningData&);
    ScanningData(ScanningData&&);
    ScanningData(const Data&);
    ScanningData& operator=(const ScanningData&);
    ScanningData& operator=(ScanningData&&);
    ~ScanningData();

    const Data& internal() const;

    eScanning_Type get_type() const;
    void set_type(eScanning_Type type);
    
    bool get_pbc() const;
    void set_pbc(bool pbc);
    
    bool get_spxs() const;
    void set_spxs(bool spxs);
    
    int get_ns() const;
    void set_ns(int ns);

    T get_x0() const;
    void set_x0(T x0);
    
    T get_y0() const;
    void set_y0(T y0);
    
    T get_xe() const;
    void set_xe(T xe);
    
    T get_ye() const;
    void set_ye(T ye);

    void set_grid();

  private:
    std::unique_ptr<Data> impl_;
  };
  
  /****************************************************************************
   * The Scanning interface
   ***************************************************************************/

  /**
   * A proxy for the Scanning class
   */
  template <typename T>
  class DetectorData {
  public:
    struct Data;
		
    using value_type = T;
		using size_type = std::size_t;

    DetectorData(const DetectorData&);
    DetectorData(DetectorData&&);
    DetectorData(const Data&);
    DetectorData& operator=(const DetectorData&);
    DetectorData& operator=(DetectorData&&);
    ~DetectorData();

    const Data& internal() const;

		size_type size() const;
		void clear();
		void resize(const size_type &new_size);
		bool is_detector_circular() const;
		bool is_detector_radial() const;
		bool is_detector_matrix() const;
		
    eDetector_Type get_type() const;
    void set_type(eDetector_Type type);

    std::vector<T> get_g_inner() const;
    void set_g_inner(const std::vector<T>& g_inner);

    std::vector<T> get_g_outer() const;
    void set_g_outer(const std::vector<T>& g_outer);

    std::vector<std::vector<T>> get_fx() const;
    void set_fx(const std::vector<std::vector<T>>& fx);
    
    std::vector<std::vector<T>> get_fR() const;
    void set_fR(const std::vector<std::vector<T>>& fR);
		
    std::vector<Grid_2d<T>> get_grid_1d() const;
		void set_grid_1d(const std::vector<Grid_2d<T>>& grid_1d);
		
    std::vector<Grid_2d<T>> get_grid_2d() const;
		void set_grid_2d(const std::vector<Grid_2d<T>>& grid_2d);
		
    std::vector<std::string> get_fn() const;
		void set_fn(const std::vector<std::string>& fn);

    std::vector<T> get_inner_ang() const;
    void set_inner_ang(const std::vector<T>& inner_ang);

    std::vector<T> get_outer_ang() const;
    void set_outer_ang(const std::vector<T>& outer_ang);

  private:
    std::unique_ptr<Data> impl_;
  };
  
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

  /****************************************************************************
   * The Input interface
   ***************************************************************************/

  /**
   * A proxy for the Input_Multislice class
   */
  template <typename T>
  class Input {
  public:
    struct Data;

    Input();
    Input(const Input&);
    Input(Input&&);
    Input(const Data&);
    Input& operator=(const Input&);
    Input& operator=(Input&&);
    ~Input();

    Data& internal();

    System_Configuration& get_system_conf() const;
    void set_system_conf(const System_Configuration& system_conf);

    eElec_Spec_Int_Model get_interaction_model() const;
    void set_interaction_model(eElec_Spec_Int_Model interaction_model);

    ePotential_Type get_potential_type() const;
    void set_potential_type(ePotential_Type potential_type);

    ePhonon_Model get_pn_model() const;
    void set_pn_model(ePhonon_Model pn_model);

    bool get_pn_coh_contrib() const;
    void set_pn_coh_contrib(bool pn_coh_contrib);

    bool get_pn_single_conf() const;
    void set_pn_single_conf(bool pn_single_conf);

    FP_Dim& get_pn_dim() const;
    void set_pn_dim(const FP_Dim &pn_dim);

    int get_fp_dist() const;
    void set_fp_dist(int fp_dist);

    int get_pn_seed() const;
    void set_pn_seed(int pn_seed);

    int get_pn_nconf() const;
    void set_pn_nconf(int pn_nconf);

    int get_fp_iconf_0() const;
    void set_fp_iconf_0(int fp_iconf_0);

    AtomData<T>& get_atoms() const;
    void set_atoms(const AtomData<T>& atoms);

    bool get_is_crystal() const;
    void set_is_crystal(bool is_crystal);

    double get_spec_rot_theta() const;
    void set_spec_rot_theta(double spec_rot_theta);

    r3d<T>& get_spec_rot_u0() const;
    void set_spec_rot_u0(const r3d<T>& spec_rot_u0);

    eRot_Point_Type get_spec_rot_center_type() const;
    void set_spec_rot_center_type(eRot_Point_Type spec_rot_center_type);

    r3d<T>& get_spec_rot_center_p() const;
    void set_spec_rot_center_p(const r3d<T>& spec_rot_center_p);

    eThick_Type get_thick_type() const;
    void set_thick_type(eThick_Type thick_type);

    std::vector<T> get_thick() const;
    void set_thick(const std::vector<T>& thick);

    ePotential_Slicing get_potential_slicing() const;
    void set_potential_slicing(ePotential_Slicing potential_slicing);

    Grid_2d<T>& get_grid_2d() const;
    void set_grid_2d(const Grid_2d<T>& grid_2d);

    Range_2d& get_output_area() const;
    void set_output_area(const Range_2d& output_area);

    eTEM_Sim_Type get_simulation_type() const;
    void set_simulation_type(eTEM_Sim_Type simulation_type);

    eIncident_Wave_Type get_iw_type() const;
    void set_iw_type(eIncident_Wave_Type iw_type);

    std::vector<std::complex<T>> get_iw_psi() const;
    void set_iw_psi(const std::vector<std::complex<T>>& iw_psi);

    std::vector<T> get_iw_x() const;
    void set_iw_x(const std::vector<T>& iw_x);

    std::vector<T> get_iw_y() const;
    void set_iw_y(const std::vector<T>& iw_y);

    double get_E_0() const;
    void set_E_0(double E_0);

    double get_lambda() const;
    void set_lambda(double lambda);

    double get_theta() const;
    void set_theta(double theta);

    double get_phi() const;
    void set_phi(double phi);

    eIllumination_Model get_illumination_model() const;
    void set_illumination_model(eIllumination_Model illumination_model);

    eTemporal_Spatial_Incoh get_temporal_spatial_incoh() const;
    void set_temporal_spatial_incoh(eTemporal_Spatial_Incoh temporal_spatial_incoh);

    Lens<T>& get_cond_lens() const;
    void set_cond_lens(const Lens<T>& cond_lens);

    Lens<T>& get_obj_lens() const;
    void set_obj_lens(const Lens<T>& obj_lens);

    ScanningData<T>& get_scanning() const;
    void set_scanning(const ScanningData<T>& scanning);

    DetectorData<T>& get_detector() const;
    void set_detector(const DetectorData<T>& detector);

    EELS<T>& get_eels_fr() const;
    void set_eels_fr(const EELS<T>& eels_fr);

    eOperation_Mode get_operation_mode() const;
    void set_operation_mode(eOperation_Mode operation_mode);

    bool get_slice_storage() const;
    void set_slice_storage(bool slice_storage);

    bool get_reverse_multislice() const;
    void set_reverse_multislice(bool reverse_multislice);

    int get_mul_sign() const;
    void set_mul_sign(int mul_sign);

    double get_Vrl() const;
    void set_Vrl(double Vrl);

    int get_nR() const;
    void set_nR(int nR);

    int get_nrot() const;
    void set_nrot(int nrot);

    eLens_Var_Type get_cdl_var_type() const;
    void set_cdl_var_type(eLens_Var_Type cdl_var_type);

    std::vector<T> get_cdl_var() const;
    void set_cdl_var(const std::vector<T>& cdl_var);

    std::vector<int> get_iscan() const;
    void set_iscan(const std::vector<int>& iscan);

    std::vector<T> get_beam_x() const;
    void set_beam_x(const std::vector<T>& beam_x);

    std::vector<T> get_beam_y() const;
    void set_beam_y(const std::vector<T>& beam_y);

    int get_islice() const;
    void set_islice(int islice);

    bool get_dp_Shift() const;
    void set_dp_Shift(bool dp_Shift);

    void set_incident_wave_type(eIncident_Wave_Type iw_type);

    void validate_parameters();

		bool is_multislice() const;
		bool is_phase_object() const;
		bool is_weak_phase_object() const;
		bool is_still_atom() const;
		bool is_absorptive_model() const;
		bool is_frozen_phonon() const;
		bool is_frozen_phonon_single_conf() const;
		bool is_whole_spec() const;
		bool is_through_slices() const;
		bool is_through_thick() const;
		bool is_slicing_by_planes() const;
		bool is_slicing_by_dz() const;
		bool is_subslicing() const;
		bool is_subslicing_whole_spec() const;
		bool is_plane_wave() const;
		bool is_convergent_wave() const;
		bool is_user_define_wave() const;
		bool is_STEM() const;
		bool is_ISTEM() const;
		bool is_CBED() const;
		bool is_CBEI() const;
		bool is_ED() const;
		bool is_HRTEM() const;
		bool is_PED() const;
		bool is_HCTEM() const;
		bool is_EWFS() const;
		bool is_EWRS() const;
		bool is_EWFS_SC() const;
		bool is_EWRS_SC() const;
		bool is_EELS() const;
		bool is_EFTEM() const;
		bool is_IWFS() const;
		bool is_IWRS() const;
		bool is_PPFS() const;
		bool is_PPRS() const;
		bool is_TFFS() const;
		bool is_TFRS() const;
		bool is_PropFS() const;
		bool is_PropRS() const;
		bool is_STEM_ISTEM() const;
		bool is_CBED_CBEI() const;
		bool is_ED_HRTEM() const;
		bool is_PED_HCTEM() const;
		bool is_EWFS_EWRS() const;
		bool is_EWFS_EWRS_SC() const;
		bool is_EELS_EFTEM() const;
		bool is_IWFS_IWRS() const;
		bool is_PPFS_PPRS() const;
		bool is_TFFS_TFRS() const;
		bool is_PropFS_PropRS() const;
		bool is_grid_FS() const;
		bool is_grid_RS() const;
		bool is_simulation_type_FS() const;
		bool is_simulation_type_RS() const;
		bool is_specimen_required() const;
		bool is_ISTEM_CBEI_HRTEM_HCTEM_EFTEM() const;
		bool is_CBED_ED_EWFS_PED() const;
		bool is_obj_lens_temp_spat() const;
		bool is_cond_lens_temp_spat() const;
		bool is_scanning() const;
		bool is_illu_mod_coherent() const;
		bool is_illu_mod_partial_coherent() const;
		bool is_illu_mod_trans_cross_coef() const;
		bool is_illu_mod_full_integration() const;
		bool is_incoh_temporal_spatial() const;
		bool is_incoh_temporal() const;
		bool is_incoh_spatial() const;
		bool is_detector_circular() const;
		bool is_detector_radial() const;
		bool is_detector_matrix() const;
		bool is_slice_storage() const;
		bool is_operation_mode_normal() const;
		bool is_operation_mode_advanced() const;
		bool is_lvt_off() const;
		bool is_lvt_m() const;
		bool is_lvt_Cs3() const;
		bool is_lvt_Cs5() const;
		bool is_lvt_mfa2() const;
		bool is_lvt_afa2() const;
		bool is_lvt_mfa3() const;
		bool is_lvt_afa3() const;
		bool is_lvt_inner_aper_ang() const;
		bool is_lvt_outer_aper_ang() const;

  private:
    std::unique_ptr<Data> impl_;
  };
  
  /****************************************************************************
   * The Output interface
   ***************************************************************************/

  /**
   * A proxy for the Output_Multislice class
   */
  template <typename T>
  class Output : public Input<T> {
  public:
    struct Data;

    typedef std::vector<T> vector_type;
    typedef std::vector<std::complex<T>> complex_vector_type;

    Output();
    Output(Output&);
    Output(Output&&);
    Output(Data&);
    Output& operator=(Output&);
    Output& operator=(Output&&);
    ~Output();

    const Data& internal() const;
    Data& internal();

    void set_input_data(Input<T> &input);
		
    eTEM_Output_Type get_output_type() const;
    void set_output_type(eTEM_Output_Type output_type);
		
    int get_ndetector() const;
		void set_ndetector(int ndetector);
    
    int get_nx() const;
    void set_nx(int nx);
    
    int get_ny() const;
    void set_ny(int ny);
    
    T get_dx() const;
    void set_dx(T dx);
    
    T get_dy() const;
    void set_dy(T dy);
		
    T get_dr() const;
		void set_dr(T dr);

    std::vector<T> get_x() const;
    void set_x(const std::vector<T>& x);
    
    std::vector<T> get_y() const;
    void set_y(const std::vector<T>& y);
    
    std::vector<T> get_r() const;
    void set_r(const std::vector<T>& r);
    
    std::vector<DetInt<vector_type>> get_image_tot() const;
    void set_image_tot(const std::vector<DetInt<vector_type>>& image_tot);
    
    std::vector<DetInt<vector_type>> get_image_coh() const;
    void set_image_coh(const std::vector<DetInt<vector_type>>& image_coh);
    
    std::vector<vector_type> get_m2psi_tot() const;
    void set_m2psi_tot(const std::vector<vector_type>& m2psi_tot);
    
    std::vector<vector_type> get_m2psi_coh() const;
    void set_m2psi_coh(const std::vector<vector_type>& m2psi_coh);
    
    std::vector<complex_vector_type> get_psi_coh() const;
    void set_psi_coh(const std::vector<complex_vector_type>& psi_coh);
    
    std::vector<vector_type> get_V() const;
    void set_V(const std::vector<vector_type>& V);
    
    std::vector<complex_vector_type> get_trans() const;
    void set_trans(const std::vector<complex_vector_type>& trans);
    
    std::vector<complex_vector_type> get_psi_0() const;
    void set_psi_0(const std::vector<complex_vector_type>& psi_0);

		std::vector<bool> get_thk_gpu() const;
		void set_thk_gpu(const std::vector<bool>& thk_gpu);

    void gather();
    void clean_temporal();
  
  private:
    std::unique_ptr<Data> impl_;

  };
  
  /****************************************************************************
   * The Multislice interface
   ***************************************************************************/

  /**
   * A proxy for the Output_Multislice class
   */
  /* template <typename T> */
  /* class MultisliceData { */
  /* public: */
  /*   struct Data; */

  /*   MultisliceData(); */
  /*   ~MultisliceData(); */

  /*   const Data& internal() const; */

		/* void set_input_data(Input<T> &input); */
		/* void operator()(Output<T> &output_multislice); */
    
  /* private: */
  /*   std::unique_ptr<Data> impl_; */

  /* }; */
  
  /****************************************************************************
   * Simulation functions
   ***************************************************************************/

  template <typename T>
  mt::Output<T> tem_simulation(Input<T>&);

}  // namespace mt

#endif
