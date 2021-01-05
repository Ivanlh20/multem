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

#include <iostream>
#include <multem.h>
#include <types.cuh>
#include <input_multislice.cuh>

namespace mt {

  /****************************************************************************
   * The SystemConfiguration interface
   ***************************************************************************/

  struct SystemConfiguration::Data {
    System_Configuration data;
    Data() {}
    Data(const System_Configuration& d)
        : data(d) {}
  };

  SystemConfiguration::SystemConfiguration()
      : impl_(std::make_unique<Data>()) {}

  SystemConfiguration::SystemConfiguration(const SystemConfiguration& other)
      : impl_(std::make_unique<Data>(*other.impl_)) {}

  SystemConfiguration::SystemConfiguration(SystemConfiguration&& other) = default;

  SystemConfiguration& SystemConfiguration::operator=(const SystemConfiguration& other) {
    *impl_ = *other.impl_;
    return *this;
  }

  SystemConfiguration::SystemConfiguration(const Data& other)
      : impl_(std::make_unique<Data>(other)) {}

  SystemConfiguration& SystemConfiguration::operator=(SystemConfiguration&&) = default;

  SystemConfiguration::~SystemConfiguration() = default;

  const SystemConfiguration::Data& SystemConfiguration::internal() const {
    return *impl_;
  }

  void SystemConfiguration::set_precision(ePrecision precision) {
    impl_->data.precision = precision;
  }

  ePrecision SystemConfiguration::get_precision() const {
    return impl_->data.precision;
  }

  void SystemConfiguration::set_device(eDevice device) {
    impl_->data.device = device;
  }

  eDevice SystemConfiguration::get_device() const {
    return impl_->data.device;
  }

  void SystemConfiguration::set_cpu_ncores(int cpu_ncores) {
    impl_->data.cpu_ncores = cpu_ncores;
  }

  int SystemConfiguration::get_cpu_ncores() const {
    return impl_->data.cpu_ncores;
  }

  void SystemConfiguration::set_cpu_nthread(int cpu_nthread) {
    impl_->data.cpu_nthread = cpu_nthread;
  }

  int SystemConfiguration::get_cpu_nthread() const {
    return impl_->data.cpu_nthread;
  }

  void SystemConfiguration::set_gpu_device(int gpu_device) {
    impl_->data.gpu_device = gpu_device;
  }

  int SystemConfiguration::get_gpu_device() const {
    return impl_->data.gpu_device;
  }

  void SystemConfiguration::set_gpu_nstream(int gpu_nstream) {
    impl_->data.gpu_nstream = gpu_nstream;
  }

  int SystemConfiguration::get_gpu_nstream() const {
    return impl_->data.gpu_nstream;
  }

  void SystemConfiguration::set_nstream(int nstream) {
    impl_->data.nstream = nstream;
  }

  int SystemConfiguration::get_nstream() const {
    return impl_->data.nstream;
  }

  void SystemConfiguration::set_active(bool active) {
    impl_->data.active = active;
  }

  bool SystemConfiguration::get_active() const {
    return impl_->data.active;
  }

  bool SystemConfiguration::is_host() const {
    return impl_->data.is_host();
  }

  bool SystemConfiguration::is_device() const {
    return impl_->data.is_device();
  }

  bool SystemConfiguration::is_float() const {
    return impl_->data.is_float();
  }

  bool SystemConfiguration::is_double() const {
    return impl_->data.is_double();
  }

  bool SystemConfiguration::is_float_host() const {
    return impl_->data.is_float_host();
  }

  bool SystemConfiguration::is_double_host() const {
    return impl_->data.is_double_host();
  }

  bool SystemConfiguration::is_float_device() const {
    return impl_->data.is_float_device();
  }

  bool SystemConfiguration::is_double_device() const {
    return impl_->data.is_double_device();
  }

  /****************************************************************************
   * The Input interface
   ***************************************************************************/

  template <typename T>
  struct Input<T>::Data {
    Input_Multislice<T> data;
  };

  template <typename T>
  Input<T>::Input()
      : impl_(std::make_unique<Data>()) {}

  template <typename T>
  Input<T>::Input(const Input& other)
      : impl_(std::make_unique<Data>(*other.impl_)) {}

  template <typename T>
  Input<T>::Input(Input&& other) = default;

  template <typename T>
  Input<T>::Input(const Input<T>::Data& other)
      : impl_(std::make_unique<Data>(other)) {}

  template <typename T>
  Input<T>& Input<T>::operator=(const Input<T>& other) {
    *impl_ = *other.impl_;
    return *this;
  }

  template <typename T>
  Input<T>& Input<T>::operator=(Input<T>&&) = default;

  template <typename T>
  Input<T>::~Input<T>() = default;

  template <typename T>
  const Input<T>::Data& Input<T>::internal() const {
    return *impl_;
  }

  template <typename T>
  SystemConfiguration Input<T>::get_system_conf() const {
    return SystemConfiguration(SystemConfiguration::Data(impl_->data.system_conf));
  }

  template <typename T>
  void Input<T>::set_system_conf(const SystemConfiguration& system_conf) {
    impl_->data.system_conf = system_conf.internal().data;
  }

  template <typename T>
  eElec_Spec_Int_Model Input<T>::get_interaction_model() const {
    return impl_->data.interaction_model;
  }

  template <typename T>
  void Input<T>::set_interaction_model(eElec_Spec_Int_Model interaction_model) {
    impl_->data.interaction_model = interaction_model;
  }

  template <typename T>
  ePotential_Type Input<T>::get_potential_type() const {
    return impl_->data.potential_type;
  }

  template <typename T>
  void Input<T>::set_potential_type(ePotential_Type potential_type) {
    impl_->data.potential_type = potential_type;
  }

  template <typename T>
  ePhonon_Model Input<T>::get_pn_model() const {
    return impl_->data.pn_model;
  }

  template <typename T>
  void Input<T>::set_pn_model(ePhonon_Model pn_model) {
    impl_->data.pn_model = pn_model;
  }

  template <typename T>
  bool Input<T>::get_pn_coh_contrib() const {
    return impl_->data.pn_coh_contrib;
  }

  template <typename T>
  void Input<T>::set_pn_coh_contrib(bool pn_coh_contrib) {
    impl_->data.pn_coh_contrib = pn_coh_contrib;
  }

  template <typename T>
  bool Input<T>::get_pn_single_conf() const {
    return impl_->data.pn_single_conf;
  }

  template <typename T>
  void Input<T>::set_pn_single_conf(bool pn_single_conf) {
    impl_->data.pn_single_conf = pn_single_conf;
  }

  /* const Input<T>::FP_Dim& get_pn_dim() const { */
  /*   return impl_->data.pn_dim; */
  /* } */

  /* void Input<T>::set_pn_dim(const FP_Dim &pn_dim) { */
  /*   impl_->data.pn_dim = pn_dim; */
  /* } */

  template <typename T>
  int Input<T>::get_fp_dist() const {
    return impl_->data.fp_dist;
  }

  template <typename T>
  void Input<T>::set_fp_dist(int fp_dist) {
    impl_->data.fp_dist = fp_dist;
  }

  template <typename T>
  int Input<T>::get_pn_seed() const {
    return impl_->data.pn_seed;
  }

  template <typename T>
  void Input<T>::set_pn_seed(int pn_seed) {
    impl_->data.pn_seed = pn_seed;
  }

  template <typename T>
  int Input<T>::get_pn_nconf() const {
    return impl_->data.pn_nconf;
  }

  template <typename T>
  void Input<T>::set_pn_nconf(int pn_nconf) {
    impl_->data.pn_nconf = pn_nconf;
  }

  template <typename T>
  int Input<T>::get_fp_iconf_0() const {
    return impl_->data.fp_iconf_0;
  }

  template <typename T>
  void Input<T>::set_fp_iconf_0(int fp_iconf_0) {
    impl_->data.fp_iconf_0 = fp_iconf_0;
  }

  /* const Atom_Data<T>& get_atoms() const; */
  /* void set_atoms(const Atom_Data<T>& atoms); */

  template <typename T>
  bool Input<T>::get_is_crystal() const {
    return impl_->data.is_crystal;
  }

  template <typename T>
  void Input<T>::set_is_crystal(bool is_crystal) {
    impl_->data.is_crystal = is_crystal;
  }

  template <typename T>
  double Input<T>::get_spec_rot_theta() const {
    return impl_->data.spec_rot_theta;
  }

  template <typename T>
  void Input<T>::set_spec_rot_theta(double spec_rot_theta) {
    impl_->data.spec_rot_theta = spec_rot_theta;
  }

  /* const r3d<T>& get_spec_rot_u0() const; */
  /* void set_spec_rot_u0(const r3d<T>& spec_rot_u0); */

  template <typename T>
  eRot_Point_Type Input<T>::get_spec_rot_center_type() const {
    return impl_->data.spec_rot_center_type;
  }

  template <typename T>
  void Input<T>::set_spec_rot_center_type(eRot_Point_Type spec_rot_center_type) {
    impl_->data.spec_rot_center_type = spec_rot_center_type;
  }

  /* const r3d<T>& get_spec_rot_center_p() const; */
  /* void set_spec_rot_center_p(const r3d<T>& spec_rot_center_p); */

  template <typename T>
  eThick_Type Input<T>::get_thick_type() const {
    return impl_->data.thick_type;
  }

  template <typename T>
  void Input<T>::set_thick_type(eThick_Type thick_type) {
    impl_->data.thick_type = thick_type;
  }

  /* const host_vector<T>& get_thick() const; */
  /* void set_thick(const host_vector<T>& thick); */

  template <typename T>
  ePotential_Slicing Input<T>::get_potential_slicing() const {
    return impl_->data.potential_slicing;
  }

  template <typename T>
  void Input<T>::set_potential_slicing(ePotential_Slicing potential_slicing) {
    impl_->data.potential_slicing = potential_slicing;
  }

  /* const Grid_2d<T>& get_grid_2d() const; */
  /* void set_grid_2d(const Grid_2d<T>& grid_2d); */

  /* const Range_2d& get_output_area() const; */
  /* void set_output_area(const Range_2d& output_area); */

  template <typename T>
  eTEM_Sim_Type Input<T>::get_simulation_type() const {
    return impl_->data.simulation_type;
  }

  template <typename T>
  void Input<T>::set_simulation_type(eTEM_Sim_Type simulation_type) {
    impl_->data.simulation_type = simulation_type;
  }

  template <typename T>
  eIncident_Wave_Type Input<T>::get_iw_type() const {
    return impl_->data.iw_type;
  }

  template <typename T>
  void Input<T>::set_iw_type(eIncident_Wave_Type iw_type) {
    impl_->data.iw_type = iw_type;
  }

  /* const host_vector<complex<T>>& get_iw_psi() const; */
  /* void set_iw_psi(const host_vector<complex<T>>& iw_psi); */

  /* const host_vector<T>& get_iw_x() const; */
  /* void set_iw_x(const host_vector<T>& iw_x); */

  /* const host_vector<T>& get_iw_y() const; */
  /* void set_iw_y(const host_vector<T>& iw_y); */

  template <typename T>
  double Input<T>::get_E_0() const {
    return impl_->data.E_0;
  }

  template <typename T>
  void Input<T>::set_E_0(double E_0) {
    impl_->data.E_0 = E_0;
  }

  template <typename T>
  double Input<T>::get_lambda() const {
    return impl_->data.lambda;
  }

  template <typename T>
  void Input<T>::set_lambda(double lambda) {
    impl_->data.lambda = lambda;
  }

  template <typename T>
  double Input<T>::get_theta() const {
    return impl_->data.theta;
  }

  template <typename T>
  void Input<T>::set_theta(double theta) {
    impl_->data.theta = theta;
  }

  template <typename T>
  double Input<T>::get_phi() const {
    return impl_->data.phi;
  }

  template <typename T>
  void Input<T>::set_phi(double phi) {
    impl_->data.phi = phi;
  }

  template <typename T>
  eIllumination_Model Input<T>::get_illumination_model() const {
    return impl_->data.illumination_model;
  }

  template <typename T>
  void Input<T>::set_illumination_model(eIllumination_Model illumination_model) {
    impl_->data.illumination_model = illumination_model;
  }

  template <typename T>
  eTemporal_Spatial_Incoh Input<T>::get_temporal_spatial_incoh() const {
    return impl_->data.temporal_spatial_incoh;
  }

  template <typename T>
  void Input<T>::set_temporal_spatial_incoh(eTemporal_Spatial_Incoh temporal_spatial_incoh) {
    impl_->data.temporal_spatial_incoh = temporal_spatial_incoh;
  }

  /* const Lens<T>& get_cond_lens() const; */
  /* void set_cond_lens(const Lens<T>& cond_lens); */

  /* const Lens<T>& get_obj_lens() const; */
  /* void set_obj_lens(const Lens<T>& obj_lens); */

  /* const Scanning<T>& get_scanning() const; */
  /* void set_scanning(const Scanning<T>& scanning); */

  /* const Detector<T, e_host>& get_detector() const; */
  /* void set_detector(const Detector<T, e_host>& detector); */

  /* const EELS<T>& get_eels_fr() const; */
  /* void set_eels_fr(const EELS<T>& eels_fr); */

  template <typename T>
  eOperation_Mode Input<T>::get_operation_mode() const {
    return impl_->data.operation_mode;
  }

  template <typename T>
  void Input<T>::set_operation_mode(eOperation_Mode operation_mode) {
    impl_->data.operation_mode = operation_mode;
  }

  template <typename T>
  bool Input<T>::get_slice_storage() const {
    return impl_->data.slice_storage;
  }

  template <typename T>
  void Input<T>::set_slice_storage(bool slice_storage) {
    impl_->data.slice_storage = slice_storage;
  }

  template <typename T>
  bool Input<T>::get_reverse_multislice() const {
    return impl_->data.reverse_multislice;
  }

  template <typename T>
  void Input<T>::set_reverse_multislice(bool reverse_multislice) {
    impl_->data.reverse_multislice = reverse_multislice;
  }

  template <typename T>
  int Input<T>::get_mul_sign() const {
    return impl_->data.mul_sign;
  }

  template <typename T>
  void Input<T>::set_mul_sign(int mul_sign) {
    impl_->data.mul_sign = mul_sign;
  }

  template <typename T>
  double Input<T>::get_Vrl() const {
    return impl_->data.Vrl;
  }

  template <typename T>
  void Input<T>::set_Vrl(double Vrl) {
    impl_->data.Vrl = Vrl;
  }

  template <typename T>
  int Input<T>::get_nR() const {
    return impl_->data.nR;
  }

  template <typename T>
  void Input<T>::set_nR(int nR) {
    impl_->data.nR = nR;
  }

  template <typename T>
  int Input<T>::get_nrot() const {
    return impl_->data.nrot;
  }

  template <typename T>
  void Input<T>::set_nrot(int nrot) {
    impl_->data.nrot = nrot;
  }

  template <typename T>
  eLens_Var_Type Input<T>::get_cdl_var_type() const {
    return impl_->data.cdl_var_type;
  }

  template <typename T>
  void Input<T>::set_cdl_var_type(eLens_Var_Type cdl_var_type) {
    impl_->data.cdl_var_type = cdl_var_type;
  }

  /* const host_vector<T>& get_cdl_var() const; */
  /* void set_cdl_var(const host_vector<T>& cdl_var); */

  /* const host_vector<int>& get_iscan() const; */
  /* void set_iscan(const host_vector<int>& iscan); */

  /* const host_vector<T>& get_beam_x() const; */
  /* void set_beam_x(const host_vector<T>& beam_x); */

  /* const host_vector<T>& get_beam_y() const; */
  /* void set_beam_y(const host_vector<T>& beam_y); */

  template <typename T>
  int Input<T>::get_islice() const {
    return impl_->data.islice;
  }

  template <typename T>
  void Input<T>::set_islice(int islice) {
    impl_->data.islice = islice;
  }

  template <typename T>
  bool Input<T>::get_dp_Shift() const {
    return impl_->data.dp_Shift;
  }

  template <typename T>
  void Input<T>::set_dp_Shift(bool dp_Shift) {
    impl_->data.dp_Shift = dp_Shift;
  }

  /****************************************************************************
   * Misc function calls
   ***************************************************************************/

  template <typename T>
  void test(const Input<T>& a) {
    std::cout << a.get_system_conf().get_precision() << std::endl;
  }

  /**
   * Explict instantiation of template functions
   */
  template void test<float>(const Input<float>&);
  template void test<double>(const Input<double>&);
}  // namespace mt
