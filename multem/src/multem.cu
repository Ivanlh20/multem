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
#include <output_multislice.hpp>
#include <tem_simulation.cuh>

namespace mt {
  
  /****************************************************************************
   * The FFT interface
   ***************************************************************************/
  
  template <typename T, eDevice dev>
  struct FFTData<T, dev>::Data {
    FFT<T, dev> data;
  };

  template <typename T, eDevice dev>
  FFTData<T, dev>::FFTData()
      : impl_(std::make_unique<Data>()) {}

  template <typename T, eDevice dev>
  FFTData<T, dev>::FFTData(const FFTData& other)
      : impl_(std::make_unique<Data>(*other.impl_)) {}

  template <typename T, eDevice dev>
  FFTData<T, dev>::FFTData(FFTData&& other) = default;

  template <typename T, eDevice dev>
  FFTData<T, dev>& FFTData<T, dev>::operator=(const FFTData& other) {
    *impl_ = *other.impl_;
    return *this;
  }

  template <typename T, eDevice dev>
  FFTData<T, dev>::FFTData(const Data& other)
      : impl_(std::make_unique<Data>(other)) {}

  template <typename T, eDevice dev>
  FFTData<T, dev>& FFTData<T, dev>::operator=(FFTData&&) = default;

  template <typename T, eDevice dev>
  FFTData<T, dev>::~FFTData() = default;

  template <typename T, eDevice dev>
  const FFTData<T, dev>::Data& FFTData<T, dev>::internal() const {
    return *impl_;
  }
  
  template <typename T, eDevice dev>
  FFTData<T, dev>::Data& FFTData<T, dev>::internal() {
    return *impl_;
  }

  template <typename T, eDevice dev>
  void FFTData<T, dev>::cleanup() {
    impl_->data.cleanup();
  }
  
  template <typename T, eDevice dev>
  void FFTData<T, dev>::destroy_plan() {
    impl_->data.destroy_plan();
  }
  
  template <typename T, eDevice dev>
  void FFTData<T, dev>::create_plan_1d(const int &nx, int nThread) {
    impl_->data.create_plan_1d(nx, nThread);
  }
  
  template <typename T, eDevice dev>
  void FFTData<T, dev>::create_plan_1d_batch(const int &ny, const int &nx, int nThread) {
    impl_->data.create_plan_1d_batch(ny, nx, nThread);
  }
  
  template <typename T, eDevice dev>
  void FFTData<T, dev>::create_plan_2d(const int &ny, const int &nx, int nThread) {
    impl_->data.create_plan_2d(ny, nx, nThread);
  }
  
  template <typename T, eDevice dev>
  void FFTData<T, dev>::create_plan_2d_batch(const int &ny, const int &nx, const int &nz, int nThread) {
    impl_->data.create_plan_2d_batch(ny, nx, nz, nThread);
  }
  
  /****************************************************************************
   * The Stream interface
   ***************************************************************************/
  
  template <eDevice dev>
  struct StreamIface<dev>::Data {
    Stream<dev> data;
    Data() = default;
    Data(int new_stream)
      : data(new_stream) {}
  };

  template <eDevice dev>
  StreamIface<dev>::StreamIface()
      : impl_(std::make_unique<Data>()) {}

  template <eDevice dev>
  StreamIface<dev>::~StreamIface() = default;

  template <eDevice dev>
  const StreamIface<dev>::Data& StreamIface<dev>::internal() const {
    return *impl_;
  }
  
  template <eDevice dev>
  StreamIface<dev>::Data& StreamIface<dev>::internal() {
    return *impl_;
  }

  template <eDevice dev>
  StreamIface<dev>::StreamIface(int new_stream)
    : impl_(std::make_unique<Data>(new_stream)) {}

  template <eDevice dev>
  int StreamIface<dev>::size() const {
    return impl_->data.size();
  }

  template <eDevice dev>
  void StreamIface<dev>::resize(int new_nstream) {
    impl_->data.resize(new_nstream);
  }

  template <eDevice dev>
  void StreamIface<dev>::synchronize() {
    impl_->data.synchronize();
  }

  template <eDevice dev>
  void StreamIface<dev>::set_n_act_stream(const int &new_n_act_stream) {
    impl_->data.set_n_act_stream(new_n_act_stream);
  }

  template <eDevice dev>
  void StreamIface<dev>::set_grid(const int &nx_i, const int &ny_i) {
    impl_->data.set_grid(nx_i, ny_i);
  }

  template <eDevice dev>
  Range_2d StreamIface<dev>::get_range(const int &istream) {
    return impl_->data.get_range(istream);
  }

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
   * The AtomData interface
   ***************************************************************************/

  template <typename T>
  struct AtomData<T>::Data {
    Atom_Data<T> data;
    Data() {}
    Data(const Atom_Data<T>& d)
        : data(d) {}
  };

  template <typename T>
  AtomData<T>::AtomData()
      : impl_(std::make_unique<Data>()) {}

  template <typename T>
  AtomData<T>::AtomData(const AtomData& other)
      : impl_(std::make_unique<Data>(*other.impl_)) {}

  template <typename T>
  AtomData<T>::AtomData(AtomData&& other) = default;

  template <typename T>
  AtomData<T>::AtomData(const AtomData<T>::Data& other)
      : impl_(std::make_unique<Data>(other)) {}

  template <typename T>
  AtomData<T>& AtomData<T>::operator=(const AtomData<T>& other) {
    *impl_ = *other.impl_;
    return *this;
  }

  template <typename T>
  AtomData<T>& AtomData<T>::operator=(AtomData<T>&&) = default;

  template <typename T>
  AtomData<T>::~AtomData<T>() = default;

  template <typename T>
  const AtomData<T>::Data& AtomData<T>::internal() const {
    return *impl_;
  }
  
  /****************************************************************************
   * The Scanning interface
   ***************************************************************************/

  template <typename T>
  struct ScanningData<T>::Data {
    Scanning<T> data;
    Data() {}
    Data(const Scanning<T>& d)
        : data(d) {}
  };

  template <typename T>
  ScanningData<T>::ScanningData()
      : impl_(std::make_unique<Data>()) {}

  template <typename T>
  ScanningData<T>::ScanningData(const ScanningData& other)
      : impl_(std::make_unique<Data>(*other.impl_)) {}

  template <typename T>
  ScanningData<T>::ScanningData(ScanningData&& other) = default;

  template <typename T>
  ScanningData<T>::ScanningData(const ScanningData<T>::Data& other)
      : impl_(std::make_unique<Data>(other)) {}

  template <typename T>
  ScanningData<T>& ScanningData<T>::operator=(const ScanningData<T>& other) {
    *impl_ = *other.impl_;
    return *this;
  }

  template <typename T>
  ScanningData<T>& ScanningData<T>::operator=(ScanningData<T>&&) = default;

  template <typename T>
  ScanningData<T>::~ScanningData<T>() = default;

  template <typename T>
  const ScanningData<T>::Data& ScanningData<T>::internal() const {
    return *impl_;
  }
  
  /****************************************************************************
   * The Detector interface
   ***************************************************************************/

  template <typename T>
  struct DetectorData<T>::Data {
    Detector<T, e_host> data;
    Data() {}
    Data(const Detector<T, e_host>& d)
        : data(d) {}
  };

  template <typename T>
  DetectorData<T>::DetectorData()
      : impl_(std::make_unique<Data>()) {}

  template <typename T>
  DetectorData<T>::DetectorData(const DetectorData& other)
      : impl_(std::make_unique<Data>(*other.impl_)) {}

  template <typename T>
  DetectorData<T>::DetectorData(DetectorData&& other) = default;

  template <typename T>
  DetectorData<T>::DetectorData(const DetectorData<T>::Data& other)
      : impl_(std::make_unique<Data>(other)) {}

  template <typename T>
  DetectorData<T>& DetectorData<T>::operator=(const DetectorData<T>& other) {
    *impl_ = *other.impl_;
    return *this;
  }

  template <typename T>
  DetectorData<T>& DetectorData<T>::operator=(DetectorData<T>&&) = default;

  template <typename T>
  DetectorData<T>::~DetectorData<T>() = default;

  template <typename T>
  const DetectorData<T>::Data& DetectorData<T>::internal() const {
    return *impl_;
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
  Input<T>::Data& Input<T>::internal() {
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

  template <typename T>
  const FP_Dim& Input<T>::get_pn_dim() const {
    return impl_->data.pn_dim;
  }

  template <typename T>
  void Input<T>::set_pn_dim(const FP_Dim &pn_dim) {
    impl_->data.pn_dim = pn_dim;
  }

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

  template <typename T>
  AtomData<T> Input<T>::get_atoms() const {
    return AtomData<T>(typename AtomData<T>::Data(impl_->data.atoms));
  }

  template <typename T>
  void Input<T>::set_atoms(const AtomData<T>& atoms) {
    impl_->data.atoms = atoms.internal().data;
  }

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

  template <typename T>
  const r3d<T>& Input<T>::get_spec_rot_u0() const {
    return impl_->data.spec_rot_u0;
  }

  template <typename T>
  void Input<T>::set_spec_rot_u0(const r3d<T>& spec_rot_u0) {
    impl_->data.spec_rot_u0 = spec_rot_u0;
  }

  template <typename T>
  eRot_Point_Type Input<T>::get_spec_rot_center_type() const {
    return impl_->data.spec_rot_center_type;
  }

  template <typename T>
  void Input<T>::set_spec_rot_center_type(eRot_Point_Type spec_rot_center_type) {
    impl_->data.spec_rot_center_type = spec_rot_center_type;
  }

  template <typename T>
  const r3d<T>& Input<T>::get_spec_rot_center_p() const {
    return impl_->data.spec_rot_center_p;
  }

  template <typename T>
  void Input<T>::set_spec_rot_center_p(const r3d<T>& spec_rot_center_p) {
    impl_->data.spec_rot_center_p = spec_rot_center_p;
  }

  template <typename T>
  eThick_Type Input<T>::get_thick_type() const {
    return impl_->data.thick_type;
  }

  template <typename T>
  void Input<T>::set_thick_type(eThick_Type thick_type) {
    impl_->data.thick_type = thick_type;
  }

  template <typename T>
  std::vector<T> Input<T>::get_thick() const {
    return std::vector<T>(impl_->data.thick.begin(), impl_->data.thick.end());
  }

  template <typename T>
  void Input<T>::set_thick(const std::vector<T>& thick) {
    impl_->data.thick.assign(thick.begin(), thick.end());
  }

  template <typename T>
  ePotential_Slicing Input<T>::get_potential_slicing() const {
    return impl_->data.potential_slicing;
  }

  template <typename T>
  void Input<T>::set_potential_slicing(ePotential_Slicing potential_slicing) {
    impl_->data.potential_slicing = potential_slicing;
  }

  template <typename T>
  const Grid_2d<T>& Input<T>::get_grid_2d() const {
    return impl_->data.grid_2d;
  }

  template <typename T>
  void Input<T>::set_grid_2d(const Grid_2d<T>& grid_2d) {
    impl_->data.grid_2d = grid_2d;
  }

  template <typename T>
  const Range_2d& Input<T>::get_output_area() const {
    return impl_->data.output_area;
  }

  template <typename T>
  void Input<T>::set_output_area(const Range_2d& output_area) {
    impl_->data.output_area = output_area;
  }

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

  template <typename T>
  std::vector<complex<T>> Input<T>::get_iw_psi() const {
    return std::vector<complex<T>>(impl_->data.iw_psi.begin(), impl_->data.iw_psi.end());
  }

  template <typename T>
  void Input<T>::set_iw_psi(const std::vector<complex<T>>& iw_psi) {
    impl_->data.iw_psi.assign(iw_psi.begin(), iw_psi.end());
  }

  template <typename T>
  std::vector<T> Input<T>::get_iw_x() const {
    return std::vector<T>(impl_->data.iw_x.begin(), impl_->data.iw_x.end());
  }

  template <typename T>
  void Input<T>::set_iw_x(const std::vector<T>& iw_x) {
    impl_->data.iw_x.assign(iw_x.begin(), iw_x.end());
  }

  template <typename T>
  std::vector<T> Input<T>::get_iw_y() const {
    return std::vector<T>(impl_->data.iw_y.begin(), impl_->data.iw_y.end());
  }

  template <typename T>
  void Input<T>::set_iw_y(const std::vector<T>& iw_y) {
    impl_->data.iw_y.assign(iw_y.begin(), iw_y.end());
  }

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

  template <typename T>
  const Lens<T>& Input<T>::get_cond_lens() const {
    return impl_->data.cond_lens;
  }
  
  template <typename T>
  void Input<T>::set_cond_lens(const Lens<T>& cond_lens) {
    impl_->data.cond_lens = cond_lens;
  }

  template <typename T>
  const Lens<T>& Input<T>::get_obj_lens() const {
    return impl_->data.obj_lens;
  }

  template <typename T>
  void Input<T>::set_obj_lens(const Lens<T>& obj_lens) {
    impl_->data.obj_lens = obj_lens;
  }

  template <typename T>
  ScanningData<T> Input<T>::get_scanning() const {
    return ScanningData<T>(typename ScanningData<T>::Data(impl_->data.scanning));
  }
  
  template <typename T>
  void Input<T>::set_scanning(const ScanningData<T>& scanning) {
    impl_->data.scanning = scanning.internal().data;
  }

  template <typename T>
  DetectorData<T> Input<T>::get_detector() const {
    return DetectorData<T>(typename DetectorData<T>::Data(impl_->data.detector));
  }
  
  template <typename T>
  void Input<T>::set_detector(const DetectorData<T>& detector) {
    impl_->data.detector = detector.internal().data;
  }

  template <typename T>
  const EELS<T>& Input<T>::get_eels_fr() const {
    return impl_->data.eels_fr;
  }
  
  template <typename T>
  void Input<T>::set_eels_fr(const EELS<T>& eels_fr) {
    impl_->data.eels_fr = eels_fr;
  }

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

  template <typename T>
  std::vector<T> Input<T>::get_cdl_var() const {
    return std::vector<T>(impl_->data.cdl_var.begin(), impl_->data.cdl_var.end());
  }

  template <typename T>
  void Input<T>::set_cdl_var(const std::vector<T>& cdl_var) {
    impl_->data.cdl_var.assign(cdl_var.begin(), cdl_var.end());
  }

  template <typename T>
  std::vector<int> Input<T>::get_iscan() const {
    return std::vector<int>(impl_->data.iscan.begin(), impl_->data.iscan.end());
  }

  template <typename T>
  void Input<T>::set_iscan(const std::vector<int>& iscan) {
    impl_->data.iscan.assign(iscan.begin(), iscan.end());
  }

  template <typename T>
  std::vector<T> Input<T>::get_beam_x() const {
    return std::vector<T>(impl_->data.beam_x.begin(), impl_->data.beam_x.end());
  }

  template <typename T>
  void Input<T>::set_beam_x(const std::vector<T>& beam_x) {
    impl_->data.beam_x.assign(beam_x.begin(), beam_x.end());
  }

  template <typename T>
  std::vector<T> Input<T>::get_beam_y() const {
    return std::vector<T>(impl_->data.beam_y.begin(), impl_->data.beam_y.end());
  }

  template <typename T>
  void Input<T>::set_beam_y(const std::vector<T>& beam_y) {
    impl_->data.beam_y.assign(beam_y.begin(), beam_y.end());
  }

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

  template <typename T>
  void Input<T>::validate_parameters() {
    impl_->data.validate_parameters();
  }

  /****************************************************************************
   * The Output interface
   ***************************************************************************/

  template <typename T>
  struct Output<T>::Data {
    Output_Multislice<T> data;

    // Constructors required due to behaviour of Output_Multislice. Would
    // probably be better if we could fix Output_Multislice so that the normal
    // copy constructors work.
    Data() = default;
    Data(Data&) = default;
    Data& operator=(Data& other) {
      data = other.data;
      return *this;
    }
  };

  template <typename T>
  Output<T>::Output()
      : impl_(std::make_unique<Data>()) {}

  template <typename T>
  Output<T>::Output(Output& other)
      : impl_(std::make_unique<Data>(*other.impl_)) {}

  template <typename T>
  Output<T>::Output(Output&& other) = default;

  template <typename T>
  Output<T>::Output(Output<T>::Data& other)
      : impl_(std::make_unique<Data>(other)) {}

  template <typename T>
  Output<T>& Output<T>::operator=(Output<T>& other) {
    *impl_ = *other.impl_;
    return *this;
  }

  template <typename T>
  Output<T>& Output<T>::operator=(Output<T>&&) = default;

  template <typename T>
  Output<T>::~Output<T>() = default;

  template <typename T>
  const Output<T>::Data& Output<T>::internal() const {
    return *impl_;
  }
  
  template <typename T>
  Output<T>::Data& Output<T>::internal() {
    return *impl_;
  }

  template <typename T>
  eTEM_Output_Type Output<T>::get_output_type() const {
    return impl_->data.output_type;
  }

  template <typename T>
  int Output<T>::get_ndetector() const {
    return impl_->data.ndetector;
  }

  template <typename T>
  int Output<T>::get_nx() const {
    return impl_->data.nx;
  }

  template <typename T>
  int Output<T>::get_ny() const {
    return impl_->data.ny;
  }

  template <typename T>
  T Output<T>::get_dx() const {
    return impl_->data.dx;
  }

  template <typename T>
  T Output<T>::get_dy() const {
    return impl_->data.dy;
  }

  template <typename T>
  T Output<T>::get_dr() const {
    return impl_->data.dr;
  }

  template <typename T>
  std::vector<T> Output<T>::get_x() const {
    return std::vector<T>(impl_->data.x.begin(), impl_->data.x.end());
  }

  template <typename T>
  std::vector<T> Output<T>::get_y() const {
    return std::vector<T>(impl_->data.y.begin(), impl_->data.y.end());
  }

  template <typename T>
  std::vector<T> Output<T>::get_r() const {
    return std::vector<T>(impl_->data.r.begin(), impl_->data.r.end());
  }
    
  template <typename T>
  std::vector<DetInt<typename Output<T>::vector_type>> Output<T>::get_image_tot() const {
    std::vector<DetInt<vector_type>> result;
    for (auto &x : impl_->data.image_tot) {
      DetInt<vector_type> d;
      for (auto &y : x.image) {
        d.image.push_back(vector_type(y.begin(), y.end()));
      }
      result.push_back(d);
    }
    return result;
  }

  template <typename T>
  std::vector<DetInt<typename Output<T>::vector_type>> Output<T>::get_image_coh() const {
    std::vector<DetInt<vector_type>> result;
    for (auto &x : impl_->data.image_coh) {
      DetInt<vector_type> d;
      for (auto &y : x.image) {
        d.image.push_back(vector_type(y.begin(), y.end()));
      }
      result.push_back(d);
    }
    return result;
  }

  template <typename T>
  std::vector<typename Output<T>::vector_type> Output<T>::get_m2psi_tot() const {
    std::vector<vector_type> result;
    for (auto &x : impl_->data.m2psi_tot) {
      result.push_back(vector_type(x.begin(), x.end()));
    }
    return result;
  }

  template <typename T>
  std::vector<typename Output<T>::vector_type> Output<T>::get_m2psi_coh() const {
    std::vector<vector_type> result;
    for (auto &x : impl_->data.m2psi_coh) {
      result.push_back(vector_type(x.begin(), x.end()));
    }
    return result;
  }

  template <typename T>
  std::vector<typename Output<T>::complex_vector_type> Output<T>::get_psi_coh() const {
    std::vector<complex_vector_type> result;
    for (auto &x : impl_->data.psi_coh) {
      result.push_back(complex_vector_type(x.begin(), x.end()));
    }
    return result;
  }

  template <typename T>
  std::vector<typename Output<T>::vector_type> Output<T>::get_V() const {
    std::vector<vector_type> result;
    for (auto &x : impl_->data.V) {
      result.push_back(vector_type(x.begin(), x.end()));
    }
    return result;
  }

  template <typename T>
  std::vector<typename Output<T>::complex_vector_type> Output<T>::get_trans() const {
    std::vector<complex_vector_type> result;
    for (auto &x : impl_->data.trans) {
      result.push_back(complex_vector_type(x.begin(), x.end()));
    }
    return result;
  }

  template <typename T>
  std::vector<typename Output<T>::complex_vector_type> Output<T>::get_psi_0() const {
    std::vector<complex_vector_type> result;
    for (auto &x : impl_->data.psi_0) {
      result.push_back(complex_vector_type(x.begin(), x.end()));
    }
    return result;
  }

  template <typename T>
  std::vector<bool> Output<T>::get_thk_gpu() const {
    return std::vector<bool>(impl_->data.thk_gpu.begin(), impl_->data.thk_gpu.end());
  }

  template <typename T>
  void Output<T>::gather() {
    impl_->data.gather();
  }
  
  template <typename T>
  void Output<T>::clean_temporal() {
    impl_->data.clean_temporal();
  }
  
  /****************************************************************************
   * The Multislice interface
   ***************************************************************************/

  template <typename T, eDevice dev>
  struct MultisliceData<T, dev>::Data {
    Multislice<T, dev> data;
  };

  template <typename T, eDevice dev>
  MultisliceData<T, dev>::MultisliceData()
      : impl_(std::make_unique<Data>()) {}

  template <typename T, eDevice dev>
  const MultisliceData<T, dev>::Data& MultisliceData<T, dev>::internal() const {
    return *impl_;
  }
		
  template <typename T, eDevice dev>
  void MultisliceData<T, dev>::set_input_data(Input<T> &input, StreamIface<dev> &stream_i, FFTData<T, dev> &fft2_i) {
    impl_->data.set_input_data(
        &input.internal().data,
        &stream_i.internal().data,
        &fft2_i.internal().data);
  }
		
  template <typename T, eDevice dev>
  void MultisliceData<T, dev>::operator()(Output<T> &output) {
    impl_->data(output.internal().data);
  }

  /****************************************************************************
   * Misc function calls
   ***************************************************************************/

  template <typename T>
  void test(const Input<T>& a) {
    std::cout << a.get_system_conf().get_precision() << std::endl;
  }
  
  /****************************************************************************
   * Explict instantiation of template classes and functions
   ***************************************************************************/
  
  /**
   * Explict instantiation of template classes
   */
  template class FFTData<float, e_host>;
  template class FFTData<float, e_device>;
  template class FFTData<double, e_host>;
  template class FFTData<double, e_device>;

  template class StreamIface<e_host>;
  template class StreamIface<e_device>;

  template class Input<float>;
  template class Input<double>;

  template class Output<float>;
  template class Output<double>;
  
  template class MultisliceData<float, e_host>;
  template class MultisliceData<float, e_device>;
  template class MultisliceData<double, e_host>;
  template class MultisliceData<double, e_device>;

  /**
   * Explict instantiation of template functions
   */
  template void test<float>(const Input<float>&);
  template void test<double>(const Input<double>&);

}  // namespace mt
