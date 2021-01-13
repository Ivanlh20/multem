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
#include <multem/error.h>
#include <types.cuh>
#include <input_multislice.cuh>
#include <output_multislice.hpp>
#include <tem_simulation.cuh>

namespace mt {
  

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
  void Output<T>::set_input_data(Input_Multislice<T> &input) {
    impl_->data.set_input_data(&input);
  }

  template <typename T>
  eTEM_Output_Type Output<T>::get_output_type() const {
    return impl_->data.output_type;
  }
  
  template <typename T>
  void Output<T>::set_output_type(eTEM_Output_Type output_type) {
    impl_->data.output_type = output_type;
  }

  template <typename T>
  int Output<T>::get_ndetector() const {
    return impl_->data.ndetector;
  }
		
  template <typename T>
  void Output<T>::set_ndetector(int ndetector) {
    impl_->data.ndetector = ndetector;
  }

  template <typename T>
  int Output<T>::get_nx() const {
    return impl_->data.nx;
  }

  template <typename T>
  void Output<T>::set_nx(int nx) {
    impl_->data.nx = nx;
  }

  template <typename T>
  int Output<T>::get_ny() const {
    return impl_->data.ny;
  }

  template <typename T>
  void Output<T>::set_ny(int ny) {
    impl_->data.ny = ny;
  }

  template <typename T>
  T Output<T>::get_dx() const {
    return impl_->data.dx;
  }

  template <typename T>
  void Output<T>::set_dx(T dx) {
    impl_->data.dx = dx;
  }

  template <typename T>
  T Output<T>::get_dy() const {
    return impl_->data.dy;
  }

  template <typename T>
  void Output<T>::set_dy(T dy) {
    impl_->data.dy = dy;
  }

  template <typename T>
  T Output<T>::get_dr() const {
    return impl_->data.dr;
  }

  template <typename T>
  void Output<T>::set_dr(T dr) {
    impl_->data.dr = dr;
  }

  template <typename T>
  std::vector<T> Output<T>::get_x() const {
    return std::vector<T>(impl_->data.x.begin(), impl_->data.x.end());
  }

  template <typename T>
  void Output<T>::set_x(const std::vector<T>& x) {
    impl_->data.x.assign(x.begin(), x.end());
  }

  template <typename T>
  std::vector<T> Output<T>::get_y() const {
    return std::vector<T>(impl_->data.y.begin(), impl_->data.y.end());
  }

  template <typename T>
  void Output<T>::set_y(const std::vector<T>& y) {
    impl_->data.y.assign(y.begin(), y.end());
  }

  template <typename T>
  std::vector<T> Output<T>::get_r() const {
    return std::vector<T>(impl_->data.r.begin(), impl_->data.r.end());
  }

  template <typename T>
  void Output<T>::set_r(const std::vector<T>& r) {
    impl_->data.r.assign(r.begin(), r.end());
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
  void Output<T>::set_image_tot(const std::vector<DetInt<vector_type>>& image_tot) {
    typedef typename Output_Multislice<T>::TVector_hr TVector_hr;
    host_vector<Det_Int<TVector_hr>> result;
    for (auto &x : image_tot) {
      Det_Int<TVector_hr> d;
      for (auto &y : x.image) {
        d.image.push_back(TVector_hr(y.begin(), y.end()));
      }
      result.push_back(d);
    }
    impl_->data.image_tot = result;
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
  void Output<T>::set_image_coh(const std::vector<DetInt<vector_type>>& image_coh) {
    typedef typename Output_Multislice<T>::TVector_hr vector_type;
    host_vector<Det_Int<vector_type>> result;
    for (auto &x : image_coh) {
      Det_Int<vector_type> d;
      for (auto &y : x.image) {
        d.image.push_back(vector_type(y.begin(), y.end()));
      }
      result.push_back(d);
    }
    impl_->data.image_coh = result;
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
  void Output<T>::set_m2psi_tot(const std::vector<vector_type>& m2psi_tot) {
    typedef typename Output_Multislice<T>::TVector_hr TVector_hr;
    host_vector<TVector_hr> result;
    for (auto &x : m2psi_tot) {
      result.push_back(TVector_hr(x.begin(), x.end()));
    }
    impl_->data.m2psi_tot = result;
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
  void Output<T>::set_m2psi_coh(const std::vector<vector_type>& m2psi_coh) {
    typedef typename Output_Multislice<T>::TVector_hr TVector_hr;
    host_vector<TVector_hr> result;
    for (auto &x : m2psi_coh) {
      result.push_back(TVector_hr(x.begin(), x.end()));
    }
    impl_->data.m2psi_coh = result;
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
  void Output<T>::set_psi_coh(const std::vector<complex_vector_type>& psi_coh) {
    typedef typename Output_Multislice<T>::TVector_hc TVector_hc;
    host_vector<TVector_hc> result;
    for (auto &x : psi_coh) {
      result.push_back(TVector_hc(x.begin(), x.end()));
    }
    impl_->data.psi_coh = result;
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
  void Output<T>::set_V(const std::vector<vector_type>& V) {
    typedef typename Output_Multislice<T>::TVector_hr TVector_hr;
    host_vector<TVector_hr> result;
    for (auto &x : V) {
      result.push_back(TVector_hr(x.begin(), x.end()));
    }
    impl_->data.V = result;
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
  void Output<T>::set_trans(const std::vector<complex_vector_type>& trans) {
    typedef typename Output_Multislice<T>::TVector_hc TVector_hc;
    host_vector<TVector_hc> result;
    for (auto &x : trans) {
      result.push_back(TVector_hc(x.begin(), x.end()));
    }
    impl_->data.trans = result;
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
  void Output<T>::set_psi_0(const std::vector<complex_vector_type>& psi_0) {
    typedef typename Output_Multislice<T>::TVector_hc TVector_hc;
    host_vector<TVector_hc> result;
    for (auto &x : psi_0) {
      result.push_back(TVector_hc(x.begin(), x.end()));
    }
    impl_->data.psi_0 = result;
  }

  template <typename T>
  std::vector<bool> Output<T>::get_thk_gpu() const {
    return std::vector<bool>(impl_->data.thk_gpu.begin(), impl_->data.thk_gpu.end());
  }

  template <typename T>
  void Output<T>::set_thk_gpu(const std::vector<bool>& thk_gpu) {
    impl_->data.thk_gpu.assign(thk_gpu.begin(), thk_gpu.end());
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
   * Logging functions
   ***************************************************************************/

  std::string to_string(const System_Configuration &self, const std::string &prefix) {
    std::ostringstream msg;
    msg << prefix << "precision:   " << self.precision << "\n";
    msg << prefix << "device:      " << self.device << "\n";
    msg << prefix << "cpu_ncores:  " << self.precision << "\n";
    msg << prefix << "cpu_nthread: " << self.precision << "\n";
    msg << prefix << "gpu_device:  " << self.precision << "\n";
    msg << prefix << "gpu_nstream: " << self.precision << "\n";
    return msg.str();
  }
  
  template <typename T>
  std::string to_string(const Lens<T> &self, const std::string &prefix) {
    std::ostringstream msg;
    msg << prefix << "m: " << self.m << "\n";
    msg << prefix << "c_10: " << self.c_10 << "\n";
    msg << prefix << "c_12: " << self.c_12 << "\n";
    msg << prefix << "phi_12: " << self.phi_12 << "\n";
    msg << prefix << "c_21: " << self.c_21 << "\n";
    msg << prefix << "phi_21: " << self.phi_21 << "\n";
    msg << prefix << "c_23: " << self.c_23 << "\n";
    msg << prefix << "phi_23: " << self.phi_23 << "\n";
    msg << prefix << "c_30: " << self.c_30 << "\n";
    msg << prefix << "c_32: " << self.c_32 << "\n";
    msg << prefix << "phi_32: " << self.phi_32 << "\n";
    msg << prefix << "c_34: " << self.c_34 << "\n";
    msg << prefix << "phi_34: " << self.phi_34 << "\n";
    msg << prefix << "c_41: " << self.c_41 << "\n";
    msg << prefix << "phi_41: " << self.phi_41 << "\n";
    msg << prefix << "c_43: " << self.c_43 << "\n";
    msg << prefix << "phi_43: " << self.phi_43 << "\n";
    msg << prefix << "c_45: " << self.c_45 << "\n";
    msg << prefix << "phi_45: " << self.phi_45 << "\n";
    msg << prefix << "c_50: " << self.c_50 << "\n";
    msg << prefix << "c_52: " << self.c_52 << "\n";
    msg << prefix << "phi_52: " << self.phi_52 << "\n";
    msg << prefix << "c_54: " << self.c_54 << "\n";
    msg << prefix << "phi_54: " << self.phi_54 << "\n";
    msg << prefix << "c_56: " << self.c_56 << "\n";
    msg << prefix << "phi_56: " << self.phi_56 << "\n";
    msg << prefix << "inner_aper_ang: " << self.inner_aper_ang << "\n";
    msg << prefix << "outer_aper_ang: " << self.outer_aper_ang << "\n";
    msg << prefix << "ti_a: " << self.ti_a << "\n";
    msg << prefix << "ti_sigma: " << self.ti_sigma << "\n";
    msg << prefix << "ti_beta: " << self.ti_beta << "\n";
    msg << prefix << "ti_npts: " << self.ti_npts << "\n";
    msg << prefix << "ti_iehwgd: " << self.ti_iehwgd << "\n";
    msg << prefix << "si_a: " << self.si_a << "\n";
    msg << prefix << "si_sigma: " << self.si_sigma << "\n";
    msg << prefix << "si_beta: " << self.si_beta << "\n";
    msg << prefix << "si_rad_npts: " << self.si_rad_npts << "\n";
    msg << prefix << "si_azm_npts: " << self.si_azm_npts << "\n";
    msg << prefix << "si_iehwgd: " << self.si_iehwgd << "\n";
    msg << prefix << "si_theta_c: " << self.si_theta_c << "\n";
    msg << prefix << "zero_defocus_type: " << self.zero_defocus_type << "\n";
    msg << prefix << "zero_defocus_plane: " << self.zero_defocus_plane << "\n";
    msg << prefix << "lambda: " << self.lambda << "\n";
    return msg.str();
  }
  
  template <typename T>
  std::string to_string(const Input_Multislice<T> &self, const std::string &prefix) {
    std::ostringstream msg;
    msg << prefix << "system_conf" << "\n";
    msg << prefix << to_string(self.system_conf, " -");
    msg << prefix << "interaction_model: " << self.interaction_model << "\n";
    msg << prefix << "potential_type: " << self.potential_type << "\n";
    msg << prefix << "pn_model: " << self.pn_model << "\n";
    msg << prefix << "pn_coh_contrib: " << self.pn_coh_contrib << "\n";
    msg << prefix << "pn_single_conf: " << self.pn_single_conf << "\n";
    /* msg << prefix << "pn_dim: " << self.pn_dim << "\n"; */
    msg << prefix << "fp_dist: " << self.fp_dist << "\n";
    msg << prefix << "pn_seed: " << self.pn_seed << "\n";
    msg << prefix << "pn_nconf: " << self.pn_nconf << "\n";
    msg << prefix << "fp_iconf_0: " << self.fp_iconf_0 << "\n";
    /* msg << prefix << Helpers<mt::AtomData<T>>::tate(self.get_atoms); */
    msg << prefix << "is_crystal: " << self.is_crystal << "\n";
    msg << prefix << "spec_rot_theta: " << self.spec_rot_theta << "\n";
    /* msg << prefix << "spec_rot_u0: " << self.spec_rot_u0 << "\n"; */
    msg << prefix << "spec_rot_center_type: " << self.spec_rot_center_type << "\n";
    /* msg << prefix << "spec_rot_center_p: " << self.spec_rot_center_p << "\n"; */
    msg << prefix << "thick_type: " << self.thick_type << "\n";
    /* msg << prefix << "thick: " << self.thick << "\n"; */
    msg << prefix << "potential_slicing: " << self.potential_slicing << "\n";
    /* msg << prefix << "grid_2d: " << self.grid_2d << "\n"; */
    /* msg << prefix << "output_area: " << self.output_area << "\n"; */
    msg << prefix << "simulation_type: " << self.simulation_type << "\n";
    msg << prefix << "iw_type: " << self.iw_type << "\n";
    /* msg << prefix << "iw_psi: " << self.iw_psi << "\n"; */
    /* msg << prefix << "iw_x: " << self.iw_x << "\n"; */
    /* msg << prefix << "iw_y: " << self.iw_y << "\n"; */
    msg << prefix << "E_0: " << self.E_0 << "\n";
    msg << prefix << "lambda: " << self.lambda << "\n";
    msg << prefix << "theta: " << self.theta << "\n";
    msg << prefix << "phi: " << self.phi << "\n";
    msg << prefix << "illumination_model: " << self.illumination_model << "\n";
    msg << prefix << "temporal_spatial_incoh: " << self.temporal_spatial_incoh << "\n";
    msg << prefix << "cond_lens" << "\n";
    msg << prefix << to_string(self.cond_lens, " -");
    msg << prefix << "obj_lens" << "\n";
    msg << prefix << to_string(self.obj_lens, " -");
    /* msg << prefix << self.scanning; */
    /* msg << prefix << self.detector; */
    /* msg << prefix << self.eels_fr; */
    msg << prefix << "operation_mode: " << self.operation_mode << "\n";
    msg << prefix << "slice_storage: " << self.slice_storage << "\n";
    msg << prefix << "reverse_multislice: " << self.reverse_multislice << "\n";
    msg << prefix << "mul_sign: " << self.mul_sign << "\n";
    msg << prefix << "Vrl: " << self.Vrl << "\n";
    msg << prefix << "nR: " << self.nR << "\n";
    msg << prefix << "nrot: " << self.nrot << "\n";
    msg << prefix << "cdl_var_type: " << self.cdl_var_type << "\n";
    /* msg << prefix << "cdl_var: " << self.cdl_var << "\n"; */
    /* msg << prefix << "iscan: " << self.iscan << "\n"; */
    /* msg << prefix << "beam_x: " << self.beam_x << "\n"; */
    /* msg << prefix << "beam_y: " << self.beam_y << "\n"; */
    msg << prefix << "islice: " << self.islice << "\n";
    msg << prefix << "dp_Shift: " << self.dp_Shift << "\n";
    return msg.str();
  }
  
  /****************************************************************************
   * Misc function calls
   ***************************************************************************/

  namespace detail {
  
    template <typename T, eDevice dev>
    mt::Output<T> tem_simulation_internal(Input_Multislice<T>& input_multislice) {

      // Ensure we have the correct function
      MULTEM_ASSERT(input_multislice.system_conf.device == dev);

      // Initialise the stream and plan
      mt::Stream<dev> stream(input_multislice.system_conf.nstream);
      mt::FFT<T, dev> fft_2d;
      fft_2d.create_plan_2d(
          input_multislice.grid_2d.ny, 
          input_multislice.grid_2d.nx, 
          input_multislice.system_conf.nstream);

      // Create the multislice structure
      mt::Multislice<T, dev> simulate;
      simulate.set_input_data(&input_multislice, &stream, &fft_2d);

      // Initialise the output data
      mt::Output<T> output_multislice;
      output_multislice.set_input_data(input_multislice);
      simulate(output_multislice.internal().data);
      stream.synchronize();

      // Finalise the output data
      output_multislice.gather();
      output_multislice.clean_temporal();
      fft_2d.cleanup();

      // If there is an error then throw
      auto err = cudaGetLastError();
      if (err != cudaSuccess) {
        std::ostringstream msg;
        msg << "CUDA error: %s\n";
        msg << cudaGetErrorString(err);
        throw std::runtime_error(msg.str());
      }

      // Return the output multislice
      return output_multislice;
    }

  }

  template <typename T>
  mt::Output<T> tem_simulation(Input_Multislice<T>& input_multislice) {
    eDevice dev = input_multislice.system_conf.device;
    return (dev == mt::e_device 
        ? detail::tem_simulation_internal<T, mt::e_device>(input_multislice)
        : detail::tem_simulation_internal<T, mt::e_host>(input_multislice));
  }
  
  /****************************************************************************
   * Explict instantiation of template classes and functions
   ***************************************************************************/
  template class Lens<float>;
  template class Lens<double>;
  
  template class EELS<float>;
  template class EELS<double>;

  template class Atom_Data<float>;
  template class Atom_Data<double>;
  
  template class Input_Multislice<float>;
  template class Input_Multislice<double>;

  template class Output<float>;
  template class Output<double>;

  template mt::Output<float> tem_simulation<float>(Input_Multislice<float>&);
  template mt::Output<double> tem_simulation<double>(Input_Multislice<double>&);

  template std::string mt::to_string<float>(const Lens<float>&, const std::string&);
  template std::string mt::to_string<double>(const Lens<double>&, const std::string&);
  template std::string mt::to_string<float>(const Input_Multislice<float>&, const std::string&);
  template std::string mt::to_string<double>(const Input_Multislice<double>&, const std::string&);

}  // namespace mt

namespace std {
  
  /****************************************************************************
   * Explict instantiation of template classes and functions from std
   ***************************************************************************/

  template class complex<float>;
  template class complex<double>;

  template class vector<float>;
  template class vector<double>;
  template class vector<complex<float>>;
  template class vector<complex<double>>;

}
