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
#include "amorp_spec.hpp"
#include "fxeg_data.hpp"
#include "atomic_data.hpp"

namespace mt {
  
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

	template<class T>
	void rdf_3d(const Atom_Data<T> &atoms, T r_max, int nr, std::vector<T> &r, std::vector<T> &rdf) { 
    Vector<T, e_host> r_temp(r.size());
    Vector<T, e_host> rdf_temp(r.size());
	  rdf_3d(atoms, r_max, nr, r_temp, rdf_temp);
    r.assign(r_temp.begin(), r_temp.end());
    rdf.assign(rdf_temp.begin(), rdf_temp.end());
  }

  namespace detail {
    
    template <typename T, eDevice dev>
    mt::Output_Multislice<T> incident_wave_internal(Input_Multislice<T>& input_multislice) {

      // Check input
      MULTEM_ASSERT(input_multislice.system_conf.device == dev);
      MULTEM_ASSERT(input_multislice.simulation_type == mt::eTEMST_IWRS); 
      
      // Initialise the stream and plan
      mt::Stream<dev> stream(input_multislice.system_conf.nstream);
      mt::FFT<T, dev> fft_2d;
      fft_2d.create_plan_2d(
          input_multislice.grid_2d.ny, 
          input_multislice.grid_2d.nx, 
          input_multislice.system_conf.nstream);

      // Create the incident wave structure
      mt::Incident_Wave<T, dev> incident_wave;
      incident_wave.set_input_data(&input_multislice, &stream, &fft_2d);

      // Initialise the output data
      mt::Output_Multislice<T> output_multislice;
      output_multislice.set_input_data(&input_multislice);
      incident_wave(mt::eS_Real, output_multislice);
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


    template <typename T, eDevice dev>
    mt::Output_Multislice<T> microscope_abberations_internal(Input_Multislice<T>& input_multislice) {

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
	    mt::Microscope_Effects<T, dev> microscope_effects;
	    microscope_effects.set_input_data(&input_multislice, &stream, &fft_2d);

      // Initialise the output data
      mt::Output_Multislice<T> output_multislice;
      output_multislice.set_input_data(&input_multislice);
	    microscope_effects(output_multislice);
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

    template <typename T, eDevice dev>
    mt::Output_Multislice<T> projected_potential_internal(Input_Multislice<T>& input_multislice) {

      // Ensure we have the correct function
      MULTEM_ASSERT(input_multislice.system_conf.device == dev);

      // Initialise the stream and plan
      mt::Stream<dev> stream(input_multislice.system_conf.nstream);

      // Create the multislice structure
	    mt::Projected_Potential<T, dev> projected_potential;
	    projected_potential.set_input_data(&input_multislice, &stream);

      // Initialise the output data
      mt::Output_Multislice<T> output_multislice;
      output_multislice.set_input_data(&input_multislice);
	    projected_potential.move_atoms(input_multislice.pn_nconf);
	    projected_potential(input_multislice.islice, output_multislice);
      stream.synchronize();

      // Finalise the output data
      output_multislice.gather();
      output_multislice.clean_temporal();

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
    
    template <typename T, eDevice dev>
    mt::Output_Multislice<T> propagate_internal(Input_Multislice<T>& input_multislice) {

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
	    mt::Propagator<T, dev> propagator;
	    propagator.set_input_data(&input_multislice, &stream, &fft_2d);

      // Initialise the output data
      mt::Output_Multislice<T> output_multislice;
      output_multislice.set_input_data(&input_multislice);
	    propagator(mt::eS_Real, input_multislice.gx_0(), input_multislice.gy_0(), input_multislice.obj_lens.c_10 , output_multislice);
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
    
    template <typename T>
    std::vector<T> spec_planes_internal(Input_Multislice<T>& input_multislice) {
	    mt::Slicing<T> slicing;
	    slicing.set_input_data(&input_multislice, &(input_multislice.atoms));
      return std::vector<T>(slicing.z_plane.begin(), slicing.z_plane.end());
    }
    
    template <typename T>
    std::tuple<Atom_Data<T>, vector<Slice<T>>> spec_slicing_internal(Input_Multislice<T>& input_multislice) {
	    mt::Spec<T> spec;
	    spec.set_input_data(&input_multislice);
	    spec.move_atoms(input_multislice.pn_nconf);
      std::vector<Slice<T>> slices(spec.slicing.slice.begin(), spec.slicing.slice.end());
      return std::make_tuple(spec.atoms, slices);
    }
    
    template <typename T, eDevice dev>
    mt::Output_Multislice<T> tem_simulation_internal(Input_Multislice<T>& input_multislice) {

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
      mt::Output_Multislice<T> output_multislice;
      output_multislice.set_input_data(&input_multislice);
      simulate(output_multislice);
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

    template <typename T, eDevice dev>
    mt::Output_Multislice<T> transmission_function_internal(Input_Multislice<T>& input_multislice) {

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
	    mt::Transmission_Function<T, dev> transmission_function;
	    transmission_function.set_input_data(&input_multislice, &stream, &fft_2d);

      // Initialise the output data
      mt::Output_Multislice<T> output_multislice;
      output_multislice.set_input_data(&input_multislice);
      transmission_function.move_atoms(input_multislice.pn_nconf);
      transmission_function.trans(input_multislice.islice, output_multislice);
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

    template <typename T, eDevice dev>
    mt::Output_Multislice<T> wave_function_internal(Input_Multislice<T>& input_multislice) {

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
	    mt::Wave_Function<T, dev> wave_function;
      wave_function.set_input_data(&input_multislice, &stream, &fft_2d);

      // Initialise the output data
      mt::Output_Multislice<T> output_multislice;
      output_multislice.set_input_data(&input_multislice);
      wave_function.move_atoms(input_multislice.pn_nconf);
      wave_function.set_incident_wave(wave_function.psi_z);
      wave_function.psi(1.0, wave_function.psi_z, output_multislice);
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
  mt::Output_Multislice<T> incident_wave(Input_Multislice<T>& input_multislice) {
    eDevice dev = input_multislice.system_conf.device;
    return (dev == mt::e_device 
        ? detail::incident_wave_internal<T, mt::e_device>(input_multislice)
        : detail::incident_wave_internal<T, mt::e_host>(input_multislice));
  }
  
  template <typename T>
  mt::Output_Multislice<T> microscope_abberations(Input_Multislice<T>& input_multislice) {
    eDevice dev = input_multislice.system_conf.device;
    return (dev == mt::e_device 
        ? detail::microscope_abberations_internal<T, mt::e_device>(input_multislice)
        : detail::microscope_abberations_internal<T, mt::e_host>(input_multislice));
  }
  
  template <typename T>
  mt::Output_Multislice<T> projected_potential(Input_Multislice<T>& input_multislice) {
    eDevice dev = input_multislice.system_conf.device;
    return (dev == mt::e_device 
        ? detail::projected_potential_internal<T, mt::e_device>(input_multislice)
        : detail::projected_potential_internal<T, mt::e_host>(input_multislice));
  }
  
  template <typename T>
  mt::Output_Multislice<T> propagate(Input_Multislice<T>& input_multislice) {
    eDevice dev = input_multislice.system_conf.device;
    return (dev == mt::e_device 
        ? detail::propagate_internal<T, mt::e_device>(input_multislice)
        : detail::propagate_internal<T, mt::e_host>(input_multislice));
  }
  
  template <typename T>
  std::vector<T> spec_planes(Input_Multislice<T>& input_multislice) {
    return detail::spec_planes_internal(input_multislice);
  }

  template <typename T>
  std::tuple<Atom_Data<T>, vector<Slice<T>>> spec_slicing(Input_Multislice<T>& input_multislice) {
    return detail::spec_slicing_internal(input_multislice);
  }

  template <typename T>
  mt::Output_Multislice<T> tem_simulation(Input_Multislice<T>& input_multislice) {
    eDevice dev = input_multislice.system_conf.device;
    return (dev == mt::e_device 
        ? detail::tem_simulation_internal<T, mt::e_device>(input_multislice)
        : detail::tem_simulation_internal<T, mt::e_host>(input_multislice));
  }
  
  template <typename T>
  mt::Output_Multislice<T> transmission_function(Input_Multislice<T>& input_multislice) {
    eDevice dev = input_multislice.system_conf.device;
    return (dev == mt::e_device 
        ? detail::transmission_function_internal<T, mt::e_device>(input_multislice)
        : detail::transmission_function_internal<T, mt::e_host>(input_multislice));
  }
  
  template <typename T>
  mt::Output_Multislice<T> wave_function(Input_Multislice<T>& input_multislice) {
    eDevice dev = input_multislice.system_conf.device;
    return (dev == mt::e_device 
        ? detail::wave_function_internal<T, mt::e_device>(input_multislice)
        : detail::wave_function_internal<T, mt::e_host>(input_multislice));
  }
  
  /****************************************************************************
   * Explict instantiation of template classes and functions
   ***************************************************************************/

  template class Lens<float>;
  template class Lens<double>;
  
  template class EELS<float>;
  template class EELS<double>;

  template class Atom_Cal<float>;
  template class Atom_Cal<double>;
  
  template class Amorp_Spec<float>;
  template class Amorp_Spec<double>;

  template class Atom_Data<float>;
  template class Atom_Data<double>;
  
  template class Input_Multislice<float>;
  template class Input_Multislice<double>;

  template class Output_Multislice<float>;
  template class Output_Multislice<double>;

  template void Atomic_Data::To_atom_type_CPU<float, e_host>(
      int, double, int, double, Atom_Type<float, e_host>&);
  template void Atomic_Data::To_atom_type_CPU<double, e_host>(
      int, double, int, double, Atom_Type<double, e_host>&);
  
  template Output_Multislice<float> incident_wave<float>(Input_Multislice<float>&);
  template Output_Multislice<double> incident_wave<double>(Input_Multislice<double>&);
  
  template Output_Multislice<float> microscope_abberations<float>(Input_Multislice<float>&);
  template Output_Multislice<double> microscope_abberations<double>(Input_Multislice<double>&);
  
  template Output_Multislice<float> projected_potential<float>(Input_Multislice<float>&);
  template Output_Multislice<double> projected_potential<double>(Input_Multislice<double>&);
  
  template Output_Multislice<float> propagate<float>(Input_Multislice<float>&);
  template Output_Multislice<double> propagate<double>(Input_Multislice<double>&);
  
  template std::vector<float> spec_planes<float>(Input_Multislice<float>&);
  template std::vector<double> spec_planes<double>(Input_Multislice<double>&);
  
  template std::tuple<Atom_Data<float>, std::vector<Slice<float>>> spec_slicing<float>(
      Input_Multislice<float>&);
  template std::tuple<Atom_Data<double>, std::vector<Slice<double>>> spec_slicing<double>(
      Input_Multislice<double>&);
  
  template Output_Multislice<float> tem_simulation<float>(Input_Multislice<float>&);
  template Output_Multislice<double> tem_simulation<double>(Input_Multislice<double>&);
  
  template Output_Multislice<float> transmission_function<float>(Input_Multislice<float>&);
  template Output_Multislice<double> transmission_function<double>(Input_Multislice<double>&);

  template Output_Multislice<float> wave_function<float>(Input_Multislice<float>&);
  template Output_Multislice<double> wave_function<double>(Input_Multislice<double>&);

	template float get_lambda<float>(const float &E_0);
	template double get_lambda<double>(const double &E_0);

	template float get_sigma<float>(const float &E_0);
	template double get_sigma<double>(const double &E_0);

	template float get_gamma<float>(const float &E_0);
	template double get_gamma<double>(const double &E_0);

	template float rad_2_rAngs<float>(const float &E_0, const float &theta);
	template double rad_2_rAngs<double>(const double &E_0, const double &theta);

	template float hwhm_2_sigma<float>(const float &v);
	template double hwhm_2_sigma<double>(const double &v);

	template float fwhm_2_sigma<float>(const float &v);
	template double fwhm_2_sigma<double>(const double &v);

	template float iehwgd_2_sigma<float>(const float &v);
	template double iehwgd_2_sigma<double>(const double &v);

	template float rad_2_sigma<float>(const float &E_0, const float &theta);
	template double rad_2_sigma<double>(const double &E_0, const double &theta);

	template float get_Vr_factor<float>(const float &E_0, const float &theta);
	template double get_Vr_factor<double>(const double &E_0, const double &theta);

	template float get_Scherzer_defocus<float>(const float &E_0, const float &c_30);
	template double get_Scherzer_defocus<double>(const double &E_0, const double &c_30);

	template float get_Scherzer_aperture<float>(const float &E_0, const float &c_30);
	template double get_Scherzer_aperture<double>(const double &E_0, const double &c_30);

	template void get_Scherzer_conditions<float>(const float &E_0, const float &c_30, float &defocus, float &aperture);
	template void get_Scherzer_conditions<double>(const double &E_0, const double &c_30, double &defocus, double &aperture);

	template void rdf_3d<float>(const Atom_Data<float> &atoms, float r_max, int nr, std::vector<float> &r, std::vector<float> &rdf);
	template void rdf_3d<double>(const Atom_Data<double> &atoms, double r_max, int nr, std::vector<double> &r, std::vector<double> &rdf);
	
  template void rotate_atoms<float>(Atom_Data<float> &atoms, float theta, r3d<float> u0, r3d<float> p0);
  template void rotate_atoms<double>(Atom_Data<double> &atoms, double theta, r3d<double> u0, r3d<double> p0);
	
  template void remove_atoms_outside_z_range<float>(Atom_Data<float> &atoms, float z_0, float z_e);
  template void remove_atoms_outside_z_range<double>(Atom_Data<double> &atoms, double z_0, double z_e);

  template std::string to_string<float>(const Lens<float>&, const std::string&);
  template std::string to_string<double>(const Lens<double>&, const std::string&);
  template std::string to_string<float>(const Input_Multislice<float>&, const std::string&);
  template std::string to_string<double>(const Input_Multislice<double>&, const std::string&);

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
  template class vector<vector<float>>;
  template class vector<vector<double>>;
  template class vector<vector<complex<float>>>;
  template class vector<vector<complex<double>>>;
  template class vector<mt::Det_Int<vector<float>>>;
  template class vector<mt::Det_Int<vector<double>>>;

}
