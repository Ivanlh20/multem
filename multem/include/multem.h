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

namespace mt {

  /**
   * A class to represent atom data
   */
  class Atom {
  public:
      int element;
      double x;
      double y;
      double z;
      double sigma;
      double occupancy;
      int region;
      int charge;
      
      Atom()
        : element(0),
          x(0),
          y(0),
          z(0),
          sigma(0),
          occupancy(0),
          region(0),
          charge(0) {}

      Atom(
            int element_, 
            double x_, 
            double y_, 
            double z_, 
            double sigma_, 
            double occupancy_, 
            int region_, 
            int charge_)
        : element(element_),
          x(x_),
          y(y_),
          z(z_),
          sigma(sigma_),
          occupancy(occupancy_),
          region(region_),
          charge(charge_) {}
  };

  /**
   * A class to hold crystal parameters
   */
  class CrystalParameters {
  public:
    
    typedef std::vector<Atom> Layer;

    int na;
    int nb;
    int nc;
    double a;
    double b;
    double c;
    std::vector<Layer> layers;

    CrystalParameters()
      : na(0),
        nb(0),
        nc(0),
        a(0),
        b(0),
        c(0) {}

  };

  /**
   * A class to hold amorphous layer information
   */
  class AmorphousLayer {
  public:

    double z_0;
    double z_e;
    double dz;

    AmorphousLayer()
      : z_0(0),
        z_e(0),
        dz(0) {}

  };

  /**
   * A class to hold STEM detector information
   */
  class STEMDetector {
  public:

    class Angles {
    public:
      double inner_ang;
      double outer_ang;
      Angles()
        : inner_ang(0),
          outer_ang(0) {}
    };

    class Radial {
    public:
      double x;
      std::vector<double> fx;
      Radial()
        : x(0) {}
    };

    class Matrix {
    public:
      double R;
      std::vector<double> fR;
      Matrix()
        : R(0) {}
    };
  
    std::string type;
    std::vector<Angles> cir;
    std::vector<Radial> radial;
    std::vector<Matrix> matrix;
  };

  /**
   * A class to hold the input specification for the simulation
   */
  class Input {
  public:

    // Electron-specimen interaction model 
    std::string interaction_model; 
    std::string potential_type;  
    std::string operation_mode;
    std::size_t memory_size;
    bool reverse_multislice;
    
    // Electron-Phonon interaction model
    std::string pn_model;
    bool pn_coh_contrib;
    bool pn_single_conf;
    int pn_nconf;
    int pn_dim;
    int pn_seed;

    // Specimen information
    std::vector<Atom> spec_atoms;
    double spec_dz;
    double spec_lx;
    double spec_ly;
    double spec_lz;
    int spec_cryst_na;
    int spec_cryst_nb;
    int spec_cryst_nc;
    double spec_cryst_a;
    double spec_cryst_b;
    double spec_cryst_c;
    double spec_cryst_x0;
    double spec_cryst_y0;
    std::vector<AmorphousLayer> spec_amorp;

    // Specimen rotation
    double spec_rot_theta;
    std::array<double, 3> spec_rot_u0;
    std::string spec_rot_center_type;
    std::array<double, 3> spec_rot_center_p;

    // Specimen thickness
    std::string thick_type;
    std::vector<double> thick;

    // Potential slicing
    std::string potential_slicing;

    // X-Y sampling
    int nx;
    int ny;
    bool bwl;

    // Simulation type
    std::string simulation_type;

    // Incident wave
    std::string iw_type;
    std::vector< std::complex<double> > iw_psi;
    std::vector<double> iw_x;
    std::vector<double> iw_y;

    // Microscope parameters
    double E_0;
    double theta;
    double phi;

    // Illumination model
    std::string illumination_model;
    std::string temporal_spatial_incoh;

    // Condenser lens
    int cond_lens_m;
    double cond_lens_c_10;
    double cond_lens_c_12;
    double cond_lens_phi_12;
    double cond_lens_c_21;
    double cond_lens_phi_21;
    double cond_lens_c_23;
    double cond_lens_phi_23;
    double cond_lens_c_30;
    double cond_lens_c_32;
    double cond_lens_phi_32;
    double cond_lens_c_34;
    double cond_lens_phi_34;
    double cond_lens_c_41;
    double cond_lens_phi_41;
    double cond_lens_c_43;
    double cond_lens_phi_43;
    double cond_lens_c_45;
    double cond_lens_phi_45;
    double cond_lens_c_50;
    double cond_lens_c_52;
    double cond_lens_phi_52;
    double cond_lens_c_54;
    double cond_lens_phi_54;
    double cond_lens_c_56;
    double cond_lens_phi_56;
    double cond_lens_inner_aper_ang;
    double cond_lens_outer_aper_ang;

    // Source spread function
    double cond_lens_ssf_sigma;
    int cond_lens_ssf_npoints;
    
    // Defocus spread function
    double cond_lens_dsf_sigma;
    int cond_lens_dsf_npoints;

    // Zero defocus reference
    std::string cond_lens_zero_defocus_type;
    double cond_lens_zero_defocus_plane;

    // Objective lens
    int obj_lens_m;
    double obj_lens_c_10;
    double obj_lens_c_12;
    double obj_lens_phi_12;
    double obj_lens_c_21;
    double obj_lens_phi_21;
    double obj_lens_c_23;
    double obj_lens_phi_23;
    double obj_lens_c_30;
    double obj_lens_c_32;
    double obj_lens_phi_32;
    double obj_lens_c_34;
    double obj_lens_phi_34;
    double obj_lens_c_41;
    double obj_lens_phi_41;
    double obj_lens_c_43;
    double obj_lens_phi_43;
    double obj_lens_c_45;
    double obj_lens_phi_45;
    double obj_lens_c_50;
    double obj_lens_c_52;
    double obj_lens_phi_52;
    double obj_lens_c_54;
    double obj_lens_phi_54;
    double obj_lens_c_56;
    double obj_lens_phi_56;
    double obj_lens_inner_aper_ang;
    double obj_lens_outer_aper_ang;

    // Defocus spread function
    double obj_lens_dsf_sigma;
    int obj_lens_dsf_npoints;

    // The phase shift
    double phase_shift;

    // Zero defocus reference
    std::string obj_lens_zero_defocus_type;
    double obj_lens_zero_defocus_plane;

    // STEM detector
    STEMDetector detector;

    // Scanning area for ISTEM/STEM/EELS
    std::string scanning_type;
    bool scanning_periodic;
    int scanning_ns;
    double scanning_x0;
    double scanning_y0;
    double scanning_xe;
    double scanning_ye;

    // PED
    double ped_nrot;
    double ped_theta;

    // HCI
    double hci_nrot;
    double hci_theta;

    // EELS
    int eels_Z;
    double eels_E_loss;
    double eels_collection_angle;
    int eels_m_selection;
    std::string eels_channelling_type;

    // EFTEM
    int eftem_Z;
    double eftem_E_loss;
    double eftem_collection_angle;
    int eftem_m_selection;
    std::string eftem_channelling_type;

    // Output region
    int output_area_ix_0;
    int output_area_iy_0;
    int output_area_ix_e;
    int output_area_iy_e;

    /**
     * Set the default parameters
     */
    Input()
      : interaction_model("Multislice"), 
        potential_type("Lobato_0_12"),
        operation_mode("Normal"),
        memory_size(0),
        reverse_multislice(false),
        pn_model("Still_Atom"),
        pn_coh_contrib(false),
        pn_single_conf(false),
        pn_nconf(1),
        pn_dim(110),
        pn_seed(300183),
        spec_dz(0.25),
        spec_lx(10),
        spec_ly(10),
        spec_lz(10),
        spec_cryst_na(1),
        spec_cryst_nb(1),
        spec_cryst_nc(1),
        spec_cryst_a(0),
        spec_cryst_b(0),
        spec_cryst_c(0),
        spec_cryst_x0(0),
        spec_cryst_y0(0),
        spec_rot_theta(0),
        spec_rot_u0({0, 0, 1}),
        spec_rot_center_type("geometric_center"),
        spec_rot_center_p({0, 0, 0}),
        thick_type("Whole_Spec"),
        thick(0),
        potential_slicing("Planes"),
        nx(256),
        ny(256),
        bwl(false),
        simulation_type("EWRS"),
        iw_type("Auto"),
        iw_psi(0),
        iw_x(0.0),
        iw_y(0.0),
        E_0(300.0),
        theta(0.0),
        phi(0.0),
        illumination_model("Partial_Coherent"),
        temporal_spatial_incoh("Temporal_Spatial"),
        cond_lens_m(0),
        cond_lens_c_10(14.0312),
        cond_lens_c_12(0.0),
        cond_lens_phi_12(0.0),
        cond_lens_c_21(0.0),
        cond_lens_phi_21(0.0),
        cond_lens_c_23(0.0),
        cond_lens_phi_23(0.0),
        cond_lens_c_30(1e-3),
        cond_lens_c_32(0.0),
        cond_lens_phi_32(0.0),
        cond_lens_c_34(0.0),
        cond_lens_phi_34(0.0),
        cond_lens_c_41(0.0),
        cond_lens_phi_41(0.0),
        cond_lens_c_43(0.0),
        cond_lens_phi_43(0.0),
        cond_lens_c_45(0.0),
        cond_lens_phi_45(0.0),
        cond_lens_c_50(0.0),
        cond_lens_c_52(0.0),
        cond_lens_phi_52(0.0),
        cond_lens_c_54(0.0),
        cond_lens_phi_54(0.0),
        cond_lens_c_56(0.0),
        cond_lens_phi_56(0.0),
        cond_lens_inner_aper_ang(0.0),
        cond_lens_outer_aper_ang(21.0),
        cond_lens_ssf_sigma(0.0072),
        cond_lens_ssf_npoints(8),
        cond_lens_dsf_sigma(32),
        cond_lens_dsf_npoints(10),
        cond_lens_zero_defocus_type("Last"),
        cond_lens_zero_defocus_plane(0),
        obj_lens_m(0),
        obj_lens_c_10(14.0312),
        obj_lens_c_12(0.0),
        obj_lens_phi_12(0.0),
        obj_lens_c_21(0.0),
        obj_lens_phi_21(0.0),
        obj_lens_c_23(0.0),
        obj_lens_phi_23(0.0),
        obj_lens_c_30(1e-3),
        obj_lens_c_32(0.0),
        obj_lens_phi_32(0.0),
        obj_lens_c_34(0.0),
        obj_lens_phi_34(0.0),
        obj_lens_c_41(0.0),
        obj_lens_phi_41(0.0),
        obj_lens_c_43(0.0),
        obj_lens_phi_43(0.0),
        obj_lens_c_45(0.0),
        obj_lens_phi_45(0.0),
        obj_lens_c_50(0.0),
        obj_lens_c_52(0.0),
        obj_lens_phi_52(0.0),
        obj_lens_c_54(0.0),
        obj_lens_phi_54(0.0),
        obj_lens_c_56(0.0),
        obj_lens_phi_56(0.0),
        obj_lens_inner_aper_ang(0.0),
        obj_lens_outer_aper_ang(24.0),
        obj_lens_dsf_sigma(32),
        obj_lens_dsf_npoints(10),
        obj_lens_zero_defocus_type("Last"),
        obj_lens_zero_defocus_plane(0),
        phase_shift(0),
        scanning_type("Line"),
        scanning_periodic(true),
        scanning_ns(10),
        scanning_x0(0.0),
        scanning_y0(0.0),
        scanning_xe(4.078),
        scanning_ye(4.078),
        ped_nrot(360),
        ped_theta(3.0),
        hci_nrot(360),
        hci_theta(3.0),
        eels_Z(79),
        eels_E_loss(80),
        eels_collection_angle(100),
        eels_m_selection(3),
        eels_channelling_type("Single_Channelling"),
        eftem_Z(79),
        eftem_E_loss(80),
        eftem_collection_angle(100),
        eftem_m_selection(3),
        eftem_channelling_type("Single_Channelling"),
        output_area_ix_0(1),
        output_area_iy_0(1),
        output_area_ix_e(1),
        output_area_iy_e(1) {}
  };

  /**
   * A class to hold the system configuration
   */
  class SystemConfiguration {
  public:

    std::string device;
    std::string precision;
    std::size_t cpu_ncores;
    std::size_t cpu_nthread;
    std::size_t gpu_device;
    std::size_t gpu_nstream;
    
    SystemConfiguration()
      : device("device"),
        precision("float"),
        cpu_ncores(1),
        cpu_nthread(1),
        gpu_device(0),
        gpu_nstream(1) {}

  };

  void test(const Input &a);

}

#endif

