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

  /**
   * A class to represent atom data
   */
  /*class Atom {*/
  /*public:*/
  /*    int element;*/
  /*    double x;*/
  /*    double y;*/
  /*    double z;*/
  /*    double sigma;*/
  /*    double occupancy;*/
  /*    int region;*/
  /*    int charge;*/

  /*    Atom()*/
  /*      : element(0),*/
  /*        x(0),*/
  /*        y(0),*/
  /*        z(0),*/
  /*        sigma(0),*/
  /*        occupancy(0),*/
  /*        region(0),*/
  /*        charge(0) {}*/

  /*    Atom(*/
  /*          int element_,*/
  /*          double x_,*/
  /*          double y_,*/
  /*          double z_,*/
  /*          double sigma_,*/
  /*          double occupancy_,*/
  /*          int region_,*/
  /*          int charge_)*/
  /*      : element(element_),*/
  /*        x(x_),*/
  /*        y(y_),*/
  /*        z(z_),*/
  /*        sigma(sigma_),*/
  /*        occupancy(occupancy_),*/
  /*        region(region_),*/
  /*        charge(charge_) {}*/
  /*};*/

  /**
  * A class to hold crystal parameters
  */
  /*class CrystalParameters {*/
  /*public:*/

  /*  typedef std::vector<Atom> Layer;*/

  /*  int na;*/
  /*  int nb;*/
  /*  int nc;*/
  /*  double a;*/
  /*  double b;*/
  /*  double c;*/
  /*  std::vector<Layer> layers;*/

  /*  CrystalParameters()*/
  /*    : na(0),*/
  /*      nb(0),*/
  /*      nc(0),*/
  /*      a(0),*/
  /*      b(0),*/
  /*      c(0) {}*/

  /*};*/

  /**
   * A class to hold amorphous layer information
   */
  /*class AmorphousLayer {*/
  /*public:*/

  /*  double z_0;*/
  /*  double z_e;*/
  /*  double dz;*/

  /*  AmorphousLayer()*/
  /*    : z_0(0),*/
  /*      z_e(0),*/
  /*      dz(0) {}*/

  /*};*/

  /**
   * A class to hold STEM detector information
   */
  /*class STEMDetector {*/
  /*public:*/

  /*  class Angles {*/
  /*  public:*/
  /*    double inner_ang;*/
  /*    double outer_ang;*/
  /*    Angles()*/
  /*      : inner_ang(0),*/
  /*        outer_ang(0) {}*/
  /*  };*/

  /*  class Radial {*/
  /*  public:*/
  /*    double x;*/
  /*    std::vector<double> fx;*/
  /*    Radial()*/
  /*      : x(0) {}*/
  /*  };*/

  /*  class Matrix {*/
  /*  public:*/
  /*    double R;*/
  /*    std::vector<double> fR;*/
  /*    Matrix()*/
  /*      : R(0) {}*/
  /*  };*/

  /*  std::string type;*/
  /*  std::vector<Angles> cir;*/
  /*  std::vector<Radial> radial;*/
  /*  std::vector<Matrix> matrix;*/
  /*};*/

  /**
   * A class to hold the input specification for the simulation
   */
  /* class Input { */
  /* public: */

  /*// Electron-specimen interaction model*/
  /*std::string interaction_model;*/
  /*std::string potential_type;*/
  /*std::string operation_mode;*/
  /*std::size_t memory_size;*/
  /*bool reverse_multislice;*/

  /*// Electron-Phonon interaction model*/
  /*std::string pn_model;*/
  /*bool pn_coh_contrib;*/
  /*bool pn_single_conf;*/
  /*int pn_nconf;*/
  /*int pn_dim;*/
  /*int pn_seed;*/

  /*// Specimen information*/
  /*std::vector<Atom> spec_atoms;*/
  /*double spec_dz;*/
  /*double spec_lx;*/
  /*double spec_ly;*/
  /*double spec_lz;*/
  /*int spec_cryst_na;*/
  /*int spec_cryst_nb;*/
  /*int spec_cryst_nc;*/
  /*double spec_cryst_a;*/
  /*double spec_cryst_b;*/
  /*double spec_cryst_c;*/
  /*double spec_cryst_x0;*/
  /*double spec_cryst_y0;*/
  /*std::vector<AmorphousLayer> spec_amorp;*/

  /*// Specimen rotation*/
  /*double spec_rot_theta;*/
  /*std::array<double, 3> spec_rot_u0;*/
  /*std::string spec_rot_center_type;*/
  /*std::array<double, 3> spec_rot_center_p;*/

  /*// Specimen thickness*/
  /*std::string thick_type;*/
  /*std::vector<double> thick;*/

  /*// Potential slicing*/
  /*std::string potential_slicing;*/

  /*// X-Y sampling*/
  /*int nx;*/
  /*int ny;*/
  /*bool bwl;*/

  /*// Simulation type*/
  /*std::string simulation_type;*/

  /*// Incident wave*/
  /*std::string iw_type;*/
  /*std::vector< std::complex<double> > iw_psi;*/
  /*std::vector<double> iw_x;*/
  /*std::vector<double> iw_y;*/

  /*// Microscope parameters*/
  /*double E_0;*/
  /*double theta;*/
  /*double phi;*/

  /*// Illumination model*/
  /*std::string illumination_model;*/
  /*std::string temporal_spatial_incoh;*/

  /*// Condenser lens*/
  /*int cond_lens_m;*/
  /*double cond_lens_c_10;*/
  /*double cond_lens_c_12;*/
  /*double cond_lens_phi_12;*/
  /*double cond_lens_c_21;*/
  /*double cond_lens_phi_21;*/
  /*double cond_lens_c_23;*/
  /*double cond_lens_phi_23;*/
  /*double cond_lens_c_30;*/
  /*double cond_lens_c_32;*/
  /*double cond_lens_phi_32;*/
  /*double cond_lens_c_34;*/
  /*double cond_lens_phi_34;*/
  /*double cond_lens_c_41;*/
  /*double cond_lens_phi_41;*/
  /*double cond_lens_c_43;*/
  /*double cond_lens_phi_43;*/
  /*double cond_lens_c_45;*/
  /*double cond_lens_phi_45;*/
  /*double cond_lens_c_50;*/
  /*double cond_lens_c_52;*/
  /*double cond_lens_phi_52;*/
  /*double cond_lens_c_54;*/
  /*double cond_lens_phi_54;*/
  /*double cond_lens_c_56;*/
  /*double cond_lens_phi_56;*/
  /*double cond_lens_inner_aper_ang;*/
  /*double cond_lens_outer_aper_ang;*/

  /*// Source spread function*/
  /*double cond_lens_ssf_sigma;*/
  /*int cond_lens_ssf_npoints;*/

  /*// Defocus spread function*/
  /*double cond_lens_dsf_sigma;*/
  /*int cond_lens_dsf_npoints;*/

  /*// Zero defocus reference*/
  /*std::string cond_lens_zero_defocus_type;*/
  /*double cond_lens_zero_defocus_plane;*/

  /*// Objective lens*/
  /*int obj_lens_m;*/
  /*double obj_lens_c_10;*/
  /*double obj_lens_c_12;*/
  /*double obj_lens_phi_12;*/
  /*double obj_lens_c_21;*/
  /*double obj_lens_phi_21;*/
  /*double obj_lens_c_23;*/
  /*double obj_lens_phi_23;*/
  /*double obj_lens_c_30;*/
  /*double obj_lens_c_32;*/
  /*double obj_lens_phi_32;*/
  /*double obj_lens_c_34;*/
  /*double obj_lens_phi_34;*/
  /*double obj_lens_c_41;*/
  /*double obj_lens_phi_41;*/
  /*double obj_lens_c_43;*/
  /*double obj_lens_phi_43;*/
  /*double obj_lens_c_45;*/
  /*double obj_lens_phi_45;*/
  /*double obj_lens_c_50;*/
  /*double obj_lens_c_52;*/
  /*double obj_lens_phi_52;*/
  /*double obj_lens_c_54;*/
  /*double obj_lens_phi_54;*/
  /*double obj_lens_c_56;*/
  /*double obj_lens_phi_56;*/
  /*double obj_lens_inner_aper_ang;*/
  /*double obj_lens_outer_aper_ang;*/

  /*// Defocus spread function*/
  /*double obj_lens_dsf_sigma;*/
  /*int obj_lens_dsf_npoints;*/

  /*// The phase shift*/
  /*double phase_shift;*/

  /*// Zero defocus reference*/
  /*std::string obj_lens_zero_defocus_type;*/
  /*double obj_lens_zero_defocus_plane;*/

  /*// STEM detector*/
  /*STEMDetector detector;*/

  /*// Scanning area for ISTEM/STEM/EELS*/
  /*std::string scanning_type;*/
  /*bool scanning_periodic;*/
  /*int scanning_ns;*/
  /*double scanning_x0;*/
  /*double scanning_y0;*/
  /*double scanning_xe;*/
  /*double scanning_ye;*/

  /*// PED*/
  /*double ped_nrot;*/
  /*double ped_theta;*/

  /*// HCI*/
  /*double hci_nrot;*/
  /*double hci_theta;*/

  /*// EELS*/
  /*int eels_Z;*/
  /*double eels_E_loss;*/
  /*double eels_collection_angle;*/
  /*int eels_m_selection;*/
  /*std::string eels_channelling_type;*/

  /*// EFTEM*/
  /*int eftem_Z;*/
  /*double eftem_E_loss;*/
  /*double eftem_collection_angle;*/
  /*int eftem_m_selection;*/
  /*std::string eftem_channelling_type;*/

  /*// Output region*/
  /*int output_area_ix_0;*/
  /*int output_area_iy_0;*/
  /*int output_area_ix_e;*/
  /*int output_area_iy_e;*/

  /**
     * Set the default parameters
     */
  /*Input()*/
  /*  : interaction_model("Multislice"),*/
  /*    potential_type("Lobato_0_12"),*/
  /*    operation_mode("Normal"),*/
  /*    memory_size(0),*/
  /*    reverse_multislice(false),*/
  /*    pn_model("Still_Atom"),*/
  /*    pn_coh_contrib(false),*/
  /*    pn_single_conf(false),*/
  /*    pn_nconf(1),*/
  /*    pn_dim(110),*/
  /*    pn_seed(300183),*/
  /*    spec_dz(0.25),*/
  /*    spec_lx(10),*/
  /*    spec_ly(10),*/
  /*    spec_lz(10),*/
  /*    spec_cryst_na(1),*/
  /*    spec_cryst_nb(1),*/
  /*    spec_cryst_nc(1),*/
  /*    spec_cryst_a(0),*/
  /*    spec_cryst_b(0),*/
  /*    spec_cryst_c(0),*/
  /*    spec_cryst_x0(0),*/
  /*    spec_cryst_y0(0),*/
  /*    spec_rot_theta(0),*/
  /*    spec_rot_u0({0, 0, 1}),*/
  /*    spec_rot_center_type("geometric_center"),*/
  /*    spec_rot_center_p({0, 0, 0}),*/
  /*    thick_type("Whole_Spec"),*/
  /*    thick(0),*/
  /*    potential_slicing("Planes"),*/
  /*    nx(256),*/
  /*    ny(256),*/
  /*    bwl(false),*/
  /*    simulation_type("EWRS"),*/
  /*    iw_type("Auto"),*/
  /*    iw_psi(0),*/
  /*    iw_x(0.0),*/
  /*    iw_y(0.0),*/
  /*    E_0(300.0),*/
  /*    theta(0.0),*/
  /*    phi(0.0),*/
  /*    illumination_model("Partial_Coherent"),*/
  /*    temporal_spatial_incoh("Temporal_Spatial"),*/
  /*    cond_lens_m(0),*/
  /*    cond_lens_c_10(14.0312),*/
  /*    cond_lens_c_12(0.0),*/
  /*    cond_lens_phi_12(0.0),*/
  /*    cond_lens_c_21(0.0),*/
  /*    cond_lens_phi_21(0.0),*/
  /*    cond_lens_c_23(0.0),*/
  /*    cond_lens_phi_23(0.0),*/
  /*    cond_lens_c_30(1e-3),*/
  /*    cond_lens_c_32(0.0),*/
  /*    cond_lens_phi_32(0.0),*/
  /*    cond_lens_c_34(0.0),*/
  /*    cond_lens_phi_34(0.0),*/
  /*    cond_lens_c_41(0.0),*/
  /*    cond_lens_phi_41(0.0),*/
  /*    cond_lens_c_43(0.0),*/
  /*    cond_lens_phi_43(0.0),*/
  /*    cond_lens_c_45(0.0),*/
  /*    cond_lens_phi_45(0.0),*/
  /*    cond_lens_c_50(0.0),*/
  /*    cond_lens_c_52(0.0),*/
  /*    cond_lens_phi_52(0.0),*/
  /*    cond_lens_c_54(0.0),*/
  /*    cond_lens_phi_54(0.0),*/
  /*    cond_lens_c_56(0.0),*/
  /*    cond_lens_phi_56(0.0),*/
  /*    cond_lens_inner_aper_ang(0.0),*/
  /*    cond_lens_outer_aper_ang(21.0),*/
  /*    cond_lens_ssf_sigma(0.0072),*/
  /*    cond_lens_ssf_npoints(8),*/
  /*    cond_lens_dsf_sigma(32),*/
  /*    cond_lens_dsf_npoints(10),*/
  /*    cond_lens_zero_defocus_type("Last"),*/
  /*    cond_lens_zero_defocus_plane(0),*/
  /*    obj_lens_m(0),*/
  /*    obj_lens_c_10(14.0312),*/
  /*    obj_lens_c_12(0.0),*/
  /*    obj_lens_phi_12(0.0),*/
  /*    obj_lens_c_21(0.0),*/
  /*    obj_lens_phi_21(0.0),*/
  /*    obj_lens_c_23(0.0),*/
  /*    obj_lens_phi_23(0.0),*/
  /*    obj_lens_c_30(1e-3),*/
  /*    obj_lens_c_32(0.0),*/
  /*    obj_lens_phi_32(0.0),*/
  /*    obj_lens_c_34(0.0),*/
  /*    obj_lens_phi_34(0.0),*/
  /*    obj_lens_c_41(0.0),*/
  /*    obj_lens_phi_41(0.0),*/
  /*    obj_lens_c_43(0.0),*/
  /*    obj_lens_phi_43(0.0),*/
  /*    obj_lens_c_45(0.0),*/
  /*    obj_lens_phi_45(0.0),*/
  /*    obj_lens_c_50(0.0),*/
  /*    obj_lens_c_52(0.0),*/
  /*    obj_lens_phi_52(0.0),*/
  /*    obj_lens_c_54(0.0),*/
  /*    obj_lens_phi_54(0.0),*/
  /*    obj_lens_c_56(0.0),*/
  /*    obj_lens_phi_56(0.0),*/
  /*    obj_lens_inner_aper_ang(0.0),*/
  /*    obj_lens_outer_aper_ang(24.0),*/
  /*    obj_lens_dsf_sigma(32),*/
  /*    obj_lens_dsf_npoints(10),*/
  /*    obj_lens_zero_defocus_type("Last"),*/
  /*    obj_lens_zero_defocus_plane(0),*/
  /*    phase_shift(0),*/
  /*    scanning_type("Line"),*/
  /*    scanning_periodic(true),*/
  /*    scanning_ns(10),*/
  /*    scanning_x0(0.0),*/
  /*    scanning_y0(0.0),*/
  /*    scanning_xe(4.078),*/
  /*    scanning_ye(4.078),*/
  /*    ped_nrot(360),*/
  /*    ped_theta(3.0),*/
  /*    hci_nrot(360),*/
  /*    hci_theta(3.0),*/
  /*    eels_Z(79),*/
  /*    eels_E_loss(80),*/
  /*    eels_collection_angle(100),*/
  /*    eels_m_selection(3),*/
  /*    eels_channelling_type("Single_Channelling"),*/
  /*    eftem_Z(79),*/
  /*    eftem_E_loss(80),*/
  /*    eftem_collection_angle(100),*/
  /*    eftem_m_selection(3),*/
  /*    eftem_channelling_type("Single_Channelling"),*/
  /*    output_area_ix_0(1),*/
  /*    output_area_iy_0(1),*/
  /*    output_area_ix_e(1),*/
  /*    output_area_iy_e(1) {}*/
  /*};*/
  
  /****************************************************************************
   * The FFT interface
   ***************************************************************************/

  template <typename T, eDevice dev>
  class FFTData {
  public:
    struct Data;

    FFTData();
    FFTData(const FFTData&);
    FFTData(FFTData&&);
    FFTData(const Data&);
    FFTData& operator=(const FFTData&);
    FFTData& operator=(FFTData&&);
    ~FFTData();
    
    const FFTData::Data& internal() const;
    FFTData::Data& internal();

    void cleanup();
    void destroy_plan();
		void create_plan_1d(const int &nx, int nThread=1);
		void create_plan_1d_batch(const int &ny, const int &nx, int nThread=1);
		void create_plan_2d(const int &ny, const int &nx, int nThread);
		void create_plan_2d_batch(const int &ny, const int &nx, const int &nz, int nThread);
  
  private:
    std::unique_ptr<Data> impl_;
  };
  
  /****************************************************************************
   * The Stream interface
   ***************************************************************************/
  
  template <eDevice dev>
  class StreamIface {
  public:
    struct Data;
    
    StreamIface();
    ~StreamIface();
    
    const StreamIface::Data& internal() const;
    StreamIface::Data& internal();

    StreamIface(int new_stream);
    int size() const;
		void resize(int new_nstream);
		void synchronize();
		void set_n_act_stream(const int &new_n_act_stream);
		void set_grid(const int &nx_i, const int &ny_i);
		Range_2d get_range(const int &istream);

  private:
    std::unique_ptr<Data> impl_;
  };

  /****************************************************************************
   * The SystemConfiguration interface
   ***************************************************************************/

  /**
   * A proxy for the System_Configuration class
   */
  class SystemConfiguration {
  public:
    struct Data;

    SystemConfiguration();
    SystemConfiguration(const SystemConfiguration&);
    SystemConfiguration(SystemConfiguration&&);
    SystemConfiguration(const Data&);
    SystemConfiguration& operator=(const SystemConfiguration&);
    SystemConfiguration& operator=(SystemConfiguration&&);
    ~SystemConfiguration();

    const SystemConfiguration::Data& internal() const;

    void set_precision(ePrecision precision);
    ePrecision get_precision() const;

    void set_device(eDevice device);
    eDevice get_device() const;

    void set_cpu_ncores(int cpu_ncores);
    int get_cpu_ncores() const;

    void set_cpu_nthread(int cpu_nthread);
    int get_cpu_nthread() const;

    void set_gpu_device(int gpu_device);
    int get_gpu_device() const;

    void set_gpu_nstream(int gpu_nstream);
    int get_gpu_nstream() const;

    void set_nstream(int nstream);
    int get_nstream() const;

    void set_active(bool active);
    bool get_active() const;

    bool is_host() const;
    bool is_device() const;
    bool is_float() const;
    bool is_double() const;
    bool is_float_host() const;
    bool is_double_host() const;
    bool is_float_device() const;
    bool is_double_device() const;

  private:
    std::unique_ptr<Data> impl_;
  };
  
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

    AtomData();
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

    ScanningData();
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

  private:
    std::unique_ptr<Data> impl_;
  };
  
  /****************************************************************************
   * The Detector interface
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

    DetectorData();
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

  private:
    std::unique_ptr<Data> impl_;
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

    const Data& internal() const;
    Data& internal();

    SystemConfiguration get_system_conf() const;
    void set_system_conf(const SystemConfiguration& system_conf);

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

    AtomData<T> get_atoms() const;
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

    std::vector<complex<T>> get_iw_psi() const;
    void set_iw_psi(const std::vector<complex<T>>& iw_psi);

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

    ScanningData<T> get_scanning() const;
    void set_scanning(const ScanningData<T>& scanning);

    DetectorData<T> get_detector() const;
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
    typedef std::vector<complex<T>> complex_vector_type;

    Output();
    Output(Output&);
    Output(Output&&);
    Output(Data&);
    Output& operator=(Output&);
    Output& operator=(Output&&);
    ~Output();

    const Data& internal() const;
    Data& internal();
		
    eTEM_Output_Type get_output_type() const;
		int get_ndetector() const;
    int get_nx() const;
    int get_ny() const;
    T get_dx() const;
    T get_dy() const;
		T get_dr() const;
    std::vector<T> get_x() const;
    std::vector<T> get_y() const;
    std::vector<T> get_r() const;
    std::vector<DetInt<vector_type>> get_image_tot() const;
    std::vector<DetInt<vector_type>> get_image_coh() const;
    std::vector<vector_type> get_m2psi_tot() const;
    std::vector<vector_type> get_m2psi_coh() const;
    std::vector<complex_vector_type> get_psi_coh() const;
    std::vector<vector_type> get_V() const;
    std::vector<complex_vector_type> get_trans() const;
    std::vector<complex_vector_type> get_psi_0() const;
		std::vector<bool> get_thk_gpu() const;

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
  template <typename T, eDevice dev>
  class MultisliceData {
  public:
    struct Data;

    typedef std::vector<T> vector_type;
    typedef std::vector<complex<T>> complex_vector_type;

    MultisliceData();
    ~MultisliceData();

    const Data& internal() const;

		void set_input_data(Input<T> &input, StreamIface<dev> &stream_i, FFTData<T, dev> &fft2_i);
		void operator()(Output<T> &output_multislice);
    
  private:
    std::unique_ptr<Data> impl_;

  };

  template <typename T>
  void test(const Input<T>&);

}  // namespace mt

#endif
