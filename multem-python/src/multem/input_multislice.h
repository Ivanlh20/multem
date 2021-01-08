/*
 *  input_multislice.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */
#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {
  
  /**
   * Define helper function for the mt::Input class
   */
  template <>
  template <typename T>
  struct Helpers <mt::Input<T>> {
   
    /**
     * Get the state
     */
    static py::tuple getstate(const mt::Input<T> &self) {
      return py::make_tuple(
        self.get_system_conf(),
        self.get_interaction_model(),
        self.get_potential_type(),
        self.get_pn_model(),
        self.get_pn_coh_contrib(),
        self.get_pn_single_conf(),
        self.get_pn_dim(),
        self.get_fp_dist(),
        self.get_pn_seed(),
        self.get_pn_nconf(),
        self.get_fp_iconf_0(),
        self.get_atoms(),
        self.get_is_crystal(),
        self.get_spec_rot_theta(),
        self.get_spec_rot_u0(),
        self.get_spec_rot_center_type(),
        self.get_spec_rot_center_p(),
        self.get_thick_type(),
        self.get_thick(),
        self.get_potential_slicing(),
        self.get_grid_2d(),
        self.get_output_area(),
        self.get_simulation_type(),
        self.get_iw_type(),
        self.get_iw_psi(),
        self.get_iw_x(),
        self.get_iw_y(),
        self.get_E_0(),
        self.get_lambda(),
        self.get_theta(),
        self.get_phi(),
        self.get_illumination_model(),
        self.get_temporal_spatial_incoh(),
        self.get_cond_lens(),
        self.get_obj_lens(),
        self.get_scanning(),
        self.get_detector(),
        self.get_eels_fr(),
        self.get_operation_mode(),
        self.get_slice_storage(),
        self.get_reverse_multislice(),
        self.get_mul_sign(),
        self.get_Vrl(),
        self.get_nR(),
        self.get_nrot(),
        self.get_cdl_var_type(),
        self.get_cdl_var(),
        self.get_iscan(),
        self.get_beam_x(),
        self.get_beam_y(),
        self.get_islice(),
        self.get_dp_Shift());
    }

    /**
     * Set the state
     */
    static mt::Input<T> setstate(py::tuple obj) {
      mt::Input<T> self;
      self.set_system_conf(obj[0].cast<mt::SystemConfiguration>());
      self.set_interaction_model(obj[1].cast<mt::eElec_Spec_Int_Model>());
      self.set_potential_type(obj[2].cast<mt::ePotential_Type>());
      self.set_pn_model(obj[3].cast<mt::ePhonon_Model>());
      self.set_pn_coh_contrib(obj[4].cast<bool>());
      self.set_pn_single_conf(obj[5].cast<bool>());
      self.set_pn_dim(obj[6].cast<mt::FP_Dim >());
      self.set_fp_dist(obj[7].cast<int>());
      self.set_pn_seed(obj[8].cast<int>());
      self.set_pn_nconf(obj[9].cast<int>());
      self.set_fp_iconf_0(obj[10].cast<int>());
      self.set_atoms(obj[11].cast<mt::AtomData<T>>());
      self.set_is_crystal(obj[12].cast<bool>());
      self.set_spec_rot_theta(obj[13].cast<double>());
      self.set_spec_rot_u0(obj[14].cast<mt::r3d<T>>());
      self.set_spec_rot_center_type(obj[15].cast<mt::eRot_Point_Type>());
      self.set_spec_rot_center_p(obj[16].cast<mt::r3d<T>>());
      self.set_thick_type(obj[17].cast<mt::eThick_Type>());
      self.set_thick(obj[18].cast<std::vector<T>>());
      self.set_potential_slicing(obj[19].cast<mt::ePotential_Slicing>());
      self.set_grid_2d(obj[20].cast<mt::Grid_2d<T>>());
      self.set_output_area(obj[21].cast<mt::Range_2d>());
      self.set_simulation_type(obj[22].cast<mt::eTEM_Sim_Type>());
      self.set_iw_type(obj[23].cast<mt::eIncident_Wave_Type>());
      self.set_iw_psi(obj[24].cast<std::vector<complex<T>>>());
      self.set_iw_x(obj[25].cast<std::vector<T>>());
      self.set_iw_y(obj[26].cast<std::vector<T>>());
      self.set_E_0(obj[27].cast<double>());
      self.set_lambda(obj[28].cast<double>());
      self.set_theta(obj[29].cast<double>());
      self.set_phi(obj[30].cast<double>());
      self.set_illumination_model(obj[31].cast<mt::eIllumination_Model>());
      self.set_temporal_spatial_incoh(obj[32].cast<mt::eTemporal_Spatial_Incoh>());
      self.set_cond_lens(obj[33].cast<mt::Lens<T>>());
      self.set_obj_lens(obj[34].cast<mt::Lens<T>>());
      self.set_scanning(obj[35].cast<mt::ScanningData<T>>());
      self.set_detector(obj[36].cast<mt::DetectorData<T>>());
      self.set_eels_fr(obj[37].cast<mt::EELS<T>>());
      self.set_operation_mode(obj[38].cast<mt::eOperation_Mode>());
      self.set_slice_storage(obj[39].cast<bool>());
      self.set_reverse_multislice(obj[40].cast<bool>());
      self.set_mul_sign(obj[41].cast<int>());
      self.set_Vrl(obj[42].cast<double>());
      self.set_nR(obj[43].cast<int>());
      self.set_nrot(obj[44].cast<int>());
      self.set_cdl_var_type(obj[45].cast<mt::eLens_Var_Type>());
      self.set_cdl_var(obj[46].cast<std::vector<T>>());
      self.set_iscan(obj[47].cast<std::vector<int>>());
      self.set_beam_x(obj[48].cast<std::vector<T>>());
      self.set_beam_y(obj[49].cast<std::vector<T>>());
      self.set_islice(obj[50].cast<int>());
      self.set_dp_Shift(obj[51].cast<bool>());
      return self;
    }
    
    static
    mt::eIncident_Wave_Type get_iw_type(const mt::Input<T>& self) {
      return self.get_iw_type();
    }
    
    static
    void set_iw_type(mt::Input<T>& self, mt::eIncident_Wave_Type iw_type) {
      self.set_incident_wave_type(iw_type);
    }
    
    static
    T get_theta(const mt::Input<T>& self) {
      return self.get_theta() / mt::c_deg_2_rad;
    }
    
    static
    void set_theta(mt::Input<T>& self, T theta) {
      self.set_theta(theta * mt::c_deg_2_rad);
    }
    
    static
    T get_phi(const mt::Input<T>& self) {
      return self.get_phi() / mt::c_deg_2_rad;
    }
    
    static
    void set_phi(mt::Input<T>& self, T phi) {
      self.set_phi(phi * mt::c_deg_2_rad);
    }
    
    static
    T get_spec_rot_theta(const mt::Input<T>& self) {
      return self.get_spec_rot_theta() / mt::c_deg_2_rad;
    }
    
    static
    void set_spec_rot_theta(mt::Input<T>& self, T spec_rot_theta) {
      self.set_spec_rot_theta(spec_rot_theta * mt::c_deg_2_rad);
    }
    
    static
    mt::r3d<T> get_spec_rot_u0(const mt::Input<T>& self) {
      return self.get_spec_rot_u0();
    }
    
    static
    void set_spec_rot_u0(mt::Input<T>& self, mt::r3d<T> spec_rot_u0) {
      spec_rot_u0.normalized();
      self.set_spec_rot_u0(spec_rot_u0);
    }

    static
    std::vector<mt::Atom<T>> get_spec_atoms(const mt::Input<T>& self) {
      return self.get_atoms().get_spec_atoms();
    }
    
    static
    void set_spec_atoms(mt::Input<T>& self, const std::vector<mt::Atom<T>>& spec_atoms) {
      self.get_atoms().set_spec_atoms(spec_atoms);
    }

    static
    T get_spec_dz(const mt::Input<T>& self) {
      return self.get_atoms().get_dz();
    }

    static
    void set_spec_dz(mt::Input<T>& self, T spec_dz) {
      self.get_atoms().set_dz(spec_dz);
      self.get_grid_2d().dz = spec_dz;
    }

    static
    T get_spec_lx(const mt::Input<T>& self) {
      return self.get_atoms().get_l_x();
    }

    static
    void set_spec_lx(mt::Input<T>& self, T spec_lx) {
      self.get_atoms().set_l_x(spec_lx);
      self.get_grid_2d().lx = spec_lx;
    }

    static
    T get_spec_ly(const mt::Input<T>& self) {
      return self.get_atoms().get_l_y();
    }

    static
    void set_spec_ly(mt::Input<T>& self, T spec_ly) {
      self.get_atoms().set_l_y(spec_ly);
      self.get_grid_2d().ly = spec_ly;
    }

    static
    T get_spec_lz(const mt::Input<T>& self) {
      return self.get_atoms().get_l_z();
    }

    static
    void set_spec_lz(mt::Input<T>& self, T spec_lz) {
      self.get_atoms().set_l_z(spec_lz);
    }

    static
    int get_spec_cryst_na(const mt::Input<T>& self) {
      return self.get_atoms().get_ct_na();
    }

    static
    void set_spec_cryst_na(mt::Input<T>& self, int spec_cryst_na) {
      self.get_atoms().set_ct_na(std::max(1, spec_cryst_na));
    }

    static
    int get_spec_cryst_nb(const mt::Input<T>& self) {
      return self.get_atoms().get_ct_nb();
    }

    static
    void set_spec_cryst_nb(mt::Input<T>& self, int spec_cryst_nb) {
      self.get_atoms().set_ct_nb(std::max(1, spec_cryst_nb));
    }

    static
    int get_spec_cryst_nc(const mt::Input<T>& self) {
      return self.get_atoms().get_ct_nc();
    }

    static
    void set_spec_cryst_nc(mt::Input<T>& self, int spec_cryst_nc) {
      self.get_atoms().set_ct_nc(std::max(1, spec_cryst_nc));
    }

    static
    T get_spec_cryst_a(const mt::Input<T>& self) {
      return self.get_atoms().get_ct_a();
    }

    static
    void set_spec_cryst_a(mt::Input<T>& self, T spec_cryst_a) {
      self.get_atoms().set_ct_a(std::max(T(0), spec_cryst_a));
    }

    static
    T get_spec_cryst_b(const mt::Input<T>& self) {
      return self.get_atoms().get_ct_b();
    }

    static
    void set_spec_cryst_b(mt::Input<T>& self, T spec_cryst_b) {
      self.get_atoms().set_ct_b(std::max(T(0), spec_cryst_b));
    }

    static
    T get_spec_cryst_c(const mt::Input<T>& self) {
      return self.get_atoms().get_ct_c();
    }

    static
    void set_spec_cryst_c(mt::Input<T>& self, T spec_cryst_c) {
      self.get_atoms().set_ct_c(std::max(T(0), spec_cryst_c));
    }

    static
    T get_spec_cryst_x0(const mt::Input<T>& self) {
      return self.get_atoms().get_ct_x0();
    }

    static
    void set_spec_cryst_x0(mt::Input<T>& self, T spec_cryst_x0) {
      self.get_atoms().set_ct_x0(std::max(T(0), spec_cryst_x0));
    }

    static
    T get_spec_cryst_y0(const mt::Input<T>& self) {
      return self.get_atoms().get_ct_y0();
    }

    static
    void set_spec_cryst_y0(mt::Input<T>& self, T spec_cryst_y0) {
      self.get_atoms().set_ct_y0(std::max(T(0), spec_cryst_y0));
    }

    static
    std::vector<mt::Amorp_Lay_Info<T>> get_spec_amorp(const mt::Input<T>& self) {
      return self.get_atoms().get_amorphous_parameters();
    }

    static
    void set_spec_amorp(mt::Input<T>& self, const std::vector<mt::Amorp_Lay_Info<T>>& spec_amorp) {
      self.get_atoms().set_amorphous_parameters(spec_amorp);
    }
    
    static
    mt::FP_Dim& get_pn_dim(const mt::Input<T>& self) {
      return self.get_pn_dim();
    }

    static
    void set_pn_dim(mt::Input<T>& self, int pn_dim) {
      self.get_pn_dim().set(pn_dim);
    }

    static
    int get_nx(const mt::Input<T>& self) {
      return self.get_grid_2d().nx;
    }

    static
    void set_nx(mt::Input<T>& self, int nx) {
      self.get_grid_2d().nx = nx;
    }

    static
    int get_ny(const mt::Input<T>& self) {
      return self.get_grid_2d().ny;
    }

    static
    void set_ny(mt::Input<T>& self, int ny) {
      self.get_grid_2d().ny = ny;
    }

    static
    bool get_bwl(const mt::Input<T>& self) {
      return self.get_grid_2d().bwl;
    }

    static
    void set_bwl(mt::Input<T>& self, bool bwl) {
      self.get_grid_2d().bwl = bwl;
    }

    static
    int get_cond_lens_m(const mt::Input<T>& self) {
      return self.get_cond_lens().m;
    }

    static
    void set_cond_lens_m(mt::Input<T>& self, int m) {
      self.get_cond_lens().m = m;
    }

    static
    T get_cond_lens_c_10(const mt::Input<T>& self) {
      return self.get_cond_lens().c_10;
    }

    static
    void set_cond_lens_c_10(mt::Input<T>& self, T c_10) {
      self.get_cond_lens().c_10 = c_10;
    }

    static
    T get_cond_lens_c_12(const mt::Input<T>& self) {
      return self.get_cond_lens().c_12;
    }

    static
    void set_cond_lens_c_12(mt::Input<T>& self, T c_12) {
      self.get_cond_lens().c_12 = c_12;
    }

    static
    T get_cond_lens_phi_12(const mt::Input<T>& self) {
      return self.get_cond_lens().phi_12 / mt::c_deg_2_rad;
    }

    static
    void set_cond_lens_phi_12(mt::Input<T>& self, T phi_12) {
      self.get_cond_lens().phi_12 = phi_12 * mt::c_deg_2_rad;
    }

    static
    T get_cond_lens_c_21(const mt::Input<T>& self) {
      return self.get_cond_lens().c_21;
    }

    static
    void set_cond_lens_c_21(mt::Input<T>& self, T c_21) {
      self.get_cond_lens().c_21 = c_21;
    }

    static
    T get_cond_lens_phi_21(const mt::Input<T>& self) {
      return self.get_cond_lens().phi_21 / mt::c_deg_2_rad;
    }

    static
    void set_cond_lens_phi_21(mt::Input<T>& self, T phi_21) {
      self.get_cond_lens().phi_21 = phi_21 * mt::c_deg_2_rad;
    }

    static
    T get_cond_lens_c_23(const mt::Input<T>& self) {
      return self.get_cond_lens().c_23;
    }

    static
    void set_cond_lens_c_23(mt::Input<T>& self, T c_23) {
      self.get_cond_lens().c_23 = c_23;
    }

    static
    T get_cond_lens_phi_23(const mt::Input<T>& self) {
      return self.get_cond_lens().phi_23 / mt::c_deg_2_rad;
    }

    static
    void set_cond_lens_phi_23(mt::Input<T>& self, T phi_23) {
      self.get_cond_lens().phi_23 = phi_23 * mt::c_deg_2_rad;
    }

    static
    T get_cond_lens_c_30(const mt::Input<T>& self) {
      return self.get_cond_lens().c_30 / mt::c_mm_2_Angs;
    }

    static
    void set_cond_lens_c_30(mt::Input<T>& self, T c_30) {
      self.get_cond_lens().c_30 = c_30 * mt::c_mm_2_Angs;
    }

    static
    T get_cond_lens_c_32(const mt::Input<T>& self) {
      return self.get_cond_lens().c_32;
    }

    static
    void set_cond_lens_c_32(mt::Input<T>& self, T c_32) {
      self.get_cond_lens().c_32 = c_32;
    }

    static
    T get_cond_lens_phi_32(const mt::Input<T>& self) {
      return self.get_cond_lens().phi_32 / mt::c_deg_2_rad;
    }

    static
    void set_cond_lens_phi_32(mt::Input<T>& self, T phi_32) {
      self.get_cond_lens().phi_32 = phi_32 * mt::c_deg_2_rad;
    }

    static
    T get_cond_lens_c_34(const mt::Input<T>& self) {
      return self.get_cond_lens().c_34;
    }

    static
    void set_cond_lens_c_34(mt::Input<T>& self, T c_34) {
      self.get_cond_lens().c_34 = c_34;
    }

    static
    T get_cond_lens_phi_34(const mt::Input<T>& self) {
      return self.get_cond_lens().phi_34 / mt::c_deg_2_rad;
    }

    static
    void set_cond_lens_phi_34(mt::Input<T>& self, T phi_34) {
      self.get_cond_lens().phi_34 = phi_34 * mt::c_deg_2_rad;
    }

    static
    T get_cond_lens_c_41(const mt::Input<T>& self) {
      return self.get_cond_lens().c_41;
    }

    static
    void set_cond_lens_c_41(mt::Input<T>& self, T c_41) {
      self.get_cond_lens().c_41 = c_41;
    }

    static
    T get_cond_lens_phi_41(const mt::Input<T>& self) {
      return self.get_cond_lens().phi_41 / mt::c_deg_2_rad;
    }

    static
    void set_cond_lens_phi_41(mt::Input<T>& self, T phi_41) {
      self.get_cond_lens().phi_41 = phi_41 * mt::c_deg_2_rad;
    }

    static
    T get_cond_lens_c_43(const mt::Input<T>& self) {
      return self.get_cond_lens().c_43;
    }

    static
    void set_cond_lens_c_43(mt::Input<T>& self, T c_43) {
      self.get_cond_lens().c_43 = c_43;
    }

    static
    T get_cond_lens_phi_43(const mt::Input<T>& self) {
      return self.get_cond_lens().phi_43 / mt::c_deg_2_rad;
    }

    static
    void set_cond_lens_phi_43(mt::Input<T>& self, T phi_43) {
      self.get_cond_lens().phi_43 = phi_43 * mt::c_deg_2_rad;
    }

    static
    T get_cond_lens_c_45(const mt::Input<T>& self) {
      return self.get_cond_lens().c_45;
    }

    static
    void set_cond_lens_c_45(mt::Input<T>& self, T c_45) {
      self.get_cond_lens().c_45 = c_45;
    }

    static
    T get_cond_lens_phi_45(const mt::Input<T>& self) {
      return self.get_cond_lens().phi_45 / mt::c_deg_2_rad;
    }

    static
    void set_cond_lens_phi_45(mt::Input<T>& self, T phi_45) {
      self.get_cond_lens().phi_45 = phi_45 * mt::c_deg_2_rad;
    }

    static
    T get_cond_lens_c_50(const mt::Input<T>& self) {
      return self.get_cond_lens().c_50 / mt::c_mm_2_Angs;
    }

    static
    void set_cond_lens_c_50(mt::Input<T>& self, T c_50) {
      self.get_cond_lens().c_50 = c_50 * mt::c_mm_2_Angs;
    }

    static
    T get_cond_lens_c_52(const mt::Input<T>& self) {
      return self.get_cond_lens().c_52;
    }

    static
    void set_cond_lens_c_52(mt::Input<T>& self, T c_52) {
      self.get_cond_lens().c_52 = c_52;
    }

    static
    T get_cond_lens_phi_52(const mt::Input<T>& self) {
      return self.get_cond_lens().phi_52 / mt::c_deg_2_rad;
    }

    static
    void set_cond_lens_phi_52(mt::Input<T>& self, T phi_52) {
      self.get_cond_lens().phi_52 = phi_52 * mt::c_deg_2_rad;
    }

    static
    T get_cond_lens_c_54(const mt::Input<T>& self) {
      return self.get_cond_lens().c_54;
    }

    static
    void set_cond_lens_c_54(mt::Input<T>& self, T c_54) {
      self.get_cond_lens().c_54 = c_54;
    }

    static
    T get_cond_lens_phi_54(const mt::Input<T>& self) {
      return self.get_cond_lens().phi_54 / mt::c_deg_2_rad;
    }

    static
    void set_cond_lens_phi_54(mt::Input<T>& self, T phi_54) {
      self.get_cond_lens().phi_54 = phi_54 * mt::c_deg_2_rad;
    }

    static
    T get_cond_lens_c_56(const mt::Input<T>& self) {
      return self.get_cond_lens().c_56;
    }

    static
    void set_cond_lens_c_56(mt::Input<T>& self, T c_56) {
      self.get_cond_lens().c_56 = c_56;
    }

    static
    T get_cond_lens_phi_56(const mt::Input<T>& self) {
      return self.get_cond_lens().phi_56 / mt::c_deg_2_rad;
    }

    static
    void set_cond_lens_phi_56(mt::Input<T>& self, T phi_56) {
      self.get_cond_lens().phi_56 = phi_56 * mt::c_deg_2_rad;
    }

    static
    T get_cond_lens_inner_aper_ang(const mt::Input<T>& self) {
      return self.get_cond_lens().inner_aper_ang / mt::c_mrad_2_rad;
    }

    static
    void set_cond_lens_inner_aper_ang(mt::Input<T>& self, T inner_aper_ang) {
      self.get_cond_lens().inner_aper_ang = inner_aper_ang * mt::c_mrad_2_rad;
    }

    static
    T get_cond_lens_outer_aper_ang(const mt::Input<T>& self) {
      return self.get_cond_lens().outer_aper_ang / mt::c_mrad_2_rad;
    }

    static
    void set_cond_lens_outer_aper_ang(mt::Input<T>& self, T outer_aper_ang) {
      self.get_cond_lens().outer_aper_ang = outer_aper_ang * mt::c_mrad_2_rad;
    }

    static
    T get_cond_lens_ti_a(const mt::Input<T>& self) {
      return self.get_cond_lens().ti_a;
    }

    static
    void set_cond_lens_ti_a(mt::Input<T>& self, T ti_a) {
      self.get_cond_lens().ti_a = ti_a;
    }

    static
    T get_cond_lens_ti_sigma(const mt::Input<T>& self) {
      return self.get_cond_lens().ti_sigma;
    }

    static
    void set_cond_lens_ti_sigma(mt::Input<T>& self, T ti_sigma) {
      self.get_cond_lens().ti_sigma = ti_sigma;
    }

    static
    T get_cond_lens_ti_beta(const mt::Input<T>& self) {
      return self.get_cond_lens().ti_beta;
    }

    static
    void set_cond_lens_ti_beta(mt::Input<T>& self, T ti_beta) {
      self.get_cond_lens().ti_beta = ti_beta;
    }

    static
    int get_cond_lens_ti_npts(const mt::Input<T>& self) {
      return self.get_cond_lens().ti_npts;
    }

    static
    void set_cond_lens_ti_npts(mt::Input<T>& self, int ti_npts) {
      self.get_cond_lens().ti_npts = ti_npts;
    }
    
    static
    T get_cond_lens_si_a(const mt::Input<T>& self) {
      return self.get_cond_lens().si_a;
    }

    static
    void set_cond_lens_si_a(mt::Input<T>& self, T si_a) {
      self.get_cond_lens().si_a = si_a;
      self.get_obj_lens().si_a = si_a;
    }

    static
    T get_cond_lens_si_sigma(const mt::Input<T>& self) {
      return self.get_cond_lens().si_sigma;
    }

    static
    void set_cond_lens_si_sigma(mt::Input<T>& self, T si_sigma) {
      self.get_cond_lens().si_sigma = si_sigma;
      self.get_obj_lens().si_sigma = si_sigma;
    }

    static
    T get_cond_lens_si_beta(const mt::Input<T>& self) {
      return self.get_cond_lens().si_beta;
    }
    
    static
    void set_cond_lens_si_beta(mt::Input<T>& self, T si_beta) {
      self.get_cond_lens().si_beta = si_beta;
      self.get_obj_lens().si_beta = si_beta;
    }
    
    static
    int get_cond_lens_si_rad_npts(const mt::Input<T>& self) {
      return self.get_cond_lens().si_rad_npts;
    }

    static
    void set_cond_lens_si_rad_npts(mt::Input<T>& self, int si_rad_npts) {
      self.get_cond_lens().si_rad_npts = si_rad_npts;
      self.get_obj_lens().si_rad_npts = si_rad_npts;
    }
    
    static
    int get_cond_lens_si_azm_npts(const mt::Input<T>& self) {
      return self.get_cond_lens().si_azm_npts;
    }

    static
    void set_cond_lens_si_azm_npts(mt::Input<T>& self, int si_azm_npts) {
      self.get_cond_lens().si_azm_npts = si_azm_npts;
      self.get_obj_lens().si_azm_npts = si_azm_npts;
    }

    static
    mt::eZero_Defocus_Type get_cond_lens_zero_defocus_type(const mt::Input<T>& self) {
      return self.get_cond_lens().zero_defocus_type;
    }

    static
    void set_cond_lens_zero_defocus_type(mt::Input<T>& self, mt::eZero_Defocus_Type zero_defocus_type) {
      self.get_cond_lens().zero_defocus_type = zero_defocus_type;
    }

    static
    T get_cond_lens_zero_defocus_plane(const mt::Input<T>& self) {
      return self.get_cond_lens().zero_defocus_plane;
    }

    static
    void set_cond_lens_zero_defocus_plane(mt::Input<T>& self, T zero_defocus_plane) {
      self.get_cond_lens().zero_defocus_plane = zero_defocus_plane;
    }

    static
    int get_obj_lens_m(const mt::Input<T>& self) {
      return self.get_obj_lens().m;
    }

    static
    void set_obj_lens_m(mt::Input<T>& self, int m) {
      self.get_obj_lens().m = m;
    }

    static
    T get_obj_lens_c_10(const mt::Input<T>& self) {
      return self.get_obj_lens().c_10;
    }

    static
    void set_obj_lens_c_10(mt::Input<T>& self, T c_10) {
      self.get_obj_lens().c_10 = c_10;
    }

    static
    T get_obj_lens_c_12(const mt::Input<T>& self) {
      return self.get_obj_lens().c_12;
    }

    static
    void set_obj_lens_c_12(mt::Input<T>& self, T c_12) {
      self.get_obj_lens().c_12 = c_12;
    }

    static
    T get_obj_lens_phi_12(const mt::Input<T>& self) {
      return self.get_obj_lens().phi_12 / mt::c_deg_2_rad;
    }

    static
    void set_obj_lens_phi_12(mt::Input<T>& self, T phi_12) {
      self.get_obj_lens().phi_12 = phi_12 * mt::c_deg_2_rad;
    }

    static
    T get_obj_lens_c_21(const mt::Input<T>& self) {
      return self.get_obj_lens().c_21;
    }

    static
    void set_obj_lens_c_21(mt::Input<T>& self, T c_21) {
      self.get_obj_lens().c_21 = c_21;
    }

    static
    T get_obj_lens_phi_21(const mt::Input<T>& self) {
      return self.get_obj_lens().phi_21 / mt::c_deg_2_rad;
    }

    static
    void set_obj_lens_phi_21(mt::Input<T>& self, T phi_21) {
      self.get_obj_lens().phi_21 = phi_21 * mt::c_deg_2_rad;
    }

    static
    T get_obj_lens_c_23(const mt::Input<T>& self) {
      return self.get_obj_lens().c_23;
    }

    static
    void set_obj_lens_c_23(mt::Input<T>& self, T c_23) {
      self.get_obj_lens().c_23 = c_23;
    }

    static
    T get_obj_lens_phi_23(const mt::Input<T>& self) {
      return self.get_obj_lens().phi_23 / mt::c_deg_2_rad;
    }

    static
    void set_obj_lens_phi_23(mt::Input<T>& self, T phi_23) {
      self.get_obj_lens().phi_23 = phi_23 * mt::c_deg_2_rad;
    }

    static
    T get_obj_lens_c_30(const mt::Input<T>& self) {
      return self.get_obj_lens().c_30 / mt::c_mm_2_Angs;
    }

    static
    void set_obj_lens_c_30(mt::Input<T>& self, T c_30) {
      self.get_obj_lens().c_30 = c_30 * mt::c_mm_2_Angs;
    }

    static
    T get_obj_lens_c_32(const mt::Input<T>& self) {
      return self.get_obj_lens().c_32;
    }

    static
    void set_obj_lens_c_32(mt::Input<T>& self, T c_32) {
      self.get_obj_lens().c_32 = c_32;
    }

    static
    T get_obj_lens_phi_32(const mt::Input<T>& self) {
      return self.get_obj_lens().phi_32 / mt::c_deg_2_rad;
    }

    static
    void set_obj_lens_phi_32(mt::Input<T>& self, T phi_32) {
      self.get_obj_lens().phi_32 = phi_32 * mt::c_deg_2_rad;
    }

    static
    T get_obj_lens_c_34(const mt::Input<T>& self) {
      return self.get_obj_lens().c_34;
    }

    static
    void set_obj_lens_c_34(mt::Input<T>& self, T c_34) {
      self.get_obj_lens().c_34 = c_34;
    }

    static
    T get_obj_lens_phi_34(const mt::Input<T>& self) {
      return self.get_obj_lens().phi_34 / mt::c_deg_2_rad;
    }

    static
    void set_obj_lens_phi_34(mt::Input<T>& self, T phi_34) {
      self.get_obj_lens().phi_34 = phi_34 * mt::c_deg_2_rad;
    }

    static
    T get_obj_lens_c_41(const mt::Input<T>& self) {
      return self.get_obj_lens().c_41;
    }

    static
    void set_obj_lens_c_41(mt::Input<T>& self, T c_41) {
      self.get_obj_lens().c_41 = c_41;
    }

    static
    T get_obj_lens_phi_41(const mt::Input<T>& self) {
      return self.get_obj_lens().phi_41 / mt::c_deg_2_rad;
    }

    static
    void set_obj_lens_phi_41(mt::Input<T>& self, T phi_41) {
      self.get_obj_lens().phi_41 = phi_41 * mt::c_deg_2_rad;
    }

    static
    T get_obj_lens_c_43(const mt::Input<T>& self) {
      return self.get_obj_lens().c_43;
    }

    static
    void set_obj_lens_c_43(mt::Input<T>& self, T c_43) {
      self.get_obj_lens().c_43 = c_43;
    }

    static
    T get_obj_lens_phi_43(const mt::Input<T>& self) {
      return self.get_obj_lens().phi_43 / mt::c_deg_2_rad;
    }

    static
    void set_obj_lens_phi_43(mt::Input<T>& self, T phi_43) {
      self.get_obj_lens().phi_43 = phi_43 * mt::c_deg_2_rad;
    }

    static
    T get_obj_lens_c_45(const mt::Input<T>& self) {
      return self.get_obj_lens().c_45;
    }

    static
    void set_obj_lens_c_45(mt::Input<T>& self, T c_45) {
      self.get_obj_lens().c_45 = c_45;
    }

    static
    T get_obj_lens_phi_45(const mt::Input<T>& self) {
      return self.get_obj_lens().phi_45 / mt::c_deg_2_rad;
    }

    static
    void set_obj_lens_phi_45(mt::Input<T>& self, T phi_45) {
      self.get_obj_lens().phi_45 = phi_45 * mt::c_deg_2_rad;
    }

    static
    T get_obj_lens_c_50(const mt::Input<T>& self) {
      return self.get_obj_lens().c_50 / mt::c_mm_2_Angs;
    }

    static
    void set_obj_lens_c_50(mt::Input<T>& self, T c_50) {
      self.get_obj_lens().c_50 = c_50 * mt::c_mm_2_Angs;
    }

    static
    T get_obj_lens_c_52(const mt::Input<T>& self) {
      return self.get_obj_lens().c_52;
    }

    static
    void set_obj_lens_c_52(mt::Input<T>& self, T c_52) {
      self.get_obj_lens().c_52 = c_52;
    }

    static
    T get_obj_lens_phi_52(const mt::Input<T>& self) {
      return self.get_obj_lens().phi_52 / mt::c_deg_2_rad;
    }

    static
    void set_obj_lens_phi_52(mt::Input<T>& self, T phi_52) {
      self.get_obj_lens().phi_52 = phi_52 * mt::c_deg_2_rad;
    }

    static
    T get_obj_lens_c_54(const mt::Input<T>& self) {
      return self.get_obj_lens().c_54;
    }

    static
    void set_obj_lens_c_54(mt::Input<T>& self, T c_54) {
      self.get_obj_lens().c_54 = c_54;
    }

    static
    T get_obj_lens_phi_54(const mt::Input<T>& self) {
      return self.get_obj_lens().phi_54 / mt::c_deg_2_rad;
    }

    static
    void set_obj_lens_phi_54(mt::Input<T>& self, T phi_54) {
      self.get_obj_lens().phi_54 = phi_54 * mt::c_deg_2_rad;
    }

    static
    T get_obj_lens_c_56(const mt::Input<T>& self) {
      return self.get_obj_lens().c_56;
    }

    static
    void set_obj_lens_c_56(mt::Input<T>& self, T c_56) {
      self.get_obj_lens().c_56 = c_56;
    }

    static
    T get_obj_lens_phi_56(const mt::Input<T>& self) {
      return self.get_obj_lens().phi_56 / mt::c_deg_2_rad;
    }

    static
    void set_obj_lens_phi_56(mt::Input<T>& self, T phi_56) {
      self.get_obj_lens().phi_56 = phi_56 * mt::c_deg_2_rad;
    }

    static
    T get_obj_lens_inner_aper_ang(const mt::Input<T>& self) {
      return self.get_obj_lens().inner_aper_ang / mt::c_mrad_2_rad;
    }

    static
    void set_obj_lens_inner_aper_ang(mt::Input<T>& self, T inner_aper_ang) {
      self.get_obj_lens().inner_aper_ang = inner_aper_ang * mt::c_mrad_2_rad;
    }

    static
    T get_obj_lens_outer_aper_ang(const mt::Input<T>& self) {
      return self.get_obj_lens().outer_aper_ang / mt::c_mrad_2_rad;
    }

    static
    void set_obj_lens_outer_aper_ang(mt::Input<T>& self, T outer_aper_ang) {
      self.get_obj_lens().outer_aper_ang = outer_aper_ang * mt::c_mrad_2_rad;
    }

    static
    T get_obj_lens_ti_a(const mt::Input<T>& self) {
      return self.get_obj_lens().ti_a;
    }

    static
    void set_obj_lens_ti_a(mt::Input<T>& self, T ti_a) {
      self.get_obj_lens().ti_a = ti_a;
    }

    static
    T get_obj_lens_ti_sigma(const mt::Input<T>& self) {
      return self.get_obj_lens().ti_sigma;
    }

    static
    void set_obj_lens_ti_sigma(mt::Input<T>& self, T ti_sigma) {
      self.get_obj_lens().ti_sigma = ti_sigma;
    }

    static
    T get_obj_lens_ti_beta(const mt::Input<T>& self) {
      return self.get_obj_lens().ti_beta;
    }

    static
    void set_obj_lens_ti_beta(mt::Input<T>& self, T ti_beta) {
      self.get_obj_lens().ti_beta = ti_beta;
    }

    static
    int get_obj_lens_ti_npts(const mt::Input<T>& self) {
      return self.get_obj_lens().ti_npts;
    }

    static
    void set_obj_lens_ti_npts(mt::Input<T>& self, int ti_npts) {
      self.get_obj_lens().ti_npts = ti_npts;
    }
    
    static
    mt::eZero_Defocus_Type get_obj_lens_zero_defocus_type(const mt::Input<T>& self) {
      return self.get_obj_lens().zero_defocus_type;
    }

    static
    void set_obj_lens_zero_defocus_type(mt::Input<T>& self, mt::eZero_Defocus_Type zero_defocus_type) {
      self.get_obj_lens().zero_defocus_type = zero_defocus_type;
    }

    static
    T get_obj_lens_zero_defocus_plane(const mt::Input<T>& self) {
      return self.get_obj_lens().zero_defocus_plane;
    }

    static
    void set_obj_lens_zero_defocus_plane(mt::Input<T>& self, T zero_defocus_plane) {
      self.get_obj_lens().zero_defocus_plane = zero_defocus_plane;
    }

    static
    mt::eScanning_Type get_scanning_type(const mt::Input<T>& self) {
      return self.get_scanning().get_type();
    }

    static
    void set_scanning_type(mt::Input<T>& self, mt::eScanning_Type scanning_type) {
      self.get_scanning().set_type(scanning_type);
    }

    static
    bool get_scanning_periodic(const mt::Input<T>& self) {
      return self.get_scanning().get_pbc();
    }

    static
    void set_scanning_periodic(mt::Input<T>& self, bool scanning_periodic) {
      self.get_scanning().set_pbc(scanning_periodic);
    }

    static
    bool get_scanning_square(const mt::Input<T>& self) {
      return self.get_scanning().get_spxs();
    }

    static
    void set_scanning_square(mt::Input<T>& self, bool scanning_square) {
      self.get_scanning().set_spxs(scanning_square);
    }

    static
    int get_scanning_ns(const mt::Input<T>& self) {
      return self.get_scanning().get_ns();
    }

    static
    void set_scanning_ns(mt::Input<T>& self, int scanning_ns) {
      self.get_scanning().set_ns(scanning_ns);
    }

    static
    T get_scanning_x0(const mt::Input<T>& self) {
      return self.get_scanning().get_x0();
    }

    static
    void set_scanning_x0(mt::Input<T>& self, T scanning_x0) {
      self.get_scanning().set_x0(scanning_x0);
    }

    static
    T get_scanning_y0(const mt::Input<T>& self) {
      return self.get_scanning().get_y0();
    }

    static
    void set_scanning_y0(mt::Input<T>& self, T scanning_y0) {
      self.get_scanning().set_y0(scanning_y0);
    }

    static
    T get_scanning_xe(const mt::Input<T>& self) {
      return self.get_scanning().get_xe();
    }

    static
    void set_scanning_xe(mt::Input<T>& self, T scanning_xe) {
      self.get_scanning().set_xe(scanning_xe);
    }

    static
    T get_scanning_ye(const mt::Input<T>& self) {
      return self.get_scanning().get_ye();
    }

    static
    void set_scanning_ye(mt::Input<T>& self, T scanning_ye) {
      self.get_scanning().set_ye(scanning_ye);
    }

    static
    T get_ped_nrot(const mt::Input<T>& self) {
      return self.get_nrot();
    }

    // FIXME Could simplify and remove
    static
    void set_ped_nrot(mt::Input<T>& self, T ped_nrot) {
      self.set_nrot(ped_nrot);
    }

    // FIXME Could simplify and remove
    static
    T get_ped_theta(const mt::Input<T>& self) {
      return self.get_theta() / mt::c_deg_2_rad;
    }

    // FIXME Could simplify and remove
    static
    void set_ped_theta(mt::Input<T>& self, T ped_theta) {
      self.set_theta(ped_theta * mt::c_deg_2_rad);
    }

    // FIXME Could simplify and remove
    static
    int get_hci_nrot(const mt::Input<T>& self) {
      return self.get_nrot();
    }

    // FIXME Could simplify and remove
    static
    void set_hci_nrot(mt::Input<T>& self, T hci_nrot) {
      self.set_nrot(hci_nrot);
    }

    // FIXME Could simplify and remove
    static
    T get_hci_theta(const mt::Input<T>& self) {
      return self.get_theta() / mt::c_deg_2_rad;
    }

    // FIXME Could simplify and remove
    static
    void set_hci_theta(mt::Input<T>& self, T hci_theta) {
      self.set_theta(hci_theta * mt::c_deg_2_rad);
    }

    static
    int get_eels_Z(const mt::Input<T>& self) {
      return self.get_eels_fr().Z;
    }

    static
    void set_eels_Z(mt::Input<T>& self, int eels_Z) {
      self.get_eels_fr().Z = eels_Z;
    }

    static
    T get_eels_E_loss(const mt::Input<T>& self) {
      return self.get_eels_fr().E_loss / mt::c_eV_2_keV;
    }

    static
    void set_eels_E_loss(mt::Input<T>& self, T eels_E_loss) {
      self.get_eels_fr().E_loss = eels_E_loss * mt::c_eV_2_keV;
    }

    static
    T get_eels_collection_angle(const mt::Input<T>& self) {
      return self.get_eels_fr().collection_angle / mt::c_mrad_2_rad;
    }

    static
    void set_eels_collection_angle(mt::Input<T>& self, T eels_collection_angle) {
      self.get_eels_fr().collection_angle = eels_collection_angle * mt::c_mrad_2_rad;
    }

    static
    int get_eels_m_selection(const mt::Input<T>& self) {
      return self.get_eels_fr().m_selection;
    }

    static
    void set_eels_m_selection(mt::Input<T>& self, int eels_m_selection) {
      self.get_eels_fr().m_selection = eels_m_selection;
    }

    static
    mt::eChannelling_Type get_eels_channelling_type(const mt::Input<T>& self) {
      return self.get_eels_fr().channelling_type;
    }

    static
    void set_eels_channelling_type(mt::Input<T>& self, mt::eChannelling_Type eels_channelling_type) {
      self.get_eels_fr().channelling_type = eels_channelling_type;
    }

    static
    int get_eftem_Z(const mt::Input<T>& self) {
      return self.get_eels_fr().Z;
    }

    static
    void set_eftem_Z(mt::Input<T>& self, int eftem_Z) {
      self.get_eels_fr().Z = eftem_Z;
    }

    static
    T get_eftem_E_loss(const mt::Input<T>& self) {
      return self.get_eels_fr().E_loss / mt::c_eV_2_keV;
    }

    static
    void set_eftem_E_loss(mt::Input<T>& self, T eftem_E_loss) {
      self.get_eels_fr().E_loss = eftem_E_loss * mt::c_eV_2_keV;
    }

    static
    T get_eftem_collection_angle(const mt::Input<T>& self) {
      return self.get_eels_fr().collection_angle / mt::c_mrad_2_rad;
    }

    static
    void set_eftem_collection_angle(mt::Input<T>& self, T eftem_collection_angle) {
      self.get_eels_fr().collection_angle = eftem_collection_angle * mt::c_mrad_2_rad;
    }

    static
    int get_eftem_m_selection(const mt::Input<T>& self) {
      return self.get_eels_fr().m_selection;
    }

    static
    void set_eftem_m_selection(mt::Input<T>& self, int eftem_m_selection) {
      self.get_eels_fr().m_selection = eftem_m_selection;
    }

    static
    mt::eChannelling_Type get_eftem_channelling_type(const mt::Input<T>& self) {
      return self.get_eels_fr().channelling_type;
    }

    static
    void set_eftem_channelling_type(mt::Input<T>& self, mt::eChannelling_Type eftem_channelling_type) {
      self.get_eels_fr().channelling_type = eftem_channelling_type;
    }

    static
    int get_output_area_ix_0(const mt::Input<T>& self) {
      return self.get_output_area().ix_0 + 1;
    }

    static
    void set_output_area_ix_0(mt::Input<T>& self, int ix_0) {
      self.get_output_area().ix_0 = ix_0 - 1;
    }

    static
    int get_output_area_iy_0(const mt::Input<T>& self) {
      return self.get_output_area().iy_0 + 1;
    }

    static
    void set_output_area_iy_0(mt::Input<T>& self, int iy_0) {
      self.get_output_area().iy_0 = iy_0 - 1;
    }

    static
    int get_output_area_ix_e(const mt::Input<T>& self) {
      return self.get_output_area().ix_e + 1;
    }

    static
    void set_output_area_ix_e(mt::Input<T>& self, int ix_e) {
      self.get_output_area().ix_0 = ix_e - 1;
    }

    static
    int get_output_area_iy_e(const mt::Input<T>& self) {
      return self.get_output_area().iy_e + 1;
    }

    static
    void set_output_area_iy_e(mt::Input<T>& self, int iy_e) {
      self.get_output_area().iy_0 = iy_e - 1;
    }

  };

}}


template <typename Module>
void export_input_multislice(Module m)
{
  typedef double T;

  // Wrap the mt::Input class
  py::class_<mt::Input<T>>(m, "Input")
    .def(py::init<>())
    .def_property("system_conf", 
        &mt::Input<T>::get_system_conf, 
        &mt::Input<T>::set_system_conf) 
    .def_property("interaction_model", 
        &mt::Input<T>::get_interaction_model, 
        &mt::Input<T>::set_interaction_model) 
    .def_property("potential_type", 
        &mt::Input<T>::get_potential_type,  
        &mt::Input<T>::set_potential_type)  
    .def_property("pn_model", 
        &mt::Input<T>::get_pn_model,
        &mt::Input<T>::set_pn_model)
    .def_property("pn_coh_contrib", 
        &mt::Input<T>::get_pn_coh_contrib,
        &mt::Input<T>::set_pn_coh_contrib)
    .def_property("pn_single_conf", 
        &mt::Input<T>::get_pn_single_conf,
        &mt::Input<T>::set_pn_single_conf)
    .def_property("fp_dist", 
        &mt::Input<T>::get_fp_dist,
        &mt::Input<T>::set_fp_dist)
    .def_property("pn_seed", 
        &mt::Input<T>::get_pn_seed,
        &mt::Input<T>::set_pn_seed)
    .def_property("pn_nconf", 
        &mt::Input<T>::get_pn_nconf,
        &mt::Input<T>::set_pn_nconf)
    .def_property("fp_iconf_0", 
        &mt::Input<T>::get_fp_iconf_0,
        &mt::Input<T>::set_fp_iconf_0)
    .def_property("atoms",
        &mt::Input<T>::get_atoms,
        &mt::Input<T>::set_atoms)
    .def_property("is_crystal",
        &mt::Input<T>::get_is_crystal,
        &mt::Input<T>::set_is_crystal)
    .def_property("spec_rot_center_type",
        &mt::Input<T>::get_spec_rot_center_type,
        &mt::Input<T>::set_spec_rot_center_type)
    .def_property("spec_rot_center_p",
        &mt::Input<T>::get_spec_rot_center_p,
        &mt::Input<T>::set_spec_rot_center_p)
    .def_property("thick_type",
        &mt::Input<T>::get_thick_type,
        &mt::Input<T>::set_thick_type)
    .def_property("thick",
        &mt::Input<T>::get_thick,
        &mt::Input<T>::set_thick)
    .def_property("potential_slicing",
        &mt::Input<T>::get_potential_slicing,
        &mt::Input<T>::set_potential_slicing)
    .def_property("grid_2d",
        &mt::Input<T>::get_grid_2d,
        &mt::Input<T>::set_grid_2d)
    .def_property("output_area",
        &mt::Input<T>::get_output_area,
        &mt::Input<T>::set_output_area)
    .def_property("simulation_type",
        &mt::Input<T>::get_simulation_type,
        &mt::Input<T>::set_simulation_type)
    .def_property("iw_psi",
        &mt::Input<T>::get_iw_psi,
        &mt::Input<T>::set_iw_psi)
    .def_property("iw_x",
        &mt::Input<T>::get_iw_x,
        &mt::Input<T>::set_iw_x)
    .def_property("iw_y",
        &mt::Input<T>::get_iw_y,
        &mt::Input<T>::set_iw_y)
    .def_property("E_0",
        &mt::Input<T>::get_E_0,
        &mt::Input<T>::set_E_0)
    .def_property("lambda",
        &mt::Input<T>::get_lambda,
        &mt::Input<T>::set_lambda)
    .def_property("illumination_model",
        &mt::Input<T>::get_illumination_model,
        &mt::Input<T>::set_illumination_model)
    .def_property("temporal_spatial_incoh",
        &mt::Input<T>::get_temporal_spatial_incoh,
        &mt::Input<T>::set_temporal_spatial_incoh)
    .def_property("cond_lens",
        &mt::Input<T>::get_cond_lens,
        &mt::Input<T>::set_cond_lens)
    .def_property("obj_lens",
        &mt::Input<T>::get_obj_lens,
        &mt::Input<T>::set_obj_lens)
    .def_property("scanning",
        &mt::Input<T>::get_scanning,
        &mt::Input<T>::set_scanning)
    .def_property("detector",
        &mt::Input<T>::get_detector,
        &mt::Input<T>::set_detector)
    .def_property("eels_fr",
        &mt::Input<T>::get_eels_fr,
        &mt::Input<T>::set_eels_fr)
    .def_property("operation_mode",
        &mt::Input<T>::get_operation_mode,
        &mt::Input<T>::set_operation_mode)
    .def_property("slice_storage",
        &mt::Input<T>::get_slice_storage,
        &mt::Input<T>::set_slice_storage)
    .def_property("reverse_multislice",
        &mt::Input<T>::get_reverse_multislice,
        &mt::Input<T>::set_reverse_multislice)
    .def_property("mul_sign",
        &mt::Input<T>::get_mul_sign,
        &mt::Input<T>::set_mul_sign)
    .def_property("Vrl",
        &mt::Input<T>::get_Vrl,
        &mt::Input<T>::set_Vrl)
    .def_property("nR",
        &mt::Input<T>::get_nR,
        &mt::Input<T>::set_nR)
    .def_property("nrot",
        &mt::Input<T>::get_nrot,
        &mt::Input<T>::set_nrot)
    .def_property("cdl_var_type",
        &mt::Input<T>::get_cdl_var_type,
        &mt::Input<T>::set_cdl_var_type)
    .def_property("cdl_var",
        &mt::Input<T>::get_cdl_var,
        &mt::Input<T>::set_cdl_var)
    .def_property("iscan",
        &mt::Input<T>::get_iscan,
        &mt::Input<T>::set_iscan)
    .def_property("beam_x",
        &mt::Input<T>::get_beam_x,
        &mt::Input<T>::set_beam_x)
    .def_property("beam_y",
        &mt::Input<T>::get_beam_y,
        &mt::Input<T>::set_beam_y)
    .def_property("islice",
        &mt::Input<T>::get_islice,
        &mt::Input<T>::set_islice)
    .def_property("dp_Shift",
        &mt::Input<T>::get_dp_Shift,
        &mt::Input<T>::set_dp_Shift)
    .def("validate_parameters", &mt::Input<T>::validate_parameters)
    .def(py::pickle(
        &py::detail::Helpers<mt::Input<T>>::getstate,
        &py::detail::Helpers<mt::Input<T>>::setstate))
    
    // The matlab interface defines some of these differently
    .def_property("theta", 
        &py::detail::Helpers<mt::Input<T>>::get_theta,
        &py::detail::Helpers<mt::Input<T>>::set_theta)
    .def_property("phi", 
        &py::detail::Helpers<mt::Input<T>>::get_phi,
        &py::detail::Helpers<mt::Input<T>>::set_phi)
    .def_property("spec_rot_theta", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_rot_theta,
        &py::detail::Helpers<mt::Input<T>>::set_spec_rot_theta)
    .def_property("spec_rot_u0", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_rot_u0,
        &py::detail::Helpers<mt::Input<T>>::set_spec_rot_u0)
    .def_property("iw_type", 
        &py::detail::Helpers<mt::Input<T>>::get_iw_type,
        &py::detail::Helpers<mt::Input<T>>::set_iw_type)
    
    // Additional properties to be compatible with matlab interface
    .def_property("spec_atoms", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_atoms,
        &py::detail::Helpers<mt::Input<T>>::set_spec_atoms)
    .def_property("spec_dz", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_dz,
        &py::detail::Helpers<mt::Input<T>>::set_spec_dz)
    .def_property("spec_lx", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_lx,
        &py::detail::Helpers<mt::Input<T>>::set_spec_lx)
    .def_property("spec_ly", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_ly,
        &py::detail::Helpers<mt::Input<T>>::set_spec_ly)
    .def_property("spec_lz", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_lz,
        &py::detail::Helpers<mt::Input<T>>::set_spec_lz)
    .def_property("spec_cryst_na", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_cryst_na,
        &py::detail::Helpers<mt::Input<T>>::set_spec_cryst_na)
    .def_property("spec_cryst_nb", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_cryst_nb,
        &py::detail::Helpers<mt::Input<T>>::set_spec_cryst_nb)
    .def_property("spec_cryst_nc", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_cryst_nc,
        &py::detail::Helpers<mt::Input<T>>::set_spec_cryst_nc)
    .def_property("spec_cryst_a", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_cryst_a,
        &py::detail::Helpers<mt::Input<T>>::set_spec_cryst_a)
    .def_property("spec_cryst_b", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_cryst_b,
        &py::detail::Helpers<mt::Input<T>>::set_spec_cryst_b)
    .def_property("spec_cryst_c", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_cryst_c,
        &py::detail::Helpers<mt::Input<T>>::set_spec_cryst_c)
    .def_property("spec_cryst_x0", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_cryst_x0,
        &py::detail::Helpers<mt::Input<T>>::set_spec_cryst_x0)
    .def_property("spec_cryst_y0", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_cryst_y0,
        &py::detail::Helpers<mt::Input<T>>::set_spec_cryst_y0)
    .def_property("spec_amorp", 
        &py::detail::Helpers<mt::Input<T>>::get_spec_amorp,
        &py::detail::Helpers<mt::Input<T>>::set_spec_amorp)
    .def_property("pn_dim", 
        &py::detail::Helpers<mt::Input<T>>::get_pn_dim,
        &py::detail::Helpers<mt::Input<T>>::set_pn_dim)
    .def_property("nx", 
        &py::detail::Helpers<mt::Input<T>>::get_nx,
        &py::detail::Helpers<mt::Input<T>>::set_nx)
    .def_property("nx", 
        &py::detail::Helpers<mt::Input<T>>::get_nx,
        &py::detail::Helpers<mt::Input<T>>::set_nx)
    .def_property("ny", 
        &py::detail::Helpers<mt::Input<T>>::get_ny,
        &py::detail::Helpers<mt::Input<T>>::set_ny)
    .def_property("bwl", 
        &py::detail::Helpers<mt::Input<T>>::get_bwl,
        &py::detail::Helpers<mt::Input<T>>::set_bwl)
    .def_property("cond_lens_m", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_m,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_m)
    .def_property("cond_lens_c_10", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_c_10,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_c_10)
    .def_property("cond_lens_c_12", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_c_12,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_c_12)
    .def_property("cond_lens_phi_12", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_phi_12,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_phi_12)
    .def_property("cond_lens_c_21", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_c_21,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_c_21)
    .def_property("cond_lens_phi_21", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_phi_21,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_phi_21)
    .def_property("cond_lens_c_23", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_c_23,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_c_23)
    .def_property("cond_lens_phi_23", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_phi_23,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_phi_23)
    .def_property("cond_lens_c_30", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_c_30,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_c_30)
    .def_property("cond_lens_c_32", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_c_32,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_c_32)
    .def_property("cond_lens_phi_32", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_phi_32,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_phi_32)
    .def_property("cond_lens_c_34", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_c_34,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_c_34)
    .def_property("cond_lens_phi_34", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_phi_34,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_phi_34)
    .def_property("cond_lens_c_41", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_c_41,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_c_41)
    .def_property("cond_lens_phi_41", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_phi_41,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_phi_41)
    .def_property("cond_lens_c_43", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_c_43,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_c_43)
    .def_property("cond_lens_phi_43", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_phi_43,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_phi_43)
    .def_property("cond_lens_c_45", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_c_45,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_c_45)
    .def_property("cond_lens_phi_45", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_phi_45,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_phi_45)
    .def_property("cond_lens_c_50", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_c_50,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_c_50)
    .def_property("cond_lens_c_52", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_c_52,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_c_52)
    .def_property("cond_lens_phi_52", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_phi_52,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_phi_52)
    .def_property("cond_lens_c_54", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_c_54,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_c_54)
    .def_property("cond_lens_phi_54", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_phi_54,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_phi_54)
    .def_property("cond_lens_c_56", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_c_56,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_c_56)
    .def_property("cond_lens_phi_56", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_phi_56,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_phi_56)
    .def_property("cond_lens_inner_aper_ang", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_inner_aper_ang,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_inner_aper_ang)
    .def_property("cond_lens_outer_aper_ang", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_outer_aper_ang,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_outer_aper_ang)
    .def_property("cond_lens_ti_a", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_ti_a,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_ti_a)
    .def_property("cond_lens_ti_sigma", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_ti_sigma,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_ti_sigma)
    .def_property("cond_lens_ti_beta", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_ti_beta,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_ti_beta)
    .def_property("cond_lens_ti_npts", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_ti_npts,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_ti_npts)
    .def_property("cond_lens_si_a", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_si_a,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_si_a)
    .def_property("cond_lens_si_sigma", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_si_sigma,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_si_sigma)
    .def_property("cond_lens_si_beta", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_si_beta,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_si_beta)
    .def_property("cond_lens_si_rad_npts", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_si_rad_npts,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_si_rad_npts)
    .def_property("cond_lens_si_azm_npts", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_si_azm_npts,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_si_azm_npts)
    .def_property("cond_lens_zero_defocus_type", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_zero_defocus_type,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_zero_defocus_type)
    .def_property("cond_lens_zero_defocus_plane", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_zero_defocus_plane,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_zero_defocus_plane)
    .def_property("obj_lens_m", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_m,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_m)
    .def_property("obj_lens_c_10", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_c_10,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_c_10)
    .def_property("obj_lens_c_12", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_c_12,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_c_12)
    .def_property("obj_lens_phi_12", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_phi_12,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_phi_12)
    .def_property("obj_lens_c_21", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_c_21,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_c_21)
    .def_property("obj_lens_phi_21", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_phi_21,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_phi_21)
    .def_property("obj_lens_c_23", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_c_23,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_c_23)
    .def_property("obj_lens_phi_23", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_phi_23,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_phi_23)
    .def_property("obj_lens_c_30", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_c_30,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_c_30)
    .def_property("obj_lens_c_32", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_c_32,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_c_32)
    .def_property("obj_lens_phi_32", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_phi_32,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_phi_32)
    .def_property("obj_lens_c_34", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_c_34,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_c_34)
    .def_property("obj_lens_phi_34", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_phi_34,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_phi_34)
    .def_property("obj_lens_c_41", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_c_41,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_c_41)
    .def_property("obj_lens_phi_41", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_phi_41,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_phi_41)
    .def_property("obj_lens_c_43", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_c_43,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_c_43)
    .def_property("obj_lens_phi_43", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_phi_43,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_phi_43)
    .def_property("obj_lens_c_45", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_c_45,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_c_45)
    .def_property("obj_lens_phi_45", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_phi_45,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_phi_45)
    .def_property("obj_lens_c_50", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_c_50,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_c_50)
    .def_property("obj_lens_c_52", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_c_52,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_c_52)
    .def_property("obj_lens_phi_52", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_phi_52,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_phi_52)
    .def_property("obj_lens_c_54", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_c_54,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_c_54)
    .def_property("obj_lens_phi_54", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_phi_54,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_phi_54)
    .def_property("obj_lens_c_56", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_c_56,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_c_56)
    .def_property("obj_lens_phi_56", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_phi_56,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_phi_56)
    .def_property("obj_lens_inner_aper_ang", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_inner_aper_ang,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_inner_aper_ang)
    .def_property("obj_lens_outer_aper_ang", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_outer_aper_ang,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_outer_aper_ang)
    .def_property("cond_lens_ti_a", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_ti_a,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_ti_a)
    .def_property("cond_lens_ti_sigma", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_ti_sigma,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_ti_sigma)
    .def_property("cond_lens_ti_beta", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_ti_beta,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_ti_beta)
    .def_property("cond_lens_ti_npts", 
        &py::detail::Helpers<mt::Input<T>>::get_cond_lens_ti_npts,
        &py::detail::Helpers<mt::Input<T>>::set_cond_lens_ti_npts)
    .def_property("obj_lens_zero_defocus_type", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_zero_defocus_type,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_zero_defocus_type)
    .def_property("obj_lens_zero_defocus_plane", 
        &py::detail::Helpers<mt::Input<T>>::get_obj_lens_zero_defocus_plane,
        &py::detail::Helpers<mt::Input<T>>::set_obj_lens_zero_defocus_plane)
    .def_property("scanning_type", 
        &py::detail::Helpers<mt::Input<T>>::get_scanning_type,
        &py::detail::Helpers<mt::Input<T>>::set_scanning_type)
    .def_property("scanning_periodic", 
        &py::detail::Helpers<mt::Input<T>>::get_scanning_periodic,
        &py::detail::Helpers<mt::Input<T>>::set_scanning_periodic)
    .def_property("scanning_square", 
        &py::detail::Helpers<mt::Input<T>>::get_scanning_square,
        &py::detail::Helpers<mt::Input<T>>::set_scanning_square)
    .def_property("scanning_ns", 
        &py::detail::Helpers<mt::Input<T>>::get_scanning_ns,
        &py::detail::Helpers<mt::Input<T>>::set_scanning_ns)
    .def_property("scanning_x0", 
        &py::detail::Helpers<mt::Input<T>>::get_scanning_x0,
        &py::detail::Helpers<mt::Input<T>>::set_scanning_x0)
    .def_property("scanning_y0", 
        &py::detail::Helpers<mt::Input<T>>::get_scanning_y0,
        &py::detail::Helpers<mt::Input<T>>::set_scanning_y0)
    .def_property("scanning_xe", 
        &py::detail::Helpers<mt::Input<T>>::get_scanning_xe,
        &py::detail::Helpers<mt::Input<T>>::set_scanning_xe)
    .def_property("scanning_ye", 
        &py::detail::Helpers<mt::Input<T>>::get_scanning_ye,
        &py::detail::Helpers<mt::Input<T>>::set_scanning_ye)
    .def_property("ped_nrot", 
        &py::detail::Helpers<mt::Input<T>>::get_ped_nrot,
        &py::detail::Helpers<mt::Input<T>>::set_ped_nrot)
    .def_property("ped_theta", 
        &py::detail::Helpers<mt::Input<T>>::get_ped_theta,
        &py::detail::Helpers<mt::Input<T>>::set_ped_theta)
    .def_property("hci_nrot", 
        &py::detail::Helpers<mt::Input<T>>::get_hci_nrot,
        &py::detail::Helpers<mt::Input<T>>::set_hci_nrot)
    .def_property("hci_theta", 
        &py::detail::Helpers<mt::Input<T>>::get_hci_theta,
        &py::detail::Helpers<mt::Input<T>>::set_hci_theta)
    .def_property("eels_Z", 
        &py::detail::Helpers<mt::Input<T>>::get_eels_Z,
        &py::detail::Helpers<mt::Input<T>>::set_eels_Z)
    .def_property("eels_E_loss", 
        &py::detail::Helpers<mt::Input<T>>::get_eels_E_loss,
        &py::detail::Helpers<mt::Input<T>>::set_eels_E_loss)
    .def_property("eels_collection_angle", 
        &py::detail::Helpers<mt::Input<T>>::get_eels_collection_angle,
        &py::detail::Helpers<mt::Input<T>>::set_eels_collection_angle)
    .def_property("eels_m_selection", 
        &py::detail::Helpers<mt::Input<T>>::get_eels_m_selection,
        &py::detail::Helpers<mt::Input<T>>::set_eels_m_selection)
    .def_property("eels_channelling_type", 
        &py::detail::Helpers<mt::Input<T>>::get_eels_channelling_type,
        &py::detail::Helpers<mt::Input<T>>::set_eels_channelling_type)
    .def_property("eftem_Z", 
        &py::detail::Helpers<mt::Input<T>>::get_eftem_Z,
        &py::detail::Helpers<mt::Input<T>>::set_eftem_Z)
    .def_property("eftem_E_loss", 
        &py::detail::Helpers<mt::Input<T>>::get_eftem_E_loss,
        &py::detail::Helpers<mt::Input<T>>::set_eftem_E_loss)
    .def_property("eftem_collection_angle", 
        &py::detail::Helpers<mt::Input<T>>::get_eftem_collection_angle,
        &py::detail::Helpers<mt::Input<T>>::set_eftem_collection_angle)
    .def_property("eftem_m_selection", 
        &py::detail::Helpers<mt::Input<T>>::get_eftem_m_selection,
        &py::detail::Helpers<mt::Input<T>>::set_eftem_m_selection)
    .def_property("eftem_channelling_type", 
        &py::detail::Helpers<mt::Input<T>>::get_eftem_channelling_type,
        &py::detail::Helpers<mt::Input<T>>::set_eftem_channelling_type)
    .def_property("output_area_ix_0", 
        &py::detail::Helpers<mt::Input<T>>::get_output_area_ix_0,
        &py::detail::Helpers<mt::Input<T>>::set_output_area_ix_0)
    .def_property("output_area_iy_0", 
        &py::detail::Helpers<mt::Input<T>>::get_output_area_iy_0,
        &py::detail::Helpers<mt::Input<T>>::set_output_area_iy_0)
    .def_property("output_area_ix_e", 
        &py::detail::Helpers<mt::Input<T>>::get_output_area_ix_e,
        &py::detail::Helpers<mt::Input<T>>::set_output_area_ix_e)
    .def_property("output_area_iy_e", 
        &py::detail::Helpers<mt::Input<T>>::get_output_area_iy_e,
        &py::detail::Helpers<mt::Input<T>>::set_output_area_iy_e)
    ;
}



