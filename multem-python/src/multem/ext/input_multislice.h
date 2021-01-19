/*
 * This file is part of multem-python.
 * Copyright 2021 Diamond Light Source
 * Copyright 2021 Rosalind Franklin Institute
 *
 * Author: James Parkhurst
 *
 * multem-python is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * multem-python is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with multem-python. If not, see <http:// www.gnu.org/licenses/>.
 */
#ifndef MULTEM_PYTHON_INPUT_MULTISLICE_H
#define MULTEM_PYTHON_INPUT_MULTISLICE_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {
  
  template <typename T>
  mt::Input_Multislice<T> Input_Multislice_constructor_internal(const mt::System_Configuration &system_conf) {
    mt::Input_Multislice<T> result;
    result.system_conf = system_conf;
    return result;
  }

  inline
  py::object Input_Multislice_constructor(const mt::System_Configuration &system_conf) {
    return system_conf.precision == mt::eP_float
      ? py::cast(Input_Multislice_constructor_internal<float>(system_conf))
      : py::cast(Input_Multislice_constructor_internal<double>(system_conf));
  }

  /**
   * Define helper function for the mt::Input_Multislice class
   */
  template <>
  template <typename T>
  struct Helpers<mt::Input_Multislice<T>> {
    /**
     * Get the state
     */
    static py::tuple getstate(const mt::Input_Multislice<T>& self) {
      return py::make_tuple(self.system_conf,
                            self.interaction_model,
                            self.potential_type,
                            self.pn_model,
                            self.pn_coh_contrib,
                            self.pn_single_conf,
                            self.pn_dim,
                            self.fp_dist,
                            self.pn_seed,
                            self.pn_nconf,
                            self.fp_iconf_0,
                            self.atoms,
                            self.is_crystal,
                            self.spec_rot_theta,
                            self.spec_rot_u0,
                            self.spec_rot_center_type,
                            self.spec_rot_center_p,
                            self.thick_type,
                            self.thick,
                            self.potential_slicing,
                            self.grid_2d,
                            self.output_area,
                            self.simulation_type,
                            self.iw_type,
                            self.iw_psi,
                            self.iw_x,
                            self.iw_y,
                            self.E_0,
                            self.lambda,
                            self.theta,
                            self.phi,
                            self.illumination_model,
                            self.temporal_spatial_incoh,
                            self.cond_lens,
                            self.obj_lens,
                            self.scanning,
                            self.detector,
                            self.eels_fr,
                            self.operation_mode,
                            self.slice_storage,
                            self.reverse_multislice,
                            self.mul_sign,
                            self.Vrl,
                            self.nR,
                            self.nrot,
                            self.cdl_var_type,
                            self.cdl_var,
                            self.iscan,
                            self.beam_x,
                            self.beam_y,
                            self.islice,
                            self.dp_Shift);
    }

    /**
     * Set the state
     */
    static mt::Input_Multislice<T> setstate(py::tuple obj) {
      mt::Input_Multislice<T> self;
      self.system_conf = obj[0].cast<mt::System_Configuration>();
      self.interaction_model = obj[1].cast<mt::eElec_Spec_Int_Model>();
      self.potential_type = obj[2].cast<mt::ePotential_Type>();
      self.pn_model = obj[3].cast<mt::ePhonon_Model>();
      self.pn_coh_contrib = obj[4].cast<bool>();
      self.pn_single_conf = obj[5].cast<bool>();
      self.pn_dim = obj[6].cast<mt::FP_Dim>();
      self.fp_dist = obj[7].cast<int>();
      self.pn_seed = obj[8].cast<int>();
      self.pn_nconf = obj[9].cast<int>();
      self.fp_iconf_0 = obj[10].cast<int>();
      self.atoms = obj[11].cast<mt::Atom_Data<T>>();
      self.is_crystal = obj[12].cast<bool>();
      self.spec_rot_theta = obj[13].cast<double>();
      self.spec_rot_u0 = obj[14].cast<mt::r3d<T>>();
      self.spec_rot_center_type = obj[15].cast<mt::eRot_Point_Type>();
      self.spec_rot_center_p = obj[16].cast<mt::r3d<T>>();
      self.thick_type = obj[17].cast<mt::eThick_Type>();
      self.thick = obj[18].cast<std::vector<T>>();
      self.potential_slicing = obj[19].cast<mt::ePotential_Slicing>();
      self.grid_2d = obj[20].cast<mt::Grid_2d<T>>();
      self.output_area = obj[21].cast<mt::Range_2d>();
      self.simulation_type = obj[22].cast<mt::eTEM_Sim_Type>();
      self.iw_type = obj[23].cast<mt::eIncident_Wave_Type>();
      self.iw_psi = obj[24].cast<std::vector<complex<T>>>();
      self.iw_x = obj[25].cast<std::vector<T>>();
      self.iw_y = obj[26].cast<std::vector<T>>();
      self.E_0 = obj[27].cast<double>();
      self.lambda = obj[28].cast<double>();
      self.theta = obj[29].cast<double>();
      self.phi = obj[30].cast<double>();
      self.illumination_model = obj[31].cast<mt::eIllumination_Model>();
      self.temporal_spatial_incoh = obj[32].cast<mt::eTemporal_Spatial_Incoh>();
      self.cond_lens = obj[33].cast<mt::Lens<T>>();
      self.obj_lens = obj[34].cast<mt::Lens<T>>();
      self.scanning = obj[35].cast<mt::Scanning<T>>();
      self.detector = obj[36].cast<mt::Detector<T, mt::e_host>>();
      self.eels_fr = obj[37].cast<mt::EELS<T>>();
      self.operation_mode = obj[38].cast<mt::eOperation_Mode>();
      self.slice_storage = obj[39].cast<bool>();
      self.reverse_multislice = obj[40].cast<bool>();
      self.mul_sign = obj[41].cast<int>();
      self.Vrl = obj[42].cast<double>();
      self.nR = obj[43].cast<int>();
      self.nrot = obj[44].cast<int>();
      self.cdl_var_type = obj[45].cast<mt::eLens_Var_Type>();
      self.cdl_var = obj[46].cast<std::vector<T>>();
      self.iscan = obj[47].cast<std::vector<int>>();
      self.beam_x = obj[48].cast<std::vector<T>>();
      self.beam_y = obj[49].cast<std::vector<T>>();
      self.islice = obj[50].cast<int>();
      self.dp_Shift = obj[51].cast<bool>();
      return self;
    }

    static mt::eIncident_Wave_Type get_iw_type(const mt::Input_Multislice<T>& self) {
      return self.iw_type;
    }

    static void set_iw_type(mt::Input_Multislice<T>& self,
                            mt::eIncident_Wave_Type iw_type) {
      self.set_incident_wave_type(iw_type);
    }

    static T get_theta(const mt::Input_Multislice<T>& self) {
      return self.theta / mt::c_deg_2_rad;
    }

    static void set_theta(mt::Input_Multislice<T>& self, T theta) {
      self.theta = theta * mt::c_deg_2_rad;
    }

    static T get_phi(const mt::Input_Multislice<T>& self) {
      return self.phi / mt::c_deg_2_rad;
    }

    static void set_phi(mt::Input_Multislice<T>& self, T phi) {
      self.phi = phi * mt::c_deg_2_rad;
    }

    static T get_spec_rot_theta(const mt::Input_Multislice<T>& self) {
      return self.spec_rot_theta / mt::c_deg_2_rad;
    }

    static void set_spec_rot_theta(mt::Input_Multislice<T>& self, T spec_rot_theta) {
      self.spec_rot_theta = spec_rot_theta * mt::c_deg_2_rad;
    }

    static mt::r3d<T> get_spec_rot_u0(const mt::Input_Multislice<T>& self) {
      return self.spec_rot_u0;
    }

    static void set_spec_rot_u0(mt::Input_Multislice<T>& self, mt::r3d<T> spec_rot_u0) {
      spec_rot_u0.normalized();
      self.spec_rot_u0 = spec_rot_u0;
    }

    static std::vector<mt::Atom<T>> get_spec_atoms(
      const mt::Input_Multislice<T>& self) {
      return Helpers<mt::Atom_Data<T>>::get_spec_atoms(self.atoms);
    }
    
    static void set_spec_atoms(mt::Input_Multislice<T>& self,
                               const py::object &spec_atoms) {
        Helpers<mt::Atom_Data<T>>::set_spec_atoms(self.atoms, spec_atoms);
    }

    static T get_spec_dz(const mt::Input_Multislice<T>& self) {
      return self.atoms.dz;
    }

    static void set_spec_dz(mt::Input_Multislice<T>& self, T spec_dz) {
      self.atoms.dz = spec_dz;
      self.grid_2d.dz = spec_dz;
    }

    static T get_spec_lx(const mt::Input_Multislice<T>& self) {
      return self.atoms.l_x;
    }

    static void set_spec_lx(mt::Input_Multislice<T>& self, T spec_lx) {
      self.atoms.l_x = spec_lx;
      self.grid_2d.lx = spec_lx;
    }

    static T get_spec_ly(const mt::Input_Multislice<T>& self) {
      return self.atoms.l_y;
    }

    static void set_spec_ly(mt::Input_Multislice<T>& self, T spec_ly) {
      self.atoms.l_y = spec_ly;
      self.grid_2d.ly = spec_ly;
    }

    static T get_spec_lz(const mt::Input_Multislice<T>& self) {
      return self.atoms.l_z;
    }

    static void set_spec_lz(mt::Input_Multislice<T>& self, T spec_lz) {
      self.atoms.l_z = spec_lz;
    }

    static int get_spec_cryst_na(const mt::Input_Multislice<T>& self) {
      return self.atoms.ct_na;
    }

    static void set_spec_cryst_na(mt::Input_Multislice<T>& self, int spec_cryst_na) {
      self.atoms.ct_na = std::max(1, spec_cryst_na);
    }

    static int get_spec_cryst_nb(const mt::Input_Multislice<T>& self) {
      return self.atoms.ct_nb;
    }

    static void set_spec_cryst_nb(mt::Input_Multislice<T>& self, int spec_cryst_nb) {
      self.atoms.ct_nb = std::max(1, spec_cryst_nb);
    }

    static int get_spec_cryst_nc(const mt::Input_Multislice<T>& self) {
      return self.atoms.ct_nc;
    }

    static void set_spec_cryst_nc(mt::Input_Multislice<T>& self, int spec_cryst_nc) {
      self.atoms.ct_nc = std::max(1, spec_cryst_nc);
    }

    static T get_spec_cryst_a(const mt::Input_Multislice<T>& self) {
      return self.atoms.ct_a;
    }

    static void set_spec_cryst_a(mt::Input_Multislice<T>& self, T spec_cryst_a) {
      self.atoms.ct_a = std::max(T(0), spec_cryst_a);
    }

    static T get_spec_cryst_b(const mt::Input_Multislice<T>& self) {
      return self.atoms.ct_b;
    }

    static void set_spec_cryst_b(mt::Input_Multislice<T>& self, T spec_cryst_b) {
      self.atoms.ct_b = std::max(T(0), spec_cryst_b);
    }

    static T get_spec_cryst_c(const mt::Input_Multislice<T>& self) {
      return self.atoms.ct_c;
    }

    static void set_spec_cryst_c(mt::Input_Multislice<T>& self, T spec_cryst_c) {
      self.atoms.ct_c = std::max(T(0), spec_cryst_c);
    }

    static T get_spec_cryst_x0(const mt::Input_Multislice<T>& self) {
      return self.atoms.ct_x0;
    }

    static void set_spec_cryst_x0(mt::Input_Multislice<T>& self, T spec_cryst_x0) {
      self.atoms.ct_x0 = std::max(T(0), spec_cryst_x0);
    }

    static T get_spec_cryst_y0(const mt::Input_Multislice<T>& self) {
      return self.atoms.ct_y0;
    }

    static void set_spec_cryst_y0(mt::Input_Multislice<T>& self, T spec_cryst_y0) {
      self.atoms.ct_y0 = std::max(T(0), spec_cryst_y0);
    }

    static std::vector<mt::Amorp_Lay_Info<T>> get_spec_amorp(
      const mt::Input_Multislice<T>& self) {
      return self.atoms.amorp_lay_info;
    }

    static void set_spec_amorp(mt::Input_Multislice<T>& self,
                               const std::vector<mt::Amorp_Lay_Info<T>>& spec_amorp) {
      self.atoms.amorp_lay_info = spec_amorp;
    }

    static int get_nx(const mt::Input_Multislice<T>& self) {
      return self.grid_2d.nx;
    }

    static void set_nx(mt::Input_Multislice<T>& self, int nx) {
      self.grid_2d.nx = nx;
    }

    static int get_ny(const mt::Input_Multislice<T>& self) {
      return self.grid_2d.ny;
    }

    static void set_ny(mt::Input_Multislice<T>& self, int ny) {
      self.grid_2d.ny = ny;
    }

    static bool get_bwl(const mt::Input_Multislice<T>& self) {
      return self.grid_2d.bwl;
    }

    static void set_bwl(mt::Input_Multislice<T>& self, bool bwl) {
      self.grid_2d.bwl = bwl;
    }

    static int get_cond_lens_m(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.m;
    }

    static void set_cond_lens_m(mt::Input_Multislice<T>& self, int m) {
      self.cond_lens.m = m;
    }

    static T get_cond_lens_c_10(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.c_10;
    }

    static void set_cond_lens_c_10(mt::Input_Multislice<T>& self, T c_10) {
      self.cond_lens.c_10 = c_10;
    }

    static T get_cond_lens_c_12(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.c_12;
    }

    static void set_cond_lens_c_12(mt::Input_Multislice<T>& self, T c_12) {
      self.cond_lens.c_12 = c_12;
    }

    static T get_cond_lens_phi_12(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.phi_12 / mt::c_deg_2_rad;
    }

    static void set_cond_lens_phi_12(mt::Input_Multislice<T>& self, T phi_12) {
      self.cond_lens.phi_12 = phi_12 * mt::c_deg_2_rad;
    }

    static T get_cond_lens_c_21(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.c_21;
    }

    static void set_cond_lens_c_21(mt::Input_Multislice<T>& self, T c_21) {
      self.cond_lens.c_21 = c_21;
    }

    static T get_cond_lens_phi_21(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.phi_21 / mt::c_deg_2_rad;
    }

    static void set_cond_lens_phi_21(mt::Input_Multislice<T>& self, T phi_21) {
      self.cond_lens.phi_21 = phi_21 * mt::c_deg_2_rad;
    }

    static T get_cond_lens_c_23(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.c_23;
    }

    static void set_cond_lens_c_23(mt::Input_Multislice<T>& self, T c_23) {
      self.cond_lens.c_23 = c_23;
    }

    static T get_cond_lens_phi_23(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.phi_23 / mt::c_deg_2_rad;
    }

    static void set_cond_lens_phi_23(mt::Input_Multislice<T>& self, T phi_23) {
      self.cond_lens.phi_23 = phi_23 * mt::c_deg_2_rad;
    }

    static T get_cond_lens_c_30(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.c_30 / mt::c_mm_2_Angs;
    }

    static void set_cond_lens_c_30(mt::Input_Multislice<T>& self, T c_30) {
      self.cond_lens.c_30 = c_30 * mt::c_mm_2_Angs;
    }

    static T get_cond_lens_c_32(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.c_32;
    }

    static void set_cond_lens_c_32(mt::Input_Multislice<T>& self, T c_32) {
      self.cond_lens.c_32 = c_32;
    }

    static T get_cond_lens_phi_32(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.phi_32 / mt::c_deg_2_rad;
    }

    static void set_cond_lens_phi_32(mt::Input_Multislice<T>& self, T phi_32) {
      self.cond_lens.phi_32 = phi_32 * mt::c_deg_2_rad;
    }

    static T get_cond_lens_c_34(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.c_34;
    }

    static void set_cond_lens_c_34(mt::Input_Multislice<T>& self, T c_34) {
      self.cond_lens.c_34 = c_34;
    }

    static T get_cond_lens_phi_34(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.phi_34 / mt::c_deg_2_rad;
    }

    static void set_cond_lens_phi_34(mt::Input_Multislice<T>& self, T phi_34) {
      self.cond_lens.phi_34 = phi_34 * mt::c_deg_2_rad;
    }

    static T get_cond_lens_c_41(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.c_41;
    }

    static void set_cond_lens_c_41(mt::Input_Multislice<T>& self, T c_41) {
      self.cond_lens.c_41 = c_41;
    }

    static T get_cond_lens_phi_41(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.phi_41 / mt::c_deg_2_rad;
    }

    static void set_cond_lens_phi_41(mt::Input_Multislice<T>& self, T phi_41) {
      self.cond_lens.phi_41 = phi_41 * mt::c_deg_2_rad;
    }

    static T get_cond_lens_c_43(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.c_43;
    }

    static void set_cond_lens_c_43(mt::Input_Multislice<T>& self, T c_43) {
      self.cond_lens.c_43 = c_43;
    }

    static T get_cond_lens_phi_43(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.phi_43 / mt::c_deg_2_rad;
    }

    static void set_cond_lens_phi_43(mt::Input_Multislice<T>& self, T phi_43) {
      self.cond_lens.phi_43 = phi_43 * mt::c_deg_2_rad;
    }

    static T get_cond_lens_c_45(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.c_45;
    }

    static void set_cond_lens_c_45(mt::Input_Multislice<T>& self, T c_45) {
      self.cond_lens.c_45 = c_45;
    }

    static T get_cond_lens_phi_45(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.phi_45 / mt::c_deg_2_rad;
    }

    static void set_cond_lens_phi_45(mt::Input_Multislice<T>& self, T phi_45) {
      self.cond_lens.phi_45 = phi_45 * mt::c_deg_2_rad;
    }

    static T get_cond_lens_c_50(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.c_50 / mt::c_mm_2_Angs;
    }

    static void set_cond_lens_c_50(mt::Input_Multislice<T>& self, T c_50) {
      self.cond_lens.c_50 = c_50 * mt::c_mm_2_Angs;
    }

    static T get_cond_lens_c_52(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.c_52;
    }

    static void set_cond_lens_c_52(mt::Input_Multislice<T>& self, T c_52) {
      self.cond_lens.c_52 = c_52;
    }

    static T get_cond_lens_phi_52(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.phi_52 / mt::c_deg_2_rad;
    }

    static void set_cond_lens_phi_52(mt::Input_Multislice<T>& self, T phi_52) {
      self.cond_lens.phi_52 = phi_52 * mt::c_deg_2_rad;
    }

    static T get_cond_lens_c_54(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.c_54;
    }

    static void set_cond_lens_c_54(mt::Input_Multislice<T>& self, T c_54) {
      self.cond_lens.c_54 = c_54;
    }

    static T get_cond_lens_phi_54(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.phi_54 / mt::c_deg_2_rad;
    }

    static void set_cond_lens_phi_54(mt::Input_Multislice<T>& self, T phi_54) {
      self.cond_lens.phi_54 = phi_54 * mt::c_deg_2_rad;
    }

    static T get_cond_lens_c_56(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.c_56;
    }

    static void set_cond_lens_c_56(mt::Input_Multislice<T>& self, T c_56) {
      self.cond_lens.c_56 = c_56;
    }

    static T get_cond_lens_phi_56(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.phi_56 / mt::c_deg_2_rad;
    }

    static void set_cond_lens_phi_56(mt::Input_Multislice<T>& self, T phi_56) {
      self.cond_lens.phi_56 = phi_56 * mt::c_deg_2_rad;
    }

    static T get_cond_lens_inner_aper_ang(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.inner_aper_ang / mt::c_mrad_2_rad;
    }

    static void set_cond_lens_inner_aper_ang(mt::Input_Multislice<T>& self,
                                             T inner_aper_ang) {
      self.cond_lens.inner_aper_ang = inner_aper_ang * mt::c_mrad_2_rad;
    }

    static T get_cond_lens_outer_aper_ang(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.outer_aper_ang / mt::c_mrad_2_rad;
    }

    static void set_cond_lens_outer_aper_ang(mt::Input_Multislice<T>& self,
                                             T outer_aper_ang) {
      self.cond_lens.outer_aper_ang = outer_aper_ang * mt::c_mrad_2_rad;
    }

    static T get_cond_lens_ti_a(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.ti_a;
    }

    static void set_cond_lens_ti_a(mt::Input_Multislice<T>& self, T ti_a) {
      self.cond_lens.ti_a = ti_a;
    }

    static T get_cond_lens_ti_sigma(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.ti_sigma;
    }

    static void set_cond_lens_ti_sigma(mt::Input_Multislice<T>& self, T ti_sigma) {
      self.cond_lens.ti_sigma = ti_sigma;
    }

    static T get_cond_lens_ti_beta(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.ti_beta;
    }

    static void set_cond_lens_ti_beta(mt::Input_Multislice<T>& self, T ti_beta) {
      self.cond_lens.ti_beta = ti_beta;
    }

    static int get_cond_lens_ti_npts(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.ti_npts;
    }

    static void set_cond_lens_ti_npts(mt::Input_Multislice<T>& self, int ti_npts) {
      self.cond_lens.ti_npts = ti_npts;
    }

    static T get_cond_lens_si_a(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.si_a;
    }

    static void set_cond_lens_si_a(mt::Input_Multislice<T>& self, T si_a) {
      self.cond_lens.si_a = si_a;
      self.obj_lens.si_a = si_a;
    }

    static T get_cond_lens_si_sigma(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.si_sigma;
    }

    static void set_cond_lens_si_sigma(mt::Input_Multislice<T>& self, T si_sigma) {
      self.cond_lens.si_sigma = si_sigma;
      self.obj_lens.si_sigma = si_sigma;
    }

    static T get_cond_lens_si_beta(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.si_beta;
    }

    static void set_cond_lens_si_beta(mt::Input_Multislice<T>& self, T si_beta) {
      self.cond_lens.si_beta = si_beta;
      self.obj_lens.si_beta = si_beta;
    }

    static int get_cond_lens_si_rad_npts(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.si_rad_npts;
    }

    static void set_cond_lens_si_rad_npts(mt::Input_Multislice<T>& self,
                                          int si_rad_npts) {
      self.cond_lens.si_rad_npts = si_rad_npts;
      self.obj_lens.si_rad_npts = si_rad_npts;
    }

    static int get_cond_lens_si_azm_npts(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.si_azm_npts;
    }

    static void set_cond_lens_si_azm_npts(mt::Input_Multislice<T>& self,
                                          int si_azm_npts) {
      self.cond_lens.si_azm_npts = si_azm_npts;
      self.obj_lens.si_azm_npts = si_azm_npts;
    }

    static mt::eZero_Defocus_Type get_cond_lens_zero_defocus_type(
      const mt::Input_Multislice<T>& self) {
      return self.cond_lens.zero_defocus_type;
    }

    static void set_cond_lens_zero_defocus_type(
      mt::Input_Multislice<T>& self,
      mt::eZero_Defocus_Type zero_defocus_type) {
      self.cond_lens.zero_defocus_type = zero_defocus_type;
    }

    static T get_cond_lens_zero_defocus_plane(const mt::Input_Multislice<T>& self) {
      return self.cond_lens.zero_defocus_plane;
    }

    static void set_cond_lens_zero_defocus_plane(mt::Input_Multislice<T>& self,
                                                 T zero_defocus_plane) {
      self.cond_lens.zero_defocus_plane = zero_defocus_plane;
    }

    static int get_obj_lens_m(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.m;
    }

    static void set_obj_lens_m(mt::Input_Multislice<T>& self, int m) {
      self.obj_lens.m = m;
    }

    static T get_obj_lens_c_10(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.c_10;
    }

    static void set_obj_lens_c_10(mt::Input_Multislice<T>& self, T c_10) {
      self.obj_lens.c_10 = c_10;
    }

    static T get_obj_lens_c_12(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.c_12;
    }

    static void set_obj_lens_c_12(mt::Input_Multislice<T>& self, T c_12) {
      self.obj_lens.c_12 = c_12;
    }

    static T get_obj_lens_phi_12(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.phi_12 / mt::c_deg_2_rad;
    }

    static void set_obj_lens_phi_12(mt::Input_Multislice<T>& self, T phi_12) {
      self.obj_lens.phi_12 = phi_12 * mt::c_deg_2_rad;
    }

    static T get_obj_lens_c_21(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.c_21;
    }

    static void set_obj_lens_c_21(mt::Input_Multislice<T>& self, T c_21) {
      self.obj_lens.c_21 = c_21;
    }

    static T get_obj_lens_phi_21(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.phi_21 / mt::c_deg_2_rad;
    }

    static void set_obj_lens_phi_21(mt::Input_Multislice<T>& self, T phi_21) {
      self.obj_lens.phi_21 = phi_21 * mt::c_deg_2_rad;
    }

    static T get_obj_lens_c_23(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.c_23;
    }

    static void set_obj_lens_c_23(mt::Input_Multislice<T>& self, T c_23) {
      self.obj_lens.c_23 = c_23;
    }

    static T get_obj_lens_phi_23(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.phi_23 / mt::c_deg_2_rad;
    }

    static void set_obj_lens_phi_23(mt::Input_Multislice<T>& self, T phi_23) {
      self.obj_lens.phi_23 = phi_23 * mt::c_deg_2_rad;
    }

    static T get_obj_lens_c_30(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.c_30 / mt::c_mm_2_Angs;
    }

    static void set_obj_lens_c_30(mt::Input_Multislice<T>& self, T c_30) {
      self.obj_lens.c_30 = c_30 * mt::c_mm_2_Angs;
    }

    static T get_obj_lens_c_32(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.c_32;
    }

    static void set_obj_lens_c_32(mt::Input_Multislice<T>& self, T c_32) {
      self.obj_lens.c_32 = c_32;
    }

    static T get_obj_lens_phi_32(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.phi_32 / mt::c_deg_2_rad;
    }

    static void set_obj_lens_phi_32(mt::Input_Multislice<T>& self, T phi_32) {
      self.obj_lens.phi_32 = phi_32 * mt::c_deg_2_rad;
    }

    static T get_obj_lens_c_34(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.c_34;
    }

    static void set_obj_lens_c_34(mt::Input_Multislice<T>& self, T c_34) {
      self.obj_lens.c_34 = c_34;
    }

    static T get_obj_lens_phi_34(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.phi_34 / mt::c_deg_2_rad;
    }

    static void set_obj_lens_phi_34(mt::Input_Multislice<T>& self, T phi_34) {
      self.obj_lens.phi_34 = phi_34 * mt::c_deg_2_rad;
    }

    static T get_obj_lens_c_41(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.c_41;
    }

    static void set_obj_lens_c_41(mt::Input_Multislice<T>& self, T c_41) {
      self.obj_lens.c_41 = c_41;
    }

    static T get_obj_lens_phi_41(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.phi_41 / mt::c_deg_2_rad;
    }

    static void set_obj_lens_phi_41(mt::Input_Multislice<T>& self, T phi_41) {
      self.obj_lens.phi_41 = phi_41 * mt::c_deg_2_rad;
    }

    static T get_obj_lens_c_43(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.c_43;
    }

    static void set_obj_lens_c_43(mt::Input_Multislice<T>& self, T c_43) {
      self.obj_lens.c_43 = c_43;
    }

    static T get_obj_lens_phi_43(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.phi_43 / mt::c_deg_2_rad;
    }

    static void set_obj_lens_phi_43(mt::Input_Multislice<T>& self, T phi_43) {
      self.obj_lens.phi_43 = phi_43 * mt::c_deg_2_rad;
    }

    static T get_obj_lens_c_45(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.c_45;
    }

    static void set_obj_lens_c_45(mt::Input_Multislice<T>& self, T c_45) {
      self.obj_lens.c_45 = c_45;
    }

    static T get_obj_lens_phi_45(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.phi_45 / mt::c_deg_2_rad;
    }

    static void set_obj_lens_phi_45(mt::Input_Multislice<T>& self, T phi_45) {
      self.obj_lens.phi_45 = phi_45 * mt::c_deg_2_rad;
    }

    static T get_obj_lens_c_50(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.c_50 / mt::c_mm_2_Angs;
    }

    static void set_obj_lens_c_50(mt::Input_Multislice<T>& self, T c_50) {
      self.obj_lens.c_50 = c_50 * mt::c_mm_2_Angs;
    }

    static T get_obj_lens_c_52(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.c_52;
    }

    static void set_obj_lens_c_52(mt::Input_Multislice<T>& self, T c_52) {
      self.obj_lens.c_52 = c_52;
    }

    static T get_obj_lens_phi_52(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.phi_52 / mt::c_deg_2_rad;
    }

    static void set_obj_lens_phi_52(mt::Input_Multislice<T>& self, T phi_52) {
      self.obj_lens.phi_52 = phi_52 * mt::c_deg_2_rad;
    }

    static T get_obj_lens_c_54(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.c_54;
    }

    static void set_obj_lens_c_54(mt::Input_Multislice<T>& self, T c_54) {
      self.obj_lens.c_54 = c_54;
    }

    static T get_obj_lens_phi_54(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.phi_54 / mt::c_deg_2_rad;
    }

    static void set_obj_lens_phi_54(mt::Input_Multislice<T>& self, T phi_54) {
      self.obj_lens.phi_54 = phi_54 * mt::c_deg_2_rad;
    }

    static T get_obj_lens_c_56(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.c_56;
    }

    static void set_obj_lens_c_56(mt::Input_Multislice<T>& self, T c_56) {
      self.obj_lens.c_56 = c_56;
    }

    static T get_obj_lens_phi_56(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.phi_56 / mt::c_deg_2_rad;
    }

    static void set_obj_lens_phi_56(mt::Input_Multislice<T>& self, T phi_56) {
      self.obj_lens.phi_56 = phi_56 * mt::c_deg_2_rad;
    }

    static T get_obj_lens_inner_aper_ang(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.inner_aper_ang / mt::c_mrad_2_rad;
    }

    static void set_obj_lens_inner_aper_ang(mt::Input_Multislice<T>& self,
                                            T inner_aper_ang) {
      self.obj_lens.inner_aper_ang = inner_aper_ang * mt::c_mrad_2_rad;
    }

    static T get_obj_lens_outer_aper_ang(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.outer_aper_ang / mt::c_mrad_2_rad;
    }

    static void set_obj_lens_outer_aper_ang(mt::Input_Multislice<T>& self,
                                            T outer_aper_ang) {
      self.obj_lens.outer_aper_ang = outer_aper_ang * mt::c_mrad_2_rad;
    }

    static T get_obj_lens_ti_a(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.ti_a;
    }

    static void set_obj_lens_ti_a(mt::Input_Multislice<T>& self, T ti_a) {
      self.obj_lens.ti_a = ti_a;
    }

    static T get_obj_lens_ti_sigma(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.ti_sigma;
    }

    static void set_obj_lens_ti_sigma(mt::Input_Multislice<T>& self, T ti_sigma) {
      self.obj_lens.ti_sigma = ti_sigma;
    }

    static T get_obj_lens_ti_beta(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.ti_beta;
    }

    static void set_obj_lens_ti_beta(mt::Input_Multislice<T>& self, T ti_beta) {
      self.obj_lens.ti_beta = ti_beta;
    }

    static int get_obj_lens_ti_npts(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.ti_npts;
    }

    static void set_obj_lens_ti_npts(mt::Input_Multislice<T>& self, int ti_npts) {
      self.obj_lens.ti_npts = ti_npts;
    }

    static mt::eZero_Defocus_Type get_obj_lens_zero_defocus_type(
      const mt::Input_Multislice<T>& self) {
      return self.obj_lens.zero_defocus_type;
    }

    static void set_obj_lens_zero_defocus_type(
      mt::Input_Multislice<T>& self,
      mt::eZero_Defocus_Type zero_defocus_type) {
      self.obj_lens.zero_defocus_type = zero_defocus_type;
    }

    static T get_obj_lens_zero_defocus_plane(const mt::Input_Multislice<T>& self) {
      return self.obj_lens.zero_defocus_plane;
    }

    static void set_obj_lens_zero_defocus_plane(mt::Input_Multislice<T>& self,
                                                T zero_defocus_plane) {
      self.obj_lens.zero_defocus_plane = zero_defocus_plane;
    }

    static mt::eScanning_Type get_scanning_type(const mt::Input_Multislice<T>& self) {
      return self.scanning.type;
    }

    static void set_scanning_type(mt::Input_Multislice<T>& self,
                                  mt::eScanning_Type scanning_type) {
      self.scanning.type = scanning_type;
    }

    static bool get_scanning_periodic(const mt::Input_Multislice<T>& self) {
      return self.scanning.pbc;
    }

    static void set_scanning_periodic(mt::Input_Multislice<T>& self,
                                      bool scanning_periodic) {
      self.scanning.pbc = scanning_periodic;
    }

    static bool get_scanning_square_pxs(const mt::Input_Multislice<T>& self) {
      return self.scanning.spxs;
    }

    static void set_scanning_square_pxs(mt::Input_Multislice<T>& self,
                                    bool scanning_square_pxs) {
      self.scanning.spxs = scanning_square_pxs;
    }

    static int get_scanning_ns(const mt::Input_Multislice<T>& self) {
      return self.scanning.ns;
    }

    static void set_scanning_ns(mt::Input_Multislice<T>& self, int scanning_ns) {
      self.scanning.ns = scanning_ns;
    }

    static T get_scanning_x0(const mt::Input_Multislice<T>& self) {
      return self.scanning.x0;
    }

    static void set_scanning_x0(mt::Input_Multislice<T>& self, T scanning_x0) {
      self.scanning.x0 = scanning_x0;
    }

    static T get_scanning_y0(const mt::Input_Multislice<T>& self) {
      return self.scanning.y0;
    }

    static void set_scanning_y0(mt::Input_Multislice<T>& self, T scanning_y0) {
      self.scanning.y0 = scanning_y0;
    }

    static T get_scanning_xe(const mt::Input_Multislice<T>& self) {
      return self.scanning.xe;
    }

    static void set_scanning_xe(mt::Input_Multislice<T>& self, T scanning_xe) {
      self.scanning.xe = scanning_xe;
    }

    static T get_scanning_ye(const mt::Input_Multislice<T>& self) {
      return self.scanning.ye;
    }

    static void set_scanning_ye(mt::Input_Multislice<T>& self, T scanning_ye) {
      self.scanning.ye = scanning_ye;
    }

    static int get_ped_nrot(const mt::Input_Multislice<T>& self) {
      return self.nrot;
    }

    // FIXME Could simplify and remove
    static void set_ped_nrot(mt::Input_Multislice<T>& self, int ped_nrot) {
      self.nrot = ped_nrot;
    }

    // FIXME Could simplify and remove
    static T get_ped_theta(const mt::Input_Multislice<T>& self) {
      return self.theta / mt::c_deg_2_rad;
    }

    // FIXME Could simplify and remove
    static void set_ped_theta(mt::Input_Multislice<T>& self, T ped_theta) {
      self.theta = ped_theta * mt::c_deg_2_rad;
    }

    // FIXME Could simplify and remove
    static int get_hci_nrot(const mt::Input_Multislice<T>& self) {
      return self.nrot;
    }

    // FIXME Could simplify and remove
    static void set_hci_nrot(mt::Input_Multislice<T>& self, int hci_nrot) {
      self.nrot = hci_nrot;
    }

    // FIXME Could simplify and remove
    static T get_hci_theta(const mt::Input_Multislice<T>& self) {
      return self.theta / mt::c_deg_2_rad;
    }

    // FIXME Could simplify and remove
    static void set_hci_theta(mt::Input_Multislice<T>& self, T hci_theta) {
      self.theta = hci_theta * mt::c_deg_2_rad;
    }

    static int get_eels_Z(const mt::Input_Multislice<T>& self) {
      return self.eels_fr.Z;
    }

    static void set_eels_Z(mt::Input_Multislice<T>& self, int eels_Z) {
      self.eels_fr.Z = eels_Z;
    }

    static T get_eels_E_loss(const mt::Input_Multislice<T>& self) {
      return self.eels_fr.E_loss / mt::c_eV_2_keV;
    }

    static void set_eels_E_loss(mt::Input_Multislice<T>& self, T eels_E_loss) {
      self.eels_fr.E_loss = eels_E_loss * mt::c_eV_2_keV;
    }

    static T get_eels_collection_angle(const mt::Input_Multislice<T>& self) {
      return self.eels_fr.collection_angle / mt::c_mrad_2_rad;
    }

    static void set_eels_collection_angle(mt::Input_Multislice<T>& self,
                                          T eels_collection_angle) {
      self.eels_fr.collection_angle = eels_collection_angle * mt::c_mrad_2_rad;
    }

    static int get_eels_m_selection(const mt::Input_Multislice<T>& self) {
      return self.eels_fr.m_selection;
    }

    static void set_eels_m_selection(mt::Input_Multislice<T>& self,
                                     int eels_m_selection) {
      self.eels_fr.m_selection = eels_m_selection;
    }

    static mt::eChannelling_Type get_eels_channelling_type(
      const mt::Input_Multislice<T>& self) {
      return self.eels_fr.channelling_type;
    }

    static void set_eels_channelling_type(mt::Input_Multislice<T>& self,
                                          mt::eChannelling_Type eels_channelling_type) {
      self.eels_fr.channelling_type = eels_channelling_type;
    }

    static int get_eftem_Z(const mt::Input_Multislice<T>& self) {
      return self.eels_fr.Z;
    }

    static void set_eftem_Z(mt::Input_Multislice<T>& self, int eftem_Z) {
      self.eels_fr.Z = eftem_Z;
    }

    static T get_eftem_E_loss(const mt::Input_Multislice<T>& self) {
      return self.eels_fr.E_loss / mt::c_eV_2_keV;
    }

    static void set_eftem_E_loss(mt::Input_Multislice<T>& self, T eftem_E_loss) {
      self.eels_fr.E_loss = eftem_E_loss * mt::c_eV_2_keV;
    }

    static T get_eftem_collection_angle(const mt::Input_Multislice<T>& self) {
      return self.eels_fr.collection_angle / mt::c_mrad_2_rad;
    }

    static void set_eftem_collection_angle(mt::Input_Multislice<T>& self,
                                           T eftem_collection_angle) {
      self.eels_fr.collection_angle = eftem_collection_angle * mt::c_mrad_2_rad;
    }

    static int get_eftem_m_selection(const mt::Input_Multislice<T>& self) {
      return self.eels_fr.m_selection;
    }

    static void set_eftem_m_selection(mt::Input_Multislice<T>& self,
                                      int eftem_m_selection) {
      self.eels_fr.m_selection = eftem_m_selection;
    }

    static mt::eChannelling_Type get_eftem_channelling_type(
      const mt::Input_Multislice<T>& self) {
      return self.eels_fr.channelling_type;
    }

    static void set_eftem_channelling_type(
      mt::Input_Multislice<T>& self,
      mt::eChannelling_Type eftem_channelling_type) {
      self.eels_fr.channelling_type = eftem_channelling_type;
    }

    static int get_output_area_ix_0(const mt::Input_Multislice<T>& self) {
      return self.output_area.ix_0 + 1;
    }

    static void set_output_area_ix_0(mt::Input_Multislice<T>& self, int ix_0) {
      self.output_area.ix_0 = ix_0 - 1;
    }

    static int get_output_area_iy_0(const mt::Input_Multislice<T>& self) {
      return self.output_area.iy_0 + 1;
    }

    static void set_output_area_iy_0(mt::Input_Multislice<T>& self, int iy_0) {
      self.output_area.iy_0 = iy_0 - 1;
    }

    static int get_output_area_ix_e(const mt::Input_Multislice<T>& self) {
      return self.output_area.ix_e + 1;
    }

    static void set_output_area_ix_e(mt::Input_Multislice<T>& self, int ix_e) {
      self.output_area.ix_e = ix_e - 1;
    }

    static int get_output_area_iy_e(const mt::Input_Multislice<T>& self) {
      return self.output_area.iy_e + 1;
    }

    static void set_output_area_iy_e(mt::Input_Multislice<T>& self, int iy_e) {
      self.output_area.iy_e = iy_e - 1;
    }
  };

}}  // namespace pybind11::detail

template <typename T>
void wrap_input_multislice(py::module_ m, const char *name) {
  typedef mt::Input_Multislice<T> Type;

  // Wrap the mt::Input_Multislice class
  py::class_<Type>(m, name)
    .def_readonly("system_conf", &Type::system_conf)
    .def_readwrite("interaction_model", &Type::interaction_model)
    .def_readwrite("potential_type", &Type::potential_type)
    .def_readwrite("pn_model", &Type::pn_model)
    .def_readwrite("pn_coh_contrib", &Type::pn_coh_contrib)
    .def_readwrite("pn_single_conf", &Type::pn_single_conf)
    .def_readwrite("fp_dist", &Type::fp_dist)
    .def_readwrite("pn_seed", &Type::pn_seed)
    .def_readwrite("pn_dim", &Type::pn_dim)
    .def_readwrite("pn_nconf", &Type::pn_nconf)
    .def_readwrite("fp_iconf_0", &Type::fp_iconf_0)
    .def_readwrite("atoms", &Type::atoms)
    .def_readwrite("is_crystal", &Type::is_crystal)
    .def_readwrite("spec_rot_center_type", &Type::spec_rot_center_type)
    .def_readwrite("spec_rot_center_p", &Type::spec_rot_center_p)
    .def_readwrite("thick_type", &Type::thick_type)
    .def_readwrite("thick", &Type::thick)
    .def_readwrite("potential_slicing", &Type::potential_slicing)
    .def_readwrite("grid_2d", &Type::grid_2d)
    .def_readwrite("output_area", &Type::output_area)
    .def_readwrite("simulation_type", &Type::simulation_type)
    .def_readwrite("iw_psi", &Type::iw_psi)
    .def_readwrite("iw_x", &Type::iw_x)
    .def_readwrite("iw_y", &Type::iw_y)
    .def_readwrite("E_0", &Type::E_0)
    .def_readwrite("lambda_", &Type::lambda)
    .def_readwrite("illumination_model", &Type::illumination_model)
    .def_readwrite("temporal_spatial_incoh", &Type::temporal_spatial_incoh)
    .def_readwrite("cond_lens", &Type::cond_lens)
    .def_readwrite("obj_lens", &Type::obj_lens)
    .def_readwrite("scanning", &Type::scanning)
    .def_readwrite("detector", &Type::detector)
    .def_readwrite("eels_fr", &Type::eels_fr)
    .def_readwrite("operation_mode", &Type::operation_mode)
    .def_readwrite("slice_storage", &Type::slice_storage)
    .def_readwrite("reverse_multislice", &Type::reverse_multislice)
    .def_readwrite("mul_sign", &Type::mul_sign)
    .def_readwrite("Vrl", &Type::Vrl)
    .def_readwrite("nR", &Type::nR)
    .def_readwrite("nrot", &Type::nrot)
    .def_readwrite("cdl_var_type", &Type::cdl_var_type)
    .def_readwrite("cdl_var", &Type::cdl_var)
    .def_readwrite("iscan", &Type::iscan)
    .def_readwrite("beam_x", &Type::beam_x)
    .def_readwrite("beam_y", &Type::beam_y)
    .def_readwrite("islice", &Type::islice)
    .def_readwrite("dp_Shift", &Type::dp_Shift)
    .def("validate_parameters", &Type::validate_parameters)
    .def(py::pickle(&py::detail::Helpers<Type>::getstate,
                    &py::detail::Helpers<Type>::setstate))

    // The matlab interface defines some of these differently
    .def_property("theta",
                  &py::detail::Helpers<Type>::get_theta,
                  &py::detail::Helpers<Type>::set_theta)
    .def_property(
      "phi", &py::detail::Helpers<Type>::get_phi, &py::detail::Helpers<Type>::set_phi)
    .def_property("spec_rot_theta",
                  &py::detail::Helpers<Type>::get_spec_rot_theta,
                  &py::detail::Helpers<Type>::set_spec_rot_theta)
    .def_property("spec_rot_u0",
                  &py::detail::Helpers<Type>::get_spec_rot_u0,
                  &py::detail::Helpers<Type>::set_spec_rot_u0)
    .def_property("iw_type",
                  &py::detail::Helpers<Type>::get_iw_type,
                  &py::detail::Helpers<Type>::set_iw_type)

    // Additional properties to be compatible with matlab interface
    .def_property("spec_atoms",
                  &py::detail::Helpers<Type>::get_spec_atoms,
                  &py::detail::Helpers<Type>::set_spec_atoms)
    .def_property("spec_dz",
                  &py::detail::Helpers<Type>::get_spec_dz,
                  &py::detail::Helpers<Type>::set_spec_dz)
    .def_property("spec_lx",
                  &py::detail::Helpers<Type>::get_spec_lx,
                  &py::detail::Helpers<Type>::set_spec_lx)
    .def_property("spec_ly",
                  &py::detail::Helpers<Type>::get_spec_ly,
                  &py::detail::Helpers<Type>::set_spec_ly)
    .def_property("spec_lz",
                  &py::detail::Helpers<Type>::get_spec_lz,
                  &py::detail::Helpers<Type>::set_spec_lz)
    .def_property("spec_cryst_na",
                  &py::detail::Helpers<Type>::get_spec_cryst_na,
                  &py::detail::Helpers<Type>::set_spec_cryst_na)
    .def_property("spec_cryst_nb",
                  &py::detail::Helpers<Type>::get_spec_cryst_nb,
                  &py::detail::Helpers<Type>::set_spec_cryst_nb)
    .def_property("spec_cryst_nc",
                  &py::detail::Helpers<Type>::get_spec_cryst_nc,
                  &py::detail::Helpers<Type>::set_spec_cryst_nc)
    .def_property("spec_cryst_a",
                  &py::detail::Helpers<Type>::get_spec_cryst_a,
                  &py::detail::Helpers<Type>::set_spec_cryst_a)
    .def_property("spec_cryst_b",
                  &py::detail::Helpers<Type>::get_spec_cryst_b,
                  &py::detail::Helpers<Type>::set_spec_cryst_b)
    .def_property("spec_cryst_c",
                  &py::detail::Helpers<Type>::get_spec_cryst_c,
                  &py::detail::Helpers<Type>::set_spec_cryst_c)
    .def_property("spec_cryst_x0",
                  &py::detail::Helpers<Type>::get_spec_cryst_x0,
                  &py::detail::Helpers<Type>::set_spec_cryst_x0)
    .def_property("spec_cryst_y0",
                  &py::detail::Helpers<Type>::get_spec_cryst_y0,
                  &py::detail::Helpers<Type>::set_spec_cryst_y0)
    .def_property("spec_amorp",
                  &py::detail::Helpers<Type>::get_spec_amorp,
                  &py::detail::Helpers<Type>::set_spec_amorp)
    .def_property(
      "nx", &py::detail::Helpers<Type>::get_nx, &py::detail::Helpers<Type>::set_nx)
    .def_property(
      "nx", &py::detail::Helpers<Type>::get_nx, &py::detail::Helpers<Type>::set_nx)
    .def_property(
      "ny", &py::detail::Helpers<Type>::get_ny, &py::detail::Helpers<Type>::set_ny)
    .def_property(
      "bwl", &py::detail::Helpers<Type>::get_bwl, &py::detail::Helpers<Type>::set_bwl)
    .def_property("cond_lens_m",
                  &py::detail::Helpers<Type>::get_cond_lens_m,
                  &py::detail::Helpers<Type>::set_cond_lens_m)
    .def_property("cond_lens_c_10",
                  &py::detail::Helpers<Type>::get_cond_lens_c_10,
                  &py::detail::Helpers<Type>::set_cond_lens_c_10)
    .def_property("cond_lens_c_12",
                  &py::detail::Helpers<Type>::get_cond_lens_c_12,
                  &py::detail::Helpers<Type>::set_cond_lens_c_12)
    .def_property("cond_lens_phi_12",
                  &py::detail::Helpers<Type>::get_cond_lens_phi_12,
                  &py::detail::Helpers<Type>::set_cond_lens_phi_12)
    .def_property("cond_lens_c_21",
                  &py::detail::Helpers<Type>::get_cond_lens_c_21,
                  &py::detail::Helpers<Type>::set_cond_lens_c_21)
    .def_property("cond_lens_phi_21",
                  &py::detail::Helpers<Type>::get_cond_lens_phi_21,
                  &py::detail::Helpers<Type>::set_cond_lens_phi_21)
    .def_property("cond_lens_c_23",
                  &py::detail::Helpers<Type>::get_cond_lens_c_23,
                  &py::detail::Helpers<Type>::set_cond_lens_c_23)
    .def_property("cond_lens_phi_23",
                  &py::detail::Helpers<Type>::get_cond_lens_phi_23,
                  &py::detail::Helpers<Type>::set_cond_lens_phi_23)
    .def_property("cond_lens_c_30",
                  &py::detail::Helpers<Type>::get_cond_lens_c_30,
                  &py::detail::Helpers<Type>::set_cond_lens_c_30)
    .def_property("cond_lens_c_32",
                  &py::detail::Helpers<Type>::get_cond_lens_c_32,
                  &py::detail::Helpers<Type>::set_cond_lens_c_32)
    .def_property("cond_lens_phi_32",
                  &py::detail::Helpers<Type>::get_cond_lens_phi_32,
                  &py::detail::Helpers<Type>::set_cond_lens_phi_32)
    .def_property("cond_lens_c_34",
                  &py::detail::Helpers<Type>::get_cond_lens_c_34,
                  &py::detail::Helpers<Type>::set_cond_lens_c_34)
    .def_property("cond_lens_phi_34",
                  &py::detail::Helpers<Type>::get_cond_lens_phi_34,
                  &py::detail::Helpers<Type>::set_cond_lens_phi_34)
    .def_property("cond_lens_c_41",
                  &py::detail::Helpers<Type>::get_cond_lens_c_41,
                  &py::detail::Helpers<Type>::set_cond_lens_c_41)
    .def_property("cond_lens_phi_41",
                  &py::detail::Helpers<Type>::get_cond_lens_phi_41,
                  &py::detail::Helpers<Type>::set_cond_lens_phi_41)
    .def_property("cond_lens_c_43",
                  &py::detail::Helpers<Type>::get_cond_lens_c_43,
                  &py::detail::Helpers<Type>::set_cond_lens_c_43)
    .def_property("cond_lens_phi_43",
                  &py::detail::Helpers<Type>::get_cond_lens_phi_43,
                  &py::detail::Helpers<Type>::set_cond_lens_phi_43)
    .def_property("cond_lens_c_45",
                  &py::detail::Helpers<Type>::get_cond_lens_c_45,
                  &py::detail::Helpers<Type>::set_cond_lens_c_45)
    .def_property("cond_lens_phi_45",
                  &py::detail::Helpers<Type>::get_cond_lens_phi_45,
                  &py::detail::Helpers<Type>::set_cond_lens_phi_45)
    .def_property("cond_lens_c_50",
                  &py::detail::Helpers<Type>::get_cond_lens_c_50,
                  &py::detail::Helpers<Type>::set_cond_lens_c_50)
    .def_property("cond_lens_c_52",
                  &py::detail::Helpers<Type>::get_cond_lens_c_52,
                  &py::detail::Helpers<Type>::set_cond_lens_c_52)
    .def_property("cond_lens_phi_52",
                  &py::detail::Helpers<Type>::get_cond_lens_phi_52,
                  &py::detail::Helpers<Type>::set_cond_lens_phi_52)
    .def_property("cond_lens_c_54",
                  &py::detail::Helpers<Type>::get_cond_lens_c_54,
                  &py::detail::Helpers<Type>::set_cond_lens_c_54)
    .def_property("cond_lens_phi_54",
                  &py::detail::Helpers<Type>::get_cond_lens_phi_54,
                  &py::detail::Helpers<Type>::set_cond_lens_phi_54)
    .def_property("cond_lens_c_56",
                  &py::detail::Helpers<Type>::get_cond_lens_c_56,
                  &py::detail::Helpers<Type>::set_cond_lens_c_56)
    .def_property("cond_lens_phi_56",
                  &py::detail::Helpers<Type>::get_cond_lens_phi_56,
                  &py::detail::Helpers<Type>::set_cond_lens_phi_56)
    .def_property("cond_lens_inner_aper_ang",
                  &py::detail::Helpers<Type>::get_cond_lens_inner_aper_ang,
                  &py::detail::Helpers<Type>::set_cond_lens_inner_aper_ang)
    .def_property("cond_lens_outer_aper_ang",
                  &py::detail::Helpers<Type>::get_cond_lens_outer_aper_ang,
                  &py::detail::Helpers<Type>::set_cond_lens_outer_aper_ang)
    .def_property("cond_lens_ti_a",
                  &py::detail::Helpers<Type>::get_cond_lens_ti_a,
                  &py::detail::Helpers<Type>::set_cond_lens_ti_a)
    .def_property("cond_lens_ti_sigma",
                  &py::detail::Helpers<Type>::get_cond_lens_ti_sigma,
                  &py::detail::Helpers<Type>::set_cond_lens_ti_sigma)
    .def_property("cond_lens_ti_beta",
                  &py::detail::Helpers<Type>::get_cond_lens_ti_beta,
                  &py::detail::Helpers<Type>::set_cond_lens_ti_beta)
    .def_property("cond_lens_ti_npts",
                  &py::detail::Helpers<Type>::get_cond_lens_ti_npts,
                  &py::detail::Helpers<Type>::set_cond_lens_ti_npts)
    .def_property("cond_lens_si_a",
                  &py::detail::Helpers<Type>::get_cond_lens_si_a,
                  &py::detail::Helpers<Type>::set_cond_lens_si_a)
    .def_property("cond_lens_si_sigma",
                  &py::detail::Helpers<Type>::get_cond_lens_si_sigma,
                  &py::detail::Helpers<Type>::set_cond_lens_si_sigma)
    .def_property("cond_lens_si_beta",
                  &py::detail::Helpers<Type>::get_cond_lens_si_beta,
                  &py::detail::Helpers<Type>::set_cond_lens_si_beta)
    .def_property("cond_lens_si_rad_npts",
                  &py::detail::Helpers<Type>::get_cond_lens_si_rad_npts,
                  &py::detail::Helpers<Type>::set_cond_lens_si_rad_npts)
    .def_property("cond_lens_si_azm_npts",
                  &py::detail::Helpers<Type>::get_cond_lens_si_azm_npts,
                  &py::detail::Helpers<Type>::set_cond_lens_si_azm_npts)
    .def_property("cond_lens_zero_defocus_type",
                  &py::detail::Helpers<Type>::get_cond_lens_zero_defocus_type,
                  &py::detail::Helpers<Type>::set_cond_lens_zero_defocus_type)
    .def_property("cond_lens_zero_defocus_plane",
                  &py::detail::Helpers<Type>::get_cond_lens_zero_defocus_plane,
                  &py::detail::Helpers<Type>::set_cond_lens_zero_defocus_plane)
    .def_property("obj_lens_m",
                  &py::detail::Helpers<Type>::get_obj_lens_m,
                  &py::detail::Helpers<Type>::set_obj_lens_m)
    .def_property("obj_lens_c_10",
                  &py::detail::Helpers<Type>::get_obj_lens_c_10,
                  &py::detail::Helpers<Type>::set_obj_lens_c_10)
    .def_property("obj_lens_c_12",
                  &py::detail::Helpers<Type>::get_obj_lens_c_12,
                  &py::detail::Helpers<Type>::set_obj_lens_c_12)
    .def_property("obj_lens_phi_12",
                  &py::detail::Helpers<Type>::get_obj_lens_phi_12,
                  &py::detail::Helpers<Type>::set_obj_lens_phi_12)
    .def_property("obj_lens_c_21",
                  &py::detail::Helpers<Type>::get_obj_lens_c_21,
                  &py::detail::Helpers<Type>::set_obj_lens_c_21)
    .def_property("obj_lens_phi_21",
                  &py::detail::Helpers<Type>::get_obj_lens_phi_21,
                  &py::detail::Helpers<Type>::set_obj_lens_phi_21)
    .def_property("obj_lens_c_23",
                  &py::detail::Helpers<Type>::get_obj_lens_c_23,
                  &py::detail::Helpers<Type>::set_obj_lens_c_23)
    .def_property("obj_lens_phi_23",
                  &py::detail::Helpers<Type>::get_obj_lens_phi_23,
                  &py::detail::Helpers<Type>::set_obj_lens_phi_23)
    .def_property("obj_lens_c_30",
                  &py::detail::Helpers<Type>::get_obj_lens_c_30,
                  &py::detail::Helpers<Type>::set_obj_lens_c_30)
    .def_property("obj_lens_c_32",
                  &py::detail::Helpers<Type>::get_obj_lens_c_32,
                  &py::detail::Helpers<Type>::set_obj_lens_c_32)
    .def_property("obj_lens_phi_32",
                  &py::detail::Helpers<Type>::get_obj_lens_phi_32,
                  &py::detail::Helpers<Type>::set_obj_lens_phi_32)
    .def_property("obj_lens_c_34",
                  &py::detail::Helpers<Type>::get_obj_lens_c_34,
                  &py::detail::Helpers<Type>::set_obj_lens_c_34)
    .def_property("obj_lens_phi_34",
                  &py::detail::Helpers<Type>::get_obj_lens_phi_34,
                  &py::detail::Helpers<Type>::set_obj_lens_phi_34)
    .def_property("obj_lens_c_41",
                  &py::detail::Helpers<Type>::get_obj_lens_c_41,
                  &py::detail::Helpers<Type>::set_obj_lens_c_41)
    .def_property("obj_lens_phi_41",
                  &py::detail::Helpers<Type>::get_obj_lens_phi_41,
                  &py::detail::Helpers<Type>::set_obj_lens_phi_41)
    .def_property("obj_lens_c_43",
                  &py::detail::Helpers<Type>::get_obj_lens_c_43,
                  &py::detail::Helpers<Type>::set_obj_lens_c_43)
    .def_property("obj_lens_phi_43",
                  &py::detail::Helpers<Type>::get_obj_lens_phi_43,
                  &py::detail::Helpers<Type>::set_obj_lens_phi_43)
    .def_property("obj_lens_c_45",
                  &py::detail::Helpers<Type>::get_obj_lens_c_45,
                  &py::detail::Helpers<Type>::set_obj_lens_c_45)
    .def_property("obj_lens_phi_45",
                  &py::detail::Helpers<Type>::get_obj_lens_phi_45,
                  &py::detail::Helpers<Type>::set_obj_lens_phi_45)
    .def_property("obj_lens_c_50",
                  &py::detail::Helpers<Type>::get_obj_lens_c_50,
                  &py::detail::Helpers<Type>::set_obj_lens_c_50)
    .def_property("obj_lens_c_52",
                  &py::detail::Helpers<Type>::get_obj_lens_c_52,
                  &py::detail::Helpers<Type>::set_obj_lens_c_52)
    .def_property("obj_lens_phi_52",
                  &py::detail::Helpers<Type>::get_obj_lens_phi_52,
                  &py::detail::Helpers<Type>::set_obj_lens_phi_52)
    .def_property("obj_lens_c_54",
                  &py::detail::Helpers<Type>::get_obj_lens_c_54,
                  &py::detail::Helpers<Type>::set_obj_lens_c_54)
    .def_property("obj_lens_phi_54",
                  &py::detail::Helpers<Type>::get_obj_lens_phi_54,
                  &py::detail::Helpers<Type>::set_obj_lens_phi_54)
    .def_property("obj_lens_c_56",
                  &py::detail::Helpers<Type>::get_obj_lens_c_56,
                  &py::detail::Helpers<Type>::set_obj_lens_c_56)
    .def_property("obj_lens_phi_56",
                  &py::detail::Helpers<Type>::get_obj_lens_phi_56,
                  &py::detail::Helpers<Type>::set_obj_lens_phi_56)
    .def_property("obj_lens_inner_aper_ang",
                  &py::detail::Helpers<Type>::get_obj_lens_inner_aper_ang,
                  &py::detail::Helpers<Type>::set_obj_lens_inner_aper_ang)
    .def_property("obj_lens_outer_aper_ang",
                  &py::detail::Helpers<Type>::get_obj_lens_outer_aper_ang,
                  &py::detail::Helpers<Type>::set_obj_lens_outer_aper_ang)
    .def_property("obj_lens_ti_a",
                  &py::detail::Helpers<Type>::get_obj_lens_ti_a,
                  &py::detail::Helpers<Type>::set_obj_lens_ti_a)
    .def_property("obj_lens_ti_sigma",
                  &py::detail::Helpers<Type>::get_obj_lens_ti_sigma,
                  &py::detail::Helpers<Type>::set_obj_lens_ti_sigma)
    .def_property("obj_lens_ti_beta",
                  &py::detail::Helpers<Type>::get_obj_lens_ti_beta,
                  &py::detail::Helpers<Type>::set_obj_lens_ti_beta)
    .def_property("obj_lens_ti_npts",
                  &py::detail::Helpers<Type>::get_obj_lens_ti_npts,
                  &py::detail::Helpers<Type>::set_obj_lens_ti_npts)
    .def_property("obj_lens_zero_defocus_type",
                  &py::detail::Helpers<Type>::get_obj_lens_zero_defocus_type,
                  &py::detail::Helpers<Type>::set_obj_lens_zero_defocus_type)
    .def_property("obj_lens_zero_defocus_plane",
                  &py::detail::Helpers<Type>::get_obj_lens_zero_defocus_plane,
                  &py::detail::Helpers<Type>::set_obj_lens_zero_defocus_plane)
    .def_property("scanning_type",
                  &py::detail::Helpers<Type>::get_scanning_type,
                  &py::detail::Helpers<Type>::set_scanning_type)
    .def_property("scanning_periodic",
                  &py::detail::Helpers<Type>::get_scanning_periodic,
                  &py::detail::Helpers<Type>::set_scanning_periodic)
    .def_property("scanning_square_pxs",
                  &py::detail::Helpers<Type>::get_scanning_square_pxs,
                  &py::detail::Helpers<Type>::set_scanning_square_pxs)
    .def_property("scanning_ns",
                  &py::detail::Helpers<Type>::get_scanning_ns,
                  &py::detail::Helpers<Type>::set_scanning_ns)
    .def_property("scanning_x0",
                  &py::detail::Helpers<Type>::get_scanning_x0,
                  &py::detail::Helpers<Type>::set_scanning_x0)
    .def_property("scanning_y0",
                  &py::detail::Helpers<Type>::get_scanning_y0,
                  &py::detail::Helpers<Type>::set_scanning_y0)
    .def_property("scanning_xe",
                  &py::detail::Helpers<Type>::get_scanning_xe,
                  &py::detail::Helpers<Type>::set_scanning_xe)
    .def_property("scanning_ye",
                  &py::detail::Helpers<Type>::get_scanning_ye,
                  &py::detail::Helpers<Type>::set_scanning_ye)
    .def_property("ped_nrot",
                  &py::detail::Helpers<Type>::get_ped_nrot,
                  &py::detail::Helpers<Type>::set_ped_nrot)
    .def_property("ped_theta",
                  &py::detail::Helpers<Type>::get_ped_theta,
                  &py::detail::Helpers<Type>::set_ped_theta)
    .def_property("hci_nrot",
                  &py::detail::Helpers<Type>::get_hci_nrot,
                  &py::detail::Helpers<Type>::set_hci_nrot)
    .def_property("hci_theta",
                  &py::detail::Helpers<Type>::get_hci_theta,
                  &py::detail::Helpers<Type>::set_hci_theta)
    .def_property("eels_Z",
                  &py::detail::Helpers<Type>::get_eels_Z,
                  &py::detail::Helpers<Type>::set_eels_Z)
    .def_property("eels_E_loss",
                  &py::detail::Helpers<Type>::get_eels_E_loss,
                  &py::detail::Helpers<Type>::set_eels_E_loss)
    .def_property("eels_collection_angle",
                  &py::detail::Helpers<Type>::get_eels_collection_angle,
                  &py::detail::Helpers<Type>::set_eels_collection_angle)
    .def_property("eels_m_selection",
                  &py::detail::Helpers<Type>::get_eels_m_selection,
                  &py::detail::Helpers<Type>::set_eels_m_selection)
    .def_property("eels_channelling_type",
                  &py::detail::Helpers<Type>::get_eels_channelling_type,
                  &py::detail::Helpers<Type>::set_eels_channelling_type)
    .def_property("eftem_Z",
                  &py::detail::Helpers<Type>::get_eftem_Z,
                  &py::detail::Helpers<Type>::set_eftem_Z)
    .def_property("eftem_E_loss",
                  &py::detail::Helpers<Type>::get_eftem_E_loss,
                  &py::detail::Helpers<Type>::set_eftem_E_loss)
    .def_property("eftem_collection_angle",
                  &py::detail::Helpers<Type>::get_eftem_collection_angle,
                  &py::detail::Helpers<Type>::set_eftem_collection_angle)
    .def_property("eftem_m_selection",
                  &py::detail::Helpers<Type>::get_eftem_m_selection,
                  &py::detail::Helpers<Type>::set_eftem_m_selection)
    .def_property("eftem_channelling_type",
                  &py::detail::Helpers<Type>::get_eftem_channelling_type,
                  &py::detail::Helpers<Type>::set_eftem_channelling_type)
    .def_property("output_area_ix_0",
                  &py::detail::Helpers<Type>::get_output_area_ix_0,
                  &py::detail::Helpers<Type>::set_output_area_ix_0)
    .def_property("output_area_iy_0",
                  &py::detail::Helpers<Type>::get_output_area_iy_0,
                  &py::detail::Helpers<Type>::set_output_area_iy_0)
    .def_property("output_area_ix_e",
                  &py::detail::Helpers<Type>::get_output_area_ix_e,
                  &py::detail::Helpers<Type>::set_output_area_ix_e)
    .def_property("output_area_iy_e",
                  &py::detail::Helpers<Type>::get_output_area_iy_e,
                  &py::detail::Helpers<Type>::set_output_area_iy_e);
}

void export_input_multislice(py::module_ m) {

  m.def("Input_Multislice", &py::detail::Input_Multislice_constructor);

  wrap_input_multislice<float>(m, "Input_Multislice_f");
  wrap_input_multislice<double>(m, "Input_Multislice_d");
}

#endif
