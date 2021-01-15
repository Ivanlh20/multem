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
#ifndef MULTEM_PYTHON_ENUMS_H
#define MULTEM_PYTHON_ENUMS_H

#include <pybind11/pybind11.h>
#include <multem.h>
#include <multem/ext/serialization.h>

namespace py = pybind11;

namespace pybind11 { namespace detail {

  /**
   * Allow conversion from string to enum
   */
  template <typename T>
  enum_<T> enum_wrapper(const handle& scope, const char* name) {
    py::enum_<T> obj(scope, name);
    obj.def("__init__", [](T& v, const std::string& str) {
      object obj = py::cast(v);
      py::dict entries = obj.get_type().attr("__entries");
      for (auto item : entries) {
        if (str == item.first.cast<std::string>()) {
          v = item.second.cast<py::tuple>()[0].cast<T>();
          return;
        }
      }
      std::string tname = typeid(T).name();
      detail::clean_type_id(tname);
      throw value_error("\"" + std::string(str)
                        + "\" is not a valid value for enum type " + tname);
    });

    // Allow implicit conversions from string to enum
    implicitly_convertible<str, T>();
    return obj;
  }

}}  // namespace pybind11::detail

void export_enums(py::module m) {
  py::detail::enum_wrapper<mt::eDevice>(m, "eDevice")
    .value("host", mt::e_host)
    .value("device", mt::e_device)
    .value("host_device", mt::e_host_device)
    .export_values();

  py::detail::enum_wrapper<mt::ePrecision>(m, "ePrecision")
    .value("float", mt::eP_float)
    .value("double", mt::eP_double)
    .export_values();

  py::detail::enum_wrapper<mt::eIllumination_Model>(m, "eIllumination_Model")
    .value("Coherent", mt::eIM_Coherent)
    .value("Partial_Coherent", mt::eIM_Partial_Coherent)
    .value("Trans_Cross_Coef", mt::eIM_Trans_Cross_Coef)
    .value("Full_Integration", mt::eIM_Full_Integration)
    .value("none", mt::eIM_none);

  py::detail::enum_wrapper<mt::eTemporal_Spatial_Incoh>(m, "eTemporal_Spatial_Incoh")
    .value("Temporal_Spatial", mt::eTSI_Temporal_Spatial)
    .value("Temporal", mt::eTSI_Temporal)
    .value("Spatial", mt::eTSI_Spatial)
    .value("none", mt::eTSI_none);

  py::detail::enum_wrapper<mt::eOperation_Mode>(m, "eOperation_Mode")
    .value("Normal", mt::eOM_Normal)
    .value("Advanced", mt::eOM_Advanced);

  py::detail::enum_wrapper<mt::eLens_Var_Type>(m, "eLens_Var_Type")
    .value("off", mt::eLVT_off)
    .value("m", mt::eLVT_m)
    .value("f", mt::eLVT_f)
    .value("Cs3", mt::eLVT_Cs3)
    .value("Cs5", mt::eLVT_Cs5)
    .value("mfa2", mt::eLVT_mfa2)
    .value("afa2", mt::eLVT_afa2)
    .value("mfa3", mt::eLVT_mfa3)
    .value("afa3", mt::eLVT_afa3)
    .value("inner_aper_ang", mt::eLVT_inner_aper_ang)
    .value("outer_aper_ang", mt::eLVT_outer_aper_ang);

  py::detail::enum_wrapper<mt::eTEM_Sim_Type>(m, "eTEM_Sim_Type")
    .value("STEM", mt::eTEMST_STEM)
    .value("ISTEM", mt::eTEMST_ISTEM)
    .value("CBED", mt::eTEMST_CBED)
    .value("CBEI", mt::eTEMST_CBEI)
    .value("ED", mt::eTEMST_ED)
    .value("HRTEM", mt::eTEMST_HRTEM)
    .value("PED", mt::eTEMST_PED)
    .value("HCTEM", mt::eTEMST_HCTEM)
    .value("EWFS", mt::eTEMST_EWFS)
    .value("EWRS", mt::eTEMST_EWRS)
    .value("EELS", mt::eTEMST_EELS)
    .value("EFTEM", mt::eTEMST_EFTEM)
    .value("IWFS", mt::eTEMST_IWFS)
    .value("IWRS", mt::eTEMST_IWRS)
    .value("PPFS", mt::eTEMST_PPFS)
    .value("PPRS", mt::eTEMST_PPRS)
    .value("TFFS", mt::eTEMST_TFFS)
    .value("TFRS", mt::eTEMST_TFRS)
    .value("PropFS", mt::eTEMST_PropFS)
    .value("PropRS", mt::eTEMST_PropRS);

  py::detail::enum_wrapper<mt::eTEM_Output_Type>(m, "eTEM_Output_Type")
    .value("image_tot_coh", mt::eTEMOT_image_tot_coh)
    .value("image_tot", mt::eTEMOT_image_tot)
    .value("m2psi_tot_coh", mt::eTEMOT_m2psi_tot_coh)
    .value("m2psi_tot", mt::eTEMOT_m2psi_tot)
    .value("m2psi_tot_psi_coh", mt::eTEMOT_m2psi_tot_psi_coh)
    .value("psi_coh", mt::eTEMOT_psi_coh)
    .value("psi_0", mt::eTEMOT_psi_0)
    .value("V", mt::eTEMOT_V)
    .value("trans", mt::eTEMOT_trans);

  py::detail::enum_wrapper<mt::eElec_Spec_Int_Model>(m, "eElec_Spec_Int_Model")
    .value("Multislice", mt::eESIM_Multislice)
    .value("Phase_Object", mt::eESIM_Phase_Object)
    .value("Weak_Phase_Object", mt::eESIM_Weak_Phase_Object);

  py::detail::enum_wrapper<mt::ePhonon_Model>(m, "ePhonon_Model")
    .value("Still_Atom", mt::ePM_Still_Atom)
    .value("Absorptive_Model", mt::ePM_Absorptive_Model)
    .value("Frozen_Phonon", mt::ePM_Frozen_Phonon);

  py::detail::enum_wrapper<mt::ePhonon_Model_Output>(m, "ePhonon_Model_Output")
    .value("Total", mt::eFMO_Total)
    .value("Coherent", mt::eFMO_Coherent);

  py::detail::enum_wrapper<mt::ePotential_Slicing>(m, "ePotential_Slicing")
    .value("Planes", mt::ePS_Planes)
    .value("dz_Proj", mt::ePS_dz_Proj)
    .value("dz_Sub", mt::ePS_dz_Sub)
    .value("Auto", mt::ePS_Auto);

  py::detail::enum_wrapper<mt::ePotential_Type>(m, "ePotential_Type")
    .value("Doyle_0_4", mt::ePT_Doyle_0_4)
    .value("Peng_0_4", mt::ePT_Peng_0_4)
    .value("Peng_0_12", mt::ePT_Peng_0_12)
    .value("Kirkland_0_12", mt::ePT_Kirkland_0_12)
    .value("Weickenmeier_0_12", mt::ePT_Weickenmeier_0_12)
    .value("Lobato_0_12", mt::ePT_Lobato_0_12)
    .value("none", mt::ePT_none);

  py::detail::enum_wrapper<mt::eIncident_Wave_Type>(m, "eIncident_Wave_Type")
    .value("Plane_Wave", mt::eIWT_Plane_Wave)
    .value("Convergent_Wave", mt::eIWT_Convergent_Wave)
    .value("User_Define_Wave", mt::eIWT_User_Define_Wave)
    .value("Auto", mt::eIWT_Auto);

  py::detail::enum_wrapper<mt::eRot_Point_Type>(m, "eRot_Point_Type")
    .value("geometric_center", mt::eRPT_geometric_center)
    .value("User_Define", mt::eRPT_User_Define);

  py::detail::enum_wrapper<mt::eSpace>(m, "eSpace")
    .value("Real", mt::eS_Real)
    .value("Reciprocal", mt::eS_Reciprocal);

  py::detail::enum_wrapper<mt::eMatch_Border>(m, "eMatch_Border")
    .value("Min", mt::eMB_Min)
    .value("Max", mt::eMB_Max)
    .value("MinMax", mt::eMB_MinMax);

  py::detail::enum_wrapper<mt::eAmorp_Lay_Type>(m, "eAmorp_Lay_Type")
    .value("Top", mt::eALT_Top)
    .value("Bottom", mt::eALT_Bottom)
    .value("Middle", mt::eALT_Middle)
    .value("none", mt::eALT_none);

  py::detail::enum_wrapper<mt::eZero_Defocus_Type>(m, "eZero_Defocus_Type")
    .value("First", mt::eZDT_First)
    .value("Middle", mt::eZDT_Middle)
    .value("Last", mt::eZDT_Last)
    .value("User_Define", mt::eZDT_User_Define);

  py::detail::enum_wrapper<mt::eThick_Type>(m, "eThick_Type")
    .value("Whole_Spec", mt::eTT_Whole_Spec)
    .value("Through_Thick", mt::eTT_Through_Thick)
    .value("Through_Slices", mt::eTT_Through_Slices);

  py::detail::enum_wrapper<mt::eScanning_Type>(m, "eScanning_Type")
    .value("Line", mt::eST_Line)
    .value("Area", mt::eST_Area);

  py::detail::enum_wrapper<mt::eGrid_Type>(m, "eGrid_Type")
    .value("Regular", mt::eGT_Regular)
    .value("Quadratic", mt::eGT_Quadratic);

  py::detail::enum_wrapper<mt::eDetector_Type>(m, "eDetector_Type")
    .value("Circular", mt::eDT_Circular)
    .value("Radial", mt::eDT_Radial)
    .value("Matrix", mt::eDT_Matrix);

  py::detail::enum_wrapper<mt::eChannelling_Type>(m, "eChannelling_Type")
    .value("Single_Channelling", mt::eCT_Single_Channelling)
    .value("Mixed_Channelling", mt::eCT_Mixed_Channelling)
    .value("Double_Channelling", mt::eCT_Double_Channelling);
}

#endif
