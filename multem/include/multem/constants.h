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

#ifndef MULTEM_CONSTANTS_H
#define MULTEM_CONSTANTS_H

namespace mt {
  const double c_Ha = 27.2113850656389;        // Hartree to electron-Volt
  const double c_a0 = 0.52917721077817892;     // Bohr radius
  const double c_Potf = 47.877645145863056;    //
  const double c_2Pi2a0 = 10.445539456905012;  // 2*pi^2*a0

  const double c_H = 6.62606876e-34;    // Planck's constant - J s
  const double c_bH = 1.054571596e-34;  // h/(2*pi) - J s
  const double c_C = 2.99792458e+8;     // Velocity of light - m s^-1
  const double c_Qe = 1.602176462e-19;  // Elementary charge
  const double c_me = 9.10938291e-31;   // Electron rest mass [kg]
  const double c_mp = 1.672621637e-27;  // Proton rest mass [kg]
  const double c_KB = 1.3806504e-23;    // Boltzmann's constant - J K^-1
  const double c_Na = 6.0221415e+23;    // Avogadro's Number - mol^-1

  const double c_E = 2.7182818284590452354;  // e (base of natural log)

  const double c_Pi = 3.141592653589793238463;     // pi
  const double c_iPi = 0.3183098861837906715378;   // 1.0/pi
  const double c_i2Pi = 1.570796326794896619231;   // pi/2
  const double c_i3Pi = 1.047197551196597746154;   // pi/3
  const double c_i4Pi = 0.7853981633974483096157;  // pi/4
  const double c_2Pi = 6.283185307179586476925;    // 2*pi
  const double c_3Pi = 9.424777960769379715388;    // 3*pi
  const double c_4Pi = 12.56637061435917295385;    // 4*pi
  const double c_Pi2 = 9.869604401089358618834;    // pi^2
  const double c_Pi3 = 31.00627668029982017548;    // pi^3
  const double c_Pi4 = 97.4090910340024372364;     // pi^4
  const double c_Pii2 = 1.772453850905516027298;   // pi^(1/2)
  const double c_Pii3 = 1.46459188756152326302;    // pi^(1/3)
  const double c_Pii4 = 1.331335363800389712798;   // pi^(1/4)

  const double c_2i2 = 1.414213562373095048802;  // 2^(1/2)
  const double c_3i2 = 1.732050807568877293527;  // 3^(1/2)
  const double c_5i2 = 2.236067977499789696409;  // 5^(1/2)
  const double c_7i2 = 2.645751311064590590502;  // 7^(1/2)

  const double c_hwhm_2_sigma = 0.84932180028801907;    // hwhm to sigma 1/(sqrt(2*log(2)))
  const double c_fwhm_2_sigma = 0.42466090014400953;    // fwhm to sigma 1/(2*sqrt(2*log(2)))
  const double c_iehwgd_2_sigma = 0.70710678118654746;  // iehwgd to sigma 1/sqrt(2)

  const double c_mrad_2_rad = 1.0e-03;                   // mrad-->rad
  const double c_deg_2_rad = 0.01745329251994329576924;  // degrees-->rad
  const double c_mm_2_Angs = 1.0e+07;                    // mm-->Angstrom
  const double c_eV_2_keV = 1e-03;                       // ev-->keV

  const int c_cSynCPU = 5;

  const int c_nAtomsTypes = 103;
  const int c_nAtomsIons = 15;
  const int c_nqz = 128;
  const int c_nR = 128;

  const int c_thrnxny = 16;
  const int c_thrnxy = 256;
  const double c_Vrl = 0.015;

  const int cSizeofI = sizeof(int);
  const int cSizeofRD = sizeof(double);
  const int cSizeofRF = sizeof(float);
  const int cSizeofCD = 2 * cSizeofRD;

  /******************************modify vector******************************/
  enum eInput_Atoms {
    eIA_yes = 1,
    eIA_no = 2
  };

  /******************************modify vector******************************/
  enum eModify_Vector {
    eMV_yes = 1,
    eMV_no = 2
  };

  /******************************e_device type******************************/
  enum eDevice {
    e_host = 1,
    e_device = 2,
    e_host_device = 3
  };

  /******************************Slice memory type******************************/
  enum eSlice_Memory_Type {
    eSMT_Transmission = 1,
    eSMT_Potential = 2,
    eSMT_none = 3
  };

  /******************************Microscope effects*****************************/
  enum eIllumination_Model {
    eIM_Coherent = 1,
    eIM_Partial_Coherent = 2,
    eIM_Trans_Cross_Coef = 3,
    eIM_Full_Integration = 4,
    eIM_none = 5
  };

  /******************************Spatial and temporal***************************/
  enum eTemporal_Spatial_Incoh {
    eTSI_Temporal_Spatial = 1,
    eTSI_Temporal = 2,
    eTSI_Spatial = 3,
    eTSI_none = 4
  };

  /********************************MULTEM type**********************************/
  enum ePrecision {
    eP_float = 1,
    eP_double = 2
  };

  /*************************************data type******************************/
  enum eData_Type {
    eDT_float = 1,
    eDT_double = 2,
    eDT_cfloat = 3,
    eDT_cdouble = 4
  };

  /*****************************Show Data Type**********************************/
  enum eShow_CData {
    eSCD_CReal = 1,
    eSCD_CImag = 2,
    eSCD_CMod = 3,
    eSCD_CPhase = 4
  };

  /*********************************Operation mode******************************/
  enum eOperation_Mode {
    eOM_Normal = 1,
    eOM_Advanced = 2
  };

  /****************************lens variable type******************************/
  enum eLens_Var_Type {
    eLVT_off = 0,
    eLVT_m = 1,
    eLVT_f = 2,
    eLVT_Cs3 = 3,
    eLVT_Cs5 = 4,
    eLVT_mfa2 = 5,
    eLVT_afa2 = 6,
    eLVT_mfa3 = 7,
    eLVT_afa3 = 8,
    eLVT_inner_aper_ang = 9,
    eLVT_outer_aper_ang = 10
  };

  /*****************************simulation type********************************/
  enum eTEM_Sim_Type {
    eTEMST_STEM = 11,
    eTEMST_ISTEM = 12,
    eTEMST_CBED = 21,
    eTEMST_CBEI = 22,
    eTEMST_ED = 31,
    eTEMST_HRTEM = 32,
    eTEMST_PED = 41,
    eTEMST_HCTEM = 42,
    eTEMST_EWFS = 51,
    eTEMST_EWRS = 52,
    eTEMST_EELS = 61,
    eTEMST_EFTEM = 62,
    eTEMST_IWFS = 71,
    eTEMST_IWRS = 72,
    eTEMST_PPFS = 81,
    eTEMST_PPRS = 82,  // projected potential
    eTEMST_TFFS = 91,
    eTEMST_TFRS = 92,  // transmission function
    eTEMST_PropFS = 101,
    eTEMST_PropRS = 102  // propagate
  };

  /*************************simulation data output*****************************/
  enum eTEM_Output_Type {
    eTEMOT_image_tot_coh = 1,
    eTEMOT_image_tot = 2,
    eTEMOT_m2psi_tot_coh = 3,
    eTEMOT_m2psi_tot = 4,
    eTEMOT_m2psi_tot_psi_coh = 5,
    eTEMOT_psi_coh = 6,
    eTEMOT_psi_0 = 7,
    eTEMOT_V = 8,
    eTEMOT_trans = 9
  };

  /******************Electron specimen interaction model**********************/
  enum eElec_Spec_Int_Model {
    eESIM_Multislice = 1,
    eESIM_Phase_Object = 2,
    eESIM_Weak_Phase_Object = 3
  };

  /*****************************Frozen lattice model**************************/
  enum ePhonon_Model {
    ePM_Still_Atom = 1,
    ePM_Absorptive_Model = 2,
    ePM_Frozen_Phonon = 3
  };

  /*******************************Extract data********************************/
  enum ePhonon_Model_Output {
    eFMO_Total = 1,
    eFMO_Coherent = 2
  };

  /*********************Projected_Potential Slicing Type**********************/
  enum ePotential_Slicing {
    ePS_Planes = 1,
    ePS_dz_Proj = 2,
    ePS_dz_Sub = 3,
    ePS_Auto = 4
  };

  /********************Projected_Potential parameterization******************/
  enum ePotential_Type {
    ePT_Doyle_0_4 = 1,
    ePT_Peng_0_4 = 2,
    ePT_Peng_0_12 = 3,
    ePT_Kirkland_0_12 = 4,
    ePT_Weickenmeier_0_12 = 5,
    ePT_Lobato_0_12 = 6,
    ePT_none = 0
  };

  /***************************Incident Wave Type******************************/
  enum eIncident_Wave_Type {
    eIWT_Plane_Wave = 1,
    eIWT_Convergent_Wave = 2,
    eIWT_User_Define_Wave = 3,
    eIWT_Auto = 4
  };

  enum eRot_Point_Type {
    eRPT_geometric_center = 1,
    eRPT_User_Define = 2
  };

  /*****************************Real or Fourier space**************************/
  enum eSpace {
    eS_Real = 1,
    eS_Reciprocal = 2
  };

  /****************************Defocus plane type*****************************/
  enum eMatch_Border {
    eMB_Min = 1,
    eMB_Max = 2,
    eMB_MinMax = 3
  };

  /****************************Amorphous layer Type***************************/
  enum eAmorp_Lay_Type {
    eALT_Top = 1,
    eALT_Bottom = 2,
    eALT_Middle = 3,
    eALT_none = 4
  };

  /****************************Defocus plane type*****************************/
  enum eZero_Defocus_Type {
    eZDT_First = 1,
    eZDT_Middle = 2,
    eZDT_Last = 3,
    eZDT_User_Define = 4
  };

  /*******************************thick Type*********************************/
  enum eThick_Type {
    eTT_Whole_Spec = 1,
    eTT_Through_Thick = 2,
    eTT_Through_Slices = 3
  };

  /******************************Scanning Type*******************************/
  enum eScanning_Type {
    eST_Line = 1,
    eST_Area = 2
  };
  /******************************grid_2d Type*******************************/
  enum eGrid_Type {
    eGT_Regular = 1,
    eGT_Quadratic = 2
  };

  /******************************Detector type*******************************/
  enum eDetector_Type {
    eDT_Circular = 1,
    eDT_Radial = 2,
    eDT_Matrix = 3
  };

  /********************************Channelling type*************************/
  enum eChannelling_Type {
    eCT_Single_Channelling = 1,
    eCT_Mixed_Channelling = 2,
    eCT_Double_Channelling = 3
  };

  /*******************************Output type*******************************/
  enum eOutput_Type {
    eOT_Matlab = 1,
    eOT_Vector = 2
  };

  /****************************Data selection type**************************/
  enum eDat_Sel_Type {
    eDST_Closest = 1,
    eDST_Less_Than = 2,
    eDST_Greater_Than = 3
  };

  /**************************structuring element****************************/
  enum eStr_Ele {
    eSE_Disk = 1,
    eSE_Square = 2
  };

  /******************************operation**********************************/
  enum eOP {
    eOP_N = 1,
    eOP_T = 2,
    eOP_C = 3
  };

  /********************************Exec type********************************/
  enum eET {
    eET_Matrix = 1,
    eET_Vector = 2
  };

}  // namespace mt

#endif
