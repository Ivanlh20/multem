/*
* This file is part of Multem.
* Copyright 2022 Ivan Lobato <Ivanlh20@gmail.com>
*
* Multem is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version of the License, or
* (at your option) any later version.
*
* Multem is distributed in the hope that it will be useful, 
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Multem. If not, see <http:// www.gnu.org/licenses/>.
*/

#pragma once

#include <algorithm>
#include <type_traits>
#include <cmath>

#include "types_mt.cuh"
#include "type_traits_gen.h"
#include "r_2d.h"
#include "r_3d.h"
#include "rot_in_parm.hpp"
#include "atomic_vib.hpp"
#include "beam_pos_2d.hpp"
#include "scan_pat.hpp"
#include "spec_slic_in_parm.hpp"
#include "lens.cuh"
#include "particles.cuh"

#include <mex.h>
#include "matlab_mex.h"

/************************ read beam positions *************************/
template <class T>
void mex_read_beam_pos(const mxArray* mex_array, mt::Beam_Pos_2d<T>& beam_pos)
{
	auto pbeam_pos = mex_get_pvctr_from_field<pMLD>(mex_array, "beam_pos");
	beam_pos.set_in_data(pbeam_pos);
}

/************************ read user define wave ************************/
template <class T>
void mex_read_user_define_wave(const mxArray* mex_array, mt::Vctr_cpu<T>& psi)
{
	auto ppsi = mex_get_pvctr_from_field<T>(mex_array, "iw_psi");
	psi.assign(ppsi);
}

/*********************** read specimen thickness ***********************/
template <class T>
void mex_read_spec_thick(const mxArray* mex_array, mt::Vctr_cpu<T>& thick)
{
	auto pthick = mex_get_pvctr_from_field<T>(mex_array, "thick");
	thick.assign(pthick);
}

/************************** read output area ***************************/
void mex_read_output_area(const mxArray* mex_array, mt::iThread_Rect_2d& output_area)
{
	const auto ip_0 = mex_get_vctr_from_field<dt_int32>(mex_array, "output_area_ip_0");
	const auto ip_e = mex_get_vctr_from_field<dt_int32>(mex_array, "output_area_ip_e");

	output_area.ix_0 = static_cast<dt_int32>(ip_0[0])-1;
	output_area.iy_0 = static_cast<dt_int32>(ip_0[1])-1;
	output_area.ix_e = static_cast<dt_int32>(ip_e[0])-1;
	output_area.iy_e = static_cast<dt_int32>(ip_e[1])-1;
}

/**********************  atomic vibration model ************************/
void mex_read_atomic_vib(const mxArray* mex_array, mt::Atomic_Vib& atomic_vib)
{
	atomic_vib.model = mex_get_enum_from_field<mt::eAtomic_Vib_Mod>(mex_array, "atomic_vib_mod");
	atomic_vib.coh_contrib = mex_get_bool_from_field(mex_array, "atomic_vib_coh_contrib");
	atomic_vib.sgl_conf = mex_get_bool_from_field(mex_array, "atomic_vib_sgl_conf");
	atomic_vib.nconf = mex_get_num_from_field<dt_int32>(mex_array, "atomic_vib_nconf");
	atomic_vib.dim = mex_get_r_3d_from_field<dt_bool>(mex_array, "atomic_vib_dim", mt::R_3d<dt_bool>(1, 1, 0));
	atomic_vib.seed = mex_get_num_from_field<dt_int32>(mex_array, "atomic_vib_seed");

	atomic_vib.set_dep_var();
}

/****************************** read atoms *****************************/
template <class T>
void mex_read_atoms(const mxArray* mex_array, mt::R_3d<T> bs, dt_bool pbc_xy, dt_bool b_statistic, mt::Ptc_Atom<T>& atoms)
{
	auto patoms = mex_get_pvctr_from_field<T>(mex_array, "spec_atoms");

	atoms.set_ptc(patoms, bs, pbc_xy, b_statistic);
}

//template <class T>
//void mex_read_atoms(const mxArray* mex_array, T bs_x, T bs_y, T bs_z, T sli_thick, mt::Ptc_Atom<T>& atoms)
//{
//	auto patoms = mex_get_pvctr_from_field<T>(mex_array, "spec_atoms");
//
//	auto ct_na = mex_get_num_from_field<dt_int32>(mex_array, "spec_cryst_na");
//	auto ct_nb = mex_get_num_from_field<dt_int32>(mex_array, "spec_cryst_nb");
//	auto ct_nc = mex_get_num_from_field<dt_int32>(mex_array, "spec_cryst_nc");
//	auto ct_a = mex_get_num_from_field<T>(mex_array, "spec_cryst_a");
//	auto ct_b = mex_get_num_from_field<T>(mex_array, "spec_cryst_b");
//	auto ct_c = mex_get_num_from_field<T>(mex_array, "spec_cryst_c");
//	auto ct_x0 = mex_get_num_from_field<T>(mex_array, "spec_cryst_x0");
//	auto ct_y0 = mex_get_num_from_field<T>(mex_array, "spec_cryst_y0");
//
//	auto mex_spec_amorp = mxGetField(mex_array, 0, "spec_amorp");
//	mt::Vctr_Spec_Lay_Info<T> spec_lay_info(mxGetN(mex_spec_amorp));
//	for(auto ik = 0; ik < spec_lay_info.size(); ik++)
//	{
//		spec_lay_info[ik].z_0 = mex_get_num_from_field<T>(mex_spec_amorp, ik, "z_0");
//		spec_lay_info[ik].z_e = mex_get_num_from_field<T>(mex_spec_amorp, ik, "z_e");
//		spec_lay_info[ik].sli_thick = mex_get_num_from_field<T>(mex_spec_amorp, ik, "sli_thick");
//	}
//
//	atoms.set_xtl_parameters(ct_na, ct_nb, ct_nc, ct_a, ct_b, ct_c, ct_x0, ct_y0);
//	atoms.set_amorphous_parameters(spec_lay_info);
//	atoms.set_ptc(patoms, bs, sli_thick);
//}

/*********************** read rotation parameters *********************/
template <class T>
void mex_read_rot_in_parm(const mxArray* mex_array, mt::Rot_In_Parm<T>& rot_in_parm)
{
	rot_in_parm.theta = mex_get_num_from_field<T>(mex_array, "spec_rot_theta")*mt::c_deg_2_rad<T>;
	rot_in_parm.u_0 = mex_get_r_3d_from_field<T>(mex_array, "spec_rot_u_0", mt::R_3d<T>(0, 0, 1));
	rot_in_parm.ctr_type = mex_get_enum_from_field<mt::eRot_Ctr_Typ>(mex_array, "spec_rot_ctr_typ");
	rot_in_parm.ctr_p = mex_get_r_3d_from_field<T>(mex_array, "spec_rot_ctr_p", mt::R_3d<T>(0, 0, 0));

	rot_in_parm.set_dep_var();
}
	
/*************************** specimen slicing *************************/
template <class T>
void mex_read_spec_sli_in_parm(const mxArray* mex_array, mt::Vctr_Spec_Slic_In_Parm<T>& spec_slic_in_parm)
{
	auto mex_spec_slic = mxGetField(mex_array, 0, "spec_slic");

	spec_slic_in_parm.resize(mex_get_s0(mex_spec_slic));
	for(auto ik = 0; ik < spec_slic_in_parm.size(); ik++)
	{
		spec_slic_in_parm[ik].typ = mex_get_enum_from_field<mt::eSpec_Slic_Typ>(mex_spec_slic, ik, "typ");
		spec_slic_in_parm[ik].sli_thick = mex_get_num_from_field<T>(mex_spec_slic, ik, "sli_thick");
		spec_slic_in_parm[ik].sel_typ = mex_get_enum_from_field<mt::eSpec_Slic_Sel_Typ>(mex_spec_slic, ik, "sel_typ");
		spec_slic_in_parm[ik].sel_tag = mex_get_num_from_field<dt_int32>(mex_spec_slic, ik, "sel_tag");
		spec_slic_in_parm[ik].sel_Z = mex_get_num_from_field<dt_int32>(mex_spec_slic, ik, "sel_Z");
		spec_slic_in_parm[ik].sel_z_lim = mex_get_r_2d_from_field<T>(mex_spec_slic, ik, "sel_z_lim");
		auto pz_plns = mex_get_pvctr_from_field<T>(mex_array, ik, "z_plns");
		spec_slic_in_parm[ik].z_plns.assign(pz_plns);
	}
}

/************************ read condenser lens *************************/
template <class T>
void mex_read_cond_lens(const mxArray* mex_array, mt::Lens<T>& cond_lens)
{
	cond_lens.m = mex_get_num_from_field<dt_int32>(mex_array, "cond_lens_m");		// momentum of the vortex
	cond_lens.c_10 = mex_get_num_from_field<T>(mex_array, "cond_lens_c_10");		// defocus (Angstrom)
	cond_lens.c_12 = mex_get_num_from_field<T>(mex_array, "cond_lens_c_12");		// 2-fold astigmatism (Angstrom)
	cond_lens.phi_12 = mex_get_num_from_field<T>(mex_array, "cond_lens_phi_12")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 2-fold astigmatism (degrees-->rad)

	cond_lens.c_21 = mex_get_num_from_field<T>(mex_array, "cond_lens_c_21");		// Axial coma (Angstrom)
	cond_lens.phi_21 = mex_get_num_from_field<T>(mex_array, "cond_lens_phi_21")*mt::c_deg_2_rad<T>;		// Azimuthal angle of axial coma (degrees-->rad)
	cond_lens.c_23 = mex_get_num_from_field<T>(mex_array, "cond_lens_c_23");		// 3-fold astigmatism (Angstrom)
	cond_lens.phi_23 = mex_get_num_from_field<T>(mex_array, "cond_lens_phi_23")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 3-fold astigmatism (degrees-->rad)

	cond_lens.c_30 = mex_get_num_from_field<T>(mex_array, "cond_lens_c_30")*mt::c_mm_2_angs<T>;		// 3rd order spherical aberration (mm-->Angstrom)
	cond_lens.c_32 = mex_get_num_from_field<T>(mex_array, "cond_lens_c_32");		// Axial star aberration (Angstrom)
	cond_lens.phi_32 = mex_get_num_from_field<T>(mex_array, "cond_lens_phi_32")*mt::c_deg_2_rad<T>;		// Azimuthal angle of axial star aberration (degrees-->rad)
	cond_lens.c_34 = mex_get_num_from_field<T>(mex_array, "cond_lens_c_34");		// 4-fold astigmatism (Angstrom)
	cond_lens.phi_34 = mex_get_num_from_field<T>(mex_array, "cond_lens_phi_34")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 4-fold astigmatism (degrees-->rad)

	cond_lens.c_41 = mex_get_num_from_field<T>(mex_array, "cond_lens_c_41");		// 4th order axial coma (Angstrom)
	cond_lens.phi_41 = mex_get_num_from_field<T>(mex_array, "cond_lens_phi_41")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 4th order axial coma (degrees-->rad)
	cond_lens.c_43 = mex_get_num_from_field<T>(mex_array, "cond_lens_c_43");		// 3-lobe aberration (Angstrom)
	cond_lens.phi_43 = mex_get_num_from_field<T>(mex_array, "cond_lens_phi_43")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 3-lobe aberration (degrees-->rad)
	cond_lens.c_45 = mex_get_num_from_field<T>(mex_array, "cond_lens_c_45");		// 5-fold astigmatism (Angstrom)
	cond_lens.phi_45 = mex_get_num_from_field<T>(mex_array, "cond_lens_phi_45")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 5-fold astigmatism (degrees-->rad)

	cond_lens.c_50 = mex_get_num_from_field<T>(mex_array, "cond_lens_c_50")*mt::c_mm_2_angs<T>;		// 5th order spherical aberration (mm-->Angstrom)
	cond_lens.c_52 = mex_get_num_from_field<T>(mex_array, "cond_lens_c_52");		// 5th order axial star aberration (Angstrom)
	cond_lens.phi_52 = mex_get_num_from_field<T>(mex_array, "cond_lens_phi_52")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 5th order axial star aberration (degrees-->rad)
	cond_lens.c_54 = mex_get_num_from_field<T>(mex_array, "cond_lens_c_54");		// 5th order rosette aberration (Angstrom)
	cond_lens.phi_54 = mex_get_num_from_field<T>(mex_array, "cond_lens_phi_54")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 5th order rosette aberration (degrees-->rad)
	cond_lens.c_56 = mex_get_num_from_field<T>(mex_array, "cond_lens_c_56");		// 6-fold astigmatism (Angstrom)
	cond_lens.phi_56 = mex_get_num_from_field<T>(mex_array, "cond_lens_phi_56")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 6-fold astigmatism (degrees-->rad)

	cond_lens.inner_aper_ang = mex_get_num_from_field<T>(mex_array, "cond_lens_inner_aper_ang")*mt::c_mrad_2_rad<T>;		// inner aperture (mrad-->rad)
	cond_lens.outer_aper_ang = mex_get_num_from_field<T>(mex_array, "cond_lens_outer_aper_ang")*mt::c_mrad_2_rad<T>;		// outer aperture (mrad-->rad)

	/********************* defocus spread function ********************/
	cond_lens.tp_inc_a = mex_get_num_from_field<T>(mex_array, "cond_lens_tp_inc_a");				// Height proportion of a normalized Gaussian [0, 1]
	cond_lens.tp_inc_sigma = mex_get_num_from_field<T>(mex_array, "cond_lens_tp_inc_sigma");		// standard deviation of the source spread function for the Gaussian component: For parallel ilumination(Å^-1);otherwise (Å)
	cond_lens.tp_inc_beta = mex_get_num_from_field<T>(mex_array, "cond_lens_tp_inc_beta");			// standard deviation of the source spread function for the exponential component: For parallel ilumination(Å^-1);otherwise (Å)
	cond_lens.tp_inc_npts = mex_get_num_from_field<dt_int32>(mex_array, "cond_lens_tp_inc_npts");	// number of integration points

	/***************** source size broadening function *****************/
	cond_lens.spt_inc_a = mex_get_num_from_field<T>(mex_array, "cond_lens_spt_inc_a");							// Height proportion of a normalized Gaussian [0, 1]
	cond_lens.spt_inc_sigma = mex_get_num_from_field<T>(mex_array, "cond_lens_spt_inc_sigma");					// standard deviation of the source spread function for the Gaussian component: For parallel ilumination(Å^-1);otherwise (Å)
	cond_lens.spt_inc_beta = mex_get_num_from_field<T>(mex_array, "cond_lens_spt_inc_beta");					// standard deviation of the source spread function for the exponential component: For parallel ilumination(Å^-1);otherwise (Å)
	cond_lens.spt_inc_rad_npts = mex_get_num_from_field<dt_int32>(mex_array, "cond_lens_spt_inc_rad_npts");		// number of radial integration points
	cond_lens.spt_inc_azm_npts = mex_get_num_from_field<dt_int32>(mex_array, "cond_lens_spt_inc_azm_npts");		// number of azimuth integration points

	/********************* zero defocus reference ********************/																																				/********************* zero defocus reference ********************/
	cond_lens.zero_def_typ = mex_get_enum_from_field<mt::eZero_Def_Typ>(mex_array, "cond_lens_zero_def_typ");		// Zero defocus type
	cond_lens.zero_def_plane = mex_get_num_from_field<T>(mex_array, "cond_lens_zero_def_plane");					// Zero defocus position
}

/************************* read objective lens ************************/
template <class T>
void mex_read_obj_lens(const mxArray* mex_array, mt::Lens<T>& obj_lens)
{
	obj_lens.m = mex_get_num_from_field<dt_int32>(mex_array, "obj_lens_m");		// momentum of the vortex
	obj_lens.c_10 = mex_get_num_from_field<T>(mex_array, "obj_lens_c_10");		// defocus (Angstrom)
	obj_lens.c_12 = mex_get_num_from_field<T>(mex_array, "obj_lens_c_12");		// 2-fold astigmatism (Angstrom)
	obj_lens.phi_12 = mex_get_num_from_field<T>(mex_array, "obj_lens_phi_12")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 2-fold astigmatism (degrees-->rad)

	obj_lens.c_21 = mex_get_num_from_field<T>(mex_array, "obj_lens_c_21");		// Axial coma (Angstrom)
	obj_lens.phi_21 = mex_get_num_from_field<T>(mex_array, "obj_lens_phi_21")*mt::c_deg_2_rad<T>;		// Azimuthal angle of axial coma (degrees-->rad)
	obj_lens.c_23 = mex_get_num_from_field<T>(mex_array, "obj_lens_c_23");		// 3-fold astigmatism (Angstrom)
	obj_lens.phi_23 = mex_get_num_from_field<T>(mex_array, "obj_lens_phi_23")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 3-fold astigmatism (degrees-->rad)

	obj_lens.c_30 = mex_get_num_from_field<T>(mex_array, "obj_lens_c_30")*mt::c_mm_2_angs<T>;		// 3rd order spherical aberration (mm-->Angstrom)
	obj_lens.c_32 = mex_get_num_from_field<T>(mex_array, "obj_lens_c_32");		// Axial star aberration (Angstrom)
	obj_lens.phi_32 = mex_get_num_from_field<T>(mex_array, "obj_lens_phi_32")*mt::c_deg_2_rad<T>;		// Azimuthal angle of axial star aberration (degrees-->rad)
	obj_lens.c_34 = mex_get_num_from_field<T>(mex_array, "obj_lens_c_34");		// 4-fold astigmatism (Angstrom)
	obj_lens.phi_34 = mex_get_num_from_field<T>(mex_array, "obj_lens_phi_34")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 4-fold astigmatism (degrees-->rad)

	obj_lens.c_41 = mex_get_num_from_field<T>(mex_array, "obj_lens_c_41");		// 4th order axial coma (Angstrom)
	obj_lens.phi_41 = mex_get_num_from_field<T>(mex_array, "obj_lens_phi_41")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 4th order axial coma (degrees-->rad)
	obj_lens.c_43 = mex_get_num_from_field<T>(mex_array, "obj_lens_c_43");		// 3-lobe aberration (Angstrom)
	obj_lens.phi_43 = mex_get_num_from_field<T>(mex_array, "obj_lens_phi_43")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 3-lobe aberration (degrees-->rad)
	obj_lens.c_45 = mex_get_num_from_field<T>(mex_array, "obj_lens_c_45");		// 5-fold astigmatism (Angstrom)
	obj_lens.phi_45 = mex_get_num_from_field<T>(mex_array, "obj_lens_phi_45")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 5-fold astigmatism (degrees-->rad)

	obj_lens.c_50 = mex_get_num_from_field<T>(mex_array, "obj_lens_c_50")*mt::c_mm_2_angs<T>;		// 5th order spherical aberration (mm-->Angstrom)
	obj_lens.c_52 = mex_get_num_from_field<T>(mex_array, "obj_lens_c_52");		// 5th order axial star aberration (Angstrom)
	obj_lens.phi_52 = mex_get_num_from_field<T>(mex_array, "obj_lens_phi_52")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 5th order axial star aberration (degrees-->rad)
	obj_lens.c_54 = mex_get_num_from_field<T>(mex_array, "obj_lens_c_54");		// 5th order rosette aberration (Angstrom)
	obj_lens.phi_54 = mex_get_num_from_field<T>(mex_array, "obj_lens_phi_54")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 5th order rosette aberration(degrees-->rad)
	obj_lens.c_56 = mex_get_num_from_field<T>(mex_array, "obj_lens_c_56");		// 6-fold astigmatism (Angstrom)
	obj_lens.phi_56 = mex_get_num_from_field<T>(mex_array, "obj_lens_phi_56")*mt::c_deg_2_rad<T>;		// Azimuthal angle of 6-fold astigmatism (degrees-->rad)

	obj_lens.inner_aper_ang = mex_get_num_from_field<T>(mex_array, "obj_lens_inner_aper_ang")*mt::c_mrad_2_rad<T>;		// inner aperture (mrad-->rad)
	obj_lens.outer_aper_ang = mex_get_num_from_field<T>(mex_array, "obj_lens_outer_aper_ang")*mt::c_mrad_2_rad<T>;		// outer aperture (mrad-->rad)

	/********************* defocus spread function ********************/
	obj_lens.tp_inc_a = mex_get_num_from_field<T>(mex_array, "obj_lens_tp_inc_a");					// Height proportion of a normalized Gaussian [0, 1] 
	obj_lens.tp_inc_sigma = mex_get_num_from_field<T>(mex_array, "obj_lens_tp_inc_sigma");			// standard deviation of the source spread function for the Gaussian component: For parallel ilumination(Å^-1);otherwise (Å)
	obj_lens.tp_inc_beta = mex_get_num_from_field<T>(mex_array, "obj_lens_tp_inc_beta");			// standard deviation of the source spread function for the exponential component: For parallel ilumination(Å^-1);otherwise (Å)
	obj_lens.tp_inc_npts = mex_get_num_from_field<dt_int32>(mex_array, "obj_lens_tp_inc_npts");		// number of integration points

	/********************* source spread function *********************/
	obj_lens.spt_inc_a = mex_get_num_from_field<T>(mex_array, "cond_lens_spt_inc_a");							// Height proportion of a normalized Gaussian [0, 1] 
	obj_lens.spt_inc_sigma = mex_get_num_from_field<T>(mex_array, "cond_lens_spt_inc_sigma");					// standard deviation of the source spread function for the Gaussian component: For parallel ilumination(Å^-1);otherwise (Å)
	obj_lens.spt_inc_beta = mex_get_num_from_field<T>(mex_array, "cond_lens_spt_inc_beta");						// standard deviation of the source spread function for the exponential component: For parallel ilumination(Å^-1);otherwise (Å)
	obj_lens.spt_inc_rad_npts = mex_get_num_from_field<dt_int32>(mex_array, "cond_lens_spt_inc_rad_npts");		// number of radial integration points
	obj_lens.spt_inc_azm_npts = mex_get_num_from_field<dt_int32>(mex_array, "cond_lens_spt_inc_azm_npts");		// number of azimuth integration points

	/********************* zero defocus reference ********************/
	obj_lens.zero_def_typ = mex_get_num_from_field<mt::eZero_Def_Typ>(mex_array, "obj_lens_zero_def_typ");		// Zero defocus type
	obj_lens.zero_def_plane = mex_get_num_from_field<T>(mex_array, "obj_lens_zero_def_plane");					// Zero defocus position
}

/*************************** read scanning ****************************/
template <class T>
void mex_read_scan_pat(const mxArray* mex_array, mt::Scan_Pat<T>& scan_pat)
{
	scan_pat.typ = mex_get_enum_from_field<mt::eScan_Pat_Typ>(mex_array, "scan_pat_typ");
	scan_pat.pbc = mex_get_bool_from_field(mex_array, "scan_pat_pbc");
	scan_pat.spxs = mex_get_bool_from_field(mex_array, "scan_pat_spxs");
	scan_pat.nsp = mex_get_r_2d_rep_from_field<dt_int32>(mex_array, "scan_pat_nsp");
	scan_pat.r_0 = mex_get_r_2d_from_field<T>(mex_array, "scan_pat_r_0");
	scan_pat.r_e = mex_get_r_2d_from_field<T>(mex_array, "scan_pat_r_e");

	if (scan_pat.is_scan_pat_user_def())
	{
		scan_pat.r = mex_get_vctr_r_2d_from_field<T>(mex_array, "scan_pat_r");
	}

	scan_pat.set_dep_var();
}

/***************************** read detector ***************************/
template <class T>
void mex_read_detector(const mxArray* mex_array, T E_0, mt::Detector<T, mt::edev_cpu>& detector)
{
	T lambda = mt::fcn_lambda(E_0);
	mxArray* mex_detector = mxGetField(mex_array, 0, "detector");
	detector.type = mex_get_enum_from_field<mt::eDetector_Typ>(mex_detector, "type");

	switch (detector.type)
	{
		case mt::edt_circular:
		{
			mex_detector = mxGetField(mex_detector, 0, "cir");
			dt_int32 ndetector = mex_get_MxN(mex_detector);
			if (ndetector > 0)
			{
				detector.resize(ndetector);
				for(auto i = 0; i < detector.size(); i++)
				{
					auto inner_ang = mex_get_num_from_field<T>(mex_detector, i, "inner_ang")*mt::c_mrad_2_rad<T>;
					detector.g_inner[i] = sin(inner_ang)/lambda;
					auto outer_ang = mex_get_num_from_field<T>(mex_detector, i, "outer_ang")*mt::c_mrad_2_rad<T>;
					detector.g_outer[i] = sin(outer_ang)/lambda;
				}
			}
		}
		break;
		case mt::edt_radial:
		{
			mex_detector = mxGetField(mex_detector, 0, "radial");
			dt_int32 ndetector = mxGetN(mex_detector);
			if (ndetector > 0)
			{
				detector.resize(ndetector);
				for(auto i = 0; i < detector.size(); i++)
				{
					// auto x = mex_get_pvctr_from_field<pMLD>(mex_detector, i, "x");
					// mt::assign(x, detector.x[i]);
					// mt::fcn_scale(detector.x[i], 1.0/lambda);

					auto fx = mex_get_pvctr_from_field<pMLD>(mex_detector, i, "fx");
					mt::assign(fx, detector.fx[i]);
				}
			}
		}
		break;
		case mt::edt_matrix:
		{
			mex_detector = mxGetField(mex_detector, 0, "matrix");
			dt_int32 ndetector = mxGetN(mex_detector);
			if (ndetector > 0)
			{
				detector.resize(ndetector);
				for(auto i = 0; i < detector.size(); i++)
				{
					// auto R = mex_get_pvctr_from_field<pMLD>(mex_detector, i, "x");
					// mt::assign(R, detector.R[i]);
					// mt::fcn_scale(detector.R[i], 1.0/lambda);
					// mt::fcn_fftsft_2d(grid_2d, detector.R[i]);

					auto fR = mex_get_pvctr_from_field<pMLD>(mex_detector, i, "fR");
					mt::assign(fR, detector.fR[i]);
				}
			}
		}
		break;
	}
}