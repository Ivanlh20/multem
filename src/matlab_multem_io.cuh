/*
 * This file is part of Multem.
 * Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
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

#ifndef MATLAB_MULTEM_IO_H
#define MATLAB_MULTEM_IO_H

#include <algorithm>
#include <type_traits>
#include <cmath>

#include "types.cuh"
#include "type_traits_gen.cuh"
#include "particles.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::R_2d;
using mt::R_3d;

/************************ read beam positions *************************/
template <class T>
void mex_read_beam_pos(const mxArray *mex_array, mt::Beam_Pos_2d<T>& beam_pos)
{
	auto beam_pos_r = mex_get_pvctr_from_field<pMLD>(mex_array, "beam_pos");
	beam_pos.set_in_data(beam_pos_r);
}

/************************ read user define wave ************************/
template <class T>
void mex_read_user_define_wave(const mxArray *mex_array, host_vector<complex<T>>& psi)
{
	auto psi_r = mex_get_pvctr_from_field<pMx_c>(mex_array, "iw_psi");
	mt::assign(psi_r, psi);
}

/*********************** read specimen thickness ***********************/
template <class T>
void mex_read_spec_thick(const mxArray *mex_array, host_vector<T>& thick)
{
	auto thick_r = mex_get_pvctr_from_field<pMLD>(mex_array, "thick");
	mt::assign(thick_r, thick);
}

/************************* Read output area *************************/
void mex_read_output_area(const mxArray *mex_array, mt::iThread_Rect_2d& output_area)
{
	auto ip_0 = mex_get_pvctr_from_field<pMLD>(mex_array, "output_area_ip_0");
	auto ip_e = mex_get_pvctr_from_field<pMLD>(mex_array, "output_area_ip_e");

	output_area.ix_0 = static_cast<dt_int32>(ip_0[0])-1;
	output_area.iy_0 = static_cast<dt_int32>(ip_0[1])-1;
	output_area.ix_e = static_cast<dt_int32>(ip_e[0])-1;
	output_area.iy_e = static_cast<dt_int32>(ip_e[1])-1;
}

/************** Electron-Phonon_Par interaction model **************/
void mex_read_phonon_par(const mxArray *mex_array, mt::Phonon_Par &phonon_par)
{
	phonon_par.model = mex_get_num_from_field<mt::ePhonon_Model>(mex_array, "pn_model");
	phonon_par.coh_contrib = mex_get_num_from_field<dt_bool>(mex_array, "pn_coh_contrib");
	phonon_par.single_conf = mex_get_num_from_field<dt_bool>(mex_array, "pn_float32_conf");
	phonon_par.nconf = mex_get_num_from_field<dt_int32>(mex_array, "pn_nconf");
	phonon_par.set_dim(mex_get_num_from_field<dt_int32>(mex_array, "pn_dim"));
	phonon_par.seed = mex_get_num_from_field<dt_int32>(mex_array, "pn_seed");
}

/****************************** Read atoms *****************************/
template <class T>
void mex_read_atoms(const mxArray *mex_array, T bs_x, T bs_y, T bs_z, T sli_thk, mt::Ptc_Atom<T>& atoms)
{
	auto atoms_r = mex_get_pvctr_from_field<pMLD>(mex_array, "spec_atoms");

	auto ct_na = mex_get_num_from_field<dt_int32>(mex_array, "spec_cryst_na");
	auto ct_nb = mex_get_num_from_field<dt_int32>(mex_array, "spec_cryst_nb");
	auto ct_nc = mex_get_num_from_field<dt_int32>(mex_array, "spec_cryst_nc");
	auto ct_a = mex_get_num_from_field<T>(mex_array, "spec_cryst_a");
	auto ct_b = mex_get_num_from_field<T>(mex_array, "spec_cryst_b");
	auto ct_c = mex_get_num_from_field<T>(mex_array, "spec_cryst_c");
	auto ct_x0 = mex_get_num_from_field<T>(mex_array, "spec_cryst_x0");
	auto ct_y0 = mex_get_num_from_field<T>(mex_array, "spec_cryst_y0");

	auto mex_spec_amorp = mxGetField(mex_array, 0, "spec_amorp");
	mt::Vctr<mt::Spec_Lay_Info<T>, mt::edev_cpu> spec_lay_info(mxGetN(mex_spec_amorp));
	for(auto ik = 0; ik < spec_lay_info.size(); ik++)
	{
		spec_lay_info[ik].z_0 = mex_get_num_from_field<T>(mex_spec_amorp, ik, "z_0");
		spec_lay_info[ik].z_e = mex_get_num_from_field<T>(mex_spec_amorp, ik, "z_e");
		spec_lay_info[ik].sli_thk = mex_get_num_from_field<T>(mex_spec_amorp, ik, "sli_thk");
	}

	atoms.set_xtl_parameters(ct_na, ct_nb, ct_nc, ct_a, ct_b, ct_c, ct_x0, ct_y0);
	atoms.set_amorphous_parameters(spec_lay_info);
	atoms.set_ptc(atoms_r.rows, atoms_r.cols, atoms_r.real, bs_x, bs_y, bs_z, sli_thk);
}

/************************** rotation parameter *************************/
template <class T>
void mex_read_rot_par(const mxArray *mex_array, mt::Rot_Par<T>& rot_par)
{
	rot_par.theta = mex_get_num_from_field<T>(mex_array, "spec_rot_theta")*mt::c_deg_2_rad<T>;
	rot_par.u0 = mex_get_r_3d_from_field<T>(mex_array, "spec_rot_u0", R_3d<T>(0, 0, 1));
	rot_par.center_type = mex_get_num_from_field<mt::eRot_Point_Typ>(mex_array, "spec_rot_center_type");
	rot_par.center_p = mex_get_r_3d_from_field<T>(mex_array, "spec_rot_center_p", R_3d<T>(0, 0, 0));
}

/*********************** Read Condenser lens ***************************/
template <class T>
void mex_read_cond_lens(const mxArray *mex_array, mt::Lens<T>& cond_lens)
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
	cond_lens.tp_inc_a = mex_get_num_from_field<T>(mex_array, "cond_lens_tp_inc_a");		// Height proportion of a normalized Gaussian [0, 1]
	cond_lens.tp_inc_sigma = mex_get_num_from_field<T>(mex_array, "cond_lens_tp_inc_sigma");		// standard deviation of the source spread function for the Gaussian component: For parallel ilumination(Å^-1);otherwise (Å)
	cond_lens.tp_inc_beta = mex_get_num_from_field<T>(mex_array, "cond_lens_tp_inc_beta");		// standard deviation of the source spread function for the exponential component: For parallel ilumination(Å^-1);otherwise (Å)
	cond_lens.tp_inc_npts = mex_get_num_from_field<dt_int32>(mex_array, "cond_lens_tp_inc_npts");		// number of integration points

	/***************** source size broadening function *****************/
	cond_lens.spt_inc_a = mex_get_num_from_field<T>(mex_array, "cond_lens_spt_inc_a");		// Height proportion of a normalized Gaussian [0, 1]
	cond_lens.spt_inc_sigma = mex_get_num_from_field<T>(mex_array, "cond_lens_spt_inc_sigma");		// standard deviation of the source spread function for the Gaussian component: For parallel ilumination(Å^-1);otherwise (Å)
	cond_lens.spt_inc_beta = mex_get_num_from_field<T>(mex_array, "cond_lens_spt_inc_beta");		// standard deviation of the source spread function for the exponential component: For parallel ilumination(Å^-1);otherwise (Å)
	cond_lens.spt_inc_rad_npts = mex_get_num_from_field<dt_int32>(mex_array, "cond_lens_spt_inc_rad_npts");		// number of radial integration points
	cond_lens.spt_inc_azm_npts = mex_get_num_from_field<dt_int32>(mex_array, "cond_lens_spt_inc_azm_npts");

	/********************* zero defocus reference ********************/																																				/********************* zero defocus reference ********************/
	cond_lens.zero_defocus_type = mex_get_num_from_field<mt::eZero_Defocus_Typ>(mex_array, "cond_lens_zero_defocus_type");		// Zero defocus type
	cond_lens.zero_defocus_plane = mex_get_num_from_field<T>(mex_array, "cond_lens_zero_defocus_plane");		// Zero defocus position
}

/*********************** Read Objective lens ***************************/
template <class T>
void mex_read_obj_lens(const mxArray *mex_array, mt::Lens<T>& obj_lens)
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
	obj_lens.tp_inc_a = mex_get_num_from_field<T>(mex_array, "obj_lens_tp_inc_a");		// Height proportion of a normalized Gaussian [0, 1] 
	obj_lens.tp_inc_sigma = mex_get_num_from_field<T>(mex_array, "obj_lens_tp_inc_sigma");		// standard deviation of the source spread function for the Gaussian component: For parallel ilumination(Å^-1);otherwise (Å)
	obj_lens.tp_inc_beta = mex_get_num_from_field<T>(mex_array, "obj_lens_tp_inc_beta");		// standard deviation of the source spread function for the exponential component: For parallel ilumination(Å^-1);otherwise (Å)
	obj_lens.tp_inc_npts = mex_get_num_from_field<dt_int32>(mex_array, "obj_lens_tp_inc_npts");		// number of integration points

	/********************* source spread function *********************/
	obj_lens.spt_inc_a = mex_get_num_from_field<T>(mex_array, "cond_lens_spt_inc_a");		// Height proportion of a normalized Gaussian [0, 1] 
	obj_lens.spt_inc_sigma = mex_get_num_from_field<T>(mex_array, "cond_lens_spt_inc_sigma");		// standard deviation of the source spread function for the Gaussian component: For parallel ilumination(Å^-1);otherwise (Å)
	obj_lens.spt_inc_beta = mex_get_num_from_field<T>(mex_array, "cond_lens_spt_inc_beta");		// standard deviation of the source spread function for the exponential component: For parallel ilumination(Å^-1);otherwise (Å)
	obj_lens.spt_inc_rad_npts = mex_get_num_from_field<dt_int32>(mex_array, "cond_lens_spt_inc_rad_npts");		// number of radial integration points
	obj_lens.spt_inc_azm_npts = mex_get_num_from_field<dt_int32>(mex_array, "cond_lens_spt_inc_azm_npts");		// number of azimuth integration points

	/********************* zero defocus reference ********************/
	obj_lens.zero_defocus_type = mex_get_num_from_field<mt::eZero_Defocus_Typ>(mex_array, "obj_lens_zero_defocus_type");		// Zero defocus type
	obj_lens.zero_defocus_plane = mex_get_num_from_field<T>(mex_array, "obj_lens_zero_defocus_plane");		// Zero defocus position
}

/*************************** Read scanning ****************************/
template <class T>
void mex_read_scanning(const mxArray *mex_array, mt::Scan<T>& scanning)
{
	scanning.type = mex_get_num_from_field<mt::eScan_Typ>(mex_array, "scanning_type");
	scanning.pbc = mex_get_num_from_field<dt_bool>(mex_array, "scanning_periodic");
	scanning.ns = mex_get_num_from_field<dt_int32>(mex_array, "scanning_ns");

	auto R_0 = mex_get_pvctr_from_field<pMLD>(mex_array, "scanning_R_0");
	auto R_e = mex_get_pvctr_from_field<pMLD>(mex_array, "scanning_R_e");

 	scanning.R_0 = mt::R_2d<T>(R_0);
	scanning.R_e = mt::R_2d<T>(R_e);
}

/*************************** Read detector ****************************/
template <class T>
void mex_read_detector(const mxArray *mex_array, T E_0, mt::Detector<T, mt::edev_cpu>& detector)
{
	T lambda = mt::fcn_lambda(E_0);
	mxArray *mex_detector = mxGetField(mex_array, 0, "detector");
	detector.type = mex_get_num_from_field<mt::eDetector_Typ>(mex_detector, "type");

	switch (detector.type)
	{
		case mt::edt_Circular:
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
		case mt::edt_Radial:
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
		case mt::edt_Matrix:
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


#endif