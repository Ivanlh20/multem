/*
 * This file is part of MULTEM.
 * Copyright 2016 Ivan Lobato <Ivanlh20@gmail.com>
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

#include <algorithm>

#include "math.cuh"
#include "types.cuh"
#include "matlab_types.cuh"
#include "traits.cuh"
#include "stream.cuh"
#include "atomic_data.hpp"
#include "host_device_functions.cuh"
#include "input_multislice.cuh"
#include "output_multislice.hpp"
#include "atom_data.hpp"
#include "atomic_cross_section.cuh"

#include <mex.h>
#include "matlab_mex.cuh"

using mt::Vector;
using mt::rmatrix_r;
using mt::rmatrix_c;
using mt::e_host;

class Output_Cross_section
{
	public:
		Output_Cross_section(): nr(0){}

		template<class TInput_Multislice>
		void set_input_data(TInput_Multislice *input_multislice)
		{ 
			nr = input_multislice->scanning.ns;

			stream.resize(input_multislice->cpu_nthread);
		}

		void init()
		{ 
			mt::fill(stream, r, 0);
			mt::fill(stream, fr, 0);
		}

		int nr;
		rmatrix_r r;
		rmatrix_r fr;

		mt::Stream<e_host> stream;
};

template<class TInput_Multislice>
void read_input_multislice(const mxArray *mx_input_multislice, TInput_Multislice &input_multislice)
{
	using value_type_r = mt::Value_type<TInput_Multislice>;

	input_multislice.precision = mt::eP_float;
	input_multislice.device = mt::e_device; 
	input_multislice.cpu_ncores = 1; 
	input_multislice.cpu_nthread = 4; 
	input_multislice.gpu_device = 0;
	input_multislice.gpu_nstream = 1;
	input_multislice.set_device();

	input_multislice.simulation_type = mt::eTEMST_STEM ;
	input_multislice.phonon_model = mt::ePM_Still_Atom;
	input_multislice.interaction_model = mt::eESIM_Multislice;
	input_multislice.potential_slicing = mt::ePS_dz_Sub;
	input_multislice.potential_type = mt::ePT_Lobato_0_12;

	input_multislice.E_0 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "E_0");
	int Z = mx_get_scalar_field<int>(mx_input_multislice, "Z");
	double rms3d = mx_get_scalar_field<double>(mx_input_multislice, "rms3d");
	double fwsig = mx_get_scalar_field<double>(mx_input_multislice, "fwhm")*mt::c_fwhm2sigma;

	bool bwl = false;
	bool pbc_xy = false;

	int nx = 1024;
	int ny = 1024;
	double lx = 20;
	double ly = 20;
	double dz = 0.4; 				

	/******************************** set atom *********************************/
	int natoms = 1;
	double atoms[6];
	atoms[0] = Z; 
	atoms[1] = 0.5*lx; 
	atoms[2] = 0.5*ly; 
	atoms[3] = 0.0; 
	atoms[4] = sqrt(rms3d*rms3d+fwsig*fwsig); 
	atoms[5] = 1.0;
	input_multislice.atoms.set_Atoms(natoms, atoms, lx, ly);
	input_multislice.grid.set_input_data(nx, ny, lx, ly, dz, bwl, pbc_xy);

	/****************************** Condenser lens ********************************/
	input_multislice.cond_lens.m = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_m"); 											// momentum of the vortex
	input_multislice.cond_lens.f = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_f"); 									// defocus(Angstrom)
	input_multislice.cond_lens.Cs3 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_Cs3")*mt::c_mm_2_Ags; 			// third order spherical aberration(mm-->Angstrom)
	input_multislice.cond_lens.Cs5 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_Cs5")*mt::c_mm_2_Ags; 			// fifth order aberration(mm-->Angstrom)
	input_multislice.cond_lens.mfa2 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_mfa2"); 							// magnitude 2-fold astigmatism(Angstrom)
	input_multislice.cond_lens.afa2 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_afa2")*mt::c_deg_2_rad; 		// angle 2-fold astigmatism(degrees-->rad)
	input_multislice.cond_lens.mfa3 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_mfa3"); 							// magnitude 3-fold astigmatism(Angstrom)
	input_multislice.cond_lens.afa3 = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_afa3")*mt::c_deg_2_rad; 		// angle 3-fold astigmatism(degrees-->rad)
	input_multislice.cond_lens.inner_aper_ang = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_inner_aper_ang")*mt::c_mrad_2_rad; 		// inner aperture(mrad-->rad)
	input_multislice.cond_lens.outer_aper_ang = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_outer_aper_ang")*mt::c_mrad_2_rad; 		// outer aperture(mrad-->rad)
	input_multislice.cond_lens.sf = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_sf"); 								// defocus spread(Angstrom)
	input_multislice.cond_lens.nsf = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_nsf"); 										// Number of integration steps for the defocus Spread
	input_multislice.cond_lens.beta = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_beta")*mt::c_mrad_2_rad; 		// divergence semi-angle(mrad-->rad)
	input_multislice.cond_lens.nbeta = mx_get_scalar_field<int>(mx_input_multislice, "cond_lens_nbeta");									// Number of integration steps for the divergence semi-angle
	input_multislice.cond_lens.zero_defocus_type = mx_get_scalar_field<mt::eZero_Defocus_Type>(mx_input_multislice, "cond_lens_zero_defocus_type");
	input_multislice.cond_lens.zero_defocus_plane = mx_get_scalar_field<value_type_r>(mx_input_multislice, "cond_lens_zero_defocus_plane");	
	input_multislice.cond_lens.set_input_data(input_multislice.E_0, input_multislice.grid);

	/********************************* Detectors ********************************/
	value_type_r lambda = mt::get_lambda(input_multislice.E_0);
	mxArray *mx_detector = mxGetField(mx_input_multislice, 0, "detector");
	input_multislice.detector.type = mt::eDT_Circular;
	mx_detector = mxGetField(mx_detector, 0, "cir");
	input_multislice.detector.resize(1);
	auto inner_ang = mx_get_scalar_field<value_type_r>(mx_detector, 0, "inner_ang")*mt::c_mrad_2_rad;
	input_multislice.detector.g_inner[0] = sin(inner_ang)/lambda;
	auto outer_ang = mx_get_scalar_field<value_type_r>(mx_detector, 0, "outer_ang")*mt::c_mrad_2_rad;
	input_multislice.detector.g_outer[0] = sin(outer_ang)/lambda;

	/********************************* Scanning ********************************/
	mt::Atom_Cal<double> atom_cal;
	mt::Atomic_Data atomic_data;
	atomic_data.Load_Data(input_multislice.potential_type);
	mt::Atom_Type<double, mt::e_host> atom_type;

	atomic_data.To_atom_type_CPU(Z, mt::c_Vrl, mt::c_nR, 0.0, atom_type);
	atom_cal.Set_Atom_Type(input_multislice.potential_type, &atom_type);
	auto rmax = atom_cal.AtomicRadius_Cutoff(3, 0.005);

	input_multislice.scanning.type = mt::eST_Line;
	input_multislice.scanning.grid_type = mt::eGT_Regular;
	input_multislice.scanning.ns = mt::c_nR;
	input_multislice.scanning.x0 = 0.5*lx;
	input_multislice.scanning.y0 = 0.5*ly;
	input_multislice.scanning.xe = 0.5*lx;
	input_multislice.scanning.ye = 0.5*ly+rmax;
	input_multislice.scanning.set_grid();

	input_multislice.validate_parameters();
 }

template<class TOutput_Cross_Section>
void set_output_cross_section(const mxArray *mx_input_cross_section, mxArray *&mx_output_cross_section, TOutput_Cross_Section &output_cross_section)
{
	mt::Input_Multislice<double> input_multislice;
	read_input_multislice(mx_input_cross_section, input_multislice);
	output_cross_section.set_input_data(&input_multislice);

	const char *field_names_output_multislice[] = {"r", "fr"};
	int number_of_fields_output_multislice = 2;
	mwSize dims_output_multislice[2] = {1, 1};

	mx_output_cross_section = mxCreateStructArray(2, dims_output_multislice, number_of_fields_output_multislice, field_names_output_multislice);
	output_cross_section.r = mx_create_matrix_field<rmatrix_r>(mx_output_cross_section, "r", 1, output_cross_section.nr);
	output_cross_section.fr = mx_create_matrix_field<rmatrix_r>(mx_output_cross_section, "fr", 1, output_cross_section.nr);
}

template<class T, mt::eDevice dev, class TOutput_Cross_Section>
void get_cross_section(const mxArray *mxB, TOutput_Cross_Section &output_cross_section)
{
	/**************************multislice calculation*******************************/
	mt::Input_Multislice<T> input_multislice;
	read_input_multislice(mxB, input_multislice);

	mt::Atomic_Cross_Section<T, dev> atomic_cross_section;
	atomic_cross_section.set_input_data(&input_multislice);

	atomic_cross_section.get(output_cross_section.r, output_cross_section.fr);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	Output_Cross_section output_cross_section;
	set_output_cross_section(prhs[0], plhs[0], output_cross_section);

	get_cross_section<float, mt::e_device>(prhs[0], output_cross_section);
}