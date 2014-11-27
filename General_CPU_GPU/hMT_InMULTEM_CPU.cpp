/*
 * This file is part of MULTEM.
 * Copyright 2014 Ivan Lobato <Ivanlh20@gmail.com>
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
 * along with MULTEM. If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstring>
#include "hConstTypes.h"
#include "hMT_InMULTEM_CPU.h"

void cMT_InMULTEM_CPU::freeMemory(){
	gpu = 0;
	SimType = 0;
	nConfFP = 0;
	DimFP = 0;
	SeedFP = 0;
	PotPar = 0;
	MEffect = 0;		
	STEffect = 0;	
	ZeroDefTyp = 0;
	ZeroDefPlane = 0;
	ApproxModel = 0;
	BandwidthLimit = 0;
	ThicknessTyp = 0;
	nThickness = 0;
	delete [] Thickness; Thickness = 0;

	E0 = 0;
	theta = 0;
	phi = 0;
	nx = 0;
	ny = 0;
	lx = 0;
	ly = 0;
	dz = 0;

	MC_m = 0;
	MC_f = 0;
	MC_Cs3 = 0;
	MC_Cs5 = 0;
	MC_mfa2 = 0;
	MC_afa2 = 0;
	MC_mfa3 = 0;
	MC_afa3 = 0;
	MC_aobjl = 0;
	MC_aobju = 0;
	MC_sf = 0;
	MC_nsf = 0;
	MC_beta = 0;
	MC_nbeta = 0;

	nAtomsM = 0;
	AtomsM = 0;

	STEM_line = 0;
	STEM_FastCal = false;
	STEM_ns = 0;
	STEM_x1u = 0;
	STEM_y1u = 0;
	STEM_x2u = 0;
	STEM_y2u = 0;
	
	STEM_nDet = 0;
	delete [] STEM_DetCir; STEM_DetCir = 0;

	CBED_x0 = 0;
	CBED_y0 = 0;

	CBEI_x0 = 0;
	CBEI_y0 = 0;

	PED_nrot = 0;
	PED_theta = 0;

	HCI_nrot = 0;
	HCI_theta = 0;
}

cMT_InMULTEM_CPU::cMT_InMULTEM_CPU(){
	gpu = 0;
	SimType = 0;
	nConfFP = 0;
	DimFP = 0;
	SeedFP = 0;
	PotPar = 0;
	MEffect = 0;		
	STEffect = 0;	
	ZeroDefTyp = 0;
	ZeroDefPlane = 0;
	ApproxModel = 0;
	BandwidthLimit = 0;
	ThicknessTyp = 0;
	nThickness = 0;
	Thickness = 0;

	E0 = 0;
	theta = 0;
	phi = 0;
	nx = 0;
	ny = 0;
	lx = 0;
	ly = 0;
	dz = 0;

	MC_m = 0;
	MC_f = 0;
	MC_Cs3 = 0;
	MC_Cs5 = 0;
	MC_mfa2 = 0;
	MC_afa2 = 0;
	MC_mfa3 = 0;
	MC_afa3 = 0;
	MC_aobjl = 0;
	MC_aobju = 0;
	MC_sf = 0;
	MC_nsf = 0;
	MC_beta = 0;
	MC_nbeta = 0;

	nAtomsM = 0;
	AtomsM = 0;

	STEM_line = 0;
	STEM_FastCal = false;
	STEM_ns = 0;
	STEM_x1u = 0;
	STEM_y1u = 0;
	STEM_x2u = 0;
	STEM_y2u = 0;
	
	STEM_nDet = 0;
	STEM_DetCir = 0;

	CBED_x0 = 0;
	CBED_y0 = 0;

	CBEI_x0 = 0;
	CBEI_y0 = 0;

	PED_nrot = 0;
	PED_theta = 0;

	HCI_nrot = 0;
	HCI_theta = 0;
}

cMT_InMULTEM_CPU::~cMT_InMULTEM_CPU(){
	freeMemory();
}