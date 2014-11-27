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
#include "hMT_MGP_CPU.h"

void cMT_MGP_CPU::freeMemory(){
	SimType = 0;
	MulOrder = 0;
	nConfFP = 0;
	DimFP = 0;
	DistFP = 0;
	SeedFP = 0;
	PotPar = 0;
	MEffect = 0;
	STEffect = 0;
	ZeroDefTyp = 0;
	ZeroDefPlane = 0;
	ApproxModel = 0;
	ThkTyp = 0;
	nThk = 0;
	delete [] Thk; Thk = 0;
	BWL = false;
	PBC_xy = false;
	Vrl = 0;
	E0 = 0;
	theta = 0;	
	phi = 0;
	lx = 0;
	ly = 0;	
	dz = 0;
	nx = 0;
	ny = 0;
}

cMT_MGP_CPU::cMT_MGP_CPU(){
	gpu = 0;
	SimType = 52;
	MulOrder = 2;
	nConfFP = 0;
	DimFP = 110;
	DistFP = 1;
	SeedFP = 1983;
	PotPar = 6;
	MEffect = 1;
	STEffect = 1;
	ZeroDefTyp = 3;
	ZeroDefPlane = 0.0;
	ApproxModel = 1;
	ThkTyp = 1;
	nThk = 0;
	Thk = 0;
	BWL = true;
	PBC_xy = true;
	Vrl = stVrl;
	E0 = 300;
	theta = 0;	
	phi = 0;
	lx = 4.078;
	ly = 4.078;	
	dz = 0.25;
	nx = 1024;
	ny = 1024;
}

cMT_MGP_CPU::~cMT_MGP_CPU(){
	freeMemory();
}

cMT_MGP_CPU& cMT_MGP_CPU::operator= (const cMT_MGP_CPU &MT_MGP_CPU){
	freeMemory();

	gpu = MT_MGP_CPU.gpu;
	SimType = MT_MGP_CPU.SimType;
	MulOrder = MT_MGP_CPU.MulOrder;
	nConfFP = MT_MGP_CPU.nConfFP;
	DimFP = MT_MGP_CPU.DimFP;
	DistFP = MT_MGP_CPU.DistFP;
	SeedFP = MT_MGP_CPU.SeedFP;
	PotPar = MT_MGP_CPU.PotPar;
	MEffect = MT_MGP_CPU.MEffect;
	STEffect = MT_MGP_CPU.STEffect;
	ZeroDefTyp = MT_MGP_CPU.ZeroDefTyp;
	ZeroDefPlane = MT_MGP_CPU.ZeroDefPlane;
	ApproxModel = MT_MGP_CPU.ApproxModel;
	ThkTyp = MT_MGP_CPU.ThkTyp;
	nThk = MT_MGP_CPU.nThk;
	Thk = new double[nThk];
	memcpy(Thk, MT_MGP_CPU.Thk, nThk*cSizeofRD);
	BWL = MT_MGP_CPU.BWL;
	PBC_xy = MT_MGP_CPU.PBC_xy;
	Vrl = MT_MGP_CPU.Vrl;
	E0 = MT_MGP_CPU.E0;
	theta = MT_MGP_CPU.theta;	
	phi = MT_MGP_CPU.phi;
	lx = MT_MGP_CPU.lx;
	ly = MT_MGP_CPU.ly;	
	dz = MT_MGP_CPU.dz;
	nx = MT_MGP_CPU.nx;
	ny = MT_MGP_CPU.ny;

	return *this;
}

void cMT_MGP_CPU::SetInputData(cMT_InMULTEM_CPU &MT_InMULTEM_CPU){
	freeMemory();

	gpu = MT_InMULTEM_CPU.gpu;
	SimType = MT_InMULTEM_CPU.SimType;	
	MulOrder = 2;
	nConfFP = MT_InMULTEM_CPU.nConfFP;		
	DimFP = MT_InMULTEM_CPU.DimFP;	
	DistFP = 1;
	SeedFP = MT_InMULTEM_CPU.SeedFP;
	PotPar = MT_InMULTEM_CPU.PotPar;
	MEffect = MT_InMULTEM_CPU.MEffect;
	STEffect = MT_InMULTEM_CPU.STEffect;
	ZeroDefTyp = MT_InMULTEM_CPU.ZeroDefTyp;
	ZeroDefPlane = MT_InMULTEM_CPU.ZeroDefPlane;
	ApproxModel = MT_InMULTEM_CPU.ApproxModel;
	ThkTyp = MT_InMULTEM_CPU.ThicknessTyp;
	nThk = (ThkTyp==1)?1:MT_InMULTEM_CPU.nThickness;
	Thk = new double[nThk];
	memcpy(Thk, MT_InMULTEM_CPU.Thickness, nThk*cSizeofRD); // change
	BWL = (MT_InMULTEM_CPU.BandwidthLimit==1)?true:false;
	Vrl = stVrl;
	E0 = MT_InMULTEM_CPU.E0;	
	theta = MT_InMULTEM_CPU.theta;	
	phi = MT_InMULTEM_CPU.phi;
	lx = MT_InMULTEM_CPU.lx;
	ly = MT_InMULTEM_CPU.ly;
	dz = MT_InMULTEM_CPU.dz;
	nx = MT_InMULTEM_CPU.nx;
	ny = MT_InMULTEM_CPU.ny;
	if(ApproxModel>1){
		MulOrder = 1;
		DimFP = DimFP - DimFP%10;
	}
}