% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>

clear; clc;
addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])

xtl_parm.na = 10;
xtl_parm.nb = 10;
xtl_parm.nc = 10;
xtl_parm.a = 4.0780;
xtl_parm.b = 4.0780;
xtl_parm.c = 4.0780;
xtl_parm.nuLayer = 2;
occ = 1;
region = 0;
charge = 0;
% Au = 79
%Z x y z sigma occupancy
rmsd_3d = 0.085;
xtl_parm.uLayer(1).atoms = [79, 0.0, 0.0, 0.0, rmsd_3d, occ, region, charge; 79, 0.5, 0.5, 0.0, rmsd_3d, occ, region, charge]; 
xtl_parm.uLayer(2).atoms = [79, 0.0, 0.5, 0.5, rmsd_3d, occ, region, charge; 79, 0.5, 0.0, 0.5, rmsd_3d, occ, region, charge];

tic;
Crys3D = ilc_crystal_by_lays(xtl_parm);
toc;

na = 1;
nb = 1;
nc = 1;
[Crys3D, lx, ly, lz, a, b, c, dz] = SrTiO3001_xtl(na, nb, nc, 2, 0.085);
[lx, ly]
% show crystal
clf;
ilm_show_crystal(1, Crys3D);