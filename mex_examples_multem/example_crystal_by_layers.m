clc;
clear all;
CrysPar.na = 10;
CrysPar.nb = 10;
CrysPar.nc = 10;
CrysPar.a = 4.0780;
CrysPar.b = 4.0780;
CrysPar.c = 4.0780;
CrysPar.nuLayer = 2;
occ = 1;
region = 0;
charge = 0;
% Au = 79
%Z x y z sigma occupancy
rmsd_3d = 0.085;
CrysPar.uLayer(1).atoms = [79, 0.0, 0.0, 0.0, rmsd_3d, occ, region, charge; 79, 0.5, 0.5, 0.0, rmsd_3d, occ, region, charge]; 
CrysPar.uLayer(2).atoms = [79, 0.0, 0.5, 0.5, rmsd_3d, occ, region, charge; 79, 0.5, 0.0, 0.5, rmsd_3d, occ, region, charge];

tic;
Crys3D = il_crystal_by_lays(CrysPar);
toc;

na = 1;
nb = 1;
nc = 1;
[Crys3D, lx, ly, lz, a, b, c, dz] = SrTiO3001Crystal(na, nb, nc, 2, 0.085);
[lx, ly]
% show crystal
clf;
show_crystal(1, Crys3D);