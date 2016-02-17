clc;
clear all;
CrysPar.na = 8;
CrysPar.nb = 8;
CrysPar.nc = 20;
CrysPar.a = 4.0780;
CrysPar.b = 4.0780;
CrysPar.c = 4.0780;
CrysPar.nuLayer = 2;
charge = 0;
% Au = 79
% x y z Z sigma occupancy
rmsd_3d = 0.085;
CrysPar.uLayer(1).atoms = [79, 0.0, 0.0, 0.0, rmsd_3d, 1, charge; 79, 0.5, 0.5, 0.0, rmsd_3d, 1, charge]; 
CrysPar.uLayer(2).atoms = [79, 0.0, 0.5, 0.5, rmsd_3d, 1, charge; 79, 0.5, 0.0, 0.5, rmsd_3d, 1, charge];

tic;
Crys3D = il_crystal_by_layers(CrysPar);
toc;

% show crystal
show_crystal(1, Crys3D);