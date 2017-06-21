function [Crys3D, lx, ly, lz, a, b, c, dz] = GaAs001Crystal(na, nb, nc, ncu, rms3d)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 5.6537; 
b = 5.6537; 
c = 5.6537;
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 4;
occ = 1;region = 0;charge = 0;
% Ga = 31, As = 33;
% Z charge x y z rms3d occupancy region charge
CrysPar.uLayer(1).atoms = [31, 0.00, 0.00, 0.00, rms3d, occ, region, charge; 31, 0.50, 0.50, 0.00, rms3d, occ, region, charge]; %Ga
CrysPar.uLayer(2).atoms = [33, 0.25, 0.25, 0.25, rms3d, occ, region, charge; 33, 0.75, 0.75, 0.25, rms3d, occ, region, charge]; %As
CrysPar.uLayer(3).atoms = [31, 0.00, 0.50, 0.50, rms3d, occ, region, charge; 31, 0.50, 0.00, 0.50, rms3d, occ, region, charge]; %Ga
CrysPar.uLayer(4).atoms = [33, 0.75, 0.25, 0.75, rms3d, occ, region, charge; 33, 0.25, 0.75, 0.75, rms3d, occ, region, charge]; %As

Crys3D = il_crystal_by_lays(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;