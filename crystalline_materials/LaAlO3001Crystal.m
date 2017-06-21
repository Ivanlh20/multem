function [Crys3D, lx, ly, lz, a, b, c, dz] = LaAlO3001Crystal(na, nb, nc, ncu, rms3d)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 3.7900; 
b = 3.79; 
c = 3.79;
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 2;
occ = 1;region = 0;charge = 0;
% Sr = 38, Ti = 22; O = 8
% La = 57, Al = 13; O = 8
% Z charge x y z rms3d occupancy region charge
CrysPar.uLayer(1).atoms = [57, 0.0, 0.0, 0.0, rms3d, occ, region, charge; 8, 0.5, 0.5, 0.0, rms3d, occ, region, charge];
CrysPar.uLayer(2).atoms = [8, 0.0, 0.5, 0.5, rms3d, occ, region, charge; 8, 0.5, 0.0, 0.5, rms3d, occ, region, charge; 13, 0.5, 0.5, 0.5, rms3d, occ, region, charge];
Crys3D = il_crystal_by_lays(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;