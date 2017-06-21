function [Crys3D, lx, ly, lz, a, b, c, dz] = PdH001Crystal(na, nb, nc, ncu, rms3d)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 4.025; 
b = 4.025; 
c = 4.025;
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 2;
occ = 1;
region = 0;
charge = 0;
% Pd = 46, H = 1
% Z charge x y z rms3d occupancy region charge
CrysPar.uLayer(1).atoms = [46, 0.0, 0.0, 0.0, rms3d, occ, region, charge; 46, 0.5, 0.5, 0.0, rms3d, occ, region, charge;...
                           1, 0.5, 0.0, 0.0, rms3d, occ, region, charge; 1, 0.0, 0.5, 0.0, rms3d, occ, region, charge];
CrysPar.uLayer(2).atoms = [46, 0.0, 0.5, 0.5, rms3d, occ, region, charge; 46, 0.5, 0.0, 0.5, rms3d, occ, region, charge;...
                           1, 0.0, 0.0, 0.5, rms3d, occ, region, charge; 1, 0.5, 0.5, 0.5, rms3d, occ, region, charge];
Crys3D = il_crystal_by_lays(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;