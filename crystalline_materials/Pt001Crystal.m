function [Crys3D, lx, ly, lz, a, b, c, dz] = Pt001Crystal(na, nb, nc, ncu, rms3d)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 3.9242; 
b = 3.9242; 
c = 3.9242;
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 2;
occ = 1;
region = 0;
charge = 0;
% Au = 79
% Z x y z rms3d occupancy charge
CrysPar.uLayer(1).atoms = [79, 0.0, 0.0, 0.0, rms3d, occ, region, charge; 79, 0.5, 0.5, 0.0, rms3d, occ, region, charge];
CrysPar.uLayer(2).atoms = [79, 0.0, 0.5, 0.5, rms3d, occ, region, charge; 79, 0.5, 0.0, 0.5, rms3d, occ, region, charge];
Crys3D = il_crystal_by_lays(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;