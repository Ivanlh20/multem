function [Crys3D, lx, ly, lz, a, b, c, dz] = Au001Crystal(na, nb, nc, ncu, rms3d)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 4.0780; 
b = 4.0780; 
c = 4.0780;
% a = 4.07003925431053; 
% b = a; 
% c = a;
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 2;
charge = 0;
% Au = 79
% Z charge x y z rms3d occupancy
CrysPar.uLayer(1).atoms = [79, charge, 0.0, 0.0, 0.0, rms3d, 1; 79, charge, 0.5, 0.5, 0.0, rms3d, 1];
CrysPar.uLayer(2).atoms = [79, charge, 0.0, 0.5, 0.5, rms3d, 1; 79, charge, 0.5, 0.0, 0.5, rms3d, 1];
Crys3D = il_crystal_by_layers(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;