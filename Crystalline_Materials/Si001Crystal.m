function [Crys3D, lx, ly, lz, a, b, c, dz] = Si001Crystal(na, nb, nc, ncu, rms3d)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 5.4307; 
b = 5.4307; 
c = 5.4307;
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 4;
charge = 0;
% Si = 14
% Z charge x y z rms3d occupancy
CrysPar.uLayer(1).atoms = [14, 0.00, 0.00, 0.00, rms3d, 1, charge; 14, 0.50, 0.50, 0.00, rms3d, 1, charge];
CrysPar.uLayer(2).atoms = [14, 0.25, 0.25, 0.25, rms3d, 1, charge; 14, 0.75, 0.75, 0.25, rms3d, 1, charge];
CrysPar.uLayer(3).atoms = [14, 0.00, 0.50, 0.50, rms3d, 1, charge; 14, 0.50, 0.00, 0.50, rms3d, 1, charge];
CrysPar.uLayer(4).atoms = [14, 0.25, 0.75, 0.75, rms3d, 1, charge; 14, 0.75, 0.25, 0.75, rms3d, 1, charge];
Crys3D = il_crystal_by_layers(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;