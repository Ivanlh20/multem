function [Crys3D, lx, ly, lz, a, b, c, dz] = Diamond110Crystal(na, nb, nc, ncu, rms3d)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 3.56/sqrt(2); 
b = 3.56; 
c = sqrt(2)*3.56/2;
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 2;
occ = 1;region = 0;charge = 0;
% C = 6
% Z x y z rms3d occupancy charge 
CrysPar.uLayer(1).atoms = [6, 0.00, 0.00, 0.00, rms3d, occ, region, charge; 6, 0.50, 0.75, 0.00, rms3d, occ, region, charge];
CrysPar.uLayer(2).atoms = [6, 0.00, 0.25, 0.50, rms3d, occ, region, charge; 6, 0.50, 0.50, 0.50, rms3d, occ, region, charge];
Crys3D = il_crystal_by_lays(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;