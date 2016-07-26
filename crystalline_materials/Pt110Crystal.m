function [Crys3D, lx, ly, lz, a, b, c, dz] = Pt110Crystal(na, nb, nc, ncu, rms3d)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 3.9242/sqrt(2); 
b = 3.9242; 
c = 3.9242/sqrt(2);
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 2;
charge = 0;
% Pt = 78
% Z charge x y z rms3d occupancy
CrysPar.uLayer(1).atoms = [78, 0.00, 0.00, 0.00, rms3d, 1.0, charge];
CrysPar.uLayer(2).atoms = [78, 0.50, 0.50, 0.50, rms3d, 1.0, charge];
Crys3D = il_crystal_by_layers(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;