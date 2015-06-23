function [Crys3D, lx, ly, lz, a, b, c, dz] = Cu001Crystal(na, nb, nc, ncu, sigma)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 3.6150; 
b = 3.6150; 
c = 3.6150;
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 2;
% x y z Z sigma occupancy
CrysPar.uLayer(1).atoms = [29, 0.0, 0.0, 0.0, sigma, 1; 29, 0.5, 0.5, 0.0, sigma, 1];
CrysPar.uLayer(2).atoms = [29, 0.0, 0.5, 0.5, sigma, 1; 29, 0.5, 0.0, 0.5, sigma, 1];
Crys3D = get_crystal_by_layers(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;