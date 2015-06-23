function [Crys3D, lx, ly, lz, a, b, c, dz] = Au001Crystal(na, nb, nc, ncu, sigma)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 4.0780; 
b = 4.0780; 
c = 4.0780;
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 2;
% x y z Z sigma occupancy
CrysPar.uLayer(1).atoms = [79, 0.0, 0.0, 0.0, sigma, 1; 79, 0.5, 0.5, 0.0, sigma, 1];
CrysPar.uLayer(2).atoms = [79, 0.0, 0.5, 0.5, sigma, 1; 79, 0.5, 0.0, 0.5, sigma, 1];
Crys3D = get_crystal_by_layers(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;