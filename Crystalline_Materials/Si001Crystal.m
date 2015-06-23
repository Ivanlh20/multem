function [Crys3D, lx, ly, lz, a, b, c, dz] = Si001Crystal(na, nb, nc, ncu, sigma)
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
% x y z Z sigma occupancy
CrysPar.uLayer(1).atoms = [14, 0.00, 0.00, 0.00, sigma, 1; 14, 0.50, 0.50, 0.00, sigma, 1];
CrysPar.uLayer(2).atoms = [14, 0.25, 0.25, 0.25, sigma, 1; 14, 0.75, 0.75, 0.25, sigma, 1];
CrysPar.uLayer(3).atoms = [14, 0.00, 0.50, 0.50, sigma, 1; 14, 0.50, 0.00, 0.50, sigma, 1];
CrysPar.uLayer(4).atoms = [14, 0.25, 0.75, 0.75, sigma, 1; 14, 0.75, 0.25, 0.75, sigma, 1];
Crys3D = get_crystal_by_layers(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;