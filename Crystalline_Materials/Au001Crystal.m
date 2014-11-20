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
CrysPar.uLayer(1).Atoms = [0.0, 0.0, 0.0, 79, sigma, 1; 0.5, 0.5, 0.0, 79, sigma, 1];
CrysPar.uLayer(2).Atoms = [0.0, 0.5, 0.5, 79, sigma, 1; 0.5, 0.0, 0.5, 79, sigma, 1];
Crys3D = getCrystalbyLayers(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;