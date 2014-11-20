function [Crys3D, lx, ly, lz, a, b, c, dz] = SrTiO3001Crystal(na, nb, nc, ncu, sigma)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 3.9050; 
b = 3.9050; 
c = 3.9050;
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 2;
% Sr = 38, Ti = 22; O = 8
% x y z Z sigma occupancy
CrysPar.uLayer(1).Atoms = [0.0, 0.0, 0.0, 38, sigma, 1; 0.5, 0.5, 0.0, 8, sigma, 1];
CrysPar.uLayer(2).Atoms = [0.0, 0.5, 0.5, 8, sigma, 1; 0.5, 0.0, 0.5, 8, sigma, 1; 0.5, 0.5, 0.5, 22, sigma, 1];
Crys3D = getCrystalbyLayers(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;