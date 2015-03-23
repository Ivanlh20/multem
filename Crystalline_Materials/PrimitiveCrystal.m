function [Crys3D, lx, ly, lz, a, b, c, dz] = PrimitiveCrystal(na, nb, nc, ncu, sigma)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 4.0780; 
b = 4.0780; 
c = 4.0780;
a = 6.50; 
b = 6.50;
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 1;
% x y z Z sigma occupancy
CrysPar.uLayer(1).Atoms = [0.5, 0.5, 0.0, 79, sigma, 1];
Crys3D = get_CrystalbyLayers_CPU(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;