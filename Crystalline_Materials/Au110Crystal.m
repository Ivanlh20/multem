function [Crys3D, lx, ly, lz, a, b, c, dz] = Au110Crystal(na, nb, nc, ncu, sigma)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 4.0780/sqrt(2); 
b = 4.0780; 
c = 4.0780/sqrt(2);
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 2;
% x y z Z sigma occupancy
CrysPar.uLayer(1).Atoms = [0.00, 0.00, 0.00, 79, sigma, 1.0];
CrysPar.uLayer(2).Atoms = [0.50, 0.50, 0.50, 79, sigma, 1.0];
Crys3D = get_CrystalbyLayers_CPU(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;