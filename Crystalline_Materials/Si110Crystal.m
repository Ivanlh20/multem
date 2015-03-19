function [Crys3D, lx, ly, lz, a, b, c, dz] = Si110Crystal(na, nb, nc, ncu, sigma)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 5.4307/sqrt(2); 
b = 5.4307; 
c = sqrt(2)*5.4307/4;
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 2;
% x y z Z sigma occupancy
CrysPar.uLayer(1).Atoms = [0.00, 0.00, 0.00, 14, sigma, 1; 0.50, 0.75, 0.00, 14, sigma, 1];
CrysPar.uLayer(2).Atoms = [0.00, 0.25, 0.50, 14, sigma, 1; 0.50, 0.50, 0.50, 14, sigma, 1];
Crys3D = getCrystalbyLayers(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;