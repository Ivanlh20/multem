function [Crys3D, lx, ly, lz, a, b, c, dz] = SrTiO3110Crystal(na, nb, nc, ncu, sigma)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 3.9050; 
b = 3.9050*sqrt(2); 
c = 3.9050*sqrt(2);
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 4;
% x y z Z sigma occupancy
CrysPar.uLayer(1).Atoms = [0.50, 0.00, 0.00, 38, sigma, 1; 0.50, 0.50, 0.00, 8, sigma, 1; 0.0, 0.5, 0.0, 22, sigma, 1];
CrysPar.uLayer(2).Atoms = [0.00, 0.25, 0.25, 8, sigma, 1; 0.00, 0.75, 0.25, 8, sigma, 1];
CrysPar.uLayer(3).Atoms = [0.50, 0.00, 0.50, 8, sigma, 1; 0.50, 0.50, 0.50, 38, sigma, 1; 0.00, 0.00, 0.50, 22, sigma, 1];
CrysPar.uLayer(4).Atoms = [0.00, 0.25, 0.75, 8, sigma, 1; 0.00, 0.75, 0.75, 8, sigma, 1];

Crys3D = get_CrystalbyLayers_CPU(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;