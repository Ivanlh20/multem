function [Crys3D, lx, ly, lz, a, b, c, dz] = GaAs001Crystal(na, nb, nc, ncu, sigma)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 5.6537; 
b = 5.6537; 
c = 5.6537;
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 4;
% Ga = 31, As = 33;
% x y z Z sigma occupancy
CrysPar.uLayer(1).Atoms = [0.00, 0.00, 0.00, 31, sigma, 1; 0.50, 0.50, 0.00, 31, sigma, 1]; %Ga
CrysPar.uLayer(2).Atoms = [0.25, 0.25, 0.25, 33, sigma, 1; 0.75, 0.75, 0.25, 33, sigma, 1]; %As
CrysPar.uLayer(3).Atoms = [0.00, 0.50, 0.50, 31, sigma, 1; 0.50, 0.00, 0.50, 31, sigma, 1]; %Ga
CrysPar.uLayer(4).Atoms = [0.75, 0.25, 0.75, 33, sigma, 1; 0.25, 0.75, 0.75, 33, sigma, 1]; %As

Crys3D = get_CrystalbyLayers_CPU(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;