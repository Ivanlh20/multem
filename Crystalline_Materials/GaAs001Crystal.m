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
CrysPar.uLayer(1).atoms = [31, 0.00, 0.00, 0.00, sigma, 1; 31, 0.50, 0.50, 0.00, sigma, 1]; %Ga
CrysPar.uLayer(2).atoms = [33, 0.25, 0.25, 0.25, sigma, 1; 33, 0.75, 0.75, 0.25, sigma, 1]; %As
CrysPar.uLayer(3).atoms = [31, 0.00, 0.50, 0.50, sigma, 1; 31, 0.50, 0.00, 0.50, sigma, 1]; %Ga
CrysPar.uLayer(4).atoms = [33, 0.75, 0.25, 0.75, sigma, 1; 33, 0.25, 0.75, 0.75, sigma, 1]; %As

Crys3D = get_crystal_by_layers(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;