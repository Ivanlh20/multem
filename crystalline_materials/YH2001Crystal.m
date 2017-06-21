function [Crys3D, lx, ly, lz, a, b, c, dz] = YH2001Crystal(na, nb, nc, ncu, rms3d_Y, rms3d_H)
CrysPar.na = na;
CrysPar.nb = nb;
CrysPar.nc = nc;
a = 5.2000; 
b = 5.2000; 
c = 5.2000; 
CrysPar.a = a;
CrysPar.b = b;
CrysPar.c = c;
CrysPar.nuLayer = 4;
occ = 1;
region = 0;
charge = 0;
% Au = 79
% Z x y z rms3d occupancy charge
CrysPar.uLayer(1).atoms = [39, 0.0, 0.0, 0.0, rms3d_Y, occ, region, charge; 39, 0.5, 0.5, 0.0, rms3d_Y, occ, region, charge];
CrysPar.uLayer(2).atoms = [1, 0.25, 0.25, 0.25, rms3d_H, occ, region, charge; 1, 0.75, 0.25, 0.25, rms3d_H, occ, region, charge; 1, 0.25, 0.75, 0.25, rms3d_H, occ, region, charge; 1, 0.75, 0.75, 0.25, rms3d_H, occ, region, charge];
CrysPar.uLayer(3).atoms = [39, 0.0, 0.5, 0.5, rms3d_Y, occ, region, charge; 39, 0.5, 0.0, 0.5, rms3d_Y, occ, region, charge];
CrysPar.uLayer(4).atoms = [1, 0.25, 0.25, 0.75, rms3d_H, occ, region, charge; 1, 0.75, 0.25, 0.75, rms3d_H, occ, region, charge; 1, 0.25, 0.75, 0.75, rms3d_H, occ, region, charge; 1, 0.75, 0.75, 0.75, rms3d_H, occ, region, charge];

Crys3D = il_crystal_by_lays(CrysPar);

dz = CrysPar.c/ncu;
lx = na*CrysPar.a; ly = nb*CrysPar.b; lz = nc*CrysPar.c;