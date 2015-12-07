clc;
clear all;
CrysPar.na = 8;
CrysPar.nb = 8;
CrysPar.nc = 20;
CrysPar.a = 4.0780;
CrysPar.b = 4.0780;
CrysPar.c = 4.0780;
CrysPar.nuLayer = 2;
charge = 0;
% Au = 79
% x y z Z sigma occupancy
CrysPar.uLayer(1).atoms = [79, charge, 0.0, 0.0, 0.0, 0.084, 1; 79, charge, 0.5, 0.5, 0.0, 0.084, 1]; 
CrysPar.uLayer(2).atoms = [79, charge, 0.0, 0.5, 0.5, 0.084, 1; 79, charge, 0.5, 0.0, 0.5, 0.084, 1];

tic;
Crys3D = il_crystal_by_layers(CrysPar);
toc;

figure(1);
plot3(Crys3D(:, 3), Crys3D(:, 4), Crys3D(:, 5), '*r');
axis equal;
length(Crys3D)
