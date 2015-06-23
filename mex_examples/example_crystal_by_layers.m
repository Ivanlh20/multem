clc;
clear all;
CrysPar.na = 1;
CrysPar.nb = 2;
CrysPar.nc = 1;
CrysPar.a = 4.0780;
CrysPar.b = 4.0780;
CrysPar.c = 4.0780;
CrysPar.nuLayer = 2;
% x y z Z sigma occupancy
CrysPar.uLayer(1).atoms = [79, 0.0, 0.0, 0.0, 0.084, 1; 79, 0.5, 0.5, 0.0, 0.084, 1]; 
CrysPar.uLayer(2).atoms = [79, 0.0, 0.5, 0.5, 0.084, 1; 79, 0.5, 0.0, 0.5, 0.084, 1];

tic;
Crys3D = get_crystal_by_layers(CrysPar);
toc;

figure(1);
plot3(Crys3D(:, 2), Crys3D(:, 3), Crys3D(:, 4), '*r');
axis equal;
length(Crys3D)
