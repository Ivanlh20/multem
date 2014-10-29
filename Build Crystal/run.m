clc;
clear all;
CrysPar.na = 1;
CrysPar.nb = 1;
CrysPar.nc = 1;
CrysPar.a = 4.0780;
CrysPar.b = 4.0780;
CrysPar.c = 4.0780;
CrysPar.nuLayer = 2;
% x y z Z sigma occupancy
CrysPar.uLayer(1).Atoms = [0.0, 0.0, 0.0, 79, 0.084, 1; 0.5, 0.5, 0.0, 79, 0.084, 1]; 
CrysPar.uLayer(2).Atoms = [0.0, 0.5, 0.5, 79, 0.084, 1; 0.5, 0.0, 0.5, 79, 0.084, 1];

tic;
Crys3D = CreateCrystalbyLayers(CrysPar);
toc;

figure(1);
plot3(Crys3D(:, 1), Crys3D(:, 2), Crys3D(:, 3), '*r');
axis equal;
length(Crys3D)
