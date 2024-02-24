% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>

clear; clc;
root = '/media/ivan/Data_4TB/multem_gh';
addpath([root filesep 'mex_bin'])
addpath([root filesep 'crystalline_materials'])
addpath([root filesep 'matlab_functions'])

lx = 50;
ly = 50;
lz = 20;
Z = 6;
rms_3d = 0.09;
d_min = 1.4;
seed = 1983;
rho = 2.2;

tic;
atoms = ilc_amorp_spec(lx, ly, lz, d_min, Z, rms_3d, rho, seed);
toc;

% path = strcat('input_spec\Si_',num2str(lx), 'x', num2str(ly), 'x', num2str(lz), '_', num2str(seed), '.mat');
% save(path, 'atoms');
% disp([iseed, lz]);
figure(1); clf;
plot3(atoms(:, 2), atoms(:, 3), atoms(:, 4),'.r');
axis equal;

NA = 6.022140857e+23;
d = 0;
lx0 = lx-d;
ly0 = ly-d;
lz0 = lz;
ii = find((0.5*d<=atoms(:,2))&(atoms(:,2)<=lx-0.5*d)&(0.5*d<=atoms(:,3))&(atoms(:,3)<=ly-0.5*d));
density = length(atoms)*12.011/(lx0*ly0*lz0*NA*(1e-8)^3);
disp(['density = ', num2str(density)])

tic;
[r, rdf] = ilc_rdf_3d(atoms(1:5e+4, :), 8, 200);
toc;

figure(2); clf;
plot(r, rdf,'-+r');