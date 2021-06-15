clear; clc;
addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bs = [50, 50, 50];

Z = 6;
rms_3d = 0.09;
occ = 1.0;
region = 0;
d_min = 1.4;
rho = 2.2;
seed = 1983;

tic;
atoms = ilc_amorp_build(Z, rms_3d, occ, region, bs, d_min, rho, seed);
toc;

disp(size(atoms))
figure(1); clf;
plot3(atoms(:, 2), atoms(:, 3), atoms(:, 4), '.r');
axis equal;

NA = 6.022140857e+23;

d = 0;
lx0 = bs(1)-d;
ly0 = bs(2)-d;
lz0 = bs(3);
ii = find((0.5*d<=atoms(:, 2)) & (atoms(:, 2)<=bs(1)-0.5*d) & (0.5*d<=atoms(:, 3)) & (atoms(:, 3)<=bs(2)-0.5*d));
% atoms = atoms(ii, :);
density = length(atoms)*12.011/(lx0*ly0*lz0*NA*(1e-8)^3);

disp([rho, density])
tic;
[r, rdf] = ilc_rdf(atoms(:, 2:4), 8, 200);
toc;


figure(2);clf;
plot(r, rdf, '-+r');