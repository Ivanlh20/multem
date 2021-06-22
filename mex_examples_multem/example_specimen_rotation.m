% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
clear;clc;
addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])

% create specimen
lx = 100;
ly = 100;
lz = 100;

na = 8; nb = 8; nc = 8; ncu = 2; rmsd_3d = 0.085;
[atoms, ~] = SrTiO3001_xtl(na, nb, nc, ncu, rmsd_3d);
atoms = ilm_spec_recenter(atoms, lx, ly, lz);
atoms = atoms(:, 1:5);

ilm_show_xtl(1, atoms);

theta = 45; % angle (º)
u_0 = [1 0 0]; % unitary vector			
rot_point_type = 1; % 1: geometric center, 2: User define		
p_0 = [0 0 0]; % rotation point

% rotate specimen
tic;
atoms_r = ilc_spec_rot(atoms, theta, u_0, rot_point_type, p_0);
toc;

figure(1); clf;
subplot(1, 2, 1);
ilm_show_xtl(1, atoms, false);
view([1 0 1]);
title('Raw')

subplot(1, 2, 2);
ilm_show_xtl(1, atoms_r, false);
view([0 0 1]);
title('Rotated')