% output_multislice = input_multem.ilc_multem perform TEM simulation
% 
% Exit wave real space (EWRS) simulation
% 
% All parameters of the input_multem structure are explained in ilm_dflt_input_multem()
% 
% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>

clear; clc;
addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])

% create specimen
lx = 100;
ly = 100;
lz = 100;

na = 8; nb = 8; nc = 8; ncu = 2; rmsd_3d = 0.085;

[atoms, ~] = Au001_xtl(na, nb, nc, ncu, rmsd_3d);
atoms = ilm_center_spec(atoms, lx, ly, lz);

theta = 45;                                 % angle (�)
u0 = [1 1 0];                               % unitary vector			
rot_point_type = 1;                         % 1: geometric center, 2: User define		
p0 = [0 0 0];                               % rotation point

% rotate specimen
atoms_r = ilc_spec_rot(atoms, theta, u0, rot_point_type, p0);
figure(1); clf;

subplot(1, 2, 1);
plot3(atoms(:, 2), atoms(:, 3), atoms(:, 4), 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'auto');
axis equal;
axis([0 lx 0 ly 0 lz]);
view([1 0 1]);

subplot(1, 2, 2);
plot3(atoms_r(:, 2), atoms_r(:, 3), atoms_r(:, 4), 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'auto');
axis equal;
axis([0 lx 0 ly 0 lz]);
view([0 0 1]);