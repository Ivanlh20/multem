% output_multislice = il_MULTEM(system_conf, input_multislice) perform TEM simulation
% 
% Exit wave real space (EWRS) simulation
% 
% All parameters of the input_multislice structure are explained in multem_default_values()
% 
% Copyright 2017 Ivan Lobato <Ivanlh20@gmail.com>

clear all; clc;

% create specimen
lx = 100;
ly = 100;
lz = 100;

na = 8; nb = 8; nc = 8; ncu = 2; rms3d = 0.085;

[atoms, ~] = Au001Crystal(na, nb, nc, ncu, rms3d);
atoms = center_spec(atoms, lx, ly, lz);

theta = 45;                                 % angle (º)
u0 = [1 0 0];                               % unitary vector			
rot_point_type = 1;                         % 1: geometric center, 2: User define		
p0 = [0 0 0];                               % rotation point

% rotate specimen
atoms_r = il_spec_rot(atoms, theta, u0, rot_point_type, p0);
figure(1); clf;

subplot(1, 2, 1);
plot3(atoms(:, 2), atoms(:, 3), atoms(:, 4), 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'auto');
axis equal;
axis([0 lx 0 ly 0 lz]);
view([1 0 0]);

subplot(1, 2, 2);
plot3(atoms_r(:, 2), atoms_r(:, 3), atoms_r(:, 4), 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'auto');
axis equal;
axis([0 lx 0 ly 0 lz]);
view([1 0 0]);