% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>

clear; clc;
addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])

na = 20; nb = 20; nc = 30; ncu = 2; rmsd_3d = 0.085;

[atoms, lx, ly, ~, a, b, c, input_multem.spec_dz] = Au001_xtl(na, nb, nc, ncu, rmsd_3d);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lz = 20;
Z = 6;
rms_3d = 0.09;
d_min = 1.4;
seed = 1983;
rho = 2.2;
lay_pos = 1; %1: top, 2: bottom

tic;
atoms = ilc_add_amorp_lay(atoms, lx, ly, lz, d_min, Z, rms_3d, rho, lay_pos, seed);
toc;

ilm_show_crystal(1, atoms)
view([1 0 0])
zlim([min(atoms(:, 4)), max(atoms(:, 4))])

disp([lx, ly])
disp([min(atoms(:, 2)), max(atoms(:, 2))])
disp([min(atoms(:, 3)), max(atoms(:, 3))])
disp([min(atoms(:, 4)), max(atoms(:,4))])

% atoms(:, 4) = atoms(:, 4)-min(atoms(:, 4));
% save_atomic_position_pdb('amorphous.pdb', atoms, a, b, c, 90, 90, 90);

figure(2); clf;
tic;
[r, rdf] = ilc_rdf_3d(atoms, 8, 200);
toc;
plot(r, rdf,'-+r');