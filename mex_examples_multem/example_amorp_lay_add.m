clear; clc;
addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])

na = 10; nb = 10; nc = 20; ncu = 2;rmsd_3d = 0.085;

[atoms, lx, ly, lz, a, b, c, input_multem.spec_dz] = Au001_xtl(na, nb, nc, ncu, rmsd_3d);
% ilm_show_xtl(1, atoms)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bs = [lx, ly, 30];

Z = 6;
rms_3d = 0.09;
occ = 1.0;
region = 0;
d_min = 1.4;
rho = 2.2;
lay_pos = 21; %1: top, 2: bottom
seed = 1983;

tic;
atoms = ilc_amorp_lay_add(atoms(:, 1:5), Z, rms_3d, occ, region, bs, d_min, rho, lay_pos, seed);
toc;

ilm_show_xtl(1, atoms)
view([1 0 0])
zlim([min(atoms(:, 4)), max(atoms(:, 4))])

disp([lx, ly])
disp([min(atoms(:, 2)), max(atoms(:, 2))])
disp([min(atoms(:, 3)), max(atoms(:, 3))])
disp([min(atoms(:, 4)), max(atoms(:, 4))])

nbins = round((max(atoms(:, 4))-min(atoms(:, 4)))/0.1);
tic;
[x, y] = ilc_hist(atoms(:, 4), nbins);
toc;

atoms(:, 4) = atoms(:, 4)-min(atoms(:, 4));
% ilm_write_ap_pdb('amorphous.pdb', atoms, a, b, c, 90, 90, 90);

figure(2); clf;
subplot(1, 2, 1);
plot(x, y, '-r');
tic;
[r, rdf] = ilc_rdf(atoms(:, 2:4), 8, 200);
toc;
subplot(1, 2, 2);
plot(r, rdf, '-+r');