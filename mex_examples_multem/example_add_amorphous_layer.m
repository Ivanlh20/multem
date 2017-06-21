clear all; clc;

na = 20; nb = 20; nc = 30; ncu = 2; rms3d = 0.085;

[atoms, lx, ly, ~, a, b, c, input_multislice.spec_dz] = Au001Crystal(na, nb, nc, ncu, rms3d);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lz = 20;
Z = 6;
rms_3d = 0.09;
d_min = 1.4;
seed = 1983;
rho = 2.2;
lay_pos = 1; %1: top, 2: bottom

tic;
atoms = il_add_amorp_lay(atoms, lx, ly, lz, d_min, Z, rms_3d, rho, lay_pos, seed);
toc;

show_crystal(1, atoms)
view([1 0 0])
zlim([min(atoms(:, 4)), max(atoms(:, 4))])

disp([lx, ly])
disp([min(atoms(:, 2)), max(atoms(:, 2))])
disp([min(atoms(:, 3)), max(atoms(:, 3))])
disp([min(atoms(:, 4)), max(atoms(:,4))])

nbins = round((max(atoms(:, 4))-min(atoms(:, 4)))/0.1);
tic;
[x, y] = il_hist(atoms(:, 4), nbins);
toc;

atoms(:, 4) = atoms(:, 4)-min(atoms(:, 4));
save_atomic_position_pdb('amorphous.pdb', atoms, a, b, c, 90, 90, 90);

figure(2); clf;
subplot(1, 2, 1);
plot(x, y, '-r');
tic;
[r, rdf] = il_rdf_3d(atoms, 8, 200);
toc;
subplot(1, 2, 2);
plot(r, rdf,'-+r');