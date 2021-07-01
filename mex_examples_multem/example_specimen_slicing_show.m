clear;clc;

input_multem = ilm_dflt_input_multem(); % Load default values;

input_multem.atomic_vib_mod = 1; % eavm_still_atom = 1, eavm_absorptive_pot = 2, eavm_frozen_phonon = 3, eavm_user_def = 4
input_multem.elec_spec_interac_mod = 1; % eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
input_multem.spec_slic(1).typ = 1; % esst_plns_proj = 1, esst_dz_proj = 2, esst_plns_sub = 3, esst_dz_sub = 4, esst_user_def = 5, esst_auto = 6
input_multem.atomic_vib_dim = 111;
input_multem.atomic_vib_seed = 300183;
input_multem.atomic_vib_nconf = 1;

input_multem.spec_rot_theta = 0; % final angle
input_multem.spec_rot_u_0 = [1 0 0]; % unitary vector			
input_multem.spec_rot_ctr_typ = 1; % 1: geometric center, 2: User define		
input_multem.spec_rot_ctr_p = [0 0 0]; % rotation point

na = 4;nb = 4;nc = 10;ncu = 4;rmsd_3d = 0.08;

[input_multem.spec_atoms, input_multem.spec_bs_x...
, input_multem.spec_bs_y, input_multem.spec_bs_z...
, a, b, c, input_multem.spec_slic(1).sli_thick] = GaAs001_xtl(na, nb, nc, ncu, rmsd_3d);

input_multem.spec_slic(1).sli_thick = 5;
% get spec slicing
tic;
input_multem.atomic_vib_mod = 1;
[atoms0, Slice0] = ilc_spec_slicing(input_multem);
toc;

[nslice0, ~] = size(Slice0);

tic;
input_multem.atomic_vib_mod = 3;
[atoms, Slice] = ilc_spec_slicing(input_multem);
toc;

[nslice, ~] = size(Slice);
x_min = min(atoms(:, 2))-2;
x_max = max(atoms(:, 2))+2;
z_min = min(atoms(:, 4))-5;
z_max = max(atoms(:, 4))+5;

figure(1); clf;
plot(atoms(:, 3), atoms(:, 4), 'ok');
set(gca, 'FontSize', 12, 'LineWidth', 1, 'PlotBoxAspectRatio', [1.25 1 1]);
title('Atomic positions');
ylabel('y', 'FontSize', 14);
xlabel('x', 'FontSize', 12);
axis equal;
axis([x_min x_max z_min z_max]);

for i = 1:nslice
    hold on;
    plot([x_min x_max], [Slice(i, 1) Slice(i, 1)], '-b', [x_min x_max], [Slice(i, 2) Slice(i, 2)], '-b');
    axis equal;
    axis([x_min x_max z_min z_max]);
end

for i = 1:nslice0
    hold on;
    plot([x_min x_max], [Slice0(i, 1) Slice0(i, 1)], '-r', [x_min x_max], [Slice0(i, 2) Slice0(i, 2)], '-r');
    axis equal;
    axis([x_min x_max z_min z_max]);
end