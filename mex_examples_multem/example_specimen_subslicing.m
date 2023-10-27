clear;clc;

input_multem = ilm_dflt_input_multem(); % Load default values;

input_multem.atomic_vib_mod = 1; % eavm_still_atom = 1, eavm_absorptive_pot = 2, eavm_frozen_phonon = 3, eavm_user_def = 4
input_multem.elec_spec_interact_mod = 1; % eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
input_multem.spec_slic(1).typ = 3; % esst_plns_proj = 1, esst_dz_proj = 2, esst_plns_sub = 3, esst_dz_sub = 4, esst_user_def = 5, esst_auto = 6
input_multem.atomic_vib_dim = [true, true, false];
input_multem.atomic_vib_seed = 300183;
input_multem.atomic_vib_nconf = 1;

input_multem.spec_rot_theta = 0; % final angle
input_multem.spec_rot_u_0 = [0 1 1]; % unitary vector			
input_multem.spec_rot_ctr_typ = 1; % 1: geometric center, 2: User define		
input_multem.spec_rot_ctr_p = [0 0 0]; % rotation point

input_multem.spec_bs_x = 10;
input_multem.spec_bs_y = 10;
input_multem.spec_bs_z = 10;
input_multem.spec_slic(1).sli_thick = 0.5;

occ = 1;
tag = 0;
charge = 0;
input_multem.spec_atoms = [29, 2, 2, 0.0, 0.8, 1.0, charge;29, 6, 2, 0.0, 0.8, 1.0, charge];
[input_multem.spec_atoms, input_multem.spec_bs_x, input_multem.spec_bs_y, lz] = graphene(1, 1.42, sqrt(0.5/(8*pi^2)));
input_multem.spec_slic(1).sli_thick = 0.5;

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

figure(1); clf;
plot(atoms(:, 2), atoms(:, 4), '*k');
set(gca, 'FontSize', 12, 'LineWidth', 1, 'PlotBoxAspectRatio', [1.25 1 1]);
title('Atomic positions');
ylabel('y', 'FontSize', 14);
xlabel('x', 'FontSize', 12);
axis equal;
axis([-2 input_multem.spec_bs_x+2 -5 input_multem.spec_bs_z + 5]);


for i = 1:nslice
    hold on;
    plot([-2 18], [Slice(i, 1) Slice(i, 1)], '-b', [-2 18], [Slice(i, 2) Slice(i, 2)], '-b');
    axis equal;
    axis([-2 input_multem.spec_bs_x+2 -5 input_multem.spec_bs_z + 5]);
end

for i = 1:nslice0
    hold on;
    plot([-2 input_multem.spec_bs_x+2], [Slice0(i, 1) Slice0(i, 1)], '-r', [-2 input_multem.spec_bs_x+2], [Slice0(i, 2) Slice0(i, 2)], '-r');
    axis equal;
    axis([-2 input_multem.spec_bs_x+2 -5 input_multem.spec_bs_z + 5]);
end