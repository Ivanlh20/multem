clear;clc;

input_multem = ilm_dflt_input_multem(); % Load default values;

input_multem.atomic_vib_mod = 3; % eavm_still_atom = 1, eavm_absorptive_pot = 2, eavm_frozen_phonon = 3, eavm_user_def = 4
input_multem.elec_spec_interact_mod = 1; % eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
input_multem.spec_slic(1).typ = 1; % esst_plns_proj = 1, esst_dz_proj = 2, esst_plns_sub = 3, esst_dz_sub = 4, esst_user_def = 5, esst_auto = 6
input_multem.atomic_vib_dim = 111;
input_multem.atomic_vib_seed = 300183;
input_multem.atomic_vib_nconf = 3;

input_multem.spec_rot_theta = 0; % final angle
input_multem.spec_rot_u_0 = [1 0 0]; % unitary vector			
input_multem.spec_rot_ctr_typ = 1; % 1: geometric center, 2: User define		
input_multem.spec_rot_ctr_p = [0 0 0]; % rotation point

na = 4;nb = 4;nc = 20;ncu = 2;rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_bs_x...
, input_multem.spec_bs_y, input_multem.spec_bs_z...
, a, b, c, input_multem.spec_slic(1).sli_thick] = Au001_xtl(na, nb, nc, ncu, rmsd_3d);

% ilm_show_xtl(1, input_multem.spec_atoms);

input_multem.spec_slic(1).sli_thick = 5.0;

view
% get spec slicing
tic;
[atoms, Slice] = ilc_spec_slicing(input_multem);
toc;
[natoms, ~] = size(atoms);[nslice, ~] = size(Slice);

for i = 1:nslice
    figure(1); clf;
    i1 = Slice(i, 5);i2 = Slice(i, 6);ii = i1:1:i2;
    plot3(atoms(:, 2), atoms(:, 3), atoms(:, 4), '.k', atoms(ii, 2), atoms(ii, 3), atoms(ii, 4), 'or');
    set(gca, 'FontSize', 12, 'LineWidth', 1, 'PlotBoxAspectRatio', [1.25 1 1]);
    title('Atomic positions');
    ylabel('y', 'FontSize', 14);
    xlabel('x', 'FontSize', 12);
    axis equal;
    i2-i1+1
    view([1 0 0]);
    pause(0.1);
end

[size(input_multem.spec_atoms, 1), natoms, nslice]
[input_multem.spec_bs_x, input_multem.spec_bs_y, input_multem.spec_bs_z]