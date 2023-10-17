% output_multem = ilc_multem(system_config, input_multem) perform TEM simulation
% 
% Exit wave real space (EWRS) simulation
% 
% All parameters of the input_multem structure are explained in ilm_dflt_input_multem()
% 
% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>

clear;clc;

input_multem = ilm_dflt_input_multem(); % Load default values;

system_config.precision = 1; % eP_Float = 1, eP_double = 2
system_config.device = 2; % eD_CPU = 1, eD_GPU = 2
system_config.cpu_n_thread = 4;
system_config.gpu_device = 0;

% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_STEM_EELS=61, eTEMST_ISTEM_EELS=62, eTEMST_EFTEMFS=71, eTEMST_EFTEMRS=72, eTEMST_ProbeFS=81, eTEMST_ProbeRS=82, eTEMST_PPFS=91, eTEMST_PPRS=92, eTEMST_TFFS=101, eTEMST_TFRS=102
input_multem.em_sim_typ = 52;
input_multem.atomic_vib_mod = 1; % eavm_still_atom = 1, eavm_absorptive_pot = 2, eavm_frozen_phonon = 3, eavm_user_def = 4
input_multem.elec_spec_interact_mod = 1; % eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
input_multem.spec_slic(1).typ = 1; % esst_plns_proj = 1, esst_dz_proj = 2, esst_plns_sub = 3, esst_dz_sub = 4, esst_user_def = 5, esst_auto = 6
input_multem.atomic_pot_parm_typ = 6; % eappt_doyle_0_4 = 1, eappt_peng_0_4 = 2, eappt_peng_0_12 = 3, eappt_kirkland_0_12 = 4, eappt_weickenmeier_0_12 = 5, eappt_lobato_0_12 = 6

input_multem.atomic_vib_dim = [true, true, false]; % phonon dimensions (xyz)
input_multem.atomic_vib_seed = 300183; % Random seed(frozen phonon)
input_multem.atomic_vib_sgl_conf = 0; % 1: true, 0:false (extract single configuration)
input_multem.atomic_vib_nconf = 100; % true: phonon configuration, false: number of frozen phonon configurations

input_multem.spec_rot_theta = 45; % angle (º)
input_multem.spec_rot_u_0 = [1 0 0]; % unitary vector			
input_multem.spec_rot_ctr_typ = 1; % 1: geometric center, 2: User define		
input_multem.spec_rot_ctr_p = [0 0 0]; % rotation point

na = 8;nb = 8;nc = 8;ncu = 2;rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_bs_x...
, input_multem.spec_bs_y, input_multem.spec_bs_z...
, a, b, c, input_multem.spec_slic(1).sli_thick] = Au001_xtl(na, nb, nc, ncu, rmsd_3d);

input_multem.spec_bs_x = 100;
input_multem.spec_bs_y = 100;
input_multem.spec_bs_z = 100;

[input_multem.spec_atoms] = ilm_spec_recenter(input_multem.spec_atoms, input_multem.spec_bs_x, input_multem.spec_bs_y, input_multem.spec_bs_z);

% get spec slicing
[atoms, Slice] = ilc_spec_slicing(input_multem);

ilm_show_xtl(1, atoms);

[natoms, ~] = size(atoms);
[nslice, ~] = size(Slice);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
natoms = size(atoms, 1);
bb = zeros(natoms, 1);
d = 0.1;
ic = 0;
xy = [];
for ia=1:natoms
    if(bb(ia)<0.1)
        x = atoms(ia, 2);
        y = atoms(ia, 3);
        ii = find(sqrt((atoms(:, 2)-x).^2+(atoms(:, 3)-y).^2)<d);
        bb(ii) = 1;
        
        xy = [xy;[mean(atoms(ii, 2)), mean(atoms(ii, 3))]];
    end 
end

save('xy_projected.mat', 'xy');
figure(2);
plot(xy(:, 1), xy(:, 2), '*r');