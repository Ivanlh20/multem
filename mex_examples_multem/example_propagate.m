% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
clear;clc;
addpath(['..', filesep, 'mex_bin'])

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

input_multem.atomic_vib_dim = [true, true, false];
input_multem.atomic_vib_seed = 300183;
input_multem.atomic_vib_sgl_conf = 0; % 1: true, 0:false
input_multem.atomic_vib_nconf = 5;

input_multem.illum_mod = 2; % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
input_multem.illum_inc = 1; % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

input_multem.bwl = 0;

input_multem.E_0 = 300; % Acceleration Voltage (keV)
input_multem.theta = 0.0; % Till ilumination (ยบ)
input_multem.phi = 0.0; % Till ilumination (ยบ)

na = 4;nb = 4;nc = 10;ncu = 2;rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_bs_x...
, input_multem.spec_bs_y, input_multem.spec_bs_z...
, a, b, c, input_multem.spec_slic(1).sli_thick] = Cu001_xtl(na, nb, nc, ncu, rmsd_3d);

input_multem.nx = 1024;
input_multem.ny = 1024;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.incdt_wav_typ = 4; % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.incdt_wav_psi = 0; % user define incident wave

%%%%%%%%%%%%%%%%%%%%%%%%%%% beam position %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.beam_pos = [input_multem.spec_bs_x/2;input_multem.spec_bs_y/2];

input_multem.em_sim_typ = 52; % eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, eTEMST_STEM_EELS=61, eTEMST_ISTEM_EELS=62, eTEMST_EFTEMFS=71, eTEMST_EFTEMRS=72
clear ilc_multem;
tic;
output_multem = ilc_multem(system_config, input_multem);
toc;

input_multem.obj_lens_c_10 = 10; %Angs

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.incdt_wav_typ = 3; % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.incdt_wav_psi = output_multem.data.psi_coh; % user define incident wave

tic;
output_propagate = ilc_propagate(system_config, input_multem);
toc;

figure(1);
subplot(2, 2, 1);
imagesc(abs(output_multem.data.psi_coh).^2);
title('wave intensity');
axis image;
colormap gray;

subplot(2, 2, 2);
imagesc(angle(output_multem.data.psi_coh));
title('Total intensity');
axis image;
colormap gray;

subplot(2, 2, 3);
imagesc(abs(output_propagate.psi).^2);
title('wave intensity');
axis image;
colormap gray;

subplot(2, 2, 4);
imagesc(angle(output_propagate.psi));
title('Total intensity');
axis image;
colormap gray;