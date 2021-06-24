clear;clc;
addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])

input_multem = ilm_dflt_input_multem(); % Load default values;

system_conf.precision = 1; % eP_Float = 1, eP_double = 2
system_conf.device = 2; % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_n_thread = 4;
system_conf.gpu_device = 0;

% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_STEM_EELS=61, eTEMST_ISTEM_EELS=62, eTEMST_EFTEMFS=71, eTEMST_EFTEMRS=72, eTEMST_ProbeFS=81, eTEMST_ProbeRS=82, eTEMST_PPFS=91, eTEMST_PPRS=92, eTEMST_TFFS=101, eTEMST_TFRS=102
input_multem.em_sim_typ = 52;
input_multem.atomic_vib_mod = 1; % eavm_still_atom = 1, eavm_absorptive_pot = 2, eavm_frozen_phonon = 3, eavm_user_def = 4
input_multem.elec_spec_interac_mod = 1; % eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
input_multem.spec_slic(1).typ = 1; % esst_planes = 1, esst_dz_proj = 2, esst_planes_sub = 3, esst_dz_sub = 4, esst_auto = 5
input_multem.atomic_pot_parm_typ = 6; % eappt_doyle_0_4 = 1, eappt_peng_0_4 = 2, eappt_peng_0_12 = 3, eappt_kirkland_0_12 = 4, eappt_weickenmeier_0_12 = 5, eappt_lobato_0_12 = 6

input_multem.atomic_vib_dim = [true, true, false];
input_multem.atomic_vib_seed = 300183;
input_multem.atomic_vib_sgl_conf = 0; % 1: true, 0:false
input_multem.atomic_vib_nconf = 5;

input_multem.illumination_model = 2; % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
input_multem.temporal_spatial_incoh = 1; % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

input_multem.bwl = 0;

input_multem.E_0 = 300; % Acceleration Voltage (keV)
input_multem.theta = 0.0; % Till ilumination (º)
input_multem.phi = 0.0; % Till ilumination (º)

na = 4;nb = 4;nc = 10;ncu = 2;rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_bs_x...
, input_multem.spec_bs_y, input_multem.spec_bs_z...
, a, b, c, input_multem.spec_slic(1).dz] = Cu001_xtl(na, nb, nc, ncu, rmsd_3d);

input_multem.nx = 1024;
input_multem.ny = 1024;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.incdt_wav_typ = 4; % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.incdt_wav_psi = 0; % user define incident wave

%%%%%%%%%%%%%%%%%%%%%%%%%%% beam position %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.beam_pos = [input_multem.spec_bs_x/2;input_multem.spec_bs_y/2];

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_spt_inc_sigma = 0.2; % standard deviation: For parallel ilumination(Å^-1);otherwise (Å)

%%%%%%%%%%%%%%%%%%%%%%%% Objective lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.obj_lens_m = 0; % Vortex momentum
input_multem.obj_lens_c_10 = 15.836; % Defocus (Å)
input_multem.obj_lens_c_30 = 1e-03; % Third order spherical aberration (mm)
input_multem.obj_lens_c_50 = 0.00; % Fifth order spherical aberration (mm)
input_multem.obj_lens_c_12 = 0.0; % Twofold astigmatism (Å)
input_multem.obj_lens_phi_12 = 0.0; % Azimuthal angle of the twofold astigmatism (º)
input_multem.obj_lens_c_23 = 0.0; % Threefold astigmatism (Å)
input_multem.obj_lens_phi_23 = 0.0; % Azimuthal angle of the threefold astigmatism (º)
input_multem.obj_lens_inner_aper_ang = 0.0; % Inner aperture (mrad) 
input_multem.obj_lens_outer_aper_ang = 24.0; % Outer aperture (mrad)
input_multem.obj_lens_tp_inc_sigma = 32; % standard deviation (Å)
input_multem.obj_lens_tp_inc_npts = 10; % # integration steps for the defocus Spread. It will be only used if illumination_model=4
input_multem.obj_lens_zero_def_typ = 3; % eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
input_multem.obj_lens_zero_def_plane = 0; % It will be only used if obj_lens_zero_def_typ = eZDT_User_Define

% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_STEM_EELS=61, eTEMST_ISTEM_EELS=62, eTEMST_EFTEMFS=71, eTEMST_EFTEMRS=72, eTEMST_ProbeFS=81, eTEMST_ProbeRS=82, eTEMST_PPFS=91, eTEMST_PPRS=92, eTEMST_TFFS=101, eTEMST_TFRS=102
input_multem.em_sim_typ = 52;
clear ilc_multem;
tic;
output_multem_0 = ilc_multem(system_conf, input_multem);
toc;

% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_STEM_EELS=61, eTEMST_ISTEM_EELS=62, eTEMST_EFTEMFS=71, eTEMST_EFTEMRS=72, eTEMST_ProbeFS=81, eTEMST_ProbeRS=82, eTEMST_PPFS=91, eTEMST_PPRS=92, eTEMST_TFFS=101, eTEMST_TFRS=102
input_multem.em_sim_typ = 32;
clear ilc_multem;
tic;
output_multem_1 = ilc_multem(system_conf, input_multem);
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.incdt_wav_typ = 3; % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.incdt_wav_psi = output_multem_0.data.psi_coh; % user define incident wave
clear ilc_multem;
tic;
output_multem_2 = ilc_microscope_aberrations(system_conf, input_multem);
toc;

figure(1);
subplot(1, 3, 1);
imagesc(abs(output_multem_0.data.psi_coh).^2);
title('Total intensity');
axis image;
colormap gray;

subplot(1, 3, 2);
imagesc(output_multem_1.data.m2psi_tot);
title('Total intensity');
axis image;
colormap gray;

subplot(1, 3, 3);
imagesc(output_multem_2.m2psi);
title('Total intensity');
axis image;
colormap gray;