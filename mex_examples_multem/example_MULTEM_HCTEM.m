% output_multem = ilc_multem(system_config, input_multem) perform TEM simulation
% 
% Hollow cone transmission electron microscopy (HCTEM) simulation
% 
% All parameters of the input_multem structure are explained in ilm_dflt_input_multem()
% 
% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>

clear;clc;
addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])

%%%%%%%%%%%%%%%%%% Load multem default parameter %%%%%%%%$$%%%%%%%%%
input_multem = ilm_dflt_input_multem(); % Load default values;

%%%%%%%%%%%%%%%%%%%%% Set system configuration %%%%%%%%%%%%%%%%%%%%%
system_config.precision = 1; % eP_Float = 1, eP_double = 2
system_config.device = 2; % eD_CPU = 1, eD_GPU = 2
system_config.cpu_n_thread = 1;
system_config.gpu_device = 0;

%%%%%%%%%%%%%%%%%%%% Set simulation experiment %%%%%%%%%%%%%%%%%%%%%
% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_STEM_EELS=61, eTEMST_ISTEM_EELS=62, eTEMST_EFTEMFS=71, eTEMST_EFTEMRS=72, eTEMST_ProbeFS=81, eTEMST_ProbeRS=82, eTEMST_PPFS=91, eTEMST_PPRS=92, eTEMST_TFFS=101, eTEMST_TFRS=102
input_multem.em_sim_typ = 42;

%%%%%%%%%%%%%% Electron-Specimen interaction model %%%%%%%%%%%%%%%%%
input_multem.elec_spec_interact_mod = 1; % eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
input_multem.atomic_pot_parm_typ = 6; % eappt_doyle_0_4 = 1, eappt_peng_0_4 = 2, eappt_peng_0_12 = 3, eappt_kirkland_0_12 = 4, eappt_weickenmeier_0_12 = 5, eappt_lobato_0_12 = 6

%%%%%%%%%%%%%%%%%%%%%%% specimen slicing %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.spec_slic(1).typ = 1; % esst_plns_proj = 1, esst_dz_proj = 2, esst_plns_sub = 3, esst_dz_sub = 4, esst_user_def = 5, esst_auto = 6

%%%%%%%%%%%%%%% atomic vibrations model %%%%%%%%%%%%%%%%%%
input_multem.atomic_vib_mod = 3; % eavm_still_atom = 1, eavm_absorptive_pot = 2, eavm_frozen_phonon = 3, eavm_user_def = 4
input_multem.atomic_vib_coh_contrib = 0;
input_multem.atomic_vib_sgl_conf = 0; % 1: true, 0:false (extract single configuration)
input_multem.atomic_vib_nconf = 10; % true: specific phonon configuration, false: number of frozen phonon configurations
input_multem.atomic_vib_dim = [true, true, false]; % phonon dimensions (xyz)
input_multem.atomic_vib_seed = 300183; % Random seed(frozen phonon)

%%%%%%%%%%%%%%%%%%%%%%% specimen information %%%%%%%%%%%%%%%%%%%%%%%
na = 16;nb = 16;nc = 20;ncu = 2;rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_bs_x...
, input_multem.spec_bs_y, input_multem.spec_bs_z...
, a, b, c, input_multem.spec_slic(1).sli_thick] = Cu001_xtl(na, nb, nc, ncu, rmsd_3d);

%%%%%%%%%%%%%%%%%%%%%% specimen thickness %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.thick_typ = 1; % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
input_multem.thick = c:c:1000; % Array of thickes (Å)

%%%%%%%%%%%%%%%%%%%%%% x-y sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.nx = 1024;
input_multem.ny = 1024;
input_multem.bwl = 0; % Band-width limit, 1: true, 0:false

%%%%%%%%%%%%%%%%%%%% microscope parameters %%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.E_0 = 300; % Acceleration Voltage (keV)
input_multem.theta = 0.0; % Till ilumination (º)
input_multem.phi = 0.0; % Till ilumination (º)

%%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.illum_mod = 1; % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
input_multem.illum_inc = 1; % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% source spread function %%%%%%%%%%%%
ssf_sigma = ilc_mrad_2_sigma(input_multem.E_0, 0.02); % mrad to standard deviation
input_multem.obj_lens_ssf_sigma = ssf_sigma; % standard deviation: For parallel ilumination(Å^-1);otherwise (Å)
input_multem.obj_lens_ssf_npoints = 4; % # of integration points. It will be only used if illum_mod=4

%%%%%%%%%%%%%%%%%%%%%%%% Objective lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.obj_lens_m = 0; % Vortex momentum
input_multem.obj_lens_c_10 = 20; % Defocus (Å)
input_multem.obj_lens_c_30 = 0.04; % Third order spherical aberration (mm)
input_multem.obj_lens_c_50 = 0.00; % Fifth order spherical aberration (mm)
input_multem.obj_lens_c_12 = 0.0; % Twofold astigmatism (Å)
input_multem.obj_lens_phi_12 = 0.0; % Azimuthal angle of the twofold astigmatism (º)
input_multem.obj_lens_c_23 = 0.0; % Threefold astigmatism (Å)
input_multem.obj_lens_phi_23 = 0.0; % Azimuthal angle of the threefold astigmatism (º)
input_multem.obj_lens_inner_aper_ang = 0.0; % Inner aperture (mrad) 
input_multem.obj_lens_outer_aper_ang = 25.0; % Outer aperture (mrad)

%%%%%%%%% defocus spread function %%%%%%%%%%%%
dsf_sigma = ilc_iehwgd_2_sigma(32); % from defocus spread to standard deviation
input_multem.obj_lens_tp_inc_sigma = dsf_sigma; % standard deviation (Å)
input_multem.obj_lens_tp_inc_npts = 5; % # of integration points. It will be only used if illum_mod=4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HCI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.hci_nrot = 30; % number of orientations
input_multem.hci_theta = 3.0; % Precession angle (degrees)

clear ilc_multem;
tic;
output_multem = ilc_multem(system_config, input_multem);
toc;

figure(1);
for i=1:length(output_multem.data)
    imagesc(output_multem.data(i).m2psi_tot);
    title(strcat('Total intensity -  Thick = ', num2str(output_multem.thick(i))));
    axis image;
    colormap gray;
    pause(0.25);
end