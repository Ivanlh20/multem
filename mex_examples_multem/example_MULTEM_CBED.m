% output_multem = ilc_multem(system_conf, input_multem) perform TEM simulation
% 
% Convergent beam electron diffraction (CBED) simulation
% 
% All parameters of the input_multem structure are explained in ilm_dflt_input_multem()
% 
% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>

clear;clc;
addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])

%%%%%%%%%%%%%%%%%% Load multem default parameter %%%%%%%%$$%%%%%%%%%
input_multem = ilm_dflt_input_multem(); % Load default values;

%%%%%%%%%%%%%%%%%%%%% Set system configuration %%%%%%%%%%%%%%%%%%%%%
system_conf.precision = 1; % eP_Float = 1, eP_double = 2
system_conf.device = 2; % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_n_thread = 1;
system_conf.gpu_device = 0;

%%%%%%%%%%%%%%%%%%%% Set simulation experiment %%%%%%%%%%%%%%%%%%%%%
% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_STEM_EELS=61, eTEMST_ISTEM_EELS=62, eTEMST_EFTEMFS=71, eTEMST_EFTEMRS=72, eTEMST_ProbeFS=81, eTEMST_ProbeRS=82, eTEMST_PPFS=91, eTEMST_PPRS=92, eTEMST_TFFS=101, eTEMST_TFRS=102
input_multem.em_sim_typ = 21;

%%%%%%%%%%%%%% Electron-Specimen interaction model %%%%%%%%%%%%%%%%%
input_multem.interaction_model = 1; % eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
input_multem.atomic_pot_parm_typ = 6; % eappt_doyle_0_4 = 1, eappt_peng_0_4 = 2, eappt_peng_0_12 = 3, eappt_kirkland_0_12 = 4, eappt_weickenmeier_0_12 = 5, eappt_lobato_0_12 = 6

%%%%%%%%%%%%%%%%%%%%%%% specimen slicing %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.spec_slic(1).typ = 1; % esst_planes = 1, esst_dz_proj = 2, esst_planes_sub = 3, esst_dz_sub = 4, esst_auto = 5

%%%%%%%%%%%%%%% atomic vibrations model %%%%%%%%%%%%%%%%%%
input_multem.atomic_vib_model = 3; % eavm_still_atom = 1, eavm_absorptive_pot = 2, eavm_frozen_phonon = 3, eavm_user_def = 4
input_multem.atomic_vib_coh_contrib = 0;
input_multem.atomic_vib_sgl_conf = 0; % 1: true, 0:false (extract single configuration)
input_multem.atomic_vib_nconf = 10; % true: specific phonon configuration, false: number of frozen phonon configurations
input_multem.atomic_vib_dim = [true, true, false]; % phonon dimensions (xyz)
input_multem.atomic_vib_seed = 300183; % Random seed(frozen phonon)

%%%%%%%%%%%%%%%%%%%%%%% specimen information %%%%%%%%%%%%%%%%%%%%%%%
na = 8;nb = 8;nc = 40;ncu = 2;rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_bs_x...
, input_multem.spec_bs_y, input_multem.spec_bs_z...
, a, b, c, input_multem.spec_slic(1).dz] = Si001_xtl(na, nb, nc, ncu, rmsd_3d);

%%%%%%%%%%%%%%%%%%%%%% specimen thickness %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.thick_typ = 1; % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
input_multem.thick = 0; % Array of thickes (Å)

%%%%%%%%%%%%%%%%%%%%%% x-y sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.nx = 1024;
input_multem.ny = 1024;
input_multem.bwl = 0; % Band-width limit, 1: true, 0:false

%%%%%%%%%%%%%%%%%%%% microscope parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.E_0 = 100; % Acceleration Voltage (keV)
input_multem.theta = 0.0; % Till ilumination (º)
input_multem.phi = 0.0; % Till ilumination (º)

%%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.illumination_model = 1; % 1: coherente mode, 4: Numerical integration
input_multem.temporal_spatial_incoh = 1; % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.iw_type = 4; % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.iw_psi = read_psi_0_multem(input_multem.nx, input_multem.ny); % user define incident wave

%%%%%%%%%%%%%%%%%%%%%%%%%%% beam position %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.beam_pos = [input_multem.spec_bs_x/2;input_multem.spec_bs_y/2];

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_m = 0; % Vortex momentum
input_multem.cond_lens_c_10 = 1110; % Defocus (Å)
input_multem.cond_lens_c_30 = 3.3; % Third order spherical aberration (mm)
input_multem.cond_lens_c_50 = 0.00; % Fifth order spherical aberration (mm)
input_multem.cond_lens_c_12 = 0.0; % Twofold astigmatism (Å)
input_multem.cond_lens_phi_12 = 0.0; % Azimuthal angle of the twofold astigmatism (º)
input_multem.cond_lens_c_23 = 0.0; % Threefold astigmatism (Å)
input_multem.cond_lens_phi_23 = 0.0; % Azimuthal angle of the threefold astigmatism (º)
input_multem.cond_lens_inner_aper_ang = 0.0; % Inner aperture (mrad) 
input_multem.cond_lens_outer_aper_ang = 7; % Outer aperture (mrad)

%%%%%%%%% defocus spread function %%%%%%%%%%%%
dsf_sigma = ilc_iehwgd_2_sigma(32); % from defocus spread to standard deviation
input_multem.cond_lens_tp_inc_sigma = dsf_sigma; % standard deviation (Å)
input_multem.cond_lens_tp_inc_npts = 5; % # of integration points. It will be only used if illumination_model=4

%%%%%%%%%% source spread function %%%%%%%%%%%%
ssf_sigma = ilc_hwhm_2_sigma(0.45); % half width at half maximum to standard deviation
input_multem.cond_lens_spt_inc_sigma = ssf_sigma; % standard deviation: For parallel ilumination(Å^-1);otherwise (Å)
input_multem.cond_lens_spt_inc_rad_npts = 4; % # of integration points. It will be only used if illumination_model=4

%%%%%%%%% zero defocus reference %%%%%%%%%%%%
input_multem.cond_lens_zero_def_typ = 1; % eZDT_First = 1, eZDT_User_Define = 2
input_multem.cond_lens_zero_def_plane = 0;

clear ilc_multem;
tic;
output_multem = ilc_multem(system_conf, input_multem);
toc;

c = 1e5;
figure(1); clf;
for i=1:length(output_multem.data)   
    m2psi_tot = output_multem.data(i).m2psi_tot;
    m2psi_tot = log(1+c*m2psi_tot/max(m2psi_tot(:)));

    I_min = min(m2psi_tot(:));
    I_max = max(m2psi_tot(:));

    imagesc(output_multem.x, output_multem.y, m2psi_tot, [I_min I_max]);
    title(strcat('Total intensity -  Thick = ', num2str(output_multem.thick(i))));
    axis image;
    colormap hot;
    pause(0.5);
end