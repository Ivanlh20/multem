% output_multem = ilc_multem(system_conf, input_multem) perform TEM simulation
% 
% STEM electron energy loss spectroscopy (EELS) simulation
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
input_multem.em_sim_typ = 61;

%%%%%%%%%%%%%% Electron-Specimen interaction model %%%%%%%%%%%%%%%%%
input_multem.elec_spec_interac_mod = 1; % eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
input_multem.atomic_pot_parm_typ = 6; % eappt_doyle_0_4 = 1, eappt_peng_0_4 = 2, eappt_peng_0_12 = 3, eappt_kirkland_0_12 = 4, eappt_weickenmeier_0_12 = 5, eappt_lobato_0_12 = 6

%%%%%%%%%%%%%%%%%%%%%%% specimen slicing %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.spec_slic(1).typ = 1; % esst_plns_proj = 1, esst_dz_proj = 2, esst_plns_sub = 3, esst_dz_sub = 4, esst_user_def = 5, esst_auto = 6

%%%%%%%%%%%%%%% atomic vibrations model %%%%%%%%%%%%%%%%%%
input_multem.atomic_vib_mod = 3; % eavm_still_atom = 1, eavm_absorptive_pot = 2, eavm_frozen_phonon = 3, eavm_user_def = 4
input_multem.atomic_vib_coh_contrib = 0;
input_multem.atomic_vib_sgl_conf = 0; % 1: true, 0:false (extract single configuration)
input_multem.atomic_vib_nconf = 5; % true: specific phonon configuration, false: number of frozen phonon configurations
input_multem.atomic_vib_dim = [true, true, false]; % phonon dimensions (xyz)
input_multem.atomic_vib_seed = 300183; % Random seed(frozen phonon)

%%%%%%%%%%%%%%%%%%%%%%% specimen information %%%%%%%%%%%%%%%%%%%%%%%
na = 8;nb = 8;nc = 5;ncu = 2;rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_bs_x...
, input_multem.spec_bs_y, input_multem.spec_bs_z...
, a, b, c, input_multem.spec_slic(1).sli_thick] = SrTiO3001_xtl(na, nb, nc, ncu, rmsd_3d);

%%%%%%%%%%%%%%%%%%%%%% specimen thickness %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.thick_typ = 1; % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
input_multem.thick = 0; % Array of thickes (Å)

%%%%%%%%%%%%%%%%%%%%%% x-y sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.nx = 1024;
input_multem.ny = 1024;
input_multem.bwl = 0; % Band-width limit, 1: true, 0:false

%%%%%%%%%%%%%%%%%%%% microscope parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.E_0 = 300; % Acceleration Voltage (keV)
input_multem.theta = 0.0; % Till ilumination (º)
input_multem.phi = 0.0; % Till ilumination (º)

%%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.illumination_model = 1; % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
input_multem.temporal_spatial_incoh = 1; % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.incdt_wav_typ = 4; % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.incdt_wav_psi = 0; % user define incident wave

%%%%%%%%%%%%%%%%%%%%%%%%%%% beam position %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.beam_pos = [input_multem.spec_bs_x/2;input_multem.spec_bs_y/2];

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_m = 0; % Vortex momentum
input_multem.cond_lens_c_10 = 88.7414; % Defocus (Å)
input_multem.cond_lens_c_30 = 0.04; % Third order spherical aberration (mm)
input_multem.cond_lens_c_50 = 0.00; % Fifth order spherical aberration (mm)
input_multem.cond_lens_c_12 = 0.0; % Twofold astigmatism (Å)
input_multem.cond_lens_phi_12 = 0.0; % Azimuthal angle of the twofold astigmatism (º)
input_multem.cond_lens_c_23 = 0.0; % Threefold astigmatism (Å)
input_multem.cond_lens_phi_23 = 0.0; % Azimuthal angle of the threefold astigmatism (º)
input_multem.cond_lens_inner_aper_ang = 0.0; % Inner aperture (mrad) 
input_multem.cond_lens_outer_aper_ang = 21.0; % Outer aperture (mrad)

%%%%%%%%% defocus spread function %%%%%%%%%%%%
dsf_sigma = ilc_iehwgd_2_sigma(32); % from defocus spread to standard deviation
input_multem.cond_lens_tp_inc_sigma = dsf_sigma; % standard deviation (Å)
input_multem.cond_lens_tp_inc_npts = 5; % # of integration points. It will be only used if illumination_model=4

%%%%%%%%%% source spread function %%%%%%%%%%%%
ssf_sigma = ilc_hwhm_2_sigma(0.45); % half width at half maximum to standard deviation
input_multem.obj_lens_ssf_sigma = ssf_sigma; % standard deviation: For parallel ilumination(Å^-1);otherwise (Å)
input_multem.obj_lens_ssf_npoints = 4; % # of integration points. It will be only used if illumination_model=4

%%%%%%%%% zero defocus reference %%%%%%%%%%%%
input_multem.cond_lens_zero_def_typ = 1; % eZDT_First = 1, eZDT_User_Define = 2
input_multem.cond_lens_zero_def_plane = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%% scanning area %%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.scan_pat_typ = 1; % eST_Line = 1, eST_Area = 2
input_multem.scan_pat_pbc = 1; % 1: true, 0:false (periodic boundary conditions)
input_multem.scan_pat_nsp = 10; % number of sampling points
input_multem.scan_pat_r_0 = [2*a;2.5*b]; % starting point (Å)
input_multem.scan_pat_r_e = [3*a;2.5*b]; % final point (Å)

input_multem.eels_E_loss = 532; % Energy loss (eV)
input_multem.eels_m_selection = 3; % selection rule
input_multem.eels_channelling_type = 1; % eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 
input_multem.eels_collection_angle = 100; % Collection half angle (mrad)
input_multem.eels_Z = 8; % atomic type

input_multem.eels_E_loss = 456; % Energy loss (eV)
input_multem.eels_m_selection = 3; % selection rule
input_multem.eels_channelling_type = 1; % eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 
input_multem.eels_collection_angle = 100; % Collection half angle (mrad)
input_multem.eels_Z = 22; % atomic type

input_multem.eels_E_loss = 1940; % Energy loss (eV)
input_multem.eels_m_selection = 3; % selection rule
input_multem.eels_channelling_type = 1; % eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 
input_multem.eels_collection_angle = 100; % Collection half angle (mrad)
input_multem.eels_Z = 38; % atomic type

clear ilc_multem;
tic;
output_multem = ilc_multem(system_conf, input_multem);
toc;

figure(1); clf;
for i=1:length(output_multem.data)
    plot(output_multem.data(i).image_tot(1).image);
    title(strcat('Thk = ', num2str(output_multem.thick(i))));
    pause(0.25);
end