clear;clc;
addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])

input_multem = ilm_dflt_input_multem(); % Load default values;

system_config.precision = 1; % eP_Float = 1, eP_double = 2
system_config.device = 2; % eD_CPU = 1, eD_GPU = 2
system_config.cpu_n_thread = 4;
system_config.gpu_device = 0;

% eST_STEM=11, eST_ISTEM=12, eST_CBED=21, eST_CBEI=22, eST_ED=31, eST_HRTEM=32, eST_PED=41, eST_HCI=42, eST_EWFS=51, eST_EWRS=52, 
% eST_EELS=61, eST_EFTEM=62, eST_ProbeFS=71, eST_ProbeRS=72, eST_PPFS=81, eST_PPRS=82, eST_TFFS=91, eST_TFRS=92
input_multem.em_sim_typ = 61;
input_multem.atomic_vib_mod = 1; % eavm_still_atom = 1, eavm_absorptive_pot = 2, eavm_frozen_phonon = 3, eavm_user_def = 4
input_multem.elec_spec_interact_mod = 1; % eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
input_multem.spec_slic(1).typ = 1; % esst_plns_proj = 1, esst_dz_proj = 2, esst_plns_sub = 3, esst_dz_sub = 4, esst_user_def = 5, esst_auto = 6
input_multem.atomic_pot_parm_typ = 6; % eappt_doyle_0_4 = 1, eappt_peng_0_4 = 2, eappt_peng_0_12 = 3, eappt_kirkland_0_12 = 4, eappt_weickenmeier_0_12 = 5, eappt_lobato_0_12 = 6

input_multem.atomic_vib_dim = [true, true, false]; % phonon dimensions (xyz)
input_multem.atomic_vib_seed = 300183; % Random seed(frozen phonon)
input_multem.atomic_vib_sgl_conf = 0; % 1: true, 0:false (extract single configuration)
input_multem.atomic_vib_nconf = 5; % true: phonon configuration, false: number of frozen phonon configurations

input_multem.bwl = 0; % Band-width limit, 1: true, 0:false

input_multem.E_0 = 300; % Acceleration Voltage (keV)
input_multem.theta = 0.0; % Till ilumination (º)
input_multem.phi = 0.0; % Till ilumination (º)

na = 4;nb = 4;nc = 10;ncu = 2;rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_bs_x...
, input_multem.spec_bs_y, input_multem.spec_bs_z...
, a, b, c, input_multem.spec_slic(1).sli_thick] = SrTiO3001_xtl(na, nb, nc, ncu, rmsd_3d);

input_multem.nx = 512;
input_multem.ny = 512;

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
input_multem.cond_lens_tp_inc_sigma = 32; % standard deviation (Å)
input_multem.cond_lens_tp_inc_npts = 10; % # of integration points. It will be only used if illum_mod=4
input_multem.cond_lens_spt_inc_sigma = 0.2; % standard deviation: For parallel ilumination(Å^-1);otherwise (Å)
input_multem.cond_lens_spt_inc_rad_npts = 8; % # of integration points. It will be only used if illum_mod=4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%STEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
output_multem = ilc_multem(system_config, input_multem);
toc;
clear ilc_multem;

figure(1);
for i=1:length(output_multem.data)
    imagesc(output_multem.data(i).image_tot(1).image);
    title(strcat('Thk = ', num2str(i), ', det = ', num2str(j)));
    axis image;
    colormap gray;
    pause(0.25);
end

% cc = [1 0 1;1 0 0; 0 0 1; 0 0 0];
% for i=1:4
%  input_multem.eels_channelling_type = i;
%  clear ilc_multem;
%  tic;
%  [eels] = ilc_multem(system_config, input_multem);
%  toc;
%  figure(1);
%  hold on;
%  plot(eels, 'color', cc(i, :));
% end