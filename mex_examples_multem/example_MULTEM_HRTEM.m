% output_multislice = input_multem.ilc_multem perform TEM simulation
% High resolution transmission electron microscopy (HRTEM) simulation
% All parameters of the input_multem structure are explained in ilm_dflt_input_multem()
% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>

clear; clc;
addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])

%%%%%%%%%%%%%%%%%% Load multem default parameter %%%%%%%%$$%%%%%%%%%
input_multem = multem_input.parameters;          % Load default values;

%%%%%%%%%%%%%%%%%%%%% Set system configuration %%%%%%%%%%%%%%%%%%%%%
input_multem.system_conf.precision = 1;                           % eP_Float = 1, eP_double = 2
input_multem.system_conf.device = 2;                              % eD_CPU = 1, eD_GPU = 2
input_multem.system_conf.cpu_nthread = 1; 
input_multem.system_conf.gpu_device = 0;

%%%%%%%%%%%%%%%%%%%% Set simulation experiment %%%%%%%%%%%%%%%%%%%%%
% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
input_multem.simulation_type = 32;

%%%%%%%%%%%%%% Electron-Specimen interaction model %%%%%%%%%%%%%%%%%
input_multem.interaction_model = 1;              % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multem.potential_type = 6;                 % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

%%%%%%%%%%%%%%%%%%%%%%% Potential slicing %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.potential_slicing = 1;              % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4

%%%%%%%%%%%%%%% Electron-Phonon interaction model %%%%%%%%%%%%%%%%%%
input_multem.pn_model = 3;                       % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multem.pn_coh_contrib = 0;
input_multem.pn_single_conf = 0;                 % 1: true, 0:false (extract single configuration)
input_multem.pn_nconf = 10;                      % true: specific phonon configuration, false: number of frozen phonon configurations
input_multem.pn_dim = 110;                       % phonon dimensions (xyz)
input_multem.pn_seed = 300183;                   % Random seed(frozen phonon)

%%%%%%%%%%%%%%%%%%%%%%% Specimen information %%%%%%%%%%%%%%%%%%%%%%%
na = 16; nb = 16; nc = 5; ncu = 2; rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_lx...
, input_multem.spec_ly, input_multem.spec_lz...
, a, b, c, input_multem.spec_dz] = Cu001_xtl(na, nb, nc, ncu, rmsd_3d);

%%%%%%%%%%%%%%%%%%%%%% Specimen thickness %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.thick_type = 2;                     % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
input_multem.thick = c:c:1000;                   % Array of thickes (�)

%%%%%%%%%%%%%%%%%%%%%% x-y sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.nx = 1024;
input_multem.ny = 1024;
input_multem.bwl = 0;                            % Band-width limit, 1: true, 0:false

%%%%%%%%%%%%%%%%%%%% Microscope parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.E_0 = 300;                          % Acceleration Voltage (keV)
input_multem.theta = 0.0;                        % Till ilumination (�)
input_multem.phi = 0.0;                          % Till ilumination (�)

%%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.illumination_model = 1;             % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
input_multem.temporal_spatial_incoh = 1;         % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

%%%%%%%%%%%%%%%%%%%%%%%% Objective lens %%%%%%%%%%%%%%%%%%%%%%%%

input_multem.obj_lens_m = 0;                   % Vortex momentum
input_multem.obj_lens_c_10 = 20;               % Defocus (�)
input_multem.obj_lens_c_30 = 0.04;             % Third order spherical aberration (mm)
input_multem.obj_lens_c_50 = 0.00;             % Fifth order spherical aberration (mm)
input_multem.obj_lens_c_12 = 0.0;              % Twofold astigmatism (�)
input_multem.obj_lens_phi_12 = 0.0;            % Azimuthal angle of the twofold astigmatism (�)
input_multem.obj_lens_c_23 = 0.0;              % Threefold astigmatism (�)
input_multem.obj_lens_phi_23 = 0.0;            % Azimuthal angle of the threefold astigmatism (�)
input_multem.obj_lens_inner_aper_ang = 0.0;    % Inner aperture (mrad) 
input_multem.obj_lens_outer_aper_ang = 0.0;    % Outer aperture (mrad)

%%%%%%%%% defocus spread function %%%%%%%%%%%%
dsf_sigma = ilc_iehwgd_2_sigma(32); % from defocus spread to standard deviation
input_multem.obj_lens_ti_sigma = dsf_sigma;   % standard deviation (�)
input_multem.obj_lens_ti_npts = 5;         % # of integration points. It will be only used if illumination_model=4

%%%%%%%%%% source spread function %%%%%%%%%%%%
ssf_sigma = ilc_mrad_2_sigma(input_multem.E_0, 0.02);  % mrad to standard deviation
input_multem.obj_lens_si_sigma = ssf_sigma;          % standard deviation: For parallel ilumination(�^-1); otherwise (�)
input_multem.obj_lens_si_rad_npts = 4;                % # of integration points. It will be only used if illumination_model=4

%%%%%%%%% zero defocus reference %%%%%%%%%%%%
input_multem.obj_lens_zero_defocus_type = 1;   % eZDT_First = 1, eZDT_User_Define = 4
input_multem.obj_lens_zero_defocus_plane = 0;

clear ilc_multem;
tic;
output_multislice = input_multem.ilc_multem; 
toc;

figure(1);
for i=1:length(output_multislice.data)
    imagesc(output_multislice.data(i).m2psi_tot);
    title(strcat('Total intensity -  Thick = ', num2str(output_multislice.thick(i))));
    axis image;
    colormap gray;
    pause(0.25);
end