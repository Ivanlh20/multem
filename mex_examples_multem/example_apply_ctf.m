% output_multislice = input_multem.ilc_multem perform TEM simulation
% inclusion of TEM aberrerations into a exitwve simulation
% All parameters of the input_multem structure are explained in ilm_dflt_input_multem()
% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>

clear; clc;
addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])

input_multem = multem_input.parameters;         % Load default values;

input_multem.system_conf.precision = 1;                     % eP_Float = 1, eP_double = 2
input_multem.system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
input_multem.system_conf.cpu_nthread = 4; 
input_multem.system_conf.gpu_device = 0;

% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
input_multem.simulation_type = 52;
input_multem.pn_model = 1;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multem.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multem.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multem.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multem.pn_dim = 110; 
input_multem.pn_seed = 300183; 
input_multem.pn_single_conf = 0;                % 1: true, 0:false
input_multem.pn_nconf = 5;

input_multem.illumination_model = 1;            % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
input_multem.temporal_spatial_incoh = 1;        % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

input_multem.bwl = 0;

input_multem.E_0 = 300;                         % Acceleration Voltage (keV)
input_multem.theta = 0.0;                       % Till ilumination (�)
input_multem.phi = 0.0;                         % Till ilumination (�)

na = 4; nb = 4; nc = 10; ncu = 2; rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_lx...
, input_multem.spec_ly, input_multem.spec_lz...
, a, b, c, input_multem.spec_dz] = Cu001_xtl(na, nb, nc, ncu, rmsd_3d);

input_multem.nx = 1024; 
input_multem.ny = 1024;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.iw_type = 4;                       % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.iw_psi = 0;                        % user define incident wave
input_multem.iw_x = input_multem.spec_lx/2;     % x position 
input_multem.iw_y = input_multem.spec_ly/2;     % y position

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_si_sigma = 0.2;          % standard deviation: For parallel ilumination(�^-1); otherwise (�)

%%%%%%%%%%%%%%%%%%%%%%%% Objective lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.obj_lens_m = 0;                    % Vortex momentum
input_multem.obj_lens_c_10 = 15.836;            % Defocus (�)
input_multem.obj_lens_c_30 = 1e-03;             % Third order spherical aberration (mm)
input_multem.obj_lens_c_50 = 0.00;              % Fifth order spherical aberration (mm)
input_multem.obj_lens_c_12 = 0.0;               % Twofold astigmatism (�)
input_multem.obj_lens_phi_12 = 0.0;             % Azimuthal angle of the twofold astigmatism (�)
input_multem.obj_lens_c_23 = 0.0;               % Threefold astigmatism (�)
input_multem.obj_lens_phi_23 = 0.0;             % Azimuthal angle of the threefold astigmatism (�)
input_multem.obj_lens_inner_aper_ang = 0.0;     % Inner aperture (mrad) 
input_multem.obj_lens_outer_aper_ang = 24.0;    % Outer aperture (mrad)
input_multem.obj_lens_ti_sigma = 32;            % standard deviation (�)
input_multem.obj_lens_ti_npts = 10;             % # integration steps for the defocus Spread. It will be only used if illumination_model=4
input_multem.obj_lens_zero_defocus_type = 3;    % eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
input_multem.obj_lens_zero_defocus_plane = 0;   % It will be only used if obj_lens_zero_defocus_type = eZDT_User_Define

% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
input_multem.simulation_type = 52;
clear ilc_multem;
tic;
output_multislice_0 = input_multem.ilc_multem; 
toc;

% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
input_multem.simulation_type = 32;
clear ilc_multem;
tic;
output_multislice_1 = input_multem.ilc_multem; 
toc;

%%%%%%%%%%%%%%%%%%%%% apply contrast transfer function %%%%%%%%%%%%%%%%%%%%
input_multem.iw_type = 3;                      % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.iw_psi = output_multislice_0.data.psi_coh;  % user define incident wave
clear ilc_multem;
tic;
output_multislice_2 = input_multem.ilc_apply_ctf; 
toc;

figure(1);
subplot(1, 3, 1);
imagesc(abs(output_multislice_0.data.psi_coh).^2);
title('Total intensity');
axis image;
colormap gray;

subplot(1, 3, 2);
imagesc(abs(output_multislice_1.data.m2psi_tot));
title('Total intensity');
axis image;
colormap gray;

subplot(1, 3, 3);
imagesc(abs(output_multislice_2.psi).^2);
title('Total intensity');
axis image;
colormap gray;