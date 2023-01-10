% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>

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
input_multem.theta = 0.0;                       % Till ilumination (º)
input_multem.phi = 0.0;                         % Till ilumination (º)

na = 4; nb = 4; nc = 10; ncu = 2; rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_lx...
, input_multem.spec_ly, input_multem.spec_lz...
, a, b, c, input_multem.spec_dz] = Cu001_xtl(na, nb, nc, ncu, rmsd_3d);

input_multem.nx = 1024; 
input_multem.ny = 1024;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.iw_type = 4;                           % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.iw_psi = 0;                            % user define incident wave
input_multem.iw_x = input_multem.spec_lx/2;     % x position 
input_multem.iw_y = input_multem.spec_ly/2;     % y position

input_multem.simulation_type = 52;                  % eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, eTEMST_EELS=61, eTEMST_EFTEM=62
clear ilc_multem;
tic;
output_multislice = input_multem.ilc_multem; 
toc;

input_multem.obj_lens_c_10 = 10;                      %Angs

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.iw_type = 3;                      % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.iw_psi = output_multislice.data.psi_coh;  % user define incident wave

input_multem.system_conf.precision = 2;                     % eP_Float = 1, eP_double = 2
input_multem.system_conf.device = 1;                        % eD_CPU = 1, eD_GPU = 2
input_multem.system_conf.cpu_nthread = 2; 
input_multem.system_conf.gpu_device = 0;
tic;
output_propagate = input_multem.ilc_propagate; 
toc;

figure(1);
subplot(2, 2, 1);
imagesc(abs(output_multislice.data.psi_coh).^2);
title('wave intensity');
axis image;
colormap gray;

subplot(2, 2, 2);
imagesc(angle(output_multislice.data.psi_coh));
title('phase');
axis image;
colormap gray;

subplot(2, 2, 3);
imagesc(abs(output_propagate.psi).^2);
title('wave intensity');
axis image;
colormap gray;

subplot(2, 2, 4);
imagesc(angle(output_propagate.psi));
title('angle');
axis image;
colormap gray;