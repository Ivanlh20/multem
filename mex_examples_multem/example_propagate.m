clear all; clc;

input_multislice = multem_default_values();         % Load default values;

input_multislice.precision = 1;                     % eP_Float = 1, eP_double = 2
input_multislice.device = 1;                        % eD_CPU = 1, eD_GPU = 2
input_multislice.cpu_ncores = 1; 
input_multislice.cpu_nthread = 4; 
input_multislice.gpu_device = 0;
input_multislice.gpu_nstream = 8;

% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCI=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
input_multislice.simulation_type = 52;
input_multislice.phonon_model = 1;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multislice.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multislice.fp_dim = 110; 
input_multislice.fp_seed = 300183; 
input_multislice.fp_single_conf = 0;                % 1: true, 0:false
input_multislice.fp_nconf = 5;

input_multislice.illumination_model = 2;             % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient
input_multislice.temporal_spatial_incoh = 1;       % 1: Spatial and temporal, 2: Temporal, 3: Spatial

input_multislice.zero_defocus_type = 4;             % eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
input_multislice.zero_defocus_plane = 0;

input_multislice.bwl = 0;

input_multislice.E_0 = 300;                         % Acceleration Voltage (keV)
input_multislice.theta = 0.0;                       % Till ilumination (º)
input_multislice.phi = 0.0;                         % Till ilumination (º)

na = 4; nb = 4; nc = 10; ncu = 2; rms3d = 0.085;

[input_multislice.atoms, input_multislice.lx...
, input_multislice.ly, input_multislice.lz...
, a, b, c, input_multislice.dz] = Cu001Crystal(na, nb, nc, ncu, rms3d);

input_multislice.nx = 1024; 
input_multislice.ny = 1024;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.iw_type = 4;                      % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multislice.iw_psi = 0;                       % user define incident wave
input_multislice.iw_x = input_multislice.lx/2;     % x position 
input_multislice.iw_y = input_multislice.ly/2;     % y position

input_multislice.simulation_type = 52;              % eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCI=42, eTEMST_EWFS=51, eTEMST_EWRS=52, eTEMST_EELS=61, eTEMST_EFTEM=62
clear il_MULTEM;
tic;
output_multislice = il_MULTEM(input_multislice); 
toc;

input_multislice.obj_lens_f = +15;                      %Angs

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.iw_type = 3;                      % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multislice.iw_psi = output_multislice.data.psi_coh;  % user define incident wave

clear il_MULTEM;
tic;
output_propagate = il_propagate(input_multislice); 
toc;

figure(1);
subplot(2, 2, 1);
imagesc(abs(output_multislice.data.psi_coh).^2);
title('wave intensity');
axis image;
colormap gray;

subplot(2, 2, 2);
imagesc(angle(output_multislice.data.psi_coh));
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