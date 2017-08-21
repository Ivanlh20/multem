clear all; clc;

input_multislice = multem_default_values();         % Load default values;

system_conf.precision = 1;                     % eP_Float = 1, eP_double = 2
system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_ncores = 1; 
system_conf.cpu_nthread = 4; 
system_conf.gpu_device = 0;
system_conf.gpu_nstream = 1;

% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCI=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
input_multislice.simulation_type = 52;
input_multislice.pn_model = 1;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multislice.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multislice.pn_dim = 110; 
input_multislice.pn_seed = 300183; 
input_multislice.pn_single_conf = 0;                % 1: true, 0:false
input_multislice.pn_nconf = 5;

input_multislice.illumination_model = 2;            % 1: coherent mode, 2: Partial coherent mode, 3: transmission cross coefficient, 4: Numerical integration
input_multislice.temporal_spatial_incoh = 1;        % 1: Spatial and temporal, 2: Temporal, 3: Spatial

input_multislice.bwl = 0;

input_multislice.E_0 = 300;                         % Acceleration Voltage (keV)
input_multislice.theta = 0.0;                       % Tilt illumination (º)
input_multislice.phi = 0.0;                         % Tilt illumination (º)

na = 4; nb = 4; nc = 10; ncu = 2; rms3d = 0.085;

[input_multislice.spec_atoms, input_multislice.spec_lx...
, input_multislice.spec_ly, input_multislice.spec_lz...
, a, b, c, input_multislice.spec_dz] = Cu001Crystal(na, nb, nc, ncu, rms3d);

input_multislice.nx = 1024; 
input_multislice.ny = 1024;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.iw_type = 4;                           % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multislice.iw_psi = 0;                            % user define incident wave
input_multislice.iw_x = input_multislice.spec_lx/2;     % x position 
input_multislice.iw_y = input_multislice.spec_ly/2;     % y position

input_multislice.simulation_type = 52;                  % eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCI=42, eTEMST_EWFS=51, eTEMST_EWRS=52, eTEMST_EELS=61, eTEMST_EFTEM=62
clear il_MULTEM;
tic;
output_multislice = il_MULTEM(system_conf, input_multislice); 
toc;

input_multislice.obj_lens_c_10 = 10;                      %Angs

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.iw_type = 3;                      % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multislice.iw_psi = output_multislice.data.psi_coh;  % user define incident wave

tic;
output_propagate = il_propagate(system_conf, input_multislice); 
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