% output_multislice = il_MULTEM(system_conf, input_multislice) perform TEM simulation
% 
% Energy filtered transmission electron microscopy (EFTEM) simulation
% 
% All parameters of the input_multislice structure are explained in multem_default_values()
% 
% Copyright 2017 Ivan Lobato <Ivanlh20@gmail.com>

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
input_multislice.simulation_type = 62;              % eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCI=42, eTEMST_EWFS=51, eTEMST_EWRS=52, eTEMST_EELS=61, eTEMST_EFTEM=62	
input_multislice.pn_model = 1;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multislice.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multislice.pn_dim = 110;                      % phonon dimensions (xyz)
input_multislice.pn_seed = 300183;                  % Random seed(frozen phonon)
input_multislice.pn_single_conf = 0;                % 1: true, 0:false (extract single configuration)
input_multislice.pn_nconf = 5;                      % true: phonon configuration, false: number of frozen phonon configurations

input_multislice.bwl = 0;

input_multislice.E_0 = 300;                         % Acceleration Voltage (keV)
input_multislice.theta = 0.0;                       % Tilt illumination (°)
input_multislice.phi = 0.0;                         % Tilt illumination (°)

na = 8; nb = 8; nc = 3; ncu = 2; rms3d = 0.085;

[input_multislice.spec_atoms, input_multislice.spec_lx...
, input_multislice.spec_ly, input_multislice.spec_lz...
, a, b, c, input_multislice.spec_dz] = SrTiO3001Crystal(na, nb, nc, ncu, rms3d);

input_multislice.nx = 1024; 
input_multislice.ny = 1024;

%%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.illumination_model = 1;            % 1: coherent mode, 2: Partial coherent mode, 3: transmission cross coefficient, 4: Numerical integration
input_multislice.temporal_spatial_incoh = 1;        % 1: Spatial and temporal, 2: Temporal, 3: Spatial

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.iw_type = 4;                       % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multislice.iw_psi = 0;                        % user define incident wave
input_multislice.iw_x = input_multislice.spec_lx/2;      % x position 
input_multislice.iw_y = input_multislice.spec_ly/2;      % y position

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.cond_lens_m = 0;                   % Vortex momentum
input_multislice.cond_lens_c_10 = 88.7414;             % Defocus (Å)
input_multislice.cond_lens_c_30 = 0.04;              % Third order spherical aberration (mm)
input_multislice.cond_lens_c_50 = 0.00;              % Fifth order spherical aberration (mm)
input_multislice.cond_lens_c_12 = 0.0;              % Twofold astigmatism (Å)
input_multislice.cond_lens_phi_12 = 0.0;              % Azimuthal angle of the twofold astigmatism (°)
input_multislice.cond_lens_c_23 = 0.0;              % Threefold astigmatism (Å)
input_multislice.cond_lens_phi_23 = 0.0;              %Azimuthal angle of the threefold astigmatism (°)
input_multislice.cond_lens_inner_aper_ang = 0.0;    % Inner aperture (mrad) 
input_multislice.cond_lens_outer_aper_ang = 21.0;   % Outer aperture (mrad)
input_multislice.cond_lens_sf = 32;                 % Defocus Spread (Å)
input_multislice.cond_lens_nsf = 10;                % Number of integration steps for the defocus Spread
input_multislice.cond_lens_beta = 0.2;              % Divergence semi-angle (mrad)
input_multislice.cond_lens_nbeta = 10;              % Number of integration steps for the divergence semi-angle
input_multislice.cond_lens_zero_defocus_type = 1;   % eZDT_First = 1, eZDT_User_Define = 2
input_multislice.cond_lens_zero_defocus_plane = 0;

%%%%%%%%%%%%%%%%%%%%%%%% Objective lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.obj_lens_m = 0;                    % Vortex momentum
input_multislice.obj_lens_c_10 = 0;                    % Defocus (Å)
input_multislice.obj_lens_c_30 = 0;                  % Third order spherical aberration (mm)
input_multislice.obj_lens_c_50 = 0.00;               % Fifth order spherical aberration (mm)
input_multislice.obj_lens_c_12 = 0.0;               % Twofold astigmatism (Å)
input_multislice.obj_lens_phi_12 = 0.0;               % Azimuthal angle of the twofold astigmatism (°)
input_multislice.obj_lens_c_23 = 0.0;               % Threefold astigmatism (Å)
input_multislice.obj_lens_phi_23 = 0.0;               %Azimuthal angle of the threefold astigmatism (°)
input_multislice.obj_lens_inner_aper_ang = 0.0;     % Inner aperture (mrad) 
input_multislice.obj_lens_outer_aper_ang = 0.0;     % Outer aperture (mrad)
input_multislice.obj_lens_sf = 32;                  % Defocus Spread (Å)
input_multislice.obj_lens_nsf = 10;                 % # integration steps for the defocus Spread. It will be only used if illumination_model=4
input_multislice.obj_lens_zero_defocus_type = 3;    % eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
input_multislice.obj_lens_zero_defocus_plane = 0;   % It will be only used if obj_lens_zero_defocus_type = eZDT_User_Define

input_multislice.eftem_E_loss = 532;                % Energy loss (eV)
input_multislice.eftem_m_selection = 3;             % selection rule
input_multislice.eftem_channelling_type = 1;        % eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 
input_multislice.eftem_collection_angle = 100;      % Collection half angle (mrad)
input_multislice.eftem_Z = 8;                       % atomic type

input_multislice.eftem_E_loss = 456;                % Energy loss (eV)
input_multislice.eftem_m_selection = 3;             % selection rule
input_multislice.eftem_channelling_type = 1;        % eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 
input_multislice.eftem_collection_angle = 100;      % Collection half angle (mrad)
input_multislice.eftem_Z = 22;                      % atomic type

input_multislice.eftem_E_loss = 1940;               % Energy loss (eV)
input_multislice.eftem_m_selection = 3;             % selection rule
input_multislice.eftem_channelling_type = 1;        % eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 
input_multislice.eftem_collection_angle = 50;       % Collection half angle (mrad)
input_multislice.eftem_Z = 38;                      % atomic type

clear il_MULTEM;
tic;
output_multislice = il_MULTEM(system_conf, input_multislice); 
toc;
clear il_MULTEM;

figure(1);
for i=1:length(output_multislice.data)
    imagesc(output_multislice.data(i).m2psi_tot);
    title(strcat('Total intensity -  Thick = ', num2str(i)));
    axis image;
    colormap gray;
    pause(0.25);
end

% cc = [1 0 1; 1 0 0; 0 0 1; 0 0 0];
% for i=1:4
%     input_multislice.eftem_channelling_type = i;
%     clear il_MULTEM;
%     tic;
%     [m2psi_tot] = il_MULTEM(system_conf, input_multislice); 
%     toc;
%     figure(1);
%     subplot(2, 1, 1);
%     hold on;
%     plot(m2psi_tot(129, :), 'color', cc(i, :));
%     subplot(2, 1, 2);
%     hold on;
%     plot(m2psi_tot(257, :), 'color', cc(i, :)); 
%     
%     figure(2);
%     subplot(2, 2, i);
%     imagesc(m2psi_tot);
%     title('Total intensity');
%     axis image;
%     colormap gray;
% end