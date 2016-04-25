% output_multislice = il_MULTEM(input_multislice) perform TEM simulation
% 
% Exit wave real space (EWRS) simulation using till illumination
% 
% All parameters of the input_multislice structure are explained in multem_default_values()
% 
% Copyright 2016 Ivan Lobato <Ivanlh20@gmail.com>

clear all; clc;

input_multislice = multem_default_values();         % Load default values;

input_multislice.precision = 1;                     % eP_Float = 1, eP_double = 2
input_multislice.device = 2;                        % eD_CPU = 1, eD_GPU = 2
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

input_multislice.fp_dim = 110;                      % phonon dimensions (xyz)
input_multislice.fp_seed = 300183;                  % Random seed(frozen phonon)
input_multislice.fp_single_conf = 0;                % 1: true, 0:false (extract single configuration)
input_multislice.fp_nconf = 5;                      % true: phonon configuration, false: number of frozen phonon configurations

input_multislice.bwl = 0;

input_multislice.E_0 = 300;                         % Acceleration Voltage (keV)
input_multislice.theta = 3.0;                       % Till ilumination (º)
input_multislice.phi = 0.0;                         % Till ilumination (º)

na = 4; nb = 4; nc = 15; ncu = 2; rms3d = 0.085;

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

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.cond_lens_m = 0;                  % Momentum of the vortex
input_multislice.cond_lens_f = 1110;               % Defocus (Å)
input_multislice.cond_lens_Cs3 = 3.3;              % Third order spherical aberration (mm)
input_multislice.cond_lens_Cs5 = 0.00;             % Fifth order spherical aberration (mm)
input_multislice.cond_lens_mfa2 = 0.0;             % Twofold astigmatism (Å)
input_multislice.cond_lens_afa2 = 0.0;             % Azimuthal angle of the twofold astigmatism (º)
input_multislice.cond_lens_mfa3 = 0.0;             % Threefold astigmatism (Å)
input_multislice.cond_lens_afa3 = 0.0;             % Azimuthal angle of the threefold astigmatism (º)
input_multislice.cond_lens_inner_aper_ang = 0.0;   % Inner aperture (mrad) 
input_multislice.cond_lens_outer_aper_ang = 7.50;  % Outer aperture (mrad)
input_multislice.cond_lens_sf = 32;                % Defocus Spread (Å)
input_multislice.cond_lens_nsf = 10;               % Number of integration steps for the defocus Spread
input_multislice.cond_lens_beta = 0.2;             % Divergence semi-angle (mrad)
input_multislice.cond_lens_nbeta = 10;             % Number of integration steps for the divergence semi-angle
input_multislice.cond_lens_zero_defocus_type = 1;  % eZDT_First = 1, eZDT_User_Define = 2
input_multislice.cond_lens_zero_defocus_plane = 0;

input_multislice.obj_lens_zero_defocus_type = 3;  % eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
input_multislice.obj_lens_zero_defocus_plane = 0;

input_multislice.thickness_type = 3;             % eTT_Whole_Specimen = 1, eTT_Through_Thickness = 2, eTT_Through_Slices = 3
input_multislice.thickness = 0:2*c:1000;           % Array of thicknesses

clear il_MULTEM;
tic;
output_multislice = il_MULTEM(input_multislice); 
toc;

figure(1);
for i=1:length(output_multislice.data)
    m2psi_coh = abs(output_multislice.data(i).psi_coh).^2;
    if(isfield(output_multislice.data, 'm2psi_tot'))
        m2psi_tot = output_multislice.data(i).m2psi_tot;  
    else
        m2psi_tot = m2psi_coh;
    end;    
    subplot(1, 2, 1);
    imagesc(m2psi_tot);
    title(strcat('Total intensity -  Thickness = ', num2str(i)));
    axis image;
    colormap gray;

    subplot(1, 2, 2);
    imagesc(m2psi_coh);
    title(strcat('Coherent intensity -  Thickness = ', num2str(i)));
    axis image;
    colormap gray;
    pause(0.25);
end;