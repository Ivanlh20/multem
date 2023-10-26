% output_multislice = input_multem.ilc_multem perform TEM simulation
% STEM electron energy loss spectroscopy (EELS) simulation
% All parameters of the input_multem structure are explained in ilm_dflt_input_multem()
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

% eST_STEM=11, eST_ISTEM=12, eST_CBED=21, eST_CBEI=22, eST_ED=31, eST_HRTEM=32, eST_PED=41, eST_HCI=42, eST_EWFS=51, eST_EWRS=52, 
% eST_EELS=61, eST_EFTEM=62, eST_ProbeFS=71, eST_ProbeRS=72, eST_PPFS=81, eST_PPRS=82,eST_TFFS=91, eST_TFRS=92
input_multem.simulation_type = 61;
input_multem.pn_model = 1;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multem.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multem.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multem.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multem.pn_dim = 110;                      % phonon dimensions (xyz)
input_multem.pn_seed = 300183;                  % Random seed(frozen phonon)
input_multem.pn_single_conf = 0;                % 1: true, 0:false (extract single configuration)
input_multem.pn_nconf = 5;                      % true: phonon configuration, false: number of frozen phonon configurations

input_multem.bwl = 0;                           % Band-width limit, 1: true, 0:false

input_multem.E_0 = 300;                         % Acceleration Voltage (keV)
input_multem.theta = 0.0;                       % Till ilumination (�)
input_multem.phi = 0.0;                         % Till ilumination (�)

na = 4; nb = 4; nc = 10; ncu = 2; rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_lx...
, input_multem.spec_ly, input_multem.spec_lz...
, a, b, c, input_multem.spec_dz] = SrTiO3001_xtl(na, nb, nc, ncu, rmsd_3d);

input_multem.nx = 512; 
input_multem.ny = 512;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.iw_type = 4;                      % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.iw_psi = 0;                       % user define incident wave
input_multem.iw_x = input_multem.spec_lx/2;     % x position 
input_multem.iw_y = input_multem.spec_ly/2;     % y position

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_m = 0;                  % Vortex momentum
input_multem.cond_lens_c_10 = 88.7414;            % Defocus (�)
input_multem.cond_lens_c_30 = 0.04;             % Third order spherical aberration (mm)
input_multem.cond_lens_c_50 = 0.00;             % Fifth order spherical aberration (mm)
input_multem.cond_lens_c_12 = 0.0;             % Twofold astigmatism (�)
input_multem.cond_lens_phi_12 = 0.0;             % Azimuthal angle of the twofold astigmatism (�)
input_multem.cond_lens_c_23 = 0.0;             % Threefold astigmatism (�)
input_multem.cond_lens_phi_23 = 0.0;             % Azimuthal angle of the threefold astigmatism (�)
input_multem.cond_lens_inner_aper_ang = 0.0;   % Inner aperture (mrad) 
input_multem.cond_lens_outer_aper_ang = 21.0;  % Outer aperture (mrad)
input_multem.cond_lens_ti_sigma = 32;                % standard deviation (�)
input_multem.cond_lens_ti_npts = 10;               % # of integration points. It will be only used if illumination_model=4
input_multem.cond_lens_si_sigma = 0.2;             % standard deviation: For parallel ilumination(�^-1); otherwise (�)
input_multem.cond_lens_si_rad_npts = 8;             % # of integration points. It will be only used if illumination_model=4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%STEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.scanning_type = 1;             % eST_Line = 1, eST_Area = 2
input_multem.scanning_periodic = 1;         % 1: true, 0:false (periodic boundary conditions)
input_multem.scanning_ns = 10;              % number of sampling points
input_multem.scanning_x0 = 2*a;             % x-starting point (�) 
input_multem.scanning_y0 = 2.5*b;           % y-starting point (�)
input_multem.scanning_xe = 3*a;             % x-final point (�)
input_multem.scanning_ye = 2.5*b;           % y-final point (�)

input_multem.eels_E_loss = 532;             % Energy loss (eV)
input_multem.eels_m_selection = 3;          % selection rule
input_multem.eels_channelling_type = 1;     % eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 
input_multem.eels_collection_angle = 100;	% Collection half angle (mrad)
input_multem.eels_Z = 8;                    % atomic type

input_multem.eels_E_loss = 456;             % Energy loss (eV)
input_multem.eels_m_selection = 3;          % selection rule
input_multem.eels_channelling_type = 1;     % eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 
input_multem.eels_collection_angle = 100;	% Collection half angle (mrad)
input_multem.eels_Z = 22;                   % atomic type

input_multem.eels_E_loss = 1940;            % Energy loss (eV)
input_multem.eels_m_selection = 3;          % selection rule
input_multem.eels_channelling_type = 1;     % eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 
input_multem.eels_collection_angle = 100;   % Collection half angle (mrad)
input_multem.eels_Z = 38;                   % atomic type

clear ilc_multem;
tic;
output_multislice = input_multem.ilc_multem; 
toc;
clear ilc_multem;

figure(1);
for i=1:length(output_multislice.data)
    imagesc(output_multislice.data(i).image_tot(1).image);
    title(strcat('Thk = ', num2str(i), ', det = ', num2str(j)));
    axis image;
    colormap gray;
    pause(0.25);
end

% cc = [1 0 1; 1 0 0; 0 0 1; 0 0 0];
% for i=1:4
%     input_multem.eels_channelling_type = i;
%     clear ilc_multem;
%     tic;
%     [eels] = input_multem.ilc_multem; 
%     toc;
%     figure(1);
%     hold on;
%     plot(eels, 'color', cc(i, :));
% end