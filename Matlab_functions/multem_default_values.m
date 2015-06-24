function [input_multislice] = multem_default_values()

input_multislice.precision = 1;                     % eP_Float = 1, eP_double = 2
input_multislice.device = 1;                        % eD_CPU = 1, eD_GPU = 2
input_multislice.cpu_ncores = 1;                    % Number of Cores CPU
input_multislice.cpu_nthread = 4;                   % Number of threads 
input_multislice.gpu_device = 0;                    % GPU device
input_multislice.gpu_nstream = 8;                   % Number of streams

input_multislice.simulation_type = 52;          % eST_STEM=11, eST_ISTEM=12, eST_CBED=21, eST_CBEI=22, eST_ED=31, eST_HRTEM=32, eST_PED=41, eST_HCI=42, eST_EWFS=51, eST_EWRS=52, eST_EELS=61, eST_EFTEM=62	
input_multislice.phonon_model = 1;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multislice.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multislice.fp_dim = 111;                      % phonon dimensions
input_multislice.fp_seed = 1983;                    % Random seed(frozen phonon)
input_multislice.fp_nconf = 1;                      % number of frozen phonon configurations
input_multislice.fp_iconf = 0;                      % phonon configuration

input_multislice.microscope_effect = 1;             % 1: Partial coherente mode, 2: transmission_fun cross coefficient
input_multislice.spatial_temporal_effect = 2;       % 1: Spatial and temporal, 2: Temporal, 3: Spatial

input_multislice.zero_defocus_type = 3;             % eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User = 4
input_multislice.zero_defocus_plane = 0;            % Zero defocus plane

input_multislice.thickness_type = 1;                % eTT_Whole_Specimen = 1, eTT_Through_Thickness = 2, eTT_Through_Planes = 3
input_multislice.thickness = 0;                     % Array of thicknesses

input_multislice.input_wave_type = 1;               % eIWT_Automatic = 1, eIWT_User_Define = 2
input_multislice.psi_0 = 0;                         % Input wave

input_multislice.fast_cal = 0;                      % 1: fast calculation(high memory consumption), 1 :normal mode(low memory consumption)
input_multislice.bwl = 0;                           % Band-width limit, 1: true, 0:false

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.E_0 = 300;                          % Acceleration Voltage (keV)
input_multislice.theta = 0.0;                       % Till ilumination (degrees)
input_multislice.phi = 0.0;                         % Till ilumination (degrees)

input_multislice.nx = 256; 
input_multislice.ny = 256;
input_multislice.lx = 10; 
input_multislice.ly = 10; 
input_multislice.dz = 0.25;

%%%%%%%%%%%%%%%%%%%%%%%% Microscope effects %%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.lens_m = 0;       %mm
input_multislice.lens_f = 0.0;     %Angs
input_multislice.lens_Cs3 = 0.00;	%mm
input_multislice.lens_Cs5 = 0.00;	%mm
input_multislice.lens_mfa2 = 0.0; 
input_multislice.lens_afa2 = 0.0; %(Angs, degrees)
input_multislice.lens_mfa3 = 0.0; 
input_multislice.lens_afa3 = 0.0; %(Angs, degrees)
input_multislice.lens_aobjl = 0.0; 
input_multislice.lens_aobju = 100000.0; %(mrad, mrad)
input_multislice.lens_sf = 32; 
input_multislice.lens_nsf = 10; % (Angs, number of steps)
input_multislice.lens_beta = 0.2; 
input_multislice.lens_nbeta = 10; %(mrad, half number of steps)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.scanning_type = 1; % eST_Line = 1, eST_Area = 2
input_multislice.scanning_ns = 10;
input_multislice.scanning_x0 = 0.0; 
input_multislice.scanning_y0 = 0.0;
input_multislice.scanning_xe = 4.078; 
input_multislice.scanning_ye = 4.078;

% Inner angle(mrad) and Outer angle(mrad)
input_multislice.det_cir(1).ang_inner = 60; 
input_multislice.det_cir(1).ang_outer = 180;
input_multislice.det_cir(2).ang_inner = 80; 
input_multislice.det_cir(2).ang_outer = 120;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CBED %%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.cbed_x0 = 0.0;      % x position 
input_multislice.cbed_y0 = 0.0;      % y position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CBED %%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.cbei_x0 = 0.0;      % x position 
input_multislice.cbei_y0 = 0.0;      % y position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HRTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.ped_nrot = 360;         % number of orientations
input_multislice.ped_theta = 3.0;        % Precession angle (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HCI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.hci_nrot = 360;         % number of orientations
input_multislice.hci_theta = 3.0;        % Precession angle (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%% EW Fourier Space %%%%%%%%%%%%%%%%%%%%%%
input_multislice.ewfs_convergent_beam = 0;     % 1: true, 0:false
input_multislice.ewfs_x0 = 0.0;                % x position 
input_multislice.ewfs_y0 = 0.0;                % y position
%%%%%%%%%%%%%%%%%%%%%%%%%% EW Real Space %%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.ewrs_convergent_beam = 0;     % 1: true, 0:false
input_multislice.ewrs_x0 = 0.0;                % x position 
input_multislice.ewrs_y0 = 0.0;                % y position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EFTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.eftem_E_loss = 0;              % Energy loss (eV)
input_multislice.eftem_m_selection = 3;         % selection rule
input_multislice.eftem_Z = 79;                  % atomic type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.eels_E_loss = 0;               % Energy loss (eV)
input_multislice.eels_m_selection = 3;          % selection rule
input_multislice.eels_collection_angle = 20;    % Collection half angle (mrad)
input_multislice.eels_Z = 79;                   % atomic type










