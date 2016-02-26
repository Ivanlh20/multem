function [input_multislice] = multem_default_values()

input_multislice.precision = 1;                     % eP_Float = 1, eP_double = 2
input_multislice.device = 1;                        % eD_CPU = 1, eD_GPU = 2
input_multislice.cpu_ncores = 1;                    % Number of Cores CPU
input_multislice.cpu_nthread = 4;                   % Number of threads 
input_multislice.gpu_device = 0;                    % GPU device
input_multislice.gpu_nstream = 8;                   % Number of streams

% eST_STEM=11, eST_ISTEM=12, eST_CBED=21, eST_CBEI=22, eST_ED=31, eST_HRTEM=32, eST_PED=41, eST_HCI=42, eST_EWFS=51, eST_EWRS=52, 
% eST_EELS=61, eST_EFTEM=62, eST_ProbeFS=71, eST_ProbeRS=72, eST_PPFS=81, eST_PPRS=82,eST_TFFS=91, eST_TFRS=92
input_multislice.simulation_type = 52;             
input_multislice.phonon_model = 1;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multislice.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multislice.fp_dim = 111;                      % phonon dimensions (xyz)
input_multislice.fp_seed = 300183;                  % Random seed(frozen phonon)
input_multislice.fp_single_conf = 0;                % 1: true, 0:false (extract single configuration)
input_multislice.fp_nconf = 1;                      % true: phonon configuration, false: number of frozen phonon configurations

input_multislice.tm_active = 0;						% 1: true, 0:false
input_multislice.tm_theta = 90; 					% angle (º)
input_multislice.tm_u0 = [0 0 1]; 					% unitary vector			
input_multislice.tm_rot_point_type = 1; 			% 1: geometric center, 2: User define		
input_multislice.tm_p0 = [0 0 0];					% rotation point

input_multislice.microscope_effect = 2;             % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient
input_multislice.spatial_temporal_effect = 1;       % 1: Spatial and temporal, 2: Temporal, 3: Spatial

input_multislice.thickness_type = 1;                % eTT_Whole_Specimen = 1, eTT_Through_Thickness = 2, eTT_Through_Slices = 3
input_multislice.thickness = 0;                     % Array of thicknesses(Å)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.bwl = 0;                           % Band-width limit, 1: true, 0:false
input_multislice.operation_mode = 1;                % eOM_Normal = 1, eOM_Advanced = 2
input_multislice.coherent_contribution = 0;         % 1: true, 0:false
input_multislice.memory_size = 0;                   % memory size to be used(Mb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.E_0 = 300;                         % Acceleration Voltage (keV)
input_multislice.theta = 0.0;                       % Till ilumination (º)
input_multislice.phi = 0.0;                         % Till ilumination (º)

input_multislice.nx = 256;                          % number of pixels in x direction
input_multislice.ny = 256;                          % number of pixels in y direction
input_multislice.lx = 10;                           % simulation box length in x direction(Å)
input_multislice.ly = 10;                           % simulation box length in y direction(Å)
input_multislice.dz = 0.25;                         % slice thickness(Å)

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.iw_type = 4;                       % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multislice.iw_psi = 0;                        % user define incident wave (it is only active for iw_type=3)
input_multislice.iw_x = 0.0;                        % x position 
input_multislice.iw_y = 0.0;                        % y position

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.cond_lens_m = 0;                  % Momentum of the vortex
input_multislice.cond_lens_f = 14.0312;            % Defocus (Å)
input_multislice.cond_lens_Cs3 = 1e-03;            % Third order spherical aberration (mm)
input_multislice.cond_lens_Cs5 = 0.00;             % Fifth order spherical aberration (mm)
input_multislice.cond_lens_mfa2 = 0.0;             % Twofold astigmatism (Å)
input_multislice.cond_lens_afa2 = 0.0;             % Azimuthal angle of the twofold astigmatism (º)
input_multislice.cond_lens_mfa3 = 0.0;             % Threefold astigmatism (Å)
input_multislice.cond_lens_afa3 = 0.0;             % Azimuthal angle of the threefold astigmatism (º)
input_multislice.cond_lens_inner_aper_ang = 0.0;   % Inner aperture (mrad) 
input_multislice.cond_lens_outer_aper_ang = 21.0;  % Outer aperture (mrad)
input_multislice.cond_lens_sf = 32;                % Defocus Spread (Å)
input_multislice.cond_lens_nsf = 10;               % Number of integration steps for the defocus Spread
input_multislice.cond_lens_beta = 0.2;             % Divergence semi-angle (mrad)
input_multislice.cond_lens_nbeta = 10;             % Number of integration steps for the divergence semi-angle
input_multislice.cond_lens_zero_defocus_type = 1;  % eZDT_First = 1, eZDT_User_Define = 2
input_multislice.cond_lens_zero_defocus_plane = 0;

%%%%%%%%%%%%%%%%%%% condenser lens variable %%%%%%%%%%%%%%%%%%%%%
% 1: Momentum of the vortex, 2: Defocus (Å), 2: Third order spherical aberration (mm)
% 3: Third order spherical aberration (mm),  4: Fifth order spherical aberration (mm)
% 5: Twofold astigmatism (Å), 2: Defocus (Å), 6: Azimuthal angle of the twofold astigmatism (º)
% 7: Threefold astigmatism (Å),  8: Azimuthal angle of the threefold astigmatism (º)
% 9: Inner aperture (mrad), 2: Defocus (Å), 10: Outer aperture (mrad)
input_multislice.cdl_var_type = 0;                   % 0:off 1: m, 2: f, 3 Cs3, 4:Cs5, 5:mfa2, 6:afa2, 7:mfa3, 8:afa3, 9:inner_aper_ang , 10:outer_aper_ang
input_multislice.cdl_var = [-2 -1 0 1 2];            % variable array

%%%%%%%%%%%%%%%%%%%%%%%% Objective lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.obj_lens_m = 0;                  % Momentum of the vortex
input_multislice.obj_lens_f = 15.836;             % Defocus (Å)
input_multislice.obj_lens_Cs3 = 1e-03;            % Third order spherical aberration (mm)
input_multislice.obj_lens_Cs5 = 0.00;             % Fifth order spherical aberration (mm)
input_multislice.obj_lens_mfa2 = 0.0;             % Twofold astigmatism (Å)
input_multislice.obj_lens_afa2 = 0.0;             % Azimuthal angle of the twofold astigmatism (º)
input_multislice.obj_lens_mfa3 = 0.0;             % Threefold astigmatism (Å)
input_multislice.obj_lens_afa3 = 0.0;             % Azimuthal angle of the threefold astigmatism (º)
input_multislice.obj_lens_inner_aper_ang = 0.0;   % Inner aperture (mrad) 
input_multislice.obj_lens_outer_aper_ang = 24.0;  % Outer aperture (mrad)
input_multislice.obj_lens_sf = 32;                % Defocus Spread (Å)
input_multislice.obj_lens_nsf = 10;               % Number of integration steps for the defocus Spread
input_multislice.obj_lens_beta = 0.2;             % Divergence semi-angle (mrad)
input_multislice.obj_lens_nbeta = 10;             % Number of integration steps for the divergence semi-angle
input_multislice.obj_lens_zero_defocus_type = 3;  % eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
input_multislice.obj_lens_zero_defocus_plane = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% iSTEM/STEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.scanning_type = 1;         % eST_Line = 1, eST_Area = 2
input_multislice.scanning_periodic = 1;     % 1: true, 0:false (periodic boundary conditions)
input_multislice.scanning_ns = 10;          % number of sampling points
input_multislice.scanning_x0 = 0.0;         % x-starting point (Å)
input_multislice.scanning_y0 = 0.0;         % y-starting point (Å)
input_multislice.scanning_xe = 4.078;       % x-final point (Å)
input_multislice.scanning_ye = 4.078;       % y-final point (Å)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Detector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.detector.type = 1;  % eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3

input_multislice.detector.cir(1).inner_ang = 60;  % Inner angle(mrad) 
input_multislice.detector.cir(1).outer_ang = 180; % Outer angle(mrad)

input_multislice.detector.radial(1).x = 0;          % radial detector angle(mrad)
input_multislice.detector.radial(1).fx = 0;         % radial sensitivity value

input_multislice.detector.matrix(1).R = 0;          % 2D detector angle(mrad)
input_multislice.detector.matrix(1).fR = 0;         % 2D sensitivity value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.ped_nrot = 360;         % number of orientations
input_multislice.ped_theta = 3.0;        % Precession angle (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HCI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.hci_nrot = 360;         % number of orientations
input_multislice.hci_theta = 3.0;        % Precession angle (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EFTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.eftem_E_loss = 80;              % Energy loss (eV)
input_multislice.eftem_m_selection = 3;         % selection rule
input_multislice.eftem_collection_angle = 100;  % Collection half angle (mrad)
input_multislice.eftem_channelling_type = 1;    % eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 
input_multislice.eftem_Z = 79;                  % atomic type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.eels_E_loss = 80;               % Energy loss (eV)
input_multislice.eels_m_selection = 3;          % selection rule
input_multislice.eels_collection_angle = 100;   % Collection half angle (mrad)
input_multislice.eels_channelling_type = 1;     % eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 
input_multislice.eels_Z = 79;                   % atomic type
