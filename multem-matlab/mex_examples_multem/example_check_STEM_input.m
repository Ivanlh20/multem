% g_max (Angstrom^-1)
g = 12; 
% STEM Pixel Size (Angstrom)
dx = 0.35; 
% Size of STEM scan field (Angstrom)
scan_field = 20;
% Offset to account for probe size (Angstrom)
probe_size = 10;

%%%%%%%%%%%%%%%%%% Load multem default parameter %%%%%%%%%%%%%%%%%%%
input_multem = multem_input.parameters;          % Load default values;

%%%%%%%%%%%%%%%%%%%%% Set system configuration %%%%%%%%%%%%%%%%%%%%%
input_multem.system_conf.precision = 1;                           % eP_Float = 1, eP_double = 2
input_multem.system_conf.device = 2;                              % eD_CPU = 1, eD_GPU = 2
input_multem.system_conf.cpu_nthread = 12; 
input_multem.system_conf.gpu_device = 0;

%%%%%%%%%%%%%%%%%%%% Set simulation experiment %%%%%%%%%%%%%%%%%%%%%
% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
input_multem.simulation_type = 11;

%%%%%%%%%%%%%% Electron-Specimen interaction model %%%%%%%%%%%%%%%%%
input_multem.interaction_model = 1;              % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multem.potential_type = 6;                 % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

%%%%%%%%%%%%%%%%%%%%%%% Potential slicing %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.potential_slicing = 1;              % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4

%%%%%%%%%%%%%%% Electron-Phonon interaction model %%%%%%%%%%%%%%%%%%
input_multem.pn_model = 3;                       % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multem.pn_coh_contrib = 0;
input_multem.pn_single_conf = 0;                 % 1: true, 0:false (extract single configuration)
input_multem.pn_nconf = 3;                      % true: specific phonon configuration, false: number of frozen phonon configurations
input_multem.pn_dim = 110;                       % phonon dimensions (xyz)
input_multem.pn_seed = 100183;                   % Random seed(frozen phonon)

%%%%%%%%%%%%%%%%%%%%%%% Specimen information %%%%%%%%%%%%%%%%%%%%%%%
rmsd_3d = 0.0697;
[atoms, lx, ly, lz] = ilm_read_ap_xyz('Pt_twin.xyz', rmsd_3d);


atoms = ilm_center_spec(atoms, lx, ly, lz);
% centre point
p0 = [lx/2 ly/2 0]; 


lx_c = scan_field + 2 * probe_size;
b_crop = atoms(:,2) < (p0(1)-(lx_c/2)) | atoms(:,2) > (p0(1)+(lx_c/2)) | atoms(:,3) < (p0(2)-(lx_c/2)) | atoms(:,3) > (p0(2)+(lx_c/2));
atoms(b_crop,:) = [];

atoms = ilm_center_spec(atoms, lx_c, lx_c, lz);

lx = lx_c; 
ly = lx_c;

input_multem.spec_lx = lx;
input_multem.spec_ly = ly;
input_multem.spec_lz = lz;
input_multem.spec_atoms = atoms;

nx = tfm_pn_fact(g*2*lx,3);
ny = tfm_pn_fact(g*2*ly,3);

%%%%%%%%%%%%%%%%%%%%%% x-y sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.nx = nx;
input_multem.ny = ny;
input_multem.bwl = 0;                            % Band-width limit, 1: true, 0:false

%%%%%%%%%%%%%%%%%%%% Microscope parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.E_0 = 300;                          % Acceleration Voltage (keV)
input_multem.theta = 0.0;                        % Till ilumination (�)
input_multem.phi = 0.0;                          % Till ilumination (�)

%%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.illumination_model = 1;             % 1: coherente mode, 4: Numerical integration
input_multem.temporal_spatial_incoh = 1;         % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_m = 0;                  % Vortex momentum
input_multem.cond_lens_c_10 = -88.741;            % Defocus (�)
input_multem.cond_lens_c_30 = 0.04;            % Third order spherical aberration (mm)
input_multem.cond_lens_c_50 = 10;             % Fifth order spherical aberration (mm)
input_multem.cond_lens_c_12 = 0.0;             % Twofold astigmatism (�)
input_multem.cond_lens_phi_12 = 0.0;             % Azimuthal angle of the twofold astigmatism (�)
input_multem.cond_lens_c_23 = 0.0;             % Threefold astigmatism (�)
input_multem.cond_lens_phi_23 = 0.0;             % Azimuthal angle of the threefold astigmatism (�)
input_multem.cond_lens_inner_aper_ang = 0.0;   % Inner aperture (mrad) 
input_multem.cond_lens_outer_aper_ang = 21.1;  % Outer aperture (mrad)

%%%%%%%%% zero defocus reference %%%%%%%%%%%%
input_multem.cond_lens_zero_defocus_type = 1;   % eZDT_First = 1, eZDT_User_Define = 4
input_multem.cond_lens_zero_defocus_plane = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_multem.scanning_type = 2; % eST_Line = 1, eST_Area = 2
input_multem.scanning_periodic = 0;     % 1: true, 0:false (periodic boundary conditions)
input_multem.scanning_ns = ceil(scan_field/dx);
input_multem.scanning_x0 = (lx-scan_field)/2; 
input_multem.scanning_y0 = (ly-scan_field)/2;
input_multem.scanning_xe = (lx+scan_field)/2;
input_multem.scanning_ye = (ly+scan_field)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Detector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.detector.type = 1;  % eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3
input_multem.detector.cir(1).inner_ang = 60;  % Inner angle(mrad) 
input_multem.detector.cir(1).outer_ang = 190; % Outer angle(mrad)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tfm_check_stem_setup(input_multem, 1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear ilc_multem;
% tic;
% output_multislice = input_multem.ilc_multem; 
% toc;


