clear; clc;
addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])

%%%%%%%%%%%%%%%%%% Load multem default parameter %%%%%%%%$$%%%%%%%%%
input_multem = multem_default_values();          % Load default values;

%%%%%%%%%%%%%%%%%%%%% Set system configuration %%%%%%%%%%%%%%%%%%%%%
system_conf.precision = 1;                           % eP_Float = 1, eP_double = 2
system_conf.device = 2;                              % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 6;
system_conf.gpu_device = 0;

%%%%%%%%%%%%%%%%%%%% Set simulation experiment %%%%%%%%%%%%%%%%%%%%%
% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52,
% eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
input_multem.em_sim_typ = 11;

%%%%%%%%%%%%%% Electron-Specimen interaction model %%%%%%%%%%%%%%%%%
input_multem.elec_spec_interac_mod = 1;              % eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
input_multem.atomic_pot_parm_typ = 6;                 % eappt_doyle_0_4 = 1, eappt_peng_0_4 = 2, eappt_peng_0_12 = 3, eappt_kirkland_0_12 = 4, eappt_weickenmeier_0_12 = 5, eappt_lobato_0_12 = 6

%%%%%%%%%%%%%%%%%%%%%%% specimen slicing %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.spec_slic(1).typ = 1;              % esst_planes = 1, esst_dz_proj = 2, esst_planes_sub = 3, esst_dz_sub = 4, esst_auto = 5

%%%%%%%%%%%%%%% atomic vibrations model %%%%%%%%%%%%%%%%%%
input_multem.atomic_vib_mod = 3;                       % eavm_still_atom = 1, eavm_absorptive_pot = 2, eavm_frozen_phonon = 3, eavm_user_def = 4
input_multem.atomic_vib_coh_contrib = 0;
input_multem.atomic_vib_sgl_conf = 0;                 % 1: true, 0:false (extract single configuration)
input_multem.atomic_vib_nconf = 10;                      % true: specific phonon configuration, false: number of frozen phonon configurations
input_multem.atomic_vib_dim = [true, true, false];                       % phonon dimensions (xyz)
input_multem.atomic_vib_seed = 300183;                   % Random seed(frozen phonon)

%%%%%%%%%%%%%%%%%%%%%%% specimen information %%%%%%%%%%%%%%%%%%%%%%%
na = 8; nb = 8; nc = 5; ncu = 2; rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_lx...
, input_multem.spec_ly, input_multem.spec_lz...
, a, b, c, input_multem.spec_slic(1).dz] = Au001Crystal(na, nb, nc, ncu, rmsd_3d);

%%%%%%%%%%%%%%%%%%%%%% specimen thickness %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.thick_typ = 2;                     % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
input_multem.thick = c/2:c:1000;                   % Array of thickes (�)

%%%%%%%%%%%%%%%%%%%%%% x-y sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.nx = 512;
input_multem.ny = 512;
input_multem.bwl = 0;                            % Band-width limit, 1: true, 0:false

%%%%%%%%%%%%%%%%%%%% microscope parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.E_0 = 300;                          % Acceleration Voltage (keV)
input_multem.theta = 0.0;                        % Till ilumination (�)
input_multem.phi = 0.0;                          % Till ilumination (�)

%%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.illumination_model = 1;             % 1: coherente mode, 4: Numerical integration
input_multem.temporal_spatial_incoh = 1;         % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_m = 0;                  % Vortex momentum
input_multem.cond_lens_c_10 = 14.0312;            % Defocus (�)
input_multem.cond_lens_c_30 = 1e-03;            % Third order spherical aberration (mm)
input_multem.cond_lens_c_50 = 0.00;             % Fifth order spherical aberration (mm)
input_multem.cond_lens_c_12 = 0.0;             % Twofold astigmatism (�)
input_multem.cond_lens_phi_12 = 0.0;             % Azimuthal angle of the twofold astigmatism (�)
input_multem.cond_lens_c_23 = 0.0;             % Threefold astigmatism (�)
input_multem.cond_lens_phi_23 = 0.0;             % Azimuthal angle of the threefold astigmatism (�)
input_multem.cond_lens_inner_aper_ang = 0.0;   % Inner aperture (mrad)
input_multem.cond_lens_outer_aper_ang = 21.0;  % Outer aperture (mrad)

%%%%%%%%% defocus spread function %%%%%%%%%%%%
dsf_sigma = ilc_iehwgd_2_sigma(32); % from defocus spread to standard deviation
input_multem.cond_lens_ti_sigma = dsf_sigma;   % standard deviation (�)
input_multem.cond_lens_ti_npts = 5;         % # of integration points. It will be only used if illumination_model=4

%%%%%%%%%% source spread function %%%%%%%%%%%%
ssf_sigma = ilc_hwhm_2_sigma(0.45); % half width at half maximum to standard deviation
input_multem.cond_lens_si_sigma = ssf_sigma;  	% standard deviation: For parallel ilumination(�^-1); otherwise (�)
input_multem.cond_lens_si_rad_npts = 4;         % # of integration points. It will be only used if illumination_model=4

%%%%%%%%% zero defocus reference %%%%%%%%%%%%
input_multem.cond_lens_zero_def_typ = 1;   % eZDT_First = 1, eZDT_User_Define = 2
input_multem.cond_lens_zero_def_plane = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%STEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.scan_pat_typ = 2; % eST_Line = 1, eST_Area = 2
input_multem.scan_pat_pbc = 1;     % 1: true, 0:false (periodic boundary conditions)
input_multem.scan_pat_nsp = 16;
input_multem.scanning_x0 = 3*a;
input_multem.scanning_y0 = 3*b;
input_multem.scanning_xe = 4*a;
input_multem.scanning_ye = 4*b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Circular Detector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.detector.type = 1;  % eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3
input_multem.detector.cir(1).inner_ang = 40;  % Inner angle(mrad)
input_multem.detector.cir(1).outer_ang = 160; % Outer angle(mrad)

clear ilc_multem;
tic;
output_radial_detector = ilc_multem(system_conf, input_multem);
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Matrix Detector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This yields the same detector as the radial detector above, to check
% the consistency of the results

nxh = input_multem.nx/2;
nyh = input_multem.ny/2;
dgx = 1/input_multem.spec_lx;
dgy = 1/input_multem.spec_ly;

detector = zeros(input_multem.ny, input_multem.nx);
[gx, gy] = meshgrid((-nxh:1:(nxh-1))*dgx, (-nyh:1:(nyh-1))*dgy);
g = sqrt(gx.^2+gy.^2);

g_min = ilm_mrad_2_rAng(input_multem.E_0, input_multem.detector.cir(1).inner_ang);
g_max = ilm_mrad_2_rAng(input_multem.E_0, input_multem.detector.cir(1).outer_ang);
detector((g_min<=g)&(g<g_max)) = 1;

input_multem.detector.type = 3;  % eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3
input_multem.detector.matrix(1).fR = detector;

% Plot Matrix detector
figure(1); clf;
imagesc(detector);
caxis([0 1]);
axis image;
colormap gray;

%%%%%%%%%%%%%%%%%%%%%%%%%% run multem %%%%%%%%%%%%%%%%%%%%%%%%%%
clear ilc_multem;
tic;
output_matrix_detector = ilc_multem(system_conf, input_multem);
toc;

figure(2); clf;
for i=1:length(output_radial_detector.data)
    suptitle(['Thickness = ' num2str(i)])
    subplot(1, 3, 1);
    imagesc(output_radial_detector.data(i).image_tot(1).image);
    title('Radial Detector');
    axis image;
    colormap jet;
    colorbar;
    subplot(1, 3, 2);
    imagesc(output_matrix_detector.data(i).image_tot(1).image);
    title('Matrix Detector');
    axis image;
    colormap jet;
    colorbar;
    subplot(1, 3, 3);
    dI = output_radial_detector.data(i).image_tot(1).image - output_matrix_detector.data(i).image_tot(1).image;
    dI(dI<1e-6) = 0;
    imagesc(dI);
    title('Difference');
    axis image;
    colormap jet;
    colorbar;
    pause(0.5);
end
