% output_multislice = input_multem.ilc_multem perform TEM simulation
% All parameters of the input_multem structure are explained in ilm_dflt_input_multem()
% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>clear; clc;

addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])

input_multem = multem_input.parameters;         % Load default values;

input_multem.system_conf.precision = 1;                     % eP_Float = 1, eP_double = 2
input_multem.system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
input_multem.system_conf.cpu_nthread = 4;
input_multem.system_conf.gpu_device = 0;

% eTEMST_EWFS=51, eTEMST_EWRS=52
input_multem.simulation_type = 52;
input_multem.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multem.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4

input_multem.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multem.pn_model = 3;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multem.pn_dim = 110;                      % phonon dimensions
input_multem.pn_seed = 300183;                  % Random seed(frozen phonon)
input_multem.pn_single_conf = 1;                % 1: true, 0:false (extract single configuration)
input_multem.pn_nconf = 3;                      % true: phonon configuration, false: number of frozen phonon configurations

input_multem.thick_type = 1;                % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
input_multem.thick = 0;                     % Array of thickes

input_multem.bwl = 0;

input_multem.E_0 = 100;
input_multem.theta = 0.01;
input_multem.phi = 0.0;

na = 8; nb = 8; nc = 3; ncu = 2; rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_lx...
, input_multem.spec_ly, input_multem.spec_lz...
, a, b, c, input_multem.spec_dz] = Au001_xtl(na, nb, nc, ncu, rmsd_3d);

input_multem.thick_type = 1;             % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
input_multem.thick = 0:2*c:1000;         % Array of thickes

input_multem.nx = 2048;
input_multem.ny = 2048;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.iw_type = 1;       % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.iw_psi = read_psi_0_multem(input_multem.nx, input_multem.ny);    % user define incident wave
input_multem.iw_x = 0.5*input_multem.spec_lx;          % x position
input_multem.iw_y = 0.5*input_multem.spec_ly;          % y position

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_m = 0;                  % Vortex momentum
input_multem.cond_lens_c_10 = 1110;               % Defocus (�)
input_multem.cond_lens_c_30 = 3.3;              % Third order spherical aberration (mm)
input_multem.cond_lens_c_50 = 0.00;             % Fifth order spherical aberration (mm)
input_multem.cond_lens_c_12 = 0.0;             % Twofold astigmatism (�)
input_multem.cond_lens_phi_12 = 0.0;             % Azimuthal angle of the twofold astigmatism (�)
input_multem.cond_lens_c_23 = 0.0;             % Threefold astigmatism (�)
input_multem.cond_lens_phi_23 = 0.0;             % Azimuthal angle of the threefold astigmatism (�)
input_multem.cond_lens_inner_aper_ang = 0.0;   % Inner aperture (mrad)
input_multem.cond_lens_outer_aper_ang = 7.50;  % Outer aperture (mrad)
input_multem.cond_lens_ti_sigma = 32;                % standard deviation (�)
input_multem.cond_lens_ti_npts = 10;               % # of integration points. It will be only used if illumination_model=4
input_multem.cond_lens_si_sigma = 0.2;             % standard deviation: For parallel ilumination(�^-1); otherwise (�)
input_multem.cond_lens_si_rad_npts = 8;             % # of integration points. It will be only used if illumination_model=4
input_multem.cond_lens_zero_defocus_type = 1;  % eZDT_First = 1, eZDT_User_Define = 4
input_multem.cond_lens_zero_defocus_plane = 0;

%%%%%%%%%%%%%%%%%%%%%%%% Objective lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.obj_lens_m = 0;                  % Vortex momentum
input_multem.obj_lens_c_10 = 15.836;             % Defocus (�)
input_multem.obj_lens_c_30 = 1e-03;            % Third order spherical aberration (mm)
input_multem.obj_lens_c_50 = 0.00;             % Fifth order spherical aberration (mm)
input_multem.obj_lens_c_12 = 0.0;             % Twofold astigmatism (�)
input_multem.obj_lens_phi_12 = 0.0;             % Azimuthal angle of the twofold astigmatism (�)
input_multem.obj_lens_c_23 = 0.0;             % Threefold astigmatism (�)
input_multem.obj_lens_phi_23 = 0.0;             % Azimuthal angle of the threefold astigmatism (�)
input_multem.obj_lens_inner_aper_ang = 0.0;   % Inner aperture (mrad)
input_multem.obj_lens_outer_aper_ang = 24.0;  % Outer aperture (mrad)
input_multem.obj_lens_ti_sigma = 32;                % standard deviation (�)
input_multem.obj_lens_ti_npts = 10;               % # of integration points. It will be only used if illumination_model=4
input_multem.obj_lens_zero_defocus_type = 4;  % eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
input_multem.obj_lens_zero_defocus_plane = 5;

input_multem.output_area_ix_0 = 1;                             % x-starting pixel
input_multem.output_area_iy_0 = 1;                             % y-starting pixel
input_multem.output_area_ix_e = 1;                             % x-final pixel
input_multem.output_area_iy_e = 1;                             % y-final pixel

clear ilc_wave_function;
tic;
ouput_multislice = input_multem.ilc_wave_function;
toc;

figure(1);
for ithk=1:length(ouput_multislice.thick)
    psi_coh = flipud(ouput_multislice.data(ithk).psi_coh);

    subplot(1, 2, 1);
    imagesc(abs(psi_coh).^2);
    colormap gray;
    axis image;
    title(strcat('Intensity, thick = ', num2str(ithk)));
    subplot(1, 2, 2);
    imagesc(angle(psi_coh));
    colormap gray;
    axis image;
    title(strcat('Phase, thick = ', num2str(ithk)));
    pause(0.25);
end